#!/usr/bin/env Rscript
# =============================================================================
# 47_full_expressed_ranking_FH.R — Full-expressed GSEA preranking inputs, broad cell types, frontal+hippo rebuild. Faithful copy of scripts/47_full_expressed_ranking.R
# broad section only (new_annotation x {Frontal,Hippo}); the neocortical-subclass and
# hippo-de-novo-subtype caches are deferred to the neuron-figure stage, where the
# subtype taxonomy is re-derived on the NEW clustering (old neuron_subtype_map.rds
# does not match the re-clustered atlas). Recipe (FoldChange + midrank Cliff's delta,
# no MAST GLM, no p-values) unchanged.
# Output: NHD_frontal_hippo_rebuild/tables/mast_dual/cliffs_delta_full_<ct>_<rg>.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr); library(future) })
set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)
plan("sequential")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
DIR <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
OUT <- file.path(DIR, "tables/mast_dual"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

PCT_FLOOR <- 0.10; MIN_PER_COND <- 10; MIN_CELLS <- 30; GENE_CHUNK <- 1500
BROAD_TYPES <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
REGIONS     <- c("Frontal","Hippo")     # OCC removed
FORCE       <- identical(Sys.getenv("FORCE_CACHE"), "1")

old_caches <- list.files(OUT, pattern = "^cliffs_delta_full_.*\\.csv$", full.names = TRUE)
if (length(old_caches)) {
  bak_dir <- file.path(OUT, sprintf("_bak_pre_full_%s", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(bak_dir, showWarnings = FALSE, recursive = TRUE)
  ok <- file.copy(old_caches, bak_dir, overwrite = FALSE)
  msg("backed up %d existing cache(s) -> %s", sum(ok), bak_dir)
}

ATLAS <- file.path(DIR, "atlas/NHD_FH_harmony.rds")
stopifnot("MISSING atlas RDS" = file.exists(ATLAS))
msg("=== loading F+Hippo atlas %s ===", format(Sys.time(), "%H:%M:%S"))
obj <- readRDS(ATLAS)

# microglia-only primary: attach the split myeloid classes so the GSEA
# prerank can be built on microglia without PVM. run_one() already takes the metadata
# column as an argument, so this needs no change to the ranking logic itself. The
# pooled "Micro-PVM" rankings are still produced (cached) as the sensitivity tier.
.cls_f <- file.path(DIR, "tables", "micro_states", "immune_umap_classes_FH.csv")
if (!file.exists(.cls_f)) stop("MISSING ", .cls_f, " — run 2N_immune_UMAP_FH.R first")
.cls <- read.csv(.cls_f, stringsAsFactors = FALSE)
obj$myeloid_class <- .cls$class[match(colnames(obj), .cls$barcode)]
msg("myeloid classes attached: %s",
    paste(sprintf("%s=%d", names(table(obj$myeloid_class)),
                  as.integer(table(obj$myeloid_class))), collapse = " "))
obj$Condition <- factor(as.character(obj$Condition), levels = c("CON","NHD"))
DefaultAssay(obj) <- "SCT"
msg("PrepSCTFindMarkers ... %s", format(Sys.time(), "%H:%M:%S"))
obj <- PrepSCTFindMarkers(obj, verbose = FALSE); invisible(gc())
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"; obj <- NormalizeData(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"
genes_keep <- grep("^(MT-|RPL|RPS|MRPL|MRPS)", rownames(obj), value = TRUE, invert = TRUE)
msg("features kept (MT/ribo dropped): %d / %d", length(genes_keep), nrow(obj))

delta_chunk <- function(mat, isN) {
  n1 <- sum(isN); n2 <- sum(!isN)
  apply(mat, 1, function(x) { r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
}
cliffs_delta_full_vec <- function(sub, genes) {
  if (!length(genes)) return(setNames(numeric(0), character(0)))
  # Rank on RNA, not SCT. SCT corrected counts fabricate expression for
  # genes with zero raw UMIs (verified: KCNIP4/RBFOX1/CSMD1/SYT1/OPCML, max raw count
  # 0 across all nuclei, were top-ranked). A fabricated gene at the head of the
  # prerank drives the GSEA result, so this must be RNA.
  dat <- GetAssayData(sub, assay = "RNA", layer = "data")
  isN <- sub$Condition == "NHD"
  out <- numeric(length(genes)); names(out) <- genes
  idx <- split(seq_along(genes), ceiling(seq_along(genes) / GENE_CHUNK))
  for (ch in idx) { g <- genes[ch]; mm <- as.matrix(dat[g, , drop = FALSE]); out[g] <- delta_chunk(mm, isN); rm(mm) }
  invisible(gc()); out
}

run_one <- function(col, group_val, label) {
  cf <- file.path(OUT, sprintf("cliffs_delta_full_%s.csv", label))
  if (!FORCE && file.exists(cf)) { msg("  REUSE %-22s", label); return(invisible(NULL)) }
  cells <- rownames(obj@meta.data)[ obj@meta.data[[col]] == group_val & obj@meta.data$Region == group_region ]
  cells <- cells[!is.na(cells)]
  if (length(cells) < MIN_CELLS) { msg("  SKIP %-22s (<%d nuclei)", label, MIN_CELLS); return(invisible(NULL)) }
  sub  <- subset(obj, cells = cells)
  tabc <- table(sub$Condition)
  if (any(tabc[c("CON","NHD")] < MIN_PER_COND) || length(unique(sub$Condition)) < 2) {
    msg("  SKIP %-22s (per-cond nuclei %s)", label, paste(tabc, collapse = "/")); return(invisible(NULL)) }
  Idents(sub) <- "Condition"
  fc <- FoldChange(sub, ident.1 = "NHD", ident.2 = "CON", assay = "RNA")
  fc$gene <- rownames(fc)
  fc <- fc[fc$gene %in% genes_keep, , drop = FALSE]
  keep <- pmax(fc$pct.1, fc$pct.2) >= PCT_FLOOR
  fc   <- fc[keep, , drop = FALSE]
  if (!nrow(fc)) { msg("  SKIP %-22s (0 genes past pct floor)", label); return(invisible(NULL)) }
  d <- cliffs_delta_full_vec(sub, fc$gene)
  out <- data.frame(gene = fc$gene, cliffs_delta_full = unname(d[fc$gene]),
                    avg_log2FC = fc$avg_log2FC, pct.1 = fc$pct.1, pct.2 = fc$pct.2, stringsAsFactors = FALSE)
  tmp <- paste0(cf, ".tmp"); write.csv(out, tmp, row.names = FALSE)
  if (!file.rename(tmp, cf)) { file.copy(tmp, cf, overwrite = TRUE); unlink(tmp) }
  msg("  %-22s genes=%5d (NHD/CON %d/%d) | delta[min/med/max]=%.2f/%.2f/%.2f",
      label, nrow(out), tabc["NHD"], tabc["CON"],
      min(out$cliffs_delta_full), median(out$cliffs_delta_full), max(out$cliffs_delta_full))
  invisible(nrow(out))
}

counts <- list()
msg("\n=== BROAD comparisons (new_annotation x region) %s ===", format(Sys.time(), "%H:%M:%S"))
for (rg in REGIONS) { group_region <- rg
  for (ct in BROAD_TYPES) { label <- paste0(ct, "_", rg); n <- run_one("new_annotation", ct, label)
    if (!is.null(n)) counts[[label]] <- n } }

# Split myeloid prerank (microglia-only primary + PVM). Cached comparisons above are
# REUSEd, so this only computes the new ones.
msg("\n=== SPLIT myeloid comparisons (myeloid_class x region) %s ===",
    format(Sys.time(), "%H:%M:%S"))
for (rg in REGIONS) { group_region <- rg
  for (ct in c("Microglia","PVM")) { label <- paste0(ct, "_", rg)
    n <- run_one("myeloid_class", ct, label)
    if (!is.null(n)) counts[[label]] <- n } }

msg("\n=== FULL-EXPRESSED gene counts per comparison ===")
cnt_df <- data.frame(comparison = names(counts), n_genes_full = unlist(counts), row.names = NULL)
cnt_df <- cnt_df[order(cnt_df$comparison), ]
print(cnt_df, row.names = FALSE)
write.csv(cnt_df, file.path(OUT, "cliffs_delta_full_gene_counts.csv"), row.names = FALSE)
msg("\nwrote %d expanded caches | gene-count range %d..%d (median %d)",
    nrow(cnt_df), min(cnt_df$n_genes_full), max(cnt_df$n_genes_full), as.integer(median(cnt_df$n_genes_full)))
cat("\n=== DONE (broad-only; subclass+hippo-subtype deferred to neuron stage) ===\n", file = stderr())
