#!/usr/bin/env Rscript
# =============================================================================
# 47b_full_expressed_subclass_FH.R — Full-expressed GSEA preranking inputs for the neuron subclass / subtype layer, frontal+hippo rebuild. Faithful port of
#   the subclass + hippo-subtype sections of scripts/47_full_expressed_ranking.R.
# -----------------------------------------------------------------------------
# The broad-cell-type caches are already built by 47_full_expressed_ranking_FH.R;
# this script adds the neuron subclass/subtype comparisons that Fig 4 needs:
#
#   cortical subclass x Frontal  (label from scripts/_neuron_subclass_FH.R =
#     data/neuron_subtype_map_FH.rds; bona-fide neurons only)
# Frontal excitatory labels = Jorstad 2023 DLPFC transfer (116b),
#     which adds L4 IT; inhibitory = Azimuth. never atlas `predicted.subclass`.
#     L2/3 IT, L4 IT, L5 IT, L6 IT, L6 IT Car3, L6 CT, L6b, L5/6 NP, Pvalb, Sst, Vip,
#     Lamp5, Sncg.                        EXCLUDED (unpowered): L5 ET, Sst Chodl.
#
#   hippo DE-NOVO subtype x Hippo (neuron_subtype from data/neuron_subtype_map_FH.rds)
#     CA3, Subiculum, Inh-CGE (VIP), Inh-CGE (LAMP5), Inh-MGE.
#                                         excluded: Unresolved (not a subtype; CON=4).
#
# RECIPE (identical to the broad full-expressed ranking): FoldChange (all genes,
# no logfc.threshold) for avg_log2FC + pct.1 + pct.2, then min.pct>=0.10 floor,
# then the exact midrank Cliff's-delta over the kept genes (chunked so the dense
# gene x cell submatrix is bounded). no MAST GLM, no p-values — pure effect sizes,
# so pseudoreplication-safe for the single donor.
#
# Output:
#   tables/mast_dual/cliffs_delta_full_<label>.csv   (schema: gene,
#       cliffs_delta_full, avg_log2FC, pct.1, pct.2; label matches 45_FH's keys)
#   tables/mast_dual/hippo_subtype_meta.csv          (comparison, cell_type,
#       region, n_CON, n_NHD, family="hippo_subclass") -> powers 44b's tier gate.
#
# Heavy (one atlas load). Run detached; log tee'd.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr); library(future) })
set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)
plan("sequential")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot("Could not resolve PROJ root" = !is.na(PROJ), dir.exists(PROJ))
DIR <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
OUT <- file.path(DIR, "tables/mast_dual"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

PCT_FLOOR <- 0.10; MIN_PER_COND <- 10; MIN_CELLS <- 30; GENE_CHUNK <- 1500
FORCE     <- identical(Sys.getenv("FORCE_CACHE"), "1")
sanitize  <- function(x) gsub("[ /]", "_", x)

# Neocortical subclass census (Frontal); excluded unpowered: L5 ET, Sst Chodl.
# From the single-source helper (adds L4 IT) — 13 = 8 Ex + 5 Inh.
source(file.path(DIR, "scripts", "_neuron_subclass_FH.R"))   # SUBCLASS_MAP, CTX_SUBCLASSES, subclass_of()
SUBCLASSES  <- setdiff(CTX_SUBCLASSES, c("L5 ET","Sst Chodl"))
stopifnot("expected 13 subclasses incl. L4 IT" = length(SUBCLASSES) == 13 && "L4 IT" %in% SUBCLASSES)
SUB_REGIONS <- c("Frontal")
# Hippo de-novo subtypes (Hippo); excluded: Unresolved (not a subtype; CON=4).
HIPPO_SUBTYPES <- c("CA3","Subiculum","Inh-CGE (VIP)","Inh-CGE (LAMP5)","Inh-MGE")
# filename-safe keys (no spaces/parens/slashes; region suffix appended).
HIPPO_KEY <- c("CA3"="CA3","Subiculum"="Subiculum",
               "Inh-CGE (VIP)"="INCGE_VIP","Inh-CGE (LAMP5)"="INCGE_LAMP5",
               "Inh-MGE"="INMGE")

# Back up existing subclass/subtype caches (not the broad ones) before overwrite.
sub_labels_glob <- c(sprintf("cliffs_delta_full_%s_Frontal.csv", sanitize(SUBCLASSES)),
                     sprintf("cliffs_delta_full_%s_Hippo.csv", unname(HIPPO_KEY)))
old_caches <- file.path(OUT, sub_labels_glob); old_caches <- old_caches[file.exists(old_caches)]
if (length(old_caches)) {
  bak_dir <- file.path(OUT, sprintf("_bak_pre_subclass_%s", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(bak_dir, showWarnings = FALSE, recursive = TRUE)
  ok <- file.copy(old_caches, bak_dir, overwrite = FALSE)
  msg("backed up %d existing subclass/subtype cache(s) -> %s", sum(ok), bak_dir)
}

ATLAS <- file.path(DIR, "atlas/NHD_FH_harmony.rds")
stopifnot("MISSING atlas RDS" = file.exists(ATLAS))
msg("=== loading F+Hippo atlas %s ===", format(Sys.time(), "%H:%M:%S"))
obj <- readRDS(ATLAS)
stopifnot("MISSING metadata column new_annotation"     = "new_annotation"     %in% colnames(obj@meta.data),
          "MISSING metadata column Region"             = "Region"             %in% colnames(obj@meta.data),
          "MISSING metadata column Condition"          = "Condition"          %in% colnames(obj@meta.data))
obj$Condition <- factor(as.character(obj$Condition), levels = c("CON","NHD"))
# bona-fide-neuron flag (the map only carries neurons; kept explicit for the mask).
obj$is_neuron <- obj$new_annotation %in% c("Neuron_Ex", "Neuron_Inh")
# Cortical subclass label from the single-source map, by barcode.
obj$neuron_subclass <- subclass_of(colnames(obj))
stopifnot("map label on a non-neuron nucleus" = all(is.na(obj$neuron_subclass) | obj$is_neuron))
msg("subclass labels joined: %d Frontal neurons; source: %s",
    sum(!is.na(obj$neuron_subclass) & obj$Region == "Frontal"), CTX_SUBCLASS_SOURCE)

DefaultAssay(obj) <- "SCT"
msg("PrepSCTFindMarkers ... %s", format(Sys.time(), "%H:%M:%S"))
obj <- PrepSCTFindMarkers(obj, verbose = FALSE); invisible(gc())
# RNA must be joined + normalized on the full object before any subset,
# otherwise GetAssayData(layer="data") hands back raw counts. Detection, Cliff's delta
# and the pct floor are all computed on RNA (SCT corrected counts fabricate expression
# for zero-UMI genes -- see 40_mast_percell_dual_FH.R).
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"; obj <- NormalizeData(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"
genes_keep <- grep("^(MT-|RPL|RPS|MRPL|MRPS)", rownames(obj), value = TRUE, invert = TRUE)
msg("features kept (MT/ribo dropped): %d / %d", length(genes_keep), nrow(obj))

# ---------------------------------------------------------------------------
# Cliff's delta — exact midrank Mann-Whitney, chunked over genes so the dense
# (genes x cells) submatrix is bounded (matches the broad full-expressed helper).
# ---------------------------------------------------------------------------
delta_chunk <- function(mat, isN) {
  n1 <- sum(isN); n2 <- sum(!isN)
  apply(mat, 1, function(x) { r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
}
cliffs_delta_full_vec <- function(sub, genes) {
  if (!length(genes)) return(setNames(numeric(0), character(0)))
  dat <- GetAssayData(sub, assay = "RNA", layer = "data")
  isN <- sub$Condition == "NHD"
  out <- numeric(length(genes)); names(out) <- genes
  idx <- split(seq_along(genes), ceiling(seq_along(genes) / GENE_CHUNK))
  for (ch in idx) { g <- genes[ch]; mm <- as.matrix(dat[g, , drop = FALSE]); out[g] <- delta_chunk(mm, isN); rm(mm) }
  invisible(gc()); out
}

# ---------------------------------------------------------------------------
# One comparison: subset by `cell_mask` (logical over obj cells) -> gate ->
# FoldChange (all genes) -> pct floor -> full-coverage Cliff's delta -> write.
# `cell_mask` is passed explicitly (not via parent scope) so the cortical-neuron
# restriction and the hippo-subtype join are both expressed as clean masks.
# ---------------------------------------------------------------------------
run_one <- function(cell_mask, label) {
  cf <- file.path(OUT, sprintf("cliffs_delta_full_%s.csv", label))
  # Reuse only a cache newer than both the atlas and the subclass map
  if (!FORCE && file.exists(cf) && file.mtime(cf) > max(file.mtime(ATLAS), NEURON_MAP_MTIME)) { msg("  REUSE %-22s (cache newer than atlas + map; FORCE_CACHE=1 to recompute)", label); return(invisible(NULL)) }
  cells <- rownames(obj@meta.data)[cell_mask & !is.na(cell_mask)]
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

# --- cortical subclass x Frontal (map label, neurons only) -------------------
msg("\n=== SUBCLASS comparisons (neuron_subclass x Frontal, neurons only) %s ===", format(Sys.time(), "%H:%M:%S"))
msg("CENSUS excludes (unpowered, not run): L5 ET, Sst Chodl")
for (rg in SUB_REGIONS) for (sc in SUBCLASSES) {
  label <- paste0(sanitize(sc), "_", rg)
  # Key on the map label (neuron_subclass), never on predicted.subclass
  mask  <- obj$is_neuron & !is.na(obj$neuron_subclass) & obj$neuron_subclass == sc & obj$Region == rg
  n <- run_one(mask, label)
  if (!is.null(n)) counts[[label]] <- n
}

# --- hippo DE-NOVO subtype x Hippo (neuron_subtype from map) ------------------
HIP_MAP <- file.path(DIR, "data/neuron_subtype_map_FH.rds")
if (!file.exists(HIP_MAP)) {
  msg("\n=== HIPPO subtypes SKIPPED — MISSING %s ===", HIP_MAP)
} else {
  msg("\n=== HIPPO de-novo subtype comparisons (neuron_subtype x Hippo) %s ===", format(Sys.time(), "%H:%M:%S"))
  # HIPPO_SUBTYPES is an allow-list, so the unresolved cluster is excluded by
  # construction rather than by name.
  msg("EXCLUDE (not a subtype; failed the private-marker gate, CON=4): Unresolved")
  hmap <- readRDS(HIP_MAP)
  stopifnot("hippo map missing required columns" =
              all(c("barcode","Region","neuron_subtype","Condition") %in% colnames(hmap)))
  hmap <- hmap[hmap$Region == "Hippo" & hmap$neuron_subtype %in% HIPPO_SUBTYPES, , drop = FALSE]
  # join onto atlas by barcode: NA for every nucleus not in the hippo de-novo map.
  obj$neuron_subtype <- hmap$neuron_subtype[match(rownames(obj@meta.data), hmap$barcode)]
  n_join <- sum(!is.na(obj$neuron_subtype))
  msg("  joined %d hippo-subtype labels onto atlas (of %d kept map rows)", n_join, nrow(hmap))
  stopifnot("hippo-subtype join produced 0 matches — barcode format mismatch" = n_join > 0)

  hippo_meta <- list()
  for (st in HIPPO_SUBTYPES) {
    label <- paste0(HIPPO_KEY[[st]], "_Hippo")
    # per-condition counts are cheap (meta only) -> always refresh the sidecar,
    # even when the expensive delta cache is REUSEd. Powers 44b's tier gate.
    cc <- table(factor(obj$Condition[!is.na(obj$neuron_subtype) & obj$neuron_subtype == st],
                        levels = c("CON","NHD")))
    hippo_meta[[label]] <- data.frame(
      comparison = label, cell_type = st, region = "Hippo",
      n_CON = as.integer(cc["CON"]), n_NHD = as.integer(cc["NHD"]),
      family = "hippo_subclass", stringsAsFactors = FALSE)
    mask <- !is.na(obj$neuron_subtype) & obj$neuron_subtype == st
    n <- run_one(mask, label)
    if (!is.null(n)) counts[[label]] <- n
  }
  hippo_meta_df <- do.call(rbind, hippo_meta)
  hmp <- file.path(OUT, "hippo_subtype_meta.csv")
  tmp <- paste0(hmp, ".tmp"); write.csv(hippo_meta_df, tmp, row.names = FALSE)
  if (!file.rename(tmp, hmp)) { file.copy(tmp, hmp, overwrite = TRUE); unlink(tmp) }
  msg("  wrote %s (%d subtypes)", hmp, nrow(hippo_meta_df))
  print(hippo_meta_df[, c("comparison","cell_type","n_CON","n_NHD")], row.names = FALSE)
}

# ---------------------------------------------------------------------------
msg("\n=== FULL-EXPRESSED gene counts per subclass/subtype comparison ===")
cnt_df <- data.frame(comparison = names(counts), n_genes_full = unlist(counts), row.names = NULL)
cnt_df <- cnt_df[order(cnt_df$comparison), ]
print(cnt_df, row.names = FALSE)
write.csv(cnt_df, file.path(OUT, "cliffs_delta_full_subclass_gene_counts.csv"), row.names = FALSE)
if (nrow(cnt_df))
  msg("\nwrote %d subclass/subtype caches | gene-count range %d..%d (median %d)",
      nrow(cnt_df), min(cnt_df$n_genes_full), max(cnt_df$n_genes_full), as.integer(median(cnt_df$n_genes_full)))
cat("\n=== DONE ===\n", file = stderr()); print(sessionInfo())
