#!/usr/bin/env Rscript
# =============================================================================
# 40_mast_percell_dual_FH.R — Per-nucleus MAST discovery tier, frontal+hippo rebuild.
# Faithful copy of manuscript_7fig/scripts/40_mast_percell_dual.R with:
#   * input  = re-clustered F+Hippo atlas  atlas/NHD_FH_harmony.rds
#   * regions = c("Frontal","Hippo")   (OCC removed)
#   * output = NHD_frontal_hippo_rebuild/tables/mast_dual/
# Model, covariates, thresholds, Cliff's-delta, discovery call unchanged.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(MAST); library(future)
})
set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)
plan("sequential")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
DIR  <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
OUT  <- file.path(DIR, "tables/mast_dual"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

msg("=== loading F+Hippo atlas %s ===", format(Sys.time(), "%H:%M:%S"))
obj <- readRDS(file.path(DIR, "atlas/NHD_FH_harmony.rds"))

samp <- as.character(obj$Sample)
obj$matter <- factor(ifelse(grepl("_g$", samp), "grey",
                     ifelse(grepl("_w$", samp), "white", NA_character_)),
                     levels = c("grey","white"))
obj$hippo_area <- factor(ifelse(grepl("Hippo1$", samp), "area1",
                         ifelse(grepl("Hippo2$", samp), "area2", NA_character_)),
                         levels = c("area1","area2"))
obj$Condition   <- factor(as.character(obj$Condition), levels = c("CON","NHD"))
obj$log10nCount <- log10(obj$nCount_RNA)

# RNA is needed for RAW detection and effect size (see block below).
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"
msg("PrepSCTFindMarkers ... %s", format(Sys.time(), "%H:%M:%S"))
obj <- PrepSCTFindMarkers(obj, verbose = FALSE); invisible(gc())

genes_keep <- grep("^(MT-|RPL|RPS|MRPL|MRPS)", rownames(obj), value = TRUE, invert = TRUE)
msg("features kept (MT/ribo dropped): %d / %d", length(genes_keep), nrow(obj))

cliffs_delta <- function(seu, genes, grp) {
  if (!length(genes)) return(setNames(numeric(0), character(0)))
  expr <- as.matrix(GetAssayData(seu, assay = "SCT", layer = "data")[genes, , drop = FALSE])
  isN <- grp == "NHD"; n1 <- sum(isN); n2 <- sum(!isN)
  apply(expr, 1, function(x) { r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
}

cell_types <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
regions    <- c("Frontal","Hippo")     # OCC removed
MIN_PER_COND <- 10

run_one <- function(ct, rg) {
  label <- paste0(ct, "_", rg)
  cells <- WhichCells(obj, expression = new_annotation == ct & Region == rg)
  if (length(cells) < 30) { msg("  SKIP %-22s (<30 nuclei)", label); return(NULL) }
  sub <- subset(obj, cells = cells)
  tabc <- table(sub$Condition)
  if (any(tabc[c("CON","NHD")] < MIN_PER_COND) || length(unique(sub$Condition)) < 2) {
    msg("  SKIP %-22s (per-cond nuclei %s)", label, paste(tabc, collapse="/")); return(NULL) }
  Idents(sub) <- "Condition"
  block <- if (rg == "Hippo") "hippo_area" else "matter"
  has_block <- length(unique(na.omit(sub@meta.data[[block]]))) >= 2
  lv <- if (has_block) c(block, "log10nCount") else "log10nCount"
  m <- tryCatch(
    FindMarkers(sub, ident.1 = "NHD", ident.2 = "CON", test.use = "MAST",
                latent.vars = lv, assay = "SCT", features = genes_keep,
                min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi = FALSE, verbose = FALSE),
    error = function(e) { msg("  MAST err %s: %s", label, conditionMessage(e)); NULL })
  if (is.null(m) || !nrow(m)) return(NULL)
  m$gene <- rownames(m)
  # =========================================================================
  # Detection and effect size on raw RNA, not SCT.
  # SCT corrected counts fabricate detection for genes with zero raw UMIs: the
  # per-lane SCT model is aliased with Condition (NHD frontal lanes are ~48%
  # neurons vs 13% CON), so neuronal transcripts acquire a large expected value.
  # Verified: KCNIP4/RBFOX1/CSMD1/SYT1/OPCML have max raw count = 0 across all
  # nuclei yet were top discovery features, and they appeared as discovery
  # features in Astro and Oligo as well as Micro-PVM -- i.e. this affected every
  # cell type and therefore Figures 1/3/4/5, not just the microglia panels.
  # MAST still supplies the p-value on SCT (MAST-forward policy unchanged);
  # pct and Cliff's delta now come from RNA, delta is computed for every tested
  # gene (so null claims carry an effect size), and discovery additionally
  # requires real detection.
  # =========================================================================
  rna  <- GetAssayData(sub, assay = "RNA", layer = "counts")[m$gene, , drop = FALSE]
  isN  <- sub$Condition == "NHD"
  m$pct.1 <- round(Matrix::rowMeans(rna[,  isN, drop = FALSE] > 0), 3)
  m$pct.2 <- round(Matrix::rowMeans(rna[, !isN, drop = FALSE] > 0), 3)
  rnad <- GetAssayData(sub, assay = "RNA", layer = "data")[m$gene, , drop = FALSE]
  n1 <- sum(isN); n2 <- sum(!isN)
  m$cliffs_delta <- apply(as.matrix(rnad), 1, function(x) {
    r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
  m$max_pct_raw <- pmax(m$pct.1, m$pct.2)
  sig <- m$gene[which(m$p_val_adj < 0.05)]
  out <- m %>% transmute(
    comparison = label, cell_type = ct, region = rg, gene,
    avg_log2FC, p_val, p_val_adj, pct.1, pct.2, cliffs_delta, max_pct_raw,
    n_NHD = as.integer(tabc["NHD"]), n_CON = as.integer(tabc["CON"]),
    latent_vars = paste(lv, collapse = "+"),
    discovery = !is.na(p_val_adj) & p_val_adj < 0.05 &
                !is.na(cliffs_delta) & abs(cliffs_delta) >= 0.15 &
                max_pct_raw >= 0.10)
  stopifnot("discovery feature below the raw detection floor" =
              all(out$max_pct_raw[out$discovery] >= 0.10))
  write.csv(out, file.path(OUT, sprintf("MAST_%s.csv", label)), row.names = FALSE)
  msg("  %-22s nuclei=%d (NHD %d/CON %d) | sig(p<.05)=%d | discovery(|d|>=.15)=%d",
      label, ncol(sub), tabc["NHD"], tabc["CON"], length(sig), sum(out$discovery))
  invisible(gc()); out
}

all_res <- list()
for (rg in regions) for (ct in cell_types) {
  msg("=== %s_%s  START %s ===", ct, rg, format(Sys.time(), "%H:%M:%S"))
  all_res[[paste0(ct,"_",rg)]] <- run_one(ct, rg)
}
all_df <- bind_rows(all_res)
write.csv(all_df, file.path(OUT, "MAST_dual_discovery_all.csv"), row.names = FALSE)

msg("\n=== DISCOVERY COUNTS (per-nucleus MAST, |delta|>=0.15) ===")
summ <- all_df %>% filter(discovery) %>%
  group_by(comparison) %>%
  summarise(n_discovery = dplyr::n(),
            up_NHD = sum(avg_log2FC > 0), up_CON = sum(avg_log2FC < 0), .groups = "drop")
print(as.data.frame(summ), row.names = FALSE)
write.csv(summ, file.path(OUT, "MAST_discovery_counts.csv"), row.names = FALSE)
cat("\n=== DONE ===\n", file = stderr()); print(sessionInfo())
