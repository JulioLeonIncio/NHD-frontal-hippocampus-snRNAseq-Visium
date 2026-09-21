#!/usr/bin/env Rscript
# =============================================================================
# 92_micro_only_MAST_FH.R — Microglia-ONLY differential expression, run with the identical model to the pooled Micro-PVM run that produced ST1.
# -----------------------------------------------------------------------------
# Why this exists.  The manuscript states that
# pooling perivascular-macrophage-like (CD163+/F13A1+) nuclei with microglia into
# the Azimuth "Micro-PVM" class does not cost gene-level sensitivity.  ST17 shows
# that claim; this script supplies its microglia-only half.
#
# Why the EXISTING 2Q2 output could not be reused.  tables/micro_states/
# pvm_sensitivity_genes_FH.csv answers a similar question, but its
# Cliff's delta is computed on the SCT `data` layer.
# ST1's deltas are RNA-based; the 2Q2 deltas are not.
# Putting the two side by side in one table would have produced a `delta_difference`
# column that measures the assay change, not the cell-selection change.  So the
# microglia-only arm is re-run here on RNA. [[nhd-sct-fabricates-detection]]
#
# Everything else is copied from 40_mast_percell_dual_FH.R so the two arms differ
# in the cell set and nothing else:
#   * same atlas (atlas/NHD_FH_harmony.rds) and same atlas-wide PrepSCTFindMarkers
#     (the SCT rescaling is scope-dependent, so preparing a subset instead would
#     change the corrected counts and break comparability);
#   * same covariates, derived in-script the same way (matter | hippo_area, log10nCount);
#   * same MAST call: assay SCT, min.pct 0.1, logfc.threshold 0.1, recorrect_umi FALSE,
#     ident.1 = NHD, features = genes_keep (MT-/RPL/RPS/MRPL/MRPS dropped);
#   * same p-value adjustment (Seurat's Bonferroni over tested features);
#   * same RNA-based pct / Cliff's delta / discovery gate.
#
# Cell SET: the 2,824 microglia that Fig 2d/2h use, taken from the shipped assignment
# tables/micro_states/pvm_vs_micro_assignment_FH.csv rather than re-clustered, so this
# script cannot silently disagree with the figures about which nucleus is a microglion.
#
# Output: tables/micro_states/MAST_microglia_only_FH.csv (all tested genes, both regions)
#         tables/micro_states/micro_selection_steps_FH.csv (the selection chain)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(MAST); library(future)
})
set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)
plan("sequential")
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
TBL <- file.path(PROJ, "tables", "micro_states")
dir.create(TBL, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

ASSIGN <- file.path(TBL, "pvm_vs_micro_assignment_FH.csv")
stopifnot("MISSING myeloid assignment — run 2L_pvm_vs_microglia_FH.R" = file.exists(ASSIGN))
asg <- read.csv(ASSIGN, row.names = 1, check.names = FALSE)
MICRO_CELLS <- rownames(asg)[asg$myeloid_class == "Microglia"]
PVM_CELLS   <- rownames(asg)[asg$myeloid_class == "PVM"]
msg("assignment: %d microglia, %d CD163+/F13A1+ (total myeloid %d)",
    length(MICRO_CELLS), length(PVM_CELLS), nrow(asg))
stopifnot("expected 2,824 microglia in the shipped assignment" = length(MICRO_CELLS) == 2824L)

msg("=== loading F+Hippo atlas %s ===", format(Sys.time(), "%H:%M:%S"))
obj <- readRDS(file.path(PROJ, "atlas", "NHD_FH_harmony.rds"))

samp <- as.character(obj$Sample)
obj$matter <- factor(ifelse(grepl("_g$", samp), "grey",
                     ifelse(grepl("_w$", samp), "white", NA_character_)),
                     levels = c("grey","white"))
obj$hippo_area <- factor(ifelse(grepl("Hippo1$", samp), "area1",
                         ifelse(grepl("Hippo2$", samp), "area2", NA_character_)),
                         levels = c("area1","area2"))
obj$Condition   <- factor(as.character(obj$Condition), levels = c("CON","NHD"))
obj$log10nCount <- log10(obj$nCount_RNA)

if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"
msg("PrepSCTFindMarkers (atlas-wide, as in 40_mast) ... %s", format(Sys.time(), "%H:%M:%S"))
obj <- PrepSCTFindMarkers(obj, verbose = FALSE); invisible(gc())

genes_keep <- grep("^(MT-|RPL|RPS|MRPL|MRPS)", rownames(obj), value = TRUE, invert = TRUE)
msg("features kept (MT/ribo dropped): %d / %d", length(genes_keep), nrow(obj))

stopifnot("assignment barcodes are not all in the atlas" =
            all(MICRO_CELLS %in% colnames(obj)))

# ---- the selection chain, counted on the atlas itself ----------------------
md <- obj@meta.data
# NB dplyr:: is explicit throughout this block: MAST loads plyr, whose count()/mutate()
# mask dplyr's and fail on a data.frame with "Argument 'x' is not a vector: list".
sel <- md %>%
  dplyr::mutate(cell = rownames(md)) %>%
  dplyr::filter(new_annotation %in% c("Micro-PVM", "Lymphocyte")) %>%
  dplyr::mutate(step = dplyr::case_when(
    new_annotation == "Lymphocyte"        ~ "T_lymphocyte",
    cell %in% MICRO_CELLS                 ~ "microglia_final",
    cell %in% PVM_CELLS                   ~ "CD163_F13A1_removed",
    TRUE                                  ~ "oligo_doublet_removed")) %>%
  dplyr::count(Condition, Region, step, name = "n")
write.csv(sel, file.path(TBL, "micro_selection_steps_FH.csv"), row.names = FALSE)
msg("\n=== selection chain ===")
print(as.data.frame(sel %>% tidyr::pivot_wider(names_from = step, values_from = n, values_fill = 0L)))

# ---- microglia-only MAST, one run per region -------------------------------
run_one <- function(rg) {
  cells <- intersect(MICRO_CELLS, colnames(obj)[obj$Region == rg])
  sub <- subset(obj, cells = cells)
  tabc <- table(sub$Condition)
  msg("\n=== microglia-only %s: %d nuclei (NHD %d / CON %d) %s ===",
      rg, ncol(sub), tabc["NHD"], tabc["CON"], format(Sys.time(), "%H:%M:%S"))
  if (any(tabc[c("CON","NHD")] < 10)) { msg("  SKIP (per-condition nuclei too few)"); return(NULL) }
  Idents(sub) <- "Condition"
  block <- if (rg == "Hippo") "hippo_area" else "matter"
  has_block <- length(unique(na.omit(sub@meta.data[[block]]))) >= 2
  lv <- if (has_block) c(block, "log10nCount") else "log10nCount"
  msg("  latent.vars = %s", paste(lv, collapse = "+"))
  m <- FindMarkers(sub, ident.1 = "NHD", ident.2 = "CON", test.use = "MAST",
                   latent.vars = lv, assay = "SCT", features = genes_keep,
                   min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi = FALSE,
                   verbose = FALSE)
  m$gene <- rownames(m)
  # Detection and effect size on raw RNA, identical to 40_mast
  rna  <- GetAssayData(sub, assay = "RNA", layer = "counts")[m$gene, , drop = FALSE]
  isN  <- sub$Condition == "NHD"
  m$pct.1 <- round(Matrix::rowMeans(rna[,  isN, drop = FALSE] > 0), 3)
  m$pct.2 <- round(Matrix::rowMeans(rna[, !isN, drop = FALSE] > 0), 3)
  rnad <- GetAssayData(sub, assay = "RNA", layer = "data")[m$gene, , drop = FALSE]
  n1 <- sum(isN); n2 <- sum(!isN)
  m$cliffs_delta <- apply(as.matrix(rnad), 1, function(x) {
    r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
  m$max_pct_raw <- pmax(m$pct.1, m$pct.2)
  out <- m %>% transmute(
    comparison = paste0("Microglia_only_", rg), cell_type = "Microglia_only", region = rg,
    gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, cliffs_delta, max_pct_raw,
    n_NHD = as.integer(tabc["NHD"]), n_CON = as.integer(tabc["CON"]),
    latent_vars = paste(lv, collapse = "+"),
    discovery = !is.na(p_val_adj) & p_val_adj < 0.05 &
                !is.na(cliffs_delta) & abs(cliffs_delta) >= 0.15 &
                max_pct_raw >= 0.10)
  msg("  tested=%d  sig(padj<.05)=%d  discovery=%d",
      nrow(out), sum(out$p_val_adj < 0.05, na.rm = TRUE), sum(out$discovery))
  invisible(gc())
  out
}

# ---- ST17's eight genes: delta/pct for every one, tested or not ------------
# MAST returns only features passing min.pct / logfc.threshold.  A gene the pooled
# run tested and the microglia-only run did not must read as "not tested", never as
# a missing value that a referee could take for "no effect" — the same distinction
# Fig 2e draws between a hollow circle and a cross.
ST17_GENES <- c("HLA-DRA","CD74","FTH1","C1QB","SPP1","CLEC7A","ITGAX","GPNMB")
untested_stats <- function(rg) {
  cells <- intersect(MICRO_CELLS, colnames(obj)[obj$Region == rg])
  sub <- subset(obj, cells = cells)
  g <- intersect(ST17_GENES, rownames(GetAssayData(sub, assay = "RNA", layer = "counts")))
  rna  <- GetAssayData(sub, assay = "RNA", layer = "counts")[g, , drop = FALSE]
  rnad <- GetAssayData(sub, assay = "RNA", layer = "data")[g, , drop = FALSE]
  isN <- sub$Condition == "NHD"; n1 <- sum(isN); n2 <- sum(!isN)
  data.frame(gene = g, region = rg,
             pct.1 = round(Matrix::rowMeans(rna[,  isN, drop = FALSE] > 0), 3),
             pct.2 = round(Matrix::rowMeans(rna[, !isN, drop = FALSE] > 0), 3),
             cliffs_delta = apply(as.matrix(rnad), 1, function(x) {
               r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 }),
             n_NHD = n1, n_CON = n2, row.names = NULL)
}

res <- lapply(c("Frontal","Hippo"), run_one)
allres <- bind_rows(res)
write.csv(allres, file.path(TBL, "MAST_microglia_only_FH.csv"), row.names = FALSE)
msg("\nwrote MAST_microglia_only_FH.csv (%d rows)", nrow(allres))

panel <- bind_rows(lapply(c("Frontal","Hippo"), untested_stats))
write.csv(panel, file.path(TBL, "micro_only_ST17gene_stats_FH.csv"), row.names = FALSE)
msg("wrote micro_only_ST17gene_stats_FH.csv (%d rows; every ST17 gene, tested or not)", nrow(panel))

writeLines(capture.output(sessionInfo()),
           file.path(PROJ, "logs", "92_micro_only_MAST_FH_sessionInfo.txt"))
msg("\n=== DONE %s ===", format(Sys.time(), "%H:%M:%S"))
