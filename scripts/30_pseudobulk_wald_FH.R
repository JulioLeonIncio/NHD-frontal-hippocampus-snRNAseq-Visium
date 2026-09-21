#!/usr/bin/env Rscript
# =============================================================================
# 30_pseudobulk_wald_FH.R — Pseudobulk DESeq2 tier-2, frontal+hippo rebuild.
# Faithful copy of scripts/30_dge_with_wald_stat.R with:
#   * input  = atlas/NHD_FH_harmony.rds
#   * regions = c("Frontal","Hippo")   (OCC removed)
#   * output = NHD_frontal_hippo_rebuild/diagnostics/03_DGE_pseudobulk/
# Emits both unshrunken (Wald `stat` for GSEA) and apeglm-shrunken LFC (volcanoes).
# Frontal = genuine 2v2 grey/white lanes (real confirmation). Hippo NHD = 1 lane
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(DESeq2); library(apeglm); library(tibble)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
DIR  <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
OUT  <- file.path(DIR, "diagnostics", "03_DGE_pseudobulk"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
RM_PATTERN <- "^(MT-|RPL|RPS|MRPL|MRPS)"

cat("\n=== Loading F+Hippo atlas ===\n")
obj <- readRDS(file.path(DIR, "atlas", "NHD_FH_harmony.rds"))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
genes_keep <- setdiff(rownames(obj), grep(RM_PATTERN, rownames(obj), value = TRUE))
cat("Kept", length(genes_keep), "genes (dropped MT/ribo)\n")

cell_types <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
regions    <- c("Frontal","Hippo")     # OCC removed

run_pb <- function(cells, label) {
  if (length(cells) < 30) { cat("  SKIP", label, "- few cells\n"); return(NULL) }
  s <- subset(obj, cells = cells)
  md_s <- s@meta.data
  pseudo_id <- as.character(md_s$SampleID)
  cond_lane_counts <- md_s %>% group_by(SampleID, Condition) %>%
    summarise(n_cells = dplyr::n(), .groups = "drop")
  for (cond in unique(as.character(md_s$Condition))) {
    cl <- cond_lane_counts %>% filter(Condition == cond)
    if (sum(cl$n_cells >= 10) < 2) {
      idx <- which(md_s$Condition == cond)
      if (length(idx) < 20) { cat("  SKIP", label, ": cannot pseudo-split\n"); return(NULL) }
      half <- sample(idx, length(idx) %/% 2)
      pseudo_id[half]               <- paste0(cond, "_pool_a")
      pseudo_id[setdiff(idx, half)] <- paste0(cond, "_pool_b")
    }
  }
  s$pseudoSample <- pseudo_id
  uniq_p <- unique(pseudo_id)
  if (length(uniq_p) < 4) { cat("  SKIP", label, ": <4 lanes\n"); return(NULL) }

  counts_mat <- GetAssayData(s, assay = "RNA", layer = "counts")[genes_keep, ]
  pb <- sapply(uniq_p, function(sid) Matrix::rowSums(counts_mat[, pseudo_id == sid, drop = FALSE]))
  colnames(pb) <- uniq_p

  meta <- data.frame(
    pseudoSample = uniq_p,
    Condition = sapply(uniq_p, function(x) {
      if (grepl("_pool_[ab]$", x)) sub("_pool_[ab]$", "", x)
      else as.character(md_s$Condition[md_s$SampleID == x][1]) }),
    n_cells = sapply(uniq_p, function(x) sum(pseudo_id == x)),
    stringsAsFactors = FALSE)
  rownames(meta) <- meta$pseudoSample
  meta$Condition <- factor(meta$Condition, levels = c("CON","NHD"))
  if (any(table(meta$Condition) < 2)) { cat("  SKIP", label, ": uneven\n"); return(NULL) }

  pb <- pb[rowSums(pb >= 5) >= 2, ]
  if (nrow(pb) < 200) { cat("  SKIP", label, ": few genes\n"); return(NULL) }

  dds <- DESeqDataSetFromMatrix(countData = round(pb), colData = meta, design = ~ Condition)
  dds <- DESeq(dds, quiet = TRUE)
  res_un <- as.data.frame(results(dds, contrast = c("Condition","NHD","CON"))) %>% rownames_to_column("gene")
  res_sh <- as.data.frame(lfcShrink(dds, coef = "Condition_NHD_vs_CON", type = "apeglm", quiet = TRUE)) %>%
    rownames_to_column("gene") %>%
    select(gene, log2FoldChange_apeglm = log2FoldChange, lfcSE_apeglm = lfcSE)

  is_split <- any(grepl("_pool_", uniq_p))
  out <- res_un %>% left_join(res_sh, by = "gene") %>%
    mutate(comparison = label, n_cells_total = ncol(s), n_lanes = length(uniq_p),
           n_NHD_lanes = sum(meta$Condition == "NHD"), n_CON_lanes = sum(meta$Condition == "CON"),
           note = if (is_split) "pseudo-split under-represented condition" else NA_character_,
           dispersion_basis = if (is_split) "pseudo_split" else "biological_2v2")
  cat(sprintf("  %-22s genes=%d  padj<0.01&FC>2.5=%d  basis=%s\n", label, nrow(out),
              sum(out$padj < 0.01 & abs(out$log2FoldChange) > log2(2.5), na.rm = TRUE),
              out$dispersion_basis[1]))
  out
}

cat("\n=== DESeq2 with Wald stat (F+Hippo) ===\n")
all_res <- list()
for (rg in regions) for (ct in cell_types) {
  label <- paste0(ct, "_", rg)
  cells <- WhichCells(obj, expression = new_annotation == ct & Region == rg)
  all_res[[label]] <- run_pb(cells, label)
}
all_df <- bind_rows(all_res)

# Retire pseudo-split (Hippo): null inference so the merge treats it as
# discovery+cross-region only, never as within-region pseudobulk confirmation.
retired_lab <- unique(all_df$comparison[all_df$dispersion_basis == "pseudo_split"])
cat("\nRetiring pseudo-split comparisons (padj/pvalue/stat nulled):\n  ",
    paste(retired_lab, collapse = ", "), "\n")
all_df <- all_df %>% mutate(
  padj  = ifelse(dispersion_basis == "pseudo_split", NA_real_, padj),
  pvalue= ifelse(dispersion_basis == "pseudo_split", NA_real_, pvalue),
  stat  = ifelse(dispersion_basis == "pseudo_split", NA_real_, stat))

write.csv(all_df, file.path(OUT, "pseudobulk_DESeq2_with_Wald_stat.csv"), row.names = FALSE)
cat("\nSaved", nrow(all_df), "rows to pseudobulk_DESeq2_with_Wald_stat.csv\n")
cat("\n=== genuine 2v2 (Frontal) confirmation summary ===\n")
print(as.data.frame(all_df %>% filter(dispersion_basis == "biological_2v2") %>% group_by(comparison) %>%
  summarise(n=dplyr::n(), strong = sum(padj<0.01 & abs(log2FoldChange)>log2(2.5), na.rm=TRUE), .groups="drop")))
cat("\n=== DONE ===\n", file = stderr()); print(sessionInfo())
