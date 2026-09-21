# =============================================================================
# 04_zhou_cross_cohort_validation.R — Script 04 — Zhou et al. 2023 Nat Immunol cross-cohort validation
# Design: 3 NHD (DAP12 c.2T>C or c.141Gdel) vs 11 CON, occipital cortex
# Approach:
#   1. Pseudobulk DESeq2 per cell type (genuine biological replicates here)
#   2. GSEA with Wald stat (fgsea, Reactome + GO:BP + DAM/curated signatures)
#   3. Signature scoring: Zhou wound-repair + DAM + our NHD curated gene sets
#   4. Overlap of Zhou significant DEGs with our Frontal top DEGs
# Output: data/zhou_validation.rds + tables/zhou_*.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tibble)
  library(DESeq2); library(apeglm); library(fgsea); library(msigdbr)
  library(ggplot2); library(tidyr)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
DATA <- file.path(PROJ, "manuscript_7fig", "data")
TBL  <- file.path(PROJ, "manuscript_7fig", "tables")
RM_PATTERN <- "^(MT-|RPL|RPS|MRPL|MRPS)"

# Sourced for the shared DE gate constants only. GSEA here ranks on the DESeq2
# Wald `stat` (genuine biological replicates in the Zhou cohort -> a legitimate
# non-FC path). The DE gate now feeds only the display-only dge_member flag
# , never a pathway pre-filter.
source(file.path(Sys.getenv("NHD_PROJ"), "scripts", "_gsea_composite_metric.R"))   # same helper as the main analysis
DE_PADJ <- 0.01
DE_LFC  <- log2(2.5)   # Zhou's own sig-gene gate (used below for overlap too)

# ---- Load objects ----
cat("=== Loading Zhou cohort RDS ===\n")
zhou <- readRDS(file.path(PROJ, "NHD_repo", "NHD_complete_harmony+PMI.rds"))
DefaultAssay(zhou) <- "RNA"
zhou <- JoinLayers(zhou, assay = "RNA")  # Seurat v5: merge split layers

# Harmonize condition label: "control" -> "CON", "NHD" -> "NHD"
zhou$Condition_std <- ifelse(zhou$Condition == "NHD", "NHD", "CON")
zhou$Condition_std <- factor(zhou$Condition_std, levels = c("CON","NHD"))

cat("Zhou total cells:", ncol(zhou), "\n")
cat("Condition breakdown:\n"); print(table(zhou$Condition_std))
cat("Cell types:\n"); print(table(zhou$new_annotation))

# Load curated signatures
sigs <- readRDS(file.path(DATA, "curated_signatures.rds"))
# Load our NHD Frontal top DEGs for overlap test
nhd_dge <- read.csv(file.path(PROJ, "diagnostics", "03_DGE_pseudobulk",
                               "pseudobulk_DESeq2_with_Wald_stat.csv"))

genes_keep <- setdiff(rownames(zhou), grep(RM_PATTERN, rownames(zhou), value=TRUE))

# ---- 1. Pseudobulk DESeq2 per cell type in Zhou cohort ----
cat("\n=== Pseudobulk DESeq2 per cell type (Zhou cohort) ===\n")

run_pb_zhou <- function(cell_type, min_cells = 50) {
  cells <- WhichCells(zhou, expression = new_annotation == cell_type)
  if (length(cells) < min_cells) {
    cat(sprintf("  %-15s SKIP (only %d cells)\n", cell_type, length(cells))); return(NULL)
  }
  s <- subset(zhou, cells = cells)
  md <- s@meta.data

  # Sample-level pseudobulk (genuine N=14 here)
  sids <- unique(as.character(md$SampleID))
  sids <- sids[sapply(sids, function(x) sum(md$SampleID == x) >= 10)]
  if (length(sids) < 4) return(NULL)

  counts_mat <- GetAssayData(s, assay = "RNA", layer = "counts")[genes_keep, ]
  pb <- sapply(sids, function(sid)
                 Matrix::rowSums(counts_mat[, md$SampleID == sid, drop=FALSE]))
  pb <- pb[rowSums(pb >= 5) >= 2, ]
  if (nrow(pb) < 200) return(NULL)

  meta <- data.frame(
    pseudoSample = sids,
    Condition    = sapply(sids, function(x) as.character(md$Condition_std[md$SampleID == x][1])),
    n_cells      = sapply(sids, function(x) sum(md$SampleID == x)),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$pseudoSample
  meta$Condition <- factor(meta$Condition, levels = c("CON","NHD"))
  if (any(table(meta$Condition) < 2)) return(NULL)

  pb <- pb[, sids]
  dds <- tryCatch({
    d <- DESeqDataSetFromMatrix(countData = round(pb), colData = meta, design = ~ Condition)
    DESeq(d, quiet = TRUE)
  }, error = function(e) { cat("  ERROR:", e$message, "\n"); NULL })
  if (is.null(dds)) return(NULL)

  res_un <- as.data.frame(results(dds, contrast = c("Condition","NHD","CON"))) %>%
    rownames_to_column("gene")
  res_sh <- tryCatch(
    as.data.frame(lfcShrink(dds, coef = "Condition_NHD_vs_CON",
                             type = "apeglm", quiet = TRUE)) %>%
      rownames_to_column("gene") %>%
      select(gene, log2FoldChange_apeglm = log2FoldChange, lfcSE_apeglm = lfcSE),
    error = function(e) data.frame(gene = res_un$gene,
                                   log2FoldChange_apeglm = NA, lfcSE_apeglm = NA))
  out <- res_un %>%
    left_join(res_sh, by = "gene") %>%
    mutate(cell_type = cell_type,
           n_cells_total = ncol(s),
           n_NHD = sum(meta$Condition == "NHD"),
           n_CON = sum(meta$Condition == "CON"))
  cat(sprintf("  %-15s genes=%d  padj<0.01_FC>2.5=%d  (N=%d NHD vs %d CON)\n",
              cell_type, nrow(out),
              sum(out$padj < 0.01 & abs(out$log2FoldChange) > log2(2.5), na.rm=TRUE),
              sum(meta$Condition=="NHD"), sum(meta$Condition=="CON")))
  out
}

cell_types_zhou <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
zhou_dge_list <- lapply(cell_types_zhou, run_pb_zhou)
names(zhou_dge_list) <- cell_types_zhou
zhou_dge <- bind_rows(zhou_dge_list)
write.csv(zhou_dge, file.path(TBL, "zhou_pseudobulk_DESeq2.csv"), row.names=FALSE)
cat("Saved zhou_pseudobulk_DESeq2.csv:", nrow(zhou_dge), "rows\n")

# ---- 2. GSEA per cell type ----
cat("\n=== GSEA per cell type (Zhou cohort) ===\n")

# Gene sets: Reactome + GO:BP + Hallmark
gobp     <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP")
reactome <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME")
hallmark <- msigdbr(species="Homo sapiens", collection="H")
gs_all   <- bind_rows(gobp, reactome, hallmark)
gs_list  <- split(gs_all$gene_symbol, gs_all$gs_name)

run_gsea_zhou <- function(dge_ct, ct_label) {
  ranks <- dge_ct %>%
    filter(!is.na(stat)) %>%
    distinct(gene, .keep_all=TRUE) %>%
    { setNames(.$stat, .$gene) }
  ranks <- sort(ranks, decreasing=TRUE)
  # ---- FDR denominator fix -----
  # fgsea runs over the full minSize/maxSize-passing family (gs_list) so BH `padj`
  # uses the valid multiple-testing denominator. The old filter_sets_by_dge()
  # pre-filter (dropping every zero-DE set before fgsea) shrank the BH family and
  # was anti-conservative. DGE-membership is now a display-only flag (`dge_member`)
  # recorded per pathway; it never changes the correction denominator. Mirrors 44.
  dge_genes <- unique(dge_ct$gene[which(dge_ct$padj < DE_PADJ &
                                        abs(dge_ct$log2FoldChange) > DE_LFC)])
  set.seed(42)
  res <- tryCatch(
    fgsea(pathways = gs_list, stats = ranks, minSize=10, maxSize=500, nPermSimple=1000),
    error = function(e) NULL)
  if (is.null(res) || !nrow(res)) return(NULL)
  # display-only DGE-membership: >=1 Zhou-cohort DE gene in the tested set?
  res$dge_member <- if (length(dge_genes))
    vapply(gs_list[res$pathway], function(g) any(g %in% dge_genes), logical(1)) else
    rep(FALSE, nrow(res))
  res$cell_type <- ct_label
  res$FDR_tier  <- ifelse(res$padj < 0.05, "FDR<0.05", ifelse(res$padj < 0.25, "FDR<0.25", "ns"))
  res
}

zhou_gsea_list <- lapply(cell_types_zhou, function(ct) {
  sub <- zhou_dge_list[[ct]]
  if (is.null(sub) || !"stat" %in% colnames(sub)) return(NULL)
  cat("  GSEA:", ct, "\n")
  run_gsea_zhou(sub, ct)
})
zhou_gsea <- bind_rows(zhou_gsea_list) %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=";")))
write.csv(zhou_gsea, file.path(TBL, "zhou_GSEA.csv"), row.names=FALSE)
cat("Saved zhou_GSEA.csv:", nrow(zhou_gsea), "rows\n")

# ---- 3. Signature scoring per cell per Zhou cohort ----
cat("\n=== Signature scoring (Zhou cohort) ===\n")

score_signatures <- function(seurat_obj, sig_list, assay="RNA") {
  scores_list <- lapply(names(sig_list), function(nm) {
    gs <- intersect(sig_list[[nm]], rownames(seurat_obj))
    if (length(gs) < 3) return(NULL)
    seurat_obj <- AddModuleScore(seurat_obj, features=list(gs), name="sig",
                                 seed=42, assay=assay)
    data.frame(
      cell     = colnames(seurat_obj),
      score    = seurat_obj$sig1,
      signature = nm,
      Condition = seurat_obj$Condition_std,
      cell_type = seurat_obj$new_annotation,
      SampleID  = seurat_obj$SampleID,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(scores_list)
}

# Focus on microglia signatures for headline validation
zhou_micro <- subset(zhou, new_annotation == "Micro-PVM")
# Seurat v5: data layer may be empty after subset — normalize counts
zhou_micro <- NormalizeData(zhou_micro, assay="RNA", verbose=FALSE)
micro_sigs_to_score <- sigs[c("DAM_up","Homeostatic_microglia","Zhou_NHD_microglia_up","LDAM","IRM")]
micro_scores <- score_signatures(zhou_micro, micro_sigs_to_score)
write.csv(micro_scores, file.path(TBL, "zhou_micro_signature_scores.csv"), row.names=FALSE)

# Also score astrocytes
zhou_astro <- subset(zhou, new_annotation == "Astro")
zhou_astro <- NormalizeData(zhou_astro, assay="RNA", verbose=FALSE)
astro_sigs_to_score <- sigs[c("Zhou_NHD_astrocyte_up","A1_reactive","A2_reactive","DAA_Habib")]
astro_scores <- score_signatures(zhou_astro, astro_sigs_to_score)
write.csv(astro_scores, file.path(TBL, "zhou_astro_signature_scores.csv"), row.names=FALSE)
cat("Signature scores saved.\n")

# ---- 4. DEG overlap: Zhou Micro-PVM vs our NHD Frontal Micro-PVM ----
cat("\n=== DEG overlap test (Zhou vs our Micro-PVM Frontal) ===\n")

our_micro_frontal <- nhd_dge %>%
  filter(comparison == "Micro-PVM_Frontal",
         padj < 0.01, abs(log2FoldChange) > log2(2.5)) %>%
  mutate(direction = ifelse(log2FoldChange > 0, "NHD_up", "CON_up"))

zhou_micro_dge <- zhou_dge_list[["Micro-PVM"]]
if (!is.null(zhou_micro_dge)) {
  zhou_micro_sig <- zhou_micro_dge %>%
    filter(padj < 0.01, abs(log2FoldChange) > log2(2.5)) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "NHD_up", "CON_up"))

  overlap_up <- intersect(
    our_micro_frontal$gene[our_micro_frontal$direction == "NHD_up"],
    zhou_micro_sig$gene[zhou_micro_sig$direction == "NHD_up"])
  overlap_dn <- intersect(
    our_micro_frontal$gene[our_micro_frontal$direction == "CON_up"],
    zhou_micro_sig$gene[zhou_micro_sig$direction == "CON_up"])

  cat(sprintf("Our Micro-PVM_Frontal sig genes: %d up + %d down in NHD\n",
              sum(our_micro_frontal$direction=="NHD_up"),
              sum(our_micro_frontal$direction=="CON_up")))
  cat(sprintf("Zhou Micro-PVM sig genes: %d up + %d down in NHD\n",
              sum(zhou_micro_sig$direction=="NHD_up"),
              sum(zhou_micro_sig$direction=="CON_up")))
  cat(sprintf("Overlap NHD-up: %d genes\n", length(overlap_up)))
  cat(sprintf("Overlap CON-up: %d genes\n", length(overlap_dn)))
  if (length(overlap_up) > 0) cat("Shared NHD-up:", paste(head(overlap_up,20), collapse=", "), "\n")

  overlap_rows <- list()
  if (length(overlap_up) > 0) overlap_rows[["up"]]  <- data.frame(gene=overlap_up, direction="NHD_up")
  if (length(overlap_dn) > 0) overlap_rows[["down"]] <- data.frame(gene=overlap_dn, direction="CON_up")
  if (length(overlap_rows) > 0) {
    overlap_tbl <- bind_rows(overlap_rows)
    write.csv(overlap_tbl, file.path(TBL, "zhou_micro_pvm_overlap.csv"), row.names=FALSE)
  }
}

# ---- 5. Save consolidated validation object ----
zhou_validation <- list(
  dge            = zhou_dge,
  gsea           = zhou_gsea,
  micro_scores   = micro_scores,
  astro_scores   = astro_scores,
  meta           = data.frame(
    n_NHD         = 3,
    n_CON         = 11,
    tissue        = "occipital cortex",
    genotype_NHD  = "DAP12 c.2T>C / c.141Gdel",
    source        = "Zhou et al. 2023 Nat Immunol",
    stringsAsFactors = FALSE)
)
saveRDS(zhou_validation, file.path(DATA, "zhou_validation.rds"))
cat("Saved zhou_validation.rds\n")

cat("\n=== DONE ===\n")
