# =============================================================================
# 01e_prep_NHD.R — NHD frontal deconvolution — port of SCZ 01e_prep_MAIN_BANKSY.R (the proven, successful local recipe). Exports snRNA reference (8 broad new_annotation main cell types: Oligo, Neuron_Ex, Neuron_Inh, Astro, OPC, Micro-PVM, Pericytes, Endo) + integrated 4-frontal Visium query to dependency-free mtx bundles.
#   Reference label = raw new_annotation (broad), batch = SampleID.
#   Query = NHD_frontal_integrated_harmony.rds (17,976 spots @ 50/25 QC, 4 imgs,
#   carries sample_id / condition / seurat_clusters domains for per-domain export).
# Output: Visium/cell2location/c2l_MAIN/{sc_ref_*, visium_*}
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
write_mtx_bundle <- function(counts, meta, prefix) {
  Matrix::writeMM(as(counts, "CsparseMatrix"), paste0(prefix, "_matrix.mtx"))
  writeLines(rownames(counts), paste0(prefix, "_genes.tsv"))
  writeLines(colnames(counts), paste0(prefix, "_barcodes.tsv"))
  write.csv(meta, paste0(prefix, "_meta.csv"), row.names = TRUE)
  cat("  wrote bundle:", basename(prefix), "(", nrow(counts), "genes x", ncol(counts), "cells )\n")
}
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
NHD   <- dirname(Sys.getenv("NHD_PROJ"))
ATLAS <- file.path(NHD, "NHD_QC_harmony.rds")
VIS   <- file.path(NHD, "Visium/integrated_harmony/NHD_frontal_integrated_harmony.rds")
OUT   <- file.path(NHD, "Visium/cell2location/c2l_MAIN")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- reference: broad new_annotation (8 main cell types) -------------------
cat("=== Loading atlas (broad new_annotation) ===\n")
obj <- readRDS(ATLAS); DefaultAssay(obj) <- "RNA"; obj <- JoinLayers(obj, assay = "RNA")
stopifnot("SampleID" %in% colnames(obj@meta.data), "new_annotation" %in% colnames(obj@meta.data))
if (anyNA(obj$SampleID)) stop("SampleID has NA — c2l batch correction would be invalid")
obj$c2l_label <- as.character(obj$new_annotation)          # RAW broad labels (no refinement)
obj <- subset(obj, cells = colnames(obj)[!is.na(obj$c2l_label)])
cat("Broad main-cell-type labels:\n"); print(table(obj$c2l_label))

# stratified subsample (cell2location estimates signatures from a subsample)
CELLS_PER_LABEL <- 3000; set.seed(42)          # broad classes -> allow more cells/label
keep_cells <- unlist(lapply(split(colnames(obj), obj$c2l_label), function(cc)
  if (length(cc) <= CELLS_PER_LABEL) cc else sample(cc, CELLS_PER_LABEL)), use.names = FALSE)
cat(sprintf("Reference subsample: %d -> %d cells (cap %d/label)\n", ncol(obj), length(keep_cells), CELLS_PER_LABEL))
cond_col   <- if ("Condition" %in% colnames(obj@meta.data)) "Condition" else NA
ref_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[, keep_cells]
ref_meta   <- data.frame(c2l_label = obj$c2l_label[keep_cells],
                         SampleID  = as.character(obj$SampleID[keep_cells]),
                         Condition = if (!is.na(cond_col)) as.character(obj@meta.data[keep_cells, cond_col]) else NA,
                         row.names = keep_cells)
rm(obj); invisible(gc())
cat("Reference cells per label (post-subsample):\n"); print(table(ref_meta$c2l_label))
write_mtx_bundle(ref_counts, ref_meta, file.path(OUT, "sc_ref")); rm(ref_counts); invisible(gc())

## ---- Visium query: integrated 4 frontal sections --------------------------
cat("\n=== Loading Visium (integrated harmony, 4 frontal) ===\n")
vis <- readRDS(VIS); DefaultAssay(vis) <- "Spatial"; vis <- JoinLayers(vis, assay = "Spatial")
vis_counts <- GetAssayData(vis, assay = "Spatial", layer = "counts")
vis_meta_cols <- intersect(c("sample_id","condition","seurat_clusters","Spatial_snn_res.0.5",
                             "percent.mt","nCount_Spatial","nFeature_Spatial"), colnames(vis@meta.data))
vis_meta <- vis@meta.data[, vis_meta_cols, drop = FALSE]
# carry image coords per spot (for spatial proportion maps in export)
if (length(vis@images) > 0) {
  co_all <- do.call(rbind, lapply(vis@images, function(im) {
    co <- tryCatch(im@coordinates, error = function(e) NULL); if (is.null(co)) return(NULL)
    co[, intersect(c("imagerow","imagecol","row","col"), colnames(co)), drop = FALSE] }))
  if (!is.null(co_all)) { common <- intersect(rownames(vis_meta), rownames(co_all))
    for (cn in colnames(co_all)) vis_meta[[cn]] <- NA_real_
    vis_meta[common, colnames(co_all)] <- co_all[common, , drop = FALSE] }
}
cat("Visium spots:", ncol(vis_counts), "| condition:\n"); print(table(vis_meta$condition, useNA = "ifany"))
cat("seurat_clusters (domains):\n"); print(table(vis_meta$seurat_clusters, useNA = "ifany"))
write_mtx_bundle(vis_counts, vis_meta, file.path(OUT, "visium"))
cat("\n=== 01e NHD prep DONE ===\n", file = stderr())
