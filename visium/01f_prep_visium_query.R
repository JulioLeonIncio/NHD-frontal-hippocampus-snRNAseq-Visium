# =============================================================================
# 01f_prep_visium_query.R — 01f — Visium query export for the spatial mapping step: writes only the Visium query bundle
# (visium_{matrix,genes,barcodes,meta}) from the integrated object. The snRNA reference bundle (sc_ref_*) and its
# signatures (ref_signatures.csv) come from the reference-regression step and are reused by 02f.
# Also writes spot_coords.csv (x = pxl_row, y = pxl_col, both x tissue_lowres_scalef — the
# convention every Fig-6 script relies on).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(jsonlite) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
NHD <- dirname(Sys.getenv("NHD_PROJ"))
VIS <- file.path(NHD, "Visium/integrated_harmony/NHD_frontal_integrated_harmony.rds")
OUT <- file.path(NHD, "Visium/cell2location/c2l_MAIN")
source(file.path(NHD, "NHD_frontal_hippo_rebuild", "scripts", "_visium_outs_FH.R"))   # VIS_SECTIONS + vis_outs(sec) from VISIUM_OUTS_MANIFEST.json
write_mtx_bundle <- function(counts, meta, prefix) {
  Matrix::writeMM(as(counts, "CsparseMatrix"), paste0(prefix, "_matrix.mtx"))
  writeLines(rownames(counts), paste0(prefix, "_genes.tsv")); writeLines(colnames(counts), paste0(prefix, "_barcodes.tsv"))
  write.csv(meta, paste0(prefix, "_meta.csv"), row.names = TRUE)
  cat("  wrote bundle:", basename(prefix), "(", nrow(counts), "genes x", ncol(counts), "cells )\n") }
vis <- readRDS(VIS); DefaultAssay(vis) <- "Spatial"; vis <- JoinLayers(vis, assay = "Spatial")
vis_counts <- GetAssayData(vis, assay = "Spatial", layer = "counts")
vis_meta_cols <- intersect(c("sample_id","condition","seurat_clusters","Spatial_snn_res.0.5","percent.mt","nCount_Spatial","nFeature_Spatial"), colnames(vis@meta.data))
vis_meta <- vis@meta.data[, vis_meta_cols, drop = FALSE]
if (length(vis@images) > 0) {
  co_all <- do.call(rbind, lapply(vis@images, function(im) { co <- tryCatch(im@coordinates, error = function(e) NULL); if (is.null(co)) return(NULL)
    co[, intersect(c("imagerow","imagecol","row","col"), colnames(co)), drop = FALSE] }))
  if (!is.null(co_all)) { common <- intersect(rownames(vis_meta), rownames(co_all))
    for (cn in colnames(co_all)) vis_meta[[cn]] <- NA_real_
    vis_meta[common, colnames(co_all)] <- co_all[common, , drop = FALSE] } }
cat("Visium spots:", ncol(vis_counts), "| condition:\n"); print(table(vis_meta$condition, useNA = "ifany"))
cat("seurat_clusters (domains):\n"); print(table(vis_meta$seurat_clusters, useNA = "ifany"))
write_mtx_bundle(vis_counts, vis_meta, file.path(OUT, "visium"))
## spot_coords.csv from the Space Ranger positions of the same outs the object was built from
co <- do.call(rbind, lapply(VIS_SECTIONS, function(s) {   # outs per section via the manifest
  sp <- file.path(vis_outs(s), "spatial"); sf <- fromJSON(file.path(sp, "scalefactors_json.json"))$tissue_lowres_scalef
  tp <- read.csv(file.path(sp, "tissue_positions_list.csv"), header = FALSE, col.names = c("bc","it","ar","ac","pr","pc"))
  data.frame(cell = paste0(s, "_", tp$bc), x = tp$pr * sf, y = tp$pc * sf) }))
co <- co[co$cell %in% colnames(vis), ]; stopifnot(nrow(co) == ncol(vis))
write.csv(co, file.path(OUT, "spot_coords.csv"), row.names = FALSE)
cat("spot_coords.csv:", nrow(co), "spots\n=== 01f DONE ===\n")
