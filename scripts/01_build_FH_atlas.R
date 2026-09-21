#!/usr/bin/env Rscript
# 01_build_FH_atlas.R — Foundation. Re-cluster/normalize the Frontal+Hippocampus data with OCC removed and QC tightened to nFeature_RNA > 200. Mirrors Steps 5-9 of the
# original NHD_QC_Azimuth+Harmony.R exactly (SCTransform per sample -> integration
# features (MT/ribo removed) -> merge -> PCA -> Harmony(by SampleID) -> UMAP ->
# FindClusters 0.6/0.8/1). QC + doublet removal were per-lane and are already baked
# into the source object, so we start from its preserved RNA counts.
#
# Inputs : ../../NHD_QC_harmony.rds  (relaxed-QC atlas, floor nFeature 121)
# Outputs: atlas/NHD_FH_subset_RNA.rds  (checkpoint, RNA-only subset, pre-SCT)
#          atlas/NHD_FH_harmony.rds     (re-clustered F+Hippo atlas)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(harmony); library(dplyr)
})
set.seed(42)
options(future.globals.maxSize = 8000*1024^2)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
OUT  <- file.path(PROJ, "NHD_frontal_hippo_rebuild", "atlas")
KEEP_REGIONS <- c("Frontal","Hippo")

cat("==== LOAD ====\n"); t0 <- Sys.time()
obj <- readRDS(file.path(PROJ, "NHD_QC_harmony.rds"))
cat("loaded", ncol(obj), "cells in", round(difftime(Sys.time(),t0,units="mins"),1), "min\n")

# ---- Subset: Frontal+Hippo, nFeature>200, nFeature<6500, mito<10 ----
md <- obj[[]]
keep <- md$Region %in% KEEP_REGIONS &
        md$nFeature_RNA > 200 & md$nFeature_RNA < 6500 &
        md$percent_mito < 10
cat("\n==== SUBSET ====\n")
cat("cells kept:", sum(keep), "of", length(keep), "\n")
cat("dropped by region (OCC):", sum(!(md$Region %in% KEEP_REGIONS)), "\n")
cat("dropped by nFeature<=200 (within F+Hippo):",
    sum(md$Region %in% KEEP_REGIONS & md$nFeature_RNA <= 200), "\n")
sub <- subset(obj, cells = colnames(obj)[keep])
rm(obj); gc()

# strip old SCT / prediction-score assays + reductions; keep only RNA counts
DefaultAssay(sub) <- "RNA"
sub <- DietSeurat(sub, assays = "RNA", dimreducs = NULL, graphs = NULL)
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])

cat("\n-- subset verification --\n")
cat("regions:\n"); print(table(sub$Region))
cat("lanes x region:\n"); print(table(sub$SampleID, sub$Region))
cat("cell types (new_annotation):\n"); print(table(sub$new_annotation))
cat("assays:", paste(Assays(sub), collapse=", "), "| cells:", ncol(sub), "\n")
saveRDS(sub, file.path(OUT, "NHD_FH_subset_RNA.rds"))
cat("checkpoint saved: NHD_FH_subset_RNA.rds\n")

# preserve annotation map (re-attach after split/merge, as in the original)
annotation_map <- data.frame(
  cell               = colnames(sub),
  new_annotation     = as.character(sub$new_annotation),
  predicted.subclass = as.character(sub$predicted.subclass),
  stringsAsFactors = FALSE)

# ---- Step 5: SCTransform per sample ----
cat("\n==== SCTransform per SampleID ====\n"); t1 <- Sys.time()
lst <- SplitObject(sub, split.by = "SampleID")
lst <- lapply(lst, SCTransform, assay = "RNA", vst.flavor = "v2", verbose = FALSE)
cat("SCTransform done in", round(difftime(Sys.time(),t1,units="mins"),1), "min\n")

# ---- Step 6: integration features (drop MT / ribo) ----
rm_pat <- "^(MT-|RPL|RPS|MRPL|MRPS)"
all.features   <- SelectIntegrationFeatures(object.list = lst, nfeatures = 5000)
brain.features <- all.features[!grepl(rm_pat, all.features, ignore.case = TRUE)]
brain.features <- head(brain.features, 3000)
cat("integration features:", length(brain.features), "(from 5000, ribo/MT removed)\n")

# ---- Step 7: merge SCT list; re-attach annotations ----
harm <- merge(x = lst[[1]], y = lst[2:length(lst)], merge.data = TRUE)
harm$new_annotation     <- annotation_map$new_annotation[match(colnames(harm), annotation_map$cell)]
harm$predicted.subclass <- annotation_map$predicted.subclass[match(colnames(harm), annotation_map$cell)]
cat("NAs new_annotation:", sum(is.na(harm$new_annotation)), "\n")

# ---- Step 8: PCA -> Harmony -> UMAP -> clusters (exact original sequence) ----
VariableFeatures(harm) <- brain.features
DefaultAssay(harm) <- "SCT"
cat("\n==== PCA / Harmony / UMAP / clusters ====\n"); t2 <- Sys.time()
harm <- harm %>%
  RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.8, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)
harm <- RunHarmony(harm, group.by.vars = "SampleID")
harm <- harm %>%
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = c(0.6, 0.8, 1), verbose = FALSE)
cat("integration done in", round(difftime(Sys.time(),t2,units="mins"),1), "min\n")

# ---- Step 9: save + verify ----
harm$Region <- factor(as.character(harm$Region), levels = KEEP_REGIONS)
saveRDS(harm, file.path(OUT, "NHD_FH_harmony.rds"))
cat("\n==== SAVED NHD_FH_harmony.rds ====\n")
cat("cells:", ncol(harm), " | assays:", paste(Assays(harm), collapse=", "),
    " | reductions:", paste(Reductions(harm), collapse=", "), "\n")
cat("Region:\n"); print(table(harm$Region))
cat("cell types:\n"); print(table(harm$new_annotation))
cat("clusters (harmony res 0.8 -> seurat_clusters):\n"); print(table(harm$seurat_clusters))
cat("celltype x cluster crosstab (sanity):\n")
print(table(harm$new_annotation, harm$seurat_clusters))
cat("\nTOTAL runtime", round(difftime(Sys.time(),t0,units="mins"),1), "min. Done.\n")
