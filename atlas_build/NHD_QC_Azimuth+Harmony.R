# NHD_QC_Azimuth+Harmony.R — per-library quality control, doublet removal, Azimuth annotation and the first Harmony integration of the single-nucleus libraries.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
library(Seurat)
library(Azimuth)
library(scDblFinder)
library(scCustomize)
library(dplyr)
library(ggplot2)
library(future)

out_dir <- Sys.getenv("NHD_ATLAS_OUT", "~/Analysis/NHD/QC")   # where the first-pass atlas is written
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────
# STEP 1: Load H5 files, mito %, QC filter
# ─────────────────────────────────────────────
h5_files <- list.files(Sys.getenv("NHD_CELLBENDER", "~/Analysis/NHD/transfers/CellBender_results"),   # CellBender-filtered h5 per library
                       pattern = "_cellbender_filtered\\.h5$",
                       full.names = TRUE,
                       recursive = TRUE)

seurat_list <- lapply(h5_files, function(file) {
  sample_name <- basename(dirname(file))
  
  # Read CellBender h5 manually
  f        <- hdf5r::H5File$new(file, mode = "r")
  barcodes <- f[["matrix/barcodes"]][]
  genes    <- f[["matrix/features/name"]][]
  data     <- f[["matrix/data"]][]
  indices  <- f[["matrix/indices"]][]
  indptr   <- f[["matrix/indptr"]][]
  shape    <- f[["matrix/shape"]][]
  f$close_all()
  
  mat <- Matrix::sparseMatrix(
    i    = indices + 1,
    p    = indptr,
    x    = as.numeric(data),
    dims = shape,
    dimnames = list(genes, barcodes)
  )
  
  # Make gene names unique
  rownames(mat) <- make.unique(rownames(mat))
  
  obj <- CreateSeuratObject(counts = mat, project = sample_name)
  obj$SampleID <- sample_name
  
  # Mitochondrial percentage
  mito.genes   <- grep("^MT-", rownames(obj), value = TRUE)
  mito.counts  <- colSums(GetAssayData(obj, layer = "counts")[mito.genes, , drop = FALSE])
  total.counts <- colSums(GetAssayData(obj, layer = "counts"))
  obj$percent_mito <- ifelse(total.counts > 0, (mito.counts / total.counts) * 100, 0)
  
  # QC filtering
  obj <- subset(obj, subset = nFeature_RNA > 350 &
                  nFeature_RNA < 6500 &
                  percent_mito < 10)
  
  Idents(obj) <- sample_name
  obj
})

names(seurat_list) <- sapply(h5_files, function(f) basename(dirname(f)))

# ─────────────────────────────────────────────
# STEP 2: Doublet Detection
# ─────────────────────────────────────────────
combined_df <- data.frame()

for (i in seq_along(seurat_list)) {
  seurat_obj <- seurat_list[[i]]
  obj_name   <- names(seurat_list)[i]
  
  seurat_obj.sce               <- as.SingleCellExperiment(seurat_obj)
  seurat_obj.sce               <- scDblFinder(seurat_obj.sce)
  seurat_obj$scDblFinder.class <- seurat_obj.sce@colData$scDblFinder.class
  rm(seurat_obj.sce)
  
  temp_df <- data.frame(
    scDblFinder = seurat_obj$scDblFinder.class,
    SampleID    = seurat_obj$SampleID,
    Sample      = obj_name
  )
  combined_df <- rbind(combined_df, temp_df)
  
  seurat_list[[i]] <- subset(seurat_obj, subset = scDblFinder.class == "singlet")
}

# Plot singlet/doublet counts per sample
plot_df <- combined_df %>%
  group_by(SampleID, scDblFinder) %>%
  summarise(cell_number = n(), .groups = "drop")

nature_colors     <- DiscretePalette_scCustomize(num_colors = 26, palette = "alphabet")
n_samples         <- length(unique(plot_df$SampleID))
individual_colors <- nature_colors[1:n_samples]

ggplot(plot_df, aes(fill = SampleID, y = cell_number, x = scDblFinder)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = individual_colors) +
  labs(x = "", y = "Cell Number", title = "Singlet/Doublet Counts per Sample") +
  theme_bw() +
  theme(text         = element_text(size = 12),
        axis.text.x  = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"))

ggsave(file.path(out_dir, "Singlets_Doublets_per_Sample.png"), width = 10, height = 6, dpi = 300)

# ─────────────────────────────────────────────
# STEP 3: Merge into one Seurat object
# ─────────────────────────────────────────────
NHD_complete <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  merge.data = TRUE
)

# ─────────────────────────────────────────────
# STEP 4: Run Azimuth
# ─────────────────────────────────────────────
NHD_complete <- RunAzimuth(NHD_complete, reference = "humancortexref", assay = "RNA")

Idents(NHD_complete) <- NHD_complete$predicted.subclass

NHD_complete <- RenameIdents(object = NHD_complete,
                             "L5 ET"      = "Neuron_Ex",
                             "Pvalb"      = "Neuron_Inh",
                             "Vip"        = "Neuron_Inh",
                             "Lamp5"      = "Neuron_Inh",
                             "Sncg"       = "Neuron_Inh",
                             "VLMC"       = "Pericytes",
                             "L2/3 IT"    = "Neuron_Ex",
                             "L5 IT"      = "Neuron_Ex",
                             "L5/6 NP"    = "Neuron_Ex",
                             "L6 CT"      = "Neuron_Ex",
                             "L6 IT"      = "Neuron_Ex",
                             "L6 IT Car3" = "Neuron_Ex",
                             "Sst"        = "Neuron_Inh",
                             "L6b"        = "Neuron_Ex",
                             "Sst Chodl"  = "Neuron_Inh")

NHD_complete$new_annotation <- Idents(NHD_complete)

# Save annotation map before SCTransform
annotation_map <- data.frame(
  cell               = colnames(NHD_complete),
  new_annotation     = as.character(NHD_complete$new_annotation),
  predicted.subclass = as.character(NHD_complete$predicted.subclass)
)

saveRDS(NHD_complete, file = file.path(out_dir, "NHD_QC_Azimuth.rds"))

# ─────────────────────────────────────────────
# STEP 5: SCTransform per sample
# ─────────────────────────────────────────────
options(future.globals.maxSize = 4000 * 1024^2)

processed_seurat_list <- SplitObject(NHD_complete, split.by = "SampleID")
processed_seurat_list <- lapply(processed_seurat_list, SCTransform, assay = "RNA", vst.flavor = "v2")

# ─────────────────────────────────────────────
# STEP 6: Select integration features (no mito/ribo)
# ─────────────────────────────────────────────
rm_pat <- "^(MT-|RPL|RPS|MRPL|MRPS)"

all.features <- SelectIntegrationFeatures(
  object.list = processed_seurat_list,
  nfeatures = 5000
)

brain.features <- all.features[!grepl(rm_pat, all.features, ignore.case = TRUE)]
brain.features <- head(brain.features, 3000)

# ─────────────────────────────────────────────
# STEP 7: Merge SCTransformed list
# ─────────────────────────────────────────────
NHD_harmony <- merge(
  x = processed_seurat_list[[1]],
  y = processed_seurat_list[2:length(processed_seurat_list)],
  merge.data = TRUE
)

# Transfer Azimuth annotations
NHD_harmony$new_annotation     <- annotation_map$new_annotation[match(colnames(NHD_harmony), annotation_map$cell)]
NHD_harmony$predicted.subclass <- annotation_map$predicted.subclass[match(colnames(NHD_harmony), annotation_map$cell)]

cat("NAs in new_annotation:",     sum(is.na(NHD_harmony$new_annotation)),     "\n")
cat("NAs in predicted.subclass:", sum(is.na(NHD_harmony$predicted.subclass)), "\n")

# ─────────────────────────────────────────────
# STEP 8: Dimensionality reduction + Harmony
# ─────────────────────────────────────────────
VariableFeatures(NHD_harmony) <- brain.features

NHD_harmony <- NHD_harmony %>%
  RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(dims = 1:30)

DefaultAssay(NHD_harmony) <- "SCT"
NHD_harmony <- harmony::RunHarmony(NHD_harmony, group.by.vars = "SampleID")
NHD_harmony <- NHD_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = c(0.6, 0.8, 1))

# ─────────────────────────────────────────────
# STEP 9: Save
# ─────────────────────────────────────────────
saveRDS(NHD_harmony, file = file.path(out_dir, "NHD_QC_harmony.rds"))
gc()
