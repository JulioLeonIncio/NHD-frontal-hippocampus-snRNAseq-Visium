#!/usr/bin/env Rscript
# =============================================================================
# 01_zhou_QC_Azimuth_Harmony.R — Zhou et al. 2023 occipital-cortex cohort (GEO GSE190015): GEO metadata, per-sample mitochondrial and doublet filtering, merge, PMI from the paper's supplementary table, Azimuth annotation, SCTransform
#   per individual, integration features without mitochondrial/ribosomal genes, PCA and Harmony. PMI is carried as a
#   covariate because the original study used it. Writes NHD_complete_Azimuth+PMI.rds and NHD_complete_harmony+PMI.rds.
# Paper: Human early-onset dementia caused by DAP12 deficiency reveals a unique signature of dysregulated microglia
#   (https://pmc.ncbi.nlm.nih.gov/articles/PMC9992145/).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

library(GEOquery)
library(Biobase)
library(Seurat)
library(Azimuth)
library(scDblFinder)
library(scCustomize)
library(readxl)
library(purrr)
library(dplyr)
library(ggplot2)
library(future)

base_dir  <- Sys.getenv("NHD_ZHOU_DIR")   # holds Downloads/ (GEO GSE190015 objects) and QC+PMI/ (output)
if (!nzchar(base_dir)) stop("set NHD_ZHOU_DIR to the folder holding the Zhou et al. 2023 GEO downloads")
input_dir <- file.path(base_dir, "Repository_QC")
out_dir   <- file.path(base_dir, "QC+PMI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────
# STEP 1: Get metadata from GEO
# ─────────────────────────────────────────────
gse  <- getGEO("GSE190015", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[2]]))

# ─────────────────────────────────────────────
# STEP 2: Load all RDS files
# ─────────────────────────────────────────────
rds_files <- list.files(file.path(base_dir, "Downloads"),
                        pattern = "\\.rds$", full.names = TRUE)
names(rds_files) <- sub("\\.rds$", "", basename(rds_files))

gsm_ids   <- sub("_.*", "", names(rds_files))
unmatched <- names(rds_files)[!gsm_ids %in% meta$geo_accession]
if (length(unmatched) > 0) warning("Unmatched files: ", paste(unmatched, collapse = ", "))

# ─────────────────────────────────────────────
# STEP 3: Read RDS, attach GEO metadata, remove unwanted individuals
# ─────────────────────────────────────────────
individuals_to_remove <- c("Control1, rep2", "NHD1, rep2")

seu_list <- imap(rds_files, function(path, code) {
  obj <- readRDS(path)
  
  gsm <- sub("_.*", "", code)
  row <- meta[meta$geo_accession == gsm, ]
  
  obj$SampleID     <- row$title
  obj$IndividualID <- row$`sample identifier:ch1`
  obj$Genotype     <- row$`dap12 genotype:ch1`
  obj$Age          <- row$`age:ch1`
  obj$Sex          <- row$`Sex:ch1`
  obj$Tissue       <- row$`tissue:ch1`
  obj$Batch        <- row$`batch:ch1`
  obj$Condition    <- row$`disease state:ch1`
  obj$geo_code     <- row$geo_accession
  
  obj
})

seu_list <- Filter(function(obj) {
  !obj$SampleID[1] %in% individuals_to_remove
}, seu_list)

cat("Remaining samples after removal:\n")
print(sapply(seu_list, function(x) x$SampleID[1]))

# ─────────────────────────────────────────────
# STEP 4: Mitochondrial percentage
# ─────────────────────────────────────────────
for (i in seq_along(seu_list)) {
  seurat_obj <- seu_list[[i]]
  
  mito.genes   <- grep("^MT-", rownames(seurat_obj), value = TRUE)
  mito.counts  <- colSums(GetAssayData(seurat_obj, layer = "counts")[mito.genes, , drop = FALSE])
  total.counts <- colSums(GetAssayData(seurat_obj, layer = "counts"))
  
  seurat_obj$percent_mito <- ifelse(total.counts > 0, (mito.counts / total.counts) * 100, 0)
  Idents(seurat_obj) <- seurat_obj$SampleID
  
  seu_list[[i]] <- seurat_obj
}

# ─────────────────────────────────────────────
# STEP 4b: Doublet detection
# ─────────────────────────────────────────────
combined_df <- data.frame()

for (i in seq_along(seu_list)) {
  seurat_obj <- seu_list[[i]]
  obj_name   <- names(seu_list)[i]
  
  seurat_obj.sce               <- as.SingleCellExperiment(seurat_obj)
  seurat_obj.sce               <- scDblFinder(seurat_obj.sce)
  seurat_obj$scDblFinder.class <- seurat_obj.sce@colData$scDblFinder.class
  rm(seurat_obj.sce)
  
  temp_df <- data.frame(
    scDblFinder = seurat_obj$scDblFinder.class,
    SampleID    = seurat_obj$SampleID,
    Condition   = seurat_obj$Condition,
    Sample      = obj_name
  )
  combined_df <- rbind(combined_df, temp_df)
  
  seu_list[[i]] <- subset(seurat_obj, subset = scDblFinder.class == "singlet")
}

plot_df <- combined_df %>%
  group_by(SampleID, Condition, scDblFinder) %>%
  summarise(cell_number = n(), .groups = "drop")

nature_colors     <- DiscretePalette_scCustomize(num_colors = 30, palette = "glasbey")
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
# STEP 5: Merge into one Seurat object
# ─────────────────────────────────────────────
NHD_complete <- merge(
  x = seu_list[[1]],
  y = seu_list[2:length(seu_list)],
  add.cell.ids = names(seu_list),
  merge.data = TRUE
)

# ─────────────────────────────────────────────
# STEP 6: QC filtering
# ─────────────────────────────────────────────
NHD_complete <- subset(NHD_complete,
                       subset = nFeature_RNA > 350 &
                         nFeature_RNA < 6500 &
                         percent_mito < 10)

# ─────────────────────────────────────────────
# STEP 7: Patch PMI from supplementary table
# ─────────────────────────────────────────────
supp_path <- file.path(out_dir, "NIHMS1867024-supplement-Supp_Tables_1_to_6.xlsx")
demo <- read_excel(supp_path, sheet = "Table 1_Demographic", skip = 3, col_names = TRUE)
colnames(demo) <- c("SampleIdentifier", "SampleSource", "Sex", "Age",
                    "APOEGenotype", "DAP12Genotype", "Tissue", "ceradsc",
                    "Education", "Braak", "PMI", "SampleID_snRNAseq", "Batch_snRNAseq")
demo <- demo[!is.na(demo$SampleID_snRNAseq), ]
cat("Samples in supplementary table:\n")
print(demo[, c("SampleID_snRNAseq", "PMI")])

NHD_complete@meta.data$SampleID_clean <- gsub(", rep\\d+", "", NHD_complete@meta.data$SampleID)
pmi_lookup <- setNames(as.numeric(demo$PMI), demo$SampleID_snRNAseq)
NHD_complete@meta.data$PMI <- pmi_lookup[NHD_complete@meta.data$SampleID_clean]

cat("NAs in PMI:", sum(is.na(NHD_complete@meta.data$PMI)), "\n")
cat("PMI per sample:\n")
print(tapply(NHD_complete@meta.data$PMI, NHD_complete@meta.data$SampleID_clean, unique))

# ─────────────────────────────────────────────
# STEP 8: Run Azimuth
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

annotation_map <- data.frame(
  cell               = colnames(NHD_complete),
  new_annotation     = as.character(NHD_complete$new_annotation),
  predicted.subclass = as.character(NHD_complete$predicted.subclass)
)

saveRDS(NHD_complete, file = file.path(out_dir, "NHD_complete_Azimuth+PMI.rds"))

# ─────────────────────────────────────────────
# STEP 9: SCTransform per individual
# ─────────────────────────────────────────────
options(future.globals.maxSize = 4000 * 1024^2)

processed_seurat_list <- SplitObject(NHD_complete, split.by = "SampleID")
processed_seurat_list <- lapply(processed_seurat_list, SCTransform, assay = "RNA", vst.flavor = "v2")

# ─────────────────────────────────────────────
# STEP 10: Select integration features (no mito/ribo)
# ─────────────────────────────────────────────
rm_pat <- "^(MT-|RPL|RPS|MRPL|MRPS)"

all.features <- SelectIntegrationFeatures(
  object.list = processed_seurat_list,
  nfeatures = 5000
)

brain.features <- all.features[!grepl(rm_pat, all.features, ignore.case = TRUE)]
brain.features <- head(brain.features, 3000)

# ─────────────────────────────────────────────
# STEP 11: Merge SCTransformed list
# ─────────────────────────────────────────────
NHD_harmony <- merge(
  x = processed_seurat_list[[1]],
  y = processed_seurat_list[2:length(processed_seurat_list)],
  merge.data = TRUE
)

NHD_harmony$new_annotation     <- annotation_map$new_annotation[match(colnames(NHD_harmony), annotation_map$cell)]
NHD_harmony$predicted.subclass <- annotation_map$predicted.subclass[match(colnames(NHD_harmony), annotation_map$cell)]

cat("NAs in new_annotation:",     sum(is.na(NHD_harmony$new_annotation)),     "\n")
cat("NAs in predicted.subclass:", sum(is.na(NHD_harmony$predicted.subclass)), "\n")

# ─────────────────────────────────────────────
# STEP 12: Dimensionality reduction + Harmony
# ─────────────────────────────────────────────
VariableFeatures(NHD_harmony) <- brain.features

NHD_harmony <- NHD_harmony %>%
  RunPCA(assay = "SCT", verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.8) %>%
  RunUMAP(dims = 1:30)

DefaultAssay(NHD_harmony) <- "SCT"
NHD_harmony <- harmony::RunHarmony(NHD_harmony, group.by.vars = c("SampleID", "Batch"))
NHD_harmony <- NHD_harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = c(0.6, 0.8, 1))

# ─────────────────────────────────────────────
# STEP 13: Final checks and save
# ─────────────────────────────────────────────
cat("NAs in geo_code:",     sum(is.na(NHD_harmony$geo_code)),     "\n")
cat("NAs in IndividualID:", sum(is.na(NHD_harmony$IndividualID)), "\n")
cat("NAs in SampleID:",     sum(is.na(NHD_harmony$SampleID)),     "\n")
cat("NAs in PMI:",          sum(is.na(NHD_harmony$PMI)),          "\n")
cat("PMI summary:\n")
print(summary(NHD_harmony$PMI))

saveRDS(NHD_harmony, file = file.path(out_dir, "NHD_complete_harmony+PMI.rds"))
cat("Saved successfully.\n")
gc()
