#!/usr/bin/env Rscript
# =============================================================================
# 2J0_micro_states_cache_FH.R — Shared microglial-state cache for Figure 2.
# -----------------------------------------------------------------------------
# Builds one embedded + clustered bona-fide Micro-PVM object and caches it, so
# the UMAP panel (2A), the cluster-prevalence panel and the cluster-marker DGE
# all describe the same clusters. Before this, 2A re-embedded and re-clustered
# from the 2.9 GB atlas on every run, which is both slow and a drift hazard:
# any downstream panel that re-derived clusters separately could disagree with
# the UMAP it is supposed to annotate.
#
# Embedding parameters are a byte-faithful copy of 2A_micro_UMAP_FH.R
# (bona-fide Micro-PVM only -> SCT -> PCA 30 -> harmony(SampleID) ->
#  FindNeighbors 1:30 -> FindClusters res 0.4 -> UMAP 1:30, set.seed(42)).
#
# Cache is mtime-guarded against the atlas, so a rebuilt atlas invalidates it.
# Output: data/_cache_micro_states_FH.rds  (Seurat object, ~3k nuclei)
#         tables/micro_states/micro_cluster_counts_FH.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
CACHE <- file.path(PROJ, "data", "_cache_micro_states_FH.rds")
TBL   <- file.path(PROJ, "tables", "micro_states")
dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
stopifnot("MISSING atlas" = file.exists(ATLAS))

stale <- !file.exists(CACHE) || file.info(CACHE)$mtime < file.info(ATLAS)$mtime
if (stale) {
  cat("== cache miss/stale: rebuilding microglial-state embedding from atlas ==\n")
  suppressPackageStartupMessages(library(harmony))
  atl <- readRDS(ATLAS)
  o   <- subset(atl, subset = new_annotation == "Micro-PVM")   # BONA-FIDE only
  rm(atl); gc()
  cat(sprintf("bona-fide Micro-PVM nuclei: %d\n", ncol(o)))
  DefaultAssay(o) <- "SCT"
  stopifnot("SCT scale.data empty — re-run SCTransform/PrepSCT before embedding" =
              nrow(GetAssayData(o, assay = "SCT", layer = "scale.data")) > 0)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- harmony::RunHarmony(o, group.by.vars = "SampleID", verbose = FALSE)
  o <- FindNeighbors(o, reduction = "harmony", dims = 1:30, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.4, verbose = FALSE)
  o <- RunUMAP(o, reduction = "harmony", dims = 1:30, verbose = FALSE)
  saveRDS(o, CACHE)
  cat(sprintf("cached -> %s\n", CACHE))
} else {
  cat("== loading cached microglial-state object ==\n")
  o <- readRDS(CACHE)
}

md <- o@meta.data
md$cluster <- factor(as.integer(as.character(md$seurat_clusters)))
md$Condition <- factor(md$Condition, levels = c("CON", "NHD"))
md$Region    <- factor(md$Region,    levels = REGION_ORDER)

cat(sprintf("\nmicroglial nuclei: %d   clusters: %s\n",
            nrow(md), paste(levels(md$cluster), collapse = ",")))

cnt <- md %>%
  count(cluster, Condition, Region, name = "n") %>%
  tidyr::complete(cluster, Condition, Region, fill = list(n = 0L))
tot <- md %>% count(cluster, name = "n_total")

cat("\n=== nuclei per cluster (total) ===\n")
print(as.data.frame(tot))
cat("\n=== cluster x condition ===\n")
print(table(md$cluster, md$Condition))
cat("\n=== cluster x region ===\n")
print(table(md$cluster, md$Region))
cat("\n=== cluster x condition x region ===\n")
print(ftable(table(md$cluster, md$Condition, md$Region)))

out <- cnt %>% left_join(tot, by = "cluster")
write.csv(out, file.path(TBL, "micro_cluster_counts_FH.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s\n", file.path(TBL, "micro_cluster_counts_FH.csv")))

small <- tot$cluster[tot$n_total < 10]
if (length(small)) {
  cat(sprintf("\n** clusters with <10 nuclei (drop candidates): %s **\n",
              paste(small, collapse = ", ")))
} else cat("\n** no cluster has <10 nuclei **\n")

cat("\n=== DONE: 2J0_micro_states_cache_FH ===\n")
