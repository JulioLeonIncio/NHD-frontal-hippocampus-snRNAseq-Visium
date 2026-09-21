#!/usr/bin/env Rscript
# =============================================================================
# 2L0_pvm_identify_FH.R — Is there a resolvable PVM / border-associated macrophage population inside the Micro-PVM subset? diagnostic only.
# -----------------------------------------------------------------------------
# "before we had a very small cluster of PVM, that cluster was
# very different. I want to include that one and at least show PVM vs Microglia
# proportions, including cell markers."
#
# Why this is a different question from the state-subclustering we just rejected
# ([[nhd-microglia-do-not-subcluster]]): PVM/BAM vs parenchymal microglia is a
# lineage/identity distinction with a large, canonical, well-replicated marker set
# (MRC1/CD163/LYVE1/F13A1/MS4A7/STAB1...), not a transcriptional "state". So it can
# be real even though DAM-style state clusters were not. This script tests that
# rather than assuming it: it asks whether the PVM markers form a coherent,
# separable population, and whether any cluster is PVM-dominant.
#
# Prints only — no panels, no classification is written. Cut choice comes after.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
CACHE <- file.path(PROJ, "data", "_cache_micro_states_FH.rds")
TBL   <- file.path(PROJ, "tables", "micro_states")
dir.create(TBL, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(CACHE))

# Canonical human PVM/BAM and parenchymal-microglia panels.
# TYROBP is deliberately EXCLUDED (it is the mutated gene in this disease).
PVM_MK   <- c("MRC1","CD163","LYVE1","F13A1","MS4A7","STAB1","COLEC12","DAB2",
              "SIGLEC1","PF4","CD209","MAF","APOE","MSR1","CD36","EGFL7")
MICRO_MK <- c("P2RY12","P2RY13","TMEM119","CX3CR1","GPR34","SELPLG","SALL1",
              "SLC2A5","OLFML3","SPARC","CSF1R","INPP5D","DOCK8")

o <- readRDS(CACHE)
o$cl04 <- factor(as.integer(as.character(o$seurat_clusters)))
DefaultAssay(o) <- "RNA"
if (inherits(o[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(o[["RNA"]], search = "counts")) > 1)
  o[["RNA"]] <- SeuratObject::JoinLayers(o[["RNA"]])
o <- NormalizeData(o, verbose = FALSE)

PVM_MK   <- PVM_MK[PVM_MK     %in% rownames(o)]
MICRO_MK <- MICRO_MK[MICRO_MK %in% rownames(o)]
cat("PVM markers present  :", paste(PVM_MK, collapse = ", "), "\n")
cat("Micro markers present:", paste(MICRO_MK, collapse = ", "), "\n\n")

o <- AddModuleScore(o, features = list(PVM_MK, MICRO_MK),
                    name = "IDSC", seed = 42, ctrl = 50)
md <- o@meta.data %>%
  mutate(PVM_score = IDSC1, MICRO_score = IDSC2, pvm_lead = IDSC1 - IDSC2)

cat("=== per-nucleus detection rate of the core PVM markers ===\n")
core <- intersect(c("MRC1","CD163","LYVE1","F13A1","MS4A7","STAB1","COLEC12","DAB2"), rownames(o))
det  <- FetchData(o, vars = core, layer = "data")
print(round(100 * colMeans(det > 0), 2))

cat("\n=== PVM-lead score by res-0.4 cluster ===\n")
print(md %>% group_by(cl04) %>%
        summarise(n = n(), PVM = median(PVM_score), MICRO = median(MICRO_score),
                  lead = median(pvm_lead), frac_lead_pos = mean(pvm_lead > 0)) %>%
        as.data.frame(), digits = 3)

cat("\n=== core PVM marker detection (%) by res-0.4 cluster ===\n")
d2 <- cbind(det, cl04 = md$cl04)
print(d2 %>% group_by(cl04) %>% summarise(across(all_of(core), ~ round(100*mean(.x > 0), 1)),
                                          n = n()) %>% as.data.frame())

cat("\n=== pvm_lead distribution (is it bimodal / is there a tail?) ===\n")
print(round(quantile(md$pvm_lead, c(.5,.75,.9,.95,.97,.98,.99,.995,1)), 3))
for (cut in c(0.15, 0.20, 0.25, 0.30, 0.40)) {
  k <- md$pvm_lead > cut
  cat(sprintf("cut %.2f -> %4d nuclei (%.2f%%) | MRC1+ %.0f%% CD163+ %.0f%% P2RY12+ %.0f%%\n",
              cut, sum(k), 100*mean(k),
              100*mean(det$MRC1[k] > 0), 100*mean(det$CD163[k] > 0),
              100*mean(FetchData(o, "P2RY12", layer = "data")[k, 1] > 0)))
}

# Does a finer clustering isolate them? (identity, unlike states, should hold up)
GRAPH <- grep("_snn$", names(o@graphs), value = TRUE)[1]
for (r in c(0.6, 0.8, 1.0)) {
  o <- FindClusters(o, resolution = r, graph.name = GRAPH, verbose = FALSE)
  cl <- factor(as.integer(as.character(o$seurat_clusters)))
  s  <- data.frame(cl = cl, lead = md$pvm_lead, MRC1 = det$MRC1 > 0, CD163 = det$CD163 > 0)
  cat(sprintf("\n=== resolution %.1f: PVM-lead by cluster ===\n", r))
  print(s %>% group_by(cl) %>%
          summarise(n = n(), med_lead = round(median(lead), 3),
                    MRC1_pct = round(100*mean(MRC1), 1),
                    CD163_pct = round(100*mean(CD163), 1)) %>%
          arrange(desc(med_lead)) %>% head(4) %>% as.data.frame())
}

write.csv(md %>% select(Condition, Region, SampleID, cl04, PVM_score, MICRO_score, pvm_lead),
          file.path(TBL, "micro_pvm_identity_scores_FH.csv"))
cat("\nwrote micro_pvm_identity_scores_FH.csv\n=== DONE: 2L0_pvm_identify_FH ===\n")
