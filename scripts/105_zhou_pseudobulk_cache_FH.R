#!/usr/bin/env Rscript
# =============================================================================
# 105_zhou_pseudobulk_cache_FH.R — One load of the re-processed Zhou 2023 object, cached as per-cell-type pseudobulk count matrices + the per-sample metadata.
# -----------------------------------------------------------------------------
# Why. Three co-authors asked for the cross-cohort contrast without donor 861 (our NHD donor
# = Zhou's NHD3): if the direction consistency we report is carried by the same brain
# measured twice, it is not external support. Running that contrast — and its two
# companions (all three donors; 861 alone) — needs the same pseudobulk inputs three times,
# so the 4.2 GB object is loaded once here and the arms (106_…) read the cache.
#
# Faithful to manuscript_7fig/scripts/04_zhou_cross_cohort_validation.R (the shipped Zhou
# DGE): same object, same RNA counts, same MT/ribosomal exclusion, same per-sample
# pseudobulk (SampleID with >= 10 cells of the type), same cell types. The gene-detection
# filter (>= 5 counts in >= 2 samples) is applied per arm in 106, because it depends on
# which samples are in the design.
#
# Output: data/_cache_zhou_pseudobulk_FH.rds  = list(pb = list(<CT> = genes x samples),
#         samples = data.frame(SampleID, Condition, <every per-sample metadata column>),
#         cells_per_sample_ct = data.frame(SampleID, cell_type, n_cells), provenance)
# Run:  Rscript scripts/105_zhou_pseudobulk_cache_FH.R    (~10 min, ~12 GB RAM)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot(!is.na(PROJ))
RB   <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
ZOBJ <- file.path(PROJ, "NHD_repo", "NHD_complete_harmony+PMI.rds")
OUT  <- file.path(RB, "data", "_cache_zhou_pseudobulk_FH.rds")
stopifnot(file.exists(ZOBJ))
RM_PATTERN <- "^(MT-|RPL|RPS|MRPL|MRPS)"                    # as in script 04
CTS <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")   # as in script 04
say <- function(...) cat(sprintf(...), "\n")

if (file.exists(OUT) && file.mtime(OUT) > file.mtime(ZOBJ)) {
  say("cache is newer than the Zhou object — nothing to do (%s)", OUT); quit(save = "no")
}
say("loading %s (%.1f GB) ...", basename(ZOBJ), file.size(ZOBJ) / 1e9)
zhou <- readRDS(ZOBJ)
DefaultAssay(zhou) <- "RNA"
zhou <- JoinLayers(zhou, assay = "RNA")
md <- zhou@meta.data
stopifnot(all(c("SampleID","Condition","new_annotation") %in% names(md)))
say("cells %d | samples %d | conditions: %s", ncol(zhou), length(unique(md$SampleID)),
    paste(names(table(md$Condition)), table(md$Condition), collapse = "; "))

# ---- per-sample metadata: every column that is constant within a sample -------------
per_sample <- md %>% group_by(SampleID) %>%
  summarise(across(everything(), ~ if (dplyr::n_distinct(.x, na.rm = TRUE) <= 1) as.character(stats::na.omit(.x)[1]) else NA_character_),   # a stray NA must not blank a constant column
            n_cells = n(), .groups = "drop")
per_sample <- per_sample[, colSums(is.na(per_sample)) < nrow(per_sample)]   # drop per-cell columns
per_sample$Condition_std <- ifelse(per_sample$Condition == "NHD", "NHD", "CON")
say("per-sample metadata columns kept: %s", paste(names(per_sample), collapse = ", "))
print(as.data.frame(per_sample %>% select(SampleID, Condition_std, n_cells, any_of(c("Donor","donor","Age","age","Sex","sex","PMI","pmi","Genotype","genotype")))))

# ---- pseudobulk per cell type (SampleID with >= 10 cells of that type) --------------
genes_keep <- setdiff(rownames(zhou), grep(RM_PATTERN, rownames(zhou), value = TRUE))
cnt <- LayerData(zhou, assay = "RNA", layer = "counts")[genes_keep, ]
stopifnot(inherits(cnt, "dgCMatrix"), max(cnt@x %% 1) == 0)
pb <- list(); cps <- list()
for (CT in CTS) {
  cells <- which(md$new_annotation == CT)
  if (length(cells) < 50) { say("  %-10s SKIP (%d cells)", CT, length(cells)); next }
  sid <- as.character(md$SampleID[cells])
  tab <- table(sid); sids <- names(tab)[tab >= 10]
  m <- sapply(sids, function(s) Matrix::rowSums(cnt[, cells[sid == s], drop = FALSE]))
  pb[[CT]] <- m
  cps[[CT]] <- data.frame(SampleID = sids, cell_type = CT, n_cells = as.integer(tab[sids]))
  say("  %-10s %d cells -> %d samples with >= 10 cells", CT, length(cells), length(sids))
}
saveRDS(list(pb = pb, samples = as.data.frame(per_sample), cells_per_sample_ct = bind_rows(cps),
             provenance = list(object = ZOBJ, object_mtime = file.mtime(ZOBJ), built = Sys.time(),
                               rm_pattern = RM_PATTERN, min_cells_ct = 50, min_cells_sample = 10)),
        OUT)
say("wrote %s", OUT)
cat("\n== sessionInfo ==\n"); print(sessionInfo())
cat("\n=== DONE ===\n", file = stderr())
