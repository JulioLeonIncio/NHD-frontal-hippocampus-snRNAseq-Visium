#!/usr/bin/env Rscript
# =============================================================================
# 131_zhou_TYROBP_per_donor_FH.R — TYROBP, TREM2 and CSF1R detection in the microglia of every donor of the re-processed Zhou et al. 2023 cohort (GEO GSE190015; 3 NHD, 11 controls).
# -----------------------------------------------------------------------------
# Why. TYROBP is the mutated gene, so it is left out of every immune score: our NHD donor (c.2T>C,
#   start-loss) retains the transcript (Results), which therefore does not report DAP12 protein. The Zhou
#   et al. NHD donors carry two alleles (NHD1 and NHD3 = donor 861: c.2T>C; NHD2: c.141delG); this table
#   records whether their microglia retain the transcript, with TREM2 and CSF1R as within-donor references and the
#   median depth of each donor as the depth check for the detection percentages.
# What, per donor (the GEO sample title; the second libraries of Control1 and NHD1 were removed during
#   re-processing, so each donor is one library), over the nuclei whose Azimuth annotation contains
#   "Micro" (in this object the Micro-PVM class only: microglia with perivascular macrophages, 2,884 nuclei)
#   — every value from the raw RNA counts (SCT-transformed values fabricate detection):
#   * n_microglia; pct_<gene>_pos = % of those nuclei with >= 1 count of the gene;
#   * TYROBP_per10k = 1e4 * sum(TYROBP counts) / sum(total counts), POOLED over the donor's microglia
#     (a ratio of sums, not a mean of per-nucleus ratios, so shallow nuclei do not dominate);
#   * median_UMI = median total counts per microglial nucleus.
# In:  <NHD>/NHD_repo/NHD_complete_harmony+PMI.rds — the re-processed Zhou et al. 2023 object written by
#      zhou_reprocessing/01_zhou_QC_Azimuth_Harmony.R (the object 105 and 04 read).
# Out: tables/zhou_TYROBP_per_donor_FH.csv and its shipped copy
#      supplementary_tables/ST32_zhou_TYROBP_per_donor_FH.csv (Supplementary Data 9, sheet TYROBP_per_donor,
#      packaged by 100_package_supplementary_data_FH.R).
# Run:  Rscript scripts/131_zhou_TYROBP_per_donor_FH.R    (~5 min, ~12 GB RAM: one load of the 4.2 GB object)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot(!is.na(PROJ))
RB   <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
ZOBJ <- file.path(PROJ, "NHD_repo", "NHD_complete_harmony+PMI.rds")
stopifnot("MISSING Zhou object NHD_repo/NHD_complete_harmony+PMI.rds — run zhou_reprocessing/01_zhou_QC_Azimuth_Harmony.R" = file.exists(ZOBJ))
OUT_TAB  <- file.path(RB, "tables", "zhou_TYROBP_per_donor_FH.csv")
OUT_SHIP <- file.path(RB, "supplementary_tables", "ST32_zhou_TYROBP_per_donor_FH.csv")

cat(sprintf("== loading %s (%.1f GB)\n", basename(ZOBJ), file.size(ZOBJ) / 1e9))
z <- readRDS(ZOBJ); DefaultAssay(z) <- "RNA"; z <- JoinLayers(z, assay = "RNA")
md <- z@meta.data
DON <- "SampleID"   # the GEO sample title, as in 105 (one library per donor after re-processing)
stopifnot(all(c(DON, "Condition", "new_annotation") %in% names(md)),
          "Zhou Condition must be control / NHD" = setequal(unique(as.character(md$Condition)), c("control", "NHD")))
mic <- grepl("Micro", md$new_annotation, ignore.case = TRUE)
stopifnot("no microglial nuclei in the Zhou object" = any(mic))
cat(sprintf("== microglial nuclei: %d of %d (annotations: %s)\n", sum(mic), nrow(md),
            paste(unique(md$new_annotation[mic]), collapse = " | ")))

# Subset the sparse counts to the microglial nuclei before any per-gene extraction (never densify the whole matrix)
cnt <- GetAssayData(z, assay = "RNA", layer = "counts")[, mic, drop = FALSE]
m   <- md[mic, , drop = FALSE]
rm(z); invisible(gc())
tot <- Matrix::colSums(cnt)
stopifnot("a microglial nucleus with zero total counts" = all(tot > 0))   # the per-10k denominator cannot be 0
g_of <- function(g) { stopifnot("gene absent from the Zhou matrix" = g %in% rownames(cnt)); as.numeric(cnt[g, ]) }
d <- data.frame(donor = as.character(m[[DON]]), condition = as.character(m$Condition),
                tot = tot, TYROBP = g_of("TYROBP"), TREM2 = g_of("TREM2"), CSF1R = g_of("CSF1R"))
# arrange() sorts in the C locale (dplyr >= 1.1), so the row order is the same on every machine
out <- d %>% group_by(condition, donor) %>%
  summarise(n_microglia = n(),
            pct_TYROBP_pos = round(100 * mean(TYROBP > 0), 1),
            TYROBP_per10k  = round(1e4 * sum(TYROBP) / sum(tot), 2),   # pooled: ratio of sums over the donor's microglia
            pct_TREM2_pos  = round(100 * mean(TREM2 > 0), 1),
            pct_CSF1R_pos  = round(100 * mean(CSF1R > 0), 1),
            median_UMI     = round(median(tot)), .groups = "drop") %>%
  arrange(desc(condition), donor)
stopifnot(!anyNA(out), !anyDuplicated(out$donor))
print(as.data.frame(out), row.names = FALSE)
write.csv(out, OUT_TAB, row.names = FALSE)
write.csv(out, OUT_SHIP, row.names = FALSE)
cat(sprintf("== wrote %s and %s (%d donors)\n", OUT_TAB, OUT_SHIP, nrow(out)))
cat("\npooled: CON", sprintf("%.1f%%", 100 * mean(d$TYROBP[d$condition != "NHD"] > 0)),
    "| NHD", sprintf("%.1f%%", 100 * mean(d$TYROBP[d$condition == "NHD"] > 0)), "\n")
cat("\n== sessionInfo ==\n"); print(sessionInfo())
cat("=== DONE ===\n")
