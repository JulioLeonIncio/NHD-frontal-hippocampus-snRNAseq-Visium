#!/usr/bin/env Rscript
# =============================================================================
# 68_micro_TF_collectri_FH.R — Which transcription factors govern the MICROGLIAL NHD-vs-CON difference, in Frontal and Hippocampus.
# -----------------------------------------------------------------------------
# Why a re-analysis rather than A RE-PLOT (after two iterations
# on the existing table showed it could not answer the question):
#
#   1. wrong prior. `TF_main_celltype_activity_per_comparison.csv` was scored on the
#      GTRD regulon, which is ChIP-derived across many cell types. Only 43 of 386
#      scored TFs are detected in >=10% of frontal microglia. That is how BARHL1
#      (hindbrain), NKX2-2 (oligodendrocyte), FOXE1 (thyroid), DLX6 (interneuron)
#      and PAX7 (muscle) reached the top of a MICROGLIAL result — and how PSMB5
#      (a proteasome subunit), GTF2A2/GTF2E2 (basal machinery) and HJURP (a histone
#      chaperone) were reported as "transcription factors" at all.
#      Worse, every canonical myeloid TF is absent from that prior: SPI1, IRF8,
#      CEBPA/B, NFKB1, STAT1/2/3, MAF, MAFB, RUNX1, JUN, FOS, SALL1, MITF are all
#      detected in the cells but have no regulon to be scored against.
#      -> CollecTRI (curated, signed, myeloid-complete) replaces GTRD here.
#
#   2. TF filter applied after scoring, not before. Spending the multiple-testing
#      burden on factors the cell cannot express costs real power. The prior is now
#      restricted to expressed TFs before run_ulm.
#
#   3. depth-confounded input. The old run fed MAST avg_log2FC, a share. Microglial
#      depth is balanced in Frontal (CON 722 vs NHD 736 median UMI, ratio 1.02) but
#      NHD is 3.4x DEEPER in Hippocampus (717 vs 2,425) — and hippocampus is exactly
#      where the old table reported 31 hits against frontal's 1. That ordering is the
#      wrong way round for biology. Input is now a depth-matched statistic.
#
# Output is deliberately a table plus a diagnostic, not a figure: the point of the
# re-run is to find out whether a defensible answer exists before drawing one.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tidyr); library(decoupleR)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
TDIR <- file.path(PROJ, "tables"); DDIR <- file.path(PROJ, "data")
LOGD <- file.path(PROJ, "logs")
for (d in c(TDIR, DDIR, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Run this script twice:
#     DEPTH_MODE=lognorm_all  (default) -> the PRIMARY estimator
#     DEPTH_MODE=thin750                -> the sensitivity arm
# why the default changed. The 750-UMI floor-then-thin was the primary treatment and it
# should not have been. It is an estimator being used where a diagnostic was needed:
#   * it corrects nothing in frontal cortex, where myeloid depth is already matched between
#     conditions (696 vs 711 UMI, KS D = 0.14), while discarding 57% of those nuclei;
#   * in hippocampus it keeps 44% of the control molecules but only 25% of the NHD ones,
#     so it taxes the arm carrying the signal and shrinks effects by regression dilution;
#   * it truncates the feature universe, and three of the factors it removes -- PPARG,
#     STAT1 and CEBPB -- clear the 10% detection floor at native depth. PPARG is the DAM-2
#     lipid regulon, i.e. the paper's central claim, and at native depth it is a clean
#     measured null (+0.70 frontal, +0.02 hippocampus, neither significant).
# The two treatments agree on the answer (Pearson r = 0.97 in both regions, regression
# slope 1.03/1.06, 3 and 4 sign flips of 48) while differing five-fold in residual depth
# artifact by the control-only negative control in 96_depth_diagnostics_FH.R -- which is
# evidence that the scores are not carried by depth. `run_ulm` returns a t statistic and is
# invariant to uniform rescaling of its input, so the per-gene attenuation does not reach
# this panel. Log-normalization is also what ST1 and every other panel already use.
# [[depth-matching-is-a-falsification-test-not-an-estimator]]
DEPTH_MODE <- Sys.getenv("DEPTH_MODE", "lognorm_all")
stopifnot("DEPTH_MODE must be lognorm_all or thin750" =
            DEPTH_MODE %in% c("lognorm_all","thin750"))
SUFFIX <- if (DEPTH_MODE == "thin750") "_thin750" else ""
cat(sprintf("DEPTH_MODE = %s   (outputs get suffix '%s')\n", DEPTH_MODE, SUFFIX))

DEPTH_T  <- 750L    # the sensitivity arm's floor; see the sweep in 96_depth_diagnostics_FH.R
DET_FLOOR<- 0.10    # a TF must be detected in >=10% of the cells to be scored
MINSIZE  <- 10L     # was 5. A z from 5 targets cannot sit beside a z from 85; every
                    # displayed factor now needs at least 10 scored targets.
# Housekeeping / stress leave-out. Several factors
# were being carried by heat-shock and nuclear-quality transcripts rather than by
# biology -- NFAT5 Frontal fell +2.07 -> +0.54 when these were removed, and CREB1's
# leading targets are MALAT1, HSPA1A, TPT1, RPS27A. Every score now ships with its
# leave-out value so the robustness is on the record instead of being assumed.
HK_RE <- "^RP[LS]|^MT-|^HSPA|^HSPB|^HSP90|^DNAJ|^MALAT1$|^NEAT1$|^TPT1$|^EEF|^FTL$|^FTH1$|^B2M$|^ACTB$|^UB[BC]$"

# ---------------------------------------------------------------------------
# 1. CollecTRI regulon (cached — the network is data-independent)
# ---------------------------------------------------------------------------
NET <- file.path(DDIR, "TF_regulon_collectri.rds")
RAW <- file.path(DDIR, "CollecTRI_raw.csv")
if (file.exists(NET)) {
  net <- readRDS(NET); cat("CollecTRI: loaded cache\n")
} else {
  # decoupleR::get_collectri() and every OmnipathR route (collectri/transcriptional/
  # dorothea) fail on this machine with an internal map_int() error — a package-version
  # bug, not the network: the server responds fine. The table is therefore fetched
  # directly and parsed here. `weight` is CollecTRI's signed mode of regulation
  # (+1 activation / -1 repression), which is exactly what run_ulm expects as `mor`.
  stopifnot("MISSING data/CollecTRI_raw.csv — fetch it from rescued.omnipathdb.org" =
              file.exists(RAW))
  raw <- read.csv(RAW, stringsAsFactors = FALSE)
  net <- raw[, c("source","target","weight")]
  names(net)[3] <- "mor"
  net <- net[!is.na(net$mor) & nzchar(net$source) & nzchar(net$target), ]
  net <- net[!duplicated(net[, c("source","target")]), ]
  saveRDS(net, NET)
  cat("CollecTRI: parsed from the raw table\n")
}
cat(sprintf("CollecTRI: %d edges, %d TFs\n", nrow(net), dplyr::n_distinct(net$source)))
CANON <- c("SPI1","IRF8","RUNX1","MEF2C","MAF","MAFB","CEBPA","CEBPB","NFKB1","RELA",
           "STAT1","STAT2","STAT3","TFEB","NR4A1","JUN","JUNB","FOS","ATF3","IRF1",
           "IRF7","EGR1","KLF4","SALL1","SMAD3","BHLHE40","MITF","TAL1")
cat("canonical myeloid TFs present in CollecTRI:",
    sum(CANON %in% net$source), "/", length(CANON), "\n")
cat("  missing:", paste(setdiff(CANON, net$source), collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 2. microglia, depth-matched, per region
# ---------------------------------------------------------------------------
CACHE <- file.path(DDIR, "_cache_micro_states_clean_FH.rds")
stopifnot("MISSING microglia cache" = file.exists(CACHE))
o <- readRDS(CACHE)
DefaultAssay(o) <- "RNA"
if (length(SeuratObject::Layers(o, assay = "RNA")) > 1) o <- SeuratObject::JoinLayers(o)
cts <- GetAssayData(o, assay = "RNA", layer = "counts")
md  <- o@meta.data
cat(sprintf("\nmicroglia cache: %d cells | regions %s\n", ncol(o),
            paste(sort(unique(as.character(md$Region))), collapse = "/")))

# =============================================================================
# Pass 1 — per-region thinned matrices and detection rates.
# The TF set must be fixed across regions. Applying the detection floor separately in
# each region silently drops a factor from one panel and renders it as blank ink: NFKB1
# is detected in 7.2% of thinned frontal cells and 10.7% of thinned hippocampal cells,
# so a per-region 10% floor scored it in one region only, and the figure then asserted
# a region difference that was a knife-edge floor artefact. The scored set is now the
# union of factors clearing the floor in either region, so both regions are scored on
# the same factors and a genuine absence is distinguishable from a threshold accident.
# =============================================================================
REGS <- c("Frontal","Hippo")
prep <- list(); gene_delta <- list()
for (rg in REGS) {
  sel <- md$Region == rg & md$Condition %in% c("CON","NHD")
  msel <- md[sel, ]
  sub_cts <- cts[, sel, drop = FALSE]
  tot <- Matrix::colSums(sub_cts)
  if (DEPTH_MODE == "thin750") {
    keep <- tot >= DEPTH_T
    m <- msel[keep, ]
    cat(sprintf("\n=== %s ===\n  thin T=%d: CON %d/%d, NHD %d/%d retained\n", rg, DEPTH_T,
        sum(keep & msel$Condition == "CON"), sum(msel$Condition == "CON"),
        sum(keep & msel$Condition == "NHD"), sum(msel$Condition == "NHD")))
    if (min(sum(m$Condition=="CON"), sum(m$Condition=="NHD")) < 40) {
      cat("  too few cells after thinning — SKIPPED\n"); next
    }
    ds <- Seurat::SampleUMI(sub_cts[, keep, drop = FALSE],
                            max.umi = DEPTH_T, upsample = FALSE, verbose = FALSE)
    rownames(ds) <- rownames(cts)
    detf <- Matrix::rowMeans(ds > 0)
  } else {
    # PRIMARY: every nucleus, library-size normalized (log1p of counts per 10,000).
    # Detection is taken from the RAW counts, so the floor means "detected", not
    # the same convention ST1 uses.
    m  <- msel
    cat(sprintf("\n=== %s ===\n  all nuclei, log-normalized: CON %d, NHD %d\n", rg,
        sum(m$Condition == "CON"), sum(m$Condition == "NHD")))
    if (min(sum(m$Condition=="CON"), sum(m$Condition=="NHD")) < 40) {
      cat("  too few cells — SKIPPED\n"); next
    }
    detf <- Matrix::rowMeans(sub_cts > 0)
    ds <- log1p(t(t(sub_cts) / pmax(tot, 1)) * 1e4)
    rownames(ds) <- rownames(cts)
  }
  # depth-matched effect size per gene: Cliff's delta, NHD vs CON
  isN <- m$Condition == "NHD"; isC <- m$Condition == "CON"
  gkeep <- detf >= 0.05
  dsg <- ds[gkeep, , drop = FALSE]
  cd <- apply(as.matrix(dsg), 1, function(v) {
    r <- rank(v); n1 <- sum(isN); n2 <- sum(isC)
    U <- sum(r[isN]) - n1*(n1+1)/2
    2*(U/(n1*n2)) - 1
  })
  stat <- matrix(cd, ncol = 1, dimnames = list(names(cd), "cliffs_delta"))
  cat(sprintf("  genes scored: %d (depth-matched Cliff's delta)\n", nrow(stat)))
  # Export the per-gene statistic. The TF scores are a weighted roll-up
  # of exactly these numbers, so the network panel (70) has to read the same vector
  # rather than recompute it from a second thinning with a different random draw --
  # otherwise the edges and the node colours come from two different datasets.
  gene_delta[[rg]] <- data.frame(gene = rownames(stat), region = rg,
                                 cliffs_delta = as.numeric(stat[, 1]),
                                 pct_detected = as.numeric(detf[gkeep]),
                                 stringsAsFactors = FALSE)

  prep[[rg]] <- list(stat = stat, detf = detf, n_CON = sum(isC), n_NHD = sum(isN))
}

# ---- the fixed, region-invariant TF set --------------------------------------
tf_any <- unique(unlist(lapply(prep, function(p)
  intersect(net$source, names(p$detf)[p$detf >= DET_FLOOR]))))
net_f  <- net %>% filter(source %in% tf_any)
cat(sprintf("\nFIXED TF set (>=%.0f%% detection in EITHER region): %d factors, %d edges\n",
            DET_FLOOR*100, length(tf_any), nrow(net_f)))
for (rg in REGS) {
  only <- setdiff(intersect(net$source, names(prep[[rg]]$detf)[prep[[rg]]$detf >= DET_FLOOR]),
                  Reduce(intersect, lapply(prep, function(p)
                    intersect(net$source, names(p$detf)[p$detf >= DET_FLOOR]))))
  if (length(only))
    cat(sprintf("  clears the floor ONLY in %s (now scored in both): %s\n",
                rg, paste(sort(only), collapse = ", ")))
}
# Factors a reader will look for that no region can score — named so the legend cannot
# claim "every canonical identity factor" when the strongest ones were never measurable
WANTED <- c("SPI1","IRF8","MEF2C","MAF","MAFB","RUNX1","TAL1","SMAD3","SALL1",
            "CEBPA","CEBPB","RELA","STAT1","JUN","FOS","EGR1","ATF3","NR4A1")
unscorable <- setdiff(intersect(WANTED, net$source), tf_any)
cat("  BELOW the detection floor in both regions (cannot be scored):",
    paste(sort(unscorable), collapse = ", "), "\n")

# =============================================================================
# Pass 2 — score the same factors in every region, twice: as-is, and with the
# housekeeping / heat-shock / nuclear-quality transcripts removed.
# =============================================================================
res <- list(); diag <- list()
for (rg in REGS) {
  stat <- prep[[rg]]$stat
  a <- decoupleR::run_ulm(mat = stat, network = net_f,
                          .source = "source", .target = "target", .mor = "mor",
                          minsize = MINSIZE)
  hk  <- grepl(HK_RE, rownames(stat))
  a_h <- decoupleR::run_ulm(mat = stat[!hk, , drop = FALSE], network = net_f,
                            .source = "source", .target = "target", .mor = "mor",
                            minsize = MINSIZE)
  cov <- net_f %>% filter(target %in% rownames(stat)) %>%
    group_by(source) %>% summarise(n_target = n(), .groups = "drop")
  a <- a %>%
    left_join(a_h %>% select(source, score_noHK = score), by = "source") %>%
    left_join(cov, by = "source")
  a$region <- rg
  res[[rg]] <- a
  cat(sprintf("  %s: %d factors scored (minsize %d), %d HK/HSP genes removed for the leave-out\n",
              rg, nrow(a), MINSIZE, sum(hk)))
  diag[[rg]] <- data.frame(region = rg, n_tf_scored = dplyr::n_distinct(a$source),
                           n_sig = sum(a$p_value < 0.05),
                           n_cells_CON = prep[[rg]]$n_CON, n_cells_NHD = prep[[rg]]$n_NHD)
}

if (!length(res)) stop("no region produced a result")
all <- bind_rows(res) %>% arrange(region, p_value)
write.csv(all, file.path(TDIR, sprintf("micro_TF_collectri%s_FH.csv", SUFFIX)), row.names = FALSE)
cat(sprintf("\nWrote tables/micro_TF_collectri%s_FH.csv\n", SUFFIX))
gd <- bind_rows(gene_delta)
write.csv(gd, file.path(TDIR, sprintf("micro_gene_cliffsdelta%s_FH.csv", SUFFIX)), row.names = FALSE)
cat(sprintf("Wrote tables/micro_gene_cliffsdelta%s_FH.csv (%d rows)\n", SUFFIX, nrow(gd)))
print(bind_rows(diag))

for (rg in names(res)) {
  cat(sprintf("\n=== %s — top TFs (CollecTRI, expressed-only, %s) ===\n", rg, DEPTH_MODE))
  t <- res[[rg]] %>% arrange(p_value) %>% head(15)
  cat(sprintf("  %-10s %8s %10s\n", "TF", "score", "p"))
  for (i in seq_len(nrow(t)))
    cat(sprintf("  %-10s %+8.2f %10.2e%s\n", t$source[i], t$score[i], t$p_value[i],
                ifelse(t$p_value[i] < 0.05, "  *", "")))
  cat("  canonical myeloid TFs in this result:\n")
  cc <- res[[rg]] %>% filter(source %in% CANON) %>% arrange(p_value)
  if (!nrow(cc)) cat("    none scored\n") else
    for (i in seq_len(nrow(cc)))
      cat(sprintf("    %-9s %+7.2f  p %.3f%s\n", cc$source[i], cc$score[i], cc$p_value[i],
                  ifelse(cc$p_value[i] < 0.05, "  *", "")))
}
writeLines(capture.output(sessionInfo()),
           file.path(LOGD, sprintf("68_micro_TF_collectri%s_FH_sessionInfo.txt", SUFFIX)))
cat("\n=== DONE ===\n")
