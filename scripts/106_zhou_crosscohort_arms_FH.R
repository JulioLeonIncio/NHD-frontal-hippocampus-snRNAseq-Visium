#!/usr/bin/env Rscript
# =============================================================================
# 106_zhou_crosscohort_arms_FH.R — The Zhou 2023 contrast without donor 861, and its two companions, with the ST13 direction-consistency statistic recomputed per arm.
# -----------------------------------------------------------------------------
# RATIONALE (Satoru, Carlo, Piero — raised independently, 2026-09): our single NHD donor is
# Zhou's NHD3 (verified from the object metadata in 105: 59-year-old female, PMI 36 h,
# DAP12 c.2T>C = Table 1's donor 861). The cross-cohort agreement we report could therefore
# be the same brain measured twice. This script asks whether the direction consistency
# survives when Zhou's NHD side is only the two independent donors:
#
#   arm all3   : NHD1 + NHD2 + NHD3 vs 11 controls  (n = 3 v 11) — must reproduce the shipped
#                zhou_validation.rds $dge (same pipeline, same object); sanity check.
#   arm no861  : NHD1 + NHD2 vs 11 controls           (n = 2 v 11) — the primary arm.
#   arm only861: NHD3 vs 11 controls                  (n = 1 v 11) — descriptive companion:
#                what the shared donor alone contributes. DESeq2 estimates dispersion from
#                the controls; no replicate on the NHD side, so it is reported as such.
#
# PIPELINE = manuscript_7fig/scripts/04_zhou_cross_cohort_validation.R, unchanged:
#   per-sample pseudobulk (cached by 105) -> genes with >= 5 counts in >= 2 samples ->
#   DESeq2 ~ Condition -> Wald results + apeglm shrinkage. The Zhou "significant" gate stays
#   padj < 0.10 (ZHOU_SIG in 70). STATISTIC = _zhou_crosscohort_core_FH.R (verbatim from 70).
#
# What to expect and how to read IT. With two NHD donors the DESeq2 power drops, so
# n_both_sig (our discovery intersect Zhou padj < 0.10) will fall; concord_frac on that
# smaller set and, above all, rho_full over all shared genes (threshold-free, power-
# insensitive) are the quantities to compare across arms. Agreement that holds in no861 is
# external support; agreement that collapses to only861 was the shared donor.
#
# Outputs (tables/):
#   zhou_arms_pseudobulk_DESeq2_FH.csv         per-arm DGE (schema of zhou_validation$dge + arm)
#   ST13_zhou_crosscohort_arms_FH.csv          ST13 schema + arm, all cell types x regions
#   zhou_arms_<arm>_<CT>_<region>_genes.csv    per-panel scatter sources
#   zhou_arms_reproduction_check_FH.csv        all3 vs shipped: per-CT Spearman of apeglm LFC
# Run:  Rscript scripts/106_zhou_crosscohort_arms_FH.R
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble); library(DESeq2); library(apeglm) })
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot(!is.na(PROJ))
RB   <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
TB   <- file.path(RB, "tables", "mast_dual")
TDIR <- file.path(RB, "tables"); LOGS <- file.path(RB, "logs")
source(file.path(RB, "scripts", "_artifact_genes.R"))
XMIN     <- .Machine$double.xmin
ZHOU_SIG <- 0.10
REGIONS  <- c("Frontal", "Hippo")
LOGF <- file.path(LOGS, "106_zhou_crosscohort_arms_FH.log"); logcon <- file(LOGF, open = "wt")
logln <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = logcon) }
logln("== 106_zhou_crosscohort_arms_FH.R  run %s ==", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
source(file.path(RB, "scripts", "_zhou_crosscohort_core_FH.R"))   # floor_p(), one_panel()

# ---- inputs -----------------------------------------------------------------
CACHE <- file.path(RB, "data", "_cache_zhou_pseudobulk_FH.rds")
MASTF <- file.path(TB, "MAST_dual_discovery_all.csv")
SHIPF <- file.path(PROJ, "manuscript_7fig", "data", "zhou_validation.rds")
for (f in c(CACHE, MASTF, SHIPF)) if (!file.exists(f)) stop("MISSING ", f)
for (f in c(CACHE, MASTF, SHIPF)) logln("input %-36s mtime %s", basename(f), format(file.mtime(f), "%Y-%m-%d %H:%M:%S"))
zc   <- readRDS(CACHE)
mast <- read.csv(MASTF, stringsAsFactors = FALSE)
if (any(grepl("OCC|Occipital", mast$region, ignore.case = TRUE))) stop("OCC in our MAST inputs")
smp  <- zc$samples
smp$Condition_std <- factor(smp$Condition_std, levels = c("CON","NHD"))

# ---- donor 861 = Zhou NHD3, identified from the metadata, not from a note ------------
nhd <- smp[smp$Condition_std == "NHD", ]
d861 <- nhd$SampleID[nhd$Age == "59" & nhd$Sex == "female" & nhd$PMI == "36"]
stopifnot("donor 861 must match exactly one Zhou NHD sample (59 y, female, PMI 36 h)" = length(d861) == 1)
logln("donor 861 = Zhou sample '%s' (%s, age %s, PMI %s h, %s)", d861,
      nhd$Sex[nhd$SampleID == d861], nhd$Age[nhd$SampleID == d861], nhd$PMI[nhd$SampleID == d861], nhd$Genotype[nhd$SampleID == d861])
cons  <- smp$SampleID[smp$Condition_std == "CON"]
other <- setdiff(nhd$SampleID, d861)                       # NHD1 (c.2T>C, F 50) and NHD2 (c.141delG, M 50)
short <- setNames(sub("^(NHD[0-9]+).*$", "\\1", nhd$SampleID), nhd$SampleID)   # "NHD1, rep1" -> "NHD1" (the comma is part of the SampleID)
ARMS <- c(list(
  all3    = smp$SampleID,
  no861   = setdiff(smp$SampleID, d861),
  only861 = c(d861, cons)),
  setNames(lapply(other, function(x) setdiff(smp$SampleID, x)), paste0("no",   short[other])),   # noNHD1, noNHD2 (2 v 11)
  setNames(lapply(other, function(x) c(x, cons)),               paste0("only", short[other])))   # onlyNHD1, onlyNHD2 (1 v 11)
for (a in names(ARMS)) logln("arm %-9s NHD = { %s } vs %d controls", a,
  paste(short[intersect(ARMS[[a]], nhd$SampleID)], collapse = " + "), sum(ARMS[[a]] %in% cons))

# ---- pseudobulk DESeq2, exactly as script 04 -----------------------------------------
run_pb <- function(CT, sids_arm) {
  pb <- zc$pb[[CT]]; if (is.null(pb)) return(NULL)
  sids <- intersect(colnames(pb), sids_arm)
  if (length(sids) < 4) return(NULL)
  pb <- pb[, sids, drop = FALSE]
  pb <- pb[rowSums(pb >= 5) >= 2, , drop = FALSE]
  if (nrow(pb) < 200) return(NULL)
  meta <- data.frame(pseudoSample = sids,
                     Condition = smp$Condition_std[match(sids, smp$SampleID)],
                     n_cells = zc$cells_per_sample_ct$n_cells[match(paste(sids, CT), paste(zc$cells_per_sample_ct$SampleID, zc$cells_per_sample_ct$cell_type))],
                     stringsAsFactors = FALSE)
  rownames(meta) <- sids
  meta$Condition <- factor(meta$Condition, levels = c("CON","NHD"))
  # script 04 required >= 2 samples per condition; the only861 arm has 1 NHD by design and is
  # allowed through (dispersion comes from the 11 controls) — flagged in the output.
  if (sum(meta$Condition == "CON") < 2 || sum(meta$Condition == "NHD") < 1) return(NULL)
  dds <- DESeqDataSetFromMatrix(countData = round(pb), colData = meta, design = ~ Condition)
  dds <- DESeq(dds, quiet = TRUE)
  res_un <- as.data.frame(results(dds, contrast = c("Condition","NHD","CON"))) %>% rownames_to_column("gene")
  res_sh <- tryCatch(
    as.data.frame(lfcShrink(dds, coef = "Condition_NHD_vs_CON", type = "apeglm", quiet = TRUE)) %>%
      rownames_to_column("gene") %>% select(gene, log2FoldChange_apeglm = log2FoldChange, lfcSE_apeglm = lfcSE),
    error = function(e) { logln("  apeglm FAILED for %s: %s", CT, conditionMessage(e))
      data.frame(gene = res_un$gene, log2FoldChange_apeglm = NA_real_, lfcSE_apeglm = NA_real_) })
  res_un %>% left_join(res_sh, by = "gene") %>%
    mutate(cell_type = CT, n_NHD = sum(meta$Condition == "NHD"), n_CON = sum(meta$Condition == "CON"))
}

CTS <- names(zc$pb)
dge_arms <- list()
for (a in names(ARMS)) {
  logln("\n== DESeq2 arm %s ==", a)
  for (CT in CTS) {
    r <- run_pb(CT, ARMS[[a]])
    if (is.null(r)) { logln("  %-10s skipped", CT); next }
    r$arm <- a; dge_arms[[paste(a, CT)]] <- r
    logln("  %-10s genes=%-6d padj<0.10=%-5d (NHD %d vs CON %d)", CT, nrow(r), sum(r$padj < 0.10, na.rm = TRUE), r$n_NHD[1], r$n_CON[1])
  }
}
dge_arms <- bind_rows(dge_arms)
write.csv(dge_arms, file.path(TDIR, "zhou_arms_pseudobulk_DESeq2_FH.csv"), row.names = FALSE)

# ---- reproduction check: arm all3 vs the shipped Zhou DGE ----------------------------
ship <- readRDS(SHIPF)$dge
rep_chk <- lapply(CTS, function(CT) {
  a <- dge_arms %>% filter(arm == "all3", cell_type == CT) %>% select(gene, lfc_new = log2FoldChange_apeglm, padj_new = padj)
  b <- ship %>% filter(cell_type == CT) %>% select(gene, lfc_ship = log2FoldChange_apeglm, padj_ship = padj)
  m <- inner_join(a, b, by = "gene") %>% filter(!is.na(lfc_new), !is.na(lfc_ship))
  data.frame(cell_type = CT, n_genes_new = nrow(a), n_genes_shipped = nrow(b), n_common = nrow(m),
             spearman_lfc = if (nrow(m) > 5) cor(m$lfc_new, m$lfc_ship, method = "spearman") else NA,
             max_abs_lfc_diff = if (nrow(m)) max(abs(m$lfc_new - m$lfc_ship)) else NA,
             sig_new = sum(a$padj_new < 0.10, na.rm = TRUE), sig_ship = sum(b$padj_ship < 0.10, na.rm = TRUE))
})
rep_chk <- bind_rows(rep_chk)
write.csv(rep_chk, file.path(TDIR, "zhou_arms_reproduction_check_FH.csv"), row.names = FALSE)
logln("\n== reproduction check (arm all3 vs shipped zhou_validation.rds) =="); print(rep_chk)
if (any(rep_chk$spearman_lfc < 0.99, na.rm = TRUE))
  stop("arm all3 does not reproduce the shipped Zhou DGE (Spearman < 0.99) — pipeline drift; do not proceed")

# ---- ST13 statistic per arm (identical core) ----------------------------------------
all_stats <- list()
for (a in names(ARMS)) {
  zhou_all <- dge_arms %>% filter(arm == a)          # the core reads this global
  CTS_a <- intersect(sort(unique(mast$cell_type)), sort(unique(zhou_all$cell_type)))
  logln("\n== ST13 arm %s (cell types %s) ==", a, paste(CTS_a, collapse = ", "))
  for (region in REGIONS) for (CT in CTS_a) {
    res <- one_panel(CT, region); s <- res$stat; s$arm <- a
    s$n_NHD_zhou <- zhou_all$n_NHD[zhou_all$cell_type == CT][1]
    s$n_CON_zhou <- zhou_all$n_CON[zhou_all$cell_type == CT][1]          # 10 for Micro-PVM / OPC (Control4 < 10 nuclei), else 11
    s$n_zhou_sig_arm <- sum(zhou_all$padj[zhou_all$cell_type == CT] < ZHOU_SIG, na.rm = TRUE)
    # Rho on the unshrunken DESeq2 log2FC: apeglm pulls the noisiest (lowest-n) arms hardest toward
    # zero, so shrunken rho is not comparable across arms of different n.
    mle <- zhou_all %>% filter(cell_type == CT) %>% select(gene, lfc_mle = log2FoldChange, wald = stat)
    dm  <- res$d %>% inner_join(mle, by = "gene") %>% filter(!is.na(lfc_mle))
    s$rho_full_mle  <- suppressWarnings(cor(dm$our_lfc, dm$lfc_mle, method = "spearman"))
    s$rho_full_wald <- suppressWarnings(cor(dm$our_lfc, dm$wald,    method = "spearman"))
    s$genomewide_sign_agreement <- s$baseline_concord                     # same quantity, honest name (threshold-free direction agreement)
    all_stats[[length(all_stats) + 1]] <- s   # p_floored + the floor are applied once, below, to every p column
    gsafe <- gsub("[^A-Za-z0-9_.-]", "_", CT)
    write.csv(res$d, file.path(TDIR, sprintf("zhou_arms_%s_%s_%s_genes.csv", a, gsafe, region)), row.names = FALSE)
    logln("%-8s %-22s n_shared=%-5d both_sig=%-3d rho_full=%+.3f conc=%d/%d (%.0f%%) base=%.1f%% binom_p=%s",
          a, s$panel, s$n_shared, s$n_both_sig, s$rho_full, s$n_concord, s$n_both_sig - s$n_dir_zero,
          100 * s$concord_frac, 100 * s$baseline_concord,
          ifelse(is.na(s$binom_p), "NA", sprintf("%.1e", s$binom_p)))
  }
}
stats_df <- bind_rows(all_stats)
# ---- p-underflow floor on every p column -----------------------------------
# The core's floor_p() catches only p == 0; a Spearman p of 8.7e-314 is a subnormal double that
# is below the shipped floor and printed as a measured value. Rule (ST0 convention): flag any
# p < P_FLOOR first, then floor with pmax(). P_FLOOR = .Machine$double.xmin rounded up to 8
# significant digits (2.2250739e-308) so the value survives the CSV round-trip; read as
# p < 2.2e-308, never as p = 0. The count per column is logged.
P_FLOOR <- 2.2250739e-308
stopifnot(P_FLOOR >= XMIN, as.numeric(sprintf("%.15g", P_FLOOR)) >= XMIN)
P_COLS <- c("rho_full_p", "rho_sig_p", "tau_sig_p", "binom_p")
stopifnot(all(P_COLS %in% names(stats_df)))
.below <- lapply(P_COLS, function(cc) !is.na(stats_df[[cc]]) & stats_df[[cc]] < P_FLOOR); names(.below) <- P_COLS
stats_df$p_floored <- Reduce(`|`, .below)                       # TRUE where any p column of the row was below the floor
for (cc in P_COLS) {
  logln("  p-floor %-10s: %d of %d values < %.8g -> floored (report as p < 2.2e-308)", cc, sum(.below[[cc]]), sum(!is.na(stats_df[[cc]])), P_FLOOR)
  stats_df[[cc]] <- ifelse(is.na(stats_df[[cc]]), NA_real_, pmax(stats_df[[cc]], P_FLOOR)) }
stopifnot("a p value below the floor survived" = all(vapply(P_COLS, function(cc) all(is.na(stats_df[[cc]]) | stats_df[[cc]] >= P_FLOOR), logical(1))))
logln("  rows with p_floored = TRUE: %d of %d", sum(stats_df$p_floored), nrow(stats_df))
ST13_COLS <- c("panel","n_shared","n_both_sig","rho_full","rho_full_p","rho_sig","rho_sig_p","tau_sig","tau_sig_p",
               "n_concord","concord_frac","binom_p","rho_loo")
extra <- c("baseline_concord","genomewide_sign_agreement","n_baseline_genes","n_dir_zero","cell_type","region","arm",
           "n_NHD_zhou","n_CON_zhou","n_zhou_sig_arm","rho_full_mle","rho_full_wald","p_floored")
out <- stats_df[, c(ST13_COLS, extra)]
lab <- c(all3 = "all three Zhou NHD donors vs controls (3 v 11)",
         no861 = "without donor 861: NHD1 + NHD2 vs controls (2 v 11; the primary contrast)",
         only861 = "donor 861 (Zhou NHD3) vs controls (1 v 11; no NHD replicate, dispersion from controls; descriptive)")
for (x in other) {
  lab[paste0("no",   short[x])] <- sprintf("without %s: %s vs controls (2 v 11; equal-power comparator)", short[x],
                                           paste(short[setdiff(nhd$SampleID, x)], collapse = " + "))
  lab[paste0("only", short[x])] <- sprintf("%s (%s, %s) vs controls (1 v 11; descriptive)", short[x], nhd$Genotype[nhd$SampleID == x], nhd$Sex[nhd$SampleID == x])
}
out$arm_label <- lab[out$arm]
OUTF <- file.path(TDIR, "ST13_zhou_crosscohort_arms_FH.csv")
write.csv(out, OUTF, row.names = FALSE)
logln("\nwrote %s (%d rows)", OUTF, nrow(out))
# Shipped copy = ST19
ST19 <- file.path(RB, "supplementary_tables", "ST19_zhou_crosscohort_donor861_arms.csv")
write.csv(out, ST19, row.names = FALSE); logln("wrote %s", ST19)

# ---- within-Zhou donor-donor agreement: is 861 an outlier inside Zhou's own cohort? -------
# Spearman between the three 1 v 11 (apeglm and unshrunken) profiles, per cell type.
solo <- names(ARMS)[grepl("^only", names(ARMS))]
wz <- list()
for (CT in CTS) for (i in seq_along(solo)) for (j in seq_along(solo)) if (i < j) {
  a <- dge_arms %>% filter(arm == solo[i], cell_type == CT) %>% select(gene, x = log2FoldChange_apeglm, xm = log2FoldChange)
  b <- dge_arms %>% filter(arm == solo[j], cell_type == CT) %>% select(gene, y = log2FoldChange_apeglm, ym = log2FoldChange)
  m <- inner_join(a, b, by = "gene") %>% filter(!is.na(x), !is.na(y))
  wz[[length(wz) + 1]] <- data.frame(cell_type = CT, donor_a = sub("^only", "", solo[i]), donor_b = sub("^only", "", solo[j]),
    n_genes = nrow(m), rho_apeglm = cor(m$x, m$y, method = "spearman"), rho_mle = cor(m$xm, m$ym, method = "spearman"))
}
wz <- bind_rows(wz); wz$donor_a <- sub("NHD3", "NHD3 (861)", wz$donor_a); wz$donor_b <- sub("NHD3", "NHD3 (861)", wz$donor_b)
write.csv(wz, file.path(TDIR, "zhou_arms_within_cohort_donor_concordance_FH.csv"), row.names = FALSE)
write.csv(wz, file.path(RB, "supplementary_tables", "ST19b_zhou_within_cohort_donor_concordance.csv"), row.names = FALSE)
logln("\n== within-Zhou donor-donor Spearman (1 v 11 profiles) =="); print(wz %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ---- side-by-side for the three headline cell types ---------------------------------
hl <- out %>% filter(cell_type %in% c("Micro-PVM","Astro","Neuron_Ex")) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3))) %>%
  select(cell_type, region, arm, n_NHD_zhou, n_CON_zhou, n_zhou_sig_arm, rho_full, rho_full_mle, n_both_sig, concord_frac, genomewide_sign_agreement) %>%
  arrange(cell_type, region, factor(arm, levels = names(ARMS)))
logln("\n== headline cell types, three arms =="); print(as.data.frame(hl), row.names = FALSE)
capture.output(print(as.data.frame(hl), row.names = FALSE), file = logcon, append = TRUE)
cat("\n== sessionInfo ==\n", file = logcon); capture.output(sessionInfo(), file = logcon, append = TRUE); close(logcon)
cat("\n=== DONE ===\n", file = stderr())
