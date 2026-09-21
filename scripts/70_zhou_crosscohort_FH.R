#!/usr/bin/env Rscript
# =============================================================================
# 70_zhou_crosscohort_FH.R — Zhou 2023 (PMID 36658241) cross-cohort direction- concordance stats for the Frontal+Hippocampus-only rebuild. data layer only
#   (no figures rendered here).
#
# REFRAMING: the paper is rebuilt on Frontal+Hippocampus only; our
#   Occipital (OCC) is fully removed.  Zhou 2023 is an independent external cohort
#   (3 NHD + 11 CON, occipital cortex).  We therefore no longer have a same-region
#   ("replication") axis — this is now cross-region generalization to an independent
#   occipital cohort:
#       our Frontal  (NHD vs CON)  x  Zhou Occipital (NHD vs CON)
#       our Hippo    (NHD vs CON)  x  Zhou Occipital (NHD vs CON)
#   no "OCC vs Zhou OCC" row is produced (our OCC is gone).  One Zhou NHD donor is
#   thought to overlap one of ours, but with no donor-ID map to exclude it and the
#   comparison now cross-region anyway, all 3 Zhou NHD donors are used as-is.
#
# Method — faithful port of the two existing reference methodologies, merged:
#   * output SCHEMA = ST13 (10d_v2_5H_frontal_crossref.R):
#       panel, n_shared, n_both_sig, rho_full, rho_full_p, rho_sig, rho_sig_p,
#       tau_sig, tau_sig_p, n_concord, concord_frac, binom_p, rho_loo
#   * our side  = the FH MAST caches (MAST-forward, threshold-free ranking):
#       - rho_full / rho_loo / n_shared over the full transcriptome ranking cache
#         cliffs_delta_full_<CT>_<region>.csv  (our avg_log2FC vs Zhou lfc apeglm).
#       - our "significant" set = the MAST discovery set (discovery==TRUE, which
#         already implies |Cliff's delta| >= 0.15; verified min = 0.150).  This is
#         the both-sig / concordance / Kendall subset, intersect Zhou-sig.
#   * ZHOU side = zhou_validation.rds $dge  (log2FoldChange_apeglm, padj), the
#       processed occipital NHD-vs-CON per-gene DGE cache.  Zhou padj = NA
#       (DESeq2 independent filtering) is kept as not-tested -> not significant,
#       never coerced to 0/1.
#   * binom_p is tested vs the empirical genome-wide sign-agreement baseline
#       (from _zhou_validation_core.R's `gw` logic: sign-agreement over all Zhou-
#       tested genes that also have an our-region measurement) — never 0.5, because
#       the 0.5 null is inflated by Zhou's genome-wide sign imbalance.  The
#       baseline + n_baseline_genes are written as extra auditable columns.
#
# reviewer-facing numerical rigor:
#   * p = 0 underflow: cor.test / binom.test p-values that underflow to a literal 0
#     are floored to .Machine$double.xmin and the flooring is logged; a literal 0 is
#     never written (it would break -log10 downstream).
#   * Zhou padj NA kept as "not tested" (not sig), not coerced.
#   * sign()==0 (exact-zero effect) genes excluded from direction tests, counted.
#
# Reproducibility: set.seed(42); provenance (inputs + mtimes + sessionInfo) logged.
# Output:
#   tables/ST13_zhou_crosscohort_FH.csv                      (ST13 schema)
#   tables/zhou_FH_<CT>_<region>_genes.csv                   (per-panel scatter src)
#   tables/fig2d_micro_zhou_FH_<region>_genes.csv            (Fig2d micro scatter src)
#   tables/fig3f_astro_zhou_FH_<region>_genes.csv            (Fig3f astro scatter src)
#   tables/fig4f_neuronEx_zhou_FH_<region>_genes.csv         (Fig4f neuron scatter src)
#   logs/70_zhou_crosscohort_FH.log
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr) })
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot(!is.na(PROJ), dir.exists(PROJ))
RB    <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
TB    <- file.path(RB, "tables", "mast_dual")
TDIR  <- file.path(RB, "tables")
LOGS  <- file.path(RB, "logs")
for (d in c(TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
source(file.path(RB, "scripts", "_artifact_genes.R"))     # is_artifact(gene, CT)

XMIN     <- .Machine$double.xmin   # p=0 underflow floor (never print a literal 0)
ZHOU_SIG <- 0.10                   # Zhou padj cut for "independently significant"
REGIONS  <- c("Frontal", "Hippo")
REGION_FULL <- c(Frontal = "Frontal cortex", Hippo = "Hippocampus")

# ---- provenance log sink ------------------------------------------------------
LOGF   <- file.path(LOGS, "70_zhou_crosscohort_FH.log")
logcon <- file(LOGF, open = "wt")
logln  <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = logcon) }
logln("== 70_zhou_crosscohort_FH.R  run %s ==", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logln("PROJ = %s", PROJ)
logln("CROSS-REGION generalization: our Frontal / Hippo (NHD vs CON) x Zhou OCCIPITAL (NHD vs CON)")
logln("Zhou = 3 NHD vs 11 CON, occipital cortex (PMID 36658241); NO 'OCC vs Zhou' row (our OCC removed)")

# ---- inputs -------------------------------------------------------------------
MASTF <- file.path(TB, "MAST_dual_discovery_all.csv")
ZHOUF <- file.path(PROJ, "manuscript_7fig", "data", "zhou_validation.rds")
for (f in c(MASTF, ZHOUF))
  if (!file.exists(f)) stop(sprintf("MISSING %s — prep step failed", f))
logln("input %-38s  mtime %s", basename(MASTF), format(file.mtime(MASTF), "%Y-%m-%d %H:%M:%S"))
logln("input %-38s  mtime %s", basename(ZHOUF), format(file.mtime(ZHOUF), "%Y-%m-%d %H:%M:%S"))

mast <- read.csv(MASTF, stringsAsFactors = FALSE)
stopifnot(all(c("cell_type","region","gene","avg_log2FC","cliffs_delta","discovery") %in% colnames(mast)))
# Guard: no OCC may leak into our-side inputs (rebuild is Frontal+Hippo only).
if (any(grepl("OCC|Occipital", mast$region, ignore.case = TRUE)))
  stop("OCC detected in our MAST inputs — rebuild must be Frontal+Hippo only")

zhou_all <- readRDS(ZHOUF)$dge
stopifnot(all(c("cell_type","gene","log2FoldChange_apeglm","padj") %in% colnames(zhou_all)))
logln("Zhou dge: %d rows across cell types { %s }", nrow(zhou_all),
      paste(sort(unique(zhou_all$cell_type)), collapse = ", "))

# Cell types present in both our data and Zhou
CTS <- intersect(sort(unique(mast$cell_type)), sort(unique(zhou_all$cell_type)))
logln("cell types in both cohorts: %s", paste(CTS, collapse = ", "))

# per-panel computation lives in _zhou_crosscohort_core_FH.R (shared with 106, the
# donor-861 exclusion arms)
source(file.path(RB, "scripts", "_zhou_crosscohort_core_FH.R"))

# ---- run all cell types x {Frontal, Hippo} ------------------------------------
logln("\n== computing panels ==")
all_stats <- list(); k <- 0L
# scatter-source config: which CTs feed which figure-panel scatter csv
fig_stub <- c("Micro-PVM" = "fig2d_micro", "Astro" = "fig3f_astro", "Neuron_Ex" = "fig4f_neuronEx")

for (region in REGIONS) {
  for (CT in CTS) {
    res  <- one_panel(CT, region)
    d    <- res$d; s <- res$stat; k <- k + 1L
    all_stats[[k]] <- s

    # per-panel gene table (scatter source for every CT)
    gsafe <- gsub("[^A-Za-z0-9_.-]", "_", CT)
    write.csv(d, file.path(TDIR, sprintf("zhou_FH_%s_%s_genes.csv", gsafe, region)), row.names = FALSE)
    # figure-panel scatter source for the three headline CTs
    if (CT %in% names(fig_stub))
      write.csv(d, file.path(TDIR, sprintf("%s_zhou_FH_%s_genes.csv", fig_stub[[CT]], region)),
                row.names = FALSE)

    logln("%-22s n_shared=%-5d n_both_sig=%-3d rho_full=%+.3f (p %s) conc=%d/%d (%.0f%%) base=%.1f%% (n=%d) binom_p=%s rho_loo=%+.3f",
          s$panel, s$n_shared, s$n_both_sig, s$rho_full,
          ifelse(s$rho_full_p <= XMIN, sprintf("<%.1e", XMIN), sprintf("%.1e", s$rho_full_p)),
          s$n_concord, s$n_both_sig - s$n_dir_zero, 100 * s$concord_frac,
          100 * s$baseline_concord, s$n_baseline_genes,
          ifelse(is.na(s$binom_p), "NA", ifelse(s$binom_p <= XMIN, sprintf("<%.1e", XMIN), sprintf("%.1e", s$binom_p))),
          s$rho_loo)
  }
}

stats_df <- bind_rows(all_stats)

# Guard: no OCC-vs-Zhou row must exist
if (any(grepl("OCC vs Zhou OCC|Occipital vs Zhou", stats_df$panel)))
  stop("an OCC-vs-Zhou row leaked into the output — must not exist in the FH rebuild")

# ---- write ST13-schema CSV (canonical column order first, extras appended) ----
ST13_COLS <- c("panel","n_shared","n_both_sig","rho_full","rho_full_p","rho_sig",
               "rho_sig_p","tau_sig","tau_sig_p","n_concord","concord_frac",
               "binom_p","rho_loo")
extra_cols <- c("baseline_concord","n_baseline_genes","n_dir_zero","cell_type","region")
out_df <- stats_df[, c(ST13_COLS, extra_cols)]
OUT <- file.path(TDIR, "ST13_zhou_crosscohort_FH.csv")
write.csv(out_df, OUT, row.names = FALSE)
logln("\nwrote %s  (%d rows)", OUT, nrow(out_df))

logln("\n== ST13 summary (signif 3) ==")
print(out_df %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))
capture.output(print(out_df %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))),
               file = logcon, append = TRUE)

cat("\n== sessionInfo ==\n", file = logcon)
capture.output(sessionInfo(), file = logcon, append = TRUE)
close(logcon)
cat("provenance:", LOGF, "\n")
cat("=== DONE ===\n", file = stderr())
