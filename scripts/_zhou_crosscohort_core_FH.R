# =============================================================================
# _zhou_crosscohort_core_FH.R — The per-panel cross-cohort statistic (ST13 schema), Sourced by 70 and 106. Reads these objects from the caller's global environment:
#   zhou_all  (Zhou DGE: cell_type, gene, log2FoldChange_apeglm, padj)
#   mast      (our MAST dual-discovery table)   TB (mast_dual dir)   is_artifact()
#   ZHOU_SIG  XMIN  logln()
# After the extraction, 70 was re-run and tables/ST13_zhou_crosscohort_FH.csv was verified
# byte-identical (md5) to the shipped file.
# -----------------------------------------------------------------------------
# ---- p-underflow floor (never emit a literal 0) -------------------------------
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
floor_p <- function(p, tag) {
  if (is.na(p)) return(p)
  if (p <= 0) { logln("  NOTE: %s p underflowed to 0 -> floored to %.3e", tag, XMIN); return(XMIN) }
  p
}

# ---- per cell-type x region core computation ----------------------------------
# x = our region avg_log2FC (full ranking cache); y = Zhou occipital lfc apeglm.
one_panel <- function(CT, region) {
  panel <- sprintf("%s %s vs Zhou OCC", CT, region)
  fullf <- file.path(TB, sprintf("cliffs_delta_full_%s_%s.csv", CT, region))
  if (!file.exists(fullf)) stop(sprintf("MISSING %s — prep step failed", fullf))

  full <- read.csv(fullf, stringsAsFactors = FALSE)
  stopifnot(all(c("gene","avg_log2FC","cliffs_delta_full") %in% colnames(full)))
  ours <- full %>% transmute(gene, our_lfc = avg_log2FC, our_delta = cliffs_delta_full) %>%
    filter(!is.na(our_lfc)) %>% distinct(gene, .keep_all = TRUE)

  # Zhou occipital for this CT; keep padj NA = not TESTED (never coerced)
  zh <- zhou_all %>% filter(cell_type == CT, !is.na(log2FoldChange_apeglm)) %>%
    transmute(gene, zhou_lfc = log2FoldChange_apeglm, zhou_padj = padj) %>%
    distinct(gene, .keep_all = TRUE)
  if (nrow(zh) == 0) stop(sprintf("no Zhou rows for cell_type '%s'", CT))

  # our discovery set (discovery==TRUE, any region for this CT; |delta|>=0.15 baked in)
  disc <- mast %>% filter(cell_type == CT, region == !!region,
                          discovery %in% c(TRUE, "TRUE")) %>% distinct(gene)

  # SHARED = full ranking cache intersect Zhou (both axes non-NA) — threshold-free
  d <- ours %>% inner_join(zh, by = "gene")
  if (nrow(d) < 5) stop(sprintf("only %d shared genes for %s — below floor", nrow(d), panel))
  d <- d %>% mutate(
    is_disc    = gene %in% disc$gene,
    zhou_sig   = !is.na(zhou_padj) & zhou_padj < ZHOU_SIG,   # NA -> FALSE (not tested)
    both_sig   = is_disc & zhou_sig,
    artifact   = is_artifact(gene, CT),
    concordant = sign(our_lfc) == sign(zhou_lfc))

  # ---- rho_full: Spearman over all shared (our avg_log2FC vs Zhou lfc) ---------
  sp_full <- suppressWarnings(cor.test(d$our_lfc, d$zhou_lfc, method = "spearman"))
  rho_full   <- unname(sp_full$estimate)
  rho_full_p <- floor_p(sp_full$p.value, sprintf("%s rho_full", panel))

  # ---- leave-one-out rho: drop the biggest |our_lfc| (outlier robustness) ------
  rho_loo <- if (nrow(d) >= 5) {
    idx <- which.max(abs(d$our_lfc))
    suppressWarnings(cor(d$our_lfc[-idx], d$zhou_lfc[-idx],
                         method = "spearman", use = "pairwise.complete.obs"))
  } else NA_real_

  # ---- both-sig subset: our discovery intersect Zhou-sig ----------------------
  bs <- d %>% filter(both_sig)
  n_both <- nrow(bs)
  sp_sig <- if (n_both >= 5) suppressWarnings(cor.test(bs$our_lfc, bs$zhou_lfc, method = "spearman")) else NULL
  kn_sig <- if (n_both >= 5) suppressWarnings(cor.test(bs$our_lfc, bs$zhou_lfc, method = "kendall"))  else NULL

  # ---- sign concordance on both-sig (exclude exact-zero effects, counted) ------
  bsd     <- bs %>% filter(sign(our_lfc) != 0, sign(zhou_lfc) != 0)
  n_dirzero <- n_both - nrow(bsd)
  n_conc  <- sum(sign(bsd$our_lfc) == sign(bsd$zhou_lfc))
  conc_fr <- if (nrow(bsd)) n_conc / nrow(bsd) else NA_real_

  # ---- empirical genome-wide baseline (never 0.5): sign-agreement over all -----
  # Zhou-tested genes with an our-region measurement (not restricted to discovery).
  gw <- ours %>% inner_join(zh, by = "gene") %>%
    filter(sign(our_lfc) != 0, sign(zhou_lfc) != 0)
  baseline   <- mean(sign(gw$our_lfc) == sign(gw$zhou_lfc))
  n_baseline <- nrow(gw)

  binom_p <- if (nrow(bsd) >= 5 && baseline > 0 && baseline < 1)
    binom.test(n_conc, nrow(bsd), p = baseline, alternative = "greater")$p.value else NA_real_
  binom_p <- floor_p(binom_p, sprintf("%s binom", panel))

  stat <- tibble::tibble(
    panel            = panel,
    n_shared         = nrow(d),
    n_both_sig       = n_both,
    rho_full         = rho_full,
    rho_full_p       = rho_full_p,
    rho_sig          = if (!is.null(sp_sig)) unname(sp_sig$estimate) else NA_real_,
    rho_sig_p        = if (!is.null(sp_sig)) floor_p(sp_sig$p.value, sprintf("%s rho_sig", panel)) else NA_real_,
    tau_sig          = if (!is.null(kn_sig)) unname(kn_sig$estimate) else NA_real_,
    tau_sig_p        = if (!is.null(kn_sig)) floor_p(kn_sig$p.value, sprintf("%s tau_sig", panel)) else NA_real_,
    n_concord        = n_conc,
    concord_frac     = conc_fr,
    binom_p          = binom_p,
    rho_loo          = rho_loo,
    # ---- extra auditable columns (empirical-baseline provenance) ----
    baseline_concord = baseline,
    n_baseline_genes = n_baseline,
    n_dir_zero       = n_dirzero,
    cell_type        = CT,
    region           = region)

  list(d = d, stat = stat)
}
