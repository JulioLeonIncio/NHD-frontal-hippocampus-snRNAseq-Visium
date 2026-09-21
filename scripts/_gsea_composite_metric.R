# =============================================================================
# _gsea_composite_metric.R — Shared GSEA preranking metric (single source of truth, sourced by both 44_gsea_MAST_cliffs.R and _run_3C_states_MAST_GSEA.R so the metric is computed identically everywhere).
# -----------------------------------------------------------------------------
# METRIC = signed avg_log2FC, winsorized at the 1st/99th percentile within each
# cell_type x region comparison's ranked gene list (so a few extreme/soup FC
# genes cannot dominate the tails). FC-ONLY — no Cliff's delta, no p-value
# (pseudoreplication explicitly excluded; the metric is a pure effect size).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

# canonical metric tag written to every GSEA output table
GSEA_METRIC_NAME <- "log2fc_winsorized"

# ----------------------------------------------------------------------------
# Shared pathway PRE-FILTER: drop every gene set that contains
# zero differentially-expressed (DGE / discovery) genes for the comparison before
# fgsea. A set with no DE member cannot be genuinely enriched and only inflates
# the BH multiple-testing denominator. Applied identically across every GSEA in
# the paper. `sets` = named list of gene-character-vectors; `dge_genes` = the
# comparison's discovery/DE gene symbols. Returns the filtered set list.
filter_sets_by_dge <- function(sets, dge_genes, label = NULL, verbose = TRUE) {
  dge_genes <- unique(dge_genes[!is.na(dge_genes)])
  if (!length(dge_genes)) {           # no DE genes -> nothing can be enriched
    if (verbose) message(sprintf("  %sno DGE genes -> all %d pathways dropped",
                                 if (is.null(label)) "" else paste0(label, ": "),
                                 length(sets)))
    return(sets[0])
  }
  keep <- vapply(sets, function(g) any(g %in% dge_genes), logical(1))
  if (verbose)
    message(sprintf("  %spathways with >=1 DGE gene: %d / %d (dropped %d)",
                    if (is.null(label)) "" else paste0(label, ": "),
                    sum(keep), length(sets), length(sets) - sum(keep)))
  sets[keep]
}

# winsorize a numeric vector at its own 1st/99th percentile (NA-safe).
winsorize_lfc <- function(lfc) {
  q <- stats::quantile(lfc, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  pmin(pmax(lfc, q[1]), q[2])
}

# preranking metric for one comparison's ranked gene list: signed avg_log2FC
# winsorized at the 1/99th pct. `avg_log2FC` is a per-gene vector; genes with NA
# become NA and should be dropped by the caller before fgsea.
fc_rank <- function(avg_log2FC) winsorize_lfc(avg_log2FC)
