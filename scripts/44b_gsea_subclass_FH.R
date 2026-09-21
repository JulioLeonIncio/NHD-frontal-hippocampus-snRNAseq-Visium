#!/usr/bin/env Rscript
# =============================================================================
# 44b_gsea_subclass_FH.R — GSEA (Hallmark) for the neuron subclass / subtype layer, frontal+hippo rebuild. Faithful port of the subclass/Hallmark part of
#   manuscript_7fig/scripts/44_gsea_MAST_cliffs.R (family != "broad" -> collection "H").
# -----------------------------------------------------------------------------
# METRIC = signed avg_log2FC winsorized at 1/99th pct within each comparison's
# ranked list (FC-only, p-value-free -> pseudoreplication-safe). Preranking terms
# come from the full-expressed caches written by 47b_full_expressed_subclass_FH.R
# (the atlas is never loaded here).
#
# FAMILIES routed to Hallmark ("H"):
#   subclass       — cortical Azimuth subclass x Frontal (MAST_subclass_discovery_all.csv)
#   hippo_subclass — hippo de-novo subtype x Hippo        (hippo_subtype_meta.csv)
#
# fgsea: gseaParam=1, nPermSimple=1000, set.seed(42), minSize=10, maxSize=500.
# Multiple testing: BH within each comparison (one fgsea call per comparison).
# Power TIERS from n_nuclei_min=min(n_NHD,n_CON): excluded<50 / low_conf 50-149 / solid>=150.
#
# OUTPUT (does not clobber the broad fGSEA_NHD_MAST_cliffs_all.csv):
#   diagnostics/05_GSEA/fGSEA_NHD_MAST_cliffs/fGSEA_NHD_MAST_cliffs_subclass_all.csv
#   diagnostics/05_GSEA/fGSEA_NHD_MAST_cliffs/fGSEA_NHD_MAST_cliffs_subclass_summary.csv
#   ...and a merged fGSEA_NHD_MAST_cliffs_combined_all.csv (broad + subclass), if
#      the broad file is present (drop-in for downstream panel scripts).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(fgsea); library(msigdbr)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot("Could not resolve PROJ root" = !is.na(PROJ), dir.exists(PROJ))
DIR      <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
MAST_DIR <- file.path(DIR, "tables/mast_dual")
OUT_DIR  <- file.path(DIR, "diagnostics/05_GSEA/fGSEA_NHD_MAST_cliffs"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(DIR, "scripts/_gsea_composite_metric.R"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

PCT_FLOOR <- 0.10; HIDE_N <- 50; FADE_N <- 150
POWER_TIERS <- c("excluded","low_conf","solid"); FDR_THRESH <- 0.05

# ============================================================================
# STEP 1: per-comparison metadata for both families
# ============================================================================
msg("=== STEP 1: subclass + hippo-subtype metadata %s ===", format(Sys.time(), "%H:%M:%S"))
subclass_path <- file.path(MAST_DIR, "MAST_subclass_discovery_all.csv")
stopifnot("MISSING MAST_subclass_discovery_all.csv — run 45_mast_subclass_FH.R" = file.exists(subclass_path))
mast_sub <- read_csv(subclass_path, show_col_types = FALSE) %>% mutate(family = "subclass")
cmp_meta <- mast_sub %>% group_by(comparison) %>%
  summarise(cell_type = first(cell_type), region = first(region), family = first(family),
            n_NHD = first(n_NHD), n_CON = first(n_CON), .groups = "drop") %>%
  mutate(n_nuclei_min = pmin(n_NHD, n_CON),
         power_tier = cut(n_nuclei_min, breaks = c(-Inf, HIDE_N, FADE_N, Inf),
                          labels = POWER_TIERS, right = FALSE, ordered_result = TRUE),
         low_n = n_nuclei_min < FADE_N)
msg("  subclass comparisons: %d", nrow(cmp_meta))

# Discovery (DGE) genes per subclass comparison — display-only membership flag
# (not an FDR pre-filter; fgsea runs over the full Hallmark family).
DISC_GENES <- mast_sub %>% filter(discovery %in% c(TRUE, "TRUE")) %>%
  group_by(comparison) %>% summarise(g = list(unique(gene)), .groups = "drop")
disc_genes_for <- function(cmp) { i <- match(cmp, DISC_GENES$comparison); if (is.na(i)) character(0) else DISC_GENES$g[[i]] }

# append hippo de-novo subtype comparisons (family="hippo_subclass"; no MAST table,
# so dge_member stays FALSE). Only those whose prerank cache actually exists.
hippo_meta_path <- file.path(MAST_DIR, "hippo_subtype_meta.csv")
if (file.exists(hippo_meta_path)) {
  hm <- read_csv(hippo_meta_path, show_col_types = FALSE) %>%
    mutate(n_nuclei_min = pmin(n_NHD, n_CON),
           power_tier = cut(n_nuclei_min, breaks = c(-Inf, HIDE_N, FADE_N, Inf),
                            labels = POWER_TIERS, right = FALSE, ordered_result = TRUE),
           low_n = n_nuclei_min < FADE_N) %>%
    select(comparison, cell_type, region, family, n_NHD, n_CON, n_nuclei_min, power_tier, low_n)
  has_cache <- file.exists(file.path(MAST_DIR, sprintf("cliffs_delta_full_%s.csv", hm$comparison)))
  hm <- hm[has_cache, , drop = FALSE]
  if (nrow(hm)) {
    cmp_meta <- bind_rows(cmp_meta, hm)
    msg("  appended %d hippo de-novo subtype comparison(s): %s", nrow(hm), paste(hm$comparison, collapse = ", "))
    print(as.data.frame(hm %>% select(comparison, n_CON, n_NHD, n_nuclei_min, power_tier)), row.names = FALSE)
  }
} else {
  msg("  (no hippo_subtype_meta.csv — hippo de-novo subtypes not included; run 47b)")
}
comparisons <- cmp_meta$comparison

# ============================================================================
# STEP 2: full-expressed caches (47b)
# ============================================================================
msg("\n=== STEP 2: load full-expressed caches (47b) %s ===", format(Sys.time(), "%H:%M:%S"))
cache_path <- function(cmp) file.path(MAST_DIR, sprintf("cliffs_delta_full_%s.csv", cmp))
missing <- comparisons[!file.exists(vapply(comparisons, cache_path, character(1)))]
if (length(missing))
  stop(sprintf("MISSING %d full-expressed cache(s) — run 47b_full_expressed_subclass_FH.R first: %s",
               length(missing), paste(missing, collapse = ", ")))
full_delta <- bind_rows(lapply(comparisons, function(cmp)
  read_csv(cache_path(cmp), show_col_types = FALSE) %>% mutate(comparison = cmp)))
need_cols <- c("gene","cliffs_delta_full","avg_log2FC","pct.1","pct.2")
miss_cols <- setdiff(need_cols, colnames(full_delta))
if (length(miss_cols))
  stop(sprintf("caches missing column(s): %s — re-run 47b", paste(miss_cols, collapse = ", ")))
msg("  loaded %d cache rows across %d comparisons; per-cmp genes %d..%d",
    nrow(full_delta), length(comparisons), min(table(full_delta$comparison)), max(table(full_delta$comparison)))

# ============================================================================
# STEP 3: gene sets — Hallmark for both families
# ============================================================================
msg("\n=== STEP 3: Hallmark gene sets %s ===", format(Sys.time(), "%H:%M:%S"))
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H") %>% {split(.$gene_symbol, .$gs_name)}
msg("  Hallmark=%d", length(hallmark_sets))
clean_pathway <- function(x, wrap = 55) { out <- sub("^(GOBP_|REACTOME_|HALLMARK_)","",x); out <- gsub("_"," ",out); str_wrap(str_to_sentence(out), width = wrap) }

# ============================================================================
# STEP 4: fgsea per comparison (BH within comparison)
# ============================================================================
msg("\n=== STEP 4: fgsea preranked on %s (signed winsorized avg_log2FC) %s ===",
    GSEA_METRIC_NAME, format(Sys.time(), "%H:%M:%S"))
run_one_gsea <- function(cmp) {
  info <- cmp_meta[cmp_meta$comparison == cmp, ]
  tbl  <- full_delta %>% filter(comparison == cmp) %>%
    filter(pmax(pct.1, pct.2) >= PCT_FLOOR, !is.na(avg_log2FC)) %>% distinct(gene, .keep_all = TRUE)
  if (nrow(tbl) < 50) { msg("  %-22s SKIP (only %d ranked genes)", cmp, nrow(tbl)); return(NULL) }
  comp  <- fc_rank(tbl$avg_log2FC)
  ranks <- sort(setNames(comp, tbl$gene), decreasing = TRUE); ranks <- ranks[!is.na(ranks)]
  # both families -> Hallmark (family != "broad")
  sets  <- hallmark_sets
  dge_g <- disc_genes_for(cmp)
  set.seed(42)
  fg <- tryCatch(fgsea(pathways = sets, stats = ranks, minSize = 10, maxSize = 500, gseaParam = 1, nPermSimple = 1000),
                 error = function(e) { msg("  fgsea err %s: %s", cmp, conditionMessage(e)); NULL })
  if (is.null(fg) || !nrow(fg)) return(NULL)
  has_dge <- if (length(dge_g)) vapply(sets[fg$pathway], function(g) any(g %in% dge_g), logical(1)) else rep(FALSE, nrow(fg))
  fg <- as_tibble(fg) %>% mutate(
      subtype = cmp, direction = if_else(NES > 0, "NHD", "CON"),
      db = "Hallmark", pathway_clean = clean_pathway(pathway), metric = GSEA_METRIC_NAME,
      family = info$family,
      n_nuclei_min = info$n_nuclei_min, power_tier = as.character(info$power_tier), low_n = info$low_n,
      dge_member = has_dge, leadingEdge = vapply(leadingEdge, paste, collapse = "|", FUN.VALUE = character(1))) %>%
    arrange(padj, desc(abs(NES)))
  msg("  %-22s [%-14s] %4d tested | sig(FDR<%.2g)=%d | DGE-member=%d | ranked=%d%s",
      cmp, info$family, nrow(fg), FDR_THRESH, sum(fg$padj < FDR_THRESH, na.rm = TRUE),
      sum(fg$dge_member, na.rm = TRUE), length(ranks),
      if (isTRUE(info$low_n)) sprintf("  [LOW_N n=%d]", info$n_nuclei_min) else "")
  fg %>% select(subtype, pathway, pval, padj, ES, NES, size, leadingEdge, direction, db,
                pathway_clean, metric, family, n_nuclei_min, power_tier, low_n, dge_member)
}
gsea_all <- bind_rows(lapply(comparisons, run_one_gsea))
sub_csv <- file.path(OUT_DIR, "fGSEA_NHD_MAST_cliffs_subclass_all.csv")
write_csv(gsea_all, sub_csv)
msg("Saved: %s (%d rows)", sub_csv, nrow(gsea_all))

# ============================================================================
# STEP 5: merge with the existing broad GSEA (non-destructive: broad file kept)
# ============================================================================
broad_csv <- file.path(OUT_DIR, "fGSEA_NHD_MAST_cliffs_all.csv")
if (file.exists(broad_csv)) {
  broad <- read_csv(broad_csv, show_col_types = FALSE)
  if (!"family" %in% colnames(broad)) broad$family <- "broad"
  common <- intersect(colnames(broad), colnames(gsea_all))
  combined <- bind_rows(broad[, common, drop = FALSE], gsea_all[, common, drop = FALSE])
  comb_csv <- file.path(OUT_DIR, "fGSEA_NHD_MAST_cliffs_combined_all.csv")
  write_csv(combined, comb_csv)
  msg("Saved: %s (%d broad + %d subclass/subtype rows)", comb_csv, nrow(broad), nrow(gsea_all))
} else {
  msg("  (broad fGSEA_NHD_MAST_cliffs_all.csv not present — combined file skipped)")
}

# ============================================================================
# STEP 6: summary
# ============================================================================
msg("\n=== STEP 6: per-comparison Hallmark summary ===")
if (nrow(gsea_all)) {
  s2 <- gsea_all %>% group_by(subtype) %>%
    summarise(family = first(family), pathways = dplyr::n(),
              sig_FDR005 = sum(padj < FDR_THRESH, na.rm = TRUE),
              sig_up_NHD = sum(padj < FDR_THRESH & NES > 0, na.rm = TRUE),
              sig_up_CON = sum(padj < FDR_THRESH & NES < 0, na.rm = TRUE),
              n_nuclei_min = first(n_nuclei_min), power_tier = first(power_tier), low_n = first(low_n),
              .groups = "drop")
  print(as.data.frame(s2), row.names = FALSE)
  write_csv(s2, file.path(OUT_DIR, "fGSEA_NHD_MAST_cliffs_subclass_summary.csv"))
}
cat("\n=== DONE ===\n", file = stderr()); print(sessionInfo())
