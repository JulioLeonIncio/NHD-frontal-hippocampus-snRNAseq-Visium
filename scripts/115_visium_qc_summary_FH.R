#!/usr/bin/env Rscript
# =============================================================================
# 115_visium_qc_summary_FH.R — Inputs : Visium/cell2location/c2l_MAIN/visium_meta.csv (the integrated object's metadata, script 01f)
#          Visium/integrated_harmony/run_L1_20260918.log      (pre-floor spot counts, "<sec>: N -> M spots")
#          <manifest outs>/tissue_call.csv per section         (tier per barcode; folder via _visium_outs_FH.R)
#          Visium/integrated_harmony/integrated_cluster_identity.csv  (cluster -> domain, via _visium_domains_FH.R)
# Output : tables/visium_qc_summary_FH.csv           (long: level, unit, metric, value)
#          tables/visium_domain_occupancy_depth_FH.csv (per cluster + per domain: n, n_NHD, frac_NHD, median UMI)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
ROOT <- dirname(PROJ)
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))   # DOM_MAP (cluster -> domain), DOM_LEV
source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))      # VIS_SECTIONS + vis_outs(sec) from VISIUM_OUTS_MANIFEST.json
say <- function(...) cat(sprintf(...), "\n")

m <- fread(file.path(ROOT, "Visium", "cell2location", "c2l_MAIN", "visium_meta.csv")); setnames(m, 1, "spot_id")
stopifnot(all(c("sample_id","condition","seurat_clusters","percent.mt","nCount_Spatial","nFeature_Spatial") %in% names(m)))
m[, domain := unname(DOM_MAP[as.character(seurat_clusters)])]; stopifnot(!anyNA(m$domain))

# ---- spots before / after the floor: BEFORE = the in-tissue barcodes of the final tissue call (tissue_call.csv of
#      the manifest outs, script 117b); AFTER = the spots in the object. The integration log is no longer parsed: a RESUME=1 run
#      reads the checkpoint and never prints the "X -> Y spots" lines.
TIERS <- c("spaceranger", "he_mask", "expression")
tier <- rbindlist(lapply(VIS_SECTIONS, function(s) {
  f <- file.path(vis_outs(s), "tissue_call.csv"); stopifnot("MISSING tissue_call.csv (117b)" = file.exists(f))
  tc <- fread(f); stopifnot(all(c("barcode","in_tissue_analysis","tissue_call_source") %in% names(tc)))
  tc <- tc[in_tissue_analysis == 1]
  stopifnot("unknown tissue_call_source among in-tissue barcodes" = all(tc$tissue_call_source %in% TIERS))
  data.table(sample_id = s, metric = paste0("n_", TIERS), value = as.numeric(vapply(TIERS, function(t) sum(tc$tissue_call_source == t), numeric(1)))) }))
tier_tot <- tier[, .(n_before = as.integer(sum(value))), by = sample_id]
floor <- merge(tier_tot, m[, .(n_after = .N), by = sample_id], by = "sample_id")
stopifnot("object has spots that are not in the final tissue call" = all(floor$n_after <= floor$n_before),
          "object carries a section missing from the manifest" = nrow(floor) == length(VIS_SECTIONS))
floor[, n_removed := n_before - n_after]
say("floor nCount >= 50 & nFeature >= 25: %d -> %d spots (%d removed)", sum(floor$n_before), sum(floor$n_after), sum(floor$n_removed))

print(dcast(tier, sample_id ~ metric, value.var = "value"))

# ---- per-section / per-condition medians
med <- function(d, by) d[, .(n_spots = as.numeric(.N), median_umi = as.numeric(median(nCount_Spatial)), median_genes = as.numeric(median(nFeature_Spatial)),
                             median_pct_mt = median(percent.mt)), by = by]
sec <- med(m, "sample_id"); con <- med(m, "condition")
qc <- rbind(
  floor[, .(level = "section", unit = sample_id, metric = "n_spots_before_floor", value = n_before)],
  floor[, .(level = "section", unit = sample_id, metric = "n_spots_removed_by_floor", value = n_removed)],
  tier[, .(level = "section", unit = sample_id, metric, value)],   # n_spaceranger / n_he_mask / n_expression per section
  tier[, .(value = sum(value)), by = metric][, .(level = "all", unit = "all", metric, value)],
  data.table(level = "all", unit = "all", metric = c("n_spots_before_floor", "n_spots_removed_by_floor", "n_spots"),
             value = c(sum(floor$n_before), sum(floor$n_removed), nrow(m))),
  melt(sec, id.vars = "sample_id", variable.name = "metric", value.name = "value")[, .(level = "section", unit = sample_id, metric, value)],
  melt(con, id.vars = "condition", variable.name = "metric", value.name = "value")[, .(level = "condition", unit = condition, metric, value)])

# ---- occupancy vs depth, per cluster and per domain
occ_cl  <- m[, .(level = "cluster", n = .N, n_NHD = sum(condition == "NHD"), median_umi = as.numeric(median(nCount_Spatial))), by = .(unit = as.character(seurat_clusters))]
occ_dom <- m[, .(level = "domain",  n = .N, n_NHD = sum(condition == "NHD"), median_umi = as.numeric(median(nCount_Spatial))), by = .(unit = domain)]
occ <- rbind(occ_cl, occ_dom)[, frac_NHD := n_NHD / n][]
occ[, unit_order := ifelse(level == "domain", match(unit, DOM_LEV), suppressWarnings(as.integer(unit)))]; setorder(occ, level, unit_order); occ[, unit_order := NULL]
rho <- cor(occ_dom$n_NHD / occ_dom$n, occ_dom$median_umi, method = "spearman")
qc <- rbind(qc, data.table(level = c("cluster","cluster","domain","domain","domain"), unit = "all",
                           metric = c("frac_NHD_min","frac_NHD_max","frac_NHD_min","frac_NHD_max","spearman_fracNHD_vs_medianUMI"),
                           value = c(range(occ_cl$n_NHD / occ_cl$n), range(occ_dom$n_NHD / occ_dom$n), rho)))
say("NHD occupancy: clusters %.1f-%.1f %%, domains %.0f-%.0f %%; Spearman rho(frac NHD, median UMI) over %d domains = %.2f",
    100 * min(occ_cl$n_NHD / occ_cl$n), 100 * max(occ_cl$n_NHD / occ_cl$n), 100 * min(occ_dom$n_NHD / occ_dom$n), 100 * max(occ_dom$n_NHD / occ_dom$n), nrow(occ_dom), rho)
print(con); print(occ[level == "domain"])
fwrite(qc,  file.path(PROJ, "tables", "visium_qc_summary_FH.csv"))
fwrite(occ, file.path(PROJ, "tables", "visium_domain_occupancy_depth_FH.csv"))
writeLines(capture.output(sessionInfo()), file.path(PROJ, "logs", "115_visium_qc_summary_FH_sessionInfo.txt"))
cat("\n=== DONE ===\n", file = stderr())
