#!/usr/bin/env Rscript
# =============================================================================
# 111_visium_proximity_severity_FH.R — Does the microglial state and its downstream consequences co-localize in tissue?
# -----------------------------------------------------------------------------
# D1  microglia-rich vs microglia-poor SPOTS, within each section and compartment (grey matter
#     = the domains flagged compartment == "grey" in integrated_domain_levels.csv; white matter = "white"). PRIMARY ranking = cell2location microglial
#     FRACTION (microglial q05 / total q05 over the eight reference types), quartiles defined
#     within section x domain, then pooled — because absolute microglial abundance in control
#     tissue is a cellularity index (it co-varies with neuron/oligodendrocyte content), so an
#     absolute-abundance contrast in controls compares cellular with acellular spots, not
#     microglia-rich with microglia-poor. Absolute abundance ships as the sensitivity version.
#     Outcomes per spot (3,000-UMI-thinned AddModuleScore, script 65; depth is out):
#     astrocyte reactive, complement/MHC-II, antigen presentation, sterol arm, structural
#     myelin, fatty-acid/sphingomyelin, lipid uptake, glycolytic shift, FTL/HAMP (iron; the
#     module is carried by FTL and HAMP alone), SPP1 (the DAM-2 module is carried by SPP1
#     alone), housekeeping null, homeostatic microglia (P2RY12/CX3CR1/CSF1R/TMEM119/SALL1 — the
#     "more microglia vs a different microglial state" reference: it rises with the microglial
#     fraction if the rich spots simply hold more microglia); plus the abundance confound rows
#     (astrocyte, oligodendrocyte, excitatory-neuron q05) and the oligodendrocyte-decile-stratified
#     delta (table only).
#     statistic: Cliff's delta with a block-bootstrap interval (tissue blocks of 6 array rows x
#     12 array columns within section — a within-donor sampling interval, never a donor-level
#     CI) and a toroidal-shift NULL: the ranking map is wrapped by 200 random parity-preserving
#     shifts of the array coordinates within the section and delta recomputed each time; the
#     2.5-97.5 % envelope is drawn behind every row and each delta is reported as a percentile
#     of it. No P value anywhere: condition is perfectly aliased with donor.
#     neighborhood version: mean fraction of the six hex neighbors (>= 5 present), excluding
#     the spot, so "near microglia-rich tissue" is tested at one-spot distance; drawn as the open
#     symbol beside the spot-ranked delta in panel a.
#     PRE-SPECIFIED primary contrast: NHD grey matter, spot version, four predicted outcomes
#     (astrocyte reactive up, complement/MHC-II up, sterol arm down, structural myelin down near
#     microglia-rich spots). Everything else is exploratory and labelled so in the table.
#     circularity, stated: the cell2location microglial signature comes from our own atlas and
#     shares genes with the complement/MHC-II and antigen-presentation modules. The reference
#     (ref_signatures.csv, genes x cell types) is ranked by microglial SPECIFICITY (Micro-PVM
#     signature / row sum over the eight types) and, for every outcome module, the number of
#     members inside the top-200 specific microglial genes is written to
#     tables/visium_proximity_signature_overlap_FH.csv; the module membership comes from
#     tables/visium_dm_module_genes_FH.csv (script 65), so the counts describe the sets as scored.
#     Those two rows are a positive control that the state is spatially concentrated, not
#     independent evidence of activation.
#     white matter: the WM domain is carried almost entirely by the NHD sections (spot counts
#     per section x domain are read from the live CSVs, printed below and shipped in
#     tables/visium_proximity_FH.csv); the WM facet therefore shows NHD sections 1 and 2 only,
#     without the toroidal-shift band (a curved 85-spot band wrapped on the array is not a usable
#     reference; the shift null stays in the table), with the housekeeping-null drift as its floor.
#
# D2(ii) SEVERITY (molecular) is computed for the table only and not plotted: the
#     structural-myelin score tracks grey/white composition, not lesion severity (it inverts the
#     microglial programs in deep layers in both donors). Histological strata (D2(i)) wait for
#     the neuropathologists' outlines.
#
# Inputs (never recomputed): tables/visium_dm_module_spots_FH.csv (65 v5),
#   Data_deposition/Visium/spot_metadata.csv.gz, tissue_positions_list.csv, c2l ref_signatures.
# Outputs: tables/visium_proximity_FH.csv (all rows, all versions, null envelopes, percentiles)
#          tables/visium_severity_FH.csv, tables/visium_proximity_signature_overlap_FH.csv
#          figures/Supplementary/SuppFig6_microglia_proximity/S6a_proximity_dots.{png,pdf}
#          figures/Figure_6/panels/F6f_microglia_proximity_dots.{png,pdf}   (the same rows at the half-row Fig-6 size; panel g = script 114)
#          (panels b and c of Supplementary Fig. 6 — composition per domain — come from script 114)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggh4x); library(patchwork); library(ragg) })
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))   # DOM_COMP: domain -> grey | white | other
source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))      # vis_outs(sec) — outs folder per section from VISIUM_OUTS_MANIFEST.json
ROOT <- dirname(PROJ); TDIR <- file.path(PROJ, "tables")
OUT  <- file.path(PROJ, "figures", "Supplementary", "SuppFig6_microglia_proximity"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
NBOOT <- 2000L; NSHIFT <- 200L; BLOCK <- c(rows = 6L, cols = 12L)

# ---- inputs -----------------------------------------------------------------
sp <- fread(file.path(TDIR, "visium_dm_module_spots_FH.csv"))
md <- fread(file.path(PROJ, "manuscript", "Final_fig_and_tables", "Data_deposition", "Visium", "spot_metadata.csv.gz"))
CT8 <- c("Astro","Endo","Micro-PVM","Neuron_Ex","Neuron_Inh","OPC","Oligo","Pericytes")
stopifnot(all(paste0("c2l_q05_", CT8) %in% names(md)))
md[, total_q05 := rowSums(.SD), .SDcols = paste0("c2l_q05_", CT8)]
mods <- c("Astrocyte reactive","Complement / MHC-II","Antigen presentation","Cholesterol (sterol arm)","Structural myelin",
          "Fatty-acid / sphingomyelin","Lipid uptake / salvage","Glycolytic shift","Iron handling","TREM2-dep lipid (DAM-2)","Housekeeping (null)",
          "Microglia homeostatic")
stopifnot(all(mods %in% names(sp)))
dall <- merge(sp, md[, .(spot_id, spatial_domain, micro = `c2l_q05_Micro-PVM`, astro = c2l_q05_Astro, oligo = c2l_q05_Oligo, neu = c2l_q05_Neuron_Ex, total_q05)], by = "spot_id")
stopifnot(nrow(dall) == nrow(sp))
dall[, micro_frac := fifelse(total_q05 > 0, micro / total_q05, NA_real_)]   # a spot with no deconvolved content has no fraction (0/0), not a fraction of 0
say("spots with zero total q05 (fraction undefined, excluded from ranking): %d", sum(is.na(dall$micro_frac)))
stopifnot("spot_metadata carries a domain absent from integrated_domain_levels.csv" = all(unique(dall$spatial_domain) %in% names(DOM_COMP)))
dall[, compartment := fifelse(DOM_COMP[spatial_domain] == "grey", "Grey matter", fifelse(DOM_COMP[spatial_domain] == "white", "White matter", NA_character_))]
say("grey-matter domains: %s | white-matter domain(s): %s", paste(DOM_GREY, collapse = ", "), paste(DOM_WHITE, collapse = ", "))
setnames(dall, c("Iron handling", "TREM2-dep lipid (DAM-2)", "Microglia homeostatic"), c("FTL/HAMP (iron)", "SPP1 (DAM-2 arm)", "Homeostatic microglia"))   # single-gene-carried modules named for their carriers
mods <- sub("^Microglia homeostatic$", "Homeostatic microglia", sub("^Iron handling$", "FTL/HAMP (iron)", sub("^TREM2-dep lipid \\(DAM-2\\)$", "SPP1 (DAM-2 arm)", mods)))
d <- dall[!is.na(compartment)]   # the statistic runs on grey + white matter; dall keeps every spot for the maps
say("spots in grey/white compartments: %d", nrow(d)); print(d[, .N, by = .(sample_id, compartment)][order(sample_id, compartment)])

# ---- signature overlap (circularity, stated) -------------------------------------
# The cell2location reference is ranked by microglial SPECIFICITY (Micro-PVM signature divided by the
# row sum over the eight reference types; genes with a zero row sum carry no signature and are
# excluded) — raw magnitude ranks housekeeping-level genes first and says nothing about circularity.
# Every outcome module scored by script 65 is counted against the top-200 specific genes; the
# magnitude-ranked count is kept as a second column for reference.
sigf <- file.path(ROOT, "Visium", "cell2location", "c2l_MAIN", "ref_signatures.csv")
stopifnot("MISSING cell2location ref_signatures.csv" = file.exists(sigf))
modgenes <- fread(file.path(TDIR, "visium_dm_module_genes_FH.csv"))[present == TRUE]
modgenes[, module := sub("^Microglia homeostatic$", "Homeostatic microglia", sub("^Iron handling$", "FTL/HAMP (iron)", sub("^TREM2-dep lipid \\(DAM-2\\)$", "SPP1 (DAM-2 arm)", module)))]
stopifnot(all(mods %in% modgenes$module))
sg <- fread(sigf); gcol <- names(sg)[1]; setnames(sg, gcol, "gene"); stopifnot(all(CT8 %in% names(sg)))
sg[, rowsum := rowSums(.SD), .SDcols = CT8]
say("reference genes: %d; zero row sum (no signature, dropped from the specificity ranking): %d", nrow(sg), sum(sg$rowsum == 0))
sg <- sg[rowsum > 0][, specificity := `Micro-PVM` / rowsum]
NTOP <- 200L
top_spec <- sg[order(-specificity, -`Micro-PVM`)][1:NTOP]$gene     # ties broken by magnitude, deterministic
top_magn <- sg[order(-`Micro-PVM`, -specificity)][1:NTOP]$gene
ov <- modgenes[module %in% mods, .(n_genes = .N, genes = paste(gene, collapse = ","),
        n_in_top200_specific = sum(gene %in% top_spec), which_specific = paste(gene[gene %in% top_spec], collapse = ","),
        n_in_top200_magnitude = sum(gene %in% top_magn), which_magnitude = paste(gene[gene %in% top_magn], collapse = ","),
        median_specificity = median(sg$specificity[match(gene, sg$gene)], na.rm = TRUE)), by = module]
ov <- ov[order(-n_in_top200_specific, module)]
ov[, specificity_threshold_top200 := min(sg$specificity[sg$gene %in% top_spec])]
fwrite(ov, file.path(TDIR, "visium_proximity_signature_overlap_FH.csv"))
say("\n== module genes inside the top-%d microglia-SPECIFIC cell2location reference genes (threshold specificity >= %.3f) ==", NTOP, ov$specificity_threshold_top200[1])
print(ov[, .(module, n_genes, n_in_top200_specific, which_specific, n_in_top200_magnitude, median_specificity = round(median_specificity, 3))])

# ---- hex neighbors from array coordinates (>= 5 of 6 present) --------------------
pos <- rbindlist(lapply(unique(d$sample_id), function(s) {
  p <- fread(file.path(vis_outs(s), "spatial", "tissue_positions_list.csv"), header = FALSE,   # manifest outs
             col.names = c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"))
  p[, `:=`(sample_id = s, spot_id = paste0(s, "_", barcode))]; p }))
d <- merge(d, pos[, .(spot_id, array_row, array_col)], by = "spot_id"); stopifnot(!anyNA(d$array_row))
OFF <- list(c(0, 2), c(0, -2), c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
nb_mean <- function(dd, v) { key <- paste(dd$array_row, dd$array_col); idx <- setNames(seq_len(nrow(dd)), key)
  vapply(seq_len(nrow(dd)), function(i) { j <- idx[vapply(OFF, function(o) paste(dd$array_row[i] + o[1], dd$array_col[i] + o[2]), character(1))]
    j <- j[!is.na(j)]; if (length(j) < 5) NA_real_ else mean(dd[[v]][j]) }, numeric(1)) }
d[, frac_nb := nb_mean(.SD, "micro_frac"), by = sample_id, .SDcols = c("array_row","array_col","micro_frac")]
d[, micro_nb := nb_mean(.SD, "micro"),      by = sample_id, .SDcols = c("array_row","array_col","micro")]
d[, block := paste(sample_id, array_row %/% BLOCK[["rows"]], array_col %/% BLOCK[["cols"]])]

# ---- statistics ------------------------------------------------------------------
cliff <- function(a, b) { if (length(a) < 5 || length(b) < 5) return(NA_real_)
  r <- rank(c(a, b)); (2 * (sum(r[seq_along(a)]) - length(a) * (length(a) + 1) / 2) / (length(a) * length(b))) - 1 }
# Block bootstrap for all outcomes at once: the block resample is drawn once per replicate and
# Cliff's delta is read off the ranks of a numeric matrix (no per-outcome data.table copies)
cliff_mat <- function(M, ia) { na <- sum(ia); nb <- length(ia) - na; if (na < 5 || nb < 5) return(rep(NA_real_, ncol(M)))
  vapply(seq_len(ncol(M)), function(o) { r <- data.table::frank(M[, o], ties.method = "average"); (2 * (sum(r[ia]) - na * (na + 1) / 2) / (na * nb)) - 1 }, numeric(1)) }
block_ci_all <- function(M, ia, block, n = NBOOT) {   # rows of M = proximal + distal spots; returns 2 x ncol(M)
  bl <- split(seq_len(nrow(M)), block); if (length(bl) < 8) return(matrix(NA_real_, 2, ncol(M)))
  v <- matrix(NA_real_, n, ncol(M))
  for (r in seq_len(n)) { i <- unlist(bl[sample.int(length(bl), replace = TRUE)], use.names = FALSE); v[r, ] <- cliff_mat(M[i, , drop = FALSE], ia[i]) }
  apply(v, 2, quantile, probs = c(0.025, 0.975), names = FALSE, na.rm = TRUE) }
cliff_strat <- function(a, b, sa, sb) { br <- unique(quantile(c(sa, sb), seq(0, 1, 0.1))); if (length(br) < 3) return(NA_real_)
  ga <- cut(sa, br, include.lowest = TRUE); gb <- cut(sb, br, include.lowest = TRUE); w <- 0; acc <- 0
  for (lv in levels(ga)) { x <- a[ga == lv]; y <- b[gb == lv]; if (length(x) >= 5 && length(y) >= 5) { acc <- acc + cliff(x, y) * (length(x) + length(y)); w <- w + length(x) + length(y) } }
  if (w > 0) acc / w else NA_real_ }
# toroidal shift null: shift the RANKING map on the array (parity-preserving: even row and column
# offsets), wrap within the section's array extent, re-form quartiles, recompute delta
shift_null <- function(dd, rank_col, outcomes, nshift = NSHIFT) {
  out <- matrix(NA_real_, nshift, length(outcomes), dimnames = list(NULL, outcomes))

  rr <- range(dd$array_row); cc <- range(dd$array_col); nr <- diff(rr) + 1; nc <- diff(cc) + 1
  key <- paste(dd$sample_id, dd$array_row, dd$array_col); idx <- setNames(seq_len(nrow(dd)), key)
  for (k in seq_len(nshift)) {
    dr <- 2 * sample.int(nr %/% 2, 1); dc <- 2 * sample.int(nc %/% 2, 1)
    r2 <- rr[1] + (dd$array_row - rr[1] + dr) %% nr; c2 <- cc[1] + (dd$array_col - cc[1] + dc) %% nc
    j <- idx[paste(dd$sample_id, r2, c2)]                                   # the spot now sitting under the shifted map
    ok <- !is.na(j); if (sum(ok) < 40) next
    v <- rep(NA_real_, nrow(dd)); v[ok] <- dd[[rank_col]][j[ok]]              # shifted ranking variable
    ok2 <- which(ok); q <- rep(NA_real_, nrow(dd))
    q[ok2] <- ave(v[ok2], paste(dd$sample_id[ok2], dd$spatial_domain[ok2]), FUN = function(z) (data.table::frank(z) - 0.5) / length(z))
    sel <- which(!is.na(q) & (q >= 0.75 | q <= 0.25)); if (length(sel) < 10) next
    out[k, ] <- cliff_mat(as.matrix(dd[sel, outcomes, with = FALSE]), q[sel] >= 0.75)
  }
  out }

OUTCOMES <- c(mods, "astro", "oligo", "neu")
OUT_LAB  <- c(setNames(mods, mods), astro = "Astrocyte abundance (c2l)", oligo = "Oligodendrocyte abundance (c2l)", neu = "Excitatory-neuron abundance (c2l)")
PRIMARY  <- c("Astrocyte reactive", "Complement / MHC-II", "Cholesterol (sterol arm)", "Structural myelin")
# panel row order: the four primary outcomes with the homeostatic reference directly under complement / MHC-II,
# then the exploratory modules in table order, then the three abundance confounds
PANEL_ORDER <- c("Astrocyte reactive", "Complement / MHC-II", "Homeostatic microglia", "Cholesterol (sterol arm)", "Structural myelin",
                 setdiff(mods, c(PRIMARY, "Homeostatic microglia")), "astro", "oligo", "neu")
stopifnot(setequal(PANEL_ORDER, OUTCOMES))
sets <- list(NHD_Frontal1 = "NHD_Frontal1", NHD_Frontal2 = "NHD_Frontal2", `NHD (both)` = c("NHD_Frontal1","NHD_Frontal2"),
             CON_Frontal1 = "CON_Frontal1", CON_Frontal2 = "CON_Frontal2", `CON (both)` = c("CON_Frontal1","CON_Frontal2"))
VERSIONS <- c(fraction = "micro_frac", abundance = "micro", fraction_neighborhood = "frac_nb", abundance_neighborhood = "micro_nb")

# the bootstrap + shift-null loop takes ~40 min; its result is cached against the two input tables'
# mtimes so the panel section can be re-rendered without recomputing (delete the cache to force)
CACHE <- file.path(PROJ, "data", "_cache_111_proximity_FH.rds")
KEY <- paste(file.info(file.path(TDIR, "visium_dm_module_spots_FH.csv"))$mtime, file.info(file.path(PROJ, "manuscript", "Final_fig_and_tables", "Data_deposition", "Visium", "spot_metadata.csv.gz"))$mtime,
             NBOOT, NSHIFT, paste(BLOCK, collapse = ","), paste(names(VERSIONS), VERSIONS, collapse = ","), paste(names(sets), sapply(sets, paste, collapse = "+"), collapse = ","),
             paste(OUTCOMES, collapse = ","), paste(PRIMARY, collapse = ","), paste(names(OUT_LAB), OUT_LAB, collapse = ","),
             paste(DOM_GREY, collapse = ","), paste(DOM_WHITE, collapse = ","))   # every constant that changes the result, incl. the compartment definition
# one contrast = one job; the 48 jobs run in parallel (L'Ecuyer streams, seed 42 -> reproducible)
run_contrast <- function(nm, cp, ver) {
  out <- list()

  dd <- d[sample_id %in% sets[[nm]] & compartment == cp & !is.na(get(VERSIONS[[ver]]))]
  if (nrow(dd) < 60) { if (ver == "fraction") say("  skip %s / %s (n = %d)", nm, cp, nrow(dd)); return(NULL) }
  dd[, rv := get(VERSIONS[[ver]])]
  dd[, q := (rank(rv) - 0.5) / .N, by = .(sample_id, spatial_domain)]         # quartiles within section x domain, then pooled
  dd[, grp := fifelse(q >= 0.75, "prox", fifelse(q <= 0.25, "dist", NA_character_))]
  prox <- dd[grp == "prox"]; dist <- dd[grp == "dist"]
  nul <- if (ver %in% c("fraction", "abundance")) shift_null(dd, VERSIONS[[ver]], OUTCOMES) else NULL
  pd_ <- dd[grp %in% c("prox","dist")]; CI <- block_ci_all(as.matrix(pd_[, OUTCOMES, with = FALSE]), pd_$grp == "prox", pd_$block); colnames(CI) <- OUTCOMES
  for (o in OUTCOMES) { a <- prox[[o]]; b <- dist[[o]]; dl <- cliff(a, b)
    ci <- CI[, o]
    nv <- if (!is.null(nul)) nul[, o] else NA
    out[[length(out) + 1]] <- data.table(set = nm, version = ver, compartment = cp, outcome = OUT_LAB[[o]],
      primary = (nm == "NHD (both)" & ver == "fraction" & cp == "Grey matter" & o %in% PRIMARY),
      n_proximal = nrow(prox), n_distal = nrow(dist), n_blocks = length(unique(dd$block)),
      rank_median_proximal = median(prox$rv), rank_median_distal = median(dist$rv),
      delta = dl, ci_lo = ci[1], ci_hi = ci[2],
      null_lo = if (all(is.na(nv))) NA_real_ else quantile(nv, 0.025, na.rm = TRUE, names = FALSE),
      null_hi = if (all(is.na(nv))) NA_real_ else quantile(nv, 0.975, na.rm = TRUE, names = FALSE),
      null_percentile = if (all(is.na(nv))) NA_real_ else mean(nv <= dl, na.rm = TRUE),
      delta_oligo_stratified = if (o == "oligo") NA_real_ else cliff_strat(a, b, prox$oligo, dist$oligo),
      median_proximal = median(a), median_distal = median(b)) }
  say("  %-13s %-12s %-24s n = %4d / %4d  blocks %3d", nm, cp, ver, nrow(prox), nrow(dist), length(unique(dd$block)))
  rbindlist(out) }
jobs <- expand.grid(nm = names(sets), cp = c("Grey matter", "White matter"), ver = names(VERSIONS), stringsAsFactors = FALSE)
res <- list()
if (file.exists(CACHE) && identical(readRDS(CACHE)$key, KEY)) { res <- readRDS(CACHE)$res; say("proximity statistics loaded from cache (%s)", basename(CACHE)) } else {
  RNGkind("L'Ecuyer-CMRG"); set.seed(42)
  res <- parallel::mclapply(seq_len(nrow(jobs)), function(i) run_contrast(jobs$nm[i], jobs$cp[i], jobs$ver[i]),
                            mc.cores = max(1L, min(6L, parallel::detectCores() - 2L)), mc.set.seed = TRUE)
  bad <- vapply(res, inherits, logical(1), "try-error"); if (any(bad)) stop("contrast failed: ", paste(res[bad][[1]]))
}
if (is.list(res) && !is.data.frame(res)) { res <- rbindlist(res); res[, condition := ifelse(grepl("^NHD", set), "NHD", "CON")]
  saveRDS(list(key = KEY, res = res), CACHE) }
fwrite(res, file.path(TDIR, "visium_proximity_FH.csv"))
fwrite(res, file.path(PROJ, "supplementary_tables", "ST20_visium_microglia_proximity.csv"))   # shipped copy -> Supplementary Data 11, sheet Microglia_proximity (script 100)
say("\n== PRIMARY contrast (NHD grey matter, fraction ranking, pooled sections) + controls ==")
print(res[version == "fraction" & set %in% c("NHD (both)","CON (both)") & compartment == "Grey matter",
          .(set, outcome, primary, delta = round(delta, 2), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi),
            null = sprintf("[%.2f, %.2f]", null_lo, null_hi), pct = round(null_percentile, 3), strat = round(delta_oligo_stratified, 2))], nrows = 40)
say("\n== white matter (NHD sections 1 and 2; no control WM domain of >= 60 spots captured) ==")
print(res[version %in% c("fraction", "fraction_neighborhood") & set %in% c("NHD_Frontal1", "NHD_Frontal2") & compartment == "White matter",
          .(set, version, outcome, n_proximal, delta = round(delta, 2), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi), null = sprintf("[%.2f, %.2f]", null_lo, null_hi))], nrows = 60)

# ---- D2(ii) severity — table only --------------------------------------------------------
sev <- list()
MICRO_PROG <- c("Complement / MHC-II","Antigen presentation","FTL/HAMP (iron)","Glycolytic shift","Astrocyte reactive","Housekeeping (null)")
TREM2_ARM  <- c("gene_LPL","gene_CST7","gene_GPNMB","gene_ITGAX","gene_TREM2","gene_SPP1","gene_CD9")
for (nm in c("NHD_Frontal1","NHD_Frontal2","NHD (both)","CON (both)")) for (dom in list(WM = DOM_WHITE, deep = c("L5","L6"))) {
  dd <- d[sample_id %in% sets[[nm]] & spatial_domain %in% dom]; if (nrow(dd) < 60) next
  dd[, sev_rank := (rank(`Structural myelin`) - 0.5) / .N, by = sample_id]
  severe <- dd[sev_rank <= 1/3]; mild <- dd[sev_rank > 2/3]
  for (o in MICRO_PROG) sev[[length(sev) + 1]] <- data.table(set = nm, domain = paste(dom, collapse = "+"), outcome = o, kind = "module",
      n_severe = nrow(severe), n_mild = nrow(mild), delta = cliff(severe[[o]], mild[[o]]), value_severe = median(severe[[o]]), value_mild = median(mild[[o]]))
  for (g in TREM2_ARM) if (g %in% names(dd)) sev[[length(sev) + 1]] <- data.table(set = nm, domain = paste(dom, collapse = "+"), outcome = sub("^gene_", "", g), kind = "detection",
      n_severe = nrow(severe), n_mild = nrow(mild), delta = mean(severe[[g]] > 0) - mean(mild[[g]] > 0), value_severe = mean(severe[[g]] > 0), value_mild = mean(mild[[g]] > 0))
}
fwrite(rbindlist(sev), file.path(TDIR, "visium_severity_FH.csv"))

# ---- panels --------------------------------------------------------------------------------
# panel a: filled symbol = spot-ranked delta with its block-bootstrap interval; open symbol = the
# same contrast with spots ranked by the mean microglial fraction of their six hex neighbors
# (spot excluded; interval in the table). Grey matter carries the toroidal-shift envelope (widest
# across the four sections); the white-matter facet is the two NHD sections and no envelope.
pa <- res[version %in% c("fraction", "fraction_neighborhood") & !grepl("both", set) & (compartment == "Grey matter" | condition == "NHD")]
pa[, facet := factor(fifelse(compartment == "White matter", "White matter (NHD)", "Grey matter"), levels = c("Grey matter", "White matter (NHD)"))]
pa[, section := ifelse(grepl("1$", set), "section 1", "section 2")]
pa[, ranking := fifelse(version == "fraction", "spot", "neighborhood")]
pa[, key := factor(paste(section, ranking, sep = ", "), levels = c("section 1, spot", "section 1, neighborhood", "section 2, spot", "section 2, neighborhood"))]
pa[, outcome_f := factor(outcome, levels = rev(unname(OUT_LAB[PANEL_ORDER])))]
pa[, dodge := interaction(condition, section)]                 # spot and neighborhood symbols of one section share a dodge slot
say("\n== panel a rows per facet x section x ranking (n_proximal = n_distal) ==")
print(dcast(unique(pa[, .(facet, condition, section, ranking, n_proximal)]), facet + condition + section ~ ranking, value.var = "n_proximal"))
nulls <- pa[facet == "Grey matter" & ranking == "spot", .(lo = min(null_lo, na.rm = TRUE), hi = max(null_hi, na.rm = TRUE)), by = .(facet, outcome_f)]
SHAPES <- c(`section 1, spot` = 21, `section 1, neighborhood` = 1, `section 2, spot` = 24, `section 2, neighborhood` = 2)
p1 <- ggplot(pa, aes(y = outcome_f, group = dodge)) +
  geom_rect(data = nulls, aes(xmin = lo, xmax = hi, ymin = as.numeric(outcome_f) - 0.42, ymax = as.numeric(outcome_f) + 0.42), fill = "grey92", inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(data = pa[ranking == "spot"], aes(xmin = ci_lo, xmax = ci_hi, colour = condition), width = 0, linewidth = 0.35, orientation = "y",
                position = position_dodge2(width = 0.6, reverse = TRUE)) +
  geom_point(data = pa[ranking == "neighborhood"], aes(x = delta, colour = condition, shape = key), size = 1.6, stroke = 0.45,
             position = position_dodge2(width = 0.6, reverse = TRUE)) +
  geom_point(data = pa[ranking == "spot"], aes(x = delta, colour = condition, fill = condition, shape = key), size = 1.7, stroke = 0.3,
             position = position_dodge2(width = 0.6, reverse = TRUE)) +
  facet_wrap(~ facet, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = PAL_COND, name = NULL, guide = guide_legend(order = 1, override.aes = list(shape = 21, fill = unname(PAL_COND[c("CON","NHD")]), stroke = 0.3, size = 1.8))) +
  scale_fill_manual(values = PAL_COND, guide = "none") +
  scale_shape_manual(values = SHAPES, name = NULL, drop = FALSE, guide = guide_legend(order = 2, override.aes = list(fill = "grey35", colour = "grey35", size = 1.8, stroke = 0.45))) +
  labs(x = "Cliff's δ, microglia-rich vs microglia-poor spots (top vs bottom quartile)", y = NULL) +
  scale_x_continuous(breaks = function(lim) seq(-1, 1, by = if (diff(lim) > 1.3) 0.5 else 0.25), expand = expansion(mult = c(0.06, 0.12))) +   # the wider (white-matter) facet takes 0.5 steps so its tick labels stay apart
  theme_pub(base_size = 7) + theme(legend.position = "right", strip.text = element_text(face = "plain", size = 6.6, colour = "black"),
                                  panel.grid.major.y = element_blank(), axis.title.x = element_text(size = 6.5, colour = "black"),
                                  legend.spacing.y = unit(4, "pt"), legend.key.height = unit(9, "pt"), legend.text = element_text(size = 6, colour = "black"),
                                  panel.spacing.x = unit(14, "pt"))   # the two facets' end tick labels must not touch
ggsave(file.path(OUT, "S6a_proximity_dots.png"), p1, width = 6.6, height = 3.3, dpi = 600, device = ragg::agg_png)
ggsave(file.path(OUT, "S6a_proximity_dots.pdf"), p1, width = 6.6, height = 3.3, device = grDevices::quartz, type = "pdf")   # quartz: the delta glyph survives in the PDF
# ---- Figure 6 panel f: the half-row version (left of the last Fig-6 row; panel g = log2-ratio, script 114) ----
# Same statistics and rows as S6a, re-typeset for 3.45 x 3.50 in: rows grouped by family in a left strip,
# short display labels (a recode of the plotted factor only — OUT_LAB is part of the cache KEY and the
# table vocabulary, so it is not edited), and the four-level combinatorial shape key split into three
# two-key guides (condition | section | spot vs neighborhood).
F6P <- file.path(PROJ, "figures", "Figure_6", "panels")
THEME_F6_ROW <- theme_pub(base_size = 7) + theme(   # shared verbatim with 114 (panel g) so the two halves of the row match
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3), axis.line = element_blank(),
  axis.text  = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 6.4, colour = "black"),
  strip.text = element_text(size = 6.6, colour = "black", face = "plain"),
  legend.position = "bottom", legend.text = element_text(size = 6, colour = "black"),
  legend.key.size = unit(9, "pt"), legend.key.spacing.x = unit(4, "pt"),
  legend.margin = margin(t = -2), legend.box.spacing = unit(3, "pt"), plot.margin = margin(3, 3, 3, 3))
DISP_LAB <- c("Astrocyte reactive" = "Reactive astrocyte", "Homeostatic microglia" = "Homeostatic", "Cholesterol (sterol arm)" = "Cholesterol",
              "Structural myelin" = "Myelin", "Glycolytic shift" = "Glycolysis", "FTL/HAMP (iron)" = "Iron (FTL/HAMP)", "SPP1 (DAM-2 arm)" = "SPP1 / DAM-2",
              "Housekeeping (null)" = "Housekeeping", "Astrocyte abundance (c2l)" = "Astrocyte", "Oligodendrocyte abundance (c2l)" = "Oligodendrocyte",
              "Excitatory-neuron abundance (c2l)" = "Excitatory neuron")   # panel-f display names, keyed by the OUT_LAB / table vocabulary; unmapped rows keep their name
FAMILY <- c("Complement / MHC-II" = "Microglia", "Antigen presentation" = "Microglia", "Homeostatic microglia" = "Microglia", "FTL/HAMP (iron)" = "Microglia", "SPP1 (DAM-2 arm)" = "Microglia",
            "Cholesterol (sterol arm)" = "Lipid", "Structural myelin" = "Lipid", "Fatty-acid / sphingomyelin" = "Lipid", "Lipid uptake / salvage" = "Lipid", "Glycolytic shift" = "Lipid",
            "Astrocyte reactive" = "Astrocyte", "Housekeeping (null)" = "Null", "astro" = "Abundance", "oligo" = "Abundance", "neu" = "Abundance")   # keyed by OUTCOMES
FAMILY_ORDER <- c("Microglia", "Lipid", "Astrocyte", "Null", "Abundance")
stopifnot(setequal(names(FAMILY), OUTCOMES), all(FAMILY %in% FAMILY_ORDER), all(names(DISP_LAB) %in% OUT_LAB))
disp_of <- function(lab) { m <- DISP_LAB[lab]; ifelse(is.na(m), lab, m) }
DISP_LEVELS <- rev(disp_of(unname(OUT_LAB[PANEL_ORDER])))        # PANEL_ORDER kept within each family (free y drops the rows a family does not carry)
pf <- copy(pa)
pf[, disp := factor(disp_of(outcome), levels = DISP_LEVELS)]
pf[, family := factor(FAMILY[names(OUT_LAB)[match(outcome, OUT_LAB)]], levels = FAMILY_ORDER)]
pf <- pf[!family %in% c("Abundance", "Null")]; pf[, family := droplevels(family)]   # panel f = program scores only: abundance is panel g, the housekeeping null stays in the table and Supp Fig 6
stopifnot(!anyNA(pf$disp), !anyNA(pf$family))
say("\n== panel f rows per family ==")
print(pf[, .(rows = paste(unique(as.character(disp)), collapse = " | ")), by = family])
# Null band per grey-matter row: the rect y-position is the row index within the family panel (free y), not the panel-wide factor code
nulls_f <- pf[facet == "Grey matter" & ranking == "spot", .(lo = min(null_lo, na.rm = TRUE), hi = max(null_hi, na.rm = TRUE)), by = .(facet, family, disp)]
nulls_f[, row := match(as.character(disp), DISP_LEVELS[DISP_LEVELS %in% as.character(disp)]), by = family]
stopifnot(!anyNA(nulls_f$row))
p1f <- ggplot(pf, aes(y = disp, group = dodge)) +
  geom_rect(data = nulls_f, aes(xmin = lo, xmax = hi, ymin = row - 0.42, ymax = row + 0.42), fill = "grey92", inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(data = pf[ranking == "spot"], aes(xmin = ci_lo, xmax = ci_hi, colour = condition), width = 0, linewidth = 0.35, orientation = "y",
                position = position_dodge2(width = 0.6, reverse = TRUE)) +
  geom_point(data = pf[ranking == "neighborhood"], aes(x = delta, colour = condition, shape = section, alpha = ranking), fill = NA, size = 1.4, stroke = 0.45,
             position = position_dodge2(width = 0.6, reverse = TRUE)) +
  geom_point(data = pf[ranking == "spot"], aes(x = delta, colour = condition, fill = condition, shape = section, alpha = ranking), size = 1.5, stroke = 0.3,
             position = position_dodge2(width = 0.6, reverse = TRUE)) +
  facet_grid(family ~ facet, scales = "free_y", space = "free_y", switch = "y") +   # one x range for both facets so 0.68 (grey) and 0.91 (white) are read on the same axis
  scale_colour_manual(values = PAL_COND, name = NULL, guide = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21, fill = unname(PAL_COND[c("CON","NHD")]), stroke = 0.3, size = 1.8))) +
  scale_fill_manual(values = PAL_COND, guide = "none") +
  scale_shape_manual(values = c(`section 1` = 21, `section 2` = 24), name = NULL,
                     guide = guide_legend(order = 2, nrow = 1, override.aes = list(fill = "grey35", colour = "grey35", size = 1.8, stroke = 0.45))) +
  scale_alpha_manual(values = c(spot = 1, neighborhood = 1), breaks = c("spot", "neighborhood"), name = NULL,   # alpha is a no-op carrier for the filled-vs-open key
                     guide = guide_legend(order = 3, nrow = 1, override.aes = list(shape = 21, fill = c("grey35", NA), colour = "grey35", size = 1.8, stroke = 0.45))) +
  labs(x = "Cliff's δ, microglia-rich vs microglia-poor spots", y = NULL) +
  scale_x_continuous(breaks = scales::breaks_width(0.5), limits = range(c(-0.3, 1, pf$delta, pf$ci_lo, pf$ci_hi, pf$null_lo, pf$null_hi), na.rm = TRUE), expand = expansion(mult = c(0.06, 0.10))) +
  THEME_F6_ROW +
  theme(strip.placement = "outside", strip.switch.pad.grid = unit(4, "pt"),
        strip.text.y.left = element_text(angle = 0, size = 6, hjust = 1, colour = "black", margin = margin(l = 3, r = 3)),
        strip.background.y = element_rect(fill = "grey94", colour = NA),
        panel.spacing.x = unit(10, "pt"), panel.spacing.y = unit(3, "pt"),
        legend.box = "horizontal", legend.spacing.x = unit(7, "pt"),
        legend.location = "plot")   # centre the three-guide row on the PAGE (3.45 in), not on the panel area: on the panel area it ran off the right edge
WF <- 3.45; HF <- 2.850
ggsave(file.path(F6P, "F6f_microglia_proximity_dots.png"), p1f, width = WF, height = HF, dpi = 600, device = ragg::agg_png)
ggsave(file.path(F6P, "F6f_microglia_proximity_dots.pdf"), p1f, width = WF, height = HF, device = grDevices::quartz, type = "pdf")   # quartz: the delta glyph survives in the PDF
say("Figure 6 panel f: F6f_microglia_proximity_dots (%.2f x %.2f in, left half of the last row)", WF, HF)

say("wrote S6a_proximity_dots (6.6 x 3.3 in) to %s; panels b and c of Supplementary Fig. 6 are written by 114_visium_domain_composition_FH.R", OUT)
writeLines(capture.output(sessionInfo()), file.path(PROJ, "logs", "111_visium_proximity_severity_FH_sessionInfo.txt"))
cat("\n=== DONE ===\n", file = stderr())
