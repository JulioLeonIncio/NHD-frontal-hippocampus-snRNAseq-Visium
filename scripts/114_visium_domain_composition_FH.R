#!/usr/bin/env Rscript
# =============================================================================
# 114_visium_domain_composition_FH.R — Supplementary Fig. 6 panels b and c: cell2location composition of every spatial domain, CON beside NHD, and its NHD / CON ratio.
# -----------------------------------------------------------------------------
# What. Per spot, the cell2location q05 abundance of each of the eight reference types is
# turned into a proportion of that spot's total q05 abundance (a spot whose total is zero has
# no composition and is excluded, counted below). Per domain x condition the mean proportion
# per type is the composition (panel b, dot matrix); per domain the log2 ratio of the two
# condition means is the shift (panel c), with a 95 % interval from a block bootstrap over
# tissue blocks (6 array rows x 12 array columns within a section; the unit of script 111),
# blocks resampled with replacement within each condition, NBOOT draws, percentile interval.
#
# What it is not. A within-donor sampling interval: n = 2 sections per condition, one donor per
# condition, so condition is aliased with donor and no P value is given anywhere. Only domains
# with >= MIN_SPOTS spots in both conditions carry a ratio; the others are listed as excluded.
# Endothelial and pericyte abundances are inflated by the reference (the atlas over-calls the
# vascular types: 27.8 % + 32.9 % of control cortex, script 71); the two rows are kept but
# shaded grey and named in the legend.
#
# Domain vocabulary (order, palette, compartment) is read from integrated_domain_levels.csv
# through _visium_domains_FH.R; the spot -> domain assignment comes from spot_metadata.csv.gz
# (script 102), which is itself mapped through integrated_cluster_identity.csv.
#
# Inputs:  manuscript/Final_fig_and_tables/Data_deposition/Visium/spot_metadata.csv.gz (102)
#          <manifest outs>/spatial/tissue_positions_list.csv per section (array coordinates -> blocks; folder from VISIUM_OUTS_MANIFEST.json via _visium_outs_FH.R)
# Outputs: tables/visium_domain_composition_FH.csv
#          supplementary_tables/ST21_visium_domain_composition.csv   (shipped copy -> Supplementary Data 11, sheet Composition_per_domain)
#          figures/Supplementary/SuppFig6_microglia_proximity/{S6b_domain_dotmatrix, S6c_domain_log2ratio}.{png,pdf}
#          figures/Figure_6/panels/{F6h_domain_composition, F6g_domain_composition_log2ratio}.{png,pdf}
#
# Column DICTIONARY (visium_domain_composition_FH.csv / ST21):
#   cell_type            cell2location reference type (Astro | Endo | Micro-PVM | Neuron_Ex | Neuron_Inh | OPC | Oligo | Pericytes)
#   cell_type_label      display name used on the panels (Astrocyte, Endothelial, Microglia / PVM, Excitatory neuron, ...)
#   reference_inflated   TRUE for Endo and Pericytes (vascular abundance over-called by the reference; shaded grey on the panels)
#   domain               spatial domain (integrated_domain_levels.csv)
#   domain_order         superficial -> deep order of the domain (1-7)
#   compartment          grey | white | other
#   n_spots_CON          control spots in the domain with a defined composition (total q05 > 0)
#   n_spots_NHD          NHD spots in the domain with a defined composition
#   n_blocks_CON         tissue blocks (6 x 12 array units, within section) carrying those control spots
#   n_blocks_NHD         tissue blocks carrying those NHD spots
#   mean_prop_CON        mean per-spot proportion of the type (q05 of the type / total q05 over the eight types), control spots
#   mean_prop_NHD        the same, NHD spots
#   mean_prop_<section>  the same per section (CON_Frontal1, CON_Frontal2, NHD_Frontal1, NHD_Frontal2); NA when the section has no spot in the domain
#   log2_ratio           log2(mean_prop_NHD / mean_prop_CON); NA when the domain is excluded (see `included`)
#   ci_lo, ci_hi         2.5th and 97.5th percentiles of log2_ratio over NBOOT block-bootstrap draws (blocks resampled within condition)
#   included             TRUE when n_spots_CON >= MIN_SPOTS and n_spots_NHD >= MIN_SPOTS (the panel-c filter)
#   nboot, min_spots     the constants used (2000, 30)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ragg) })
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))
source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))   # vis_outs(sec) — outs folder per section from VISIUM_OUTS_MANIFEST.json
ROOT <- dirname(PROJ); TDIR <- file.path(PROJ, "tables")
OUT  <- file.path(PROJ, "figures", "Supplementary", "SuppFig6_microglia_proximity"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
NBOOT <- 2000L; MIN_SPOTS <- 30L; BLOCK <- c(rows = 6L, cols = 12L)   # block as script 111

CT8 <- c("Astro","Endo","Micro-PVM","Neuron_Ex","Neuron_Inh","OPC","Oligo","Pericytes")
CT_LAB <- c("Micro-PVM" = "Microglia / PVM", "Oligo" = "Oligodendrocyte", "OPC" = "OPC", "Astro" = "Astrocyte",
            "Neuron_Ex" = "Excitatory neuron", "Neuron_Inh" = "Inhibitory neuron", "Endo" = "Endothelial", "Pericytes" = "Pericyte")   # as Fig 6f (script 71)
CT_ORDER <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh","Endo","Pericytes")   # panel rows: glia, neurons, then the reference-inflated vascular pair
INFLATED <- c("Endo","Pericytes")
stopifnot(setequal(CT8, CT_ORDER), all(CT8 %in% names(PAL_CELLTYPE_VIVID)))

# ---- inputs --------------------------------------------------------------------------------
mdf <- file.path(PROJ, "manuscript", "Final_fig_and_tables", "Data_deposition", "Visium", "spot_metadata.csv.gz")
stopifnot("MISSING spot_metadata.csv.gz — run 102 first" = file.exists(mdf))
md <- fread(mdf)
stopifnot(all(c("spot_id","sample_id","condition","spatial_domain", paste0("c2l_q05_", CT8)) %in% names(md)),
          "spot_metadata carries a domain absent from integrated_domain_levels.csv" = all(unique(md$spatial_domain) %in% DOM_LEV))
setnames(md, paste0("c2l_q05_", CT8), CT8)
md[, total_q05 := rowSums(.SD), .SDcols = CT8]
n0 <- nrow(md); md <- md[total_q05 > 0]
say("spots: %d -> %d with a defined composition (%d with zero total q05 excluded)", n0, nrow(md), n0 - nrow(md))
for (ct in CT8) set(md, j = ct, value = md[[ct]] / md$total_q05)          # proportion of the spot's deconvolved content
stopifnot(abs(rowSums(md[, ..CT8]) - 1) < 1e-8)

pos <- rbindlist(lapply(unique(md$sample_id), function(s) {
  f <- file.path(vis_outs(s), "spatial", "tissue_positions_list.csv")   # manifest outs, as 111
  stopifnot(file.exists(f))
  p <- fread(f, header = FALSE, col.names = c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col"))
  p[, spot_id := paste0(s, "_", barcode)]; p[, .(spot_id, array_row, array_col)] }))
n1 <- nrow(md); md <- merge(md, pos, by = "spot_id")
stopifnot("array coordinates missing for some spots" = nrow(md) == n1, !anyNA(md$array_row))
md[, block := paste(sample_id, array_row %/% BLOCK[["rows"]], array_col %/% BLOCK[["cols"]])]
md[, domain := factor(spatial_domain, levels = DOM_LEV)]
md[, condition := factor(condition, levels = c("CON","NHD"))]

# ---- composition per domain x condition (panel b) --------------------------------------------
long <- melt(md[, c("spot_id","sample_id","condition","domain","block", CT8), with = FALSE], id.vars = c("spot_id","sample_id","condition","domain","block"),
             variable.name = "cell_type", value.name = "prop", variable.factor = FALSE)
comp <- long[, .(n_spots = .N, n_blocks = uniqueN(block), mean_prop = mean(prop)), by = .(domain, condition, cell_type)]
comp_sec <- long[, .(mean_prop = mean(prop)), by = .(domain, sample_id, cell_type)]
ncell <- unique(comp[, .(domain, condition, n_spots)])
say("\n== spots per domain x condition (defined composition) ==")
print(dcast(ncell, domain ~ condition, value.var = "n_spots", fill = 0L))

# ---- log2(NHD / CON) with a block-bootstrap interval (panel c) ------------------------------
# per domain x condition: per-block sums of the proportions (B x 8) and block sizes; a bootstrap
# draw is a multinomial re-weighting of the blocks, so 2,000 draws are 2,000 small matrix products
boot_means <- function(dd, n = NBOOT) {
  S  <- as.matrix(dd[, lapply(.SD, sum), by = block, .SDcols = CT8][, ..CT8])
  nb <- dd[, .N, by = block]$N
  B  <- length(nb)
  W  <- t(vapply(seq_len(n), function(i) tabulate(sample.int(B, B, replace = TRUE), nbins = B), integer(B)))   # n x B block multiplicities
  (W %*% S) / as.numeric(W %*% nb)          # n x 8 bootstrap means of the proportion
}
ratio <- list(); excluded <- character(0)
for (dm in DOM_LEV) {
  dc <- md[domain == dm & condition == "CON"]; dn <- md[domain == dm & condition == "NHD"]
  ok <- nrow(dc) >= MIN_SPOTS && nrow(dn) >= MIN_SPOTS
  if (!ok) { excluded <- c(excluded, sprintf("%s (CON %d, NHD %d)", dm, nrow(dc), nrow(dn))) }
  mC <- if (nrow(dc)) colMeans(dc[, ..CT8]) else setNames(rep(NA_real_, 8), CT8)
  mN <- if (nrow(dn)) colMeans(dn[, ..CT8]) else setNames(rep(NA_real_, 8), CT8)
  lr <- lo <- hi <- setNames(rep(NA_real_, 8), CT8)
  if (ok) {
    zero <- mC == 0 | mN == 0                                     # a zero mean would give log2 = +/-Inf; reported as NA and logged
    if (any(zero)) say("  %s: mean proportion of %s is 0 in one condition — ratio reported as NA", dm, paste(CT8[zero], collapse = ", "))
    lr[!zero] <- log2(mN[!zero] / mC[!zero])
    bC <- boot_means(dc); bN <- boot_means(dn)
    bl <- log2(bN / bC); bl[!is.finite(bl)] <- NA                # a draw with a zero mean is dropped from the percentile (counted below)
    ndrop <- colSums(is.na(bl))
    if (any(ndrop > 0)) say("  %s: %s bootstrap draws with a zero mean dropped (%s)", dm, paste(ndrop[ndrop > 0], collapse = "/"), paste(CT8[ndrop > 0], collapse = ", "))
    q <- apply(bl, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    lo[!zero] <- q[1, !zero]; hi[!zero] <- q[2, !zero]
  }
  ratio[[dm]] <- data.table(domain = dm, cell_type = CT8, n_spots_CON = nrow(dc), n_spots_NHD = nrow(dn),
                            n_blocks_CON = uniqueN(dc$block), n_blocks_NHD = uniqueN(dn$block),
                            mean_prop_CON = unname(mC), mean_prop_NHD = unname(mN),
                            log2_ratio = unname(lr), ci_lo = unname(lo), ci_hi = unname(hi), included = ok)
  say("  %-12s CON %4d spots / %3d blocks  NHD %4d spots / %3d blocks  %s", dm, nrow(dc), uniqueN(dc$block), nrow(dn), uniqueN(dn$block), if (ok) "ratio" else "EXCLUDED")
}
res <- rbindlist(ratio)
secw <- dcast(comp_sec, domain + cell_type ~ sample_id, value.var = "mean_prop")
setnames(secw, setdiff(names(secw), c("domain","cell_type")), paste0("mean_prop_", setdiff(names(secw), c("domain","cell_type"))))
res <- merge(res, secw, by = c("domain","cell_type"), all.x = TRUE)
res[, `:=`(cell_type_label = unname(CT_LAB[cell_type]), reference_inflated = cell_type %in% INFLATED,
           domain_order = match(domain, DOM_LEV), compartment = unname(DOM_COMP[domain]), nboot = NBOOT, min_spots = MIN_SPOTS)]
setcolorder(res, c("cell_type","cell_type_label","reference_inflated","domain","domain_order","compartment",
                   "n_spots_CON","n_spots_NHD","n_blocks_CON","n_blocks_NHD","mean_prop_CON","mean_prop_NHD",
                   grep("^mean_prop_(CON|NHD)_Frontal", names(res), value = TRUE),
                   "log2_ratio","ci_lo","ci_hi","included","nboot","min_spots"))
res <- res[order(domain_order, match(cell_type, CT_ORDER))]
fwrite(res, file.path(TDIR, "visium_domain_composition_FH.csv"))
dir.create(file.path(PROJ, "supplementary_tables"), showWarnings = FALSE)
fwrite(res, file.path(PROJ, "supplementary_tables", "ST21_visium_domain_composition.csv"))   # shipped copy -> Supplementary Data 11, sheet Composition_per_domain (script 100)
say("\nexcluded from panel c (< %d spots in one condition): %s", MIN_SPOTS, if (length(excluded)) paste(excluded, collapse = "; ") else "none")
say("\n== log2(NHD / CON) of the mean proportion, included domains (for the legend) ==")
print(res[included == TRUE, .(domain, cell_type, n_CON = n_spots_CON, n_NHD = n_spots_NHD, mean_CON = signif(mean_prop_CON, 3), mean_NHD = signif(mean_prop_NHD, 3),
                              log2_ratio = round(log2_ratio, 2), ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi))], nrows = 80)

# ---- panel b: dot matrix ---------------------------------------------------------------------
pb <- copy(comp)
pb[, cell_type := factor(cell_type, levels = rev(CT_ORDER))]
pb <- merge(pb, ncell[, .(domain, condition, n_col = n_spots)], by = c("domain","condition"))
pb[, col_lab := sprintf("%s\n%s", as.character(domain), format(n_col, big.mark = ",", trim = TRUE))]   # domain over its spot count
pb[, col_lab := factor(col_lab, levels = unique(pb[order(condition, domain)]$col_lab))]              # facet-local columns in level order
shade <- data.table(ymin = which(rev(CT_ORDER) %in% INFLATED) - 0.5)[, ymax := ymin + 1]     # grey band behind the reference-inflated rows
# one builder for the supplement (facets side by side, 6.6 in) and the Fig-6 copy (facets stacked, 3.3 in)
dot_matrix <- function(stacked = FALSE) {
  ggplot(pb, aes(col_lab, cell_type)) +
    geom_rect(data = shade, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax), fill = "grey92", colour = NA, inherit.aes = FALSE) +
    geom_point(aes(size = mean_prop, fill = mean_prop), shape = 21, colour = "grey30", stroke = 0.25) +
    facet_wrap(~ condition, nrow = if (stacked) 2 else 1, scales = "free_x") +
    scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, NA), name = "mean proportion of spot content",
                         guide = guide_colourbar(order = 1, title.position = "top", direction = "horizontal")) +
    scale_size_area(max_size = if (stacked) 3.6 else 4.2, limits = c(0, NA), breaks = c(0.05, 0.15, 0.3), name = NULL,
                    guide = guide_legend(order = 2, direction = "horizontal", override.aes = list(fill = "grey60"))) +
    scale_y_discrete(labels = CT_LAB) +
    labs(x = "spatial domain (spots with a composition)", y = NULL) +
    theme_pub(base_size = 7) +
    theme(axis.text.x = element_text(size = 5.6, colour = "black", lineheight = 0.9), axis.text.y = element_text(size = 6.2, colour = "black"),
          axis.title.x = element_text(size = 6.4, colour = "black"), strip.text = element_text(size = 6.8, colour = "black", face = "plain"),
          strip.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3), axis.line = element_blank(),
          panel.grid.major = element_blank(), panel.spacing.x = unit(12, "pt"), panel.spacing.y = unit(4, "pt"),
          legend.position = "bottom", legend.box = if (stacked) "vertical" else "horizontal", legend.title = element_text(size = 6, colour = "black"),
          legend.text = element_text(size = 5.8, colour = "black"), legend.key.width = unit(0.6, "cm"), legend.key.height = unit(0.2, "cm"),
          legend.spacing.x = unit(14, "pt"), legend.spacing.y = unit(0, "pt"), legend.margin = margin(t = 0, b = 0),
          legend.box.spacing = unit(3, "pt"), plot.margin = margin(3, 3, 3, 3))
}
p_b <- dot_matrix(stacked = FALSE)
WB <- 6.6; HB <- 2.95
ggsave(file.path(OUT, "S6b_domain_dotmatrix.png"), p_b, width = WB, height = HB, dpi = 600, device = ragg::agg_png)
ggsave(file.path(OUT, "S6b_domain_dotmatrix.pdf"), p_b, width = WB, height = HB, useDingbats = FALSE)

# ---- panel c: log2 ratio with block-bootstrap intervals ----------------------------------------
pc <- res[included == TRUE & !is.na(log2_ratio)]
pc[, domain := factor(domain, levels = DOM_LEV[DOM_LEV %in% unique(domain)])]
pc[, cell_type := factor(cell_type, levels = CT_ORDER)]
pal_c <- PAL_CELLTYPE_VIVID[CT_ORDER]
yl <- range(c(pc$ci_lo, pc$ci_hi, 0)); yl <- yl + c(-0.05, 0.05) * diff(yl)
# shared verbatim with 111 (panel f) so the two halves of the last Fig-6 row match
THEME_F6_ROW <- theme_pub(base_size = 7) + theme(
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3), axis.line = element_blank(),
  axis.text  = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 6.4, colour = "black"),
  strip.text = element_text(size = 6.6, colour = "black", face = "plain"),
  legend.position = "bottom", legend.text = element_text(size = 6, colour = "black"),
  legend.key.size = unit(9, "pt"), legend.key.spacing.x = unit(4, "pt"),
  legend.margin = margin(t = -2), legend.box.spacing = unit(3, "pt"), plot.margin = margin(3, 3, 3, 3))
# "right" = Supplementary Fig. 6c (6.6 in, legend at the right); "bottom" = Figure 6 panel g (3.05 x 3.50 in half row:
# wider dodge, smaller dots, two-row legend, the cell2location statement carried in the y title)
log2_panel <- function(legend = c("right","bottom"), dom_drop = character()) {
  legend <- match.arg(legend); fig6 <- legend == "bottom"
  dw <- if (fig6) 0.85 else 0.7
  p <- ggplot(pc, aes(domain, log2_ratio, colour = cell_type, group = cell_type)) +
    geom_vline(xintercept = seq_len(nlevels(pc$domain) - 1) + 0.5, colour = "grey88", linewidth = 0.25) +   # domain separators
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, linewidth = if (fig6) 0.35 else 0.4, position = position_dodge(width = dw)) +
    geom_point(size = if (fig6) 1.1 else 1.5, position = position_dodge(width = dw)) +
    scale_colour_manual(values = pal_c, labels = CT_LAB[CT_ORDER], name = NULL,
                        guide = guide_legend(nrow = if (fig6) 3 else 8, byrow = TRUE, override.aes = list(size = if (fig6) 2 else 2.2, linewidth = 0))) +   # fig6: 3 rows x 3 columns — two rows (4 columns) run 3.4 in wide at 6 pt and clip OPC / Pericyte at 3.05 in
    scale_y_continuous(limits = yl, expand = expansion(mult = 0))
  if (!fig6) return(p +
    labs(x = "spatial domain", y = "log2(NHD / CON), mean proportion") +
    theme_pub(base_size = 7) +
    theme(axis.text.x = element_text(size = 6.2, colour = "black"), axis.text.y = element_text(size = 6, colour = "black"),
          axis.title = element_text(size = 6.4, colour = "black"), panel.grid.major.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3), axis.line = element_blank(),
          legend.position = legend, legend.text = element_text(size = 6, colour = "black"), legend.key.height = unit(9, "pt"),
          legend.key.width = unit(8, "pt"), legend.margin = margin(t = -2), legend.box.spacing = unit(3, "pt"), plot.margin = margin(3, 3, 3, 3)))
  # Fig 6 panel g: narrow forest — domains as rows (L2/3 at the top), log2 ratio along x,
  # the eight cell types dodged within each row; "cell2location" stated in the panel title and the x-axis title
  pcf <- copy(pc)[!domain %in% dom_drop]; pcf[, domain := factor(domain, levels = rev(setdiff(levels(domain), dom_drop)))]   # Fig 6g shows the cortical layers, L1/pia and WM only
  ggplot(pcf, aes(x = log2_ratio, y = domain, colour = cell_type, group = cell_type)) +
    geom_hline(yintercept = seq_len(nlevels(pcf$domain) - 1) + 0.5, colour = "grey88", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0, linewidth = 0.35, orientation = "y", position = position_dodge(width = 0.85)) +
    geom_point(size = 1.1, position = position_dodge(width = 0.85)) +
    scale_colour_manual(values = pal_c, labels = CT_LAB[CT_ORDER], name = NULL,
                        guide = guide_legend(ncol = 1, override.aes = list(size = 2, linewidth = 0))) +
    scale_x_continuous(limits = yl, expand = expansion(mult = 0), breaks = scales::breaks_width(1)) +
    labs(y = NULL, x = "log2 (NHD / CON)", title = "cell2location abundance") +
    THEME_F6_ROW +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(size = 6.4, colour = "black", face = "plain", hjust = 0.5, margin = margin(b = 2)),
          legend.position = "right", legend.justification = "center", legend.key.width = unit(8, "pt"), legend.key.spacing.y = unit(2, "pt"),   # legend on the right, vertically centred
          legend.margin = margin(l = 2), legend.box.spacing = unit(3, "pt"))
}
p_c <- log2_panel("right")
WC <- 6.6; HC <- 2.35
ggsave(file.path(OUT, "S6c_domain_log2ratio.png"), p_c, width = WC, height = HC, dpi = 600, device = ragg::agg_png)
ggsave(file.path(OUT, "S6c_domain_log2ratio.pdf"), p_c, width = WC, height = HC, useDingbats = FALSE)

# ---- Figure-6 copies: F6h (dot matrix, stacked) and panel g (log2 ratio, right half of the last Fig-6 row beside F6f from 111) ----
F6P <- file.path(PROJ, "figures", "Figure_6", "panels")
p_h <- dot_matrix(stacked = TRUE); WH <- 3.3; HH <- 4.05
ggsave(file.path(F6P, "F6h_domain_composition.png"), p_h, width = WH, height = HH, dpi = 600, device = ragg::agg_png)
ggsave(file.path(F6P, "F6h_domain_composition.pdf"), p_h, width = WH, height = HH, useDingbats = FALSE)
p_g <- log2_panel("bottom", dom_drop = c("InN", "Vasc/immune")); WI <- 3.10; HI <- 2.85   # InN + Vasc/immune dropped from the panel; both stay in the table (ST21 / Supplementary Data 24)   # 10 % narrower than 3.45; same height as F6f
ggsave(file.path(F6P, "F6g_domain_composition_log2ratio.png"), p_g, width = WI, height = HI, dpi = 600, device = ragg::agg_png)
ggsave(file.path(F6P, "F6g_domain_composition_log2ratio.pdf"), p_g, width = WI, height = HI, useDingbats = FALSE)
say("Figure-6 copies: F6h_domain_composition (%.2f x %.2f in), F6g_domain_composition_log2ratio (%.2f x %.2f in)", WH, HH, WI, HI)
say("wrote S6b_domain_dotmatrix (%.2f x %.2f in) and S6c_domain_log2ratio (%.2f x %.2f in) to %s", WB, HB, WC, HC, OUT)
writeLines(capture.output(sessionInfo()), file.path(PROJ, "logs", "114_visium_domain_composition_FH_sessionInfo.txt"))
cat("\n=== DONE ===\n", file = stderr())
