#!/usr/bin/env Rscript
# =============================================================================
# 67_visium_supplychain_spatial_FH.R — Figure 6: "6_supplychain_spatial" and its quantitative companion "6_supplychain_matched".
# -----------------------------------------------------------------------------
# Why: "to avoid showing just myelin genes we can improve what
# we show in terms of module scores after learning what we learnt from the oligo
# OPC figure."
#
# The shipped spatial panels score one programme — myelin — which is the half of
# the Figure-4 argument that is not the finding. The finding is the dissociation:
# structural myelin held (or raised) per cell while every lipid arm that has to
# supply the membrane falls. This panel puts that dissociation in tissue.
#
# All four arms clear the 55 um floor on this slide (CON_Frontal1, fraction of
# spots): sterol 13/13 genes with DHCR24 68%, FDFT1 61%, MSMO1 61%, HMGCR 46%;
# SCD 78%, ELOVL5 36%; UGT8 22%, GALC 22%; SLC44A1 53%. So this is measurable
# spatial biology, not a stretch.
#
# Two controls baked in, both learned the hard way in this project:
#   * depth. NHD spots carry half the depth of CON spots (median UMI per section
#     in tables/visium_dm_qc_FH.csv, written by 65).
#     Scores come from 65_visium_prep_depthmatched_FH.R, where every spot is
#     thinned to a common 3,000 UMI first. Without it, five unrelated programmes
#     all "fall" by a similar amount and the panel measures library size.
#   * gene dominance. A summed-share module is carried by its most abundant gene
#     (FTH1 alone drove a -0.92 of a -0.69 "stress" move; PLP1 dwarfs the myelin
#     set). Scores are AddModuleScore with expression-matched control bins, and the
#     housekeeping null is plotted alongside so the reader can see the baseline
#     does not move.
#
# n = 2 CON + 2 NHD sections, one donor per condition: per-section points, ratio of
# section means, no p-values (condition is perfectly aliased with donor).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales); library(patchwork); library(ggh4x)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)

MOD <- file.path(TDIR, "visium_dm_module_spots_FH.csv")
stopifnot("MISSING per-spot module table — run 65_visium_prep_depthmatched_FH.R" =
            file.exists(MOD))
d <- read.csv(MOD, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("prep did not carry coordinates — re-run 65 (v2)" = all(c("x","y") %in% names(d)))

# The panel now opens with the glial response, then the oligodendrocyte
# supply chain ("6_supplychain_spatial is about myelination, why don't we make an
# upper section about the findings in astrocytes and microglia?"). Putting the two
# reactive limbs above the lipid arms makes the figure read in causal order -- who is
# activated, and what fails -- instead of starting in the middle of the story. It also
# lets the reader see, in the same four sections, that the microglial and astrocytic
# signals rise in the same tissue where the lipid arms fall.
GLIA <- c("Complement / MHC-II", "Reactive astrocyte")
ARMS <- c("Structural myelin", "Cholesterol (sterol arm)",
          "Fatty-acid / sphingomyelin", "Galactolipid", "Lipid uptake / salvage")
NULLM <- "Housekeeping (null)"
have <- intersect(c(ARMS, NULLM), names(d))
stopifnot("prep table lacks the supply-chain arms — re-run 65 (v2)" = length(have) >= 5)
cat("arms present:", paste(have, collapse = " | "), "\n")

SEC <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")
d$sample_id <- factor(d$sample_id, levels = SEC)
d$Condition <- factor(d$Condition, levels = c("CON","NHD"))

# the astrocyte limb has no precomputed module, so it is built here from the two
# reactive genes the prep exports, each z-scored within section first so that neither
# gene's absolute level can dominate the pair.
zc <- function(v) { s <- stats::sd(v, na.rm = TRUE); if (!is.finite(s) || s == 0) return(v*0)
                    (v - mean(v, na.rm = TRUE)) / s }
astro_genes <- intersect(c("gene_GFAP", "gene_SERPINA3"), names(d))
stopifnot("astrocyte genes missing — re-run 65 (v3)" = length(astro_genes) == 2)
d <- d %>% group_by(sample_id) %>%
  mutate(`Astrocyte: reactive` = rowMeans(cbind(zc(.data[[astro_genes[1]]]),
                                                zc(.data[[astro_genes[2]]])))) %>%
  ungroup() %>%
  mutate(`Reactive astrocyte` = `Astrocyte: reactive`)
cat(sprintf("astrocyte limb built from %s\n", paste(astro_genes, collapse = " + ")))

ROWS <- c(GLIA, ARMS)
long <- d %>%
  select(spot_id, sample_id, Condition, x, y, oligo_prop, all_of(ROWS)) %>%
  pivot_longer(all_of(ROWS), names_to = "arm", values_to = "score") %>%
  mutate(arm    = factor(arm, levels = ROWS),
         block  = factor(ifelse(arm %in% GLIA, "Glial response",
                                "Oligodendrocyte supply chain"),
                         levels = c("Glial response", "Oligodendrocyte supply chain")))

# ---------------------------------------------------------------------------
# PANEL 1 — the arms in space
# ---------------------------------------------------------------------------
# Per-section coordinates are centred and scaled by the larger axis span so all four
# sections occupy the same box with aspect preserved (the house spmap grammar).
long <- long %>% group_by(sample_id) %>%
  mutate(.sp = max(diff(range(x)), diff(range(y))),
         xr = (x - mean(range(x)))/.sp, yr = (y - mean(range(y)))/.sp) %>% ungroup()
# spots below the 3,000-UMI depth floor (nearly all white matter: control WM median ~200 UMI) carry no
# depth-matched score; they are drawn as a light-grey underlay in the same per-section frame so the
# section keeps its shape and the reader sees "not measured", not "no signal"
.norm <- long %>% group_by(sample_id) %>% summarise(.sp = max(diff(range(x)), diff(range(y))), cx = mean(range(x)), cy = mean(range(y)), .groups = "drop")
.all  <- read.csv(file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN", "spot_coords.csv"), stringsAsFactors = FALSE) %>%
  transmute(spot_id = cell, x, y, sample_id = sub("_[ACGT]+-1$", "", cell)) %>% filter(sample_id %in% unique(long$sample_id))
below_floor <- .all %>% filter(!spot_id %in% unique(long$spot_id)) %>% inner_join(.norm, by = "sample_id") %>%
  mutate(xr = (x - cx)/.sp, yr = (y - cy)/.sp, sample_id = factor(sample_id, levels = levels(long$sample_id)))

# One colour scale per ARM (free), since the arms sit at different absolute levels;
# limits are the 1st-99th percentile so a few extreme spots cannot flatten the map
# ggplot cannot give each facet row its own colour scale, and the arms sit at very
# different absolute levels — with one shared ramp only the myelin row had visible
# colour and the other four rendered uniformly black. Rescale within each arm to its
# own 1st-99th percentile so every row uses the full ramp; the legend is therefore
# relative within arm, which the axis title now states.
# per-section relative SCALING ("make each section have their own
# relative z score"). Each map is now rescaled to its own section-and-row distribution
# (1st-99th percentile), not to the row pooled across sections.
#
# Read the consequence before writing the legend. This deliberately removes the
# between-section level difference from the maps: every panel now spans low-to-high
# within itself, so a section cannot look uniformly brighter than another. What the maps
# now show is where a programme sits inside each tissue -- its spatial pattern -- and
# not how much higher it is in NHD. The CON-versus-NHD magnitude lives in
# 6_supplychain_matched (quantified at matched oligodendrocyte content) and in
# F6d_supplychain_genes (log2 NHD/CON per cortical band). The legend must not read a
# level difference off this panel.
lim <- long %>% group_by(arm, sample_id) %>%
  summarise(lo = quantile(score, 0.01), hi = quantile(score, 0.99), .groups = "drop")
long <- long %>% left_join(lim, by = c("arm", "sample_id")) %>%
  mutate(score_c = pmin(pmax((score - lo)/(hi - lo), 0), 1))

# Second scaling, matched to the gene HEATMAP ("do it in the way it
# replicates as much as possible what we see in the gene heatmaps"). The showcase panels
# below use this one; the seven-row grid keeps the per-section stretch above.
# 73_visium_gene_heatmap_FH.R centres each gene on its own mean and then divides by the
# pathway block's SD. The spatial analogue of "a gene" is the programme score itself, so
# here each programme is centred on its own mean and divided by its own SD, pooled across
# all four SECTIONS. Consequences, both intended:
#   * the unit becomes SD-within-pathway, identical to the heatmap's, so a red spot and a
#     red heatmap cell now mean the same thing and the legend may cross-reference them;
#   * nothing is clipped (the percentile squash above is not applied here);
#   * pooling across sections RESTORES the CON-versus-NHD level difference to the maps,
#     which the per-section stretch had deliberately removed.
long <- long %>% group_by(arm) %>%
  mutate(z = (score - mean(score)) / stats::sd(score)) %>% ungroup()

# No smoothing.
# The trade is visible and deliberate -- a single
# 55 um spot carries a noisy module score, so the maps are grainier than the gene heatmap
# (which averages ~1,400 spots per cell and is smooth by construction). Nothing here is
# filtered, averaged or interpolated; the grain is the measurement.

# Two single-programme PANELS ("for the spatial part we can just
# showcase two of them: Complement/MHC-II and Fatty-acid/sphingomyelin, generate two
# separate panels"). One programme per panel, four sections across, which lets each map
# be ~2.5x the area it had in the seven-row grid -- the point of showing these two
# spatially is that the reader can see the tissue, and at seven rows they could not.
#
# The programme is named in the COLOURBAR title rather than in a left-hand strip: with a
# single row a strip would cost ~0.9 in of width to say one thing, and it keeps the panel
# free of any text that is not a data label.
#
# The seven-row grid (all programmes) is still built above and remains available for the
# supplement; these two are the main-figure panels.
# Raw SPOTS only. These two panels draw the per-spot value exactly as measured -- no
# neighbour averaging, no interpolation, no percentile squash. A 6-nearest-neighbour
# smooth was tried and rejected; FNN is deliberately not imported by
# this script, so reintroducing one would be a visible change, not a silent one.
SHOWCASE <- c("Complement / MHC-II", "Fatty-acid / sphingomyelin")
YSPAN <- long %>% group_by(sample_id) %>%
  summarise(sp = diff(range(yr)), .groups = "drop") %>% pull(sp) %>% max()
cat(sprintf("tallest section spans %.2f of the cell -> row height %.2f in\n",
            YSPAN, 0.95 * YSPAN))
stopifnot("a showcase programme is missing from the scored rows" =
            all(SHOWCASE %in% levels(long$arm)))

GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
# The compression is fitted, not copied. The gene heatmap uses power 0.60 on a |z| that
# tops out at 3.72, because it plots band-AVERAGED gene means. Per-spot scores have far
# heavier tails -- |z| reaches 9.42 here -- so reusing 0.60 would push the whole tissue
# into the pale centre and a handful of outlying spots would own the ramp. The power is
# therefore solved so that the 90th percentile of |z| lands at 0.60 of the half-range,
# which is where the heatmap's own 90th percentile sits. Same construction, same feel,
# fitted to this data's spread. Still no clipping.
.zs  <- abs(long$z[long$arm %in% SHOWCASE])
ZMAX <- max(.zs)
Q90  <- as.numeric(quantile(.zs, 0.90))
POW  <- max(0.45, log(0.60) / log(Q90 / ZMAX))   # floor: below this the middle over-expands
compress <- function(v) sign(v) * (abs(v) / ZMAX)^POW
BR_Z <- c(-2, 0, 2); BR_Z <- BR_Z[abs(BR_Z) <= ZMAX]   # three ticks: the compressive mapping places ±1 and ±2 too close for two labels
cat(sprintf("showcase scaling: SD within pathway | |z| max %.2f, 90th pct %.2f -> power %.2f\n",
            ZMAX, Q90, POW))

spatial_one <- function(programme, basename) {
  dd <- long %>% filter(arm == programme) %>% mutate(pos = compress(z))
  p <- ggplot(dd, aes(xr, yr, colour = pos)) +
    geom_point(data = below_floor, aes(xr, yr), colour = "grey88", size = 0.58, stroke = 0, inherit.aes = FALSE) +   # below the depth floor
    geom_point(size = 0.58, stroke = 0) +   # sized to close the grid gaps
    # The magma ramp is KEPT ("can we just keep the colour scale we
    # had before"). Only the COLOURS revert -- the quantity, the units and the ticks stay
    # exactly as the gene heatmap defines them: centred, divided by the pathway SD, drawn
    # through the same compressive mapping, ticked in SD. So a bright spot still means
    # the same number a red heatmap cell means; it is the palette that differs, not the
    # measurement. Magma also solves the mid-value problem a diverging ramp has on a
    # spatial map, where near-zero spots rendered near-white and punched holes in the
    # tissue.
    # Dark = high, pale = low (magma reversed); one horizontal bar centred above the four maps
    scale_colour_gradientn(colours = c("#FFFFFF", rev(viridisLite::magma(9, begin = 0.08, end = 0.97))), limits = c(-1, 1),   # white = low ... black = high
                           breaks = compress(BR_Z), labels = ifelse(BR_Z == 0, "0", sprintf("%+.0f", BR_Z)),
                           name = paste0(programme, " (SD within pathway)"),
                           guide = guide_colourbar(title.position = "left", title.vjust = 0.9, direction = "horizontal")) +
    facet_wrap(~ sample_id, nrow = 1,
               labeller = labeller(sample_id = function(x) sub("_Frontal", " ", x))) +
    coord_equal() +
      # The panel row is sized to the TALLEST section's actual aspect, not to a square.
    # Each section is normalised by its larger span, so a section wider than it is tall
    # has a y-range < 1 and a square cell leaves a band of blank above and below it.
    # Sizing the row to the real y-extent removes that dead space.
    ggh4x::force_panelsizes(rows = unit(0.95 * YSPAN, "in"),
                            cols = unit(0.95, "in")) +
    theme_pub(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(), axis.line = element_blank(),
          strip.text.x = element_text(size = 6.4, colour = "black"),
          strip.background = element_blank(),
          panel.spacing = unit(0.04, "cm"),
          legend.position = "top", legend.direction = "horizontal",
          legend.justification = "centre",
          legend.key.width = unit(1.1, "cm"), legend.key.height = unit(0.18, "cm"),
          legend.title = element_text(size = 5.8, colour = "black"),
          legend.text = element_text(size = 5.8, colour = "black"),
          legend.margin = margin(b = 0), legend.box.spacing = unit(0.06, "cm"),
          plot.margin = margin(1, 1, 1, 1))
  W <- 4.35; H <- 0.95 * YSPAN + 0.60   # bar + strip above the maps
  ggsave(file.path(PANEL, paste0(basename, ".pdf")), p, width = W, height = H,
         useDingbats = FALSE)
  ggsave(file.path(PANEL, paste0(basename, ".png")), p, width = W, height = H, dpi = 600,
         device = ragg::agg_png)
  cat(sprintf("Saved %s.{png,pdf} (%.2f x %.2f in) — %s\n", basename, W, H, programme))
}
spatial_one(SHOWCASE[1], "F6e1_spatial_complement_MHCII")
spatial_one(SHOWCASE[2], "F6e2_spatial_fattyacid_sphingomyelin")

p_map <- ggplot(long, aes(xr, yr, colour = score_c)) +
  geom_point(size = 0.30, stroke = 0) +
  scale_colour_viridis_c(option = "magma", name = "relative within section",
                         breaks = c(0, 0.5, 1), labels = c("low","mid","high"),
                         guide = guide_colourbar(title.position = "top",
                                                 title.hjust = 0.5)) +
  # Strip LAYOUT (measured). A nested outer strip would cost 0.65 in of a 3.45 in panel
  # and is not text-sized (shortening "Oligodendrocyte supply chain" to "Oligo" changes
  # its width by nothing), so that level is dropped and 19% of the panel goes to the maps
  # (cell 0.49 -> 0.65 in). The two sections are still distinguished, by strip fill: the two
  # glial rows sit on a darker ground than the five supply-chain rows. Spell the section
  # names out once in the legend.
  ggh4x::facet_grid2(arm ~ sample_id, switch = "y",
                     labeller = labeller(
                       sample_id = function(x) sub("_Frontal", " ", x),
                       arm = function(x) sub(" / ", " /\n", sub(" \\(", "\n(", x))),
                     strip = ggh4x::strip_themed(
                       # Section identity now rides on the text, not on a filled box:
                       # the glial rows are set in a darker ink than the supply-chain rows.
                       background_y = ggh4x::elem_list_rect(fill = NA, colour = NA),
                       text_y = ggh4x::elem_list_text(
                         angle = rep(0, length(ROWS)),
                         size  = rep(5.4, length(ROWS)),
                         lineheight = rep(0.82, length(ROWS)),
                         hjust = rep(1, length(ROWS)),
                         colour = rep("black", length(ROWS))))) +
  coord_equal() +
  ggh4x::force_panelsizes(rows = unit(0.60, "in"), cols = unit(0.60, "in")) +
  theme_pub(base_size = 8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        # theme_pub() sets axis.title.x/.y explicitly, and a child set earlier beats a
        # parent set later — so `axis.title = element_blank()` alone left the raw column
        # names "xr" / "yr" on the panel. Blank the children by name. (Same trap the
        # oligo UMAP script documents at umap_theme_arrow().)
        axis.title = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.line = element_blank(),
        strip.text.x = element_text(size = 6.0),
        strip.placement = "outside", strip.clip = "off",
        panel.spacing = unit(0.02, "cm"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.width = unit(1.05, "cm"), legend.key.height = unit(0.18, "cm"),
        legend.margin = margin(t = -2, b = 0), legend.box.spacing = unit(0.06, "cm"),
        legend.text = element_text(size = 5.8), legend.title = element_text(size = 6.0),
        plot.margin = margin(1, 1, 1, 1))

BN1 <- "6_supplychain_spatial"; W1 <- 3.15; H1 <- 5.15
ggsave(file.path(PANEL, paste0(BN1, ".pdf")), p_map, width = W1, height = H1, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN1, ".png")), p_map, width = W1, height = H1, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s.{png,pdf} (%.2f x %.2f in)\n", BN1, W1, H1))

# ---------------------------------------------------------------------------
# PANEL 2 — the same arms at matched OLIGODENDROCYTE content
# ---------------------------------------------------------------------------
# Matching on deconvolved oligodendrocyte content is what stops cell loss from
# masquerading as a per-cell change: an arm that converges once content is matched
# was carried by abundance; one that stays separated is a per-cell deficit.
BR <- unique(quantile(long$oligo_prop, seq(0, 1, length.out = 11), na.rm = TRUE))
long <- long %>% mutate(bin = cut(oligo_prop, BR, include.lowest = TRUE, labels = FALSE))
sec <- long %>% filter(!is.na(bin)) %>%
  group_by(arm, Condition, sample_id, bin) %>%
  summarise(m = mean(score), .groups = "drop") %>%          # per bin, per section
  group_by(arm, Condition, sample_id) %>%
  summarise(score = mean(m), .groups = "drop")              # bins weighted equally
summ <- sec %>% group_by(arm, Condition) %>%
  summarise(m = mean(score), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = m) %>%
  mutate(diff = NHD - CON)
cat("\n== NHD - CON at matched oligodendrocyte content ==\n")
print(as.data.frame(summ %>% mutate(across(where(is.numeric), ~round(.x, 4)))))

p_q <- ggplot(sec, aes(x = Condition, y = score)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.26,
               colour = "grey35") +
  geom_point(aes(fill = Condition), shape = 21, colour = "grey25", stroke = 0.25,
             size = 1.9, position = position_jitter(width = 0.07, height = 0)) +
  scale_fill_manual(values = PAL_COND, guide = "none") +
  facet_wrap(~ arm, nrow = 1, scales = "free_y", labeller = label_wrap_gen(width = 15)) +
  labs(x = NULL, y = "Module score at matched depth\nand oligodendrocyte content") +
  theme_pub(base_size = 8) +
  theme(axis.text = element_text(size = 6.0), axis.title.y = element_text(size = 6.4),
        strip.text = element_text(size = 6.2, lineheight = 0.9),
        panel.spacing = unit(0.22, "cm"))

BN2 <- "6_supplychain_matched"; W2 <- 6.6; H2 <- 2.0
ggsave(file.path(PANEL, paste0(BN2, ".pdf")), p_q, width = W2, height = H2, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN2, ".png")), p_q, width = W2, height = H2, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s.{png,pdf} (%.2f x %.2f in)\n", BN2, W2, H2))

write.csv(bind_rows(
  sec  %>% transmute(part = "per section", arm, Condition, sample_id, value = score),
  summ %>% transmute(part = "NHD - CON",   arm, Condition = "diff",
                     sample_id = NA_character_, value = diff)),
  file.path(TDIR, "fig6_supplychain_FH.csv"), row.names = FALSE)
cat("Wrote tables/fig6_supplychain_FH.csv\n=== DONE ===\n")
