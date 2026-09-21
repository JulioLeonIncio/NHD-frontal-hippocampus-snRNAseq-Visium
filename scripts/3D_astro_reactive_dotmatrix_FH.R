#!/usr/bin/env Rscript
# =============================================================================
# 3D_astro_reactive_dotmatrix_FH.R — Compact companion to Figure 3d, frontal + hippocampus. Port of 63_fig3_astro_reactive_dotmatrix.R.
# -----------------------------------------------------------------------------
# Same load-bearing info as the CON/NHD reactive-state violin grid, a fraction of
# the footprint: a Cliff's-delta DOT-MATRIX (house 3e/3c style).
#   rows     = the 11 astrocyte-state signatures
#   columns  = the 2 regions (Frontal / Hippocampus, full names), pale PAL_REGION
#              banners + black plain text (strip_region_x)
#   dot fill = Cliff's delta (NHD - CON), house diverging steel -> white -> rose
#   dot size = -log10(BH q), p=0 underflow floored + finite cap (both logged)
#   overlay  = 3d's gated star convention verbatim (*/**/*** q<a & |d|>=0.15;
#              grey "q*" q<a & |d|<0.15; n.s. not drawn)
# RE-VIZ, not re-analysis: reads tables/fig3d_astro_signature_stats.csv so every
# number MATCHES the violin panel.  Atlas-light.
# Output: figures/Figure_3/panels/F3f_astro_reactive_dotmatrix.{pdf,png}
# -----------------------------------------------------------------------------
# Harmonised to figure 2. Fig 3 had a 1.5x internal spread in type
# (axis.text 5.6 in the LR/receptor dotplots vs 8.5 in the dotmatrix/gas6 panels) and none of
# it matched Fig 2's house spec. All Fig-3 panels now use the same tokens as Fig 2's
# 2D_secretome_MAST reference -- base_size 8 / axis.text 6.3 / axis.title 7 / strip 6.4-7 /
# legend 6.2 -- rendered at their placed width so 6.3 pt is 6.3 pt on the page.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(stringr); library(ggh4x)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("MISSING project root" = !is.na(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # REGION_FULL, strip_region_x, REGION_ORDER

PANEL <- file.path(PROJ, "figures", "Figure_3", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)

CSV <- file.path(PROJ, "tables", "fig3d_astro_signature_stats.csv")
if (!file.exists(CSV))
  stop("MISSING ", CSV, " — run 3D_astro_state_violins_FH.R first (it writes this table).")
d <- read.csv(CSV, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("signature","Region","cliff_d","q_BH","star","n_CON","n_NHD") %in% names(d)))

# A1/A2 renamed
# to provenance-anchored names; Metallothionein -> "(lost)". must match the strings
# written by 3D_astro_state_violins_FH.R (the setequal() guard below enforces it).
SIG_LEVELS <- c("Complement/IFN-reactive (Liddelow)","Pan-reactive",
                "Ischemic/S100A10-reactive (Zamanian)","DAA (Habib)",
                "Zhou NHD astrocyte","Interferon (Hasel)","STAT3 targets",
                "Metallothionein (lost)","Homeostatic astro","Heat-shock (HSF1)","Synaptogenic")
REG_LEVELS <- REGION_ORDER            # c("Frontal","Hippo")
stopifnot("stats signatures do not match expected 11-state set" =
            setequal(unique(as.character(d$signature)), SIG_LEVELS),
          setequal(unique(as.character(d$Region)), REG_LEVELS),
          nrow(d) == length(SIG_LEVELS) * length(REG_LEVELS))
cat("Loaded stats table:", nrow(d), "rows\n")
d <- d %>% mutate(signature = factor(as.character(signature), levels = SIG_LEVELS),
                  Region    = factor(as.character(Region),    levels = REG_LEVELS))

# ---- degenerate-value guards (reviewer-facing) -----------------------------
XMIN <- .Machine$double.xmin           # underflow floor
CAPL <- 50                             # ceiling on -log10(q): q = 1e-50 (stated)
n_under <- sum(!is.na(d$q_BH) & d$q_BH > 0 & d$q_BH < XMIN)
n_zero  <- sum(!is.na(d$q_BH) & d$q_BH == 0)
raw_mlp <- ifelse(is.na(d$q_BH), NA_real_, -log10(pmax(d$q_BH, XMIN)))
n_cap   <- sum(!is.na(raw_mlp) & raw_mlp > CAPL)
cat(sprintf(paste0("size guard: %d q == 0 and %d in (0, xmin) floored to xmin ",
                   "(-log10 = %.1f); %d cell(s) with -log10(q) > %d capped (q = 1e-%d).\n"),
            n_zero, n_under, -log10(XMIN), n_cap, CAPL, CAPL))
d$dotsize <- pmin(raw_mlp, CAPL)

L <- ceiling(max(abs(d$cliff_d), na.rm = TRUE) * 20) / 20
cat(sprintf("fill scale: symmetric limits +/- %.2f (max |delta| = %.3f)\n",
            L, max(abs(d$cliff_d), na.rm = TRUE)))
d$delta_c <- pmax(pmin(d$cliff_d, L), -L)

GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
# No per-nucleus */**/*** overlay (single-donor
# pseudoreplication). Only the grey effect-gated "q*" mark is drawn; the Cliff's
# delta fill carries the effect. (The upstream CSV now emits only "" or "q*".)
stars_qgrey <- d %>% filter(star == "q*")

# Hippocampus astro = single NHD lane -> grey banner +
# "(low n)" suffix, matching every other astro and oligo panel of Figure 3.
low_conf <- "Hippo"
strip_spec <- strip_region_x(REG_LEVELS, fontsize = 7.5, low_conf = low_conf, clip = "off")
ylab_map <- setNames(str_wrap(SIG_LEVELS, width = 18), SIG_LEVELS)

# Vertical layout: signatures down the
# y-axis, Region as two facet columns. Transposed, the region strips rotated onto the right
# edge and the first signature label clipped to "/IFN-reactive".
p <- ggplot(d, aes(x = 1, y = signature)) +
  geom_point(aes(fill = delta_c, size = dotsize), shape = 21, colour = "grey35", stroke = 0.3) +
  geom_text(data = stars_qgrey, aes(x = 1, y = signature, label = star), inherit.aes = FALSE,
            nudge_y = 0.28, size = 2.6, colour = "grey50", vjust = 0.5) +
  facet_grid2(cols = vars(Region), strip = strip_spec,
              labeller = labeller(Region = region_full_lowconf_labeller(low_conf))) +
  scale_fill_gradient2(low = GSEA_STEEL, mid = "white", high = GSEA_ROSE, midpoint = 0,
    limits = c(-L, L), breaks = c(-L, 0, L),
    # multi-line title and multi-line break labels collided once the legend went
    # INLINE (they overlapped each other and "+0.8" clipped). Same fix as the GSEA panel: keep
    # the title and ticks short; direction (+ = NHD-up) goes in the figure legend.
    labels = c(sprintf("%.1f", -L), "0", sprintf("+%.1f", L)),
    name = expression("Cliff's " * delta)) +
  scale_size_continuous(range = c(0.8, 5.2), limits = c(0, CAPL),
                        breaks = c(1.3, 10, 50), labels = c("0.05", "1e-10", "<=1e-50"),
                        name = "BH q") +
  # wider x-expansion widens each single-dot facet so the centred full-name banner
  # (esp. "Hippocampus (low n)") has room and does not overflow into the sibling strip.
  # Expansion 1.1 -> 0.75. It was sized for a long "(low n)" strip label that is
  # now retired, and it padded the facets far wider than their dot columns ("a lot of space on
  # the sides of the dots"). note x is continuous here, so this is the scale to touch -- adding
  # a scale_x_discrete only triggers "Scale for x is already present" and is dropped.
  scale_x_continuous(expand = expansion(add = 0.75)) +
  scale_y_discrete(limits = rev(SIG_LEVELS), labels = function(x) ylab_map[x]) +
  # Same legend idiom as 3e -- a vertical
  # stack of single rows with INLINE (left) titles and a horizontal colourbar. Stacked blocks
  # with top titles put each guide on its own baseline and read as scattered.
  # Legend at the bottom for the dot plots. Horizontal colourbar with an
  # inline (left) title, keys on one row -- same idiom as the GSEA panel, just relocated.
  # vertical colourbar again now that the legend is a right-hand column
  guides(fill = guide_colourbar(order = 1, barwidth = unit(0.28, "cm"),
                                barheight = unit(1.5, "cm"), title.position = "top"),
         # The size key sat on the
         # colourbar's "-0.8" label. Title on top + real vertical spacing between the two guides.
         size = guide_legend(order = 2, title.position = "top")) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8) +
  theme(
    axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 6.3, colour = "black", lineheight = 0.82),
    axis.line    = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.spacing.x = unit(1.5, "pt"), strip.clip = "off",
    # Legend on the right.
    legend.position = "right", legend.box = "vertical", legend.box.just = "left",
    legend.justification = "center",
    # separate the colourbar and size blocks: "+0.8" was touching "BH q".
    legend.spacing.x = unit(0.30, "cm"),
    legend.spacing.y = unit(0.40, "cm"),
    # Margin was l=-3, tuned for a right legend column; for a bottom legend the pull
    # needs to be on the TOP edge instead.
    legend.margin = margin(t = -2, r = 0, b = 0, l = 0), legend.key.size = unit(0.4, "lines"),
    legend.text = element_text(size = 6.2), legend.title = element_text(size = 6.2, lineheight = 0.9),
    plot.margin = margin(5, 4, 3, 3))

# Widened 3.4 -> 4.3 in — the "(low n)" suffix on the Hippo
# banner (clip="off") overflowed into the Frontal strip at 3.4 in. Extra width gives
# both full-name banners room so they never collide/clip.
W <- 3.95; H <- 2.90   # right-hand legend column
BN <- "F3f_astro_reactive_dotmatrix"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600, device = ragg::agg_png)
cat(sprintf("Saved %s.{pdf,png}  footprint %.1f x %.1f in\n", BN, W, H))
writeLines(capture.output(sessionInfo()),
           file.path(PANEL, "_3D_astro_reactive_dotmatrix_FH.sessionInfo.txt"))
cat("=== DONE ===\n", file = stderr())
