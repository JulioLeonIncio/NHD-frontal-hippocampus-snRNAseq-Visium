#!/usr/bin/env Rscript
# =============================================================================
# 108b_visium_HE_panel_FH.R — Fig. 6b rebuilt as two rows: H&E of each capture area (top) over the spatial-domain map (bottom), in identical frames.
# -----------------------------------------------------------------------------
# The domain row is script 74's panel b, reproduced line for line (same spot coordinates, same
# per-section normalisation, same palette, same panel sizes). The H&E row places the crop made
# by 108a (already transposed into the Fig-6b frame) with annotation_raster at the crop's map
# extent, so it is registered to the spots by construction — and the registration is checked:
# figures/_diagnostics/F6b_HE_registration.png overlays the spot centroids on the H&E.
# No scale bar on the H&E (submitted without; can be added at revision if asked); the
# 1-mm length in map units is still computed from spot_diameter_fullres (55 µm) via HE_frames.json
# and printed, so a bar can be added in one line.
#
# Output: figures/Figure_6/panels/F6b_visium_HE_domains.{png,pdf}   (4.35 in wide, two rows)
#         figures/_diagnostics/F6b_HE_registration.png  (WM + Vasc/immune spots, large, over the H&E)
# The single-row F6b_visium_domain_spatial stays on disk; the composite relinks 6b to this file
# and the page re-flows (c-f move down by the added row height).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(png); library(jsonlite); library(ragg); library(patchwork) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
ROOT <- dirname(PROJ); C2L <- file.path(ROOT, "Visium", "cell2location", "c2l_MAIN")
HEP  <- file.path(ROOT, "Visium", "_HE_fullres", "panels")
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels"); DIAG <- file.path(PROJ, "figures", "_diagnostics")
say <- function(...) cat(sprintf(...), "\n")
SEC <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))   # DOM_MAP / DOM_LEV / DOM_PAL from integrated_harmony/*.csv (as script 74)

# ---- spots exactly as script 74 ------------------------------------------------
meta <- read.csv(file.path(C2L, "visium_meta.csv"), stringsAsFactors = FALSE); names(meta)[1] <- "spot_id"
xy   <- read.csv(file.path(C2L, "spot_coords.csv"), stringsAsFactors = FALSE)
meta <- meta %>% mutate(domain = factor(unname(DOM_MAP[as.character(seurat_clusters)]), levels = DOM_LEV),
                        sample_id = factor(sample_id, levels = SEC)) %>%
  left_join(xy %>% transmute(spot_id = cell, x, y), by = "spot_id") %>% filter(!is.na(domain), !is.na(x))
norm <- meta %>% group_by(sample_id) %>% summarise(sp = max(diff(range(x)), diff(range(y))), cx = mean(range(x)), cy = mean(range(y)), .groups = "drop")
sp <- meta %>% left_join(norm, by = "sample_id") %>% mutate(xr = (x - cx) / sp, yr = (y - cy) / sp)
YSPAN <- max(tapply(sp$yr, sp$sample_id, function(v) diff(range(v))))
XSPAN <- max(tapply(sp$xr, sp$sample_id, function(v) diff(range(v))))

# ---- H&E rasters in the same normalised frame ------------------------------------
fr <- fromJSON(file.path(HEP, "HE_frames.json"))
he <- lapply(SEC, function(s) { f <- fr[[s]]; n <- norm[norm$sample_id == s, ]
  img <- readPNG(file.path(HEP, paste0(s, "_HE_fig6b_frame.png")))
  list(sample_id = s, img = img, xmin = (f$x_min - n$cx) / n$sp, xmax = (f$x_max - n$cx) / n$sp,
       ymin = (f$y_min - n$cy) / n$sp, ymax = (f$y_max - n$cy) / n$sp, bar = f$scalebar_1mm_map_units / n$sp) })
names(he) <- SEC
lim <- sp %>% group_by(sample_id) %>% summarise(x0 = min(xr), x1 = max(xr), y0 = min(yr), y1 = max(yr), .groups = "drop")

base_theme <- theme_pub(base_size = 8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.line = element_blank(),
        strip.text.x = element_text(size = 6.8, colour = "black"), strip.background = element_blank(),
        panel.spacing = unit(0.04, "cm"), plot.margin = margin(1, 1, 1, 1))

# annotation_raster is not facet-aware: one GeomCustomAnn layer per facet, keyed by sample_id
raster_layer <- function(h) layer(data = data.frame(sample_id = factor(h$sample_id, levels = SEC), x = 0, y = 0), mapping = aes(x = x, y = y),
  stat = "identity", geom = ggplot2::GeomCustomAnn, position = "identity", inherit.aes = FALSE,
  params = list(grob = grid::rasterGrob(h$img, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc")),
                xmin = h$xmin, xmax = h$xmax, ymin = h$ymin, ymax = h$ymax))
p_he <- ggplot(sp, aes(xr, yr)) + geom_blank() +
  facet_wrap(~ sample_id, nrow = 1, labeller = labeller(sample_id = function(x) sub("_Frontal", " ", x)))
for (s in SEC) p_he <- p_he + raster_layer(he[[s]])
p_he <- p_he + coord_equal() + ggh4x::force_panelsizes(rows = unit(0.95 * YSPAN, "in"), cols = unit(0.95, "in")) + base_theme + labs(x = NULL, y = NULL)

p_dom <- ggplot(sp, aes(xr, yr, colour = domain)) + geom_point(size = 0.6, stroke = 0) +   # closes the hex grid at 0.95 in
  scale_colour_manual(values = DOM_PAL, name = NULL, guide = guide_legend(override.aes = list(size = 1.9), nrow = 1)) +
  facet_wrap(~ sample_id, nrow = 1) + coord_equal() +
  ggh4x::force_panelsizes(rows = unit(0.95 * YSPAN, "in"), cols = unit(0.95, "in")) + base_theme + labs(x = NULL, y = NULL) +
  theme(strip.text.x = element_blank(), legend.position = "bottom", legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 6.4, colour = "black"), legend.margin = margin(t = -2), legend.box.spacing = unit(0.06, "cm"))

p <- p_he / p_dom
# canvas measured from the assembled gtable (fixed panels + strips + legend) plus a small margin —
# no blank band above, between or below the rows
g_b <- patchwork::patchworkGrob(p)
W <- sum(grid::convertWidth(g_b$widths, "in", valueOnly = TRUE)) + 0.06
H <- sum(grid::convertHeight(g_b$heights, "in", valueOnly = TRUE)) + 0.06
BN <- "F6b_visium_HE_domains"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600, device = ragg::agg_png)
say("wrote %s (%.2f x %.2f in, measured canvas); 1-mm = %s map units", BN, W, H, paste(round(sapply(he, `[[`, "bar"), 3), collapse = " / "))

# registration check: spot centroids over the H&E, coloured by domain, thin points
p_chk <- p_he + geom_point(data = sp[sp$domain %in% c(DOM_WHITE, "Vasc/immune"), ], aes(xr, yr, colour = domain), size = 0.9, alpha = 0.75, stroke = 0, inherit.aes = FALSE) +
  scale_colour_manual(values = DOM_PAL, guide = "none")
ggsave(file.path(DIAG, "F6b_HE_registration.png"), p_chk, width = 2 * W, height = 2 * (0.95 * YSPAN + 0.4), dpi = 300, device = ragg::agg_png)
say("registration check: figures/_diagnostics/F6b_HE_registration.png")
cat("\n=== DONE ===\n", file = stderr())
