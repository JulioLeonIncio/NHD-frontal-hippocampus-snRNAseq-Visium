#!/usr/bin/env Rscript
# =============================================================================
# 4Z_volcano_composite_FH.R — Figure 4: Oligodendrocyte and OPC thematic volcanoes side by side in one panel.
# -----------------------------------------------------------------------------
# Assembles the two ggplot objects stashed by 4_oligolineage_volcano_MAST_FH.R —
# a native re-plot, never a raster paste-up of the two PNGs, so the type stays
# vector and the on-page point sizes stay honest.
#
# Three things the composite must do that the standalone panels do not:
#   1. one key, not two. The plotmath encoding key (nominal/strict/floors) is
#      identical for both cell types, so printing it twice wastes ~1.2 in of width
#      and invites the reader to look for a difference between them.
#   2. one theme legend. Both panels are rebuilt with VOLC_UNION_LEVELS=1 so they
#      carry the same eight theme levels; patchwork then merges the guides instead
#      of stacking two near-identical legends. Themes absent from a cell type simply
#      draw no points there — which is itself informative.
#   3. Name the halves. Without headers the reader cannot tell which pair of facets
#      is which cell type; the region strips are identical in both.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
# The stashed ggplots carry the house region labeller, which resolves REGION_FULL lazily at
# print time — so the assembler must have the theme in scope or printing fails with
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
PANEL <- file.path(PROJ, "figures", "Figure_4", "panels")

f_oli <- file.path(PROJ, "data", "4_oligo_volcano_MAST_gg.rds")
f_opc <- file.path(PROJ, "data", "4_opc_volcano_MAST_gg.rds")
for (f in c(f_oli, f_opc))
  if (!file.exists(f))
    stop("MISSING ", basename(f), " — run 4_oligolineage_volcano_MAST_FH.R with ",
         "VOLC_SAVE_GG=1 VOLC_UNION_LEVELS=1 for both VOLC_CT values first.")
oli <- readRDS(f_oli); opc <- readRDS(f_opc)

hdr <- function(g, lab) g + ggtitle(lab) +
  theme(plot.title = element_text(size = 9.5, hjust = 0.5, face = "plain",
                                  colour = "black", margin = margin(b = 2)),
        legend.position = "none")

# ---- one hand-built theme legend -------------------------------------------
# patchwork's guides="collect" does not merge these two: even with drop = FALSE the two
# panels draw keys only for the levels their own data contain, so the guides are not
# identical and patchwork keeps both — which rendered as two legends where the absent
# themes appeared as text with no swatch. Building the legend
# by hand is deterministic: every theme gets its dot whether or not it appears in a given
# cell type, which is the point of a shared key.
PAL <- c(
  "Myelin protein (up)"               = "#C1272D",
  "Stress / reactive (up)"            = "#D81B60",
  "OPC state shift (up)"              = "#E8A33D",
  "Sterol synthesis (lost)"           = "#2E86AB",
  "Fatty-acid / sphingomyelin (lost)" = "#1B9E77",
  "Lipid import (lost)"               = "#7E57A6",
  "Iron / scaffold (lost)"            = "#8C564B",
  "OPC identity / commitment (lost)"  = "#16456B")
# 4 columns x 2 rows with a non-uniform pitch. "Fatty-acid / sphingomyelin (lost)" is by
# far the longest label, so an even pitch makes it collide with column 4 — it did at both
# canvas widths tried. Column 3 therefore gets ~1.6x the gap of the others. The positions
# are in data units on a fixed 0-4.9 axis, so they scale with the canvas and this holds if
# the panel is resized again.
leg_df <- data.frame(
  lab = names(PAL), col = unname(PAL),
  x   = rep(c(0, 1.15, 2.30, 3.95), each = 2) + 0.02,
  y   = rep(c(1, 0), times = 4))
legend_plot <- ggplot(leg_df, aes(x, y)) +
  geom_point(aes(colour = I(col)), size = 2.6) +
  geom_text(aes(label = lab), hjust = 0, nudge_x = 0.075, size = 2.35) +
  scale_x_continuous(limits = c(-0.05, 5.55)) +
  scale_y_continuous(limits = c(-0.5, 1.5)) +
  theme_void() + theme(plot.margin = margin(2, 2, 2, 2))

# The encoding key is a tall narrow block with absolute internal coordinates, so any
# short-and-wide slot squashes its lines into each other (it did, twice). Give it its own
# full-height column on the right — the shape it was designed for in the standalone panels.
left  <- (hdr(oli$p, "Oligodendrocytes") | hdr(opc$p, "OPC")) / legend_plot +
  plot_layout(heights = c(1, 0.30))
combo <- left | oli$key
combo <- combo + plot_layout(widths = c(1, 0.20))

BN <- "F4b_oligolineage_volcano_composite"
# Taller: at 3.15 in the -log10 axis was compressed and the point
# cloud read as a flat band. Height only — the width is already at page.
W  <- 8.10; H <- 3.70
ragg::agg_png(file.path(PANEL, paste0(BN, ".png")), width = W, height = H,
              units = "in", res = 600)
print(combo); invisible(dev.off())
ggsave(file.path(PANEL, paste0(BN, ".pdf")), combo, width = W, height = H,
       useDingbats = FALSE)
cat(sprintf("Wrote %s.{png,pdf}  (%.2f x %.2f in)\n", BN, W, H))
