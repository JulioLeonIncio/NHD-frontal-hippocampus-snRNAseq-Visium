#!/usr/bin/env Rscript
# =============================================================================
# 06d_figure1_umap_compact_FH.R — Figure 1 panel b (FH rebuild) — compact atlas UMAP, SCZ-resubmission style.
# -----------------------------------------------------------------------------
# Faithful port of manuscript_7fig/scripts/06d_figure1_umap_compact.R to the
# Frontal + Hippocampus atlas (NHD_FH_harmony.rds). Three equal small UMAPs
# (cell type / condition / region): near-fixed-size points (abundant pops look
# abundant, not density-shrunk), blank axes, legends at bottom.
# OCC does not exist in this atlas — only Frontal + Hippocampus are shown.
#
# Builds/uses a UMAP+metadata cache (umap_df.rds) so re-runs never reload the
# 2.9 GB Seurat object. Cache auto-rebuilds if older than the atlas (mtime guard).
# Output: figures/Figure_1/panels/F1b_atlas_umap.{png,pdf}
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tibble); library(patchwork); library(ggrepel)
})
set.seed(42)

# --- portable project root (julio.l origin/RIKEN or JulioLeon local) --------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
if (is.na(PROJ)) stop("MISSING project root — neither julio.l nor JulioLeon path exists")
DIR   <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
ATLAS <- file.path(DIR, "atlas", "NHD_FH_harmony.rds")
CACHE <- file.path(DIR, "data", "_cache_Fig1"); dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
PANEL <- file.path(DIR, "figures", "Figure_1", "panels"); dir.create(PANEL, recursive = TRUE, showWarnings = FALSE)
THEME <- file.path(DIR, "scripts", "22_publication_theme_FH.R")
if (!file.exists(ATLAS)) stop(sprintf("MISSING %s — atlas build failed", ATLAS))
if (!file.exists(THEME)) stop(sprintf("MISSING %s — FH theme missing", THEME))
source(THEME)   # PAL_CELLTYPE / PAL_COND / PAL_REGION_UMAP / theme_pub

# --- UMAP + metadata cache (mtime-guarded so a fresh atlas invalidates it) --
umap_cache <- file.path(CACHE, "umap_df.rds")
stale <- !file.exists(umap_cache) ||
         file.info(umap_cache)$mtime < file.info(ATLAS)$mtime
if (stale) {
  cat("== cache miss/stale: extracting UMAP + metadata from atlas ==\n")
  suppressPackageStartupMessages(library(Seurat))
  obj <- readRDS(ATLAS)
  umap_key <- if ("umap" %in% names(obj@reductions)) "umap" else "UMAP"
  umap_df <- as.data.frame(Embeddings(obj, reduction = umap_key)) %>%
    setNames(c("UMAP_1","UMAP_2")) %>%
    rownames_to_column("cell") %>%
    left_join(obj@meta.data %>% rownames_to_column("cell"), by = "cell") %>%
    mutate(cell_type = as.character(new_annotation),
           Condition = factor(Condition, levels = c("CON","NHD")),
           Region    = factor(Region,    levels = c("Frontal","Hippo")))
  saveRDS(umap_df, umap_cache)
  rm(obj); gc()
  cat(sprintf("   cached %d nuclei -> %s\n", nrow(umap_df), umap_cache))
} else {
  cat("== loading cached UMAP + metadata ==\n")
  umap_df <- readRDS(umap_cache)
}
u <- umap_df
stopifnot(all(c("UMAP_1","UMAP_2","cell_type","Condition","Region") %in% colnames(u)))

# Canonical 8-type overview (drop low-count Lymphocyte / Micro-PVM_doublet, as in
# the reference — those minor categories are reported in the QC supplement).
CT_LEVELS <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh","Endo","Pericytes")
n0 <- nrow(u)
u <- u[u$cell_type %in% CT_LEVELS, ]
cat(sprintf("cell-type filter: %d -> %d nuclei (dropped %d non-canonical)\n",
            n0, nrow(u), n0 - nrow(u)))
u$cell_type <- factor(u$cell_type, levels = CT_LEVELS)
u$Condition <- factor(u$Condition, levels = c("CON","NHD"))
u$Region    <- factor(u$Region,    levels = c("Frontal","Hippo"))
cat(sprintf("region split: Frontal=%d  Hippo=%d\n",
            sum(u$Region == "Frontal"), sum(u$Region == "Hippo")))

# panel sizing (matches the reference SCZ-style compact layout):
#   No aspect.ratio lock (that left side gaps); shape set purely by canvas dims.
PW <- 8.4 * 0.8 * 0.9 * 0.85   # -> 5.14 in
PH <- 2.64 * 0.8               # -> 2.11 in

# --- geometry for the shared corner-arrow axis glyph ------------------------
# the glyph used to be
# anchored inside the data range at its bottom-right, on the assumption that corner
# was empty. It is not — the Oligo cloud fills it, so the arrows and the "UMAP 2"
# label were drawn straight across the blob.
# The glyph is now placed outside the data entirely: the panel limits are padded to
# the right and below, and the L sits in that padding (gx > max(UMAP_1)), so no
# amount of re-clustering can ever make it collide with a cell again.
# The same limits are applied to all three UMAPs so the three copies of the
# embedding stay at identical scale (they are compared side by side).
# The glyph lives in a padding strip added below the data, anchored at the left.
# Bottom-left (not bottom-right) because the horizontal "UMAP 1" label is ~2x
# longer than the arrow arm at this panel size, so it has to run rightward into
# open strip; anchored at the right edge it would overflow the panel.
xr <- range(u$UMAP_1); yr <- range(u$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
PAD_B <- 0.20                                      # bottom strip that houses the glyph
XLIM <- c(xr[1] - 0.11 * dx, xr[2] + 0.02 * dx)   # left pad clears the rotated y-label
YLIM <- c(yr[1] - PAD_B * dy, yr[2] + 0.02 * dy)
gx <- xr[1] + 0.005 * dx                           # glyph origin: left edge ...
gy <- yr[1] - 0.135 * dy                           # ... and below all data
armx <- 0.11 * dx; army <- 0.11 * dy               # vertical arm stops below yr[1]

# Compact SCZ-style UMAP: near-fixed point size (abundant pops look abundant),
# blank axes, bottom legend, seeded row shuffle so no level buries another.
compact_umap <- function(df, colvar, pal, legend_nrow = 1, lead = NULL,
                         pt_alpha = 0.8, pt_size = 0.22, shuffle_seed = NULL,
                         legend_labels = ggplot2::waiver()) {
  if (!is.null(shuffle_seed)) set.seed(shuffle_seed)
  if (!is.null(lead)) {
    df <- bind_rows(df[df[[colvar]] == lead, ],
                    df[df[[colvar]] != lead, ][sample(sum(df[[colvar]] != lead)), ])
  } else df <- df[sample(nrow(df)), ]
  # Dots are drawn as solid points (shape 16), not
  # shape 21 with a grey20 ring. At this point size (0.22) the ring covered a large
  # fraction of each dot and mixed grey into every colour — side-by-side renders
  # showed it was the single biggest thing muting the palette. Removing it is what
  # makes the clusters read at journal saturation. The ring standard still applies
  # to LARGER dot marks elsewhere (effect-size dots etc.), not to dense UMAP points.
  ggplot(df, aes(UMAP_1, UMAP_2, colour = .data[[colvar]])) +
    geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
    scale_colour_manual(values = pal, name = NULL, labels = legend_labels) +
    coord_cartesian(xlim = XLIM, ylim = YLIM, clip = "off") +
    theme_pub(base_size = 8) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank(),
          legend.position = "bottom",
          legend.key.size = unit(0.22, "cm"), legend.text = element_text(size = 6.2),
          legend.box.spacing = unit(1, "pt"), legend.margin = margin(-2, 0, 0, 0),
          plot.margin = margin(1, 3, 1, 2)) +
    guides(colour = guide_legend(nrow = legend_nrow, byrow = TRUE,
                                 override.aes = list(shape = 16, size = 1.9, alpha = 1)))
}

# single corner L-arrow axis glyph (cell-type panel only) — one orientation cue
# for the whole row instead of redundant per-panel axis titles.
corner_axis <- list(
  annotate("segment", x = gx, xend = gx + armx, y = gy, yend = gy,
           arrow = arrow(length = unit(1.3, "mm"), type = "closed"),
           linewidth = 0.3, colour = "grey25"),
  annotate("segment", x = gx, xend = gx, y = gy, yend = gy + army,
           arrow = arrow(length = unit(1.3, "mm"), type = "closed"),
           linewidth = 0.3, colour = "grey25"),
  # Labels are anchored at the corner (hjust = 0) so they run away from it —
  # rightward under the x-arm, upward beside the y-arm. Centring them on the arms
  # (the old hjust = 0.5) made them overhang and cross the arrows, because the
  # words are wider than the arms are long.
  annotate("text", x = gx + 0.012 * dx, y = gy - 0.012 * dy,
           label = "UMAP 1", size = 1.9, colour = "grey25", hjust = 0, vjust = 1),
  annotate("text", x = gx - 0.055 * dx, y = gy + 0.012 * dy,
           label = "UMAP 2", size = 1.9, colour = "grey25", hjust = 0,
           vjust = 0.5, angle = 90))

# cell-type panel: in-plot centroid labels (grey leader lines) instead of legend.
# NB house rule: no bold anywhere -> centroid labels are PLAIN (the reference used
# fontface="bold"; we use "plain" per the FH house standard). These are data
# labels, not caption text.
# They
# are small smears that collide with Astro in the lower-left; they remain as
# coloured points and are named in the figure legend. Repel force is boosted so the
# remaining lower-left labels (Astro / Micro-PVM) separate cleanly.
LABEL_DROP <- c("Endo","Pericytes")
cent <- u %>% group_by(cell_type) %>%
  summarise(x = median(UMAP_1), y = median(UMAP_2), .groups = "drop") %>%
  filter(!cell_type %in% LABEL_DROP) %>%
  # house rule: no code tokens in blessed labels -> display recode (Neuron_Ex ->
  # "Neuron Ex", etc.). The coded `cell_type` is kept for the join above; only the
  # printed `disp` string reaches geom_text_repel.
  mutate(disp = ct_disp(cell_type))
cat(sprintf("centroid labels: %d types shown (dropped %s)\n",
            nrow(cent), paste(LABEL_DROP, collapse = ", ")))
p_ct <- compact_umap(u, "cell_type", PAL_CELLTYPE, legend_nrow = 3, lead = "Oligo") +
          corner_axis +
          geom_text_repel(data = cent, aes(x, y, label = disp), inherit.aes = FALSE,
                          # Micro-PVM's centroid label drifts right until it touches
                          # the Oligo cloud; seat it in the gap above its own cluster.
                          nudge_x = ifelse(cent$cell_type == "Micro-PVM", -0.05 * dx, 0),
                          nudge_y = ifelse(cent$cell_type == "Micro-PVM",  0.05 * dy, 0),
                          size = 2.2, fontface = "plain", colour = "grey10",
                          bg.color = "white", bg.r = 0.2, box.padding = 0.32, force = 6,
                          force_pull = 0.6, max.overlaps = 20, seed = 42,
                          # Min.segment.length was 0, which
                          # forced a leader line on every label — including the ones
                          # sitting on their own cluster, so each label had a short
                          # grey stub scratched across its cells. A leader is now
                          # drawn only when the label is genuinely displaced.
                          segment.color = "grey55",
                          min.segment.length = unit(12, "pt"),
                          segment.size = 0.25,
                          # keep labels inside the data area, out of the axis-glyph padding
                          xlim = c(xr[1], xr[2]), ylim = c(yr[1], yr[2])) +
          theme(legend.position = "none")

# Condition copy: blue/red kept lighter than the saturated PAL_COND so shuffled
# Paled 10% via the theme's .pale() in the same pass, to match the other palettes.
PAL_COND_PALE <- setNames(.pale(c(CON = "#3B8FD4", NHD = "#E8503C"), 0.10), c("CON","NHD"))
p_cond <- compact_umap(u, "Condition", PAL_COND_PALE, legend_nrow = 1,
                       pt_alpha = 0.70, shuffle_seed = 43)

# Region copy: DARKER region tones (PAL_REGION_UMAP, 0.32 blend) so points read as
# scattered dots; analysis-panel facet banners keep PAL_REGION_PALE (0.55).
# Plain seeded shuffle so neither region is drawn wholly on top.
# via REGION_FULL, matching panels d/e/f. Data/factor levels stay short codes.
p_reg <- compact_umap(u, "Region", PAL_REGION_UMAP, legend_nrow = 1,
                      pt_alpha = 0.70, shuffle_seed = 42,
                      legend_labels = REGION_FULL[levels(u$Region)])

pB <- p_ct + p_cond + p_reg + plot_layout(nrow = 1)

ggsave(file.path(PANEL, "F1b_atlas_umap.png"), pB, width = PW, height = PH,
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "F1b_atlas_umap.pdf"), pB, width = PW, height = PH, bg = "white")
cat(sprintf("  wrote F1b_atlas_umap.{png,pdf}  (%.2f x %.2f in)\n", PW, PH))

# provenance
si <- file.path(DIR, "logs", "1B_umap_compact_FH_sessionInfo.txt")
dir.create(dirname(si), recursive = TRUE, showWarnings = FALSE)
writeLines(c(sprintf("run: %s", format(Sys.time())), sprintf("PROJ: %s", PROJ),
             capture.output(sessionInfo())), si)
cat("\n=== DONE: F1b_atlas_umap (FH) ===\n")
