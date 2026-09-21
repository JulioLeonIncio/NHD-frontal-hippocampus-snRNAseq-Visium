# =============================================================================
# 22_publication_theme.R — Publication theme + palettes (used by other scripts via source())
# -----------------------------------------------------------------------------
# Aesthetic inspired by recent Nature/Cell single-cell figures:
# Clean white background, thin black axes, no gridlines, big bold
# titles, italic gene names, journal-style categorical palette
# (teal, magenta, steel blue, salmon, mustard, light grey).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

library(ggplot2)

# --- Color palettes (matching the example figure) -------------------
# Disease comparison
PAL_DGE <- c(
  "Up in NHD"  = "#C0392B",   # warm red
  "Up in CON"  = "#2C7FB8",   # cool blue
  "ns"         = "#BDBDBD",   # neutral grey
  "Highlight"  = "#E8B83A",   # mustard amber
  "Notable"    = "#4DBDA1"    # teal
)

# Cell types — mapped to the journal palette in the example image
PAL_CELLTYPE <- c(
  "Micro-PVM"  = "#B83C7E",   # magenta — microglia stand out
  "Astro"      = "#4DBDA1",   # teal
  "Oligo"      = "#5985C5",   # steel blue
  "OPC"        = "#E6817D",   # salmon
  "Neuron_Ex"  = "#E8B83A",   # mustard
  "Neuron_Inh" = "#C8A2C8",   # light orchid
  "Endo"       = "#7F7F7F",   # neutral
  "Pericytes"  = "#C4C4C4"    # light grey
)

# Conditions
PAL_COND <- c("CON" = "#2C7FB8", "NHD" = "#C0392B")

# Regions — neuroanatomical mnemonic palette
# Frontal = dusty aubergine (prefrontal "higher-order" violet)
# OCC     = caramel / burnt bronze (visual cortex warmth)
# Hippo   = moss green (limbic earth-tone)
# Chosen to be distinct from PAL_CELLTYPE (Astro teal, Neuron_Ex mustard, OPC
# salmon) and PAL_COND (CON blue, NHD red) — no cross-palette collisions.
PAL_REGION <- c(
  "Frontal" = "#6B4E7D",   # dusty aubergine
  "OCC"     = "#A66B2E",   # caramel / burnt bronze
  "Hippo"   = "#4F7942"    # moss green
)

# Pale region tones — PAL_REGION blended 55% toward white. These are the same
# tones as the region facet-strip BANNERS (strip_region_x below) and, from
# Defined via the same blend used by strip_region_x (not hardcoded) so the two
# stay byte-identical: Frontal #BCAFC4, OCC #D6BCA0, Hippo #AFC2A9.
# Saturated PAL_REGION is unchanged — it still feeds the banners (via this
# blend) and any other saturated use (heatmap column strips, etc.).
.region_lighten <- function(hex, f = 0.55) {
  m <- grDevices::col2rgb(hex)
  grDevices::rgb(t(m * (1 - f) + 255 * f), maxColorValue = 255)
}
PAL_REGION_PALE <- setNames(.region_lighten(PAL_REGION, 0.55), names(PAL_REGION))

# DARKER region tones for UMAP points only: the 0.55 pale
# tones were too washed out to see as scattered points, so region UMAP points
# use a 0.32 blend toward white (deeper mauve/tan/sage, still soft). The facet
# strip banners keep PAL_REGION_PALE (0.55) so their black text stays legible.
# Derived via the same .region_lighten() blend (not hand-typed) so it tracks
# PAL_REGION: Frontal #9A86A6, OCC #C29A70, Hippo #87A37E.
PAL_REGION_UMAP <- setNames(.region_lighten(PAL_REGION, 0.32), names(PAL_REGION))

# Full region display names — use everywhere a region is
# LABELLED (strip banners, axes, legends). Data/factor levels stay the short
# codes; this maps code -> printed name. Pass as `labeller = REGION_LABELLER`
# to facet_*(), or `labels = REGION_FULL[levels]` to scale_*_discrete().
REGION_FULL <- c(Frontal = "Frontal", OCC = "Occipital", Hippo = "Hippocampus")
REGION_LABELLER <- ggplot2::as_labeller(REGION_FULL)

# Helper: coloured facet strip for region-faceted panels (ggh4x).
# Returns a `strip_themed` object you pass to facet_grid2 / facet_wrap2
# so the Region strip header is filled with PAL_REGION and text is white.
#   Usage:
#     library(ggh4x)
#     p + facet_grid2(sig ~ Region, strip = strip_region_x(c("Frontal","Hippo")))
#   low_conf: optional character vector of region names that are statistically
#     low-confidence (e.g. n < 150 nuclei/group). Those strips are overridden to
#     a desaturated grey background with dark text so a whole-region power caveat
#     reads at the facet-strip level (does not touch per-point aesthetics).
strip_region_x <- function(regions = c("Frontal","OCC","Hippo"),
                           fontsize = 7, text_colour = "black",   # black labels on the colour banners
                           low_conf = character(0),
                           low_conf_fill = "grey80",
                           low_conf_text = "grey25",
                           # clip: "inherit" (default) keeps prior behaviour for every
                           # existing caller (incl. Figs 1-2). Pass clip="off" for panels
                           # whose facets are so narrow the banner would truncate the
                           # (full) region name, letting the centred label overflow the
                           # pale banner instead of being cut (theme(strip.clip=) does
                           # not reach ggh4x strip_themed).
                           clip = "inherit") {
  if (!requireNamespace("ggh4x", quietly = TRUE))
    stop("strip_region_x() needs the ggh4x package.")
  # Pale region banners: PAL_REGION blended ~55% toward white so
  # black strip labels read clearly (the full-saturation banners washed them out).
  # Uses the shared PAL_REGION_PALE constant (defined above from the identical
  # blend) so banners stay byte-identical and match the region UMAP points.
  fills <- unname(PAL_REGION_PALE[regions])
  txt   <- rep(text_colour, length(regions))
  is_lc <- regions %in% low_conf
  fills[is_lc] <- low_conf_fill
  txt[is_lc]   <- low_conf_text
  ggh4x::strip_themed(
    clip         = clip,
    background_x = ggh4x::elem_list_rect(fill = fills, colour = NA),
    text_x       = ggh4x::elem_list_text(colour = txt,
                                         size   = fontsize,
                                         face   = "plain")
  )
}

# Region facet labeller that appends a plain-ASCII low-confidence marker
# " (low n)" to the named low-confidence regions and leaves the others untouched.
# ASCII (not a † dagger): the base pdf() device cannot encode U+2020 (mbcsToSbcs
# conversion failure), and these panels save both .png (ragg) and .pdf (base).
LOWCONF_MARK <- " (low n)"
region_lowconf_labeller <- function(low_conf = character(0),
                                    mark = LOWCONF_MARK) {
  ggplot2::as_labeller(function(x)
    ifelse(x %in% low_conf, paste0(x, mark), x))
}
# Standard one-line footnote for any panel carrying >=1 low-confidence region.
LOWCONF_CAPTION <- "(low n) = n < 150 nuclei/group, lower confidence"

# --- Theme function -------------------------------------------------
# Bigger labels, no gridlines, thin black axis lines (matches example)
theme_pub <- function(base_size = 8, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.caption     = element_text(size = base_size - 1, color = "grey45",
                                      hjust = 0, margin = margin(t = 4)),
      axis.title       = element_text(face = "plain", size = base_size + 1,
                                      color = "black"),
      axis.title.x     = element_text(margin = margin(t = 3)),
      axis.title.y     = element_text(margin = margin(r = 3)),
      axis.text        = element_text(color = "black", size = base_size),
      axis.line        = element_line(color = "black", linewidth = 0.35),
      axis.ticks       = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = unit(0.10, "cm"),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text       = element_text(face = "plain", size = base_size,
                                      color = "black",
                                      margin = margin(b = 2, t = 2)),
      legend.position  = "right",
      legend.title     = element_text(face = "plain", size = base_size),
      legend.text      = element_text(size = base_size),
      legend.key.size  = unit(0.35, "cm"),
      legend.spacing.x = unit(0.2, "cm"),
      legend.margin    = margin(0, 0, 0, 0),
      panel.grid       = element_blank(),
      plot.margin      = margin(4, 4, 4, 4),
      plot.background  = element_rect(fill = "white", color = NA),
      plot.tag         = element_text(face = "bold", size = base_size + 2,
                                      color = "black")
    )
}
theme_pub_categorical <- theme_pub  # alias
