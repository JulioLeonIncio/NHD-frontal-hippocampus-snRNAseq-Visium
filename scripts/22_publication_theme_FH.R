# =============================================================================
# 22_publication_theme_FH.R — Publication theme + palettes — frontal + hippocampus rebuild (used by other scripts via source())
# -----------------------------------------------------------------------------
# Faithful adaptation of scripts/22_publication_theme.R with OCC dropped
# everywhere (this rebuild is Frontal + Hippocampus only, 2 regions).
# The only changes vs the 3-region original are region-scoped:
#   * PAL_REGION / PAL_REGION_PALE / PAL_REGION_UMAP drop OCC
#   * REGION_FULL / REGION_LABELLER drop OCC
#   * strip_region_x() default `regions` = c("Frontal","Hippo")
# Everything else (the lighten-0.55 pale-banner blend, lighten-0.32 UMAP
# tones, black strip text on pale fill, theme_pub(), low-conf machinery,
# PAL_CELLTYPE/PAL_COND/PAL_DGE) is identical to the original.
#
# Aesthetic inspired by recent Nature/Cell single-cell figures:
# Clean white background, thin black axes, no gridlines, italic gene
# names, journal-style categorical palette (teal, magenta, steel blue,
# salmon, mustard, light grey).
#
# House RULES (hard constraints, from the pi):
#   * no bold fonts anywhere in any panel (default all text plain).
#   * no caption text inside panels (all explanation goes to legends).
#   * Region facet strips = pale fill + black text, never saturated + white.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

library(ggplot2)

# Blend a colour toward white by fraction f (0 = unchanged, 1 = white).
# Single source for every "pale" tone in this theme so the derived palettes can
# never drift from their base hues.
.pale <- function(hex, f = 0.10) {
  m <- grDevices::col2rgb(hex)
  grDevices::rgb(t(m * (1 - f) + 255 * f), maxColorValue = 255)
}

# --- Color palettes (matching the example figure) -------------------
# Disease comparison
PAL_DGE <- c(
  "Up in NHD"  = "#C0392B",   # warm red
  "Up in CON"  = "#2C7FB8",   # cool blue
  "ns"         = "#BDBDBD",   # neutral grey
  "Highlight"  = "#E8B83A",   # mustard amber
  "Notable"    = "#4DBDA1"    # teal
)

# Cell types — vibrant journal palette ("pick up more vibrant
# tones ... so they look like UMAPs from Science or Nature"). Each hue keeps its
# established IDENTITY (microglia magenta, astro teal, oligo blue, OPC coral,
# Neuron_Ex gold, Neuron_Inh purple) — only saturation/value were lifted, so
# cross-figure colour meaning is unchanged. Compared by eye against the previous
# muted set and a ggsci-NPG-style alternative; the NPG navy made the 46%-of-atlas
# Oligo blob read as a heavy dark mass, so this set was chosen.
# "nice colors, we just need to pale them 10%" —
# so the vivid set below is the BASE and every consumer gets it blended 10% toward
# white. Derived via .pale() rather than hand-typed hex so base and shipped tone
# can never drift apart.
PAL_CELLTYPE_VIVID <- c(
  "Micro-PVM"  = "#D62A86",   # magenta — microglia stand out
  "Astro"      = "#00B294",   # teal
  "Oligo"      = "#3B7DDD",   # blue
  "OPC"        = "#F4623A",   # coral
  "Neuron_Ex"  = "#F2B01E",   # gold
  "Neuron_Inh" = "#9C5FD0",   # purple
  "Endo"       = "#6E6E6E",   # neutral
  "Pericytes"  = "#B6B6B6"    # light grey
)
PAL_CELLTYPE <- setNames(.pale(PAL_CELLTYPE_VIVID, 0.10), names(PAL_CELLTYPE_VIVID))

# Cell-type display names (no code tokens / underscores in any
# blessed panel label). Data/factor levels & join keys stay the coded symbols
# ("Neuron_Ex"); this maps code -> printed label. Use everywhere a cell type is
# LABELLED (facet strips, axis ticks, on-plot annotations, legends): pass as
# `labels = CT_DISPLAY[levels]` to scale_*_discrete(), `labeller = CT_LABELLER`
# to facet_*(), or `ct_disp(x)` for a plain character vector. Only Neuron_Ex /
# Neuron_Inh carry underscores; the rest are identity so this is safe to apply
# blanket across all cell types.
CT_DISPLAY <- c(
  "Micro-PVM"  = "Micro-PVM",
  "Astro"      = "Astro",
  "Oligo"      = "Oligo",
  "OPC"        = "OPC",
  "Neuron_Ex"  = "Neuron Ex",
  "Neuron_Inh" = "Neuron Inh",
  "Endo"       = "Endo",
  "Pericytes"  = "Pericytes"
)
# vectorised recode: any symbol not in the map passes through unchanged.
ct_disp <- function(x) ifelse(as.character(x) %in% names(CT_DISPLAY),
                              unname(CT_DISPLAY[as.character(x)]), as.character(x))
CT_LABELLER <- ggplot2::as_labeller(ct_disp)

# Conditions
# immune-compartment display labels — single source.
# The myeloid partition is not independently identifiable as perivascular macrophage
# (re-derived without CD163/F13A1 it does not resolve at all; Jaccard 0.012), so the
# population must be named by its markers, not by an assumed cell identity. Data keys
# stay "PVM" for continuity with the assignment CSVs; only the printed label changes.
IMM_LABEL <- c("Microglia"    = "Microglia",
               "PVM"          = "CD163+/F13A1+",
               "T-lymphocyte" = "T-lymphocyte")
IMM_LABELLER <- ggplot2::as_labeller(IMM_LABEL)

PAL_COND <- c("CON" = "#2C7FB8", "NHD" = "#C0392B")
# Role tint — single source. CON/NHD may vary in
# lightness by panel role but never in hue. Figure 2 had accumulated three
# incompatible pairs (#82B2D6/#D6837A in the violins/UMAP, #3E7CB1/#D1495B in the
# DAM/secretome bars, PAL_COND elsewhere) purely because each was hand-typed
# locally; #D1495B is a different hue (rose vs brick), which is what made the drift
# visible side by side. Panels needing a lighter fill (violins, density ribbons,
# dense UMAP points) use PAL_COND_PALE; bars/points/slopes use PAL_COND directly.
PAL_COND_PALE <- setNames(.pale(PAL_COND, 0.45), names(PAL_COND))

# Cross-cohort (Zhou 2023) concordance classes — single source. The four Zhou scatters each hand-typed green (#1B7837) / red (#B2182B),
# which is exactly the red-green contrast Nature Portfolio asks authors to recolour
# ("avoid the use of red and green for contrast"). Okabe-Ito ink black vs orange: legible
# to every common colour-vision deficiency, and neither hue is CON/NHD, a region or a
# cell type, so the pair carries no borrowed meaning.
# SHAPE_ZHOU_CLS is kept so panels use
# both tokens together (scale_colour_manual + scale_shape_manual) and a shape cue can be
# restored with one edit.
# Blue is close to the CON token, accepted.
PAL_ZHOU_CLS   <- c("Zhou-sig, concordant" = "#0072B2",
                    "Zhou-sig, discordant" = "#D55E00",
                    "n.s. in Zhou"         = "grey78")
SHAPE_ZHOU_CLS <- c("Zhou-sig, concordant" = 16, "Zhou-sig, discordant" = 16, "n.s. in Zhou" = 16)

# Regions — neuroanatomical mnemonic palette
# Frontal = dusty aubergine (prefrontal "higher-order" violet)
# Hippo   = moss green (limbic earth-tone)
# OCC dropped: this rebuild is Frontal + Hippocampus only.
# Chosen to be distinct from PAL_CELLTYPE (Astro teal, Neuron_Ex mustard, OPC
# salmon) and PAL_COND (CON blue, NHD red) — no cross-palette collisions.
PAL_REGION <- c(
  "Frontal" = "#6B4E7D",   # dusty aubergine
  "Hippo"   = "#4F7942"    # moss green
)

# Pale region tones — PAL_REGION blended 55% toward white. These are the same
# tones as the region facet-strip BANNERS (strip_region_x below) and, from
# Defined via the same blend used by strip_region_x (not hardcoded) so the two
# stay byte-identical: Frontal #BCAFC4, Hippo #AFC2A9.
# Saturated PAL_REGION is unchanged — it still feeds the banners (via this
# blend) and any other saturated use (heatmap column strips, etc.).
.region_lighten <- function(hex, f = 0.55) {
  m <- grDevices::col2rgb(hex)
  grDevices::rgb(t(m * (1 - f) + 255 * f), maxColorValue = 255)
}
PAL_REGION_PALE <- setNames(.region_lighten(PAL_REGION, 0.55), names(PAL_REGION))

# Vivid region tones for UMAP points only (vibrancy pass).
# History: points were once the 0.55 pale banner tone (invisible), then a 0.32
# white blend (mauve/sage — still washed). A white blend can only ever desaturate,
# so the UMAP tones are now set directly: the same violet/green mnemonic as
# PAL_REGION (frontal violet, hippocampal moss) at journal saturation.
# The facet-strip banners are unchanged — they still use PAL_REGION_PALE (0.55)
# with black text, so banner identity and legibility are untouched.
# Paled 10% (same pass as PAL_CELLTYPE).
PAL_REGION_UMAP <- setNames(
  .pale(c("Frontal" = "#8E5BB5",   # vivid violet  (banner #BCAFC4)
          "Hippo"   = "#5AA24A"),  # vivid moss    (banner #AFC2A9)
        0.10),
  c("Frontal", "Hippo"))

# Full region display names — use everywhere a region is
# LABELLED (strip banners, axes, legends). Data/factor levels stay the short
# codes; this maps code -> printed name. Pass as `labeller = REGION_LABELLER`
# to facet_*(), or `labels = REGION_FULL[levels]` to scale_*_discrete().
REGION_FULL <- c(Frontal = "Frontal", Hippo = "Hippocampus")
REGION_LABELLER <- ggplot2::as_labeller(REGION_FULL)

# Canonical region order for this rebuild (Frontal first, then Hippocampus).
REGION_ORDER <- c("Frontal","Hippo")

# Helper: coloured facet strip for region-faceted panels (ggh4x).
# Returns a `strip_themed` object you pass to facet_grid2 / facet_wrap2
# so the Region strip header is filled with PAL_REGION_PALE and text is black.
#   Usage:
#     library(ggh4x)
#     p + facet_grid2(sig ~ Region, strip = strip_region_x(c("Frontal","Hippo"), clip = "off"))
#   The single-NHD-Hippocampus lane caveat
#     now lives in the legend, never on the plot. This argument is kept only so
#     the many existing callers don't error; it is ignored — every region strip
#     gets the normal pale region tone with black text, including Hippocampus.
strip_region_x <- function(regions = c("Frontal","Hippo"),
                           fontsize = 7, text_colour = "black",   # black labels on the colour banners
                           low_conf = character(0),               # ignored — see note above
                           low_conf_fill = "grey80",              # IGNORED (kept for signature compat)
                           low_conf_text = "grey25",              # IGNORED (kept for signature compat)
                           # clip: "inherit" (default) keeps prior behaviour for every
                           # existing caller. Pass clip="off" for panels whose facets
                           # are so narrow the banner would truncate the (full) region
                           # name, letting the centred label overflow the pale banner
                           # instead of being cut (theme(strip.clip=) does not reach
                           # ggh4x strip_themed).
                           clip = "inherit") {
  if (!requireNamespace("ggh4x", quietly = TRUE))
    stop("strip_region_x() needs the ggh4x package.")
  # Pale region banners: PAL_REGION blended ~55% toward white so
  # black strip labels read clearly (the full-saturation banners washed them out).
  # Uses the shared PAL_REGION_PALE constant (defined above from the identical
  # blend) so banners stay byte-identical and match the region UMAP points.
  # low_conf is deliberately not applied: Hippocampus gets the normal
  # pale tone like every other region; no grey desaturation, no on-plot caveat.
  fills <- unname(PAL_REGION_PALE[regions])
  txt   <- rep(text_colour, length(regions))
  # Scalar fontsize styled only the first strip and every later one silently fell
  # back to theme_pub()'s strip.text size (8). In a 2-region panel that rendered
  # "Hippocampus" visibly larger than "Frontal". fontsize/colour are now recycled to length(regions), so the
  # requested size applies to every banner. Panels rendered before this date carry
  # the oversized second strip until they are re-rendered.
  fsz   <- rep_len(fontsize, length(regions))
  ggh4x::strip_themed(
    clip         = clip,
    background_x = ggh4x::elem_list_rect(fill = fills, colour = NA),
    text_x       = ggh4x::elem_list_text(colour = txt,
                                         size   = fsz,
                                         face   = "plain")
  )
}

# LOWCONF_MARK is now an empty string and both labellers below ignore their
# `low_conf` argument (kept only for caller-signature compatibility). Region
# strips therefore read a clean "Frontal" / "Hippocampus" with no suffix.
LOWCONF_MARK <- ""
region_lowconf_labeller <- function(low_conf = character(0),
                                    mark = LOWCONF_MARK) {
  # low_conf / mark ignored — no suffix is ever appended.
  ggplot2::as_labeller(function(x) x)
}
# full-name labeller (region facet strips use the full name
# Frontal/Hippocampus in all panels). Maps the short code to REGION_FULL. The
# `low_conf` argument is ignored — no "(low n)" suffix.
region_full_lowconf_labeller <- function(low_conf = character(0),
                                         mark = LOWCONF_MARK) {
  ggplot2::as_labeller(function(x)
    ifelse(x %in% names(REGION_FULL), unname(REGION_FULL[x]), x))
}
# Retired on-plot footnote (kept as an empty string so any caller referencing it
# prints nothing). The lane caveat now lives in the figure legend only.
LOWCONF_CAPTION <- ""

# --- Theme function -------------------------------------------------
# no gridlines, thin black axis lines, all text plain (house rule: no bold).
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

# --- sender marker for ligand-receptor panels -------------------------------
# "adding the cell type sender with an arrow really meets the
# purpose, so lets propagate this format to all the l-r panels in all the figures."
#
# An l-r dot plot names its RECEIVERS on the x axis but had no way to say who is
# sending, so the reader had to get it from the legend. This returns a patchwork
# header strip — "Microglia ->" — to sit above the plot.
#
# Why a header strip and not annotate(): a faceted plot has no single corner to
# anchor to, and an annotation placed outside the coord limits with clip = "off" is
# exactly what has been clipping labels elsewhere in this project. A header cannot
# clip and cannot overlap data.
# Why a drawn segment and not a glyph: base pdf() cannot encode U+2192 (it has
# already broken a render here). plotmath's %->% does work and self-spaces, but it
# renders a THIN arrow; the filled arrowhead below is the approved look.
#
# arrow_x defaults to ~0.102 data units per character, calibrated by eye on
# "Microglia" (9 chars -> 0.92) at the ~2.4 in panel width these panels share. Pass
# an explicit value for a panel whose plotting area is much wider or narrower, and
# always eye-check the rendered header rather than trusting the estimate.
sender_header <- function(label, size = 2.6, arrow_x = NULL, xmax = 3.2,
                          panel_width = 4.35) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("sender_header() needs patchwork.")
  # width-aware spacing. The header's x-scale is fixed at 0..xmax
  # while the label is sized in mm, so the label's width in data UNITS is inversely
  # proportional to the saved panel width: the same 0.102/char that is correct at
  # 4.35 in is far too small at 3.9 in (arrow lands on the word) and far too large at
  # 6.4 in (arrow drifts away from it). Scale the offset by the width ratio, with
  # 4.35 in as the calibration point where it was verified by eye.
  if (is.null(arrow_x))
    arrow_x <- 0.102 * nchar(label) * (4.35 / panel_width)
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = label,
                      hjust = 0, vjust = 0.5, size = size, colour = "black") +
    ggplot2::annotate("segment", x = arrow_x, xend = arrow_x + 0.28 * (4.35 / panel_width),
                      y = 0, yend = 0,
                      arrow = grid::arrow(length = grid::unit(0.045, "in"),
                                          type = "closed"),
                      linewidth = 0.4, colour = "black") +
    ggplot2::scale_x_continuous(limits = c(0, xmax)) +
    ggplot2::scale_y_continuous(limits = c(-0.5, 0.5)) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(-0.5, 0.5), clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 0, 3, 2))
}

# Stack a sender header above an l-r panel. Height fraction is small enough that the
# plot area is barely reduced; the caller should add ~0.15 in to the saved height so
# the panel body keeps its size.
# Always pass panel_width = the width this panel is saved at, or the arrow will be
# mis-spaced. Eye-check the rendered header; the estimate is a starting point.
with_sender <- function(p, label, size = 2.6, panel_height = NULL,
                        height = NULL, ...) {
  # Height now FOLLOWS the text size. The header row was a fixed
  # 0.055 of the panel height, which is ~0.15 in on a 2.85 in panel -- fine for the old
  # 2.6 mm label but shorter than the GLYPHS once the apparent-type parity pass raised
  # the sender labels to ~3.4 mm, so "Astrocytes" rendered with its top sliced off.
  # Give the row the label's own height (cap + descender + breathing room) whenever the
  # caller tells us how tall the panel is.
  if (is.null(height)) {
    hdr_in <- size / 25.4 * 1.90 + 6 / 72   # glyph box + the 3 pt top/bottom margins
    height <- if (is.null(panel_height)) 0.075 else hdr_in / (panel_height - hdr_in)
  }
  patchwork::wrap_plots(sender_header(label, size = size, ...), p, ncol = 1,
                        heights = c(height, 1))
}

# =============================================================================
# axis_arrows() — the house UMAP axis key: an L of two closed arrows in the lower
# left with "UMAP 1" / "UMAP 2" beside them, in place of axes on an embedding.
# Lifted verbatim from 64_fig4_neuron_umap_FH.R (Figure 5) so the supplementary
# UMAPs cannot drift from the main figures; scripts that already define their own
# local copy keep it (a local definition masks this one).
#
# House RULE, from the same reference: one axis key per composite, at the reading
# entry point -- the leftmost panel of a side-by-side composite, the lower panel
# of a stacked one. Never repeat it on every sub-panel.
#
# Xr / yr are the view limits (what coord_cartesian shows), not the raw data range,
# so a cropped embedding still puts the key inside the frame. head_in and linewidth
# default to the Figure-5 values at base_size 8; scale them with the panel's
# declared type (head_in = 0.055 * base/8) so the key keeps its apparent weight.
# The panel needs clip = "off" -- the labels sit just outside the data range.
# =============================================================================
axis_arrows <- function(xr, yr, frac = 0.16, lab = c("UMAP 1", "UMAP 2"),
                        lab_size = 2.5, head_in = 0.055, linewidth = 0.4) {
  rng <- max(diff(xr), diff(yr))
  x0  <- xr[1] - 0.01 * diff(xr); y0 <- yr[1] - 0.01 * diff(yr)
  a   <- grid::arrow(length = grid::unit(head_in, "in"), type = "closed")
  list(
    ggplot2::annotate("segment", x = x0, xend = x0 + frac * rng, y = y0, yend = y0,
                      arrow = a, linewidth = linewidth),
    ggplot2::annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + frac * rng,
                      arrow = a, linewidth = linewidth),
    ggplot2::annotate("text", x = x0, y = y0 - 0.03 * rng, hjust = 0, vjust = 1,
                      label = lab[1], size = lab_size, colour = "black"),
    ggplot2::annotate("text", x = x0 - 0.03 * rng, y = y0, angle = 90, hjust = 0,
                      vjust = 0, label = lab[2], size = lab_size, colour = "black"))
}
