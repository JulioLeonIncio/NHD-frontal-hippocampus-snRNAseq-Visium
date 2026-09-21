#!/usr/bin/env Rscript
# =============================================================================
# 4E_micro_to_ExInh_chord_FH.R — F5e1_micro_to_ExInh_chord_FH.R — Figure 4 (Part B): CON vs NHD chord plots of the microglia -> Ex/Inh neuron signalling programme, frontal + hippocampus.
# Faithful port of manuscript_7fig/scripts/10e_v3_5E_micro_to_ExInh_chord.R to
# the FH rebuild (3 regions -> 2; OCC dropped).
# -----------------------------------------------------------------------------
# 2x2 layout (Region rows x Condition cols) + right pathway legend column:
#        Frontal (CON)      |  Frontal (NHD)
#        Hippocampus (CON)  |  Hippocampus (NHD)
# Perimeter sectors: Micro-PVM, Ex, Inh (3 sectors).  Ribbons coloured by the
# micro->neuron ligand axis (the significant FH set only).
#
# KEY FH difference vs the reference: the FH withNeurons_v2 CellChat objects were
# fit at major-class resolution (targets are already "Neuron_Ex"/"Neuron_Inh"),
# so no subtype roll-up is needed — the reference's Azimuth/Tippani collapse maps
# are replaced by a direct Neuron_Ex->Ex / Neuron_Inh->Inh rename.
#
# Critical gotcha (project): CellChat pathway_name != ligand (GAS6's pathway is
# "GAS", ENTPD1's is "CD39", SEMA4D's is "SEMA4", GRN's is "GRN"). We select on
# ligand and label ribbons at the ligand-receptor level so "GAS6" is not misread
# as an abbreviation.
#
# FH LR REALITY (the 5 axes):
#   SPP1  -> integrins (ITGAV/ITGA4/8/9+ITGB1)   present both regions.
#   GAS6  -> TYRO3      present; up in NHD both regions (trophic ligand).
#   GRN   -> SORT1      Hippocampus-NHD-restricted (tiny).
#   SEMA4D-> PLXNB2     CON-skewed (lost in NHD).
#   ENTPD1-> ADORA1     CON-skewed (purinergic brake, lost in NHD).
# Valence SPLIT (house rule): SPP1 (phagocytic/synaptotoxic) and GAS6 (trophic,
# but uncoupled — its neuronal TAM receptor TYRO3 falls) are shown as separate
# ligand ribbons, never lumped. No causal verbs anywhere.
#
# Output: figures/Figure_5/panels/F5e1_micro_to_ExInh_chord.{pdf,png} + stats CSV.
# -----------------------------------------------------------------------------
# apparent-type parity PASS  ("adjust the font/panel size ratio
#   so everything looks more or less uniform, favouring the readability of the figure")
# Measured from the assembled Fig5.pdf: each panel is placed as a linked image with a
# known scale (parsed from the page's placement matrices). On-page type was 3.2-3.9 pt
# for most panels but 6.4-6.9 pt for the chord -- a 2x spread, and most of it below the
# print-legibility floor. Declared sizes here are therefore set to
#     declared = apparent_target / placement_scale
# with one set of apparent targets shared by every Figure-5 panel:
#     tick / legend text 5.2 pt | legend title 5.4 | axis title & facet strip 6.0
#     panel title 6.4 | in-panel gene labels & annotations 5.4
# canvas size and PIXEL DIMENSIONS are unchanged, so re-linking in Illustrator is a
# no-op: the frame keeps its box and only the type grows.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(CellChat); library(dplyr); library(tidyr)
  library(circlize); library(grid)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # REGION_FULL, PAL_COND

CC    <- file.path(PROJ, "data", "cellchat_per_region")
PANEL <- file.path(PROJ, "figures", "Figure_5", "panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

REGIONS <- c("Frontal","Hippo")   # canonical FH order (OCC dropped)
CONDS   <- c("CON","NHD")
SENDER  <- "Micro-PVM"
NEURON  <- c("Neuron_Ex","Neuron_Inh")

SECTORS <- c("Micro-PVM","Ex","Inh")
DISP    <- c("Micro-PVM" = "Micro", "Ex" = "Ex", "Inh" = "Inh")

# Curated micro->neuron ligand axes (data keys = CellChatDB pathway_name; the
# five FH-significant axes). Existence-guarded below.
TARGET_PATHS <- c("SPP1","GAS","CD39","GRN","SEMA4")
PATH_PAL <- c(
  SPP1  = "#E8A33D",   # gold            — phagocytic/synaptotoxic axis
  GAS   = "#1F6FB2",   # saturated blue  — trophic-but-uncoupled (GAS6->TYRO3)
  CD39  = "#2CA089",   # teal            — purinergic brake (ENTPD1->ADORA1)
  GRN   = "#B49AC6",   # muted lilac     — lysosomal (GRN->SORT1; Hippo-NHD only)
  SEMA4 = "#D55E00")   # orange          — SEMA4D->PLXNB2 (CON-skewed)
# Display names at the ligand-receptor level (house rule: expand GAS6 not GAS).
PATH_DISP <- c(SPP1 = "SPP1-integrin",
               GAS  = "GAS6-TYRO3",
               CD39 = "ENTPD1-ADORA1",
               GRN  = "GRN-SORT1",
               SEMA4= "SEMA4D-PLXNB2")
# Sector palette: in this micro->neuron chord the interneuron receiver (Inh) is the
# house neuron magenta (#C42E7B), which is nearly identical to the microglia
# PAL_CELLTYPE magenta (#B83C7E) — sender and one receiver would collide. So the
# Micro SENDER gets a distinct deep-violet accent here so
# Micro / Ex / Inh are three unambiguous sector colours. Ex teal + Inh magenta are
# the locked house neuron-class hues and are unchanged.
SECTOR_PAL <- c(
  "Micro-PVM" = "#5E3C99",   # deep violet — microglia sender (distinct from Inh magenta)
  "Ex"        = "#006D77",   # neuron teal (house)
  "Inh"       = "#C42E7B")   # neuron magenta (house)
FLOOR_PROB <- 1e-4

# ---- Load per-region withNeurons CellChat merged objects -------------------
load_region <- function(region) {
  path <- file.path(CC, sprintf("merged_NHDvsCON_%s_withNeurons_v2.rds", region))
  if (!file.exists(path))
    stop(sprintf("MISSING %s — CellChat withNeurons prep for %s failed", path, region))
  m  <- readRDS(path)
  df <- subsetCommunication(m, slot.name = "net")
  bind_rows(df$CON %>% mutate(Condition = "CON"),
            df$NHD %>% mutate(Condition = "NHD")) %>% mutate(Region = region)
}

raw_all <- bind_rows(lapply(REGIONS, load_region)) %>% as_tibble() %>%
  filter(source == SENDER, target %in% NEURON) %>%
  # FH major-class objects: rename Neuron_Ex->Ex, Neuron_Inh->Inh (no roll-up).
  mutate(src = "Micro-PVM",
         tgt = dplyr::recode(as.character(target),
                             Neuron_Ex = "Ex", Neuron_Inh = "Inh"))
stopifnot(nrow(raw_all) > 0)

# --- Existence guard: drop curated pathways absent from the FH objects -------
present_pw <- sort(unique(raw_all$pathway_name))
absent     <- setdiff(TARGET_PATHS, present_pw)
if (length(absent) > 0) {
  cat("WARNING dropping absent pathway(s) [not micro->neuron in FH objects]:",
      paste(absent, collapse = ", "), "\n")
  TARGET_PATHS <- intersect(TARGET_PATHS, present_pw)
  PATH_PAL     <- PATH_PAL[TARGET_PATHS]
}
stopifnot(length(TARGET_PATHS) >= 1)
cat("Chord ligand axes used:", paste(TARGET_PATHS, collapse = ", "), "\n")

raw <- raw_all %>% filter(pathway_name %in% TARGET_PATHS)
cat(sprintf("Micro-PVM -> {Ex,Inh} rows (curated axes): %d\n", nrow(raw)))

# Aggregate to (Region, Condition, src, tgt, pathway) sum-prob, log1p-compress.
agg <- raw %>%
  group_by(Region, Condition, src, tgt, pathway_name) %>%
  summarise(prob_sum = sum(prob, na.rm = TRUE), .groups = "drop") %>%
  filter(prob_sum >= FLOOR_PROB) %>%
  mutate(weight = log1p(prob_sum * 20))
write.csv(agg, file.path(TDIR, "fig4E_micro_neuron_chord_stats.csv"), row.names = FALSE)
cat(sprintf("Wrote stats: %d (region x cond x tgt x pathway) entries\n", nrow(agg)))

# Ported from 3A_micro_glia_chord_FH.R ("labels are hard to read,
# why don't you copy the style of the astrocyte chord in terms of font size").
# This panel never received it.
# circos/base graphics scale text as pointsize * cex, and both pdf() and png() default to
# pointsize 12, so these are absolute points: sector ~11.4 pt, key ~10.6 pt, title ~11 pt.
# Titles are black, not tinted by condition — the condition is spelled out in the title
# text, so the colour carried no information and was low-contrast on white.
CEX_TITLE <- 0.88; CEX_SECTOR <- 0.78; CEX_KEY <- 0.72; CEX_HDR <- 0.75

draw_chord <- function(df_rc, panel_title, title_col = "black") {
  if (nrow(df_rc) == 0) {
    plot.new(); title(panel_title, col.main = title_col, cex.main = CEX_TITLE, col.main = "black",
                      line = -0.7, font.main = 1)
    text(0.5, 0.5, "no links", col = "black", cex = 0.85); return(invisible(NULL))
  }
  order_sectors <- intersect(SECTORS, unique(c(df_rc$src, df_rc$tgt)))
  grid_col <- SECTOR_PAL[order_sectors]
  link_col <- PATH_PAL[df_rc$pathway_name]
  circos.clear()
  circos.par(start.degree = 90,
             gap.after   = setNames(rep(8, length(order_sectors)), order_sectors),
             track.margin = c(0.005, 0.005),
             # "the empty blank space in the chord panels, right in
             # the middle in between the chords"). canvas.xlim/ylim is the drawing box the circle
             # is fitted into, so a larger limit shrinks the chord inside its cell and the padding
             # shows up as a gap between the four panels. Tightening it draws each chord larger in
             # the same cell and closes the gap without changing any text size.
             canvas.xlim  = c(-1.06, 1.06), canvas.ylim = c(-1.06, 1.06),
             points.overflow.warning = FALSE)
  chordDiagram(
    x = df_rc[, c("src","tgt","weight")], order = order_sectors,
    grid.col = grid_col, col = link_col, transparency = 0.25,
    directional = 1, direction.type = "arrows", link.arr.type = "big.arrow",
    link.arr.length = 0.06, link.lwd = 0.15, link.border = NA,
    annotationTrack = "grid", annotationTrackHeight = mm_h(2.6),
    preAllocateTracks = list(track.height = 0.025),
    link.sort = TRUE, link.decreasing = TRUE)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    nm   <- get.cell.meta.data("sector.index")
    disp <- DISP[nm]; if (is.na(disp)) disp <- nm
    # facing="downward": all sector labels drawn perfectly horizontal, no per-letter
    # bending along the arc. The bending faces ("outside"/"bending.outside") re-arced
    # the glyphs so a wide "Ex" sector splayed the two letters and read as "ex" while
    # a narrow sector read "Ex" — a casing mismatch across sibling sub-panels.
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + 1.9, disp,
                # Curved along the arc, as the astro chord does ("in the
                # astrocyte chord the cell type labels follow the chord round shape").
                # "downward" renders them flat and they read as pinned on rather than part
                # of the ring; niceFacing keeps them upright where the arc would invert them.
                facing = "bending.outside", niceFacing = TRUE, adj = c(0.5, 0.5),
                cex = CEX_SECTOR, col = "black", font = 1)
  }, bg.border = NA)
  title(panel_title, col.main = "black", cex.main = CEX_TITLE, line = -0.6, font.main = 1)
}
slice_rc <- function(region, condition)
  agg %>% filter(Region == region, Condition == condition) %>%
    select(src, tgt, weight, pathway_name)

# 1x4 row layout at the astro chord's geometry (was 2x2 at 3.24 x 2.90, which left the
# labels large relative to a small panel and stretched the legend down the full height).
W <- 3.56; H <- 2.30   # 2x2 + legend; cells stay square (1.148 in); 10% smaller
out_pdf <- file.path(PANEL, "F5e1_micro_to_ExInh_chord.pdf")
out_png <- file.path(PANEL, "F5e1_micro_to_ExInh_chord.png")

render <- function(fn_open, fn_args, fn_close) {
  do.call(fn_open, fn_args)
  # 2x2 grid + full-height legend column (after trying both 1x4 and a
  # single vertical column). Region rows x condition columns puts the CON/NHD comparison
  # side by side within a region — the contrast the panel exists to show — while the
  # legend keeps a tall column of its own, which is what the long ligand-RECEPTOR names
  # need (the astro chord's short pathway names fit a 1x4 row; "SEMA4D-PLXNB2" does not).
  # Region 5 is repeated down the third column so the legend spans both rows.
  layout(matrix(c(1, 2, 5,
                  3, 4, 5), ncol = 3, byrow = TRUE),
         # Widths chosen so each chord cell is square. A chord is a circle: it fits the
         # smaller of the cell's two dimensions, so a cell wider than it is tall leaves
         # leftover width as a vertical white band between the columns
         # spotted. Row height here is H/2 = 1.275 in, so each chord column must be 1.275 in
         # too: 2 x 1.275 + a 1.40 in legend column = 3.95 in total.
         widths = c(1, 1, 1.10))
  par(mar = c(0.05, 0.05, 0.35, 0.05), xpd = NA)
  for (r in REGIONS) for (cd in CONDS) {
    ttl_col <- if (cd == "NHD") unname(PAL_COND["NHD"]) else unname(PAL_COND["CON"])
    ttl <- sprintf("%s  |  %s", REGION_FULL[r], cd)
    draw_chord(slice_rc(r, cd), ttl, title_col = ttl_col)
  }
  par(mar = c(0.2, 0.05, 0.2, 0.1)); plot.new()
  leg_y <- 0.80 - (seq_along(PATH_PAL[TARGET_PATHS]) - 1) * 0.085
  for (i in seq_along(TARGET_PATHS)) {
    p <- TARGET_PATHS[i]
    points(0.08, leg_y[i], pch = 22, cex = 1.9, bg = PATH_PAL[p], col = "black")
    disp_p <- if (!is.na(PATH_DISP[p])) PATH_DISP[[p]] else p
    text(0.24, leg_y[i], disp_p, cex = CEX_KEY, col = "black", adj = c(0, 0.5))
  }
  # "Ligand axis (ribbon colour)" is caption text and clips a narrow legend column;
  # shortened to one word, as the astro chord did.
  text(0.02, 0.95, "Ligand", cex = CEX_HDR, font = 1, col = "black", adj = c(0, 1))
  do.call(fn_close, list())
}
render(pdf, list(file = out_pdf, width = W, height = H, useDingbats = FALSE), dev.off)
render(png, list(filename = out_png, width = W, height = H, units = "in", res = 600), dev.off)

cat("\nWrote:\n  ", out_pdf, "\n  ", out_png, "\n  ",
    file.path(TDIR, "fig4E_micro_neuron_chord_stats.csv"), "\n", sep = "")
cat("=== DONE ===\n", file = stderr())
