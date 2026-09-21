#!/usr/bin/env Rscript
# =============================================================================
# 3A_micro_glia_chord_FH.R — Figure 3 panel a: CON vs NHD chord plots of the curated micro->glia signalling programme, frontal + HIPPOCAMPUS (2 regions).
# Faithful port of 09a2_v2_micro_glia_chord_by_region.R to the FH rebuild.
# -----------------------------------------------------------------------------
# 2x2 layout:
#        Frontal (CON)      |  Frontal (NHD)
#        Hippocampus (CON)  |  Hippocampus (NHD)
# Perimeter sectors: Astro, OPC, Oligo, Micro-PVM (PAL_CELLTYPE hues).
# Ribbons coloured by pathway, restricted to the curated micro->glia programme.
#
# CRITICAL (project gotcha): CellChat pathway_name != ligand.  GAS6's pathway is
# "GAS", PROS1's is "PROS", PSAP's is "PSAP", SPP1's is "SPP1".  We select on the
# curated pathway set but the ligand-level identity is what carries the biology.
#
# FH-RE-VALIDATION of the reference 4-pathway set
# {SPP1, GAS, PSAP, SEMA4} against tables/cellchat_all_LR_pairs.csv:
#   * SPP1  micro->{Astro(CD44),OPC/Oligo(integrins)}  present CON+NHD, both regions.
#   * GAS   (GAS6->MERTK[astro]/TYRO3[OPC,oligo]/AXL)   present CON+NHD, both regions;
#           GAS6->astro-MERTK up in NHD (F 0.003->0.029, H 0.008->0.013).
#   * PSAP  (PSAP->GPR37[oligo]/GPR37L1[astro,OPC])     NHD-ONLY (absent in CON).
#   * PROS  (PROS1->MERTK/TYRO3/AXL) — second TAM ligand, nhd-only.  added to the
#           curated set (it co-carries the GAS6 TAM biology and only shows in NHD,
#           so it belongs in the micro->glia clearance programme).
#   * SEMA4 (SEMA4D->PLXNB1) present; kept.
#   Held out as ambient/off-programme (documented): Glutamate (SLC1A3/GLS artefact),
#   ADGRL (TENM4 teneurin adhesion), PDGF (endothelial PDGFB/C), vista, CD39, WNT.
#
# Ribbon width = summed CellChat `prob` for all LR entries of that
# (source -> target) pair within that pathway.  Reads the aggregated
# subsetCommunication output (tables/cellchat_all_LR_pairs.csv) directly — no
# heavy merged-CellChat object load needed (atlas-light).
#
# Output: figures/Figure_3/panels/3a_micro_glia_chord.{pdf,png} + stats CSV.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(circlize); library(grid)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # REGION_FULL, PAL_CELLTYPE, PAL_COND

LR    <- file.path(PROJ, "tables", "cellchat_all_LR_pairs.csv")
PANEL <- file.path(PROJ, "figures", "Figure_3", "panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(LR)) stop(sprintf("MISSING %s — CellChat prep step failed", LR))

# ---- Config ----------------------------------------------------------------
# ---------------------------------------------------------------------------
# CHORD_SET. One parameterised script, not three forks:
#   micro_astro  (default) -> Micro-PVM <-> Astro           -> Figure_3/F3g1_micro_astro_chord
#   astro_neuron           -> Astro <-> Neuron_Ex/Inh       -> Figure_3/F3h1_astro_neuron_chord
#   micro_glia             -> the original 4-way glia chord -> Figure_3/3a_micro_glia_chord
# The two sets read different sources: the aggregated glia LR table contains only
# Astro/Micro-PVM/Oligo/OPC (no neurons), so astro_neuron must come from the
# *_withNeurons_v2 CellChat objects. Provenance of those was verified before use: built
# after the FH atlas (14:08), inside this rebuild, Frontal+
# Hippo only -- not the old 3-region manuscript_7fig analysis.
# CHORD_MODE: "state" keeps the CON | NHD panels; "delta" collapses each region to
# one chord of the NHD-CON change. Needed because a state chord is dominated by whatever pathway
# is largest, not by whatever moved: in astro->neuron, NRXN is ~70% of the probability mass in
# all four panels, so they looked near-identical while Frontal NCAM was rising +604%
# (0.148 -> 1.040) and Hippocampus CADM was falling -65%. Plotting the change puts the movers in
# front. micro_astro keeps state (GAS/PROS genuinely expand there, so CON vs NHD reads).
# Default mode is per chord SET. astro_neuron defaults to "share" because one
# pathway (NRXN, ~70% of the mass) otherwise fixes the shape of all four panels; micro_astro
# stays "state" because GAS/PROS genuinely expand there against a stable SPP1 background.
# The default mode writes the canonical filename -- no suffix -- so the improved panel is
# F3h1_astro_neuron_chord and no stale variant can be picked up at assembly.
.default_mode <- if (tolower(Sys.getenv("CHORD_SET", "micro_astro")) == "astro_neuron")
                   "share" else "state"
CHORD_MODE <- tolower(Sys.getenv("CHORD_MODE", .default_mode))
stopifnot("CHORD_MODE must be state | share | delta" =
            CHORD_MODE %in% c("state","share","delta"))
CHORD_SET <- tolower(Sys.getenv("CHORD_SET", "micro_astro"))
stopifnot("CHORD_SET must be micro_astro | astro_neuron | micro_glia" =
            CHORD_SET %in% c("micro_astro","astro_neuron","micro_glia"))
CELLS   <- switch(CHORD_SET,
                  micro_astro  = c("Micro-PVM","Astro"),
                  astro_neuron = c("Astro","Neuron_Ex","Neuron_Inh"),
                  micro_glia   = c("Astro","OPC","Oligo","Micro-PVM"))
BN_SET  <- switch(CHORD_SET,
                  micro_astro  = "F3g1_micro_astro_chord",
                  astro_neuron = "F3h1_astro_neuron_chord",
                  micro_glia   = "3a_micro_glia_chord")
cat(sprintf("CHORD_SET=%s | sectors: %s | out: %s\n",
            CHORD_SET, paste(CELLS, collapse=", "), BN_SET))
REGIONS <- c("Frontal","Hippo")                    # canonical FH order
CONDS   <- c("CON","NHD")
# Curated micro->glia programme (FH-revalidated above). pros added vs reference.
# astro_neuron uses exactly the families already vetted in the 4G astro->neuron bubble panel
# (NRXN / NCAM / CADM / NRG), so the chord and that panel show the same programme two ways and
# there is no second curation to defend. Held out and documented: LAMININ + COLLAGEN (ECM bulk
# -- they dominate on ribbon width without being a signalling programme), ADGRL (teneurin
# adhesion, already held out of the micro->glia chord as off-programme), and BMP/EPHA/PTN/
# UNC5/JAM/CDH/SLITRK/FGF/PTPR/ANGPTL.
TARGET_PATHS <- switch(CHORD_SET,
                       astro_neuron = c("NRXN","NCAM","CADM","NRG"),
                       c("SPP1","GAS","PROS","PSAP","SEMA4"))
FLOOR_PROB   <- 1e-4                               # drop links below this (readability)

# pathway ribbon palette — distinct, high-contrast; SPP1 gold + GAS emerald lead.
# Ribbon palette. Now that the sectors carry the saturated cell-type hues, the ribbons must
# avoid those hues or the two encodings collide (a teal ribbon on a teal Astro sector is
# unreadable). Palettes are therefore chosen per chord against that chord's own sectors:
#   micro_astro  sectors = magenta (Micro) + teal (Astro)   -> ribbons avoid magenta/teal
#   astro_neuron sectors = teal (Astro) + gold (Ex) + purple (Inh) -> ribbons avoid those three
# All values are colour-blind-safe (Okabe-Ito, plus one deep purple / deep red where a 5th and
# 4th distinct hue is needed).
PATH_PAL_AN <- c(NRXN = "#0072B2",   # blue
                 # NCAM gets a bright orange. It is the one pathway that
                 # actually changes in Frontal (1.3% -> 9.3% of signalling) but at Okabe-Ito
                 # vermillion #D55E00 it read as a muted rust that sat too close to NRG's deep
                 # red, so the ribbons that carry the finding did not stand out. #FF6D00 is
                 # clearly orange, distinct from NRG red and from the gold Ex SECTOR (#F3B734,
                 # which is yellower), so it cannot be confused with either.
                 NCAM = "#FF6D00",   # bright orange -- the mover, deliberately loudest
                 CADM = "#56B4E9",   # sky blue
                 NRG  = "#B2182B")   # deep red
PATH_PAL <- c(
  SPP1  = "#E69F00",   # orange       (lead)
  GAS   = "#0072B2",   # blue         (TAM lead)
  PROS  = "#56B4E9",   # sky blue     (TAM co-ligand, same family as GAS)
  PSAP  = "#D55E00",   # vermillion
  SEMA4 = "#7D3C98")   # purple
if (CHORD_SET == "astro_neuron") PATH_PAL <- PATH_PAL_AN

# perimeter cell-node palette = PAL_CELLTYPE (matches the rest of the figures).
# Sector palette. Sectors are now PAL_CELLTYPE at full strength -- byte-identical to the Fig-1
# compact UMAP (Micro-PVM #DA3F92, Astro #19B99E, Neuron_Ex #F3B734, Neuron_Inh #A56FD4), so a
# cell type is the same colour everywhere in the paper. Earlier passes used a 0.55 white blend
# (washed-out pink) and then neutral greys; this is the cross-figure-consistent choice.
CELL_PAL <- c(
  Astro       = unname(PAL_CELLTYPE["Astro"]),
  `Micro-PVM` = unname(PAL_CELLTYPE["Micro-PVM"]),
  OPC         = unname(PAL_CELLTYPE["OPC"]),
  Oligo       = unname(PAL_CELLTYPE["Oligo"]),
  Neuron_Ex   = unname(PAL_CELLTYPE["Neuron_Ex"]),
  Neuron_Inh  = unname(PAL_CELLTYPE["Neuron_Inh"]))

# ---- Load LR (aggregated subsetCommunication output) -----------------------
load_withneurons <- function() {
  suppressPackageStartupMessages(library(CellChat))
  CCDIR <- file.path(PROJ, "data", "cellchat_per_region")
  bind_rows(lapply(REGIONS, function(rg) {
    f <- file.path(CCDIR, sprintf("merged_NHDvsCON_%s_withNeurons_v2.rds", rg))
    stopifnot("MISSING withNeurons CellChat object" = file.exists(f))
    d <- subsetCommunication(readRDS(f), slot.name = "net")
    bind_rows(d$CON %>% mutate(Condition = "CON"),
              d$NHD %>% mutate(Condition = "NHD")) %>% mutate(Region = rg)
  })) %>% as_tibble()
}

raw <- if (CHORD_SET == "astro_neuron") {
  load_withneurons() %>%
    filter(source %in% CELLS, target %in% CELLS, source != target,
           pathway_name %in% TARGET_PATHS,
           Region %in% REGIONS, Condition %in% CONDS)
} else read.csv(LR, stringsAsFactors = FALSE) %>%
  # region_cond is e.g. "CON_Frontal" -> split into Condition + Region
  mutate(Condition = sub("_.*$", "", region_cond),
         Region    = sub("^[^_]*_", "", region_cond)) %>%
  filter(source %in% CELLS, target %in% CELLS,
         pathway_name %in% TARGET_PATHS,
         Region %in% REGIONS, Condition %in% CONDS) %>%
  as_tibble()
stopifnot("no rows for this CHORD_SET" = nrow(raw) > 0)
cat(sprintf("Loaded %d curated %s programme rows\n", nrow(raw), CHORD_SET))
cat("Pathways retained:", paste(sort(unique(raw$pathway_name)), collapse = ", "), "\n")

# aggregate to (Region, Condition, source, target, pathway) sum of prob; ribbon
# weight = log1p-scaled so small but real programme links stay visible next to
# the dominant SPP1 links (visualisation choice; raw prob kept in the stats CSV).
agg <- raw %>%
  group_by(Region, Condition, source, target, pathway_name) %>%
  summarise(prob_sum = sum(prob, na.rm = TRUE), .groups = "drop") %>%
  filter(prob_sum >= FLOOR_PROB) %>%
  mutate(weight = log1p(prob_sum * 20))

if (CHORD_MODE == "share") {
  # Share mode.
  # Comparable. In state mode the ribbon mass is absolute probability, so whichever pathway is
  # biggest fixes the shape of every panel -- for astro->neuron NRXN is ~70% of the mass in all
  # four, which is why they looked identical. Here each panel is normalised to its own total, so
  # the chord shows composition and a pathway that grows from 1.3% to 9.3% of the signalling
  # (Frontal NCAM) is immediately visible. The absolute total is not lost -- it is printed in
  # each panel title, so magnitude is disclosed rather than implied by ribbon size.
  .tot <- agg %>% group_by(Region, Condition) %>%
    summarise(total = sum(prob_sum), .groups = "drop")
  agg <- agg %>% left_join(.tot, by = c("Region", "Condition")) %>%
    mutate(share = prob_sum / total, weight = share)
  PANEL_TOTAL <<- setNames(.tot$total, paste(.tot$Region, .tot$Condition))
  cat("SHARE mode: pathway share per panel (%)\n")
  print(as.data.frame(agg %>% group_by(Region, Condition, pathway_name) %>%
          summarise(pct = round(100 * sum(share), 1), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = Condition, values_from = pct)), row.names = FALSE)
}
if (CHORD_MODE == "delta") {
  # signed change per (Region, source, target, pathway); ribbon WIDTH = |delta|, and the
  # direction is carried by colour (see PATH_PAL override below), because a chord ribbon has no
  # natural sign of its own.
  agg <- agg %>%
    select(Region, Condition, source, target, pathway_name, prob_sum) %>%
    tidyr::pivot_wider(names_from = Condition, values_from = prob_sum,
                       values_fill = list(prob_sum = 0)) %>%
    mutate(delta = NHD - CON,
           dir   = ifelse(delta >= 0, "up", "down"),
           weight = log1p(abs(delta) * 20)) %>%
    filter(abs(delta) >= FLOOR_PROB)
  cat(sprintf("DELTA mode: %d link(s) with |NHD-CON| >= %.0e | up %d / down %d\n",
              nrow(agg), FLOOR_PROB, sum(agg$dir == "up"), sum(agg$dir == "down")))
  print(as.data.frame(agg %>% group_by(Region, pathway_name) %>%
          summarise(delta = round(sum(delta), 3), .groups = "drop") %>%
          arrange(Region, delta)), row.names = FALSE)
}
write.csv(agg, file.path(TDIR, paste0(BN_SET, "_stats.csv")), row.names = FALSE)
cat(sprintf("Wrote stats: %d (region x cond x src x tgt x pathway) entries\n", nrow(agg)))

# ---- Chord draw ------------------------------------------------------------
# text sizes. Base graphics
# scales text as pointsize * cex; both pdf() and png() default to pointsize 12, so the values
# below are absolute points. Set explicitly rather than left implicit, and raised for the two
# things a reader actually needs to decode -- the sector labels (which cell is which) and the
# pathway KEY -- while the panel titles come down slightly so they stop dominating.
#   panel title 0.82 -> 0.72  (~8.6 pt)
#   sector      0.84 -> 0.95  (~11.4 pt)
#   legend key  0.78 -> 0.88  (~10.6 pt)
#   legend hdr  0.68 -> 0.78  (~9.4 pt)
# titles were
# 0.72 and tinted by condition (blue CON / red NHD), which is low-contrast on white and reads
# poorly at final size. The condition is already spelled out in the title text itself, so the
# colour was carrying no unique information -- now black and larger.
CEX_TITLE <- 0.92; CEX_SECTOR <- 0.95; CEX_KEY <- 0.88; CEX_HDR <- 0.78

draw_chord <- function(df_rc, panel_title, title_col = "black") {
  if (nrow(df_rc) == 0) {
    plot.new(); title(panel_title, col.main = title_col, cex.main = CEX_TITLE, line = -0.6, font.main = 1)
    text(0.5, 0.5, "no links", col = "grey50", cex = 0.9); return(invisible(NULL))
  }
  order_sectors <- intersect(CELLS, unique(c(as.character(df_rc$source),
                                             as.character(df_rc$target))))
  grid_col <- CELL_PAL[order_sectors]
  link_col <- if (CHORD_MODE == "delta")
    c(up = "#C0392B", down = "#2C7FB8")[as.character(df_rc$dir)]
  else PATH_PAL[as.character(df_rc$pathway_name)]
  circos.clear()
  circos.par(start.degree = 90,
             # 18 deg between sectors was sized for the 4-way glia
             # chord; with 2-3 sectors it opened large wedges of empty canvas. 10 deg keeps
             # the sectors visually separate without the waste.
             gap.after    = setNames(rep(10, length(order_sectors)), order_sectors),
             track.margin = c(0.005, 0.005),
             # canvas was +-1.28 around a radius-1 circle -- a 0.28 empty ring on every
             # side, ~22% of each cell. +-1.14 still clears the bending sector labels.
             canvas.xlim  = c(-1.14, 1.14), canvas.ylim = c(-1.14, 1.14),
             points.overflow.warning = FALSE)
  chordDiagram(
    x = df_rc[, c("source","target","weight")],
    order = order_sectors, grid.col = grid_col, col = link_col,
    transparency = 0.25, directional = 0, link.lwd = 0.15, link.border = NA,
    annotationTrack = "grid", annotationTrackHeight = mm_h(3),
    preAllocateTracks = list(track.height = 0.02),
    link.sort = TRUE, link.decreasing = TRUE)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    nm   <- get.cell.meta.data("sector.index")
    disp <- c("Astro"="Astro","OPC"="OPC","Oligo"="Oligo","Micro-PVM"="Micro",
              "Neuron_Ex"="Ex","Neuron_Inh"="Inh")[nm]
    if (is.na(disp)) disp <- nm
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + 1.9, disp,
                facing = "bending.outside", niceFacing = TRUE, adj = c(0.5, 0.5),
                cex = CEX_SECTOR, col = "black", font = 1)
  }, bg.border = NA)
  title(panel_title, col.main = title_col, cex.main = CEX_TITLE, line = -0.6, font.main = 1)
}
slice_rc <- function(region, condition) {
  # braces required: a bare `else` on its own line at the top level of a function body is a
  # parse error in R.
  if (CHORD_MODE == "delta") {
    agg %>% filter(Region == region) %>% select(source, target, weight, dir)
  } else {
    agg %>% filter(Region == region, Condition == condition) %>%
      select(source, target, weight, pathway_name)
  }
}

# ---- Render (2 region rows x 2 condition cols + right legend col) ----------
# Tightened by removing empty
# canvas rather than by shrinking content: sector gaps 18->10 deg, chord canvas +-1.28 ->
# +-1.14, sector labels pulled in 2.5 -> 1.9, legend keys on a fixed pitch. 3.24 x 2.90 ->
# 2.70 x 2.35 (area 9.4 -> 6.3 in2, a 33% saving) with the same ribbons and type sizes.
# 1x4 row layout. Four chord cells side by side plus the legend column; height
# is now a single chord row rather than two.
# delta legend keys ("down in NHD") are longer than the pathway names, so the 2-chord
# layout still needs ~3.7 in for the key column not to clip.
W <- if (CHORD_MODE == "delta") 3.70 else 5.40; H <- 1.55
out_pdf <- file.path(PANEL, paste0(BN_SET, if (CHORD_MODE == .default_mode) "" else paste0("_", CHORD_MODE), ".pdf"))
out_png <- file.path(PANEL, paste0(BN_SET, if (CHORD_MODE == .default_mode) "" else paste0("_", CHORD_MODE), ".png"))

render <- function(fn_open, fn_args, fn_close) {
  do.call(fn_open, fn_args)
  # One row -- one chord per column, four columns, legend last.
  # Region x condition now reads left-to-right (Frontal CON, Frontal NHD, Hippo CON,
  # Hippo NHD) instead of as a 2x2 block.
  if (CHORD_MODE == "delta") {
    layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.55))
  } else {
    layout(matrix(c(1, 2, 3, 4, 5), nrow = 1), widths = c(1, 1, 1, 1, 0.42))
  }
  par(mar = c(0.1, 0.2, 0.4, 0.2), xpd = NA)
  if (CHORD_MODE == "delta") {
    for (r in REGIONS)
      draw_chord(slice_rc(r, NA), sprintf("%s  |  NHD - CON", REGION_FULL[r]), title_col = "black")
  } else {
    for (r in REGIONS) for (cd in CONDS)
      {
        draw_chord(slice_rc(r, cd), sprintf("%s  |  %s", REGION_FULL[r], cd),
                   title_col = "black")
        # The per-panel total subtitle is dropped from the plot. share
        # mode still normalises per panel, so the totals must appear in the figure legend
        # instead: Frontal CON 11.0 / NHD 11.2; Hippocampus CON 10.6 / NHD 8.8.
      }
  }
  par(mar = c(0.2, 0.05, 0.2, 0.1)); plot.new()
  # Was seq(0.82, 0.16, length.out = n), which STRETCHED the keys across the
  # full column height whatever n was -- with 4-5 pathways that left ~0.15 of blank column
  # between every swatch. Fixed pitch from the top instead, so the key block is compact and
  # the legend column can be narrowed.
  .keys <- if (CHORD_MODE == "delta")
             c("up in NHD" = "#C0392B", "down in NHD" = "#2C7FB8") else PATH_PAL[TARGET_PATHS]
  leg_y <- 0.80 - (seq_along(.keys) - 1) * 0.085
  for (i in seq_along(.keys)) {
    points(0.08, leg_y[i], pch = 22, cex = 1.9, bg = .keys[i], col = "grey30")
    text(0.24, leg_y[i], names(.keys)[i], cex = CEX_KEY, col = "black", adj = c(0, 0.5))
  }
  # Title shortened to "Pathway". After the tighten the legend column is too
  # narrow for "(ribbon colour)" -- it clipped to "(ribbon colou" -- and that parenthetical
  # is caption text, which belongs in the figure legend, not in the panel.
  text(0.02, 0.95, if (CHORD_MODE == "delta") "Change" else "Pathway", cex = CEX_HDR, font = 1,
       col = "black", adj = c(0, 1))
  do.call(fn_close, list())
}
render(pdf, list(file = out_pdf, width = W, height = H, useDingbats = FALSE), dev.off)
render(png, list(filename = out_png, width = W, height = H, units = "in", res = 600), dev.off)

cat("\nWrote:\n  ", out_pdf, "\n  ", out_png, "\n  ",
    file.path(TDIR, paste0(BN_SET, "_stats.csv")), "\n", sep = "")
cat("=== DONE ===\n", file = stderr())
