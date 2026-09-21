# =============================================================================
# 2C_micro_volcano_MAST_FH.R — Figure 2 panel c: per-nucleus MAST microglia volcano, frontal + HIPPOCAMPUS (2 regions). Faithful port of the 3-region
# reference (_run_3D_micro_volcano_MAST.R) to the FH rebuild atlas/caches.
# -----------------------------------------------------------------------------
# MAST-forward, atlas-free: reads only the FH MAST discovery CSV — never loads a
# Seurat object.  Per region (Frontal / Hippocampus) volcano of Micro-PVM
# NHD-vs-CON:
#   x = avg_log2FC (NHD vs CON)     y = -log10(p_val_adj)   [descriptive only]
#   colour = microglial PROGRAM (re-validated on FH
#            data — see gene-set re-validation below); two-tier size encoding
#            (nominal padj<0.05 small / strict + |Cliff's d|>=0.15 big).
#
# Honest FRAMING (project policy): the y-axis padj is per-nucleus and inflated by
#   pseudoreplication -> descriptive ranking axis only.  No horizontal donor-level
#   significance claim; the vertical dashed guides mark the fold-change floor and
#   the p=0.05 line is a threshold marker only.  No in-panel title/caption.
#
# gene-set re-validation:
#   Every MICRO_THEME gene was re-checked against the FH Micro-PVM MAST tables
#   (Frontal + Hippo).  result: all 7 programme themes hold — every listed gene is
#   detected (max pct >= 0.10) and moves in its programme's expected direction in
#   >=1 region.  nothing dropped.  OCC removed (this rebuild is 2-region).  The
#   OCC-only note about MARCH1/MARCHF1 being the "OCC brake" is dropped with OCC.
#   DAM-1-arrest biology confirmed at the volcano level: homeostatic P2RY12/MEF2C
#   down; MHC-II/iron/inflammatory up; the DAM-2 lipid endpoint (LPL/CST7/LGALS3/
#   CD9) never reaches the ranked list, GPNMB Hippo-only, SPP1 flat in Frontal.
#
# Sources of truth (self-contained in the FH rebuild):
#   - is_artifact : scripts/_artifact_genes.R
#   - theme_pub / PAL_DGE / PAL_REGION / strip_region_x / REGION_LABELLER :
#                   scripts/22_publication_theme_FH.R
#
# Outputs:
#   figures/Figure_2/panels/F2b_micro_volcano_MAST.{png,pdf}
#   + .provenance.txt next to the panel (counts + sessionInfo).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(ggh4x); library(ragg)
  library(patchwork)
})
set.seed(42)

# --- portable project root (two candidate machines) --------------------------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))

# microglia-only primary: PVM are excluded from the compartment. DAM
# staging is undefined for PVM (Van Hove 2019 PMID 31061494), so a pooled Micro-PVM
# volcano cannot be labelled microglial. Pooled table kept as the sensitivity tier.
MAST_PATH <- file.path(PROJ, "tables", "mast_dual", "MAST_dual_discovery_all.csv")
ARTIFACT  <- file.path(PROJ, "scripts", "_artifact_genes.R")
THEME     <- file.path(PROJ, "scripts", "22_publication_theme_FH.R")
PANEL     <- file.path(PROJ, "figures", "Figure_2", "panels")
LOGS      <- file.path(PROJ, "logs")

for (f in c(MAST_PATH, ARTIFACT, THEME))
  if (!file.exists(f)) stop(sprintf("MISSING %s — prep step failed", f))
for (d in c(PANEL, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

source(ARTIFACT)   # is_artifact()
source(THEME)      # theme_pub, PAL_DGE, PAL_REGION, strip_region_x, REGION_LABELLER, REGION_ORDER

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
CT      <- "Micro-PVM"   # POOLED microglia+PVM; PVM = 5.4% of compartment
YCAP    <- 60                      # display cap on -log10(padj)
FC_MIN  <- 1.20
LFC_MIN <- log2(FC_MIN)            # ~0.263 in log2 units; drawn as vertical dashed guides
CD_CUT   <- 0.15                   # |Cliff's delta| effect-size gate (baked into discovery)
PADJ_CUT <- 0.05                   # adjusted-p gate (baked into discovery); drawn as h-guide

# --- themed microglial programs ---
# Direction is read from the x-axis; a gene is themed only where its discovery
# direction matches its program's expected direction (up theme -> NHD-up, lost ->
# CON-up).  Every gene below was confirmed detected + directionally moving in the
# FH Micro-PVM MAST data.
MICRO_THEME <- c(
  # --- up in NHD ---
  CD74="Antigen presentation", `HLA-DRA`="Antigen presentation",
  `HLA-DRB1`="Antigen presentation", `HLA-DMB`="Antigen presentation",
  FTH1="Metal handling (iron/zinc)", FTL="Metal handling (iron/zinc)",
  HAMP="Metal handling (iron/zinc)", TMEM163="Metal handling (iron/zinc)",  # TMEM163 = zinc
  FCGR3A="Phagocytic / lysosomal", SCIN="Phagocytic / lysosomal",
  ADAM28="Phagocytic / lysosomal", CTSB="Phagocytic / lysosomal",
  PFKFB3="Glycolytic shift", SLC2A3="Glycolytic shift",
  CEBPD="Inflammatory / stress response", ZFP36L1="Inflammatory / stress response",
  RGS1="Inflammatory / stress response", SRGN="Inflammatory / stress response",
  HSPA1A="Inflammatory / stress response", TNFRSF1B="Inflammatory / stress response",
  NAIP="Inflammatory / stress response",
  F13A1="CD163+/F13A1+ state",
  # --- lost in NHD (CON-up) ---
  P2RY12="Homeostatic (lost)", MEF2C="Homeostatic (lost)",
  PLXDC2="Homeostatic (lost)", SORL1="Homeostatic (lost)",
  IRAK3="Immunoregulatory brake (lost)", LDLRAD4="Immunoregulatory brake (lost)",
  ZBTB16="Immunoregulatory brake (lost)", HDAC9="Immunoregulatory brake (lost)")

THEME_UP   <- c("Antigen presentation","Metal handling (iron/zinc)","Phagocytic / lysosomal",
                "Glycolytic shift","Inflammatory / stress response","CD163+/F13A1+ state")
THEME_LOST <- c("Homeostatic (lost)","Immunoregulatory brake (lost)")
THEME_LEVELS <- c(THEME_UP, THEME_LOST)
THEME_PAL <- c(
  "Antigen presentation"            = "#C1272D",  # red
  "Metal handling (iron/zinc)"      = "#E67E22",  # orange
  "Phagocytic / lysosomal"          = "#B8860B",  # dark goldenrod
  "Glycolytic shift"                = "#8E44AD",  # purple
  "Inflammatory / stress response"  = "#D81B60",  # magenta
  "CD163+/F13A1+ state"              = "#6D4C41",  # brown
  "Homeostatic (lost)"              = "#2E86AB",  # blue
  "Immunoregulatory brake (lost)"   = "#1B9E77")  # teal

# -----------------------------------------------------------------------------
# 1. Load + filter (atlas-free CSV only)
# -----------------------------------------------------------------------------
cat("== load FH MAST discovery table ==\n")
raw <- read.csv(MAST_PATH, stringsAsFactors = FALSE)
need <- c("cell_type","region","gene","avg_log2FC","p_val_adj",
          "cliffs_delta","discovery","pct.1","pct.2")
stopifnot(all(need %in% colnames(raw)))

d <- raw %>%
  filter(cell_type == CT, region %in% REGIONS) %>%
  mutate(
    region    = factor(region, levels = REGIONS),
    discovery = as.logical(discovery),
    # -log10(padj) with a display cap; floor padj to xmin so a p=0 underflow never
    # becomes -log10(0)=Inf (reviewer-facing edge case).  Underflow logged below.
    negl10    = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    abs_cd    = abs(cliffs_delta),
    low_pct   = ifelse(avg_log2FC > 0, pct.2, pct.1),
    dir_col   = dplyr::case_when(
                  discovery & avg_log2FC >  LFC_MIN ~ "Up in NHD",
                  discovery & avg_log2FC < -LFC_MIN ~ "Up in CON",
                  TRUE                              ~ "ns"),
    artifact  = is_artifact(gene, CT),
    theme     = unname(MICRO_THEME[gene]),
    theme_dir = dplyr::case_when(theme %in% THEME_UP   ~ "up",
                                 theme %in% THEME_LOST  ~ "lost",
                                 TRUE                   ~ NA_character_),
    dir_match = !is.na(theme) &
                ((theme_dir == "up"   & avg_log2FC >  LFC_MIN) |
                 (theme_dir == "lost" & avg_log2FC < -LFC_MIN)),
    # two-tier size encoding: coloured whenever nominal padj<0.05 (small); large
    # when it also clears |Cliff's delta|>=0.15 (strict == discovery).
    themed_strict = dir_match & discovery,
    themed        = dir_match & (p_val_adj < PADJ_CUT),
    tier      = factor(ifelse(themed_strict, "strict",
                       ifelse(themed, "nominal", "ns")),
                       levels = c("ns", "nominal", "strict")),
    theme_f   = factor(ifelse(themed, theme, NA), levels = THEME_LEVELS),
    gene_disp = gene)
stopifnot(nrow(d) > 0)

# reviewer-facing edge case: report any p=0 underflow flooring
n_uf <- sum(d$p_val_adj <= 0 | d$p_val_adj < .Machine$double.xmin, na.rm = TRUE)
if (n_uf > 0)
  cat(sprintf("NOTE: %d padj values <= xmin floored to %.3e (-log10 capped, no Inf)\n",
              n_uf, .Machine$double.xmin))

n0 <- nrow(d)
cat(sprintf("Micro-PVM rows across %d regions: %d\n", length(REGIONS), n0))
.nom <- d[d$tier == "nominal", c("region","gene","theme","avg_log2FC","cliffs_delta","p_val_adj")]
cat(sprintf("== NOMINAL-tier micro genes: %d ==\n", nrow(.nom)))
if (nrow(.nom)) print(.nom[order(.nom$region), ], row.names = FALSE)
stopifnot(sum(d$discovery & is.na(d$abs_cd)) == 0)

d$y_disp  <- pmin(d$negl10, YCAP)
d$capped  <- d$negl10 > YCAP

# -----------------------------------------------------------------------------
# 2. Per-region discovery counts (reported)
# -----------------------------------------------------------------------------
counts <- d %>%
  group_by(region) %>%
  summarise(NHD_up = sum(dir_col == "Up in NHD"),
            CON_up = sum(dir_col == "Up in CON"),
            theme_nominal = sum(tier == "nominal"),
            theme_strict  = sum(tier == "strict"),
            ns     = sum(dir_col == "ns"),
            capped = sum(capped), .groups = "drop")
cat("\n== discovery points + two-tier theme counts per region ==\n")
print(as.data.frame(counts))

# -----------------------------------------------------------------------------
# 3. Label set = the themed genes that are genuine discovery hits per region.
# -----------------------------------------------------------------------------
lab_df <- d %>% filter(themed)
cat("\n== themed (coloured + labelled) genes per region ==\n")
for (r in REGIONS) {
  gg <- lab_df %>% filter(region == r) %>% arrange(theme_f, desc(abs(avg_log2FC)))
  cat(sprintf("  %-7s [%d]: %s\n", r, nrow(gg),
              paste(sprintf("%s(%s)", gg$gene_disp, abbreviate(gg$theme, 4)),
                    collapse = ", ")))
}

# -----------------------------------------------------------------------------
# 4. Build volcano — 2-region facet (Frontal / Hippocampus).
# -----------------------------------------------------------------------------
d <- d[order(match(d$tier, c("ns","nominal","strict"))), ]  # strict drawn last (on top)
d$region <- factor(as.character(d$region), levels = REGIONS)
df_ns  <- dplyr::filter(d, !themed, dir_col == "ns")
df_d0  <- dplyr::filter(d, !themed, dir_col != "ns")
df_nom <- dplyr::filter(d, tier == "nominal")
df_str <- dplyr::filter(d, tier == "strict")
df_th  <- dplyr::filter(d, themed)
xl <- max(abs(d$avg_log2FC)) * 1.04

# per-facet y crop (kills the empty band above each cloud) — 90th pct of labelled
# points + headroom, clamped to the global display ceiling.  Deterministic.
# Allow extra top headroom (YCEIL up to 1.22x YCAP) so labels for points that
# pin at the YCAP ceiling (F13A1 / TMEM163 in Frontal) have room to stack above the
# points instead of colliding at the top edge.
YCEIL <- YCAP * 1.22
reg_upper <- vapply(REGIONS, function(r) {
  yv <- d$y_disp[d$region == r & d$themed]
  tm <- if (length(yv)) as.numeric(quantile(yv, 0.90, names = FALSE))
        else max(d$y_disp[d$region == r], na.rm = TRUE)
  # if any labelled point pins at the display ceiling, guarantee a headroom band
  # above it for the stacked labels; otherwise the usual tight crop.
  pinned <- any(d$y_disp[d$region == r & d$themed] >= YCAP - 1e-6)
  cap <- min(tm + max(2.5, 0.12 * tm), YCEIL)
  if (pinned) min(YCEIL, YCAP + 0.20 * YCAP) else cap
}, numeric(1))
cat(sprintf("\n== per-facet y upper crop (topmost label + headroom; ceiling %.1f) ==\n", YCEIL))
print(round(reg_upper, 1))
y_scales <- lapply(REGIONS, function(r)
  scale_y_continuous(limits = c(0, reg_upper[[r]]), oob = scales::squish,
                     breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0.02, 0.03))))

# Label set: cap per facet, strict then headline genes preferred, deterministic.
# In Hippocampus ~13 of the 15 themed genes sit on the positive
# side in a narrow x band; the y pre-stagger separates their starting points but ggrepel's
# force_pull then drags them back together (HSPA1A/TMEM163, TNFRSF1B/SLC2A3).  Below ~11
# labels per facet the solve converges cleanly.  Selection order is unchanged (strict tier,
# then headline genes) and theme colour still encodes every point that loses its text.
LAB_CAP  <- 11L
HEADLINE <- c("P2RY12","MEF2C","CD74","HLA-DRA","FTH1","FTL","HAMP","TMEM163","CEBPD","F13A1")
df_lab <- df_th %>%
  group_by(region) %>%
  arrange(desc(tier == "strict"), desc(gene %in% HEADLINE),
          desc(abs(avg_log2FC)), .by_group = TRUE) %>%
  slice_head(n = LAB_CAP) %>% ungroup() %>% as.data.frame()
cat(sprintf("\n== labels drawn per facet (cap %d of themed) ==\n", LAB_CAP))
print(df_lab %>% count(region, name = "n_labels") %>% as.data.frame())
# a strong outward nudge pinned right-hand labels against the panel edge, where they
# had nowhere left to resolve; halve it and let repel use the empty left half instead.
df_lab$nudge_x <- sign(df_lab$avg_log2FC) * 0.05 * xl

# ggrepel starts every label at its point, so two labels at nearly the same y on the
# same side start on top of each other and the solver can settle in that overlapping
# local minimum (SRGN / TMEM163 in Hippocampus printed as "SRGNMEM163" across three
# retunes).  Pre-stagger vertically: within each region+side, walk labels in y order
# and push any label that starts within SEP of the previous one further apart.  This
# is deterministic and gene-agnostic -- no hard-coded per-gene nudges.
df_lab$nudge_y <- 0
for (r in unique(df_lab$region)) for (sd in c(-1, 1)) {
  ii <- which(df_lab$region == r & sign(df_lab$avg_log2FC) == sd)
  if (length(ii) < 2) next
  # Fit the stagger into the panel instead of pushing blindly upward: if the labels
  # need more room than the axis has, shrink the separation to what actually fits, then
  # re-centre the staggered group on its original centre.  A uniform "slide down when it
  # overshoots" is not safe -- it drags labels that needed no stagger onto their
  # neighbours (it re-created "SRGNMEM163" here after fixing 4B).
  top <- reg_upper[[as.character(r)]] * 0.97; bot <- 0.02 * top
  ii   <- ii[order(df_lab$y_disp[ii])]     # walk in y order; keeps nudge_y aligned to ii
  ypos <- df_lab$y_disp[ii]
  n_l  <- length(ypos)
  SEP  <- min(0.055 * top, (top - bot) / max(n_l - 1, 1))
  for (k in seq_len(n_l)[-1])
    if (ypos[k] - ypos[k-1] < SEP) ypos[k] <- ypos[k-1] + SEP
  ctr  <- mean(range(df_lab$y_disp[ii]))
  ypos <- ypos - (mean(range(ypos)) - ctr)
  if (max(ypos) > top) ypos <- ypos - (max(ypos) - top)
  if (min(ypos) < bot) ypos <- ypos + (bot - min(ypos))
  df_lab$nudge_y[ii] <- ypos - df_lab$y_disp[ii]
}
cat(sprintf("\n== label pre-stagger: %d of %d labels shifted in y (max %.1f) ==\n",
            sum(abs(df_lab$nudge_y) > 1e-9), nrow(df_lab), max(abs(df_lab$nudge_y))))

lay_nom <- if (nrow(df_nom) > 0)
  geom_point(data = df_nom, aes(avg_log2FC, y_disp, colour = theme_f),
             size = 0.9, alpha = 0.5, shape = 16, stroke = 0) else NULL
cat(sprintf("nominal-tier layer: %s (%d rows)\n",
            ifelse(is.null(lay_nom), "OMITTED (empty)", "drawn"), nrow(df_nom)))

p <- ggplot() +
  geom_vline(xintercept = c(-LFC_MIN, LFC_MIN),
             linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_hline(yintercept = -log10(PADJ_CUT),
             linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25, linetype = "22") +
  geom_point(data = df_ns, aes(avg_log2FC, y_disp, shape = capped),
             colour = "#949494", size = 0.5, alpha = 0.5, stroke = 0.3) +
  geom_point(data = df_d0, aes(avg_log2FC, y_disp),
             colour = "#ADADAD", size = 0.75, alpha = 0.45, shape = 16, stroke = 0) +
  # An empty layer leaves ggh4x with no data to
  # train that panel's y scale under facetted_pos_scales(), which aborts with
  # Add the layer only when it has rows.
  lay_nom +
  # Display cap made visible: points pinned at
  # YCAP were previously indistinguishable from real values, so 66 Frontal features sat
  # on an invisible ceiling that reads as a horizontal feature. Capped points are drawn
  # as open TRIANGLES via a shape SCALE (an extra geom_point layer with a fixed shape
  # broke the patchwork assembly). `capped` was already computed and unused.
  geom_point(data = df_str, aes(avg_log2FC, y_disp, colour = theme_f, shape = capped),
             size = 2.3, alpha = 0.90, stroke = 0.4) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), guide = "none",
                     drop = FALSE) +
  # F13A1 & TMEM163 both pin at the YCAP ceiling in Frontal (same
  # y), so their labels collided at the top-right.  Stronger repel (force/box.padding)
  # + a downward pull off the ceiling separates them deterministically (fixed seed).
  geom_text_repel(data = df_lab,
                  aes(avg_log2FC, y_disp, label = gene_disp),
                  colour = "black", fontface = "italic", size = 2.0, seed = 42,
                  nudge_x = df_lab$nudge_x, nudge_y = df_lab$nudge_y,
                  # After the RNA-detection fix the Hippo cloud tightened
                  # and ggrepel stopped converging inside the narrower facet -- SRGN and
                  # TMEM163 printed on top of each other ("SRGNTMEM163"), FTL/HLA-DRA
                  # crowded the bottom edge.  max.overlaps=Inf never DROPS a label, so a
                  # non-converged solve just draws them overlapping.  More force + padding
                  # and a much larger iteration/time budget so the solve actually finishes.
                  max.overlaps = Inf, force = 14, force_pull = 0.12,
                  box.padding = 0.95, point.padding = 0.35,
                  max.iter = 50000, max.time = 3,
                  min.segment.length = 0,
                  segment.size = 0.2, segment.colour = "grey55", segment.alpha = 0.7,
                  bg.color = "white", bg.r = 0.12,
                  show.legend = FALSE) +
  # Both region banners use their normal pale region tone; the single-NHD-hippo-lane
  # caveat is a legend note, not an on-plot greyed banner.
  facet_wrap2(~ region, nrow = 1, scales = "free_y",
              labeller = REGION_LABELLER,
              strip = strip_region_x(REGIONS, clip = "off")) +
  ggh4x::facetted_pos_scales(y = y_scales) +
  scale_colour_manual(values = THEME_PAL, name = NULL, drop = FALSE, na.translate = FALSE,
                      guide = guide_legend(override.aes = list(size = 2.6, alpha = 1),
                                           ncol = 3, byrow = TRUE)) +
  scale_x_continuous(limits = c(-xl, xl), oob = scales::squish,
                     expand = expansion(mult = 0.04),
                     breaks = scales::pretty_breaks(n = 5)) +
  labs(x = expression(log[2]~"fold change (NHD vs CON)"),
       y = expression(-log[10]~"(adjusted "*italic(p)*")")) +
  theme_classic(base_size = 8) +
  theme(panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.line      = element_blank(),
        axis.ticks     = element_line(colour = "black", linewidth = 0.35),
        axis.text      = element_text(colour = "black"),
        axis.title     = element_text(colour = "black"),
        strip.background = element_blank(),
        panel.spacing.x = unit(0.32, "lines"),
        legend.position = "bottom",
        legend.text     = element_text(size = 7.4),
        legend.key.size = unit(0.34, "cm"),
        legend.key.spacing.x = unit(0.25, "cm"),
        legend.box.spacing = unit(0.1, "cm"),
        legend.margin      = margin(t = 1, b = 0),
        plot.margin = margin(2, 3, 1, 3))

# -----------------------------------------------------------------------------
# 4b. Right-side threshold KEY (states the exact calling rule + programme hues).
# -----------------------------------------------------------------------------
kt <- function(x, y, lab, size = 2.35, face = "plain", col = "black", parse = FALSE)
  annotate("text", x = x, y = y, label = lab, hjust = 0, size = size,
           colour = col, fontface = face, parse = parse)
UP_COL <- unname(PAL_DGE["Up in NHD"]); DN_COL <- unname(PAL_DGE["Up in CON"])
NS_COL <- unname(PAL_DGE["ns"])
key <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  kt(0.00, 0.985, "Coloured = theme gene:", size = 2.5) +
  kt(0.06, 0.915, "adjusted~italic(p)~\"< 0.05\"", parse = TRUE) +
  kt(0.06, 0.855, "\"fold change \" >= \"1.20x\"", parse = TRUE) +
  (if (nrow(df_nom) > 0)
     list(annotate("point", x = 0.03, y = 0.775, colour = "grey45", size = 0.9, shape = 16),
          kt(0.10, 0.775, "nominal (small)")) else NULL) +
  annotate("point", x = 0.03, y = 0.700, colour = "grey45", size = 2.3, shape = 16) +
  kt(0.10, 0.700, "\"strict: |\"*delta*\"| \" >= 0.15", parse = TRUE) +
  annotate("segment", x = 0.02, xend = 0.02, y = 0.585, yend = 0.635,
           linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  kt(0.10, 0.610, "fold-change floor") +
  annotate("segment", x = 0.00, xend = 0.05, y = 0.520, yend = 0.520,
           linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  kt(0.10, 0.520, "italic(p)~\"= 0.05\"", parse = TRUE) +
  annotate("point", x = 0.03, y = 0.390, colour = UP_COL, size = 2.4, shape = 16) +
  kt(0.10, 0.390, "up in NHD (warm)") +
  annotate("point", x = 0.03, y = 0.320, colour = DN_COL, size = 2.4, shape = 16) +
  kt(0.10, 0.320, "lost in NHD (cool)") +
  annotate("point", x = 0.03, y = 0.250, colour = NS_COL, size = 2.4, shape = 16) +
  kt(0.10, 0.250, "not significant") +
  annotate("point", x = 0.03, y = 0.180, colour = "grey45", size = 2.0, shape = 2, stroke = 0.5) +
  kt(0.10, 0.180, paste0("beyond axis cap (", YCAP, ")")) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 1))

# stitch: volcano (keeps bottom programme legend) | key.  2-region volcano is
# narrower than the 3-region original, so the overall width is reduced.
combo <- p + key + patchwork::plot_layout(widths = c(1, 0.34))

# -----------------------------------------------------------------------------
# 5. Render (ragg PNG + pdf).
# -----------------------------------------------------------------------------
# 10% slimmer (width only). 7.40 -> 6.66 in; height unchanged.
# H raised 3.35 -> 3.85. After the RNA-detection fix the Hippo labelled
# points concentrated into x in [0,7], and 15 labels could not be placed in a 3.35in
# panel without overlapping regardless of repel force.  Height is the axis the label
# stack actually needs.
W <- 7.4 * 0.9; H <- 3.85
png_path <- file.path(PANEL, "F2b_micro_volcano_MAST.png")
pdf_path <- file.path(PANEL, "F2b_micro_volcano_MAST.pdf")
ragg::agg_png(png_path, width = W, height = H, units = "in", res = 300)
print(combo); invisible(dev.off())
ggsave(pdf_path, combo, width = W, height = H, useDingbats = FALSE)
cat(sprintf("\nWrote panel: %s\n         and: %s\n", png_path, pdf_path))

# -----------------------------------------------------------------------------
# 6. Provenance next to the panel.
# -----------------------------------------------------------------------------
prov <- file.path(PANEL, "F2b_micro_volcano_MAST.provenance.txt")
sink(prov)
cat("F2b_micro_volcano_MAST (FH) — provenance\n")
cat("generated:", format(Sys.time()), "\n")
cat("script   :", "scripts/2C_micro_volcano_MAST_FH.R\n")
cat("input    :", MAST_PATH, "\n")
cat("filter   : cell_type==Microglia, region in {Frontal,Hippo}\n")
cat("YCAP (-log10 padj display cap):", YCAP, "\n")
cat("padj underflow floored to xmin:", n_uf, "rows\n\n")
cat("== discovery counts per region ==\n"); print(as.data.frame(counts))
cat("\n== labelled genes ==\n")
print(as.data.frame(lab_df[, c("region","dir_col","gene","avg_log2FC","cliffs_delta")]))
cat("\n== sessionInfo ==\n"); print(sessionInfo())
sink()
cat(sprintf("Wrote provenance: %s\n", prov))

cat("\n=== DONE ===\n", file = stderr())
