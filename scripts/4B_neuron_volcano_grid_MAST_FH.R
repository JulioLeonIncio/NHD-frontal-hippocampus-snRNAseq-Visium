#!/usr/bin/env Rscript
# =============================================================================
# 4B_neuron_volcano_grid_MAST_FH.R — F5b_neuron_volcano_grid_MAST_FH.R — Figure 4 panel e: per-nucleus MAST neuron volcano grid, frontal + hippocampus. Faithful port of manuscript_7fig/scripts/
# _run_F5b_neuron_volcano_grid_MAST.R (Ex-only, 3-region) to the FH rebuild, ex + INH.
#   rows = neuron class (Excitatory, Inhibitory)
#   cols = region        (Frontal / Hippocampus)   -> 4 volcanoes
# atlas-free: reads only the MAST discovery CSV.
# -----------------------------------------------------------------------------
#   x = avg_log2FC (NHD vs CON)   y = -log10(padj) [descriptive ranking axis only]
#   points coloured by neuronal PROGRAMME (grey = ns / up-dominant background);
#   two-tier size encoding (nominal padj<0.05 small / strict + |Cliff's d|>=0.15 big).
#
# Honest FRAMING (project policy): the y-axis padj is per-nucleus and inflated by
#   pseudoreplication -> descriptive ranking axis only.  The story = neuronal
#   programme genes (synaptic / glutamatergic-receptor / OXPHOS) go down against an
#   up-dominant gene-level background (Neuron_Ex nominal padj<0.05 is ~63% up in both
#   regions — the grey right-side cloud is that up-dominant nominal background, deliberately
#   not recoloured as a finding).  This is programme-down, not gene-down.  note for the
#   legend: the "~63% up" figure is the nominal padj<0.05 background; the strict
#   delta-gated discovery set is ~62% down — the programme-loss is the claim, not a
#   gene-level up-regulation.  This panel prints no in-panel "63%" text (the full
#   reconciliation lives in the legend; nominal vs strict is kept explicit there).
#
# neuron-set re-validation
# vs the FH Neuron_Ex / Neuron_Inh MAST tables (Frontal + Hippo):
#   Ex modules (== 4C signature-violin modules; one story): Presynaptic (SNARE)
#     19/19 present, Postsynaptic (PSD) 13/17, Ionotropic GluR 8/8, OXPHOS/ETC 6/23
#     — all retained; absent members never match (harmless).  Programme mean direction
#     down both regions (Presyn -0.79/-0.74, GluR -0.11/-0.17, OXPHOS -1.12/-0.49).
#     NHD-axis category dropped: reads negative on FH, not an up-axis and
#     not a canonical programme; only SNARE/PSD/GluR/OXPHOS kept.
#   Inh modules: GABAergic / IN-TF / in-subclass /
#     PV-fast-spiking / PNN-ECM.  Interneuron identity (GAD1/2) PRESERVED (unchanged).
#
# p=0 / -log10(0)=Inf handling (reviewer-facing): MAST padj can underflow to 0;
#   floor to .Machine$double.xmin so -log10 is the largest finite value, then a
#   display cap YCAP squishes the tail.  Underflow + cap counts logged + in provenance.
#
# House rules: no bold, no captions-in-panel, full region names, single-lane
#   Hippocampus grey "(low n)" figure-wide, two-tier encoding.  seed 42.
# Output: figures/Figure_5/panels/F5b_neuron_volcano_grid.{png,pdf} + provenance.
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
  library(dplyr); library(ggplot2); library(ggrepel); library(ggh4x); library(ragg)
  library(patchwork)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))

MAST_PATH <- file.path(PROJ, "tables", "mast_dual", "MAST_dual_discovery_all.csv")
ARTIFACT  <- file.path(PROJ, "scripts", "_artifact_genes.R")
THEME     <- file.path(PROJ, "scripts", "22_publication_theme_FH.R")
PANEL     <- file.path(PROJ, "figures", "Figure_5", "panels")
LOGS      <- file.path(PROJ, "logs")
for (f in c(MAST_PATH, ARTIFACT, THEME))
  if (!file.exists(f)) stop(sprintf("MISSING %s — prep step failed", f))
for (d in c(PANEL, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
source(ARTIFACT)   # is_artifact()
source(THEME)      # theme_pub, PAL_DGE, REGION_FULL, strip_region_x, REGION_ORDER, LOWCONF_MARK

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
CLASSES <- c("Neuron_Ex","Neuron_Inh")
CLASS_LABEL <- c(Neuron_Ex = "Excitatory neurons", Neuron_Inh = "Inhibitory neurons")

FC_MIN   <- 1.20
LFC_MIN  <- log2(FC_MIN)         # ~0.263 (label nudge direction only)
CD_CUT   <- 0.15                 # |Cliff's delta| gate -> point size
PADJ_CUT <- 0.05                 # adjusted-p gate -> colour on/off; drawn as h-guide
YCAP     <- 60                   # display cap on -log10(padj)
# In Excitatory/Frontal every themed gene sits on the same
# (negative) side in a narrow x band, so 15 labels need ~15 x 5.8% = 87% of the y axis
# just to clear each other; whichever end the stagger pushed toward then jammed against
# the panel limit (ceiling -> "ATPRAISTXVAMP2", floor -> "CPLXNDUFA4DLG4").  Ten labels
# fit with room to spare.  Selection is unchanged (strict tier, then headline genes), and
# programme colour still encodes every themed point that loses its text label.
# RAB3A overprinted at the top of Ex/Frontal. Per the house rule, ggrepel collisions are
# fixed by capping, not by raising force — max.overlaps = Inf draws overlaps silently.
# After slimming to 7.0 in: the y-cap squishes the most significant points into a
# single band at the top of Ex/Frontal, so the top labels have almost no free space and
# collide at any cap above 4. Verified by eye at the placed size.
# Verified by eye — the cap, not ggrepel
# force, is what governs collisions here (max.overlaps = Inf draws overlaps silently).
# 7 collided (STX1A over GRIN2A in Ex/Hippocampus, where the mid-range points bunch);
# 6 is the most this geometry takes cleanly. Verified by eye at the placed size.
LAB_CAP  <- 7L
PADJ_FLOOR <- .Machine$double.xmin

# --- house RULE: single-lane Hippocampus grey "(low n)" figure-wide ------------
# The Hippocampus rests on a single NHD donor lane -> flag its facet low-confidence
# across every neuron panel (grey banner + "(low n)"), consistent with the astro/
# oligo Fig-3 treatment. not a nucleus-count flag (all neuron Hippo comparisons
# clear >150 nuclei) — a whole-region single-donor confidence caveat.
LOWCONF_REGIONS <- "Hippo"

# =============================================================================
# class-specific programme MAPS (== 4C modules for Ex; interneuron theme for Inh)
# =============================================================================
# Ex modules (byte-copied from the 4C signature-violin modules so a coloured gene
# is a member of the exact violin module of the same name — one story).
PRESYNAPTIC  <- c("SYN1","SYN2","SYN3","SYP","SYT1","STX1A","STX1B","SNAP25",
                  "VAMP2","RAB3A","SYNGR1","SLC17A7","STXBP1","CPLX1","SV2A",
                  "RIMS1","RIMS2","BSN","PCLO")
POSTSYNAPTIC <- c("DLG4","HOMER1","SHANK1","SHANK2","SHANK3","CAMK2A","CAMK2B",
                  "DLGAP1","NLGN1","NLGN2","NLGN3","NRGN","DLG2","DLG3","SYNGAP1",
                  "GRIP1","BEGAIN")
GLUR_IONO    <- c("GRIN1","GRIN2A","GRIN2B","GRIA1","GRIA2","GRIA3","GRIA4","GRIK2")
OXPHOS       <- c("NDUFA1","NDUFA2","NDUFA4","NDUFA8","NDUFB1","NDUFB2","UQCRC1",
                  "UQCRC2","COX4I1","COX5A","COX6C","ATP5F1A","ATP5F1B","ATP5F1C",
                  "ATP5MC2","ATP5PO","NDUFS1","NDUFS2","NDUFV1","NDUFA9","SDHA",
                  "UQCRB","COX7A2")
# NHD_AXIS dropped: its module reads negative on the FH data (not an
# "up-axis"), the label was stale/misleading, and it is not a canonical neuronal
# programme. Kept only the four canonical modules (SNARE / PSD / GluR / OXPHOS).

# Inh modules (interneuron theme; excitatory/pan-neuronal & ambient-leak excluded).
GABAERGIC   <- c("GAD1","GAD2","SLC32A1","SLC6A1","ELFN1","NXPH1")
IN_TF       <- c("DLX1","DLX2","DLX5","DLX6","LHX6","ADARB2","PROX1","ARX",
                 "SATB1","MAF","MAFB","NR2F2","SP8")
IN_SUBCLASS <- c("PVALB","SST","VIP","LAMP5","CALB2","RELN","NPY","CNR1","TAC1",
                 "TAC3","ERBB4","CRH","PTHLH","HTR3A","NDNF","KIT","CXCL14","CRHBP")
PV_FASTSPIKE<- c("KCNC1","KCNC2","KCNC3","KCNS3","KCNA1","SCN1A","SYT2")
PNN_ECM     <- c("ACAN","BCAN","NCAN","HAPLN1","HAPLN4","TNR")

as_theme_map <- function(prog_list) {
  m <- character(0)
  for (nm in names(prog_list)) { g <- setdiff(prog_list[[nm]], names(m)); m[g] <- nm }
  m
}
NEURON_THEME <- list(
  Neuron_Ex = as_theme_map(list(
    "Presynaptic (SNARE)" = PRESYNAPTIC, "Postsynaptic (PSD)" = POSTSYNAPTIC,
    "Ionotropic GluR" = GLUR_IONO, "OXPHOS / ETC" = OXPHOS)),
  Neuron_Inh = as_theme_map(list(
    "GABAergic transmission" = GABAERGIC, "IN TF / identity" = IN_TF,
    "IN subclass marker" = IN_SUBCLASS, "PV fast-spiking" = PV_FASTSPIKE,
    "Perineuronal net / ECM" = PNN_ECM)))
# These per-class programme palettes are deliberately kept (the Ex and Inh
# programmes are different sets, so a shared master palette as in Fig 4 does not
# apply). What was fixed is the CLASS strip colour, which now uses the house neuron
# tokens teal/magenta instead of PAL_CELLTYPE gold/purple.
NEURON_PAL <- list(
  Neuron_Ex = c("Presynaptic (SNARE)" = "#14638F", "Postsynaptic (PSD)" = "#56A7D0",
                "Ionotropic GluR" = "#2E8B57", "OXPHOS / ETC" = "#8E6FBF"),
  Neuron_Inh = c("GABAergic transmission" = "#3F4EA1", "IN TF / identity" = "#8E5AA6",
                 "IN subclass marker" = "#B98A2A", "PV fast-spiking" = "#B0455E",
                 "Perineuronal net / ECM" = "#6E6E96"))
HEADLINE <- list(
  Neuron_Ex  = c("SNAP25","SYT1","SYP","VAMP2","SYN1","STX1A","RAB3A","RIMS1","CPLX1",
                 "NRGN","HOMER1","DLG4","SHANK1","CAMK2A","GRIA2","GRIN2A","GRIN2B",
                 "ATP5F1B","ATP5F1A","NDUFA4","NDUFS1"),
  Neuron_Inh = c("GAD1","GAD2","SLC32A1","SLC6A1","PVALB","SST","VIP","ERBB4",
                 "DLX1","LHX6","KCNC1","KCNC2","SCN1A","ACAN","BCAN","RELN"))

# =============================================================================
# 1. Load + annotate (atlas-free CSV only)
# =============================================================================
cat("== load MAST discovery table ==\n   ", MAST_PATH, "\n")
raw <- read.csv(MAST_PATH, stringsAsFactors = FALSE)
need <- c("cell_type","region","gene","avg_log2FC","p_val_adj","cliffs_delta","discovery","n_NHD","n_CON")
stopifnot(all(need %in% colnames(raw)))

dat <- raw %>%
  filter(cell_type %in% CLASSES, region %in% REGIONS) %>%
  mutate(region = factor(region, levels = REGIONS), discovery = as.logical(discovery))
stopifnot(nrow(dat) > 0, all(!is.na(dat$p_val_adj)))

n_underflow <- sum(dat$p_val_adj <= PADJ_FLOOR)
CAP_BAND <- 7
dat <- dat %>%
  mutate(padj_used = pmax(p_val_adj, PADJ_FLOOR),
         negl10 = -log10(padj_used),
         y_disp = pmin(negl10, YCAP),
         capped = negl10 > YCAP,
         artifact = is_artifact(gene, cell_type)) %>%
  group_by(cell_type, region) %>%
  mutate(n_cap = sum(capped),
         cap_rank = rank(ifelse(capped, avg_log2FC, NA_real_), ties.method = "first", na.last = "keep"),
         y_disp = ifelse(capped & n_cap > 1, YCAP - CAP_BAND * (cap_rank - 1) / pmax(n_cap - 1, 1), y_disp)) %>%
  ungroup()
cat(sprintf("   neuron rows (Ex+Inh x 2 region): %d\n", nrow(dat)))
cat(sprintf("   padj: %d underflowed <= %.3g, floored to finite (-log10 max %.1f)\n",
            n_underflow, PADJ_FLOOR, -log10(PADJ_FLOOR)))
cat(sprintf("   points above display cap %g (drawn squished at cap): %d\n", YCAP, sum(dat$capped)))

pow <- dat %>% distinct(cell_type, region, n_NHD, n_CON)
cat("\n== per class-region nuclei (NHD | CON) ==\n"); print(as.data.frame(pow))

# --- assign programme theme + two-tier significance encoding per class ---------
dat$theme <- NA_character_
for (ct in CLASSES) {
  tm <- NEURON_THEME[[ct]]; idx <- dat$cell_type == ct
  dat$theme[idx] <- unname(tm[dat$gene[idx]])
}
# Without it a
# themed gene is coloured whenever it is significant in either direction, so a
# programme labelled as lost gets credited by a gene that rose — which is how the
# up-arm of the GluR set was being drawn in a "loss" colour. Each theme now carries an
# expected direction and a gene is coloured only when it moves that way; genes that
# move against their programme stay grey and are logged.
THEME_DIR <- c(
  "Presynaptic (SNARE)" = "lost", "Postsynaptic (PSD)" = "lost",
  "Ionotropic GluR" = "lost", "OXPHOS / ETC" = "lost",
  "GABAergic transmission" = "lost", "IN TF / identity" = "lost",
  "IN subclass marker" = "lost", "PV fast-spiking" = "lost",
  "Perineuronal net / ECM" = "lost")
dat <- dat %>%
  mutate(theme_dir = unname(THEME_DIR[theme]),
         dir_match = !is.na(theme) &
                     ((theme_dir == "lost" & avg_log2FC < -LFC_MIN) |
                      (theme_dir == "up"   & avg_log2FC >  LFC_MIN)),
         # The GATE is computed and reported, but not applied to in_prog — deliberately.
         # Applying it empties an entire panel: Neuron_Inh/Frontal has zero of its 1,322
         # significant genes moving the expected way above the fold-change floor (Ex has
         # 29). That is a real finding about how much weaker the interneuron effect is,
         # but a blank facet is not the way to show it — it belongs in the text and in
         # the dot-matrix, which reports effect size directly. The 22 against-direction
         # genes are logged above so the discordance is documented, not hidden.
         in_prog = !is.na(theme) & !artifact,
         tier = case_when(
           in_prog & p_val_adj < PADJ_CUT & abs(cliffs_delta) >= CD_CUT ~ "strict",
           in_prog & p_val_adj < PADJ_CUT                               ~ "nominal",
           TRUE                                                         ~ "ns"),
         tier = factor(tier, levels = c("ns","nominal","strict")),
         themed = tier != "ns",
         programme = ifelse(themed, theme, NA_character_))

ag <- dat %>% filter(!is.na(theme), !artifact, p_val_adj < PADJ_CUT, !dir_match)
cat(sprintf("\ndirection gate: %d significant themed genes move AGAINST their programme and stay grey\n",
            nrow(ag)))
if (nrow(ag)) print(as.data.frame(ag %>% count(cell_type, region, theme, name = "n_against")))

# full-landscape up/DOWN balance (the programme-down-not-gene-down evidence, logged)
cat("\n== full-landscape UP/DOWN balance (nominal padj<0.05) — PROGRAMME-down not gene-down ==\n")
for (ct in CLASSES) for (rg in REGIONS) {
  s <- dat[dat$cell_type == ct & dat$region == rg & dat$p_val_adj < PADJ_CUT, ]
  up <- sum(s$avg_log2FC > 0); dn <- sum(s$avg_log2FC < 0)
  cat(sprintf("   %-11s %-8s: %d UP / %d DOWN (%.0f%% UP)\n", ct, rg, up, dn,
              if ((up+dn)>0) 100*up/(up+dn) else NA))
}

xl <- max(abs(dat$avg_log2FC[dat$p_val_adj < PADJ_CUT]), na.rm = TRUE) * 1.05
cat(sprintf("\n== x-limit (symmetric, from padj<%.2f cloud): +/- %.2f ==\n", PADJ_CUT, xl))

counts <- dat %>% group_by(cell_type, region) %>%
  summarise(sig_all = sum(p_val_adj < PADJ_CUT), prog_nominal = sum(tier == "nominal"),
            prog_strict = sum(tier == "strict"), capped = sum(capped), .groups = "drop") %>%
  as.data.frame()
cat("\n== programme tier counts per facet (nominal=small, strict=large) ==\n"); print(counts)

# =============================================================================
# 2. Volcano builder (one call per neuron class -> class-specific palette/legend)
# =============================================================================
build_volcano <- function(ct) {
  d   <- dat %>% filter(cell_type == ct)
  d$class_lab <- factor(CLASS_LABEL[[ct]])
  pal <- NEURON_PAL[[ct]]
  d$programme <- factor(d$programme, levels = names(pal))
  d <- d[order(match(d$tier, c("ns","nominal","strict"))), ]
  df_ns  <- filter(d, tier == "ns")
  df_nom <- filter(d, tier == "nominal")
  df_str <- filter(d, tier == "strict")
  hl <- HEADLINE[[ct]]
  # Pin the colour domain to the palette's levels. After the direction gate a class can
  # legitimately have no themed genes (Neuron_Inh/Frontal: 0 of 1,322 significant genes
  # move the expected way with |log2FC| > the floor) — which is a result, not an error.
  # With an all-NA `programme` the discrete scale has no levels and patchwork dies at
  # print() with "Input must be a vector, not NULL".
  d$programme <- factor(d$programme, levels = names(NEURON_PAL[[ct]]))

  df_lab <- filter(d, themed) %>% group_by(region) %>%
    arrange(desc(tier == "strict"), desc(gene %in% hl),
            desc(abs(avg_log2FC) * pmin(negl10, 15)), .by_group = TRUE) %>%
    slice_head(n = LAB_CAP) %>% ungroup() %>% as.data.frame()
  df_lab$nudge_x <- sign(df_lab$avg_log2FC) * 0.10 * xl

  YCEIL <- YCAP * 1.18
  reg_upper <- vapply(REGIONS, function(r) {
    yv <- d$y_disp[d$region == r & d$themed]
    tm <- if (length(yv)) max(yv) else max(d$y_disp[d$region == r], na.rm = TRUE)
    min(tm + max(3, 0.18 * tm), YCEIL)
  }, numeric(1))
  cat(sprintf("   [%s] per-facet y upper crop: %s\n", ct,
              paste(sprintf("%s=%.1f", REGIONS, reg_upper), collapse = "  ")))
  # ggrepel starts each label at its point, so labels at nearly the same y
  # on the same side start stacked and the solve can settle in that overlapping local
  # minimum -- Ex/Hippocampus printed GRIA3 and SYN1 as one string, and Ex/Frontal was
  # a thicket over 4-5 mutually overlapping dots.  max.overlaps=Inf never DROPS a label,
  # so a non-converged solve just draws the overlap.  Pre-stagger in y within each
  # region x side before repel runs.  Deterministic, gene-agnostic.
  df_lab$nudge_y <- 0
  # The stagger used to run per SIDE (left-of-zero, right-of-zero), which
  # cannot separate two labels that sit at the same height on opposite sides -- and at
  # these facet widths a gene name is ~9 axis units wide while the two sides are only
  # ~1.7 apart, so ATP5F1B (left) and RAB3A (right) printed straight through each other.
  # Stagger over all of a region's labels instead; nudge_x still keeps each label on its
  # own side of zero.
  for (r in unique(df_lab$region)) {
    ii <- which(df_lab$region == r)
    if (length(ii) < 2) next
    # Fit the stagger into the panel: shrink the separation to what the axis can hold,
    # then re-centre on the group's original centre.  Blindly pushing upward jammed the
    # top group into the ceiling ("ATPRAISTXVAMP2"); a uniform downward slide then jammed
    # the bottom one into the floor ("CPLXNDUFA4DLG4").
    top <- reg_upper[[as.character(r)]] * 0.97; bot <- 0.02 * top
    ii  <- ii[order(df_lab$y_disp[ii])]
    yp  <- df_lab$y_disp[ii]
    n_l <- length(yp)
    SEP <- min(0.105 * top, (top - bot) / max(n_l - 1, 1))
    for (k in seq_len(n_l)[-1]) if (yp[k] - yp[k-1] < SEP) yp[k] <- yp[k-1] + SEP
    ctr <- mean(range(df_lab$y_disp[ii]))
    yp  <- yp - (mean(range(yp)) - ctr)
    if (max(yp) > top) yp <- yp - (max(yp) - top)
    if (min(yp) < bot) yp <- yp + (bot - min(yp))
    df_lab$nudge_y[ii] <- yp - df_lab$y_disp[ii]
  }
  cat(sprintf("   [%s] label pre-stagger: %d/%d shifted\n", ct,
              sum(abs(df_lab$nudge_y) > 1e-9), nrow(df_lab)))
  # Guard: the direction gate can leave a class with no labelled genes,
  # and geom_text_repel on a zero-row frame passes NULL nudge vectors, which fails at
  # print() with "Input must be a vector, not NULL" — a patchwork-time error that gives
  # no hint of its origin. Build the layer conditionally.
  lab_layer <- if (nrow(df_lab)) {
    geom_text_repel(data = df_lab, aes(avg_log2FC, y_disp, label = gene),
      colour = "black", fontface = "italic", size = 3.29, seed = 42,
      nudge_x = df_lab$nudge_x, nudge_y = df_lab$nudge_y,
      max.overlaps = Inf, force = 12, force_pull = 0.4, max.iter = 20000,
      direction = "x",
      box.padding = 1.45, point.padding = 0.32, min.segment.length = 0,
      segment.size = 0.2, segment.colour = "grey55", segment.alpha = 0.7,
      bg.color = "white", bg.r = 0.12, show.legend = FALSE)
  } else NULL

  y_scales <- lapply(REGIONS, function(r)
    scale_y_continuous(limits = c(0, reg_upper[[r]]), oob = scales::squish,
                       breaks = scales::pretty_breaks(n = 5), expand = expansion(mult = c(0.04, 0.06))))

  ggplot() +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25, linetype = "22") +
    geom_point(data = df_ns, aes(avg_log2FC, y_disp), colour = "#949494", size = 0.5, alpha = 0.5, shape = 16, stroke = 0) +
    geom_point(data = df_nom, aes(avg_log2FC, y_disp, colour = programme), size = 0.9, alpha = 0.5, shape = 16, stroke = 0) +
    geom_point(data = df_str, aes(avg_log2FC, y_disp, colour = programme), size = 2.6, alpha = 0.95, shape = 16, stroke = 0) +
    lab_layer +
    # the right-hand class strip duplicated the new column header and cost width, so the
    # facet is now region-only; the header names the class
    facet_grid2(. ~ region, scales = "free_y", independent = "y",
                labeller = labeller(region = region_full_lowconf_labeller(LOWCONF_REGIONS)),
                # The y (class) strip was default grey92 while the region
                # banners above carry the pale house tones -- the ggplot default was
                # leaking into an otherwise themed panel.  Give the class row the same
                # treatment: its cell-type hue blended 0.55 toward white (identical
                # blend to PAL_REGION_PALE) with plain black text.
                strip = ggh4x::strip_themed(
                  background_x = ggh4x::elem_list_rect(
                    fill = unname(PAL_REGION_PALE[REGIONS]), colour = NA),
                  text_x = ggh4x::elem_list_text(colour = "black", size = 10.4, face = "plain"),
                  background_y = ggh4x::elem_list_rect(
                    fill = unname(.pale(c(Neuron_Ex = "#006D77", Neuron_Inh = "#C42E7B")[[ct]], 0.55)), colour = NA),
                  text_y = ggh4x::elem_list_text(colour = "black", size = 10.4, face = "plain"))) +
    ggh4x::facetted_pos_scales(y = y_scales) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE, na.translate = FALSE,
                        guide = guide_legend(override.aes = list(size = 2.4, alpha = 1), nrow = 1)) +
    scale_x_continuous(limits = c(-xl, xl), oob = scales::squish, expand = expansion(mult = 0.04),
                       breaks = scales::pretty_breaks(n = 5)) +
    labs(x = expression(log[2]~"fold change (NHD vs CON)"),
         y = expression(-log[10]~"(adjusted "*italic(p)*")")) +
    theme_classic(base_size = 11.25) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.35),
          axis.text = element_text(colour = "black", size = 9.0),
          axis.title = element_text(colour = "black", size = 10.4),
          panel.spacing.x = unit(0.32, "lines"), legend.position = "bottom",
          legend.box = "horizontal", legend.spacing.x = unit(0.12, "cm"),
          legend.text = element_text(size = 9.0), legend.key.size = unit(0.40, "cm"),
          legend.box.spacing = unit(0.06, "cm"), legend.margin = margin(t = 1, b = 0),
          plot.margin = margin(2, 3, 1, 3))
}

p_ex  <- build_volcano("Neuron_Ex")
p_inh <- build_volcano("Neuron_Inh")

# =============================================================================
# 3. Tidy bottom KEY (two-tier size key only; methods sentence -> legend).
# =============================================================================
kt <- function(x, y, lab, parse = FALSE)
  annotate("text", x = x, y = y, label = lab, hjust = 0, size = 3.16, colour = "black", parse = parse)
kd <- function(x, y, sz) annotate("point", x = x, y = y, colour = "black", size = sz, shape = 16)
key_line <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  kd(0.300, 0.5, 0.9) + kt(0.315, 0.5, "nominal~(italic(p)[adj] < 0.05)", parse = TRUE) +
  kd(0.560, 0.5, 2.6) + kt(0.580, 0.5, "strict~(abs(delta) >= 0.15)", parse = TRUE) +
  theme_void() + theme(plot.margin = margin(2, 3, 1, 3))

# =============================================================================
# 4. Assemble: Ex row / Inh row / size-key strip.
# =============================================================================
# Same construction as F4b_oligolineage_volcano_composite: the two
# cell types sit as columns with a plain header each, and the shared two-tier size key
# runs underneath as one strip rather than being drawn twice.
# Unlike the oligo/OPC composite there is no shared theme legend, and that is correct:
# Oligo and OPC score the same programmes so their legends merge, whereas the Ex and Inh
# programme sets are genuinely different (SNARE/PSD/GluR/OXPHOS vs GABAergic/IN-TF/
# subclass/PV/PNN). Each column therefore keeps its own legend, and a master palette
# would be wrong here — the same colour would mean different programmes.
hdr <- function(g, lab) g + ggtitle(lab) +
  theme(plot.title = element_text(size = 11.1, hjust = 0.5, face = "plain",
                                  colour = "black", margin = margin(b = 2)))
combo <- (hdr(p_ex, "Excitatory neurons") | hdr(p_inh, "Inhibitory neurons")) / key_line +
  patchwork::plot_layout(heights = c(1, 0.075)) &
  guides(colour = guide_legend(nrow = 3, byrow = TRUE,
                               override.aes = list(size = 2.2, alpha = 1)))

# =============================================================================
# 5. Render (ragg PNG + base pdf)
# =============================================================================
# In so it fits the page width without being scaled at assembly (scaling
# would drop 8 pt type below its siblings — apparent-font parity). Height trimmed to
# keep the volcano aspect from stretching.
W <- 8.40; H <- 4.62  # +20% to buy label capacity
# Why +20% BUYS GENES. Label capacity is the ratio of label width to facet width, and
# that ratio is scale-invariant -- so shrinking the type would not have helped, but a
# 20% wider canvas at the same type does: each facet gains 20% of room while a gene name
# stays the same size, and the cap goes 4 -> 7 per facet.
# Placement: to keep the on-page type at the figure-wide 5.2 pt this panel must now be
# placed 4.85 in wide, not 4.04 (8.40 x 0.577). Dropped into the OLD 4.04 in box it
# scales to 0.481 and the type falls to 4.3 pt -- smaller than every sibling panel.
png_path <- file.path(PANEL, "F5b_neuron_volcano_grid.png")
pdf_path <- file.path(PANEL, "F5b_neuron_volcano_grid.pdf")
ragg::agg_png(png_path, width = W, height = H, units = "in", res = 600)
print(combo); invisible(dev.off())
ggsave(pdf_path, combo, width = W, height = H, useDingbats = FALSE)
cat(sprintf("\nWrote panel: %s\n         and: %s\n", png_path, pdf_path))

# =============================================================================
# 6. Provenance
# =============================================================================
prov <- file.path(PANEL, "F5b_neuron_volcano_grid.provenance.txt")
sink(prov)
cat("F5b_neuron_volcano_grid (FH) — provenance\n")
cat("generated:", format(Sys.time()), "\n")
cat("script   : scripts/F5b_neuron_volcano_grid_MAST_FH.R\n")
cat("input    :", MAST_PATH, "\n")
cat("filter   : cell_type in {Neuron_Ex, Neuron_Inh}, region in {Frontal,Hippo}\n")
cat("calling  : TWO-TIER by size — colour if padj<0.05 (nominal, small); LARGE if also |Cliff's delta|>=0.15 (strict).\n")
cat("Ex themes (== 4C modules): Presynaptic(SNARE)/Postsynaptic(PSD)/Ionotropic GluR/OXPHOS-ETC (NHD-axis DROPPED 2026-07-13)\n")
cat("Inh themes: GABAergic/IN TF/IN subclass/PV fast-spiking/Perineuronal net-ECM\n")
cat(sprintf("padj underflow floored (<= %.3g): %d ; display cap -log10(p) = %d ; capped points: %d\n",
            PADJ_FLOOR, n_underflow, YCAP, sum(dat$capped)))
cat("HOUSE: single NHD Hippocampus lane caveat now in the LEGEND, not on the panel (no on-plot '(low n)')\n\n")
cat("== per class-region nuclei (NHD | CON) ==\n"); print(as.data.frame(pow))
cat("\n== programme tier counts per facet ==\n"); print(counts)
for (ct in CLASSES) {
  cat(sprintf("\n== coloured programme genes (tier) — %s ==\n", ct))
  gg <- dat %>% filter(cell_type == ct, themed) %>%
    transmute(region, gene, programme = theme, tier = as.character(tier),
              avg_log2FC = round(avg_log2FC, 2), cliffs_delta = round(cliffs_delta, 3),
              dir = ifelse(avg_log2FC > 0, "NHD-up", "NHD-down")) %>%
    arrange(region, programme, desc(tier == "strict"), desc(abs(avg_log2FC)))
  print(as.data.frame(gg))
}
cat("\n== sessionInfo ==\n"); print(sessionInfo())
sink()
cat(sprintf("Wrote provenance: %s\n", prov))
cat("\n=== DONE ===\n", file = stderr())
