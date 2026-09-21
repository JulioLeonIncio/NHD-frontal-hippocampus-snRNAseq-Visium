#!/usr/bin/env Rscript
# =============================================================================
# 4E_LR_pairs_ExInh_FH.R — F5e2_LR_pairs_ExInh_FH.R — Figure 4 (Part B): microglia -> Ex/Inh neuron ligand-receptor pair bubble, frontal + hippocampus. Faithful port of
# manuscript_7fig/scripts/10h_v3_5E_LR_pairs_ExInh.R to the FH rebuild.
# -----------------------------------------------------------------------------
# y = curated micro->neuron LR pair (grouped by ligand family);
# x = Condition (CON/NHD); facet class(Ex/Inh) x Region (Frontal|Hippocampus).
# fill = signed delta prob (NHD - CON) steel<->rose; size = summed comm prob.
#
# FH data source: the withNeurons_v2 CellChat objects (major-class Neuron_Ex/Inh
# targets — no subtype roll-up needed; a direct Neuron_Ex->Ex / Neuron_Inh->Inh
# rename replaces the reference's Azimuth/Tippani collapse).
#
# Critical gotcha: pairs keyed by interaction_name_2 (the ligand-RECEPTOR key),
# not pathway_name (GAS6's pathway is "GAS", ENTPD1's is "CD39").
#
# Curated FH whitelist = the 5 significant micro->neuron axes:
#   SPP1-(ITGAV/ITGA4/8/9+ITGB1), GAS6-TYRO3, GRN-SORT1, SEMA4D-PLXNB2,
#   ENTPD1-ADORA1.  Only pairs present in >=1 FH context are shown; absent
#   whitelist pairs are dropped + logged (never a blank row).
# Valence SPLIT (house rule): SPP1 and GAS6 shown as separate ligand-family row
# blocks — not lumped. GAS6-TYRO3 is trophic-but-UNCOUPLED (receptor fate in the
# receptor panel). No causal verbs; the delta is the effect size.
#
# Full region names on the pale x-banner; no in-panel caption; single-donor
# Hippocampus flagged " (low n)" at the region strip (figure-wide).
#
# Output: figures/Figure_5/panels/F5e2_LR_pairs_ExInh.{pdf,png} + stats CSV.
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
  library(patchwork)
  library(CellChat); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales); library(forcats); library(ggh4x)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))

CC    <- file.path(PROJ, "data", "cellchat_per_region")
PANEL <- file.path(PROJ, "figures", "Figure_5", "panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

REGIONS <- c("Frontal","Hippo")
LOW_CONF_REGIONS <- "Hippo"      # single Hippocampus donor/lane (figure-wide caveat)
SENDER  <- "Micro-PVM"
NEURON  <- c("Neuron_Ex","Neuron_Inh")

# Curated micro->neuron LR whitelist (FH-audited; the 5 axes, receptor-resolved).
LR_WHITELIST <- c(
  "SPP1 - (ITGAV+ITGB1)", "SPP1 - (ITGA9+ITGB1)",
  "SPP1 - (ITGA8+ITGB1)", "SPP1 - (ITGA4+ITGB1)",
  "GAS6 - TYRO3",
  "GRN - SORT1",
  "SEMA4D - PLXNB2",
  "ENTPD1 - ADORA1")

load_region <- function(region) {
  path <- file.path(CC, sprintf("merged_NHDvsCON_%s_withNeurons_v2.rds", region))
  if (!file.exists(path))
    stop(sprintf("MISSING %s — CellChat withNeurons prep for %s failed", path, region))
  m  <- readRDS(path)
  df <- subsetCommunication(m, slot.name = "net")
  bind_rows(df$CON %>% mutate(Condition = "CON"),
            df$NHD %>% mutate(Condition = "NHD")) %>% mutate(Region = region)
}

raw <- bind_rows(lapply(REGIONS, load_region)) %>% as_tibble() %>%
  filter(source == SENDER, target %in% NEURON) %>%
  mutate(class = dplyr::recode(as.character(target),
                               Neuron_Ex = "Ex", Neuron_Inh = "Inh"))
stopifnot(nrow(raw) > 0)
cat(sprintf("Raw micro->neuron rows: %d\n", nrow(raw)))

# One aggregate per (Region, class, LR pair); FH objects have a single Ex/Inh
# receiver so prob_sum == prob, but sum keeps parity with the reference grammar.
agg <- raw %>%
  group_by(Region, class, interaction_name_2, pathway_name, Condition) %>%
  summarise(prob_sum = sum(prob, na.rm = TRUE),
            min_pval = suppressWarnings(min(pval, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(min_pval = ifelse(is.infinite(min_pval), 1, min_pval))

wide <- agg %>%
  pivot_wider(id_cols = c(Region, class, interaction_name_2, pathway_name),
              names_from = Condition, values_from = c(prob_sum, min_pval),
              values_fill = list(prob_sum = 0, min_pval = 1)) %>%
  mutate(min_p_either = pmin(min_pval_CON, min_pval_NHD),
         q_either     = p.adjust(min_p_either, method = "BH"),
         delta        = prob_sum_NHD - prob_sum_CON,
         peak         = pmax(prob_sum_CON, prob_sum_NHD))
# Note: q_either is kept for the audit CSV/log (transparency) but is no
# longer encoded as a "thick border" in the panel. On the curated FH whitelist that
# border was NON-DISCRIMINATING (essentially every curated pair passes BH q<0.1), so
# it read as a spurious significance claim. Dropped from the plot (Part B fix).

# Diagnostic: all present micro->neuron pairs (for audit/curation transparency).
wide %>% group_by(interaction_name_2, pathway_name) %>%
  summarise(net_delta = sum(delta), max_abs = max(abs(delta)),
            min_q = min(q_either), n_ctx = dplyr::n(), .groups = "drop") %>%
  arrange(desc(max_abs)) %>%
  write.csv(file.path(TDIR, "fig4E_present_micro_neuron_LR.csv"), row.names = FALSE)

present_lr <- sort(unique(wide$interaction_name_2))
absent_lr  <- setdiff(LR_WHITELIST, present_lr)
if (length(absent_lr) > 0)
  cat("WARNING whitelist pair(s) absent from FH objects, dropping:",
      paste(absent_lr, collapse = " | "), "\n")
LR_WHITELIST <- intersect(LR_WHITELIST, present_lr)
stopifnot(length(LR_WHITELIST) >= 1)
cat("Curated LR pairs used (", length(LR_WHITELIST), "):",
    paste(LR_WHITELIST, collapse = " | "), "\n")

# per-axis significance audit (the requested SPP1/GAS6/GRN/SEMA4D/ENTPD1 status)
cat("\n== micro->neuron per-axis status (prob CON/NHD, delta, q<0.1) ==\n")
audit <- wide %>% filter(interaction_name_2 %in% LR_WHITELIST) %>%
  transmute(Region, class, interaction_name_2, pathway_name,
            prob_CON = round(prob_sum_CON, 4), prob_NHD = round(prob_sum_NHD, 4),
            delta = round(delta, 4), q = round(q_either, 3),
            q_lt_0.1 = q_either < 0.1) %>%   # audit-only; not a panel encoding
  arrange(match(pathway_name, c("SPP1","GAS","GRN","SEMA4","CD39")),
          interaction_name_2, Region, class)
print(as.data.frame(audit), row.names = FALSE)

bubble <- wide %>%
  filter(interaction_name_2 %in% LR_WHITELIST) %>%
  pivot_longer(cols = c(prob_sum_CON, prob_sum_NHD),
               names_to = "Condition", values_to = "prob") %>%
  mutate(Condition = sub("prob_sum_", "", Condition),
         class     = factor(class,     levels = c("Ex","Inh")),
         Region    = factor(Region,    levels = REGIONS),
         Condition = factor(Condition, levels = c("CON","NHD")))

# Order LR pairs: keep ligand-family blocks together (SPP1, GAS6, GRN, SEMA4D,
# ENTPD1), within a family by net delta. rev() so the block reads top->bottom.
fam_order <- c("SPP1","GAS","GRN","SEMA4","CD39")
y_rank <- bubble %>% distinct(interaction_name_2, pathway_name, delta) %>%
  group_by(interaction_name_2, pathway_name) %>%
  summarise(net_delta = sum(delta, na.rm = TRUE), .groups = "drop") %>%
  mutate(pathway_name = factor(pathway_name, levels = fam_order)) %>%
  arrange(pathway_name, net_delta)
bubble$interaction_name_2 <- factor(bubble$interaction_name_2,
                                    levels = rev(unique(y_rank$interaction_name_2)))

# Signed-delta colour cap: data-driven 95th pct of |delta| (stated). SPP1 vs the
# tiny GRN deltas differ ~10x; winsorize colour only, raw delta kept in the CSV.
delta_cells <- bubble %>% distinct(Region, class, interaction_name_2, delta)
delta_lim   <- as.numeric(stats::quantile(abs(delta_cells$delta), 0.95, na.rm = TRUE))
if (!is.finite(delta_lim) || delta_lim == 0)
  delta_lim <- max(abs(delta_cells$delta), na.rm = TRUE)
n_wins <- sum(abs(delta_cells$delta) > delta_lim, na.rm = TRUE)
cat(sprintf("\nbubble: %d dots / %d pairs; delta colour cap |%.4f| (95th pct), %d winsorized for colour\n",
            nrow(bubble), length(unique(bubble$interaction_name_2)), delta_lim, n_wins))

write.csv(bubble, file.path(TDIR, "fig4E_micro_neuron_LR_bubble.csv"), row.names = FALSE)

# data-driven size breaks (positive comm prob), always in range
pos_prob  <- bubble$prob[bubble$prob > 0]
sz_breaks <- sort(unique(signif(stats::quantile(pos_prob, c(0.2, 0.6, 0.9), na.rm = TRUE), 1)))
sz_labels <- format(sz_breaks, trim = TRUE, scientific = FALSE)

present_regions <- REGIONS[REGIONS %in% as.character(bubble$Region)]

p <- ggplot(bubble,
            aes(x = Condition, y = interaction_name_2, size = prob,
                fill = pmax(pmin(delta, delta_lim), -delta_lim))) +
  geom_point(shape = 21, stroke = 0.25, color = "grey20", alpha = 0.95) +
  # (non-discriminating on the curated
  #  whitelist, read as a spurious significance claim. delta is the effect size.)
  facet_grid2(class ~ Region, switch = "y", scales = "free_y", space = "free_y",
              labeller = labeller(
                Region = ggplot2::as_labeller(c(Frontal = "Frontal", Hippo = "Hippo"))),
              strip = strip_region_x(present_regions, fontsize = 10.1,
                                     low_conf = LOW_CONF_REGIONS,
                                     clip = "off")) +   # let full "Hippocampus (low n)" overflow the narrow banner
  scale_fill_gradient2(low = unname(PAL_DGE["Up in CON"]), mid = "white",
                       high = unname(PAL_DGE["Up in NHD"]), midpoint = 0,
                       limits = c(-delta_lim, delta_lim), oob = scales::squish,
                       name = expression(Delta ~ "prob"),
                       breaks = c(-delta_lim, 0, delta_lim),
                       labels = function(b) ifelse(abs(b) < 1e-12, "0", sprintf("%+.3f", b)),
                       guide = guide_colorbar(order = 1, title.position = "left", title.vjust = 0.9,
                                              barheight = unit(0.26, "cm"),
                                              barwidth = unit(1.9, "cm"))) +
  scale_size_continuous(range = c(0.8, 3.6), name = "prob",
                        breaks = sz_breaks[c(1, length(sz_breaks))],
                        labels = sz_labels[c(1, length(sz_labels))],
                        guide = guide_legend(order = 2, title.position = "left", nrow = 1)) +
  # Sender declared on the y axis (propagated across every l-r panel).
  # The rows are ligand - receptor pairs, so the axis names both halves; the grey strips
  # beside them are CellChat pathway names, which are not always the ligand (TENM4 sits
  # under "ADGRL", named for the receptor; VSIR under "vista"; GAS6 under "GAS"). Naming
  # the axis for the pair — and the x axis for the receiver — lets the reader get sender
  # and receiver from the panel instead of the legend.
  # Never a unicode arrow here: base pdf() cannot encode U+2192 and it has already broken
  # a render in this project. ASCII hyphen only.
  labs(x = NULL, y = "Microglial ligand - receptor") +
  theme_pub(base_size = 12.6) +
  theme(strip.background.y = element_rect(fill = "grey92", colour = NA),
        strip.text.y.left  = element_text(colour = "black", face = "plain", size = 10.1),
        axis.text.x        = element_text(size = 8.8),
        axis.text.y        = element_text(size = 8.8, face = "italic"),
        axis.title         = element_text(size = 10.1),
        panel.border       = element_rect(colour = "grey70", fill = NA, linewidth = 0.25),
        panel.spacing.x    = unit(0.20, "lines"),
        panel.spacing.y    = unit(0.20, "lines"),
        panel.grid.major.y = element_line(linewidth = 0.10, colour = "grey94"),
        panel.grid.major.x = element_blank(),
        # Legend under the panel: at the parity type size a right-hand
        # legend left the four data columns ~0.35 in apart and CON/NHD overprinted.
        legend.position    = "bottom", legend.box = "vertical",
        legend.box.just    = "center", legend.justification = "center",
        legend.box.spacing = unit(0.10, "cm"), legend.margin = margin(t = 0, b = 0),
        legend.key.size    = unit(0.22, "cm"),
        legend.title       = element_text(size = 9.1),
        legend.text        = element_text(size = 8.8),
        plot.margin        = margin(4, 4, 4, 4),
        strip.placement    = "outside")

W <- 3.9; H <- 3.75   # +0.15 in for the sender header
# Sender marker — width-aware so the arrow neither collides with the label nor
# drifts away from it (22_publication_theme_FH.R).
p <- with_sender(p, "Microglia", size = 3.56, panel_width = W, panel_height = H, arrow_x = 1.760)
BN <- "F5e2_LR_pairs_ExInh"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nWrote %s.{pdf,png} + fig4E_micro_neuron_LR_bubble.csv\n", BN))
cat("=== DONE ===\n", file = stderr())
