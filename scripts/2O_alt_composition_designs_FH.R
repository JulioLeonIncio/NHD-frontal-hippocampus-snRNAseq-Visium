#!/usr/bin/env Rscript
# =============================================================================
# 2O_alt_composition_designs_FH.R — Alternative designs for the NHD immune proportion panel ("generate a couple of alternatives more").
# -----------------------------------------------------------------------------
# Same numbers as 2O/2O2, three different ways of showing them. All read the
# composition CSV, so no re-analysis and no atlas load — the alternatives cannot
# disagree with the shipped panel.
#
#  ALT 1  slope / dumbbell on a log10 axis — all three populations and both regions
#         in one small panel. The log axis is what makes it possible to show
#         microglia (~1000 per 10k) beside T cells (~5 per 10k) without flattening
#         either; the slope carries the message (microglia down, PVM up, T flat).
#  ALT 2  log2 fold-change of abundance vs CON — the most compact statement of the
#         same thing; one bar per population per region, zero line = no change.
#  ALT 3  paired stacked composition, CON|NHD adjacent per region, with the immune
#         compartment expressed as % — the conventional view, kept for comparison.
#
# No stars anywhere: n = 1 donor per condition, the shift is the effect size.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggh4x)
})
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
PANEL <- file.path(PROJ, "figures", "Figure_2", "panels", "_alternatives")
dir.create(PANEL, recursive = TRUE, showWarnings = FALSE)
TBL <- file.path(PROJ, "tables", "micro_states")
comp <- read.csv(file.path(TBL, "immune_composition_FH.csv"))

IMM_LEVELS <- c("Microglia", "PVM", "T-lymphocyte")
PAL_IMM <- c("Microglia" = unname(PAL_CELLTYPE[["Micro-PVM"]]),
             "PVM" = "#6D4C41", "T-lymphocyte" = "#C77B0A")
comp <- comp %>%
  mutate(immune_class = factor(immune_class, levels = IMM_LEVELS),
         Condition = factor(Condition, levels = c("CON","NHD")),
         Region = factor(Region, levels = REGION_ORDER))

# apparent-type parity. House rule:
# on-page pt = declared pt x placement scale (placed_width/saved_width in the
# assembled figure); main-figure band 4.8-5.3 pt, target 5.1.
# Assembled placements in SuppFig3 (Fig3_suppl.pdf):
#   ALT1_immune_slope_log10      3.00 -> 2.273 in (0.758): 7 pt -> 5.31 pt on page
#   ALT2_immune_log2FC_abundance 2.90 -> 2.302 in (0.794): 7 pt -> 5.56 pt on page
# ALT3 is not placed in any figure, so it keeps base_size 7 untouched.
BASE_ALT1 <- 6.7   # 6.7 x 0.758 = 5.08 pt on page
BASE_ALT2 <- 7.2

# ---- ALT 1: slope on log10 -------------------------------------------------
p1 <- ggplot(comp, aes(x = Condition, y = per_10k_total,
                       colour = immune_class, group = immune_class)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.5) +
  facet_wrap2(~ Region, nrow = 1, labeller = region_full_lowconf_labeller(),
              strip = strip_region_x(REGION_ORDER, clip = "off")) +
  # Display the CD163+/F13A1+ marker name, not bare "PVM". SuppFig4 panel c
  # (2N) already uses the IMM_LABEL recode, so d/e/f saying "PVM" made one figure carry two
  # names for the same population. Bare "PVM" also asserts a perivascular origin we cannot
  # establish (CD163 can mark a microglial state) -- only the joint-class label "Micro-PVM"
  # is sanctioned. See MANIFEST_Fig2.md -> "Pooled wording".
  scale_colour_manual(values = PAL_IMM, labels = IMM_LABEL, name = NULL) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300, 1000),
                labels = c("3","10","30","100","300","1000")) +
  labs(x = NULL, y = "nuclei per 10,000 total (log scale)") +
  theme_pub(base_size = BASE_ALT1) +
  theme(legend.position = "bottom", legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(t = -2), panel.spacing = unit(6, "pt"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.2))
ggsave(file.path(PANEL, "ALT1_immune_slope_log10.png"), p1, width = 3.0, height = 2.2,
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "ALT1_immune_slope_log10.pdf"), p1, width = 3.0, height = 2.2, bg = "white")

# ---- ALT 2: log2 fold-change of abundance ---------------------------------
fc <- comp %>% select(Region, Condition, immune_class, per_10k_total) %>%
  pivot_wider(names_from = Condition, values_from = per_10k_total) %>%
  mutate(log2fc = log2(NHD / CON))
p2 <- ggplot(fc, aes(x = immune_class, y = log2fc, fill = immune_class)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_col(width = 0.6, linewidth = 0) +
  facet_wrap2(~ Region, nrow = 1, labeller = region_full_lowconf_labeller(),
              strip = strip_region_x(REGION_ORDER, clip = "off")) +
  # Display the CD163+/F13A1+ marker name, not bare "PVM". SuppFig4 panel c
  # (2N) already uses the IMM_LABEL recode, so d/e/f saying "PVM" made one figure carry two
  # names for the same population. Bare "PVM" also asserts a perivascular origin we cannot
  # establish (CD163 can mark a microglial state) -- only the joint-class label "Micro-PVM"
  # is sanctioned. See MANIFEST_Fig2.md -> "Pooled wording".
  # The three angled class names under each facet cost ~0.45 in
  # of height; the classes are keyed by fill in a one-row legend below instead (same PAL_IMM / IMM_LABEL
  # as panels d and e of the same figure), and the canvas drops 2.9 x 2.1 -> 2.3 x 1.55 in.
  scale_fill_manual(values = PAL_IMM, labels = IMM_LABEL, name = NULL,
                    guide = guide_legend(nrow = 1, override.aes = list(width = 0.5))) +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(breaks = c(-1, 0, 1, 2), expand = expansion(mult = c(0.08, 0.06))) +   # a -1 tick under the bars that reach -0.8
  labs(x = NULL, y = expression(log[2]~"abundance NHD / CON")) +
  theme_pub(base_size = BASE_ALT2) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "bottom", legend.key.size = unit(0.20, "cm"), legend.text = element_text(size = BASE_ALT2 * 0.85),
        legend.margin = margin(t = 1), legend.box.spacing = unit(4, "pt"), legend.spacing.x = unit(0.08, "cm"),
        panel.spacing = unit(6, "pt"))
ggsave(file.path(PANEL, "ALT2_immune_log2FC_abundance.png"), p2, width = 2.3, height = 1.65,
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "ALT2_immune_log2FC_abundance.pdf"), p2, width = 2.3, height = 1.65, bg = "white")

# ---- ALT 3: paired stacked composition ------------------------------------
p3 <- ggplot(comp, aes(x = Condition, y = pct_immune, fill = immune_class)) +
  geom_col(width = 0.66, linewidth = 0) +
  geom_text(data = comp %>% filter(immune_class != "Microglia", pct_immune > 1),
            aes(label = sprintf("%.1f", pct_immune)),
            position = position_stack(vjust = 0.5), size = 1.7, colour = "white") +
  facet_wrap2(~ Region, nrow = 1, labeller = region_full_lowconf_labeller(),
              strip = strip_region_x(REGION_ORDER, clip = "off")) +
  scale_fill_manual(values = PAL_IMM, labels = IMM_LABEL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "% of immune nuclei") +
  theme_pub(base_size = 7) +
  theme(legend.position = "bottom", legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(t = -2), panel.spacing = unit(6, "pt"))
ggsave(file.path(PANEL, "ALT3_immune_stacked_labelled.png"), p3, width = 2.8, height = 2.2,
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "ALT3_immune_stacked_labelled.pdf"), p3, width = 2.8, height = 2.2, bg = "white")

cat("=== log2 abundance fold-change (NHD/CON) ===\n"); print(as.data.frame(fc), digits = 3)
cat("\nwrote ALT1/ALT2/ALT3 to figures/Figure_2/panels/_alternatives/\n")
