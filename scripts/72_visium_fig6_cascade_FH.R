#!/usr/bin/env Rscript
# =============================================================================
# 72_visium_fig6_cascade_FH.R — Figure 6, column 3: the cascade in tissue.
#   6_laminar_profile     (g)  Figure 5's deep-layer vulnerability, spatially
#   6_micro_vs_myelin     (h)  does myelin fail where the microglia are?
#   6_crossmodality       (i)  snRNA per-cell delta vs Visium spot delta
# -----------------------------------------------------------------------------
# (g) the axis is validated before it is used. There is no layer annotation on these
#     sections, so depth is derived (see _visium_depth_FH.R). Before any disease claim
#     rests on it, it has to order the layer markers correctly in the controls, and it
#     does — Spearman rho against depth, CON_Frontal1 / CON_Frontal2:
#         upper  CUX2 +0.25/+0.15   LAMP5 +0.22/+0.16   CALB1 +0.13/+0.05
#         deep   HS3ST4 -0.21/-0.07  SEMA3E -0.18/-0.08  TLE4 -0.11/-0.06
#                FOXP2 -0.11/-0.10   PCP4 -0.08/-0.19
#     Every upper marker rises toward the pia and every deep marker falls, in both
#     control sections. RORB is the one inconsistent marker (+0.05/-0.18), which is
#     expected: it is a mid-cortical L4/L5 marker with no monotone expectation.
# The GRADIENTS are SHALLOW (|rho| 0.05-0.25). This is a 55 um-resolution proxy
#     built from transcript, not a laminar segmentation. The panel therefore shows the
#     PROFILES, so the reader sees the size of the effect and the section-to-section
#     spread, and the legend must call this the weakest panel of the three.
#
# (h) is the adjacency the paper currently only implies: Figure 3 says microglia are
#     activated, Figure 4 says myelin fails; this asks whether the two coincide in the
#     same 55 um spot. Split by compartment, because myelin lives in white matter and
#     pooling the two would let a compartment shift masquerade as a local relationship.
#
# (i) is the one panel that says "two modalities, same answer". The programme
#     definitions are literally the same objects (_oligo_programmes_FH.R), so the two
#     axes are comparable by construction rather than by assertion.
#
# n = 2 sections per condition, one donor per condition. Lines are per section, never
# pooled into a single condition curve, and there are no p-values anywhere.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales); library(ggrepel)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_visium_depth_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(file.path(TDIR, "visium_dm_module_spots_FH.csv"),
              check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("layer genes missing — re-run 65 (v3)" = "gene_CUX2" %in% names(d))
cat("== deriving the cortical-depth axis ==\n")
d <- add_cortical_depth(d)
SEC <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")
d$sample_id <- factor(d$sample_id, levels = SEC)
d$Condition <- factor(d$Condition, levels = c("CON","NHD"))

# ---------------------------------------------------------------------------
# (g) laminar profiles
# ---------------------------------------------------------------------------
UPPER <- c("CUX2","LAMP5","CALB1"); DEEP <- c("PCP4","TLE4","FOXP2","SEMA3E","HS3ST4")
gm <- d %>% filter(!is_wm, !is.na(depth))
prof <- gm %>%
  select(sample_id, Condition, depth, all_of(paste0("gene_", c(UPPER, DEEP)))) %>%
  pivot_longer(-c(sample_id, Condition, depth), names_to = "gene", values_to = "e") %>%
  mutate(gene = sub("^gene_", "", gene),
         class = ifelse(gene %in% UPPER, "Upper-layer markers", "Deep-layer markers")) %>%
  # Scale each gene within each section to its own mean, so the profile is a shape and
  # a section with more counts of a gene cannot dominate the class average
  group_by(sample_id, gene) %>% mutate(e = e / mean(e[e > 0], na.rm = TRUE)) %>% ungroup() %>%
  mutate(bin = cut(depth, breaks = seq(0, 1, length.out = 11),
                   include.lowest = TRUE, labels = FALSE)) %>%
  filter(!is.na(bin)) %>%
  group_by(class, sample_id, Condition, bin) %>%
  summarise(m = mean(e, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(mid = (bin - 0.5)/10,
         class = factor(class, levels = c("Deep-layer markers", "Upper-layer markers")))

slope <- prof %>% group_by(class, sample_id, Condition) %>%
  summarise(slope = coef(lm(m ~ mid))[2], .groups = "drop")
cat("\n== laminar gradient slope per section (expression per unit depth) ==\n")
print(as.data.frame(slope %>% mutate(slope = round(slope, 3))), row.names = FALSE)

p_g <- ggplot(prof, aes(mid, m, colour = Condition, group = sample_id)) +
  geom_line(linewidth = 0.42) +
  geom_point(size = 0.7, stroke = 0) +
  scale_colour_manual(values = PAL_COND, name = NULL) +
  facet_wrap(~ class, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("WM edge", "mid", "pial")) +
  labs(x = "derived cortical depth", y = "expression\n(section mean = 1)") +
  theme_pub(base_size = 8) +
  theme(axis.text = element_text(size = 6.0, colour = "black"),
        axis.title = element_text(size = 6.6),
        strip.text = element_text(size = 6.4),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        axis.line = element_blank(), panel.spacing.x = unit(0.18, "cm"),
        legend.position = "bottom", legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 6.2), legend.margin = margin(t = -2),
        plot.margin = margin(4, 5, 3, 3))
BN <- "6_laminar_profile"; W <- 3.20; H <- 1.95
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_g, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_g, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))

# ---------------------------------------------------------------------------
# (h) myelin against microglial content, in the same spot
# ---------------------------------------------------------------------------
mm <- d %>%
  mutate(compartment = ifelse(is_wm, "White matter", "Cortex")) %>%
  group_by(sample_id) %>%
  mutate(mbin = ntile(micro_prop, 8)) %>% ungroup() %>%
  group_by(compartment, sample_id, Condition, mbin) %>%
  summarise(micro = mean(micro_prop), myelin = mean(`Structural myelin`),
            n = n(), .groups = "drop") %>%
  mutate(compartment = factor(compartment, levels = c("White matter", "Cortex")))
cat("\n== myelin score across microglial-content octiles ==\n")
print(as.data.frame(mm %>% filter(mbin %in% c(1, 8)) %>%
      mutate(across(c(micro, myelin), ~round(.x, 3)))), row.names = FALSE)

p_h <- ggplot(mm, aes(mbin, myelin, colour = Condition, group = sample_id)) +
  geom_line(linewidth = 0.42) + geom_point(size = 0.7, stroke = 0) +
  scale_colour_manual(values = PAL_COND, name = NULL) +
  facet_wrap(~ compartment, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 4, 8), labels = c("low", "mid", "high")) +
  labs(x = "microglial content (spot octile)", y = "structural myelin\nscore") +
  theme_pub(base_size = 8) +
  theme(axis.text = element_text(size = 6.0, colour = "black"),
        axis.title = element_text(size = 6.6),
        strip.text = element_text(size = 6.4),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        axis.line = element_blank(), panel.spacing.x = unit(0.18, "cm"),
        legend.position = "bottom", legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 6.2), legend.margin = margin(t = -2),
        plot.margin = margin(4, 5, 3, 3))
BN <- "6_micro_vs_myelin"; W <- 3.20; H <- 1.95
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_h, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_h, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))

# ---------------------------------------------------------------------------
# (i) cross-modality concordance
# ---------------------------------------------------------------------------
SN <- file.path(TDIR, "oligo_lineage_deltas_FH.csv")
VS <- file.path(TDIR, "fig6_supplychain_FH.csv")
if (file.exists(SN) && file.exists(VS)) {
  sn <- read.csv(SN, stringsAsFactors = FALSE) %>%
    filter(Region == "Frontal", lineage == "Oligodendrocyte") %>%
    transmute(arm = signature, sn_delta = cliff_d_depth_matched)
  vs <- read.csv(VS, stringsAsFactors = FALSE) %>%
    filter(part == "NHD - CON") %>% transmute(arm, vis_delta = value)
  cm <- inner_join(sn, vs, by = "arm")
  cat(sprintf("\ncross-modality: %d programmes shared\n", nrow(cm)))
  if (nrow(cm) >= 3) {
    rho <- suppressWarnings(cor(cm$sn_delta, cm$vis_delta, method = "spearman"))
    cat(sprintf("Spearman rho = %.2f over %d programmes\n", rho, nrow(cm)))
    print(as.data.frame(cm %>% mutate(across(where(is.numeric), ~round(.x, 3)))),
          row.names = FALSE)
    p_i <- ggplot(cm, aes(sn_delta, vis_delta)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
      geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
      geom_point(size = 1.9, shape = 21, fill = "#3E7CB1", colour = "black", stroke = 0.25) +
      ggrepel::geom_text_repel(aes(label = arm), size = 1.95, colour = "black", seed = 42,
                               box.padding = 0.42, min.segment.length = 0.2,
                               segment.size = 0.18, segment.colour = "grey60",
                               max.overlaps = 20) +
      labs(x = expression(atop("snRNA per-cell Cliff's " * delta,
                               "(oligodendrocyte, frontal)")),
           y = expression(atop("Visium spot " * Delta,
                               "(matched depth + content)"))) +
      theme_pub(base_size = 8) +
      theme(axis.text = element_text(size = 6.0, colour = "black"),
            axis.title = element_text(size = 6.6),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
            axis.line = element_blank(), plot.margin = margin(4, 6, 3, 3))
    BN <- "6_crossmodality"; W <- 2.70; H <- 2.15
    ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_i, width = W, height = H, useDingbats = FALSE)
    ggsave(file.path(PANEL, paste0(BN, ".png")), p_i, width = W, height = H, dpi = 600,
           device = ragg::agg_png)
    cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))
    write.csv(cm, file.path(TDIR, "fig6_crossmodality_FH.csv"), row.names = FALSE)
  } else cat("SKIPPED (i): fewer than 3 shared programmes\n")
} else cat("SKIPPED (i): missing", SN, "or", VS, "\n")

write.csv(prof,  file.path(TDIR, "fig6_laminar_profile_FH.csv"), row.names = FALSE)
write.csv(slope, file.path(TDIR, "fig6_laminar_slope_FH.csv"), row.names = FALSE)
write.csv(mm,    file.path(TDIR, "fig6_micro_vs_myelin_FH.csv"), row.names = FALSE)
cat("=== DONE ===\n")
