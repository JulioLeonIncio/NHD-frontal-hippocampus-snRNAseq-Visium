#!/usr/bin/env Rscript
# =============================================================================
# 66_visium_markers_panel_FH.R — Figure 6: "6_depth_matched_markers".
# Selective loss in tissue, measured as counts per matched sequencing depth.
# -----------------------------------------------------------------------------
# Replaces `panel_f_prevalence`, which is not salvageable: it plots the fraction
# of spots in which a marker is detected, and NHD spots carry half the depth, so
# the panel largely measures library size. Depth-matched, CUX2's apparent loss is
# entirely artifact (raw 0.699 -> 1.046 matched) while PCP4 is entirely real
# (0.279 -> 0.278). Both were in the same panel.
#
# What this panel earns: "the loss is selective." Deep-layer neuronal markers and
# myelin fall while upper-layer neuronal markers do not, and glial activation genes
# rise — all measured in the same spots at the same depth. The upper-layer group is
# the internal positive control that licenses the word "selective"; without it the
# honest reading is only "global degradation".
#
# n = 2 sections per condition, one donor per condition. Per-section dots, ratio of
# section means, no p-values.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales); library(forcats)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)

MARK <- file.path(TDIR, "visium_dm_marker_counts_FH.csv")
stopifnot("MISSING marker table — run 65_visium_prep_depthmatched_FH.R" = file.exists(MARK))
m <- read.csv(MARK, stringsAsFactors = FALSE)

GRP <- c("Deep-layer neuron","Upper-layer neuron","Pan-neuronal","Myelin",
         "Microglial activation","Astrocyte reactive")
# A gene has to be measurable before a ratio means anything: require a mean of at
# least 0.02 counts per 3,000 UMI in the control sections (~1 count per 50 spots).
# The floor is applied to the maximum of the two arms, so a gene is kept if either
# condition can measure it: a control-only floor would delete exactly the genes that
# are absent in control and induced in disease (e.g. SERPINA3, CON 0.013 -> NHD 0.151,
# log2 +3.54, both NHD sections concordant; AIF1 +1.68), i.e. the glial-activation
# genes this panel exists to show.
det <- m %>% group_by(gene, Condition) %>%
  summarise(mu = mean(mean_count), .groups = "drop") %>%
  group_by(gene) %>% summarise(peak = max(mu), .groups = "drop")
drop <- det$gene[det$peak < 0.02]
if (length(drop)) cat("below the measurable floor in CON (dropped):",
                      paste(drop, collapse = ", "), "\n")
m <- m %>% filter(!gene %in% drop)

sec <- m %>% group_by(gene, group, Condition) %>%
  summarise(v = mean(mean_count), .groups = "drop")
rat <- sec %>% pivot_wider(names_from = Condition, values_from = v) %>%
  mutate(log2_ratio = log2(NHD / CON))
# per-section points: each NHD section against the CON section mean, so the spread
# shown is the between-section spread that actually exists (2 v 2)
pts <- m %>% group_by(gene, group, Condition, sample_id) %>%
  summarise(v = mean(mean_count), .groups = "drop") %>%
  left_join(sec %>% filter(Condition == "CON") %>% select(gene, con = v), by = "gene") %>%
  mutate(log2_ratio = log2(v / con))

ord <- rat %>% arrange(match(group, GRP), log2_ratio) %>% pull(gene)
lv <- function(d) d %>% mutate(gene = factor(gene, levels = ord),
                               group = factor(group, levels = GRP))
rat <- lv(rat); pts <- lv(pts)

cat("\n== log2(NHD/CON) counts per 3,000 UMI ==\n")
print(as.data.frame(rat %>% arrange(group, log2_ratio) %>%
      mutate(across(c(CON, NHD, log2_ratio), ~round(.x, 3)))))

p <- ggplot(rat, aes(x = log2_ratio, y = gene)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2_ratio, yend = gene), colour = "grey75",
               linewidth = 0.7, lineend = "round") +
  # Both conditions' section values: plotting only the NHD sections against
  # the CON mean hid the control's own between-section spread, so a reader could not tell
  # whether an effect cleared the noise. CUX2's two CON sections differ 1.9-fold.
  geom_point(data = pts, aes(x = log2_ratio, y = gene, fill = Condition),
             shape = 21, size = 1.4, colour = "grey30", stroke = 0.2, alpha = 0.9) +
  scale_fill_manual(values = PAL_COND, name = NULL) +
  # The summary diamond is now neutral: its old five-colour fill used hues nearly
  # identical to PAL_COND and read as a condition encoding.
  geom_point(shape = 23, size = 2.3, fill = "grey25", colour = "grey20", stroke = 0.3) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  labs(x = "log2 (NHD / CON), counts per 3,000 UMI", y = NULL) +
  theme_pub(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.3, colour = "black", face = "italic"),
        axis.text.x = element_text(size = 6.3),
        axis.title.x = element_text(size = 7),
        strip.text.y.left = element_text(size = 6.2, angle = 0, hjust = 1, colour = "black"),
        strip.background.y = element_rect(fill = "grey96", colour = "grey85", linewidth = 0.2),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        axis.line = element_blank(), panel.spacing.y = unit(0.10, "cm"))

BN <- "6_depth_matched_markers"
W <- 3.5; H <- 3.6
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nSaved %s.{png,pdf} (%.2f x %.2f in)\n", BN, W, H))
write.csv(rat %>% left_join(
  pts %>% filter(Condition == "NHD") %>% select(gene, sample_id, section_log2 = log2_ratio),
  by = "gene"), file.path(TDIR, "fig6_depth_matched_markers_FH.csv"), row.names = FALSE)
cat("Wrote tables/fig6_depth_matched_markers_FH.csv\n=== DONE ===\n")
