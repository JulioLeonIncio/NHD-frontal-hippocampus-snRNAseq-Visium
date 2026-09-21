#!/usr/bin/env Rscript
# =============================================================================
# 71_visium_fig6_overview_FH.R — Figure 6, column 1: the tissue and its controls.
#   6_tissue_depth_map   (a)  the four sections + the derived cortical-depth axis
#   6_depth_control      (b)  why every downstream panel is depth-matched
#   F6f_c2l_composition    (c)  deconvolved composition, CON vs NHD, per section
# -----------------------------------------------------------------------------
# Figure 6 exists to close the single-nucleus claims in tissue, so this column has to
# earn the reader's trust before the biology columns spend it:
#   (a) shows there is a laminar axis and that it was derived, not assumed. There is no
#       layer annotation on these sections, so white matter is each section's own
#       high-oligodendrocyte tail and depth is distance to it, scaled per section.
#   (b) is the panel that licenses everything else. NHD spots carry half the depth of
#       CON spots (per-section median UMI is written to tables/fig6_depth_control_FH.csv
#       by this script). Un-matched, every programme "falls" in NHD and the figure would
#       be measuring library size.
#   (c) is the payoff of moving to tissue: composition measured in intact tissue, where
#       the dissociation bias that shapes a nuclei prep does not operate.
#
# (c) design. the DECONVOLUTION is MIS-CALIBRATED and the panel must not
# hide IT. In control frontal cortex cell2location calls Endothelial 27.8% + Pericytes
# 22.7% = 50.5% vascular, against a plausible ~1-3% each, and oligodendrocytes at 5.1%
# where 20-40% is expected. The proportions are closed (they sum to exactly 1.000 per
# spot), so an over-called vascular compartment compresses every other type toward zero.
# Three consequences, all handled here:
#   1. absolute abundance is now primary. Proportions are compositionally closed against
#      the 52% oligodendrocyte loss; the absolute q05 estimate is not. Read from
#      c2l_MAIN/abundance_q05.csv rather than the proportion table.
#   2. A DECONVOLUTION-FREE estimator is plotted beside IT. Marker transcript per matched
#      depth needs no reference and no closure. If the two estimators agree, the call
#      does not rest on the deconvolution at all.
#   3. the reference is our own ATLAS (NHD_QC_harmony.rds), so the snRNA composition
#      result and the cell2location result are one line of evidence, not two. Only the
#      marker series is independent. The legend must say this.
#
# n = 2 sections per condition, one donor per condition. Per-section points, ratio of
# section means, no p-values: spots are sub-sampling units of one donor and condition is
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
source(file.path(PROJ, "scripts", "_visium_depth_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)

MOD <- file.path(TDIR, "visium_dm_module_spots_FH.csv")
QC  <- file.path(TDIR, "visium_dm_qc_FH.csv")
stopifnot("MISSING per-spot table — run 65 (v3)" = file.exists(MOD),
          "MISSING QC table — run 65 (v3)"      = file.exists(QC))
d  <- read.csv(MOD, check.names = FALSE, stringsAsFactors = FALSE)
qc <- read.csv(QC, stringsAsFactors = FALSE)
stopifnot("c2l columns missing — re-run 65 (v3)" = any(grepl("^c2l_", names(d))),
          "myelin score missing — re-run 65" = "Structural myelin" %in% names(d))
# Depth is derived here, from the shared helper, overwriting whatever the prep wrote:
# the prep's first version keyed white matter off the deconvolved oligodendrocyte tail,
# which is not spatially coherent on these sections (see _visium_depth_FH.R).
cat("\n== deriving the cortical-depth axis ==\n")
d <- add_cortical_depth(d)

SEC <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")
d$sample_id <- factor(d$sample_id, levels = SEC)
d$Condition <- factor(d$Condition, levels = c("CON","NHD"))
cat(sprintf("spots %d across %d sections\n", nrow(d), nlevels(d$sample_id)))

# per-section coordinates centred and scaled by the larger span (house spmap grammar)
d <- d %>% group_by(sample_id) %>%
  mutate(.sp = max(diff(range(x)), diff(range(y))),
         xr = (x - mean(range(x)))/.sp, yr = (y - mean(range(y)))/.sp) %>% ungroup()

# ---------------------------------------------------------------------------
# (a) the sections and the derived depth axis
# ---------------------------------------------------------------------------
# WM spots are drawn in their own colour rather than as depth 0, so the reader can see
# what the axis was measured from and judge whether the tail really is white matter.
# If those spots do not form a contiguous compartment the panel says so immediately --
# which is exactly how the first (cell2location-based) definition was caught.
pa_df <- d %>% mutate(shown = ifelse(is_wm, NA_real_, depth))
p_a <- ggplot(pa_df, aes(xr, yr)) +
  geom_point(aes(colour = shown), size = 0.26, stroke = 0) +
  scale_colour_viridis_c(option = "mako", direction = -1, na.value = "#C2185B",
                         limits = c(0, 1), breaks = c(0, 0.5, 1),
                         labels = c("WM edge", "mid", "pial"),
                         name = "cortical\ndepth") +
  facet_wrap(~ sample_id, nrow = 1) +
  coord_equal() +
  theme_pub(base_size = 8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.line = element_blank(),
        strip.text = element_text(size = 6.2),
        panel.spacing = unit(0.06, "cm"),
        legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.60, "cm"),
        legend.text = element_text(size = 6.0), legend.title = element_text(size = 6.2,
                                                                            lineheight = 0.9))
BN <- "6_tissue_depth_map"; W <- 5.20; H <- 1.70
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_a, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_a, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))

# ---------------------------------------------------------------------------
# (b) the depth control
# ---------------------------------------------------------------------------
med <- d %>% group_by(sample_id, Condition) %>%
  summarise(med = median(total_umi_pre), n = n(), .groups = "drop")
cat("\n== per-section median UMI before thinning ==\n"); print(as.data.frame(med))

p_b <- ggplot(d, aes(x = total_umi_pre, y = fct_rev(sample_id), fill = Condition)) +
  geom_violin(scale = "width", width = 0.85, linewidth = 0.18, colour = "grey30") +
  geom_point(data = med, aes(x = med, y = fct_rev(sample_id)), inherit.aes = FALSE,
             shape = 23, size = 1.5, fill = "white", colour = "black", stroke = 0.3) +
  geom_vline(xintercept = 3000, linetype = "dashed", colour = "black", linewidth = 0.35) +
  scale_fill_manual(values = PAL_COND, name = NULL) +
  scale_x_continuous(trans = "log10", labels = label_number(big.mark = ","),
                     breaks = c(1000, 3000, 10000, 30000)) +
  labs(x = "UMI per spot before thinning", y = NULL) +
  theme_pub(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.2, colour = "black"),
        axis.text.x = element_text(size = 6.2), axis.title.x = element_text(size = 6.6),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        axis.line = element_blank(),
        legend.position = "bottom", legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 6.2), legend.margin = margin(t = -2),
        plot.margin = margin(4, 5, 3, 3))
BN <- "6_depth_control"; W <- 2.60; H <- 1.70
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_b, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_b, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))

# ---------------------------------------------------------------------------
# (c) deconvolved composition
# ---------------------------------------------------------------------------
CT_LAB <- c("c2l_Micro-PVM" = "Microglia / PVM", "c2l_Oligo" = "Oligodendrocyte",
            "c2l_OPC" = "OPC", "c2l_Astro" = "Astrocyte",
            "c2l_Neuron_Ex" = "Excitatory neuron", "c2l_Neuron_Inh" = "Inhibitory neuron",
            "c2l_Endo" = "Endothelial", "c2l_Pericytes" = "Pericyte")
have <- intersect(names(CT_LAB), names(d))
# ---- absolute abundance (primary) ------------------------------------------
ABU <- file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN", "abundance_q05.csv")
stopifnot("MISSING absolute abundance table" = file.exists(ABU))
ab <- read.csv(ABU, check.names = FALSE, stringsAsFactors = FALSE)
ab_ct <- intersect(sub("^c2l_", "", names(CT_LAB)), names(ab))
ab <- ab %>% filter(spot_id %in% d$spot_id) %>%
  select(spot_id, all_of(ab_ct)) %>%
  left_join(d %>% select(spot_id, sample_id, Condition), by = "spot_id")
cat(sprintf("absolute abundance joined for %d of %d QC-passed spots\n", nrow(ab), nrow(d)))
comp <- ab %>% pivot_longer(all_of(ab_ct), names_to = "ct0", values_to = "p") %>%
  mutate(ct = paste0("c2l_", ct0)) %>%
  group_by(ct, Condition, sample_id) %>% summarise(v = mean(p), .groups = "drop")
sec_mean <- comp %>% group_by(ct, Condition) %>%
  summarise(m = mean(v), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = m) %>%
  mutate(log2_ratio = log2(NHD / CON), label = unname(CT_LAB[ct]))
pts <- comp %>% left_join(sec_mean %>% select(ct, CON), by = "ct") %>%
  mutate(log2_ratio = log2(v / CON), label = unname(CT_LAB[ct]))
ordl <- sec_mean %>% arrange(log2_ratio) %>% pull(label)
sec_mean$label <- factor(sec_mean$label, levels = ordl)
pts$label      <- factor(pts$label,      levels = ordl)

cat("\n== deconvolved composition, log2(NHD/CON) of section means ==\n")
print(as.data.frame(sec_mean %>% arrange(desc(log2_ratio)) %>%
      mutate(across(c(CON, NHD, log2_ratio), ~round(.x, 3))) %>%
      select(label, CON, NHD, log2_ratio)), row.names = FALSE)

# ---- the deconvolution-free estimator (marker transcript, matched depth) -----
MK <- list("c2l_Micro-PVM" = c("AIF1","C1QB","CD74"),
           "c2l_Oligo"     = c("MBP","PLP1","MOBP"),
           # (no astrocyte row: GFAP measures reactivity, not number)
           "c2l_Neuron_Ex" = c("SLC17A7"),
           "c2l_Neuron_Inh"= c("GAD1","GAD2"))
mk <- lapply(names(MK), function(k) {
  g <- paste0("gene_", intersect(MK[[k]], sub("^gene_", "", grep("^gene_", names(d), value = TRUE))))
  if (!length(g)) return(NULL)
  d %>% select(Condition, all_of(g)) %>%
    pivot_longer(all_of(g), names_to = "gene", values_to = "e") %>%
    group_by(gene, Condition) %>% summarise(m = mean(e), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = m) %>%
    summarise(ct = k, marker_log2 = mean(log2(NHD / CON)), n_gene = n())
}) %>% bind_rows()
mk$label <- unname(CT_LAB[mk$ct])
cat("\n== deconvolution-free marker estimator, log2(NHD/CON) ==\n")
print(as.data.frame(mk %>% mutate(marker_log2 = round(marker_log2, 3))), row.names = FALSE)
mk$label <- factor(mk$label, levels = ordl)

# The two estimators go in one legend, not an axis subtitle: explanatory text under an
# axis is caption text inside a panel, which this project does not do. Condition moves
# from fill to colour so that fill is free to carry the estimator.
est <- bind_rows(
  sec_mean %>% transmute(label, x = log2_ratio, estimator = "cell2location abundance"),
  mk       %>% transmute(label, x = marker_log2, estimator = "marker transcript"))
est$estimator <- factor(est$estimator,
                        levels = c("cell2location abundance", "marker transcript"))

p_c <- ggplot(sec_mean, aes(x = log2_ratio, y = label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2_ratio, yend = label), colour = "grey75",
               linewidth = 0.7, lineend = "round") +
  geom_point(data = pts, aes(x = log2_ratio, y = label, colour = Condition),
             inherit.aes = FALSE, shape = 16, size = 1.3, alpha = 0.9) +
  scale_colour_manual(values = PAL_COND, name = NULL,
                      guide = guide_legend(order = 1)) +
  geom_point(data = est, aes(x = x, y = label, shape = estimator, fill = estimator),
             inherit.aes = FALSE, size = 2.0, colour = "black", stroke = 0.35) +
  scale_shape_manual(values = c(`cell2location abundance` = 23, `marker transcript` = 24),
                     name = NULL, guide = guide_legend(order = 2)) +
  scale_fill_manual(values = c(`cell2location abundance` = "grey25",
                               `marker transcript` = "white"),
                    name = NULL, guide = guide_legend(order = 2)) +
  scale_x_continuous(breaks = pretty_breaks(5),
                     limits = range(c(sec_mean$log2_ratio, pts$log2_ratio,
                                      mk$marker_log2), na.rm = TRUE) + c(-0.18, 0.18)) +
  labs(x = expression("log"[2]*" (NHD / CON), 5th-percentile abundance"), y = NULL) +
  theme_pub(base_size = 8) +
  theme(axis.text.y = element_text(size = 6.2, colour = "black"),
        axis.text.x = element_text(size = 6.2), axis.title.x = element_text(size = 6.6),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        axis.line = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "center", legend.justification = "center",
        legend.key.size = unit(0.26, "cm"), legend.spacing.y = unit(0.02, "cm"),
        legend.text = element_text(size = 6.2), legend.margin = margin(t = -2, b = 0),
        legend.box.spacing = unit(0.08, "cm"),
        plot.margin = margin(4, 5, 3, 3))
BN <- "F6f_c2l_composition"; W <- 3.10; H <- 2.45
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_c, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_c, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))

write.csv(sec_mean %>% select(cell_type = label, CON, NHD, log2_ratio) %>%
            left_join(mk %>% select(cell_type = label, marker_log2, n_gene), by = "cell_type"),
          file.path(TDIR, "figF6f_c2l_composition_FH.csv"), row.names = FALSE)
write.csv(med, file.path(TDIR, "fig6_depth_control_FH.csv"), row.names = FALSE)
cat("Wrote tables/fig6_{c2l_composition,depth_control}_FH.csv\n=== DONE ===\n")
