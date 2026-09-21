#!/usr/bin/env Rscript
# =============================================================================
# 74_visium_fig6_abc_FH.R — Figure 6 panels a, b, c, rebuilt inside this project.
#   F6a_visium_domain_umap      spot UMAP coloured by laminar domain
#   F6b_visium_domain_spatial   the same domains in tissue, 4 sections
#   F6c_laminar_marker_profile  layer markers across the grey-matter domains (superficial -> deep), CON vs NHD
# -----------------------------------------------------------------------------
# Why these panels live in this project. They are produced here, linked (not embedded)
# in the composite so they refresh with the data, and drawn on the same depth-matched
# prep (65_visium_prep_depthmatched_FH.R) as the rest of Figure 6: without depth
# matching, upper-layer markers such as CUX2 appear to flatten in NHD purely through
# library size (65 quantifies this as "entirely artifact").
#
# Panel c includes white matter and is therefore computed on every spot at
# native depth with a depth-free estimator: per spot, gene UMIs per 10,000 UMIs (a ratio, unbiased
# at any depth), z-scored per gene across all spots, then mean ± s.e.m. per domain and condition.
# Control white matter is shallow (median ~200 UMI per spot; tables/fig6_domain_depth_asymmetry_FH.csv),
# so its per-spot values are noisy, but the mean over ~900 spots is not; this is the same estimator
# as the panel-e maps (67b). Domains run L2/3 > L4 > L5 > L6 > WM (depth_rank order, then WM).
#
# Therefore panel c is drawn over L2/3 -> L6 only, on depth-matched spots (every spot
# thinned to a common 3,000 UMI by 65), and the WM asymmetry is reported as its own
# number for the legend rather than smuggled into a curve. That asymmetry is probably
# real biology -- gliotic, cell-rich NHD white matter against hypocellular control
# white matter -- but it is a cellularity statement, not a laminar-identity one.
#
# n = 2 sections per condition, one donor per condition. Descriptive only; no p-values.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggh4x); library(patchwork)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels")
C2L   <- file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN")

# cluster -> domain map, level order, palette and the laminar (grey-matter) axis all come from
# integrated_harmony/integrated_{cluster_identity,domain_levels}.csv via _visium_domains_FH.R
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))
LAM     <- DOM_GREY                                          # depth-matchable domains, superficial -> deep
SEC     <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")

meta <- read.csv(file.path(C2L, "visium_meta.csv"), stringsAsFactors = FALSE)
names(meta)[1] <- "spot_id"
xy   <- read.csv(file.path(C2L, "spot_coords.csv"), stringsAsFactors = FALSE)
meta <- meta %>%
  mutate(domain = factor(unname(DOM_MAP[as.character(seurat_clusters)]), levels = DOM_LEV),
         sample_id = factor(sample_id, levels = SEC)) %>%
  left_join(xy %>% transmute(spot_id = cell, x, y), by = "spot_id") %>%
  filter(!is.na(domain), !is.na(x))
cat(sprintf("spots with a domain and coordinates: %d\n", nrow(meta)))

# the WM depth asymmetry, computed once and written out for the legend
wm <- meta %>% group_by(domain, condition) %>%
  summarise(n = n(), median_umi = median(nCount_Spatial),
            n_ge3000 = sum(nCount_Spatial >= 3000), .groups = "drop")
cat("\n== per-spot depth by domain (this is why panel c stops at L6) ==\n")
print(as.data.frame(wm), row.names = FALSE)
write.csv(wm, file.path(TDIR, "fig6_domain_depth_asymmetry_FH.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# b — domains in tissue (built from coordinates + cluster labels; no Seurat object)
# ---------------------------------------------------------------------------
sp <- meta %>% group_by(sample_id) %>%
  mutate(.sp = max(diff(range(x)), diff(range(y))),
         xr = (x - mean(range(x)))/.sp, yr = (y - mean(range(y)))/.sp) %>% ungroup()
YSPAN <- max(tapply(sp$yr, sp$sample_id, function(v) diff(range(v))))

p_b <- ggplot(sp, aes(xr, yr, colour = domain)) +
  geom_point(size = 0.42, stroke = 0) +
  scale_colour_manual(values = DOM_PAL, name = NULL,
                      guide = guide_legend(override.aes = list(size = 1.9), nrow = 1)) +
  facet_wrap(~ sample_id, nrow = 1,
             labeller = labeller(sample_id = function(x) sub("_Frontal", " ", x))) +
  coord_equal() +
  ggh4x::force_panelsizes(rows = unit(0.95 * YSPAN, "in"), cols = unit(0.95, "in")) +
  theme_pub(base_size = 8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.line = element_blank(),
        strip.text.x = element_text(size = 6.4, colour = "black"),
        strip.background = element_blank(), panel.spacing = unit(0.04, "cm"),
        legend.position = "bottom", legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 5.8), legend.margin = margin(t = -2),
        legend.box.spacing = unit(0.06, "cm"), plot.margin = margin(1, 1, 1, 1))
BN <- "F6b_visium_domain_spatial"; W <- 4.35; H <- 0.95 * YSPAN + 0.58
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_b, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_b, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nSaved %s (%.2f x %.2f in)\n", BN, W, H))

# ---------------------------------------------------------------------------
# c — laminar marker profile, depth-matched, grey-matter domains in depth_rank order
# ---------------------------------------------------------------------------
GEN <- c("CUX2" = "L2/3", "RORB" = "L4", "PCP4" = "L5", "TLE4" = "L6", "MBP" = "myelin")
GEN_PAL <- setNames(unname(DOM_PAL[c("L2/3","L4","L5","L6","WM")]), c("CUX2","RORB","PCP4","TLE4","MBP")); stopifnot(!anyNA(GEN_PAL))
LAMW <- c(LAM, "WM")
# gene UMIs per 10,000 on every spot (cached on the object's mtime)
CACHE_C <- file.path(PROJ, "data", "_cache_74_panelc_umi_FH.rds")
IH_RDS <- file.path(dirname(PROJ), "Visium", "integrated_harmony", "NHD_frontal_integrated_harmony.rds")
KEY_C <- paste(file.info(IH_RDS)$mtime, paste(names(GEN), collapse = ","))
if (file.exists(CACHE_C) && identical(readRDS(CACHE_C)$key, KEY_C)) { dc <- readRDS(CACHE_C)$dc } else {
  suppressPackageStartupMessages(library(Seurat)); vv <- readRDS(IH_RDS); DefaultAssay(vv) <- "Spatial"; vv <- JoinLayers(vv)
  cnt <- GetAssayData(vv, layer = "counts"); tot <- Matrix::colSums(cnt); g <- intersect(names(GEN), rownames(cnt))
  dc <- data.frame(spot_id = colnames(vv), Condition = vv$condition, stringsAsFactors = FALSE)
  for (gg in g) dc[[paste0("gene_", gg)]] <- as.numeric(cnt[gg, ] / tot * 1e4)
  saveRDS(list(key = KEY_C, dc = dc), CACHE_C); rm(vv, cnt); invisible(gc()) }
d <- dc %>% left_join(meta %>% select(spot_id, domain), by = "spot_id") %>%
  filter(domain %in% LAMW) %>% mutate(domain = factor(domain, levels = LAMW))
have <- intersect(paste0("gene_", names(GEN)), names(d))
stopifnot("layer genes missing from the object" = length(have) >= 4)
cat(sprintf("\npanel c: %d spots (all, native depth) across %s\n", nrow(d), paste(LAMW, collapse = " > ")))
print(table(d$domain, d$Condition))

prof <- d %>% select(Condition, domain, all_of(have)) %>%
  pivot_longer(all_of(have), names_to = "gene", values_to = "e") %>%
  mutate(gene = sub("^gene_", "", gene)) %>%
  group_by(gene) %>% mutate(z = (e - mean(e)) / stats::sd(e)) %>%   # one z per gene
  group_by(gene, Condition, domain) %>%
  summarise(m = mean(z), se = stats::sd(z)/sqrt(n()), .groups = "drop") %>%
  mutate(Condition = factor(Condition, levels = c("CON","NHD")))
cat("\n== panel c, mean z by domain ==\n")
print(as.data.frame(prof %>% select(gene, Condition, domain, m) %>%
      mutate(m = round(m, 2)) %>% pivot_wider(names_from = domain, values_from = m)),
      row.names = FALSE)

p_c <- ggplot(prof, aes(domain, m, colour = gene, group = gene)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se, fill = gene), alpha = 0.18,
              colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.5) + geom_point(size = 0.9, stroke = 0) +
  # each trace wears the hue of the layer it marks (CUX2 L2/3, RORB L4, PCP4 L5, TLE4 L6, MBP WM), read from the
  # domain palette, so panel b keys onto a/c instead of borrowing Dark2 hues with the wrong meaning
  scale_colour_manual(values = GEN_PAL, name = NULL) +
  scale_fill_manual(values = GEN_PAL) +
  facet_wrap(~ Condition, ncol = 1) +   # CON over NHD
  labs(x = "spatial domain", y = "z (UMI per 10,000)") +
  theme_pub(base_size = 7.4) +          # same type-to-panel ratio as panel a, placed at 100 %
  theme(axis.text = element_text(colour = "black"), strip.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),   # seven domains on a 0.90-in axis — horizontal labels overprinted (L1/pia|L2/3, L6|WM)
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        axis.line = element_blank(), panel.spacing.y = unit(0.22, "cm"),
        legend.position = "bottom", legend.key.size = unit(0.25, "cm"), legend.margin = margin(t = -2), legend.key.spacing.x = unit(2, "pt"), legend.location = "plot",
        legend.text = element_text(face = "italic"), plot.margin = margin(3, 3, 2, 3)) +   # gene symbols italic; legend at the bottom
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +   # five genes: 3 + 2 fits the 1.79-in plot width
  ggh4x::force_panelsizes(rows = unit(0.72, "in"), cols = unit(0.90, "in"))   # same plot height as panel a
g_c <- ggplotGrob(p_c)
W <- sum(grid::convertWidth(g_c$widths, "in", valueOnly = TRUE)) + 0.22
H <- sum(grid::convertHeight(g_c$heights, "in", valueOnly = TRUE)) + 0.10
BN <- "F6c_laminar_marker_profile"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_c, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_c, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))
write.csv(prof, file.path(TDIR, "fig6c_laminar_profile_FH.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# a — composite UMAP, matching the layout of the panel it replaces:
#     TOP    spots coloured by laminar domain, legend beneath in two rows
#     bottom the same embedding coloured by condition, legend beneath in one row
#     Canvas is 3.90 x 4.70 in (2340 x 2820 px), identical to the old panel, so it
#     drops into the existing Illustrator frame without re-placing.
# The embedding is cached on first run so later edits to this panel do not have to
# reload the 528 MB integrated object.
# ---------------------------------------------------------------------------
EMB <- file.path(PROJ, "data", "_cache_visium_umap_FH.csv")
IH  <- file.path(dirname(PROJ), "Visium", "integrated_harmony",
                 "NHD_frontal_integrated_harmony.rds")
if (!file.exists(EMB) ||
    (file.exists(IH) && file.info(EMB)$mtime < file.info(IH)$mtime)) {
  stopifnot("integrated object not found and no cached embedding" = file.exists(IH))
  suppressPackageStartupMessages(library(Seurat))
  cat("\nloading the integrated object for the UMAP embedding ...\n")
  o  <- readRDS(IH)
  em <- as.data.frame(Seurat::Embeddings(o, "umap"))[, 1:2]
  names(em) <- c("UMAP_1", "UMAP_2"); em$spot_id <- rownames(em)
  rm(o); invisible(gc())
  write.csv(em, EMB, row.names = FALSE)
  cat("cached the embedding ->", basename(EMB), "\n")
} else {
  em <- read.csv(EMB, stringsAsFactors = FALSE)
  cat("\nusing the cached UMAP embedding\n")
}

u <- em %>% left_join(meta %>% select(spot_id, domain, condition), by = "spot_id") %>%
  filter(!is.na(domain))
set.seed(42); u <- u[sample(nrow(u)), ]
u$condition <- factor(u$condition, levels = c("CON","NHD"))
cat(sprintf("panel a: %d spots\n", nrow(u)))

umap_base <- function(p, base = 7.4)   # base 7.4 on 0.9 x 0.72-in plots, placed at 100 % (axis text 6, titles 6.8: the row standard)
  p + coord_equal() +
      labs(x = "UMAP 1", y = "UMAP 2") +
      theme_pub(base_size = base) +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
            axis.line = element_blank(),
            axis.text = element_text(size = base - 1.4, colour = "black"),
            axis.title = element_text(size = base - 0.6),
            legend.position = "bottom", legend.key.size = unit(0.30, "cm"), legend.location = "plot",   # centred on the whole canvas, not the panel
            legend.text = element_text(size = base - 0.6),
            legend.margin = margin(t = -2, b = 0),
            legend.box.spacing = unit(0.10, "cm"),
            plot.margin = margin(2, 3, 1, 2))

p_dom <- umap_base(
  ggplot(u, aes(UMAP_1, UMAP_2, colour = domain)) +
    geom_point(size = 0.26, stroke = 0, alpha = 0.85) +
    scale_colour_manual(values = DOM_PAL, name = NULL,
                        guide = guide_legend(ncol = 3, byrow = TRUE,   # legend cells are equal-width (widest = "Vasc/immune"): 4 per row overflows the panel width, 3 per row reads L2/3 L4 L5 / L6 InN WM / Vasc/immune
                                             override.aes = list(size = 1.9, alpha = 1))))

p_cond <- umap_base(
  ggplot(u, aes(UMAP_1, UMAP_2, colour = condition)) +
    geom_point(size = 0.26, stroke = 0, alpha = 0.75) +
    scale_colour_manual(values = PAL_COND_PALE, name = NULL,   # dense UMAP points = the pale pair
                        guide = guide_legend(nrow = 1,
                                             override.aes = list(size = 1.9, alpha = 1))))

# No dead margin. With coord_equal and an embedding 1.45x wider than tall, two stacked
# sub-panels are height-limited and every surplus inch of width becomes blank margin
# (measured at up to 21% of a 3.16 in panel otherwise; the Visium spatial maps behave the
# same way). The remedy is geometric: fix the panel size to the data's own aspect and cut
# the canvas to fit it, so the width the panel is given is the width it can actually use.
ASP  <- diff(range(u$UMAP_1)) / diff(range(u$UMAP_2))   # ~1.45, from the data
PH_A <- 0.72                                            # per-panel plot height, inches: a is placed at 100 % in a ~1.5-in column (type-to-panel ratio)
PW_A <- PH_A * ASP
p_dom  <- p_dom  + ggh4x::force_panelsizes(rows = unit(PH_A, "in"), cols = unit(PW_A, "in"))
p_cond <- p_cond + ggh4x::force_panelsizes(rows = unit(PH_A, "in"), cols = unit(PW_A, "in"))
p_a <- patchwork::wrap_plots(p_dom, p_cond, ncol = 1)
# canvas = the measured gtable (fixed panels + furniture) plus a small margin — a legend that
# gains a row can no longer clip the key or the top tick (measured, not a hand-tuned constant)
g_a <- patchwork::patchworkGrob(p_a)
W <- sum(grid::convertWidth(g_a$widths, "in", valueOnly = TRUE)) + 0.22   # + room for the 3-column key (gtable widths do not include legend overhang)
H <- sum(grid::convertHeight(g_a$heights, "in", valueOnly = TRUE)) + 0.10
BN <- "F6a_visium_domain_umap"
cat(sprintf("panel a: embedding aspect %.2f -> plot %.2f x %.2f in each; canvas %.2f x %.2f\n",
            ASP, PW_A, PH_A, W, H))
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_a, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p_a, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s (%.2f x %.2f in)\n", BN, W, H))
cat("=== DONE ===\n")
