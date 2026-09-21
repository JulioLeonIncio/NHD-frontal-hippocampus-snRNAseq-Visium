#!/usr/bin/env Rscript
# =============================================================================
# 67b_visium_fig6e_alldepth_FH.R — Fig. 6e spatial maps of the two showcase programs on every spot at native depth (no 3,000-UMI floor), so the white matter is drawn.
# -----------------------------------------------------------------------------
# The depth-matched design of Fig. 6c–d and Supplementary Fig. 6 excludes nearly all white
# matter (control WM median ~200 UMI per spot). For the maps the biology comes first: each program is scored on all spots of the integrated object after
# per-spot log-normalisation (Seurat NormalizeData + AddModuleScore, ctrl = 100, seed 42, the
# same gene sets as script 65, read from tables/visium_dm_module_genes_FH.csv), centred and
# divided by its pooled SD (SD-within-pathway, as 67/73), drawn with the same compressive
# mapping, frame, ramp (white = low, dark = high) and layout as script 67. The depth caveat
# (shallow spots score noisily; control WM is ~200 UMI) goes in the legend. Overrides
# F6e1/F6e2 written by 67; the driver runs this step after 67.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2); library(ragg); library(ggh4x); library(patchwork) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))   # vis_outs(sec) — outs folder per section from VISIUM_OUTS_MANIFEST.json
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels"); TDIR <- file.path(PROJ, "tables")
VIS <- file.path(dirname(PROJ), "Visium", "integrated_harmony", "NHD_frontal_integrated_harmony.rds")
XY  <- file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN", "spot_coords.csv")
SEC <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")
SHOWCASE <- c("Complement / MHC-II", "Fatty-acid / sphingomyelin")
say <- function(...) cat(sprintf(...), "\n")

# ---- depth-free program signal on every spot -------------------------------------------
# Score = UMIs of the program's genes per 10,000 UMIs, pooled over the spot and its six hex
# neighbours (array coordinates; ratio of sums, so a single UMI in a ~200-UMI white-matter spot
# cannot dominate), log2(x + 1), pooled across sections. Per-domain means are printed and written so the
# legend can state them (control WM carries ~7x the grey-matter sphingomyelin signal; NHD WM does not).
CACHE <- file.path(PROJ, "data", "_cache_67b_alldepth_scores_FH.rds")
mg <- read.csv(file.path(TDIR, "visium_dm_module_genes_FH.csv"), stringsAsFactors = FALSE)
KEY <- paste("v3-umi", file.info(VIS)$mtime, file.info(XY)$mtime, paste(mg$gene[mg$module %in% SHOWCASE], collapse = ","))
if (file.exists(CACHE) && identical(readRDS(CACHE)$key, KEY)) { sc <- readRDS(CACHE)$sc; say("scores loaded from cache") } else {
  v <- readRDS(VIS); DefaultAssay(v) <- "Spatial"; v <- JoinLayers(v)
  cnt <- GetAssayData(v, layer = "counts"); tot <- Matrix::colSums(cnt)
  sc <- data.frame(spot_id = colnames(v), sample_id = v$sample_id, nCount = tot, stringsAsFactors = FALSE)
  for (nm in SHOWCASE) { g <- intersect(mg$gene[mg$module == nm], rownames(cnt))
    say("program %-26s %d genes: %s", nm, length(g), paste(g, collapse = " "))
    sc[[paste0(nm, " UMI")]] <- as.numeric(Matrix::colSums(cnt[g, , drop = FALSE])) }
  saveRDS(list(key = KEY, sc = sc), CACHE); rm(v, cnt); invisible(gc()) }
xy <- read.csv(XY, stringsAsFactors = FALSE) %>% transmute(spot_id = cell, x, y)
pos <- bind_rows(lapply(SEC, function(s) { f <- file.path(vis_outs(s), "spatial", "tissue_positions_list.csv")   # manifest outs
  p <- read.csv(f, header = FALSE, col.names = c("bc","it","ar","ac","pr","pc")); data.frame(spot_id = paste0(s, "_", p$bc), ar = p$ar, ac = p$ac) }))
sc <- sc %>% inner_join(xy, by = "spot_id") %>% inner_join(pos, by = "spot_id") %>% mutate(sample_id = factor(sample_id, levels = SEC))
stopifnot(nrow(sc) > 15000)
OFF <- list(c(0, 2), c(0, -2), c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
# Neighbourhood value = ratio of sums over the spot and its present neighbours (program UMIs / total
# UMIs x 10,000): deeper spots weigh more and a single UMI in a 200-UMI spot cannot blow up the ratio
nb_ratio <- function(dd, v) { key <- paste(dd$ar, dd$ac); idx <- setNames(seq_len(nrow(dd)), key)
  vapply(seq_len(nrow(dd)), function(i) { j <- idx[vapply(OFF, function(o) paste(dd$ar[i] + o[1], dd$ac[i] + o[2]), character(1))]
    k <- c(i, j[!is.na(j)]); 1e4 * sum(dd[[v]][k]) / sum(dd$nCount[k]) }, numeric(1)) }
for (nm in SHOWCASE) sc <- sc %>% group_by(sample_id) %>% mutate(!!paste0(nm, " (nb)") := nb_ratio(pick(everything()), paste0(nm, " UMI"))) %>% ungroup()
long <- tidyr::pivot_longer(sc, all_of(paste0(SHOWCASE, " (nb)")), names_to = "arm", values_to = "per10k") %>%
  mutate(arm = sub(" \\(nb\\)$", "", arm), score = log2(per10k + 1)) %>%
  group_by(sample_id) %>% mutate(.sp = max(diff(range(x)), diff(range(y))), xr = (x - mean(range(x)))/.sp, yr = (y - mean(range(y)))/.sp) %>% ungroup()
YSPAN <- long %>% group_by(sample_id) %>% summarise(sp = diff(range(yr)), .groups = "drop") %>% pull(sp) %>% max()
md <- read.csv(file.path(PROJ, "manuscript", "Final_fig_and_tables", "Data_deposition", "Visium", "spot_metadata.csv.gz"), stringsAsFactors = FALSE)
dom <- long %>% inner_join(md %>% select(spot_id, spatial_domain), by = "spot_id") %>% mutate(condition = ifelse(grepl("^CON", sample_id), "CON", "NHD")) %>%
  group_by(arm, condition, spatial_domain) %>% summarise(mean_per10k = mean(per10k), median_per10k = median(per10k), n = dplyr::n(), .groups = "drop")
write.csv(dom, file.path(TDIR, "fig6e_program_umi_per10k_by_domain_FH.csv"), row.names = FALSE); print(as.data.frame(dom %>% filter(spatial_domain %in% c("WM","L6","L2/3"))))
write.csv(long %>% select(spot_id, sample_id, arm, nCount, per10k, score), file.path(TDIR, "fig6e_alldepth_spot_scores_FH.csv"), row.names = FALSE)

spatial_one <- function(programme, basename) {
  dd <- long %>% filter(arm == programme); cap <- as.numeric(quantile(dd$score, 0.99)); dd <- dd %>% mutate(pos = pmin(score, cap))
  br <- pretty(c(0, cap), n = 3); br <- br[br <= cap]
  p <- ggplot(dd, aes(xr, yr, colour = pos)) + geom_point(size = 0.58, stroke = 0) +
    scale_colour_gradientn(colours = c("#FFFFFF", rev(viridisLite::magma(9, begin = 0.08, end = 0.97))), limits = c(0, cap),
                           breaks = br, labels = sprintf("%g", round(2^br - 1, 1)),
                           name = paste0(programme, " (UMI per 10,000)"),
                           guide = guide_colourbar(title.position = "left", title.vjust = 0.9, direction = "horizontal")) +
    facet_wrap(~ sample_id, nrow = 1, labeller = labeller(sample_id = function(x) sub("_Frontal", " ", x))) +
    coord_equal() + ggh4x::force_panelsizes(rows = unit(0.95 * YSPAN, "in"), cols = unit(0.95, "in")) +
    theme_pub(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.line = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          strip.text.x = element_text(size = 6.4, colour = "black"), strip.background = element_blank(), panel.spacing = unit(0.04, "cm"),
          legend.position = "top", legend.direction = "horizontal", legend.justification = "centre",
          legend.key.width = unit(1.1, "cm"), legend.key.height = unit(0.18, "cm"),
          legend.title = element_text(size = 5.8, colour = "black"), legend.text = element_text(size = 5.8, colour = "black"),
          legend.margin = margin(b = 0), legend.box.spacing = unit(0.06, "cm"), plot.margin = margin(1, 1, 1, 1))
  g <- ggplotGrob(p)   # measured canvas: fixed panels + strips + bar, no blank band
  W <- sum(grid::convertWidth(g$widths, "in", valueOnly = TRUE)) + 0.06
  H <- sum(grid::convertHeight(g$heights, "in", valueOnly = TRUE)) + 0.06
  ggsave(file.path(PANEL, paste0(basename, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
  ggsave(file.path(PANEL, paste0(basename, ".png")), p, width = W, height = H, dpi = 600, device = ragg::agg_png)
  say("Saved %s.{png,pdf} (%.2f x %.2f in) — %s, all %d spots", basename, W, H, programme, nrow(dd)) }
spatial_one(SHOWCASE[1], "F6e1_spatial_complement_MHCII")
spatial_one(SHOWCASE[2], "F6e2_spatial_fattyacid_sphingomyelin")

# ---- F6e: both programs in one panel, sized to sit beside panel d at 100 % ---------------
# maps MAPW in wide (four across), section labels once (top row), one horizontal bar per program
# under its row (the two programs have different ranges), everything at 100 % on the page.
MAPW <- 0.46
map_row <- function(programme, strips) {
  dd <- long %>% filter(arm == programme); cap <- as.numeric(quantile(dd$score, 0.99)); dd <- dd %>% mutate(pos = pmin(score, cap))
  br <- pretty(c(0, cap), n = 3); br <- br[br <= cap]
  ggplot(dd, aes(xr, yr, colour = pos)) + geom_point(size = 0.5, stroke = 0) +   # closes the hex grid at 0.52 in
    scale_colour_gradientn(colours = c("#FFFFFF", rev(viridisLite::magma(9, begin = 0.08, end = 0.97))), limits = c(0, cap),
                           breaks = br, labels = sprintf("%g", round(2^br - 1, 1)), name = paste0(programme, " (UMI per 10,000)"),
                           guide = guide_colourbar(title.position = "top", title.hjust = 0, direction = "horizontal")) +   # title above a short bar: nothing clips
    facet_wrap(~ sample_id, nrow = 1, labeller = labeller(sample_id = function(x) sub("_Frontal", " ", x))) +
    coord_equal() + ggh4x::force_panelsizes(rows = unit(MAPW * YSPAN, "in"), cols = unit(MAPW, "in")) +
    theme_pub(base_size = 7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_blank(),
          strip.text.x = if (strips) element_text(size = 6.2, colour = "black") else element_blank(), strip.background = element_blank(),
          panel.spacing = unit(0.03, "cm"), legend.position = "bottom", legend.direction = "horizontal",
          legend.key.width = unit(0.55, "cm"), legend.key.height = unit(0.16, "cm"), legend.justification = "left",
          legend.title = element_text(size = 5.8, colour = "black"), legend.text = element_text(size = 5.6, colour = "black"),   # row standard: labels 6, small text 5.6
          legend.margin = margin(t = 1, b = 0), legend.box.spacing = unit(0.12, "cm"), plot.margin = margin(1, 1, 1, 1)) }   # air between the maps and the bar
# Vertical reading: one column per program, the four sections stacked (CON 1, CON 2, NHD 1, NHD 2),
# section labels once on the left column, program title on top, one horizontal bar under each column
map_col <- function(programme, ylabels) {
  dd <- long %>% filter(arm == programme); cap <- as.numeric(quantile(dd$score, 0.99)); dd <- dd %>% mutate(pos = pmin(score, cap))
  br <- pretty(c(0, cap), n = 3); br <- br[br <= cap]
  ggplot(dd, aes(xr, yr, colour = pos)) + geom_point(size = 0.5, stroke = 0) +
    scale_colour_gradientn(colours = c("#FFFFFF", rev(viridisLite::magma(9, begin = 0.08, end = 0.97))), limits = c(0, cap),
                           breaks = c(0, cap), labels = c("0", sprintf("%g", round(2^cap - 1))), name = NULL,   # ends only; unit (UMI per 10,000) in the legend text
                           guide = guide_colourbar(direction = "horizontal", theme = theme(legend.key.width = unit(0.75, "cm"), legend.key.height = unit(0.15, "cm")))) +
    facet_wrap(~ sample_id, ncol = 1, strip.position = "left", labeller = labeller(sample_id = function(x) sub("_Frontal", " ", x))) +
    coord_equal() + ggh4x::force_panelsizes(rows = unit(MAPW * YSPAN, "in"), cols = unit(MAPW, "in")) +
    labs(title = sub(" / ", " /\n", programme)) + theme_pub(base_size = 7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line = element_blank(),
          strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = if (ylabels) element_text(size = 6.2, colour = "black", angle = 0, hjust = 1) else element_blank(),
          plot.title = element_text(size = 6.2, colour = "black", face = "plain", hjust = 0.5, lineheight = 0.9, margin = margin(b = 2)),
          panel.spacing = unit(0.04, "cm"), legend.position = "bottom", legend.direction = "horizontal", legend.justification = "centre",
          legend.key.width = unit(0.30, "cm"), legend.key.height = unit(0.15, "cm"), legend.location = "panel",   # bar centred under the maps, not under maps + row labels
          legend.title = element_text(size = 5.8, colour = "black"), legend.text = element_text(size = 5.6, colour = "black"),
          legend.margin = margin(t = 1, b = 0), legend.box.spacing = unit(0.12, "cm"), plot.margin = margin(1, 4, 1, 2)) }
pe <- map_col(SHOWCASE[1], TRUE) | map_col(SHOWCASE[2], FALSE)
g_e <- patchwork::patchworkGrob(pe)
We <- sum(grid::convertWidth(g_e$widths, "in", valueOnly = TRUE)) + 0.06
He <- sum(grid::convertHeight(g_e$heights, "in", valueOnly = TRUE)) + 0.06
ggsave(file.path(PANEL, "F6e_spatial_maps.pdf"), pe, width = We, height = He, useDingbats = FALSE)
ggsave(file.path(PANEL, "F6e_spatial_maps.png"), pe, width = We, height = He, dpi = 600, device = ragg::agg_png)
say("Saved F6e_spatial_maps.{png,pdf} (%.2f x %.2f in): both programs, place at 100 %% beside d", We, He)
writeLines(capture.output(sessionInfo()), file.path(PROJ, "logs", "67b_visium_fig6e_alldepth_FH_sessionInfo.txt"))
cat("\n=== DONE ===\n", file = stderr())
