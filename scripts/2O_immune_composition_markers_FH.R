#!/usr/bin/env Rscript
# =============================================================================
# 2O_immune_composition_markers_FH.R — NHD immune compartment: proportions + three-population marker panel (Microglia / PVM / CD8 T-lymphocyte).
# -----------------------------------------------------------------------------
# "people will like to see some sort of proportion analysis
# regarding NHD immune cells", plus the go-ahead for a dedicated T-cell marker panel.
#
# Two DENOMINATORS, because they answer different questions and only reporting one
# would be misleading:
#   (a) % of the immune compartment  — how the immune compartment is re-composed.
#       This alone would overstate PVM: microglia collapse in NHD, so every other
#       immune population's share rises even if its own abundance is unchanged.
#   (b) per 10,000 total nuclei      — abundance against the whole tissue, which
#       separates "PVM rose" from "microglia fell". Still a compositional measure
#       (snRNA gives no absolute density), so it is never described as cell density.
#
# Stats: n = 1 donor per condition. The fraction shift is the effect size; no test,
# no stars [[composition-panels-effectsize-not-pseudorep]].
#
# Markers are DGE-derived for all three populations (ribosomal/mitochondrial genes
# excluded from the test universe), ranked by detection difference with a ceiling on
# the other populations — PVM and T nuclei differ in depth from microglia, so
# log2FC-ranking would surface depth-driven housekeeping genes.
#
# Outputs: figures/Figure_2/panels/2O_immune_composition.{png,pdf}
#          figures/Figure_2/panels/2P_immune_markers.{png,pdf}
#          tables/micro_states/immune_{composition,markers}_FH.csv
#          data/_cache_immune_compartment_FH.rds
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(ggh4x)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_artifact_genes.R"))
ATLAS  <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
PANEL  <- file.path(PROJ, "figures", "Figure_2", "panels")
TBL    <- file.path(PROJ, "tables", "micro_states")
CACHE  <- file.path(PROJ, "data", "_cache_immune_compartment_FH.rds")
CLS    <- file.path(TBL, "immune_umap_classes_FH.csv")
stopifnot(file.exists(ATLAS), "run 2N_immune_UMAP_FH.R first" = file.exists(CLS))

PAL_IMM <- c("Microglia"    = unname(PAL_CELLTYPE[["Micro-PVM"]]),
             "PVM"          = "#6D4C41",
             "T-lymphocyte" = "#C77B0A")
IMM_LEVELS <- names(PAL_IMM)

# The whole-tissue denominator needs the atlas, but the immune object itself can come
# from cache. Cache is mtime-guarded against the atlas so a rebuilt atlas invalidates it.
TOTCSV <- file.path(TBL, "total_nuclei_by_region_condition_FH.csv")
cache_ok <- file.exists(CACHE) && file.exists(TOTCSV) &&
            file.info(CACHE)$mtime >= file.info(ATLAS)$mtime
if (cache_ok) {
  cat("== using cached immune compartment (no atlas load) ==\n")
  o   <- readRDS(CACHE)
  tot <- read.csv(TOTCSV)
} else {
  cat("== cache miss/stale: building immune compartment from atlas ==\n")
  atl <- readRDS(ATLAS)
  if (inherits(atl[["RNA"]], "Assay5")) atl[["RNA"]] <- SeuratObject::JoinLayers(atl[["RNA"]])
  # Denominator: all nuclei per condition x region
  tot <- atl@meta.data %>% count(Region, Condition, name = "n_total")
  write.csv(tot, TOTCSV, row.names = FALSE)
  cls <- read.csv(CLS)                       # barcode, class  (from 2N)
  o <- subset(atl, cells = cls$barcode); rm(atl); gc()
  o$immune_class <- factor(cls$class[match(colnames(o), cls$barcode)], levels = IMM_LEVELS)
  o$Condition <- factor(o$Condition, levels = c("CON","NHD"))
  o$Region    <- factor(o$Region,    levels = REGION_ORDER)
  stopifnot("class assignment failed" = !any(is.na(o$immune_class)))
  saveRDS(o, CACHE)
}
tot$Region <- factor(tot$Region, levels = REGION_ORDER)
tot$Condition <- factor(tot$Condition, levels = c("CON","NHD"))
cat(sprintf("immune compartment: %d nuclei\n", ncol(o)))
print(table(o$immune_class, o$Condition, o$Region))

# ---- 1. composition --------------------------------------------------------
comp <- o@meta.data %>%
  count(Region, Condition, immune_class, name = "n") %>%
  complete(Region, Condition, immune_class, fill = list(n = 0L)) %>%
  group_by(Region, Condition) %>% mutate(n_immune = sum(n)) %>% ungroup() %>%
  left_join(tot, by = c("Region","Condition")) %>%
  mutate(pct_immune = 100 * n / n_immune,
         per_10k_total = 1e4 * n / n_total)
write.csv(comp, file.path(TBL, "immune_composition_FH.csv"), row.names = FALSE)
cat("\n=== immune composition ===\n"); print(as.data.frame(comp), digits = 3)

# (a) PRIMARY: slope on a log10 axis. All three
# populations and both regions in one small panel; the log axis is what lets
# microglia (~1000 per 10k) sit beside T cells (~5 per 10k) without flattening
# either, and the slope carries the message — microglia down, PVM up, T flat.
p_slope <- ggplot(comp, aes(x = Condition, y = per_10k_total,
                            colour = immune_class, group = immune_class)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.5) +
  facet_wrap2(~ Region, nrow = 1, labeller = region_full_lowconf_labeller(),
              strip = strip_region_x(REGION_ORDER, clip = "off")) +
  scale_colour_manual(values = PAL_IMM, labels = IMM_LABEL, name = NULL) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300, 1000),
                labels = c("3","10","30","100","300","1000")) +
  labs(x = NULL, y = "nuclei per 10,000 total (log scale)") +
  # This is the panel the assembled Fig3_suppl.pdf actually embeds (verified by
  # extracting the placed 1800x1320 raster from the PDF and diffing it against both
  # candidates). Its twin `_alternatives/ALT1_immune_slope_log10` renders the same
  # plot from `2O_alt_composition_designs_FH.R` and had already been corrected to
  # BASE_ALT1 = 6.7 -- so the correction was landing on the copy the figure does not
  # use, leaving this panel and its row-mate ALT2 at two different type sizes.
  # 6.7 x the 0.758 placement scale = 5.08 pt on the page, the house band.
  theme_pub(base_size = 6.7) +
  theme(legend.position = "bottom", legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(t = -2), panel.spacing = unit(6, "pt"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.2))
ggsave(file.path(PANEL, "2O_immune_composition.png"), p_slope,
       width = 3.0, height = 2.2, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "2O_immune_composition.pdf"), p_slope,
       width = 3.0, height = 2.2, bg = "white")
cat("wrote 2O_immune_composition (slope, primary)\n")

# (a-alt) stacked share — kept, but demoted to _alternatives now that the slope
# design is the shipped panel. On its own the stacked view OVERSTATES PVM, because
# microglial collapse inflates every other population's share.
ALT <- file.path(PANEL, "_alternatives"); dir.create(ALT, showWarnings = FALSE, recursive = TRUE)
p_share <- ggplot(comp, aes(x = Condition, y = pct_immune, fill = immune_class)) +
  geom_col(width = 0.66, linewidth = 0) +
  facet_wrap2(~ Region, nrow = 1, labeller = region_full_lowconf_labeller(),
              strip = strip_region_x(REGION_ORDER, clip = "off")) +
  scale_fill_manual(values = PAL_IMM, labels = IMM_LABEL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "% of immune nuclei") +
  theme_pub(base_size = 7) +
  theme(legend.position = "bottom", legend.key.size = unit(0.26, "cm"),
        legend.margin = margin(t = -2), panel.spacing = unit(6, "pt"))

# (b) abundance against all nuclei — free y, because microglia are ~2 orders of
# magnitude more abundant than T cells and a shared axis would flatten both.
p_abund <- ggplot(comp, aes(x = Condition, y = per_10k_total, fill = Condition)) +
  geom_col(width = 0.62, linewidth = 0) +
  facet_grid2(immune_class ~ Region, scales = "free_y", independent = "y",
              labeller = labeller(Region = REGION_FULL),
              strip = strip_region_x(REGION_ORDER, clip = "off")) +
  scale_fill_manual(values = PAL_COND, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "nuclei per 10,000 total nuclei") +
  theme_pub(base_size = 7) +
  theme(panel.spacing = unit(5, "pt"),
        strip.text.y = element_text(size = 6.2, angle = -90))

ggsave(file.path(ALT, "ALT3b_immune_stacked_share.png"), p_share,
       width = 2.6, height = 2.1, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(ALT, "ALT3b_immune_stacked_share.pdf"), p_share,
       width = 2.6, height = 2.1, bg = "white")
ggsave(file.path(PANEL, "2O2_immune_abundance_per10k.png"), p_abund,
       width = 2.7, height = 3.1, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "2O2_immune_abundance_per10k.pdf"), p_abund,
       width = 2.7, height = 3.1, bg = "white")
cat("wrote _alternatives/ALT3b_immune_stacked_share + 2O2_immune_abundance_per10k\n")

# ---- 2. three-population marker DGE ---------------------------------------
DefaultAssay(o) <- "RNA"
o <- NormalizeData(o, verbose = FALSE)
if (inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- SeuratObject::JoinLayers(o[["RNA"]])
RIBO_MITO_RE <- "^RP[LS]|^MRP[LS]|^MT-|^MTRNR|^MTATP|^MTCO|^MTND|^MTCYB"
feat_use <- setdiff(rownames(o), grep(RIBO_MITO_RE, rownames(o), value = TRUE))
cat(sprintf("\nDGE universe: %d genes (ribosomal/mitochondrial excluded)\n", length(feat_use)))

Idents(o) <- o$immune_class
am <- FindAllMarkers(o, features = feat_use, only.pos = TRUE, min.pct = 0.15,
                     logfc.threshold = 0.25, test.use = "wilcox", verbose = FALSE) %>%
  mutate(is_artifact = is_artifact(gene, "Micro-PVM"), pct_diff = pct.1 - pct.2)
write.csv(am, file.path(TBL, "immune_markers_FH.csv"), row.names = FALSE)

N_MK <- 6
top <- am %>%
  filter(!is_artifact, p_val_adj < 0.05, pct.1 >= 0.25, pct.2 <= 0.20) %>%
  group_by(cluster) %>% slice_max(pct_diff, n = N_MK, with_ties = FALSE) %>% ungroup()
cat("\n=== display markers per population (detection difference) ===\n")
print(as.data.frame(top %>% select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj)), digits = 3)

GENES <- top %>% arrange(factor(cluster, levels = IMM_LEVELS), desc(pct_diff)) %>% pull(gene)
expr  <- FetchData(o, vars = unique(GENES), layer = "data")
expr$class <- o$immune_class
dot <- expr %>%
  pivot_longer(-class, names_to = "gene", values_to = "e") %>%
  group_by(class, gene) %>%
  summarise(pct_exp = 100 * mean(e > 0), avg_exp = mean(e), .groups = "drop") %>%
  mutate(marks = factor(top$cluster[match(gene, top$gene)], levels = IMM_LEVELS),
         gene  = factor(gene, levels = unique(GENES)),
         class = factor(class, levels = rev(IMM_LEVELS)))

p_dot <- ggplot(dot, aes(x = gene, y = class)) +
  geom_point(aes(size = pct_exp, colour = avg_exp)) +
  facet_grid2(~ marks, scales = "free_x", space = "free_x") +
  scale_size_continuous(range = c(0.3, 3.2), name = "% expressing",
                        breaks = c(25, 50, 75)) +
  scale_colour_gradient(low = "grey90", high = "#B3132B", name = "mean expression") +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(size = 6.4),
        legend.position = "right", legend.key.size = unit(0.28, "cm"),
        legend.title = element_text(size = 6.2), legend.text = element_text(size = 6),
        panel.grid.major = element_line(colour = "grey94", linewidth = 0.2),
        plot.margin = margin(7, 4, 3, 3))
ggsave(file.path(PANEL, "2P_immune_markers.png"), p_dot,
       width = 5.0, height = 1.9, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "2P_immune_markers.pdf"), p_dot,
       width = 5.0, height = 1.9, bg = "white")
cat("wrote 2P_immune_markers.{png,pdf}\n")
cat("\n=== DONE: 2O_immune_composition_markers_FH ===\n")
