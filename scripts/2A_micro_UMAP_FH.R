#!/usr/bin/env Rscript
# =============================================================================
# 2A_micro_UMAP_FH.R — Microglia-subset UMAP panels for Figure 2, frontal + HIPPOCAMPUS (2 regions). Port of 48_fig2_micro_UMAP.R.
# -----------------------------------------------------------------------------
# Re-extracts the POOLED myeloid compartment (microglia + the CD163+/F13A1+ population
# labelled "PVM" in the class file; T-lymphocytes and doublets excluded) and re-embeds it.
# The CD163+/F13A1+ cells are not asserted to be peripheral macrophages -- CD163 can mark
# a microglial state, and they do not separate in this embedding (they sit inside the two
# main microglial clusters: 65/1335 in cluster 0, 96/1022 in cluster 1, 1 in cluster 2).
# (harmony over SampleID, like the atlas) so the UMAP reflects the current data.
# Renders:
#   2_micro_UMAP_condition        single UMAP, CON vs NHD (pale pair)
#   2_micro_UMAP_condition_split  CON | NHD facets, cells coloured by cluster
#   2_micro_UMAP_clusters         microglial clusters (states) 0..k
#   2_micro_UMAP_region           coloured by Region (Frontal/Hippo)
#   F2a_micro_UMAP_composite       condition | region stacked (panel a)
#
# Cache (mtime-guarded, read by 109_suppfig5_split_umaps_FH.R):
#   data/_cache_2A_myeloid_umap_FH.rds   list(df = barcode/UMAP_1/UMAP_2/Condition/Region/
#                                        seurat_clusters/myeloid_class, provenance = ...)
#
# ragg PNG + base pdf.  No Cairo/ggrastr.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(scales)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, PAL_REGION_UMAP, REGION_FULL, REGION_ORDER

# CON/NHD tints now come from the theme (PAL_COND_PALE) — see role-tint rule there.

suppressPackageStartupMessages(library(harmony))
ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
PANEL <- file.path(PROJ, "figures/Figure_2/panels")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")

stopifnot("MISSING atlas" = file.exists(ATLAS))
# Pooled myeloid: microglia and perivascular
# macrophages are analysed together in every downstream analysis, so this lead UMAP must
# show the same population the DGE was run on. Previously this panel was cut
# to microglia-only to fix a category error — the panel was titled "Microglia" while built
# on pooled cells. With pooling now the primary, the fix runs the other way: pool the
# cells and title the panel for what it contains.
# T-LYMPHOCYTES are excluded — they are a separate lineage, not part of the myeloid
# compartment the DGE covers, and they are shown in their own right in 2N_immune_UMAP.
MYELOID_CLASSES <- c("Microglia", "PVM")
.cls_f <- file.path(PROJ, "tables", "micro_states", "immune_umap_classes_FH.csv")
if (!file.exists(.cls_f)) stop("MISSING ", .cls_f, " — run 2N_immune_UMAP_FH.R first")
.cls <- read.csv(.cls_f, stringsAsFactors = FALSE)
cat("immune classes on file: ",
    paste(sprintf("%s=%d", names(table(.cls$class)), as.integer(table(.cls$class))),
          collapse = "  "), "\n", sep = "")

# ---- EMBEDDING ------------------------------------------------
# The exact embedding this panel draws (UMAP coords + cluster + myeloid class + Condition/
# Region per barcode) is cached with an MTIME guard, same pattern as 2N
# ([[atlas-single-source-cache-guards]]): the cache is reused only if it is newer than both
# inputs (the atlas and the immune-class CSV). Touch either and it rebuilds; delete the
# cache to force a rebuild. Downstream readers (109_suppfig5_split_umaps_FH.R, row b) take
# their coordinates from this cache, so the supplementary split is Fig. 2a, not a re-embed.
# On a miss, the SCT scale.data for the myeloid nuclei is taken from the immune-compartment
# cache (data/_cache_immune_compartment_FH.rds, itself mtime-guarded against the atlas and
# written by 2O as `subset(atlas, cells = <immune barcodes>)`): scale.data is only column-
# subset by subset(), so the PCA input matrix is identical to subsetting the 2.9 GB atlas,
# and the 10 GB atlas read is skipped. The atlas path is kept as the fallback.
UCACHE <- file.path(PROJ, "data", "_cache_2A_myeloid_umap_FH.rds")
IMMC   <- file.path(PROJ, "data", "_cache_immune_compartment_FH.rds")
.cache_ok <- file.exists(UCACHE) &&
  file.mtime(UCACHE) > file.mtime(ATLAS) &&
  file.mtime(UCACHE) > file.mtime(.cls_f)
if (.cache_ok) {
  .cc <- readRDS(UCACHE)
  stopifnot("bad 2A embedding cache" = is.list(.cc) && is.data.frame(.cc$df) &&
              all(c("barcode","UMAP_1","UMAP_2","Condition","Region","seurat_clusters","myeloid_class") %in% names(.cc$df)))
  cat(sprintf("== embedding CACHE HIT (%s): %d nuclei ==\n", format(file.mtime(UCACHE), "%Y-%m-%d %H:%M"), nrow(.cc$df)))
  o <- NULL
} else {
  cat("== embedding CACHE MISS -- recomputing ==\n")
  .keep_all <- .cls$barcode[.cls$class %in% MYELOID_CLASSES]
  .use_immc <- file.exists(IMMC) && file.mtime(IMMC) > file.mtime(ATLAS)
  if (.use_immc) {
    src <- readRDS(IMMC)
    .use_immc <- all(.keep_all %in% colnames(src)) &&
      "SCT" %in% Assays(src) && nrow(GetAssayData(src, assay = "SCT", layer = "scale.data")) > 0
    if (!.use_immc) { cat("   immune-compartment cache unusable for this set -> falling back to the atlas\n"); rm(src) }
  }
  if (!.use_immc) { cat("   reading the full atlas (slow, ~10 GB)\n"); src <- readRDS(ATLAS) }
  else            cat(sprintf("   source = %s (%s)\n", basename(IMMC), format(file.mtime(IMMC), "%Y-%m-%d %H:%M")))
  .keep <- intersect(.keep_all, colnames(src))
  if (length(.keep) < 100) stop("myeloid assignment matched only ", length(.keep), " cells")
  o   <- subset(src, cells = .keep); rm(src); gc()
  # carry the class through so the composition can be reported (and audited) from the panel
  o$myeloid_class <- .cls$class[match(colnames(o), .cls$barcode)]
  .n_pvm <- sum(o$myeloid_class == "PVM")
  cat(sprintf("POOLED myeloid nuclei: %d (microglia %d + PVM %d = %.1f%% PVM); T-lymphocytes excluded\n",
              ncol(o), ncol(o) - .n_pvm, .n_pvm, 100 * .n_pvm / ncol(o)))
  DefaultAssay(o) <- "SCT"
  stopifnot("SCT scale.data empty — re-run SCTransform/PrepSCT before embedding" =
              nrow(GetAssayData(o, assay = "SCT", layer = "scale.data")) > 0)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- harmony::RunHarmony(o, group.by.vars = "SampleID", verbose = FALSE)
  o <- FindNeighbors(o, reduction = "harmony", dims = 1:30, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.4, verbose = FALSE)   # 0.4: fewer, robust states for ~3k nuclei
  o <- RunUMAP(o, reduction = "harmony", dims = 1:30, verbose = FALSE)
  cat("Microglia subset re-embedded:", ncol(o), "nuclei\n")
  .emb <- as.data.frame(Embeddings(o, "umap"))[, 1:2]; colnames(.emb) <- c("UMAP_1", "UMAP_2")
  .cc <- list(
    df = data.frame(barcode = colnames(o), .emb,
                    o@meta.data[, c("Condition", "Region", "seurat_clusters", "myeloid_class")],
                    stringsAsFactors = FALSE, row.names = NULL),
    provenance = list(built = Sys.time(), script = "2A_micro_UMAP_FH.R",
                      source = if (.use_immc) basename(IMMC) else basename(ATLAS),
                      atlas_mtime = file.mtime(ATLAS), classes_mtime = file.mtime(.cls_f),
                      params = list(npcs = 30, harmony = "SampleID", dims = "1:30", resolution = 0.4, seed = 42),
                      pkgs = sapply(c("Seurat","SeuratObject","harmony","uwot"), function(p) as.character(packageVersion(p)))))
  .tmp <- paste0(UCACHE, ".tmp"); saveRDS(.cc, .tmp); file.rename(.tmp, UCACHE)
  cat(sprintf("wrote embedding cache: %s\n", UCACHE))
}
df <- .cc$df; rownames(df) <- df$barcode
.n_pvm <- sum(df$myeloid_class == "PVM")
cat(sprintf("POOLED myeloid nuclei drawn: %d (microglia %d + PVM %d = %.1f%% PVM)\n",
            nrow(df), nrow(df) - .n_pvm, .n_pvm, 100 * .n_pvm / nrow(df)))

CLUST_COL <- "seurat_clusters"
df <- df[, c("UMAP_1", "UMAP_2", "Condition", "Region", CLUST_COL)]   # same frame the pre-cache script built
df$Condition <- factor(df$Condition, levels = c("CON", "NHD"))
df$Region    <- factor(df$Region, levels = intersect(REGIONS, unique(df$Region)))
df$cluster   <- factor(df[[CLUST_COL]], levels = sort(unique(as.integer(as.character(df[[CLUST_COL]])))))
cat("microglia UMAP: ", nrow(df), " nuclei; clusters: ",
    paste(levels(df$cluster), collapse = ","), "\n", sep = "")
print(table(df$Condition)); print(table(df$Region))

CLUST_PAL <- setNames(
  c("#4E79A7","#F28E2B","#59A14F","#E15759","#B07AA1","#9C755F","#EDC948",
    "#76B7B2","#FF9DA7","#BAB0AC")[seq_along(levels(df$cluster))],
  levels(df$cluster))

umap_theme <- function(base = 8) {
  theme_pub(base_size = base) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = base - 1, hjust = 0.02),
          legend.key.height = unit(0.8, "lines"),
          legend.text = element_text(size = base - 1),
          legend.title = element_text(size = base - 0.5))
}
umap_theme_arrow <- function(base = 8) {
  umap_theme(base) +
    theme(axis.line   = element_blank(),
          axis.title  = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text   = element_blank(), axis.ticks = element_blank())
}
axis_arrows <- function(d, x = "UMAP_1", y = "UMAP_2",
                        frac = 0.16, lab = c("UMAP 1", "UMAP 2")) {
  xr <- range(d[[x]]); yr <- range(d[[y]]); rng <- max(diff(xr), diff(yr))
  x0 <- xr[1] - 0.01 * diff(xr); y0 <- yr[1] - 0.01 * diff(yr)
  a  <- grid::arrow(length = grid::unit(0.055, "in"), type = "closed")
  list(
    annotate("segment", x = x0, xend = x0 + frac*rng, y = y0, yend = y0,
             arrow = a, linewidth = 0.4),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + frac*rng,
             arrow = a, linewidth = 0.4),
    annotate("text", x = x0, y = y0 - 0.03*rng, hjust = 0, vjust = 1,
             label = lab[1], size = 2.5),
    annotate("text", x = x0 - 0.03*rng, y = y0, angle = 90, hjust = 0, vjust = 0,
             label = lab[2], size = 2.5))
}
PT  <- 0.55
ALP <- 0.85
cents <- df %>% group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

save_panel <- function(p, base, w, h) {
  ggsave(file.path(PANEL, paste0(base, ".png")), p, width = w, height = h, dpi = 600)
  ggsave(file.path(PANEL, paste0(base, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  cat("wrote", base, "\n")
}

# ---- 1. by Condition: pale-blue CON / pale-red NHD, SHUFFLED draw order ------
set.seed(43); df_sh <- df[sample(nrow(df)), ]
p_cond_base <- ggplot(df_sh, aes(UMAP_1, UMAP_2, fill = Condition)) +
  geom_point(size = 0.75, alpha = 0.80, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_COND_PALE,
    breaks = c("CON", "NHD"), labels = c("CON", "NHD"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1))) +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.04)) +
  scale_y_continuous(expand = expansion(mult = 0.04)) +
  umap_theme_arrow() +
  theme(legend.position   = c(0.82, 0.92),
        legend.background  = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key         = element_blank())
p_cond <- p_cond_base + axis_arrows(df)
save_panel(p_cond, "2_micro_UMAP_condition", 3.5, 3.4)

# ---- 2. by cluster (states) with centroid numbers --------------------------
p_clust <- ggplot(df, aes(UMAP_1, UMAP_2, fill = cluster)) +
  geom_point(size = PT, alpha = ALP, stroke = 0.04, shape = 21, colour = "grey20") +
  geom_text(data = cents, aes(x = UMAP_1, y = UMAP_2, label = cluster),
            color = "black", inherit.aes = FALSE,
            size = 3, fontface = "plain", show.legend = FALSE) +
  scale_fill_manual(values = CLUST_PAL, name = "cluster") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") + coord_equal() + umap_theme()
save_panel(p_clust, "2_micro_UMAP_clusters", 3.8, 3.2)

# ---- 3. split CON | NHD, cells coloured by cluster -------------------------
p_split <- ggplot(df, aes(UMAP_1, UMAP_2, fill = cluster)) +
  geom_point(size = PT * 0.9, alpha = ALP, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = CLUST_PAL, name = "cluster") +
  facet_wrap(~ Condition, ncol = 2) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") + coord_equal() +
  umap_theme() +
  theme(strip.text = element_text(size = 9, face = "plain"))
save_panel(p_split, "2_micro_UMAP_condition_split", 6.2, 3.2)

# ---- 4. by Region (shuffled; smaller/fainter so intermixing reads honestly) -
dsh2 <- df[sample(nrow(df)), ]
p_reg_base <- ggplot(dsh2, aes(UMAP_1, UMAP_2, fill = Region)) +
  geom_point(size = 0.60, alpha = 0.68, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_REGION_UMAP, labels = REGION_FULL, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1))) +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0.04)) +
  scale_y_continuous(expand = expansion(mult = 0.04)) +
  umap_theme_arrow() +
  theme(legend.position   = c(0.80, 0.90),
        legend.background  = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key         = element_blank())
p_reg <- p_reg_base + axis_arrows(df)
save_panel(p_reg, "2_micro_UMAP_region", 3.5, 3.4)

# ---- 5. PANEL a composite: Condition | Region (shared arrow-axis style) -----
suppressPackageStartupMessages(library(patchwork))
# Reverted to VERTICAL (stacked)
# assembled Fig2.pdf already uses. Parenthesised because `&` binds tighter than `/` in R.
panel_a <- ((p_cond_base + ggtitle("by condition")) /
            (p_reg + ggtitle("by region"))) &
  theme(plot.title = element_text(size = 8.5, hjust = 0.5, face = "plain",
                                  margin = margin(b = 1)))
panel_a <- panel_a + plot_annotation(
  # Title = "Micro-PVM", the CT_DISPLAY term already printed on Fig 1 (compact UMAP,
  # both concordance panels), the TF heatmaps and SuppFig 1, and the cell_type key in
  # every MAST table -- so Fig 2's lead UMAP does not drift from the rest of the paper.
  # "Micro-PVM" is a joint-class label (the convention when the two are not resolvable);
  # it names a class, it does not claim any individual cell is perivascular. Because the
  # abbreviation does expand to "perivascular macrophage", it must be defined once in the
  # legend + Methods -- see MANIFEST_Fig2.md "Pooled wording".
  title = "Micro-PVM",
  theme = theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain",
                                          colour = "black", family = "",
                                          margin = margin(t = 2, b = 4))))
save_panel(panel_a, "F2a_micro_UMAP_composite", 2.6, 4.35)

cat("=== DONE ===\n")
