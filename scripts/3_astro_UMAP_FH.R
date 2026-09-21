#!/usr/bin/env Rscript
# =============================================================================
# 3_astro_UMAP_FH.R — Astrocyte-subset UMAP panels for Figure 3, frontal + HIPPOCAMPUS (2 regions). Port of 60_fig3_astro_umap.R to the FH rebuild.
# -----------------------------------------------------------------------------
# Re-embeds the Astro subset of the NEW FH atlas on its OWN (the atlas UMAP is
# all-cells): SCT -> RunPCA(30) -> RunHarmony(SampleID) -> FindNeighbors ->
# RunUMAP(seed 42).  Mirrors 2A_micro_UMAP_FH.R.
#
# One continuous astrocyte population — no clustering, no invented substates
# (oligo/astro continuum lesson).  The UMAP shows a continuous reactive shift
# (CON/NHD intermix + a graded reactive + metallothionein module score), not
# discrete NHD islands.
#
# Renders:
#   3_astro_umap_condition   CON vs NHD (pale pair, shuffled)
#   3_astro_umap_region      Frontal / Hippocampus (companion)
#   3_astro_umap_reactive    continuous reactive/DAA module score
#   3_astro_umap_MT          continuous metallothionein module (MT2A/MT3/MT1G/MT1E)
#   F3a_astro_UMAP_composite   condition | region stacked (matches 2a_micro sizing)
#
# FH module-gene re-validation vs the astro MAST full-expressed
# tables: reactive genes kept = those detected + moving UP (GFAP/CD44/OSMR/TNC/...);
# MT module = MT2A/MT3/MT1G/MT1E (brain-canonical MT2A/MT3 headline; MT1F/MT1G/MT1E
# corroborate) — all detected + down in NHD in >=1 region.  Genes present/dropped
# are logged.
#
# Cache: data/_cache_astro_umap_FH.rds (mtime-invalidated vs atlas).  seed 42.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(scales)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, PAL_REGION_UMAP, REGION_FULL, REGION_ORDER

PAL_COND_PALE <- c(CON = "#82B2D6", NHD = "#D6837A")

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
DDIR  <- file.path(PROJ, "data")
CACHE <- file.path(DDIR, "_cache_astro_umap_FH.rds")
PANEL <- file.path(PROJ, "figures", "Figure_3", "panels")
LOGD  <- file.path(PROJ, "logs")
for (d in c(DDIR, PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING atlas NHD_FH_harmony.rds" = file.exists(ATLAS))

REGIONS <- REGION_ORDER                     # c("Frontal","Hippo")

# ---- module gene sets (FH-revalidated; intersect with atlas below) ----------
# Reactive = curated pan-reactive/DAA genes that are detected + up in FH astro MAST.
REACTIVE_GENES <- c("GFAP","CD44","OSMR","TNC","SERPINA3","VIM","CHI3L1","EMP1","STAT3","BAG3")
# MT module: brain-canonical MT2A/MT3 headline + MT1G/MT1E/MT1F corroborators (all down in FH).
MT_GENES <- c("MT2A","MT3","MT1G","MT1E","MT1F")

# =============================================================================
# BUILD (or reload) the embedded astro subset + module scores (mtime guard)
# =============================================================================
cache_stale <- !file.exists(CACHE) || file.mtime(CACHE) < file.mtime(ATLAS)
if (cache_stale) {
  cat("== _cache_astro_umap_FH.rds stale/missing -> re-embed astro subset ==\n")
  suppressPackageStartupMessages(library(harmony))
  atl <- readRDS(ATLAS)
  stopifnot("atlas missing new_annotation" = "new_annotation" %in% colnames(atl@meta.data),
            "atlas missing SampleID"        = "SampleID"       %in% colnames(atl@meta.data),
            "atlas missing Condition"       = "Condition"      %in% colnames(atl@meta.data),
            "atlas missing Region"          = "Region"         %in% colnames(atl@meta.data))
  o <- subset(atl, subset = new_annotation == "Astro"); rm(atl); gc()
  n_ast <- ncol(o); cat(sprintf("   Astro nuclei: %d\n", n_ast))
  stopifnot("empty astro subset" = n_ast > 0)
  print(table(Region = o$Region, Condition = o$Condition))

  DefaultAssay(o) <- "SCT"
  stopifnot("SCT scale.data empty — re-run SCTransform/PrepSCT before embedding" =
              nrow(GetAssayData(o, assay = "SCT", layer = "scale.data")) > 0)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- harmony::RunHarmony(o, group.by.vars = "SampleID", verbose = FALSE)
  o <- FindNeighbors(o, reduction = "harmony", dims = 1:30, verbose = FALSE)
  o <- RunUMAP(o, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)
  cat(sprintf("   re-embedded: %d nuclei\n", ncol(o)))

  DefaultAssay(o) <- "RNA"
  if (length(SeuratObject::Layers(o, assay = "RNA")) > 1) o <- SeuratObject::JoinLayers(o)
  o <- NormalizeData(o, assay = "RNA", verbose = FALSE)
  rg <- intersect(REACTIVE_GENES, rownames(o))
  mg <- intersect(MT_GENES,       rownames(o))
  cat(sprintf("   Reactive genes present: %d/%d (%s); dropped: %s\n",
              length(rg), length(REACTIVE_GENES), paste(rg, collapse = ", "),
              paste(setdiff(REACTIVE_GENES, rg), collapse = ", ")))
  cat(sprintf("   MT genes present: %d/%d (%s); dropped: %s\n",
              length(mg), length(MT_GENES), paste(mg, collapse = ", "),
              paste(setdiff(MT_GENES, mg), collapse = ", ")))
  stopifnot("no reactive genes present" = length(rg) > 0, "no MT genes present" = length(mg) > 0)
  o <- AddModuleScore(o, features = list(rg), name = "Reactive_score", seed = 42, assay = "RNA")
  o <- AddModuleScore(o, features = list(mg), name = "MT_score",       seed = 42, assay = "RNA")
  o$Reactive_score <- o$Reactive_score1; o$Reactive_score1 <- NULL
  o$MT_score       <- o$MT_score1;       o$MT_score1       <- NULL
  attr(o, "reactive_genes_used") <- rg; attr(o, "MT_genes_used") <- mg
  tmp <- paste0(CACHE, ".tmp"); saveRDS(o, tmp); file.rename(tmp, CACHE)
  cat(sprintf("   wrote %s\n", basename(CACHE)))
} else {
  cat("== fresh _cache_astro_umap_FH.rds (>= atlas) -> reload ==\n")
  o <- readRDS(CACHE)
}

# =============================================================================
# plotting frame
# =============================================================================
emb <- as.data.frame(Embeddings(o, "umap"))[, 1:2]
colnames(emb) <- c("UMAP_1", "UMAP_2")
df <- cbind(emb, o@meta.data[, c("Condition", "Region", "Reactive_score", "MT_score")])
df$Condition <- factor(df$Condition, levels = c("CON", "NHD"))
df$Region    <- factor(df$Region, levels = intersect(REGIONS, unique(df$Region)))

n_tot <- nrow(df)
cat(sprintf("\nAstro UMAP frame: %d nuclei\n", n_tot))
cat("By condition:\n"); print(table(df$Condition))
cat("By region:\n");    print(table(df$Region))
cat("By region x condition:\n"); print(table(df$Region, df$Condition))
cat(sprintf("Reactive score range [%.3f,%.3f] median CON=%.3f NHD=%.3f\n",
            min(df$Reactive_score), max(df$Reactive_score),
            median(df$Reactive_score[df$Condition=="CON"]), median(df$Reactive_score[df$Condition=="NHD"])))
cat(sprintf("MT score range       [%.3f,%.3f] median CON=%.3f NHD=%.3f\n",
            min(df$MT_score), max(df$MT_score),
            median(df$MT_score[df$Condition=="CON"]), median(df$MT_score[df$Condition=="NHD"])))

# =============================================================================
# house UMAP theme helpers (self-contained, from 2A_micro_UMAP_FH.R)
# =============================================================================
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
    theme(axis.line = element_blank(), axis.title = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank())
}
axis_arrows <- function(d, x = "UMAP_1", y = "UMAP_2", frac = 0.16,
                        lab = c("UMAP 1", "UMAP 2"), xr = NULL, yr = NULL) {
  if (is.null(xr)) xr <- range(d[[x]]); if (is.null(yr)) yr <- range(d[[y]])
  rng <- max(diff(xr), diff(yr))
  x0 <- xr[1] - 0.01 * diff(xr); y0 <- yr[1] - 0.01 * diff(yr)
  a  <- grid::arrow(length = grid::unit(0.055, "in"), type = "closed")
  list(
    annotate("segment", x = x0, xend = x0 + frac*rng, y = y0, yend = y0, arrow = a, linewidth = 0.4),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + frac*rng, arrow = a, linewidth = 0.4),
    annotate("text", x = x0, y = y0 - 0.03*rng, hjust = 0, vjust = 1, label = lab[1], size = 2.5),
    annotate("text", x = x0 - 0.03*rng, y = y0, angle = 90, hjust = 0, vjust = 0, label = lab[2], size = 2.5))
}

# ---- view crop (radial) — trim UMAP outlier tail so the cloud FILLS the box --
# (astro UMAPs fling a few nuclei far off; radial 97% keeps the main cloud tight,
# same limits on all panels; all nuclei retained in the embedding/cache.)
.cx <- median(df$UMAP_1); .cy <- median(df$UMAP_2)
.rad <- sqrt((df$UMAP_1 - .cx)^2 + (df$UMAP_2 - .cy)^2)
.keep <- .rad <= quantile(.rad, 0.97, na.rm = TRUE)
.qx <- range(df$UMAP_1[.keep]); .qy <- range(df$UMAP_2[.keep])
.px <- 0.02 * diff(.qx); .py <- 0.02 * diff(.qy)
XLIM <- c(.qx[1] - .px, .qx[2] + .px); YLIM <- c(.qy[1] - .py, .qy[2] + .py)
n_outside <- sum(df$UMAP_1 < XLIM[1] | df$UMAP_1 > XLIM[2] | df$UMAP_2 < YLIM[1] | df$UMAP_2 > YLIM[2])
cat(sprintf("view crop (radial 97%%): %d/%d nuclei outside the plotted frame (all retained in embedding)\n",
            n_outside, nrow(df)))
df_in <- df[df$UMAP_1 >= XLIM[1] & df$UMAP_1 <= XLIM[2] &
            df$UMAP_2 >= YLIM[1] & df$UMAP_2 <= YLIM[2], ]

save_panel <- function(p, base, w, h) {
  ggsave(file.path(PANEL, paste0(base, ".png")), p, width = w, height = h, dpi = 600)
  ggsave(file.path(PANEL, paste0(base, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  cat("wrote", base, "\n")
}

# ---- 1. by CONDITION (both coloured + shuffled: continuous shift, not islands) --
set.seed(42); df_sh <- df[sample(nrow(df)), ]
p_cond <- ggplot(df_sh, aes(UMAP_1, UMAP_2, fill = Condition)) +
  geom_point(size = 0.605, alpha = 0.72, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_COND_PALE, breaks = c("CON","NHD"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1))) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") +
  umap_theme_arrow() +
  theme(legend.position = c(0.84, 0.93),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank()) +
  axis_arrows(df_in, xr = XLIM, yr = YLIM)
save_panel(p_cond, "3_astro_umap_condition", 2.6, 2.5)

# ---- 2. by REGION (companion) ----------------------------------------------
set.seed(42); df_shr <- df[sample(nrow(df)), ]
p_reg <- ggplot(df_shr, aes(UMAP_1, UMAP_2, fill = Region)) +
  geom_point(size = 0.605, alpha = 0.70, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_REGION_UMAP, labels = REGION_FULL, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1))) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") +
  umap_theme_arrow() +
  theme(legend.position = c(0.80, 0.93),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank()) +
  axis_arrows(df_in, xr = XLIM, yr = YLIM)
save_panel(p_reg, "3_astro_umap_region", 2.6, 2.5)

# ---- 3/4. continuous module-score FeaturePlots (reactive + metallothionein) --
feature_panel <- function(d, scorecol, title, base = "", opt = "D") {
  d2 <- d[order(d[[scorecol]]), ]                       # low first, high on top
  lims <- as.numeric(quantile(d2[[scorecol]], c(0.01, 0.99), na.rm = TRUE))
  cat(sprintf("%s: colour limits capped to 1st-99th pct [%.3f, %.3f] (oob=squish)\n",
              base, lims[1], lims[2]))
  p <- ggplot(d2, aes(UMAP_1, UMAP_2, fill = .data[[scorecol]])) +
    geom_point(size = 0.605, alpha = 0.9, stroke = 0.04, shape = 21, colour = "grey20") +
    scale_fill_viridis_c(option = opt, name = title, limits = lims, oob = scales::squish,
                         guide = guide_colourbar(barwidth = unit(0.28, "cm"),
                                                 barheight = unit(1.5, "cm"))) +
    coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") +
    umap_theme_arrow() +
    theme(legend.position = "right", legend.title = element_text(size = 6.5),
          legend.text = element_text(size = 6)) +
    axis_arrows(df_in, xr = XLIM, yr = YLIM)
  save_panel(p, base, 2.9, 2.5)
}
feature_panel(df, "Reactive_score", "Reactive\n(GFAP/CD44)", "3_astro_umap_reactive", opt = "D")
feature_panel(df, "MT_score",       "Metallothionein\n(MT2A/MT3)", "3_astro_umap_MT", opt = "D")

# =============================================================================
# Primary composite — visual sibling of F2a_micro_UMAP_composite (FH).
# =============================================================================
suppressPackageStartupMessages(library(patchwork))
p_cond_comp <- ggplot(df_sh, aes(UMAP_1, UMAP_2, fill = Condition)) +
  geom_point(size = 0.682, alpha = 0.80, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_COND_PALE, breaks = c("CON","NHD"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1))) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  umap_theme_arrow() +
  theme(legend.position = c(0.82, 0.92),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank())
p_reg_comp <- ggplot(df_shr, aes(UMAP_1, UMAP_2, fill = Region)) +
  geom_point(size = 0.66, alpha = 0.68, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_REGION_UMAP, labels = REGION_FULL, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.4, alpha = 1))) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  umap_theme_arrow() +
  theme(legend.position = c(0.84, 0.90),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank()) +
  axis_arrows(df_in, xr = XLIM, yr = YLIM)
# Reflowed STACKED (/) -> SIDE-BY-SIDE (|). Stacked at 2.60 x 4.35 the
# panel wasted ~40% of its area on white and, being the tallest thing in its row, it forced
# the whole row to 4.35 in -- which is what squeezed the volcano beside it. Side-by-side at
# 2.85 x 2.25 removes the dead space and lets the row sit at the volcano's height.
# Note this deliberately breaks the visual rhyme with F2a_micro_UMAP_composite (still stacked);
# Fig 3's row geometry needs the width more than the rhyme.
# Legends collected under the composite. Floated inside the point cloud they
# kept truncating "Hippocampus" (the panel is only ~1.6 in per UMAP) and they sat on top of
# the data. A shared bottom row also stops the two sub-plots disagreeing on legend position.
# Back to VERTICAL (stacked). The side-by-side reflow saved height but
# it is also what Fig 2's Micro-PVM UMAP uses, so the
# two figures rhyme. Parenthesised because `&` binds tighter than `|`/`/` in R, which would
# otherwise apply the title theme to the lower plot only.
composite <- (((p_cond_comp + ggtitle("by condition")) /
               (p_reg_comp  + ggtitle("by region"))) &
  theme(plot.title = element_text(size = 8.5, hjust = 0.5, face = "plain",
                                  margin = margin(b = 1)))) +
  plot_annotation(title = "Astrocytes",
    theme = theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain",
                                            colour = "black", family = "",
                                            margin = margin(t = 2, b = 4))))
save_panel(composite, "F3a_astro_UMAP_composite", 2.60, 4.35)

# =============================================================================
# provenance
# =============================================================================
prov <- file.path(LOGD, "3_astro_UMAP_FH_provenance.txt")
sink(prov)
cat("3_astro_UMAP_FH.R provenance\n")
cat("run time:", format(Sys.time()), "\n")
cat("PROJ:", PROJ, "\n")
cat("atlas:", ATLAS, " mtime:", format(file.mtime(ATLAS)), "\n")
cat("cache:", CACHE, " mtime:", format(file.mtime(CACHE)), "\n")
cat("n astro total:", n_tot, "\n")
cat("recipe: SCT -> RunPCA(30) -> RunHarmony(SampleID) -> FindNeighbors(1:30) -> RunUMAP(seed 42)\n")
cat("reactive genes used:", paste(intersect(REACTIVE_GENES, rownames(o)), collapse = ", "), "\n")
cat("MT genes used:", paste(intersect(MT_GENES, rownames(o)), collapse = ", "), "\n\n")
cat("by condition:\n"); print(table(df$Condition))
cat("by region:\n");    print(table(df$Region))
cat("\n--- sessionInfo() ---\n"); print(sessionInfo())
sink()
cat("wrote provenance ->", prov, "\n")
cat("=== DONE ===\n", file = stderr())
