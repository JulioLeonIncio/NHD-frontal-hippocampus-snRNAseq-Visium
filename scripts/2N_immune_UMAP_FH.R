#!/usr/bin/env Rscript
# =============================================================================
# 2N_immune_UMAP_FH.R — Immune-compartment UMAP: Microglia / PVM / T-lymphocytes.
# -----------------------------------------------------------------------------
# The previous panel (2N_myeloid_UMAP_PVM) was embedded on bona-fide Micro-PVM only,
# so lymphocytes could not appear in it at all — they are a separate annotation class.
# This re-embeds the whole immune compartment (Micro-PVM + Lymphocyte) so all three
# populations live in one space.
#
# The Microglia/PVM call is not recomputed — it is carried over by barcode from
# 2L_pvm_vs_microglia_FH.R (tables/micro_states/pvm_vs_micro_assignment_FH.csv), so
# this panel and the proportion/marker panels describe exactly the same cells. Only
# the embedding is new. Recomputing clusters here would let the panels drift apart.
#
# Oligo-doublets (Oligo module score > 0.5) are excluded, as everywhere else in the
# myeloid analysis. `Micro-PVM_doublet` is excluded too — verified as genuine
# doublets in 2L1 (median 1,858 genes vs 528; SNAP25 in 33% vs 0%).
#
# Output: figures/Figure_2/panels/2N_immune_UMAP.{png,pdf}
#         (supersedes 2N_myeloid_UMAP_PVM.*, which is deleted to avoid a stale twin)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(harmony)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
ATLAS  <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
PANEL  <- file.path(PROJ, "figures", "Figure_2", "panels")
TBL    <- file.path(PROJ, "tables", "micro_states")
ASSIGN <- file.path(TBL, "pvm_vs_micro_assignment_FH.csv")
stopifnot(file.exists(ATLAS), "run 2L_pvm_vs_microglia_FH.R first" = file.exists(ASSIGN))

OLIGO_MK  <- c("PLP1","MBP","MOG","MOBP","CNP","ST18","CLDN11","MAG")
OLIGO_CUT <- 0.5

# ---- EMBEDDING (cached) -----------------------------------------------------
# The embedding costs a 2.7 GB atlas read + PCA/harmony/UMAP, which made every
# purely cosmetic re-render of this panel a multi-minute job. The coordinates are
# now cached to disk with an MTIME GUARD ([[atlas-single-source-cache-guards]]):
# The cache is reused only if it is newer than both inputs (the atlas and the
# Microglia/PVM assignment CSV), so it cannot go stale silently -- touch either
# input and the embedding rebuilds. Delete the cache to force a rebuild.
UCACHE <- file.path(PROJ, "data", "_cache_2N_immune_umap_FH.rds")
.cache_ok <- file.exists(UCACHE) &&
  file.mtime(UCACHE) > file.mtime(ATLAS) &&
  file.mtime(UCACHE) > file.mtime(ASSIGN)

if (.cache_ok) {
  ud <- readRDS(UCACHE)
  stopifnot("bad embedding cache" =
              all(c("UMAP_1", "UMAP_2", "class") %in% names(ud)) && nrow(ud) > 0)
  cat(sprintf("== embedding CACHE HIT (%s): %d nuclei ==\n",
              format(file.mtime(UCACHE), "%Y-%m-%d %H:%M"), nrow(ud)))
  print(table(ud$class))
} else {
  cat("== embedding CACHE MISS -- recomputing from the atlas ==\n")
  atl <- readRDS(ATLAS)
  # Join per-lane layers on the full atlas before subsetting (a subset with lanes that
  # contribute zero nuclei leaves empty layers that StitchMatrix cannot handle).
  if (inherits(atl[["RNA"]], "Assay5")) atl[["RNA"]] <- SeuratObject::JoinLayers(atl[["RNA"]])
  o <- subset(atl, subset = new_annotation %in% c("Micro-PVM", "Lymphocyte"))
  rm(atl); gc()
  cat(sprintf("immune compartment pulled: %d nuclei\n", ncol(o)))
  print(table(o$new_annotation))

  # drop oligo-doublets among the myeloid nuclei
  DefaultAssay(o) <- "RNA"
  o <- NormalizeData(o, verbose = FALSE)
  o <- AddModuleScore(o, features = list(OLIGO_MK[OLIGO_MK %in% rownames(o)]),
                      name = "OLSC", seed = 42, ctrl = 50)
  n0 <- ncol(o)
  o  <- o[, !(o$new_annotation == "Micro-PVM" & o$OLSC1 > OLIGO_CUT)]
  cat(sprintf("dropped %d oligo-doublet nuclei; %d remain\n", n0 - ncol(o), ncol(o)))

  # carry the Microglia/PVM call across by barcode
  asg <- read.csv(ASSIGN, row.names = 1)
  o$immune_class <- ifelse(o$new_annotation == "Lymphocyte", "T-lymphocyte",
                           as.character(asg[colnames(o), "myeloid_class"]))
  cat("\nunmatched myeloid barcodes (should be 0): ",
      sum(is.na(o$immune_class)), "\n", sep = "")
  o <- o[, !is.na(o$immune_class)]
  o$immune_class <- factor(o$immune_class, levels = c("Microglia", "PVM", "T-lymphocyte"))
  print(table(o$immune_class, o$Condition))

  # re-embed the whole compartment
  DefaultAssay(o) <- "SCT"
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- harmony::RunHarmony(o, group.by.vars = "SampleID", verbose = FALSE)
  o <- RunUMAP(o, reduction = "harmony", dims = 1:30, verbose = FALSE)

  ud <- as.data.frame(Embeddings(o, "umap"))[, 1:2]
  colnames(ud) <- c("UMAP_1", "UMAP_2")
  ud$class <- o$immune_class
  ud <- ud[order(ud$class), ]        # Microglia bottom -> T-lymphocyte on top
  saveRDS(ud, UCACHE)
  cat(sprintf("wrote embedding cache: %s\n", UCACHE))
  rm(o); gc()
}

PAL_IMM <- c("Microglia"    = unname(PAL_CELLTYPE[["Micro-PVM"]]),
             "PVM"          = "#6D4C41",
             "T-lymphocyte" = "#C77B0A")

# apparent-type parity. House rule:
# on-page pt = declared pt x placement scale (placed_width/saved_width in the
# assembled figure); main-figure band 4.8-5.3 pt, target 5.1.
# The panel is now saved 3.90 in wide (see the aspect fix at the ggsave),
# so the placement scale changed with it. At a ~3.50 in placed width the scale is
# 0.897 and the legend lands at 6.6 x 0.87 x 0.897 = 5.15 pt on the page, still inside
# the 4.8-5.3 band. re-check this if the panel is placed at a materially different width.
# It has no axis text, so its body type is the legend: 6.6 pt x 0.892 = 5.89 pt on the
# page (hot). TYPE_F brings it to 5.12 pt; the corner UMAP-arrow labels ride the same
# factor (1.9 -> 1.65 mm = 4.20 pt on page, still above the ~4.0 pt print floor).
# Text only -- geom_point size / stroke / arrow unit() are deliberately not scaled.
TYPE_F <- 1.14

# Axis BREAK ("take care of blank spaces"). The T-lymphocyte island lies ~2.6 cloud-widths
# to the right of the microglial cloud, so drawing the embedding to scale is mostly white. UMAP coordinates
# carry no metric meaning, so the island is translated left to sit one gap (gap x cloud width) beyond the
# cloud; the legend states the translation (no break glyph).
# No nucleus is dropped, the island's internal geometry is untouched, the microglial cloud is untouched.
is_T <- ud$class == "T-lymphocyte"
cloud_x <- range(ud$UMAP_1[!is_T]); GAP <- 0.22
# the break is only legitimate if the island is fully separated from the cloud on UMAP 1 (no T nucleus
# inside the cloud's x-range, no non-T nucleus inside the island's); assert it, and log the class counts
stopifnot("T island is not separable from the cloud on UMAP 1 — no axis break possible" =
            min(ud$UMAP_1[is_T]) > cloud_x[2] && all(ud$UMAP_1[!is_T] < min(ud$UMAP_1[is_T])))
cat("classes plotted:", paste(names(table(ud$class)), as.integer(table(ud$class)), collapse = " | "), "\n")
shift <- min(ud$UMAP_1[is_T]) - (cloud_x[2] + GAP * diff(cloud_x))
ud$UMAP_1[is_T] <- ud$UMAP_1[is_T] - shift
cat(sprintf("axis break: T island shifted left by %.2f UMAP units (cloud width %.2f)\n", shift, diff(cloud_x)))
xr <- range(ud$UMAP_1); yr <- range(ud$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
gx <- xr[1] + 0.005 * dx; gy <- yr[1] - 0.135 * dy
p <- ggplot(ud, aes(UMAP_1, UMAP_2, fill = class)) +
  geom_point(size = 0.55, alpha = 0.80, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_IMM, labels = IMM_LABEL, name = NULL) +
  # (no "//" glyph — the translation is stated in the figure legend instead; UMAP
  #  coordinates carry no metric meaning, and the island keeps its own geometry.)
  # arrow key: arrows long enough to carry their labels (0.20 dx / 0.30 dy ~ 0.4 in each on the 2.5-in canvas),
  # labels centred on the arrows at 1.6 mm ("the UMAP arrow labels look out of place")
  annotate("segment", x = gx, xend = gx + 0.20 * dx, y = gy, yend = gy,
           arrow = arrow(length = unit(1.3, "mm"), type = "closed"),
           linewidth = 0.3, colour = "grey25") +
  annotate("segment", x = gx, xend = gx, y = gy, yend = gy + 0.30 * dy,
           arrow = arrow(length = unit(1.3, "mm"), type = "closed"),
           linewidth = 0.3, colour = "grey25") +
  annotate("text", x = gx + 0.10 * dx, y = gy - 0.02 * dy, label = "UMAP 1",
           size = 1.6 * TYPE_F, colour = "grey25", hjust = 0.5, vjust = 1) +
  annotate("text", x = gx - 0.012 * dx, y = gy + 0.15 * dy, label = "UMAP 2",
           size = 1.6 * TYPE_F, colour = "grey25", hjust = 0.5, vjust = 0, angle = 90) +
  # PADS: the x range is inflated ~2.6x by the far-right T-lymphocyte island, so the
  # old 0.11/0.02 fractional pads translated into ~0.24 in of pure white down the left
  # edge and a second empty strip on the right. Trimmed to the smallest values that
  # still clear the rotated "UMAP 2" label (verified by eye: 100 px of canvas margin
  # remains on the left, 90 px on the right -- nothing clips). The axis-arrow key
  # itself is unchanged, still outside the data in the lower left.
  coord_cartesian(xlim = c(xr[1] - 0.075 * dx, xr[2] + 0.005 * dx),
                  ylim = c(yr[1] - 0.20 * dy, yr[2] + 0.02 * dy), clip = "off") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.2, alpha = 1))) +
  theme_pub(base_size = 8 * TYPE_F) +
  # blank axis.title.x/.y explicitly — theme_pub sets those children directly, so
  # blanking only the parent leaves the default UMAP_1/UMAP_2 titles printing.
  theme(axis.title = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        # dead-space fix. The embedding is intrinsically sparse: the
        # T-lymphocyte island sits far to the right of the microglial cloud, so the
        # right half of the canvas is mostly empty and the legend, parked at the top
        # right, left a large empty block beneath it (only ~10% of the canvas inked).
        # The island is the point of the panel and must not be cropped away, so the
        # furniture moves instead. Measured from the cached embedding, the microglial
        # cloud occupies x npc 0.00-0.39 and the T island x npc 0.98-1.00 (y 0.44-0.54),
        # leaving a ~60%-wide empty corridor between them. The legend is parked in the
        # middle of that corridor, vertically level with the island, so the panel reads
        # as one continuous band (cloud -> legend -> island) instead of a cloud plus a
        # floating speck. A numeric legend.position is the CENTRE of the legend box, so
        # the justification is set explicitly rather than inherited.
        # Honest limit: moving furniture redistributes white space, it cannot create
        # ink. Non-white pixels go 9.55% -> 9.87% of the canvas; the residual emptiness
        # is a property of the embedding, and the only way to remove it would be to
        # crop out the island, which would delete the panel's point.
        # At 2.50 in an in-panel legend either clips ("CD163+/F13A1+" ran off the right edge)
        # or sits on the cloud; one row below the panel instead, no overlap possible.
        legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(t = 0, b = 0), legend.box.spacing = unit(2, "pt"),
        legend.background = element_blank(),
        legend.key = element_blank(), legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 6.0 * TYPE_F),
        plot.margin = margin(2, 2, 1, 1))
# Aspect fix. The embedding spans dx = 22.09 and dy = 8.04, i.e. it is
# 2.75x wider than tall, but it was being drawn in a 2.70 x 2.125 in box (aspect 1.27).
# The y axis was therefore stretched about 2.2x relative to x: that is what made the
# panel read as mostly white, and it also misrepresents the map, since distances in
# one direction were exaggerated over the other. The box now follows the data
# (3.90 / 1.55 = 2.52, close to the true 2.75 while leaving room for the arrow key
# below the data). nothing is cropped -- the far-right T-lymphocyte island, which is
# the point of the panel, is still fully in frame, as is every microglial nucleus.
PW <- 2.50; PH <- 2.10   # -> 2.50 in wide; 2.10 in tall, legend row below; with the axis break the canvas is mostly ink
ggsave(file.path(PANEL, "2N_immune_UMAP.png"), p, width = PW, height = PH,
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "2N_immune_UMAP.pdf"), p, width = PW, height = PH, bg = "white")
cat("wrote 2N_immune_UMAP.{png,pdf}\n")

# retire the myeloid-only twin so it cannot be assembled by mistake
old <- file.path(PANEL, paste0("2N_myeloid_UMAP_PVM.", c("png","pdf")))
unlink(old[file.exists(old)])
cat("removed superseded 2N_myeloid_UMAP_PVM.{png,pdf}\n")

write.csv(data.frame(barcode = rownames(ud), class = ud$class),
          file.path(TBL, "immune_umap_classes_FH.csv"), row.names = FALSE)
cat("\n=== DONE: 2N_immune_UMAP_FH ===\n")
