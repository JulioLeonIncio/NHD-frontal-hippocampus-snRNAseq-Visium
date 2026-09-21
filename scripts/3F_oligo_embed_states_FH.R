#!/usr/bin/env Rscript
# =============================================================================
# 3F_oligo_embed_states_FH.R — Shared oligo re-embed + continuum states + Monocle 3 pseudotime cache for the Figure-3 oligo half (panels f/h/i/j).
# Frontal + hippocampus rebuild (2 regions).
# -----------------------------------------------------------------------------
# This is the one heavy step. It re-embeds the Oligo subset of the NEW FH atlas
# on its OWN (the atlas UMAP is all-cells), Louvain-clusters it, ASSIGNS states
# data-driven from marker scores, and fits a Monocle 3 trajectory. Everything is
# cached (mtime-guarded vs the atlas) so 3h/3i/3j reload instead of recomputing.
#
# Conceptual REFRAME (project-established, load-bearing — do not revert):
#   * Oligodendrocytes are one continuum, not discrete subtypes. Clusters are
#     treated as POSITIONS along a maturation AXIS (OPC-adjacent -> mature MOL),
#     not named MFOL/DAO/DAM-like substates.
#   * At most two NHD-associated states are called on top of the continuum:
#         - a high-myelin state  (high structural-myelin module; the surviving
#           heavily-myelinating pole)
#         - a stress state       (heat-shock / CRYAB / ferritin module; NHD-up)
#     These are assigned by module score, data-driven, and reported.
#   * The NHD myelin loss is COMPOSITIONAL (a shift in state proportion along the
#     continuum), not a per-cell collapse of myelin genes — the canonical myelin
#     genes read flat-to-mildly-UP per cell in the FH MAST tables. This cache
#     carries the state labels + pseudotime that let the figure panels show the
#     compositional shift; the per-cell myelin flatness is shown in panel 3g.
#
# The Oligo subset here is Oligo only (new_annotation == "Oligo"); OPC is a
# separate population in this atlas (n=3597) and is scored as a continuum anchor
# but is not re-annotated into the oligo object (no cross-population re-labeling).
#
# Outputs (idempotent, mtime-guarded):
#   data/_cache_oligo_umap.rds        embedded + clustered + scored + labeled +
#                                     pseudotime Seurat object (shared by f/h/i/j)
#   data/_cds_oligo_FH.rds            Monocle 3 CDS (shared by 3h/3j)
#   figures/Figure_3/panels/3h_oligo_umap_trajectory.{png,pdf}
#   figures/Figure_3/panels/3h_oligo_umap_condition.{png,pdf}
#   figures/Figure_3/panels/3h_oligo_umap_states.{png,pdf}
#   figures/Figure_3/panels/3h_oligo_umap_pseudotime.{png,pdf}
#   figures/Figure_3/panels/F4g_oligo_density_shift.{png,pdf}
#   tables/oligo_state_composition_FH.csv   state x region x condition fractions
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(tibble); library(scales); library(patchwork)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
DDIR  <- file.path(PROJ, "data")
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_4", "panels")   # oligo -> Figure 4
LOGD  <- file.path(PROJ, "logs")
for (d in c(DDIR, TDIR, PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING atlas NHD_FH_harmony.rds" = file.exists(ATLAS))

CACHE <- file.path(DDIR, "_cache_oligo_umap.rds")
CDS_C <- file.path(DDIR, "_cds_oligo_FH.rds")

PAL_COND_PALE <- c(CON = "#82B2D6", NHD = "#D6837A")

# ---- Continuum state palette (not subtypes) --------------------------------
# Positions along the maturation continuum + the 2 NHD states, in lineage order.
# Deliberately not the OKI7 subtype scheme: only two NHD states are named, on top
# of a graded mature continuum (early -> mature). Shared by 3h/3i.
PAL_STATE <- c(
  "OPC-proximal (early)"    = "#CCBB44",  # mustard — least-mature oligo pole
  "Mature (continuum)"      = "#4477AA",  # blue    — the mature-oligo continuum
  "High-myelin (NHD)"       = "#228833",  # green   — NHD high-myelin state
  "Stress (NHD)"            = "#AA3377"   # purple  — NHD heat-shock/ferritin state
)
STATE_ORDER <- names(PAL_STATE)

# ---- Module gene sets (continuum anchors + the 2 NHD-state modules) ---------
# Intersected with the assay below; genes present/dropped are logged. Chosen from
# the FH oligo MAST tables + canonical literature; dropout/near-absent genes are
# excluded (min-pct floor enforced at scoring in panels; here we score for state
# assignment, using detected genes only).
sigs <- readRDS(file.path(DDIR, "curated_signatures.rds"))
MOD <- list(
  # early / OPC-proximal maturation flux (low in the fully mature pole)
  early      = c("GPR17","BCAS1","TCF7L2","ENPP6","SOX6","DSCAM"),
  # structural myelin (defines the high-myelin pole)
  myelin     = c("PLP1","MOBP","MOG","MAG","CNP","CLDN11","OPALIN","PLLP"),
  # heat-shock / stress / ferritin (NHD stress state)
  stress     = c("CRYAB","HSPA1A","HSPB1","DNAJB1","FTL","FTH1","APOD","CLU")
)

# =============================================================================
# BUILD (or reload) the embedded + clustered + scored + labeled oligo object
# =============================================================================
cache_stale <- !file.exists(CACHE) || file.mtime(CACHE) < file.mtime(ATLAS)
if (cache_stale) {
  cat("== _cache_oligo_umap.rds stale/missing -> re-embed oligo subset ==\n")
  suppressPackageStartupMessages(library(harmony))
  atl <- readRDS(ATLAS)
  stopifnot("atlas missing new_annotation" = "new_annotation" %in% colnames(atl@meta.data),
            "atlas missing SampleID"        = "SampleID"       %in% colnames(atl@meta.data),
            "atlas missing Condition"       = "Condition"      %in% colnames(atl@meta.data),
            "atlas missing Region"          = "Region"         %in% colnames(atl@meta.data))
  o <- subset(atl, subset = new_annotation == "Oligo"); rm(atl); gc()   # release atlas ASAP
  n_ol <- ncol(o); cat(sprintf("   Oligo nuclei: %d\n", n_ol))
  stopifnot("empty oligo subset" = n_ol > 0)
  cat("   Region x Condition:\n"); print(table(o$Region, o$Condition))
  cat("   Hippo NHD lane(s):\n"); print(table(o$SampleID[o$Region=="Hippo" & o$Condition=="NHD"]))

  DefaultAssay(o) <- "SCT"
  stopifnot("SCT scale.data empty — re-run SCTransform/PrepSCT before embedding" =
              nrow(GetAssayData(o, assay = "SCT", layer = "scale.data")) > 0)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- harmony::RunHarmony(o, group.by.vars = "SampleID", verbose = FALSE)
  o <- FindNeighbors(o, reduction = "harmony", dims = 1:30, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.3, random.seed = 42, verbose = FALSE)
  o <- RunUMAP(o, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)
  cat(sprintf("   re-embedded + clustered: %d nuclei, %d clusters (res 0.3)\n",
              ncol(o), nlevels(o$seurat_clusters)))

  # ---- module scores on RNA (normalized) for data-driven state assignment ----
  DefaultAssay(o) <- "RNA"
  if (length(SeuratObject::Layers(o, assay = "RNA")) > 1) o <- SeuratObject::JoinLayers(o)
  o <- NormalizeData(o, assay = "RNA", verbose = FALSE)
  used <- list()
  for (m in names(MOD)) {
    g <- intersect(MOD[[m]], rownames(o))
    cat(sprintf("   module %-7s present %d/%d: %s ; dropped: %s\n",
                m, length(g), length(MOD[[m]]), paste(g, collapse=","),
                paste(setdiff(MOD[[m]], g), collapse=",")))
    stopifnot(length(g) >= 2)
    o <- AddModuleScore(o, features = list(g), name = paste0(m, "_score"),
                        seed = 42, assay = "RNA")
    o[[paste0(m, "_score")]] <- o[[paste0(m, "_score", 1)]]
    o[[paste0(m, "_score", 1)]] <- NULL
    used[[m]] <- g
  }
  attr(o, "module_genes_used") <- used

  # ---- data-driven continuum-position + 2-NHD-state assignment ---------------
  # For each Louvain cluster (a position on the continuum), compute mean module
  # scores + NHD enrichment, then map to one of 4 continuum positions:
  #   OPC-proximal (early)  : highest early score
  #   Mature (continuum)    : the mature bulk (default)
  #   High-myelin (NHD)     : top myelin-score cluster that is NHD-enriched
  #   Stress (NHD)          : top stress-score cluster that is NHD-enriched
  # "NHD-enriched" = NHD nucleus fraction in the cluster exceeds the global NHD
  # fraction. This keeps the two special states data-driven and NHD-associated,
  # matching the reframe (2 NHD states on a continuum, not invented subtypes).
  md <- o@meta.data
  glob_nhd <- mean(md$Condition == "NHD")
  cl_summ <- md %>%
    mutate(cl = as.character(seurat_clusters)) %>%
    group_by(cl) %>%
    summarise(n = dplyr::n(),
              early = mean(early_score), myelin = mean(myelin_score),
              stress = mean(stress_score),
              nhd_frac = mean(Condition == "NHD"), .groups = "drop") %>%
    mutate(nhd_enriched = nhd_frac > glob_nhd)
  cat("\n   per-cluster module means + NHD enrichment (glob NHD frac = ",
      sprintf("%.3f", glob_nhd), "):\n", sep = "")
  print(as.data.frame(cl_summ), row.names = FALSE)

  # z-score module means across clusters so "top" is comparable between modules.
  z <- function(v) if (sd(v) > 0) (v - mean(v)) / sd(v) else v * 0
  cl_summ <- cl_summ %>% mutate(z_early = z(early), z_myelin = z(myelin), z_stress = z(stress))
  # Stress state = the NHD-enriched cluster with the highest z_stress (and z_stress>0).
  stress_cl <- cl_summ %>% filter(nhd_enriched, z_stress > 0) %>%
    slice_max(z_stress, n = 1) %>% pull(cl)
  # High-myelin state = NHD-enriched cluster with highest z_myelin (excluding stress).
  hi_cl <- cl_summ %>% filter(nhd_enriched, z_myelin > 0, !cl %in% stress_cl) %>%
    slice_max(z_myelin, n = 1) %>% pull(cl)
  # Early = cluster with highest z_early (continuum root anchor), excluding the 2 states.
  early_cl <- cl_summ %>% filter(!cl %in% c(stress_cl, hi_cl)) %>%
    slice_max(z_early, n = 1) %>% pull(cl)
  cat(sprintf("\n   state clusters -> early: %s | high-myelin(NHD): %s | stress(NHD): %s | rest: Mature\n",
              paste(early_cl, collapse=","), paste(hi_cl, collapse=","),
              paste(stress_cl, collapse=",")))
  state_map <- setNames(rep("Mature (continuum)", nrow(cl_summ)), cl_summ$cl)
  state_map[early_cl]  <- "OPC-proximal (early)"
  state_map[hi_cl]     <- "High-myelin (NHD)"
  state_map[stress_cl] <- "Stress (NHD)"
  o$oligo_state <- factor(unname(state_map[as.character(o$seurat_clusters)]),
                          levels = STATE_ORDER)
  cat("\n   oligo_state counts:\n"); print(table(o$oligo_state))
  cat("   oligo_state x Condition:\n"); print(table(o$oligo_state, o$Condition))

  tmp <- paste0(CACHE, ".tmp"); saveRDS(o, tmp); file.rename(tmp, CACHE)
  cat(sprintf("   wrote %s\n", basename(CACHE)))
} else {
  cat("== fresh _cache_oligo_umap.rds (>= atlas) -> reload ==\n")
  o <- readRDS(CACHE)
}

# =============================================================================
# Monocle 3 pseudotime (idempotent, mtime-guarded vs the cache)
# =============================================================================
suppressPackageStartupMessages({ library(monocle3); library(SeuratWrappers); library(igraph) })
cds_stale <- !file.exists(CDS_C) || file.mtime(CDS_C) < file.mtime(CACHE)
if (cds_stale) {
  cat("== _cds_oligo_FH.rds stale/missing -> fit Monocle 3 ==\n")
  cds <- as.cell_data_set(o)
  cds <- cluster_cells(cds, reduction_method = "UMAP", verbose = FALSE)
  cds <- learn_graph(cds, use_partition = FALSE,
                     learn_graph_control = list(minimal_branch_len = 5),
                     verbose = FALSE)
  # Root = the OPC-proximal-early cell with the highest early-maturation score
  # (data-driven continuum root; no OPC cells are in this Oligo-only object).
  early_cells <- colnames(o)[o$oligo_state == "OPC-proximal (early)"]
  stopifnot("no early cells for root" = length(early_cells) > 0)
  root <- early_cells[which.max(o$early_score[early_cells])]
  cds <- order_cells(cds, root_cells = root)
  tmp <- paste0(CDS_C, ".tmp"); saveRDS(cds, tmp); file.rename(tmp, CDS_C)
  cat(sprintf("   wrote %s (root=%s)\n", basename(CDS_C), root))
} else {
  cat("== fresh _cds_oligo_FH.rds (>= cache) -> reload ==\n")
  cds <- readRDS(CDS_C)
}
pt <- monocle3::pseudotime(cds); pt[is.infinite(pt)] <- NA_real_
o$pseudotime <- pt[colnames(o)]
cat(sprintf("pseudotime: %.1f%% cells finite, range [%.2f, %.2f]\n",
            100*mean(!is.na(o$pseudotime)), min(o$pseudotime, na.rm=TRUE),
            max(o$pseudotime, na.rm=TRUE)))

# =============================================================================
# plotting frame + shared UMAP helpers (house style; mirror 3_astro_UMAP_FH.R)
# =============================================================================
emb <- as.data.frame(Embeddings(o, "umap"))[, 1:2]; colnames(emb) <- c("UMAP_1","UMAP_2")
df <- cbind(emb, o@meta.data[, c("Condition","Region","oligo_state","pseudotime",
                                 "myelin_score","stress_score","early_score")])
df$Condition <- factor(df$Condition, levels = c("CON","NHD"))
df$Region    <- factor(df$Region, levels = intersect(REGION_ORDER, unique(df$Region)))
df$oligo_state <- factor(df$oligo_state, levels = STATE_ORDER)

umap_theme_arrow <- function(base = 8) {
  theme_pub(base_size = base) +
    # Blank axis.title and its .x/.y children — theme_pub sets axis.title.x/.y
    # explicitly (margins), so blanking only the parent leaves the frame column
    # names (UMAP_1/UMAP_2) leaking as large axis titles. The arrows carry the axes.
    theme(axis.line = element_blank(),
          axis.title = element_blank(), axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = base - 1),
          legend.title = element_text(size = base - 0.5))
}
axis_arrows <- function(xr, yr, frac = 0.16, lab = c("UMAP 1","UMAP 2")) {
  rng <- max(diff(xr), diff(yr)); x0 <- xr[1] - 0.01*diff(xr); y0 <- yr[1] - 0.01*diff(yr)
  a <- grid::arrow(length = grid::unit(0.055, "in"), type = "closed")
  list(annotate("segment", x=x0, xend=x0+frac*rng, y=y0, yend=y0, arrow=a, linewidth=0.4),
       annotate("segment", x=x0, xend=x0, y=y0, yend=y0+frac*rng, arrow=a, linewidth=0.4),
       annotate("text", x=x0, y=y0-0.03*rng, hjust=0, vjust=1, label=lab[1], size=2.5),
       annotate("text", x=x0-0.03*rng, y=y0, angle=90, hjust=0, vjust=0, label=lab[2], size=2.5))
}
# radial 99% view crop (all nuclei retained in the embedding/cache)
.cx <- median(df$UMAP_1); .cy <- median(df$UMAP_2)
.rad <- sqrt((df$UMAP_1-.cx)^2 + (df$UMAP_2-.cy)^2)
.keep <- .rad <= quantile(.rad, 0.99, na.rm = TRUE)
.qx <- range(df$UMAP_1[.keep]); .qy <- range(df$UMAP_2[.keep])
.px <- 0.02*diff(.qx); .py <- 0.02*diff(.qy)
XLIM <- c(.qx[1]-.px, .qx[2]+.px); YLIM <- c(.qy[1]-.py, .qy[2]+.py)
cat(sprintf("view crop (radial 99%%): %d/%d nuclei outside plotted frame (all retained)\n",
            sum(df$UMAP_1<XLIM[1]|df$UMAP_1>XLIM[2]|df$UMAP_2<YLIM[1]|df$UMAP_2>YLIM[2]), nrow(df)))

save_panel <- function(p, base, w, h) {
  ggsave(file.path(PANEL, paste0(base, ".png")), p, width=w, height=h, dpi=600, device = ragg::agg_png)
  ggsave(file.path(PANEL, paste0(base, ".pdf")), p, width=w, height=h, useDingbats=FALSE)
  cat("wrote", base, "\n")
}

# ---- 3h-condition: CON vs NHD (pale pair, shuffled) — abundant => look abundant
set.seed(42); df_sh <- df[sample(nrow(df)), ]
p_cond <- ggplot(df_sh, aes(UMAP_1, UMAP_2, fill = Condition)) +
  geom_point(size = 0.605, alpha = 0.72, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_COND_PALE, breaks = c("CON","NHD"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape=21, colour="grey20", stroke=0.3, size=2.4, alpha=1))) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") + umap_theme_arrow() +
  theme(legend.position = c(0.84, 0.93),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank()) +
  axis_arrows(XLIM, YLIM)
save_panel(p_cond, "3h_oligo_umap_condition", 2.6, 2.5)

# ---- 3h-region: Frontal / Hippocampus (companion) ---------------------------
# Visual sibling of 3_astro_umap_region and the Micro-PVM region UMAP: same point
# size/alpha grammar, same PAL_REGION_UMAP tones, full region names in the key.
set.seed(42); df_shr <- df[sample(nrow(df)), ]
p_reg <- ggplot(df_shr, aes(UMAP_1, UMAP_2, fill = Region)) +
  geom_point(size = 0.605, alpha = 0.70, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_REGION_UMAP, labels = REGION_FULL, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape=21, colour="grey20", stroke=0.3, size=2.4, alpha=1))) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") + umap_theme_arrow() +
  theme(legend.position = c(0.80, 0.93),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank()) +
  axis_arrows(XLIM, YLIM)
save_panel(p_reg, "3h_oligo_umap_region", 2.6, 2.5)

# =============================================================================
# Primary composite — visual sibling of F3a_astro_UMAP_composite and
# F2a_micro_UMAP_composite ("as we did for astro and microglia").
# Stacked condition-over-region at the same 2.60 x 4.35 as its siblings, same point
# grammar, same title sizes, axis arrows on the lower panel only. Deliberately
# identical so the three cell-type figures rhyme across the paper.
# =============================================================================
suppressPackageStartupMessages(library(patchwork))
p_cond_comp <- ggplot(df_sh, aes(UMAP_1, UMAP_2, fill = Condition)) +
  geom_point(size = 0.682, alpha = 0.80, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_COND_PALE, breaks = c("CON","NHD"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape=21, colour="grey20", stroke=0.3, size=2.4, alpha=1))) +
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
  guides(fill = guide_legend(override.aes = list(shape=21, colour="grey20", stroke=0.3, size=2.4, alpha=1))) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  umap_theme_arrow() +
  theme(legend.position = c(0.84, 0.90),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank(),
        # Clipping fix. axis_arrows() places its labels
        # outside the coord limits and relies on clip="off"; with the composite's
        # expand = 0 they then fall into theme_pub()'s 4 pt plot margin and the "UMAP 1"
        # / "UMAP 2" glyphs are cut by the device edge. Widen the bottom/left margin so
        # the labels have somewhere to sit. note: F3a_astro_UMAP_composite and
        # F2a_micro_UMAP_composite carry the identical defect and need the same one-line
        # fix to stay true siblings.
        plot.margin = margin(4, 4, 12, 12)) +
  axis_arrows(XLIM, YLIM)
composite <- (((p_cond_comp + ggtitle("by condition")) /
               (p_reg_comp  + ggtitle("by region"))) &
  theme(plot.title = element_text(size = 8.5, hjust = 0.5, face = "plain",
                                  margin = margin(b = 1)))) +
  plot_annotation(title = "Oligodendrocytes",
    theme = theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain",
                                            colour = "black", family = "",
                                            margin = margin(t = 2, b = 4))))
save_panel(composite, "F4a_oligo_UMAP_composite", 2.60, 4.35)
cat(sprintf("composite: %d nuclei (CON %d / NHD %d; Frontal %d / Hippo %d)\n",
            nrow(df), sum(df$Condition=="CON"), sum(df$Condition=="NHD"),
            sum(df$Region=="Frontal"), sum(df$Region=="Hippo")))

# ---- 3h-states: continuum position + 2 NHD states --------------------------
set.seed(42); df_shs <- df[sample(nrow(df)), ]
p_state <- ggplot(df_shs, aes(UMAP_1, UMAP_2, fill = oligo_state)) +
  geom_point(size = 0.605, alpha = 0.78, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_STATE, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape=21, colour="grey20", stroke=0.3, size=2.4, alpha=1), ncol=1)) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") + umap_theme_arrow() +
  # Legend outside the panel (right), matching the pseudotime UMAP. The inside
  # placement at c(0.72,0.14) sat over the lower-right cloud and its semi-opaque
  # white box ghosted the Mature/Stress points beneath.
  theme(legend.position = "right",
        legend.key = element_blank(), legend.text = element_text(size = 6)) +
  axis_arrows(XLIM, YLIM)
save_panel(p_state, "3h_oligo_umap_states", 3.5, 2.5)

# ---- 3h-pseudotime: continuum axis (viridis) ---------------------------------
df_pt <- df[order(!is.na(df$pseudotime), df$pseudotime), ]   # NA first, low->high on top
p_pt <- ggplot(df_pt, aes(UMAP_1, UMAP_2, fill = pseudotime)) +
  geom_point(size = 0.605, alpha = 0.9, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_viridis_c(option = "viridis", direction = -1, name = "Pseudotime", na.value = "grey88",   # was "turbo" (rainbow); rainbow scales avoided; direction -1 so yellow = OPC-proximal end, matching the mustard state above
                       # ticks = FALSE: no white tick slices; siblings 4c/4d have none
                       guide = guide_colourbar(barwidth = unit(0.28,"cm"), barheight = unit(1.5,"cm"), ticks = FALSE)) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") + umap_theme_arrow() +
  theme(legend.position = "right", legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6)) +
  axis_arrows(XLIM, YLIM)
save_panel(p_pt, "3h_oligo_umap_pseudotime", 2.9, 2.5)

# =============================================================================
# Second composite: states on top, pseudotime below.
# Same stacked grammar as F4a_oligo_UMAP_composite so the two composites in Figure 4
# rhyme with each other and with the astro/micro siblings — sub-titles at 8.5 pt, a
# header at 11 pt, axis arrows on the lower panel only, and the same plot.margin
# clipping fix (axis_arrows draws outside the coord limits and is otherwise cut by the
# device edge).
# Header is "Maturation continuum" rather than "Oligodendrocytes": the sibling composite
# already carries the cell-type name, and repeating it would waste the one line each
# composite gets to say what it shows. This pair is the continuum argument — where a
# nucleus sits (state) and how far along it sits (pseudotime).
# =============================================================================
suppressPackageStartupMessages(library(patchwork))
# Top panel is rebuilt without axis_arrows rather than reusing p_state.
# p_state carries arrows because it also ships standalone; inside a stacked composite the
# arrows belong to the lower panel only — one axis key for the pair, exactly as
# F4a_oligo_UMAP_composite and the astro/micro siblings do it. Two sets read as two
# unrelated plots and duplicate the same key.
p_state_comp <- ggplot(df_shs, aes(UMAP_1, UMAP_2, fill = oligo_state)) +
  geom_point(size = 0.605, alpha = 0.78, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_STATE, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape=21, colour="grey20", stroke=0.3,
                                                 size=2.4, alpha=1), ncol=1)) +
  coord_equal(xlim = XLIM, ylim = YLIM, clip = "off") + umap_theme_arrow() +
  theme(legend.position = "right", legend.key = element_blank(),
        legend.text = element_text(size = 6),
        plot.margin = margin(4, 4, 4, 12)) +
  ggtitle("by state")
p_pt_comp <- p_pt + ggtitle("by pseudotime") +
  theme(plot.margin = margin(4, 4, 12, 12))
# Guides = "collect": the state legend ("OPC-proximal (early)" ...) is far
# wider than the pseudotime colourbar, so with per-panel legends the two UMAPs were drawn
# at different widths — the same embedding rendered as two different-sized clouds, which
# invites the reader to think they are different data. Collecting both guides into one
# shared column on the right gives the two panels identical plotting areas.
composite_ct <- ((p_state_comp / p_pt_comp + plot_layout(guides = "collect")) &
  theme(plot.title = element_text(size = 8.5, hjust = 0.5, face = "plain",
                                  margin = margin(b = 1)))) +
  plot_annotation(title = "Maturation continuum",
    theme = theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain",
                                            colour = "black", family = "",
                                            margin = margin(t = 2, b = 4))))
save_panel(composite_ct, "F4e_oligo_continuum_composite", 3.50, 4.35)

# =============================================================================
# 3h trajectory overlay — Monocle principal graph + direction-of-time arrows,
# cells coloured by continuum state. (Port of 09o, simplified; house style.)
# =============================================================================
umap_cells <- df[, c("UMAP_1","UMAP_2")]; colnames(umap_cells) <- c("UMAP1","UMAP2")
umap_cells$state <- df$oligo_state
node_coords <- t(cds@principal_graph_aux$UMAP$dp_mst)
node_df <- as.data.frame(node_coords); colnames(node_df) <- c("UMAP1","UMAP2"); node_df$node <- rownames(node_df)
# Orient early -> right? convention: time flows left->right; mirror if early on right.
.ex <- median(umap_cells$UMAP1[umap_cells$state == "OPC-proximal (early)"], na.rm=TRUE)
.mx <- median(umap_cells$UMAP1[umap_cells$state == "Mature (continuum)"], na.rm=TRUE)
FLIP <- isTRUE(.ex > .mx)
if (FLIP) { cat(sprintf("Flip UMAP1: early on RIGHT (%.2f>%.2f) -> mirror so early is LEFT\n", .ex, .mx))
  umap_cells$UMAP1 <- -umap_cells$UMAP1; node_df$UMAP1 <- -node_df$UMAP1
  XT <- rev(-XLIM) } else { cat("No flip (early already LEFT)\n"); XT <- XLIM }
pg <- principal_graph(cds)$UMAP
edges <- igraph::as_edgelist(pg)
# Qualify dplyr::rename — monocle3/SeuratWrappers load after dplyr here and
# mask rename() with an S4 generic (dispatch error otherwise).
edge_df <- data.frame(from = edges[,1], to = edges[,2], stringsAsFactors = FALSE) %>%
  left_join(node_df %>% dplyr::rename(x=UMAP1, y=UMAP2), by = c("from"="node")) %>%
  left_join(node_df %>% dplyr::rename(xend=UMAP1, yend=UMAP2), by = c("to"="node"))
# per-node avg pseudotime to orient arrows low->high
closest <- as.data.frame(cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex)
colnames(closest) <- "vertex"; closest$cell <- rownames(closest); closest$pt <- pt[closest$cell]
node_pt <- closest %>% group_by(vertex) %>% summarise(pt = mean(pt, na.rm=TRUE), .groups="drop") %>%
  mutate(node = paste0("Y_", vertex))
edge_df <- edge_df %>%
  left_join(node_pt %>% select(node, pt_from=pt), by = c("from"="node")) %>%
  left_join(node_pt %>% select(node, pt_to=pt),   by = c("to"="node"))
flip_e <- with(edge_df, !is.na(pt_from) & !is.na(pt_to) & pt_from > pt_to)
edge_df[flip_e, c("x","y","xend","yend","pt_from","pt_to")] <-
  edge_df[flip_e, c("xend","yend","x","y","pt_to","pt_from")]

p_traj <- ggplot() +
  # Point size 0.55 to match the sibling oligo UMAPs (states / condition /
  # pseudotime) so all four share point weight and the cell cloud reads equally
  # dense; the trajectory graph + arrows are drawn after, staying on top.
  geom_point(data = umap_cells, aes(UMAP1, UMAP2, colour = state),
             size = 0.55, alpha = 0.6, stroke = 0) +
  geom_segment(data = edge_df, aes(x=x, y=y, xend=xend, yend=yend),
               colour = "grey35", linewidth = 0.35, lineend = "round",
               arrow = grid::arrow(length = unit(0.03,"in"), type = "closed", angle = 20)) +
  scale_colour_manual(values = PAL_STATE, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1), ncol = 1)) +
  coord_equal(xlim = XT, ylim = YLIM, clip = "off") + umap_theme_arrow() +
  # Legend outside the panel (right), matching the states/pseudotime UMAPs. The
  # inside placement at c(0.70,0.15) sat over the lower-right cloud and its
  # semi-opaque white box ghosted the points beneath.
  theme(legend.position = "right",
        legend.key = element_blank(), legend.text = element_text(size = 6)) +
  axis_arrows(XT, YLIM)
save_panel(p_traj, "3h_oligo_umap_trajectory", 3.5, 2.5)

# =============================================================================
# 3h compositional story — NHD-vs-CON density shift along the continuum
# (pseudotime). This is the load-bearing panel: myelin loss = a shift in where
# cells sit on the continuum, not a per-cell gene collapse. Single-lane-Hippo
# honesty: Hippo NHD rests on one lane -> grey strip + labeller caveat.
# =============================================================================
HIPPO_NHD_LANES <- 1L                      # TFHS000535 is the only NHD Hippo lane
low_conf <- if (HIPPO_NHD_LANES < 2) "Hippo" else character(0)
dens_df <- df %>% filter(!is.na(pseudotime))
# Use nrd0 with adjust=1.6 for a smoother curve
# that reads as a population shift, not sampling noise. Annotate the x-axis poles
# ("OPC-proximal -> Mature -> Stress/High-myelin") so the NHD rightward shift reads
# as a redistribution along the continuum (not a per-cell collapse).
pt_rng  <- range(dens_df$pseudotime, na.rm = TRUE)
pole_lo <- pt_rng[1] + 0.02 * diff(pt_rng)
pole_hi <- pt_rng[2] - 0.02 * diff(pt_rng)
p_dens <- ggplot(dens_df, aes(x = pseudotime, fill = Condition, colour = Condition)) +
  geom_density(alpha = 0.5, linewidth = 0.4, bw = "nrd0", adjust = 1.6) +
  facet_wrap(~Region, nrow = 1,
             labeller = region_full_lowconf_labeller(low_conf)) +
  ggh4x::facetted_pos_scales() +
  scale_fill_manual(values = PAL_COND_PALE, name = NULL) +
  scale_colour_manual(values = c(CON="#5A87AB", NHD="#B36055"), guide = "none") +
  # ASCII " to " (not a unicode arrow): the base pdf() device cannot encode U+2192
  # (mbcsToSbcs conversion failure) and these panels save both .png (ragg) and .pdf.
  scale_x_continuous(name = "Pseudotime  (OPC-proximal to Mature to Stress / High-myelin)",
                     breaks = scales::pretty_breaks(4),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(y = "Density") +
  theme_pub(base_size = 8) +
  # Legend right: on top it sat between the facet strips and the
  # curves, pushing the plotting area down and stealing height from a panel that is
  # only 2.3 in tall. On the right it uses width the two facets were not using.
  theme(legend.position = "right", legend.key.size = unit(0.3,"cm"),
        legend.margin = margin(0, 0, 0, 2),
        strip.text = element_text(size = 8), panel.spacing = unit(0.5,"lines"),
        axis.title.x = element_text(size = 6.6))
# grey the Hippo strip where it rests on the single NHD lane (house honesty rule).
if (length(low_conf)) {
  suppressPackageStartupMessages(library(ggh4x))
  p_dens <- p_dens + facet_wrap2(~Region, nrow = 1,
              labeller = region_full_lowconf_labeller(low_conf),
              strip = strip_region_x(intersect(REGION_ORDER, unique(as.character(dens_df$Region))),
                                     low_conf = low_conf, clip = "off"))
}
save_panel(p_dens, "F4g_oligo_density_shift", 3.84, 1.68)   # 20% smaller

# =============================================================================
# state composition table (feeds panel-i / legend numbers)
# =============================================================================
comp <- o@meta.data %>%
  group_by(Region, Condition, oligo_state) %>% summarise(n = dplyr::n(), .groups = "drop") %>%
  group_by(Region, Condition) %>% mutate(frac = n / sum(n)) %>% ungroup() %>%
  mutate(Region = factor(Region, levels = REGION_ORDER),
         oligo_state = factor(oligo_state, levels = STATE_ORDER)) %>%
  arrange(Region, Condition, oligo_state)
write.csv(comp, file.path(TDIR, "oligo_state_composition_FH.csv"), row.names = FALSE)
cat("\nState composition (fraction within region x condition):\n")
print(as.data.frame(comp %>% mutate(frac = round(frac, 3))), row.names = FALSE)

# =============================================================================
# provenance
# =============================================================================
prov <- file.path(LOGD, "3F_oligo_embed_states_FH_provenance.txt")
sink(prov)
cat("3F_oligo_embed_states_FH.R provenance\n")
cat("run time:", format(Sys.time()), "\n"); cat("PROJ:", PROJ, "\n")
cat("atlas:", ATLAS, " mtime:", format(file.mtime(ATLAS)), "\n")
cat("cache:", CACHE, " mtime:", format(file.mtime(CACHE)), "\n")
cat("cds:",   CDS_C, " mtime:", format(file.mtime(CDS_C)), "\n")
cat("n oligo:", ncol(o), "\n")
cat("recipe: SCT->PCA(30)->Harmony(SampleID)->FindClusters(res0.3)->UMAP(seed42); Monocle3 root=early-max\n")
cat("module genes used:\n"); print(attr(o, "module_genes_used"))
cat("Hippo NHD lanes:", HIPPO_NHD_LANES, "-> low_conf:", paste(low_conf, collapse=","), "\n")
cat("state x condition:\n"); print(table(o$oligo_state, o$Condition))
cat("\n--- sessionInfo() ---\n"); print(sessionInfo())
sink()
cat("wrote provenance ->", prov, "\n")
cat("=== DONE ===\n", file = stderr())
