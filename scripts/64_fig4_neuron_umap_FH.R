#!/usr/bin/env Rscript
# =============================================================================
# 64_fig4_neuron_umap_FH.R — Neuron UMAP composites for Figure 4 (two composites), frontal + hippocampus rebuild. Faithful port of manuscript_7fig/scripts/
# 64_fig4_neuron_umap.R to the 2-region FH atlas.
# -----------------------------------------------------------------------------
# Composite A  `F5a_neuron_UMAP_composite`  header "Neurons"; single panel from the full
#   neuron re-embedding: by CONDITION (CON/NHD both-coloured + shuffled, arrows).
#   The former "by region" sub-panel was dropped: it duplicated the
#   region split already carried by `F5f1_neuron_UMAP_subtypes` (separate cortical +
#   hippocampus embeddings). Condition is the distinct view retained here.
# Composite B  `F5f1_neuron_UMAP_subtypes`  header "Neuron subtypes"; 2 stacked panels,
#   each a separate re-embedding so each group resolves its own taxonomy:
#     Cortical (Frontal)  — 15 subclasses from scripts/_neuron_subclass_FH.R
#                           (excitatory = Jorstad 2023 DLPFC transfer incl.
#                           L4 IT; inhibitory = Azimuth). Labels are applied at render
#                           time from data/neuron_subtype_map_FH.rds, so the cached
#                           embedding is label-independent and never re-embedded for a
#                           relabel; the cache is refreshed when the map is newer.
#     Hippocampus         — de-novo 6-subtype taxonomy from 60_neuron_reannotation_FH
#                           (CA3 / Subiculum / Inh-CGE VIP / Inh-CGE
#                            LAMP5 / Inh-MGE; no EC-like — did not resolve on FH).
#
# What changed vs the 3-region reference (all faithful, reported):
#   * OCC dropped everywhere (2-region rebuild). Cortical panel = frontal only.
#   * Hippo sub-embedding + labels are reused from data/_cache_hippo_neuron_umap_FH.rds
#     (built by 60_neuron_reannotation_FH.R) — the labels are already re-derived on
#     this FH clustering, so we do not re-cluster here (single source of truth).
#   * Hippo taxonomy = six subtypes (EC-like did not resolve on the FH atlas), not 7.
#
# Recipe (full + cortical embeddings): SCT -> RunPCA(30) -> RunHarmony(SampleID) ->
#   FindNeighbors(1:30) -> RunUMAP(seed 42).  Radial view-crop q=0.97.  seed 42.
#   Fixed point size PT_NEURON (abundant clouds must look abundant).
#
# Caches (mtime-invalidated vs atlas): data/_cache_neuron_umap_FH.rds (full),
#   data/_cache_neuron_ctx_umap_FH.rds (cortical).  Hippo cache from 60 (reused).
#
# House rules: no bold, no caption-in-panel, full region names, neuron palette
#   Ex #006D77 / Inh #C42E7B (class), teal/magenta subclass families.  seed 42.
# Run: OMP_NUM_THREADS=2 KMP_DUPLICATE_LIB_OK=TRUE Rscript scripts/64_fig4_neuron_umap_FH.R
# -----------------------------------------------------------------------------
# apparent-type parity. Measured from the assembled Fig5.pdf:
#   composite A is placed at scale 0.511 -> its type is raised to hit the figure-wide
#   on-page targets (tick/legend 5.2 pt, panel title 6.4).
#   composite B was placed at scale 0.397 -- far smaller than every other panel (~0.6) --
#   so matching it by font alone would need a 13 pt legend, and its 14-entry cortical
#   legend does not fit a 3.9 in canvas at that size. It is therefore re-exported on a
#   smaller canvas at the same aspect (3.90x5.60 -> 2.60x3.73): dropped into the same
#   Illustrator frame it lands at scale ~0.60 like its siblings, and its ordinary 8-9 pt
#   type then reads at the same ~5.2 pt on the page. Same aspect, so no distortion.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(scales)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")
`%||%` <- function(a, b) if (is.null(a)) b else a

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, PAL_REGION_UMAP, REGION_FULL, REGION_ORDER
source(file.path(PROJ, "scripts", "_neuron_subclass_FH.R"))       # SUBCLASS_MAP, CTX_*_SUBCLASSES, CTX_SUBCLASS_PAL, subclass_of()

ATLAS     <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
DDIR      <- file.path(PROJ, "data")
CACHE     <- file.path(DDIR, "_cache_neuron_umap_FH.rds")       # full neuron (A)
CACHE_CTX <- file.path(DDIR, "_cache_neuron_ctx_umap_FH.rds")   # cortical (B)
CACHE_HIP <- file.path(DDIR, "_cache_hippo_neuron_umap_FH.rds") # hippo (B) — from 60
PANEL     <- file.path(PROJ, "figures", "Figure_5", "panels")
LOGD      <- file.path(PROJ, "logs")
for (d in c(DDIR, PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING atlas NHD_FH_harmony.rds" = file.exists(ATLAS),
          "MISSING hippo cache (run 60_neuron_reannotation_FH.R first)" = file.exists(CACHE_HIP))

REGIONS <- REGION_ORDER                     # c("Frontal","Hippo")

# ---- cortical subclasses (15; Ex depth order then Inh) + palette: single source ----
# taken from _neuron_subclass_FH.R so 64 and 70 can never drift apart.
EX_SUBC  <- CTX_EX_SUBCLASSES
INH_SUBC <- CTX_INH_SUBCLASSES
CTX_PAL  <- CTX_SUBCLASS_PAL

# ---- hippo de-novo subtype legend order + palette (== 60_neuron_reannotation_FH) ----
# six subtypes on the FH atlas (EC-like did not resolve). Interneurons first.
# The cluster that carried it failed the
# private-marker gate in 60_neuron_reannotation_FH.R (TP73 / LHX1 / NDNF absent; the
# call rested on RELN alone), is the lowest-complexity cluster in the hippocampus
# (median nFeature 1312) and draws 186 of 190 nuclei from one library. It is reported
# as an unresolved low-complexity cluster in neutral grey, never as a subtype.
UNRESOLVED_LAB  <- "Unresolved"
HIP_LABEL_ORDER <- c("Inh-CGE (VIP)","Inh-CGE (LAMP5)","Inh-MGE",
                     "CA3","Subiculum", UNRESOLVED_LAB)
HIP_LAB_PAL <- c(
  "Inh-CGE (VIP)"   = "#8E0152",  "Inh-CGE (LAMP5)" = "#C2549D",
  "Inh-MGE"         = "#8E44AD",  "CA3"             = "#08519C",
  "Subiculum"       = "#41AB5D")
HIP_LAB_PAL[UNRESOLVED_LAB] <- "#9E9E9E"

PT_NEURON <- 0.605  # fixed point size (abundant clouds must look abundant)
PT_ALPHA  <- 0.70

# =============================================================================
# BUILD (or reload) full + cortical embeddings.  Atlas read once, only if stale.
# =============================================================================
atlas_mtime <- file.mtime(ATLAS)
stale <- function(cf) !file.exists(cf) || file.mtime(cf) < atlas_mtime
need_full <- stale(CACHE)
need_ctx  <- stale(CACHE_CTX)
cat(sprintf("cache status: full=%s ctx=%s (hip reused from 60)\n",
            ifelse(need_full,"REBUILD","fresh"), ifelse(need_ctx,"REBUILD","fresh")))

embed_subset <- function(o_sub, seed = 42) {
  DefaultAssay(o_sub) <- "SCT"
  stopifnot("SCT scale.data empty — re-run SCTransform/PrepSCT before embedding" =
              nrow(GetAssayData(o_sub, assay = "SCT", layer = "scale.data")) > 0)
  o_sub <- RunPCA(o_sub, npcs = 30, verbose = FALSE)
  o_sub <- harmony::RunHarmony(o_sub, group.by.vars = "SampleID", verbose = FALSE)
  o_sub <- FindNeighbors(o_sub, reduction = "harmony", dims = 1:30, verbose = FALSE)
  o_sub <- RunUMAP(o_sub, reduction = "harmony", dims = 1:30, seed.use = seed, verbose = FALSE)
  Embeddings(o_sub, "umap")[, 1:2]
}
save_cache <- function(obj, cf) { tmp <- paste0(cf, ".tmp"); saveRDS(obj, tmp); file.rename(tmp, cf)
  cat(sprintf("   wrote %s\n", basename(cf))) }
rt_ctx <- NA_real_

if (need_full || need_ctx) {
  cat("== loading atlas (>=1 cache stale) ==\n")
  suppressPackageStartupMessages(library(harmony))
  atl <- readRDS(ATLAS)
  stopifnot("atlas missing required meta" =
              all(c("new_annotation","SampleID","Condition","Region","predicted.subclass")
                  %in% colnames(atl@meta.data)))
  o_neu <- subset(atl, subset = new_annotation %in% c("Neuron_Ex","Neuron_Inh"))
  rm(atl); gc()
  cat(sprintf("   Neuron nuclei: %d\n", ncol(o_neu)))
  stopifnot("empty neuron subset" = ncol(o_neu) > 0)
  print(table(Region = o_neu$Region, Condition = o_neu$Condition))

  if (need_full) {
    cat("== build full neuron embedding ==\n")
    md <- o_neu@meta.data
    emb <- embed_subset(o_neu)
    keep_meta <- md[rownames(emb), c("Condition","Region","new_annotation","predicted.subclass"),
                    drop = FALSE]
    save_cache(list(embedding = emb, meta = keep_meta, atlas_mtime = atlas_mtime), CACHE)
  }

  if (need_ctx) {
    cat("== build cortical (Frontal) sub-embedding ==\n")
    ctx_cells <- colnames(o_neu)[as.character(o_neu$Region) == "Frontal"]
    o_ctx <- subset(o_neu, cells = ctx_cells)
    cat(sprintf("   cortical (Frontal) nuclei: %d\n", ncol(o_ctx)))
    stopifnot("empty cortical subset" = ncol(o_ctx) > 0)
    tt <- system.time(emb_ctx <- embed_subset(o_ctx)); rt_ctx <- unname(tt["elapsed"])
    # The label is not baked from predicted.subclass any more; it is
    # applied from the map below (label_source / map_mtime recorded in the cache).
    meta_ctx <- data.frame(Condition = as.character(o_ctx$Condition),
                           Region    = as.character(o_ctx$Region),
                           azimuth_subclass = as.character(o_ctx$predicted.subclass),
                           subtype   = NA_character_, row.names = colnames(o_ctx),
                           stringsAsFactors = FALSE)
    save_cache(list(embedding = emb_ctx[rownames(meta_ctx), , drop = FALSE],
                    meta = meta_ctx, atlas_mtime = atlas_mtime, runtime = rt_ctx), CACHE_CTX)
  }
  rm(o_neu); gc()
}

cache_obj <- readRDS(CACHE)
ctx_obj   <- readRDS(CACHE_CTX)
hip_obj   <- readRDS(CACHE_HIP)
if (is.na(rt_ctx)) rt_ctx <- ctx_obj$runtime %||% NA_real_
stopifnot("hippo cache missing meta$neuron_subtype" = "neuron_subtype" %in% colnames(hip_obj$meta))

# The embedding is label-independent; only meta$subtype is refreshed. The cache is
# rewritten whenever the map is newer than the stamp it carries, so 70_suppfig1 (which
# reads this cache) sees the same labels without re-embedding anything.
new_lab <- subclass_of(rownames(ctx_obj$meta))
stopifnot("cortical cache nuclei missing from neuron_subtype_map_FH (re-run 60)" = !anyNA(new_lab),
          "cortical label outside the 15-name taxonomy" = all(new_lab %in% c(EX_SUBC, INH_SUBC)))
ctx_stamp <- ctx_obj$map_mtime %||% as.POSIXct(0, origin = "1970-01-01")
if (is.null(ctx_obj$meta$subtype) || !identical(unname(new_lab), unname(ctx_obj$meta$subtype)) ||
    ctx_stamp < NEURON_MAP_MTIME) {
  n_chg <- if (is.null(ctx_obj$meta$subtype)) length(new_lab) else sum(new_lab != ctx_obj$meta$subtype, na.rm = TRUE)
  cat(sprintf("== cortical labels refreshed from map: %d / %d nuclei relabelled ==\n", n_chg, length(new_lab)))
  ctx_obj$meta$subtype <- unname(new_lab)
  ctx_obj$label_source <- CTX_SUBCLASS_SOURCE
  ctx_obj$map_mtime    <- NEURON_MAP_MTIME
  save_cache(ctx_obj, CACHE_CTX)
} else cat("== cortical labels in cache already match the map ==\n")

# =============================================================================
# UMAP theme helpers (self-contained, from FH sibling panels)
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
axis_arrows <- function(x = "UMAP_1", y = "UMAP_2", frac = 0.16,
                        lab = c("UMAP 1","UMAP 2"), xr, yr, lab_size = 2.5) {
  rng <- max(diff(xr), diff(yr))
  x0 <- xr[1] - 0.01 * diff(xr); y0 <- yr[1] - 0.01 * diff(yr)
  a  <- grid::arrow(length = grid::unit(0.055, "in"), type = "closed")
  list(
    annotate("segment", x = x0, xend = x0 + frac*rng, y = y0, yend = y0, arrow = a, linewidth = 0.4),
    annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + frac*rng, arrow = a, linewidth = 0.4),
    annotate("text", x = x0, y = y0 - 0.03*rng, hjust = 0, vjust = 1, label = lab[1], size = lab_size),
    annotate("text", x = x0 - 0.03*rng, y = y0, angle = 90, hjust = 0, vjust = 0, label = lab[2], size = lab_size))
}
crop_lims <- function(df, q = 0.97, pad = 0.02) {
  cx <- median(df$UMAP_1); cy <- median(df$UMAP_2)
  r  <- sqrt((df$UMAP_1 - cx)^2 + (df$UMAP_2 - cy)^2)
  keep <- r <= quantile(r, q, na.rm = TRUE)
  qx <- range(df$UMAP_1[keep]); qy <- range(df$UMAP_2[keep])
  px <- pad * diff(qx); py <- pad * diff(qy)
  list(XLIM = c(qx[1]-px, qx[2]+px), YLIM = c(qy[1]-py, qy[2]+py))
}
square_lims <- function(lim) {
  xr <- lim$XLIM; yr <- lim$YLIM
  s  <- max(diff(xr), diff(yr)); cx <- mean(xr); cy <- mean(yr)
  list(XLIM = c(cx - s/2, cx + s/2), YLIM = c(cy - s/2, cy + s/2))
}
save_panel <- function(p, base, w, h) {
  ggsave(file.path(PANEL, paste0(base, ".png")), p, width = w, height = h, dpi = 600,
         device = ragg::agg_png)
  ggsave(file.path(PANEL, paste0(base, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  cat("wrote", base, "\n")
}

# =============================================================================
# Composite A — `F5a_neuron_UMAP_composite`  by condition + region
# =============================================================================
emb <- as.data.frame(cache_obj$embedding)[, 1:2]; colnames(emb) <- c("UMAP_1","UMAP_2")
df <- cbind(emb, cache_obj$meta[, c("Condition","Region")])
df$Condition <- factor(df$Condition, levels = c("CON","NHD"))
df$Region    <- factor(df$Region, levels = intersect(REGIONS, unique(df$Region)))
n_tot <- nrow(df)
cat(sprintf("\n[A] Neuron UMAP frame: %d nuclei\n", n_tot))
cat("By condition:\n"); print(table(df$Condition))
cat("By region:\n");    print(table(df$Region))

limA <- crop_lims(df); XA <- limA$XLIM; YA <- limA$YLIM
n_out <- sum(df$UMAP_1 < XA[1] | df$UMAP_1 > XA[2] | df$UMAP_2 < YA[1] | df$UMAP_2 > YA[2])
cat(sprintf("[A] view crop (radial 0.97): %d/%d nuclei outside frame (all retained)\n", n_out, n_tot))

suppressPackageStartupMessages(library(patchwork))
set.seed(42); df_shr <- df[sample(nrow(df)), ]   # region shuffle
set.seed(43); df_shd <- df[sample(nrow(df)), ]   # condition shuffle

# UMAP." So this is now a STACKED composite in the exact grammar of
# F3a_astro_UMAP_composite / F4a_oligo_UMAP_composite: 2.60 x 4.35, an 11 pt cell-type
# header over 8.5 pt sub-titles, point size 0.682/0.66 and alpha 0.80/0.68, shuffled
# draw order, legends floated top-right, axis arrows on the lower panel only, and the
# plot.margin clipping fix.
# TOP = neuronal CLASS in the house two-colour palette (Ex #006D77 teal / Inh #C42E7B
# magenta, [[neuron-class-palette-teal-magenta]]) — the Ex/Inh split is this figure's
# vulnerability axis, so it earns the leading panel. BOTTOM = condition.
df$Class <- factor(ifelse(cache_obj$meta$new_annotation == "Neuron_Ex",
                          "Excitatory", "Inhibitory"),
                   levels = c("Excitatory","Inhibitory"))
PAL_CLASS <- c(Excitatory = "#006D77", Inhibitory = "#C42E7B")
cat("By class:\n"); print(table(df$Class))
set.seed(44); df_shc <- df[sample(nrow(df)), ]   # class shuffle

p_class <- ggplot(df_shc, aes(UMAP_1, UMAP_2, fill = Class)) +
  geom_point(size = 0.682, alpha = 0.80, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_CLASS, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.8, alpha = 1))) +
  coord_equal(xlim = XA, ylim = YA, clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  umap_theme_arrow(11.2) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.margin = margin(t = 0, b = 0), legend.box.spacing = unit(0.06, "cm"),
        legend.key.width = unit(0.42, "cm"))

p_cond_comp <- ggplot(df_shd, aes(UMAP_1, UMAP_2, fill = Condition)) +
  geom_point(size = 0.66, alpha = 0.68, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_COND_PALE, breaks = c("CON","NHD"), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.8, alpha = 1))) +
  coord_equal(xlim = XA, ylim = YA, clip = "off") +
  scale_x_continuous(expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  umap_theme_arrow(11.2) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.background = element_blank(), legend.key = element_blank(),
        legend.margin = margin(t = 0, b = 0), legend.box.spacing = unit(0.06, "cm"),
        legend.key.width = unit(0.42, "cm"),
        plot.margin = margin(4, 4, 12, 12)) +
  axis_arrows(xr = XA, yr = YA, lab_size = 3.59)

compositeA <- (((p_class     + ggtitle("by class")) /
                (p_cond_comp + ggtitle("by condition"))) &
  theme(plot.title = element_text(size = 11.7, hjust = 0.5, face = "plain",
                                  margin = margin(b = 1)))) +
  plot_annotation(title = "Neurons",
    theme = theme(plot.title = element_text(size = 12.5, hjust = 0.5, face = "plain",
                  colour = "black", family = "", margin = margin(t = 2, b = 4))))
save_panel(compositeA, "F5a_neuron_UMAP_composite", 2.60, 4.35)

# =============================================================================
# Composite B — `F5f1_neuron_UMAP_subtypes`  cortical (Frontal) + hippo (de-novo)
# =============================================================================
subtype_panel <- function(emb_mat, grp_vec, levels_order, pal,
                          arrows = FALSE, legend_title = NULL, legend_ncol = 1,
                          lab_size = 3.0, key_h = 0.62, leg_txt = 8.7) {
  d <- as.data.frame(emb_mat)[, 1:2]; colnames(d) <- c("UMAP_1","UMAP_2")
  d$grp <- as.character(grp_vec)
  present <- unique(d$grp)
  known   <- intersect(levels_order, present)
  oth     <- sort(setdiff(present, levels_order))
  lv      <- c(known, oth)
  d$grp   <- factor(d$grp, levels = lv)
  known_cols <- if (!is.null(names(pal))) unname(pal[known]) else pal[seq_along(known)]
  cols <- setNames(c(known_cols, grey.colors(length(oth), start = 0.55, end = 0.75)), lv)
  lim <- square_lims(crop_lims(d)); XL <- lim$XLIM; YL <- lim$YLIM
  n_out <- sum(d$UMAP_1 < XL[1] | d$UMAP_1 > XL[2] | d$UMAP_2 < YL[1] | d$UMAP_2 > YL[2])
  cat(sprintf("   panel: n=%d, %d levels, %d/%d outside frame\n", nrow(d), length(lv), n_out, nrow(d)))
  set.seed(42); d_sh <- d[sample(nrow(d)), ]
  p <- ggplot(d_sh, aes(UMAP_1, UMAP_2, fill = grp)) +
    geom_point(size = PT_NEURON, alpha = PT_ALPHA, stroke = 0.04, shape = 21, colour = "grey20") +
    scale_fill_manual(values = cols, name = legend_title, na.translate = FALSE) +
    guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20", stroke = 0.3, size = 1.7, alpha = 1),
                               ncol = legend_ncol, byrow = TRUE)) +
    coord_cartesian(xlim = XL, ylim = YL, clip = "off") +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    umap_theme_arrow() +
    theme(aspect.ratio = 1, legend.position = "right",
          legend.title = element_text(size = 9.0),
          legend.key.height = grid::unit(key_h, "lines"),
          legend.key.width  = grid::unit(0.40, "lines"),
          legend.text = element_text(size = leg_txt), legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = grid::unit(2, "pt"))
  if (arrows) p <- p + axis_arrows(xr = XL, yr = YL, lab_size = lab_size)
  list(plot = p, levels = lv)
}

cat("\n[B] cortical (Frontal) panel:\n")
ctxP <- subtype_panel(ctx_obj$embedding, ctx_obj$meta$subtype, c(EX_SUBC, INH_SUBC),
                      # One axis key at the reading entry point: the leftmost panel in a
                      # side-by-side composite, the lower panel in a stacked one. Now that
                      # this composite is stacked, the arrows move back to the hippocampal
                      # (lower) panel — same rule as 3h_oligo/3_astro/4_neuron composites.
                      # 15 entries (L4 IT added) — key height 0.62 -> 0.54 lines
                      # so the column no longer overprints the composite title.
                      CTX_PAL, arrows = FALSE, key_h = 0.54, leg_txt = 8.4)
cat("[B] hippocampus panel (de-novo 6 subtypes from 60):\n")
hipP <- subtype_panel(hip_obj$embedding, hip_obj$meta$neuron_subtype, HIP_LABEL_ORDER,
                      HIP_LAB_PAL, arrows = TRUE)

n_ctx <- nrow(ctx_obj$meta); n_hip <- nrow(hip_obj$meta)
cat(sprintf("[B] cortical n=%d, levels: %s\n", n_ctx, paste(ctxP$levels, collapse = ", ")))
hip_lab_n <- table(factor(hip_obj$meta$neuron_subtype, levels = HIP_LABEL_ORDER))
cat(sprintf("[B] hippocampus n=%d, %d subtype labels; per-label n:\n", n_hip, length(HIP_LABEL_ORDER)))
print(hip_lab_n)

# side-by-side: the stacked layout left a large empty vertical band
# between the Cortical and Hippocampus embeddings and split the two legends into a
# top block + bottom block with a big blank right-margin gap.
# Placing the two UMAPs as two columns puts each panel's own legend immediately to
# its right (distinct taxonomies -> distinct colour scales, never merged), closing
# both the inter-block gap and the mid-right blank. `guides="keep"` so each panel
# keeps its own legend rather than patchwork collecting them.
# Vertical / stacked. Each panel keeps its own legend immediately to
# its right — the two taxonomies are distinct (14 cortical Azimuth subclasses vs 6 de-novo
# hippocampal), so the colour scales are never merged. Stacking also puts the two
# embeddings on a common horizontal axis, which side by side they did not share.
compositeB <- ((ctxP$plot + ggtitle("Cortical (Frontal)")) /
               (hipP$plot + ggtitle("Hippocampus"))) +
  plot_layout(guides = "keep") &
  theme(plot.title = element_text(size = 10.6, hjust = 0.5, face = "plain", margin = margin(b = 1)))
compositeB <- compositeB +
  plot_annotation(title = "Neuron subtypes",
    theme = theme(plot.title = element_text(size = 11, hjust = 0.5, face = "plain",
                  colour = "black", family = "", margin = margin(t = 2, b = 4))))
# portrait: two UMAPs stacked, each with its own legend column.
save_panel(compositeB, "F5f1_neuron_UMAP_subtypes", 2.60, 3.73)

# =============================================================================
# provenance
# =============================================================================
prov <- file.path(LOGD, "64_fig4_neuron_umap_FH_provenance.txt")
sink(prov)
cat("64_fig4_neuron_umap_FH.R provenance (TWO composites)\n")
cat("run time:", format(Sys.time()), "\n")
cat("PROJ:", PROJ, "\n")
cat("atlas:", ATLAS, " mtime:", format(atlas_mtime), "\n")
cat("caches: full=", CACHE, " ctx=", CACHE_CTX, " hip(from 60)=", CACHE_HIP, "\n")
cat("cortical embedding runtime (s):", sprintf("%.1f", rt_ctx), "\n\n")
cat("A) full neuron n:", n_tot, "\n")
cat("   by condition:\n"); print(table(df$Condition))
cat("   by region:\n");    print(table(df$Region))
cat("B) cortical (Frontal) n:", n_ctx, " — ", CTX_SUBCLASS_SOURCE, "\n")
cat("   levels:", paste(ctxP$levels, collapse = ", "), "\n")
print(table(ctx_obj$meta$subtype))
cat("\nB) hippocampus n:", n_hip, " — de-novo 6-subtype taxonomy (60_neuron_reannotation_FH)\n")
cat("   EC-like NOT resolved on the FH atlas (6 labels, not 7)\n")
cat("   per-label n:\n"); print(hip_lab_n)
cat("\nNOTE: cortical (DFC Ex + Azimuth Inh) & hippocampal (de-novo) taxonomies are DISTINCT,\n")
cat("shown in SEPARATE panels/colour scales — never merged.\n")
cat("\n--- sessionInfo() ---\n"); print(sessionInfo())
sink()
cat("wrote provenance ->", prov, "\n")
cat("=== DONE ===\n", file = stderr())
