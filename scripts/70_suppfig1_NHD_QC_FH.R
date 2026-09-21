# =============================================================================
# 70_suppfig1_NHD_QC_FH.R — Supplementary Figure 1 — NHD single-nucleus atlas QC + annotation frontal + hippocampus rebuild (2 regions, 7 libraries, 55,354 nuclei)
# Faithful port of manuscript_7fig/scripts/70_suppfig1_NHD_QC.R to the NEW
# NHD_FH_harmony.rds atlas. The only substantive differences vs the original
# are data-scoped (region-scoped + cache field names), not stylistic:
#   * atlas  = atlas/NHD_FH_harmony.rds
#   * regions = c("Frontal","Hippo")  (OCC dropped — 2 regions, not 3)
#   * 7 libraries; 55,354 nuclei
#   * theme = scripts/22_publication_theme_FH.R (FH theme, OCC dropped)
#   * panels e/f re-plot from the FH Fig-4 caches (_cache_neuron_ctx_umap_FH.rds
#     / _cache_hippo_neuron_umap_FH.rds / _cache_hippo_marker_dotplot_FH.rds),
#     which have FH-specific field names (meta$subtype / meta$neuron_subtype /
#     $dp). Hippo taxonomy = 6 de-novo subtypes (CA3 / Subiculum / Inh-CGE VIP /
#     Inh-CGE LAMP5 / Inh-MGE / Cajal-Retzius). note: EC-like did not resolve in
#     the FH atlas, so it is absent everywhere (6 subtypes, not the original 7).
#     Cortical Azimuth subclass UMAP is FRONTAL-ONLY (no Occipital in this rebuild).
#
#   a  QC VIOLINS  — nCount_RNA (UMIs), nFeature_RNA (Genes), percent.mt
#      (% mito); one violin per library (7 now), CON vs NHD grouped (dashed
#      divider), white boxplot inset. Three stacked rows, per-metric 99th-pctile
#      display caps. Reflects the nFeature>200 QC floor.
#   b  COVARIATE UMAPs — harmony UMAP small-multiples: Cell type, Condition,
#      Region (Frontal/Hippo — 2 now), Library (SampleID, 7).
#   c  ANNOTATED UMAP — hero atlas UMAP coloured by cell type (clean legend).
#   d  Azimuth MAPPING-SCORE per subclass (violin/box, per-subclass n on top,
#      dashed 0.5 line). predicted.subclass.score is present in the FH atlas so
#      panel d is always the mapping-score panel (the marker-dotplot fallback of
#      the original is not built — it would only fire if the score were absent).
#   e  neuronal-subtype UMAPs — 15 cortical subclasses beside
#      the 6 de-novo hippocampal subtypes. re-plotted natively here (theme_pub
#      base 7, b/c-style points) from the FH Fig-4 caches (never the atlas, never
#      re-running Fig 4).
#   f  hippo subtype marker-validation dotplot — block-diagonal, 6 hippo subtypes
#      as ROWS x canonical private markers as COLUMNS (grouped by subtype, grey
#      strips). size = % expressing, fill = scaled mean expr.
#
# Cache: data/_cache_suppfig1_qc_FH.rds — rebuilt automatically when the atlas
# .rds is newer than the cache (guards the stale-cache trap). The atlas (~2.7 GB)
# is loaded once into the cache; re-styling reuses it.
#
# This script writes panels only (a-f, .png + .pdf).
#
# Run: OMP_NUM_THREADS=2 KMP_DUPLICATE_LIB_OK=TRUE \
#      Rscript scripts/70_suppfig1_NHD_QC_FH.R
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

# apparent-type parity. base_size per panel is set so declared x placement-
# scale lands at ~5.1 pt on the assembled Fig1_suppl page, the same band as the main
# figures. Measured scales: a .579  b .401  c .579  d .651  e .661  f .906.
# Canvas sizes are unchanged so re-linking in Illustrator stays a no-op.
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(tibble); library(patchwork); library(scales); library(forcats)
  library(grid); library(ggh4x)
})
set.seed(42)

# --- Portable project root (never a single hardcoded home) -------------------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ (rebuild) root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_neuron_subclass_FH.R"))   # CTX_*_SUBCLASSES, CTX_SUBCLASS_PAL, subclass_of()

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
DATA  <- file.path(PROJ, "data")
CACHE <- file.path(DATA, "_cache_suppfig1_qc_FH.rds")
OUT   <- file.path(PROJ, "figures", "Supplementary", "SuppFig1_NHD_QC")
PANEL <- file.path(OUT, "panels")
LOGD  <- file.path(PROJ, "logs")
for (d in c(DATA, OUT, PANEL, LOGD))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Canonical region order for this rebuild (Frontal first, then Hippocampus).
# REGION_ORDER / REGION_FULL come from 22_publication_theme_FH.R.
REGION_LVLS <- REGION_ORDER            # c("Frontal","Hippo")

cat("== Supplementary Figure 1 (NHD FH atlas QC) ==\n")
cat("PROJ :", PROJ, "\n")
cat("ATLAS:", ATLAS, "\n")

# Candidate column names for each covariate (first match wins). The FH atlas has
# canonical names, but keep the resolver so the script stays portable.
CANDS <- list(
  sample_id = c("SampleID","Sample","orig.ident","library","Library"),
  subclass  = c("predicted.subclass","predicted_subclass","predicted.id",
                "subclass","Subclass"),
  score     = c("predicted.subclass.score","predicted.id.score",
                "predicted.subclass_score","mapping.score","mapping_score",
                "prediction.score.max","predicted.score"),
  mito      = c("percent.mt","percent_mito","percent.mito","pct_mito"))

pick_col <- function(md, cand) {
  hit <- cand[cand %in% names(md)]
  if (length(hit)) return(hit[1])
  lc  <- tolower(cand); nm <- tolower(names(md))
  hit <- names(md)[match(lc, nm)]; hit <- hit[!is.na(hit)]
  if (length(hit)) hit[1] else NA_character_
}

# ============================================================================
# Build / reuse cache (rebuild when atlas is newer than cache)
# ============================================================================
stale <- !file.exists(CACHE) ||
         (file.exists(ATLAS) && file.mtime(ATLAS) > file.mtime(CACHE))
if (stale) {
  cat("== building cache (atlas load — heavy, ~2.7 GB) ==\n")
  stopifnot("MISSING atlas — NHD_FH_harmony.rds not found" = file.exists(ATLAS))
  obj <- readRDS(ATLAS)
  cat(sprintf("atlas loaded: %d nuclei x %d genes\n", ncol(obj), nrow(obj)))
  md <- obj@meta.data

  stopifnot("atlas missing new_annotation"       = "new_annotation" %in% names(md),
            "atlas missing Condition"            = "Condition"      %in% names(md),
            "atlas missing Region"               = "Region"         %in% names(md),
            "atlas missing nCount_RNA"           = "nCount_RNA"     %in% names(md),
            "atlas missing nFeature_RNA"         = "nFeature_RNA"   %in% names(md))

  cols <- lapply(CANDS, pick_col, md = md)
  cat("-- covariate column detection --\n")
  for (nm in names(cols))
    cat(sprintf("  %-9s -> %s\n", nm, ifelse(is.na(cols[[nm]]), "ABSENT",
                                             cols[[nm]])))
  stopifnot("no sample/library id column found" = !is.na(cols$sample_id),
            "no mito column found"              = !is.na(cols$mito),
            "no Azimuth subclass column found"  = !is.na(cols$subclass),
            "no Azimuth mapping-score column"   = !is.na(cols$score))

  # predicted.subclass.score is present in the FH atlas → panel d = mapping score.
  score_ok <- is.numeric(md[[cols$score]]) && sum(!is.na(md[[cols$score]])) > 0
  stopifnot("mapping score column not usable (non-numeric / all-NA)" = score_ok)
  cat("panel d = Azimuth mapping score (col ", cols$score, ")\n")

  # Guard: the FH rebuild is Frontal + Hippo only. Fail loud if OCC leaks in.
  reg_seen <- sort(unique(as.character(md$Region)))
  stopifnot("unexpected region levels (OCC should be absent in FH rebuild)" =
              all(reg_seen %in% REGION_LVLS))
  if (any(!REGION_LVLS %in% reg_seen))
    cat("NOTE region(s) with 0 nuclei:",
        paste(setdiff(REGION_LVLS, reg_seen), collapse = ", "), "\n")

  # --- Assemble a tidy per-nucleus meta -------------------------------------
  meta <- tibble(
    cell         = rownames(md),
    sample_id    = as.character(md[[cols$sample_id]]),
    Region       = factor(as.character(md$Region), levels = REGION_LVLS),
    Condition    = factor(as.character(md$Condition), levels = c("CON","NHD")),
    cell_type    = as.character(md$new_annotation),
    nCount_RNA   = as.numeric(md$nCount_RNA),
    nFeature_RNA = as.numeric(md$nFeature_RNA),
    percent.mt   = as.numeric(md[[cols$mito]]),
    subclass     = as.character(md[[cols$subclass]]),
    map_score    = as.numeric(md[[cols$score]]))
  stopifnot("NA region level after factoring" = !anyNA(meta$Region),
            "NA condition level after factoring" = !anyNA(meta$Condition))

  # --- UMAP embedding --------------------------------------------------------
  umap_key <- if ("umap" %in% names(obj@reductions)) "umap" else
              if ("UMAP" %in% names(obj@reductions)) "UMAP" else
              names(obj@reductions)[grepl("umap", names(obj@reductions),
                                          ignore.case = TRUE)][1]
  stopifnot("no UMAP reduction found in harmony atlas" = !is.na(umap_key))
  cat("UMAP reduction:", umap_key, "\n")
  emb <- as.data.frame(Embeddings(obj, reduction = umap_key))[, 1:2]
  colnames(emb) <- c("UMAP_1","UMAP_2")
  emb$cell <- rownames(emb)
  meta <- left_join(meta, emb, by = "cell")
  stopifnot("UMAP join produced NA coords" = !anyNA(meta$UMAP_1))

  qc <- list(n_total     = nrow(meta),
             n_samples   = dplyr::n_distinct(meta$sample_id),
             n_regions   = dplyr::n_distinct(as.character(meta$Region)),
             ct_counts   = sort(table(meta$cell_type), decreasing = TRUE),
             cols        = cols, score_ok = score_ok, umap_key = umap_key,
             nfeat_range = range(meta$nFeature_RNA),
             mito_max    = max(meta$percent.mt, na.rm = TRUE),
             # per-library QC medians (reproducibility record; nFeature floor)
             lib_medians = meta %>% group_by(sample_id, Condition, Region) %>%
               summarise(med_UMI  = median(nCount_RNA),
                         med_gene = median(nFeature_RNA),
                         med_mito = median(percent.mt), n = dplyr::n(),
                         .groups = "drop"))
  saveRDS(list(meta = meta, qc = qc), CACHE)
  rm(obj); gc()
  cat("== cache written:", CACHE, "==\n")
} else {
  cat("== reusing cache", basename(CACHE), "==\n")
}

CC   <- readRDS(CACHE)
meta <- CC$meta
qc   <- CC$qc
cols <- qc$cols

# Sanity guards on the loaded cache (fail loud if a stale/wrong cache is reused).
stopifnot("cache region levels != Frontal/Hippo" =
            setequal(levels(meta$Region), REGION_LVLS),
          "cache nuclei count unexpected" = qc$n_total == nrow(meta))

cat("\n-- atlas summary --\n")
cat("nuclei     :", format(qc$n_total, big.mark = ","), "\n")
cat("libraries  :", qc$n_samples, "\n")
cat("regions    :", qc$n_regions, "(", paste(levels(meta$Region), collapse=", "), ")\n")
cat("nFeature   :", qc$nfeat_range[1], "-", qc$nfeat_range[2], "\n")
cat("max %mito  :", round(qc$mito_max, 2), "\n")
cat("panel d    : Azimuth mapping score\n\n")
cat("-- per-library QC medians --\n")
print(as.data.frame(qc$lib_medians))
cat("\n")

# --- Sample ordering: Condition, then Region, then name ----------------------
samp_ord <- meta %>% distinct(sample_id, Condition, Region) %>%
  arrange(Condition, Region, sample_id) %>% pull(sample_id)
meta$sample_id <- factor(meta$sample_id, levels = samp_ord)
nCON <- meta %>% distinct(sample_id, Condition) %>%
  summarise(n = sum(Condition == "CON")) %>% pull(n)

# --- Palettes ----------------------------------------------------------------
# Region tones + condition palettes from 22_publication_theme_FH.R. The covariate
# UMAP points use PAL_REGION_UMAP (0.32 blend); banners keep PAL_REGION_PALE.
PAL_COND_PALE <- c(CON = "#82B2D6", NHD = "#D6837A")   # matches Fig-4 panel a
extra_ct <- setdiff(unique(meta$cell_type), names(PAL_CELLTYPE))
# The two QC-transparency classes (Lymphocyte, Micro-PVM_doublet) used to be filled
# from hcl.colors("Dark 3"), which handed Micro-PVM_doublet a teal indistinguishable
# from Astro once the palette went vibrant — the doublet nuclei read as
# astrocytes. They are now pinned to hues the 8-type palette does not use: brown for
# Lymphocyte, near-black for the doublet class (which also reads correctly as
# "technical, not a cell type"). Any further unforeseen class still falls back to the
# HCL ramp.
PAL_CT_FIXED_EXTRA <- c("Lymphocyte"        = "#7B4B2A",   # brown
                        "Micro-PVM_doublet" = "#2F2F2F")   # near-black
extra_fallback <- setdiff(extra_ct, names(PAL_CT_FIXED_EXTRA))
PAL_CT_EXT <- c(PAL_CELLTYPE,
                PAL_CT_FIXED_EXTRA[intersect(names(PAL_CT_FIXED_EXTRA), extra_ct)],
                if (length(extra_fallback))
                  setNames(grDevices::hcl.colors(length(extra_fallback), "Dark 3"),
                           extra_fallback))
ct_lvls <- c(intersect(names(PAL_CELLTYPE), unique(meta$cell_type)), extra_ct)

# display-only relabel of underscored raw new_annotation levels for legends.
CT_DISPLAY <- function(x) {
  m <- c("Neuron_Ex"         = "Excitatory neuron",
         "Neuron_Inh"        = "Inhibitory neuron",
         "Micro-PVM_doublet" = "Micro-PVM (doublet)")
  ifelse(x %in% names(m), unname(m[x]), x)
}

# ============================================================================
# Panel a — QC violins (3 stacked rows: UMIs / Genes / % mito), CON vs NHD
# ============================================================================
ML <- c(nCount_RNA = "UMIs", nFeature_RNA = "Genes", percent.mt = "% mito")
ql <- meta %>% select(sample_id, Condition, all_of(names(ML))) %>%
  pivot_longer(all_of(names(ML)), names_to = "metric", values_to = "value") %>%
  # cap at per-metric 99th pctile so a few extreme nuclei don't crush the violins
  group_by(metric) %>%
  mutate(value = pmin(value, quantile(value, 0.99, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(metric = factor(ML[metric], levels = ML))
grp_lab <- data.frame(
  metric = factor(ML[1], levels = ML),
  x   = c((1 + nCON) / 2, (nCON + 1 + qc$n_samples) / 2),
  lab = c("CON","NHD"))
pA <- ggplot(ql, aes(sample_id, value, fill = Condition)) +
  geom_vline(xintercept = nCON + 0.5, linetype = "22", linewidth = 0.3,
             colour = "grey55") +
  geom_violin(scale = "width", width = 0.78, linewidth = 0.18, colour = "grey30") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.18,
               fill = "white", alpha = 0.9) +
  geom_text(data = grp_lab, aes(x = x, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = -0.35, size = 3.10,
            colour = c(PAL_COND["CON"], PAL_COND["NHD"])) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = PAL_COND_PALE, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.4),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, size = 8.8),
        panel.spacing = unit(3, "pt"), plot.margin = margin(15, 4, 2, 2))

# ============================================================================
# Panel b — covariate UMAPs (small-multiples): cell type, Condition, Region,
#            Library (SampleID). No Sex/Age/Batch columns in the FH atlas.
# ============================================================================
set.seed(42)
n_draw <- min(45000, nrow(meta))
um <- meta[sample(nrow(meta), n_draw), ]
um <- um[sample(nrow(um)), ]                       # shuffle draw order

umap_mini <- function(df, colour, title, legend = TRUE, arrows = FALSE) {
  # Thin dark outline on every dot (shape 21 + grey20 border, data colour on fill).
  p <- ggplot(df, aes(UMAP_1, UMAP_2, fill = .data[[colour]])) +
    geom_point(size = 0.34, alpha = 0.85, stroke = 0.04, shape = 21, colour = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0.13, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.13, 0.02))) +
    labs(title = title, x = NULL, y = NULL, colour = NULL, fill = NULL) +
    theme_pub(base_size = 12.7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(size = 13.7, hjust = 0.5, face = "plain",
                                    margin = margin(b = 1)),
          legend.key.size = unit(0.42, "cm"),
          legend.text = element_text(size = 12.7),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"))
  if (!legend) p <- p + theme(legend.position = "none")
  # House UMAP axis key (axis_arrows(), shared theme). one key per composite, at the
  # reading entry point -- here the lower-left tile of the 2x2. lab_size in mm is the
  # panel's own declared type (12.7 pt / 2.845) and the arrow head is scaled off the
  # Figure-5 value at base 8, so the key carries the same apparent weight there.
  if (arrows) p <- p +
    coord_cartesian(clip = "off") +
    axis_arrows(xr = range(df$UMAP_1), yr = range(df$UMAP_2),
                lab_size = 12.7 / 2.845, head_in = 0.055 * 12.7 / 8,
                linewidth = 0.4 * 12.7 / 8)
  p
}
# shape-21 legend key override (grey20 border + fill from scale), panel-b parity
ov <- list(shape = 21, colour = "grey20", stroke = 0.3, size = 2.8, alpha = 1)

pB_ct <- umap_mini(um, "cell_type", "Cell type") +
  scale_fill_manual(values = PAL_CT_EXT, breaks = ct_lvls, labels = CT_DISPLAY) +
  guides(fill = guide_legend(ncol = 1, override.aes = ov))
pB_cond <- umap_mini(um, "Condition", "Condition") +
  scale_fill_manual(values = PAL_COND_PALE) +   # pale pair, matches Fig-4 panel a
  guides(fill = guide_legend(override.aes = ov))
pB_reg <- umap_mini(um, "Region", "Region", arrows = TRUE) +   # lower-left tile = the axis key
  scale_fill_manual(values = PAL_REGION_UMAP,   # darker (0.32) tones for points
                    labels = REGION_FULL[levels(um$Region)]) +
  guides(fill = guide_legend(override.aes = ov))
# Library (SampleID) — the FH atlas has no separate Donor column (one lib/donor),
# so the fourth covariate is the 7 libraries. named qualitative palette
# (setNames to levels) — an unnamed or diverging palette silently returns NA at
# its midpoint for odd n and would drop a library from the plot.
um$Library <- factor(as.character(um$sample_id), levels = levels(um$sample_id))
lib_pal <- setNames(grDevices::hcl.colors(nlevels(um$Library), "Dark 3"),
                    levels(um$Library))
stopifnot("library palette has NA colour" = !anyNA(lib_pal))
pB_lib <- umap_mini(um, "Library", "Library") +
  scale_fill_manual(values = lib_pal) +
  guides(fill = guide_legend(ncol = 1 + (nlevels(um$Library) > 8),
                             override.aes = ov))
pB <- wrap_plots(list(pB_ct, pB_cond, pB_reg, pB_lib), nrow = 2)

# ============================================================================
# Panel c — hero annotated UMAP by cell type
# ============================================================================
set.seed(42)
umc <- meta[sample(nrow(meta)), ]                  # full set, shuffled
umc$cell_type <- factor(umc$cell_type, levels = ct_lvls)
pC <- ggplot(umc, aes(UMAP_1, UMAP_2, fill = cell_type)) +
  geom_point(size = 0.36, alpha = 0.9, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_CT_EXT, breaks = ct_lvls, labels = CT_DISPLAY,
                    name = NULL) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  guides(fill = guide_legend(ncol = 1,
           override.aes = list(shape = 21, colour = "grey20", stroke = 0.3,
                               size = 2.2, alpha = 1))) +
  theme_pub(base_size = 8.8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 8.8))

# ============================================================================
# Panel d — Azimuth mapping score per subclass (violin/box), n on top, 0.5 line
# ============================================================================
cat("panel d: Azimuth mapping score per subclass\n")
ds <- meta %>% filter(!is.na(subclass), !is.na(map_score))
n_by <- ds %>% count(subclass)
ord  <- ds %>% group_by(subclass) %>%
  summarise(med = median(map_score), .groups = "drop") %>%
  arrange(med) %>% pull(subclass)          # worst -> best (left -> right)
ds$subclass <- factor(ds$subclass, levels = ord)
n_by$subclass <- factor(n_by$subclass, levels = ord)
sc_pal <- setNames(grDevices::hcl.colors(length(ord), "Dark 3"), ord)
y_top  <- 1.05
pD <- ggplot(ds, aes(subclass, map_score, fill = subclass)) +
  geom_hline(yintercept = 0.5, linetype = "22", linewidth = 0.25,
             colour = "grey60") +
  geom_violin(scale = "width", width = 0.82, linewidth = 0.15,
              colour = "grey35") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.15,
               fill = "white", alpha = 0.9) +
  geom_text(data = n_by, aes(x = subclass, y = y_top, label = comma(n)),
            inherit.aes = FALSE, size = 2.75, colour = "black",
            angle = 45, hjust = 0, vjust = 0) +
  scale_fill_manual(values = sc_pal, guide = "none") +
  scale_y_continuous(limits = c(0, 1.20), breaks = c(0, .25, .5, .75, 1),
                     expand = expansion(mult = c(0.01, 0))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Azimuth mapping score") +
  theme_pub(base_size = 7.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        axis.title.y = element_text(size = 8.8),
        plot.margin = margin(6, 3, 2, 2))

# ============================================================================
# Panels e/f — neuronal-subtype-annotation panels, re-plotted natively in the
# SuppFig1 house style from the FH Fig-4 CACHES only (never the atlas, never
# re-running Fig 4). Content identical to Fig 4; only the theme differs.
#   e  neuronal-subtype UMAPs  (cortical Azimuth subclass beside hippo de-novo)
#   f  hippo subtype marker-validation dotplot (block-diagonal)
# ============================================================================
CTX_CACHE <- file.path(DATA, "_cache_neuron_ctx_umap_FH.rds")
HIP_CACHE <- file.path(DATA, "_cache_hippo_neuron_umap_FH.rds")
DOT_CACHE <- file.path(DATA, "_cache_hippo_marker_dotplot_FH.rds")
stopifnot(
  "MISSING e source — _cache_neuron_ctx_umap_FH.rds (run Fig-4 script 64 first)"
    = file.exists(CTX_CACHE),
  "MISSING e source — _cache_hippo_neuron_umap_FH.rds (run Fig-4 script 64 first)"
    = file.exists(HIP_CACHE),
  "MISSING f source — _cache_hippo_marker_dotplot_FH.rds (run Fig-4 script 4C first)"
    = file.exists(DOT_CACHE))
ctx_obj <- readRDS(CTX_CACHE)
hip_obj <- readRDS(HIP_CACHE)
dot_obj <- readRDS(DOT_CACHE)
# FH cache field names: ctx meta$subtype ; hippo meta$neuron_subtype ; dotplot $dp.
stopifnot("ctx cache missing meta$subtype" = "subtype" %in% colnames(ctx_obj$meta),
          "hip cache missing meta$neuron_subtype" =
            "neuron_subtype" %in% colnames(hip_obj$meta),
          "hippo dotplot cache missing $dp" = !is.null(dot_obj$dp))

# --- taxonomy orders + palettes: single source shared with 64_fig4_neuron_umap_FH ---
# ex/INH order + palette come from _neuron_subclass_FH.R (15 names incl.
# L4 IT = #8DD3C7; Sst Chodl moved to #B3B3B3 — it used to collide with L4 IT here).
EX_SUBC  <- CTX_EX_SUBCLASSES
INH_SUBC <- CTX_INH_SUBCLASSES
CTX_PAL  <- CTX_SUBCLASS_PAL
# The cortical cache is written by 64 with labels applied from the map; re-apply here
# from the same map so this panel can never lag behind Fig 5 (fail loud if it would).
.lab70 <- subclass_of(rownames(ctx_obj$meta))
stopifnot("ctx cache nuclei missing from neuron_subtype_map_FH (re-run 60 then 64)" = !anyNA(.lab70),
          "cortical label outside the 15-name taxonomy" = all(.lab70 %in% c(EX_SUBC, INH_SUBC)))
if (!identical(unname(.lab70), unname(ctx_obj$meta$subtype)))
  cat("NOTE: ctx cache labels differed from the map; using the map (re-run 64 to refresh the cache)\n")
ctx_obj$meta$subtype <- unname(.lab70)
# 6 de-novo hippocampal subtypes (EC-like did not resolve in the FH atlas).
# The cluster that carried it failed the
# private-marker gate in 60_neuron_reannotation_FH.R (TP73 / LHX1 / NDNF absent; the
# call rested on RELN alone), is the lowest-complexity cluster in the hippocampus
# (median nFeature 1312) and draws 186 of 190 nuclei from one library. It is reported
# as an unresolved low-complexity cluster in neutral grey, never as a subtype.
UNRESOLVED_LAB  <- "Unresolved"
HIP_LABEL_ORDER <- c("Inh-CGE (VIP)","Inh-CGE (LAMP5)","Inh-MGE",
                     "CA3","Subiculum", UNRESOLVED_LAB)
HIP_LAB_PAL <- c(
  "Inh-CGE (VIP)"="#8E0152","Inh-CGE (LAMP5)"="#C2549D","Inh-MGE"="#8E44AD",
  "CA3"="#08519C","Subiculum"="#41AB5D")
HIP_LAB_PAL[UNRESOLVED_LAB] <- "#9E9E9E"

# radial view-crop (q=0.97) so each cloud fills its box (copied from Fig-4).
crop_lims <- function(df, q = 0.97, pad = 0.02) {
  cx <- median(df$UMAP_1); cy <- median(df$UMAP_2)
  r  <- sqrt((df$UMAP_1 - cx)^2 + (df$UMAP_2 - cy)^2)
  keep <- r <= quantile(r, q, na.rm = TRUE)
  qx <- range(df$UMAP_1[keep]); qy <- range(df$UMAP_2[keep])
  px <- pad * diff(qx); py <- pad * diff(qy)
  list(XLIM = c(qx[1]-px, qx[2]+px), YLIM = c(qy[1]-py, qy[2]+py))
}
# Square the crop so both panel-e subtype UMAPs render at the same box size.
square_lims <- function(lim) {
  xr <- lim$XLIM; yr <- lim$YLIM
  s  <- max(diff(xr), diff(yr)); cx <- mean(xr); cy <- mean(yr)
  list(XLIM = c(cx - s/2, cx + s/2), YLIM = c(cy - s/2, cy + s/2))
}
# One subtype UMAP in SuppFig1 house style. Panel identity carried by the legend
# TITLE (no plot super-title/caption). Unknown levels -> grey.
sub_umap <- function(cobj, grp_col, levels_order, pal, legend_title,
                     legend_ncol = 1, arrows = FALSE) {
  d <- as.data.frame(cobj$embedding)[, 1:2]; colnames(d) <- c("UMAP_1","UMAP_2")
  d$grp <- as.character(cobj$meta[[grp_col]])
  known <- intersect(levels_order, unique(d$grp))
  oth   <- sort(setdiff(unique(d$grp), levels_order))
  lv    <- c(known, oth)
  d$grp <- factor(d$grp, levels = lv)
  known_cols <- if (!is.null(names(pal))) unname(pal[known]) else pal[seq_along(known)]
  cols  <- setNames(c(known_cols,
                      grDevices::grey.colors(length(oth), 0.55, 0.75)), lv)
  stopifnot("subtype palette has NA colour" = !anyNA(cols))
  lim <- square_lims(crop_lims(d))                     # square view -> equal box size
  set.seed(42); d <- d[sample(nrow(d)), ]              # shuffle draw order
  p <- ggplot(d, aes(UMAP_1, UMAP_2, fill = grp)) +
    geom_point(size = 0.42, alpha = 0.85, stroke = 0.04, shape = 21, colour = "grey20") +
    scale_fill_manual(values = cols, name = legend_title, na.translate = FALSE) +
    guides(fill = guide_legend(ncol = legend_ncol,
             override.aes = list(shape = 21, colour = "grey20", stroke = 0.3,
                                 size = 2, alpha = 1))) +
    coord_cartesian(xlim = lim$XLIM, ylim = lim$YLIM, clip = "off") +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(x = NULL, y = NULL) +
    ggh4x::force_panelsizes(rows = unit(1.85, "in"), cols = unit(1.85, "in")) +
    theme_pub(base_size = 7.7) +
    theme(
          axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.title = element_text(size = 7.7),
          legend.text  = element_text(size = 7.7),
          legend.key.size = unit(0.28, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"),
          plot.margin = margin(2, 2, 2, 7))   # left pad: the rotated "UMAP 2" key
  # Same house key as Figure 5's F5f1_neuron_UMAP_subtypes, which this panel mirrors:
  # Stacked composite -> the key goes on the lower panel only. xr/yr are the view
  # limits, not the raw range, because the embedding is cropped by crop_lims().
  if (arrows) p <- p + axis_arrows(xr = lim$XLIM, yr = lim$YLIM,
                                   lab_size = 7.7 / 2.845,
                                   head_in = 0.055 * 7.7 / 8,
                                   linewidth = 0.4 * 7.7 / 8)
  p
}
cat("== re-plotting e (neuron-subtype UMAPs) in house style ==\n")
pE_ctx <- sub_umap(ctx_obj, "subtype", c(EX_SUBC, INH_SUBC), CTX_PAL,
                   "Cortical subclass")
pE_hip <- sub_umap(hip_obj, "neuron_subtype", HIP_LABEL_ORDER, HIP_LAB_PAL,
                   "Hippocampal subtype", arrows = TRUE)   # lower panel = the axis key
cat(sprintf("   e cortical n=%d (%d levels); hippo n=%d (%d levels)\n",
            nrow(ctx_obj$meta), dplyr::n_distinct(ctx_obj$meta$subtype),
            nrow(hip_obj$meta), dplyr::n_distinct(hip_obj$meta$neuron_subtype)))
# Guard the "6 hippo subtypes, no EC-like" invariant.
hip_seen <- unique(as.character(hip_obj$meta$neuron_subtype))
stopifnot("EC-like unexpectedly present in FH hippo cache" =
            !("EC-like" %in% hip_seen))
if (!setequal(hip_seen, HIP_LABEL_ORDER))
  cat("NOTE hippo cache subtypes differ from expected 6:",
      paste(sort(hip_seen), collapse = ", "), "\n")
pE <- patchwork::wrap_plots(
        pE_ctx  + ggtitle("Cortical (Frontal)"),
        pE_hip  + ggtitle("Hippocampus"),
        ncol = 1) +
      patchwork::plot_layout(guides = "keep") &
      theme(plot.title = element_text(size = 8.3, hjust = 0.5, face = "plain",
                                      margin = margin(b = 1)))

# --- panel f: hippo subtype marker-validation dotplot (block-diagonal) --------
# 6 hippocampal subtypes as ROWS x canonical private markers as columns. Marker
# blocks below are the validated FH markers (per the rebuild spec); every gene
# must exist in the dotplot cache ($dp: per gene x subtype avg.exp/pct.exp/scaled).
# The Cajal-Retzius row previously had no Cajal-Retzius marker block, so the
# one label that needed validating was the only one this panel could not validate. The
# block is added here and the row is renamed: the cluster failed the private-marker gate
# in 60_neuron_reannotation_FH.R (TP73 1.1%; the call rested on RELN, which is promiscuous)
# and is now reported as an unresolved low-complexity cluster.
SUBTYPE_ORDER <- c("CA3","Subiculum","Inh-CGE (VIP)","Inh-CGE (LAMP5)",
                   "Inh-MGE", UNRESOLVED_LAB)
SUBTYPE_DISPLAY <- c("CA3"="CA3","Subiculum"="Subiculum","Inh-CGE (VIP)"="Inh-CGE (VIP)",
                     "Inh-CGE (LAMP5)"="Inh-CGE (LAMP5)","Inh-MGE"="Inh-MGE")
SUBTYPE_DISPLAY[UNRESOLVED_LAB] <- "Unresolved"
MARKER_BLOCKS <- list(
  "CA3"             = c("NECAB1","GRIK4"),
  "Subiculum"       = c("TLE4","ROBO1"),
  "Inh-CGE (VIP)"   = c("VIP","CCK","CNR1"),
  "Inh-CGE (LAMP5)" = c("LAMP5","ID2","NDNF"),
  "Inh-MGE"         = c("LHX6","SST","PVALB"),
  "Cajal-Retzius"   = c("TP73","RELN","LHX1"),
  "pan-gate"        = c("SLC17A7","GAD1","GAD2"))
BLOCK_ORDER   <- names(MARKER_BLOCKS)
BLOCK_DISPLAY <- c(
  "CA3"="CA3","Subiculum"="Subiculum","Inh-CGE (VIP)"="Inh-CGE\n(VIP)",
  "Inh-CGE (LAMP5)"="Inh-CGE\n(LAMP5)","Inh-MGE"="Inh-MGE",
  "Cajal-Retzius"="Cajal-\nRetzius","pan-gate"="Ex / Inh")
row_tbl <- do.call(rbind, lapply(BLOCK_ORDER, function(b)
  data.frame(block = b, gene = MARKER_BLOCKS[[b]], stringsAsFactors = FALSE)))
row_tbl$key <- paste(row_tbl$block, row_tbl$gene, sep = "__")
dpf <- dot_obj$dp
# Fail loud if a requested marker is missing from the cache (would drop a column).
req_genes  <- unique(row_tbl$gene)
miss_genes <- setdiff(req_genes, unique(dpf$gene))
stopifnot("panel-f markers absent from dotplot cache" = length(miss_genes) == 0)
# Keep only the 6 hippo subtypes we plot (cache may carry more).
dpf <- dpf %>% filter(subtype %in% SUBTYPE_ORDER)
plot_df <- row_tbl %>%
  left_join(dpf, by = "gene", relationship = "many-to-many") %>%
  filter(!is.na(scaled))
stopifnot("empty f plot frame" = nrow(plot_df) > 0)
# TRANSPOSED layout (house convention, SuppFig3 panel-d): subtypes are ROWS (y),
# marker gene-keys are COLUMNS (x). y factor reversed so the first SUBTYPE_ORDER
# entry (CA3) sits at the TOP; x keys stay in block order.
plot_df$subtype <- factor(plot_df$subtype, levels = rev(SUBTYPE_ORDER))
plot_df$block   <- factor(plot_df$block,   levels = BLOCK_ORDER)
present_keys    <- row_tbl$key[row_tbl$gene %in% dpf$gene]
plot_df$key     <- factor(plot_df$key, levels = present_keys)
gene_of_key     <- setNames(row_tbl$gene, row_tbl$key)
cat(sprintf("== re-plotting f (hippo marker dotplot): %d tiles, scaled [%.2f, %.2f] ==\n",
            nrow(plot_df), min(plot_df$scaled), max(plot_df$scaled)))

GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
strip_f <- strip_nested(
  background_x = elem_list_rect(fill = "grey90", colour = NA),
  text_x = elem_list_text(colour = "black", face = "plain", size = 6,
                          angle = 0, lineheight = 0.82))
# Detection FLOOR: `scaled` is a z across six rows, so a gene detected in ~1% of nuclei
# still renders a confident red in whichever row is highest (TP73 did exactly that). Dots
# below the floor are drawn hollow; colour is reserved for values that clear detection.
DET_FLOOR <- 10
plot_lo <- plot_df[plot_df$pct.exp <  DET_FLOOR, , drop = FALSE]
plot_hi <- plot_df[plot_df$pct.exp >= DET_FLOOR, , drop = FALSE]
cat(sprintf("panel f detection floor %d%%: %d of %d dots hollow\n",
            DET_FLOOR, nrow(plot_lo), nrow(plot_df)))
pF <- ggplot(plot_df, aes(key, subtype)) +
  geom_point(data = plot_lo, aes(size = pct.exp), shape = 21, fill = "white",
             colour = "grey70", stroke = 0.25) +
  geom_point(data = plot_hi, aes(fill = scaled, size = pct.exp), shape = 21,
             colour = "grey35", stroke = 0.25) +
  facet_nested(cols = vars(block), scales = "free_x", space = "free_x",
               strip = strip_f,
               labeller = labeller(block = as_labeller(BLOCK_DISPLAY))) +
  scale_fill_gradient2(low = GSEA_STEEL, mid = "grey95", high = GSEA_ROSE,
                       midpoint = 0, limits = c(-2.5, 2.5), breaks = c(-2, 0, 2),
                       name = "scaled\nmean expr") +
  scale_size_continuous(range = c(0.3, 4.2), limits = c(0, 100),
                        breaks = c(25, 50, 75), name = "% cells\nexpressing") +
  scale_x_discrete(labels = function(k) gene_of_key[k]) +   # italic gene names
  scale_y_discrete(labels = function(x) unname(SUBTYPE_DISPLAY[x])) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(0.28, "cm"),
                                barheight = unit(1.4, "cm")),
         size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 5.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   size = 6.4, face = "italic", colour = "black"),
        axis.text.y = element_text(size = 5.6, colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x = unit(1.5, "pt"),
        strip.placement = "outside", strip.clip = "off",
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 5.6),
        legend.title = element_text(size = 5.6, lineheight = 0.9),
        plot.margin = margin(4, 4, 4, 4))

# ============================================================================
# Save individual panels (a-f; pdf + png). panels only
# composite manually, so no composite figure is written here.
# ============================================================================
save_panel <- function(p, name, w, h) {
  ggsave(file.path(PANEL, paste0(name, ".pdf")), p, width = w, height = h,
         useDingbats = FALSE, limitsize = FALSE)
  ggsave(file.path(PANEL, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, limitsize = FALSE, device = ragg::agg_png)
}
save_panel(pA, "SuppFig1_a_QC_violins",           4.6, 3.6)
save_panel(pB, "SuppFig1_b_covariateUMAP",        8.6, 5.2)
save_panel(pC, "SuppFig1_c_annotatedUMAP",        5.0, 3.6)
save_panel(pD, "SuppFig1_d_panelD",               5.4, 3.2)
save_panel(pE, "SuppFig1_e_neuron_subtype_UMAP",  3.15, 4.35)  # vertical, as Fig 5; sized to the pinned 1.85 in plots
# 40% slimmer (9.60 -> 5.76), 10% smaller overall, then 20% shorter
save_panel(pF, "SuppFig1_f_hippo_marker_validation", 5.18, 2.30)

# ============================================================================
# Manifest
# ============================================================================
man <- c(
  "# Supplementary Figure 1 — NHD single-nucleus atlas QC + annotation (Frontal + Hippocampus rebuild)",
  "",
  sprintf("_Built %s by scripts/70_suppfig1_NHD_QC_FH.R (set.seed(42))._", Sys.Date()),
  "",
  sprintf("Atlas: `atlas/NHD_FH_harmony.rds` — %s nuclei, %d libraries, %d regions (Frontal, Hippocampus).",
          format(qc$n_total, big.mark = ","), qc$n_samples, qc$n_regions),
  "",
  "Faithful port of manuscript_7fig SuppFig1 to the 7-lane / 2-region FH atlas.",
  "House style: no bold, full region names, big harmonized panel letters, legible UMAPs.",
  "",
  "## Panels",
  "| tag | panel | content |",
  "|-----|-------|---------|",
  "| a | QC violins | nCount_RNA (UMIs), nFeature_RNA (Genes), percent.mt (% mito); one violin per library (7), CON vs NHD grouped (dashed divider), white boxplot inset; 3 stacked rows; per-metric 99th-pctile display caps. |",
  "| b | Covariate UMAPs | harmony UMAP small-multiples: Cell type, Condition, Region (Frontal/Hippocampus), Library (SampleID, 7). |",
  "| c | Annotated UMAP | hero atlas UMAP coloured by cell type (new_annotation), full legend. |",
  sprintf("| d | Azimuth mapping score | per-subclass violin+box, n on top, dashed 0.5 line; subclass col = `%s`, score col = `%s`. |",
          cols$subclass, cols$score),
  "| e | Neuronal-subtype UMAPs | two side-by-side UMAPs: cortical subclasses (excitatory by prefrontal label transfer, inhibitory Azimuth; Frontal only) beside 5 de-novo hippocampal subtypes (CA3 / Subiculum / Inh-CGE VIP / Inh-CGE LAMP5 / Inh-MGE) plus one Unresolved low-complexity cluster (grey; failed the private-marker gate, TP73 1.1%); EC-like did not resolve. RE-PLOTTED NATIVELY in house style from FH Fig-4 caches `_cache_neuron_ctx_umap_FH.rds` + `_cache_hippo_neuron_umap_FH.rds`. |",
  "| f | Hippo subtype marker validation | block-diagonal dotplot: 6 hippocampal subtypes as ROWS, validated private markers as COLUMNS grouped by subtype (grey strips on top) + pan-gate SLC17A7/GAD1/GAD2; size = % expressing, fill = scaled mean expr. RE-PLOTTED NATIVELY from FH Fig-4 cache `_cache_hippo_marker_dotplot_FH.rds`. |",
  "",
  "**Panels e/f:** SAME content as the FH Figure-4 neuronal-subtype panels, RE-PLOTTED natively in the SuppFig1 house style from the FH Figure-4 caches — not raster transplants. Figure 4 + scripts 64/4C untouched.",
  "",
  "## Files (panels only — the composite is assembled by hand)",
  "- `panels/SuppFig1_a_QC_violins.{pdf,png}`",
  "- `panels/SuppFig1_b_covariateUMAP.{pdf,png}`",
  "- `panels/SuppFig1_c_annotatedUMAP.{pdf,png}`",
  "- `panels/SuppFig1_d_panelD.{pdf,png}`",
  "- `panels/SuppFig1_e_neuron_subtype_UMAP.{pdf,png}` (native re-plot, FH Fig-4 caches)",
  "- `panels/SuppFig1_f_hippo_marker_validation.{pdf,png}` (native re-plot, FH Fig-4 cache)",
  "",
  "## Reproducibility",
  "- Cache: `data/_cache_suppfig1_qc_FH.rds` (mtime-guarded vs atlas).",
  "- Seed 42; atlas loaded once; re-styling reuses the cache.")
writeLines(man, file.path(OUT, "SuppFig1_panels.md"))

cat("\n=== DONE ===\n")
cat("wrote panels/ x6 (a-f; pdf+png) + SuppFig1_panels.md\n")
cat("NO composite written — the composite is assembled by hand.\n")
cat("\n-- provenance --\n")
writeLines(capture.output(utils::sessionInfo()),
           file.path(LOGD, "70_suppfig1_NHD_QC_FH_sessionInfo.txt"))
print(utils::sessionInfo())
