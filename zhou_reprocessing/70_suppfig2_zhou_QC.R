# =============================================================================
# 70_suppfig2_zhou_QC.R — Supplementary Figure 2 — Zhou et al. 2023 external cohort: QC + annotation
# Sibling of Supp Fig 1 (NHD atlas QC). Same house style (theme_pub, PAL_*),
# but the Zhou cohort is single-region (occipital cortex), so the grouping
# dimension is DIAGNOSIS (CON vs NHD/DAP12-null) instead of Region.
#
# Cohort: Zhou et al. 2023 Nat Immunol. Occipital cortex.
#   n = 3 NHD (DAP12-null: NHD1/NHD3 c.2T>C, NHD2 c.141Gdel) vs 11 CON.
#   14,293 NHD + 49,065 CON = 63,358 nuclei; 8 annotated cell types.
#
# Panels (mirror atlas Supp Fig 1 a-e; no embedded caption banner):
#   a. Per-donor QC violins (UMIs / Genes / % mito), 3 stacked rows over one
#      shared donor x-axis, CON vs NHD grouped (dashed divider, white boxplot
#      inset), fill = PAL_COND (same tokens as Supp Fig 1) — pixel-sibling.
#   b. Covariate UMAPs: cell type, diagnosis, donor, Sex, Age, Batch.
#   c. Annotated cell-type UMAP (labelled centroids, clean legend) — the hero.
#   d. Canonical marker dot plot (cell type x marker; size = % expressing,
#      fill = scaled mean expr) — justifies new_annotation labels.
#   e. Cohort composition: cell-type fraction of nuclei, CON vs NHD (grouped
#      bars, descriptive; no pseudo-rep stars — the fraction shift is the
#      effect size, per house composition policy).
#
# Cache: data/_cache_SuppFig2_zhou_QC.rds — rebuilt when the source .rds is
# newer (stale-cache guard). Heavy 4.2 GB load runs once; plot tweaks reuse it.
#
# Reproducibility: set.seed(42); provenance (sessionInfo) at end; done sentinel.
# Numerical rigor: no p-values in a QC/annotation figure; Age coerced to numeric
#   with an explicit NA report (never silent).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(tibble); library(patchwork); library(scales); library(forcats)
})
set.seed(42)

# --- Portable project root ---------------------------------------------------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot("MISSING project root" = !is.na(PROJ))

DATA   <- file.path(PROJ, "manuscript_7fig", "data")
# Dedicated SuppFig2 folder (own assembled figure + panels/ + manifest), a
# sibling of SuppFig1_NHD_QC/ — no longer collides in the shared Supplementary/.
FIGDIR <- file.path(PROJ, "manuscript_7fig", "figures", "_proposed_5fig",
                     "Supplementary", "SuppFig2_Zhou_QC")
PANELS <- file.path(FIGDIR, "panels")
CACHE  <- file.path(DATA, "_cache_SuppFig2_zhou_QC.rds")
ZHOU   <- file.path(PROJ, "NHD_repo", "NHD_complete_harmony+PMI.rds")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PANELS, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA,   recursive = TRUE, showWarnings = FALSE)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])), "22_publication_theme.R"))   # the theme file next to this script

# --- Canonical markers (same block-diagonal set as atlas Supp Fig 1) ---------
MARKERS <- list(
  "Micro-PVM"  = c("P2RY12","CSF1R","CX3CR1","C1QB","MRC1","CD163"),
  "Astro"      = c("AQP4","GFAP","SLC1A2"),
  "Oligo"      = c("PLP1","MBP","MOBP"),
  "OPC"        = c("PDGFRA","OLIG1","VCAN"),
  "Neuron_Ex"  = c("RBFOX3","SLC17A7","SATB2"),
  "Neuron_Inh" = c("GAD1","GAD2"),
  "Endo"       = c("CLDN5","FLT1"),
  "Pericytes"  = c("PDGFRB","RGS5"))
CT_ORDER <- names(MARKERS)

# ============================================================================
# Build / reuse cache (rebuild when the Zhou .rds is newer than the cache)
# ============================================================================
stale <- !file.exists(CACHE) ||
         (file.exists(ZHOU) && file.mtime(ZHOU) > file.mtime(CACHE))
if (stale) {
  cat("== building cache (Zhou 4.2 GB load - heavy) ==\n")
  stopifnot("MISSING Zhou object - NHD_complete_harmony+PMI.rds not found" =
              file.exists(ZHOU))
  obj <- readRDS(ZHOU)
  cat(sprintf("Zhou object loaded: %d nuclei x %d genes\n", ncol(obj), nrow(obj)))
  cat("Assays:", paste(Assays(obj), collapse = ","), "\n")
  cat("Reductions:", paste(Reductions(obj), collapse = ","), "\n")

  md <- obj@meta.data
  stopifnot("missing new_annotation" = "new_annotation" %in% names(md),
            "missing Condition"      = "Condition"      %in% names(md),
            "missing SampleID"       = "SampleID"       %in% names(md))

  # mito column: percent.mt is the harmonized field in this object (0-10%,
  # already gated by Zhou); fall back to percent_mito if absent.
  mito_col <- if ("percent.mt" %in% names(md)) "percent.mt" else "percent_mito"
  cat("mito column used:", mito_col, "\n")

  # Age -> numeric with explicit NA accounting (never silent).
  age_num <- suppressWarnings(as.numeric(as.character(md$Age)))
  n_age_na <- sum(is.na(age_num)) - sum(is.na(md$Age))
  if (n_age_na > 0)
    cat(sprintf("NOTE: %d Age values failed numeric coercion -> NA\n", n_age_na))

  meta <- md %>%
    rownames_to_column("cell") %>%
    transmute(
      cell,
      Donor        = as.character(SampleID),
      Diagnosis    = factor(ifelse(as.character(Condition) == "NHD", "NHD", "CON"),
                            levels = c("CON","NHD")),
      cell_type    = as.character(new_annotation),
      nFeature_RNA = nFeature_RNA,
      nCount_RNA   = nCount_RNA,
      percent_mito = md[[mito_col]],
      Sex          = if ("Sex"      %in% names(md)) as.character(Sex)      else NA_character_,
      Age          = age_num,
      Batch        = if ("Batch"    %in% names(md)) as.character(Batch)    else NA_character_,
      Genotype     = if ("Genotype" %in% names(md)) as.character(Genotype) else NA_character_)

  qc <- list(
    n_total    = nrow(meta),
    n_donors   = dplyr::n_distinct(meta$Donor),
    n_nhd      = dplyr::n_distinct(meta$Donor[meta$Diagnosis == "NHD"]),
    n_con      = dplyr::n_distinct(meta$Donor[meta$Diagnosis == "CON"]),
    ct_counts  = sort(table(meta$cell_type), decreasing = TRUE),
    diag_cells = table(meta$Diagnosis),
    mito_col   = mito_col,
    donor_geno = meta %>% distinct(Donor, Diagnosis, Genotype) %>% arrange(Diagnosis, Donor))
  cat("n_total:", qc$n_total, " donors:", qc$n_donors,
      " (NHD:", qc$n_nhd, " CON:", qc$n_con, ")\n")
  cat("cell-type counts:\n"); print(qc$ct_counts)

  # --- Marker dot plot data (RNA lognorm) ------------------------------------
  cat("== computing marker DotPlot data ==\n")
  DefaultAssay(obj) <- "RNA"
  obj <- tryCatch(JoinLayers(obj, assay = "RNA"), error = function(e) obj)  # v5 split layers
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  obj$ct_dot <- as.character(obj$new_annotation)
  sub        <- subset(obj, cells = WhichCells(obj, expression = ct_dot %in% CT_ORDER))
  sub$ct_dot <- factor(sub$ct_dot, levels = CT_ORDER)
  mk_genes   <- intersect(unique(unlist(MARKERS)), rownames(sub))
  missing_mk <- setdiff(unique(unlist(MARKERS)), rownames(sub))
  if (length(missing_mk))
    cat("NOTE markers absent from Zhou object:", paste(missing_mk, collapse=", "), "\n")
  dot <- DotPlot(sub, features = mk_genes, group.by = "ct_dot",
                 assay = "RNA", scale = TRUE)$data %>%
    dplyr::rename(gene = features.plot, cell_type = id,
                  pct_exp = pct.exp, avg_exp = avg.exp.scaled)

  # --- UMAP embedding --------------------------------------------------------
  umap_key <- if ("umap" %in% names(obj@reductions)) "umap" else Reductions(obj)[1]
  cat("UMAP reduction used:", umap_key, "\n")
  umap <- as.data.frame(Embeddings(obj, reduction = umap_key)) %>%
    setNames(c("UMAP_1","UMAP_2")) %>%
    rownames_to_column("cell") %>%
    left_join(meta, by = "cell")

  saveRDS(list(meta = meta, qc = qc, dot = dot, umap = umap,
               markers = MARKERS, missing_mk = missing_mk), CACHE)
  rm(obj, sub); gc()
  cat("== cache written:", basename(CACHE), "==\n")
} else {
  cat("== reusing cache", basename(CACHE), "==\n")
}
CC        <- readRDS(CACHE)
meta      <- CC$meta
qc        <- CC$qc
dot       <- CC$dot
umap      <- CC$umap
missing_mk<- CC$missing_mk

# --- Donor ordering: Diagnosis, then NATURAL (numeric) donor order ----------
# Lexical sort places Control10/Control11 before Control2; extract the leading
# embedded integer so both the violin x-axis and the panel-b donor legend read
# Control1 -> Control11 (and NHD1 -> NHD3) in natural order. Non-numeric donors
# (should be none) sort last, stably.
donor_key <- function(x) {
  n <- suppressWarnings(as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1", x)))
  ifelse(is.na(n), Inf, n)
}
donor_ord <- meta %>% distinct(Donor, Diagnosis) %>%
  arrange(Diagnosis, donor_key(Donor), Donor) %>% pull(Donor)
meta$Donor <- factor(meta$Donor, levels = donor_ord)
umap$Donor <- factor(umap$Donor, levels = donor_ord)

# --- Palettes ----------------------------------------------------------------
# Diagnosis uses the house PAL_COND (CON blue / NHD red).
# Donor palette: blue family for CON, red family for NHD (reads as diagnosis
# blocks). Age continuous. Sex / Batch small discrete sets.
n_con_d <- sum(grepl("^Control", levels(meta$Donor)))
con_donors <- levels(meta$Donor)[grepl("^Control", levels(meta$Donor))]
nhd_donors <- setdiff(levels(meta$Donor), con_donors)
PAL_DONOR <- c(
  setNames(colorRampPalette(c("#9EC9E6","#1B4F72"))(length(con_donors)), con_donors),
  setNames(colorRampPalette(c("#E8998D","#7B241C"))(length(nhd_donors)), nhd_donors))
PAL_SEX   <- c(female = "#C8A2C8", male = "#5985C5")
# Diagnosis CON/NHD colour = the pale microglia module-violin pair (CON pale
# blue / NHD pale red), applied to the QC-violin fills, the CON/NHD group labels
# and the diagnosis covariate UMAP so this supp is consistent with the module
# violins used across all condition UMAPs. Overrides the
# saturated house PAL_COND from 22_publication_theme.R for this figure only.
PAL_COND <- c(CON = "#82B2D6", NHD = "#D6837A")

# cell-type levels present, in canonical order
ct_present <- intersect(CT_ORDER, unique(meta$cell_type))
PAL_CT     <- PAL_CELLTYPE[ct_present]

# display-only relabel of the underscored raw new_annotation levels so panels
# b/c/d/e print full names ("Excitatory neuron" / "Inhibitory neuron") instead
# of Neuron_Ex / Neuron_Inh. Mirrors sibling Supp Fig 1 (70_suppfig1_NHD_QC.R
# :249-254). Factor levels, order and COLOURS are unchanged — this only maps the
# printed text via scale `labels=` (legends / axes) and the centroid `label`.
# The other six names (Micro-PVM, Astro, Oligo, OPC, Endo, Pericytes) pass
# through untouched.
CT_DISPLAY <- function(x) {
  m <- c("Neuron_Ex"  = "Excitatory neuron",
         "Neuron_Inh" = "Inhibitory neuron")
  x <- as.character(x)
  ifelse(x %in% names(m), unname(m[x]), x)
}

# ============================================================================
# Panel a — per-donor QC violins (3 stacked rows: UMIs / Genes / % mito) over
#   one shared donor x-axis, CON vs NHD grouped (dashed divider, white boxplot
#   inset). Mirrors atlas Supp Fig 1 panel a exactly: same facet-stacked layout,
#   same clean strip labels, same PAL_COND fill, same 99th-pctile display cap.
# ============================================================================
ML        <- c(nCount_RNA = "UMIs", nFeature_RNA = "Genes", percent_mito = "% mito")
n_samples <- nlevels(meta$Donor)
ql <- meta %>% select(Donor, Diagnosis, all_of(names(ML))) %>%
  pivot_longer(all_of(names(ML)), names_to = "metric", values_to = "value") %>%
  # cap at per-metric 99th pctile so a few extreme nuclei don't crush the violins
  group_by(metric) %>%
  mutate(value = pmin(value, quantile(value, 0.99, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(metric = factor(ML[metric], levels = ML))
grp_lab <- data.frame(
  metric = factor(ML[1], levels = ML),
  x   = c((1 + n_con_d) / 2, (n_con_d + 1 + n_samples) / 2),
  lab = c("CON","NHD"))
pA <- ggplot(ql, aes(Donor, value, fill = Diagnosis)) +
  geom_vline(xintercept = n_con_d + 0.5, linetype = "22", linewidth = 0.3,
             colour = "grey55") +
  geom_violin(scale = "width", width = 0.78, linewidth = 0.18, colour = "grey30") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.18,
               fill = "white", alpha = 0.9) +
  geom_text(data = grp_lab, aes(x = x, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = -0.35, size = 3.6,
            colour = c(PAL_COND["CON"], PAL_COND["NHD"])) +
  facet_grid(metric ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = PAL_COND, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  # apparent-type parity: saved 4.60 in, placed 2.840 in in
  # Fig2_suppl.pdf => scale 0.6174. Axis ticks were 8 / 7.5 pt = 4.94 / 4.63 pt
  # on the page (mean 4.79); house target is 5.1 pt, so all text x 1.066.
  # Canvas dims unchanged; point/line sizes and unit() dims not scaled.
  theme_pub(base_size = 9.6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
        axis.text.y = element_text(size = 8.0),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, size = 11.2),
        panel.spacing = unit(3, "pt"), plot.margin = margin(11, 3, 2, 2))

# ============================================================================
# UMAP scaffold (downsample + shuffle draw order for fair overplotting)
# ============================================================================
set.seed(42)
idx <- sample(nrow(umap), min(45000, nrow(umap)))
umap_sub <- umap[idx, ]
umap_sub <- umap_sub[sample(nrow(umap_sub)), ]

# ============================================================================
# apparent-type parity — house rule:
#     on-page pt = declared pt x placement scale,  scale = placed_w / saved_w
# The NHD main figures sit at 4.8-5.3 pt on the page; target 5.1 pt.
# Measured from the assembled Fig2_suppl.pdf:
#     SuppFig2_b_celltype       saved 3.40" -> placed 2.698"  scale 0.794
#     SuppFig2_a_covariate_UMAP saved 8.60" -> placed 3.675"  scale 0.427
#     SuppFig2_a_QC_violins     saved 4.60" -> placed 2.840"  scale 0.617
#     SuppFig2_c_marker_dotplot saved 5.60" -> placed 3.889"  scale 0.694
# The same covariate UMAPs are saved at two canvases placed at very different
# scales, so they need different declared type sizes:
#   individual 3.40 x 2.80 in : legend 8.5 pt x 0.794 = 6.75 pt  -> x 0.756
#   2x2 composite 8.60 x 5.20 : legend 8.5 pt x 0.427 = 3.63 pt  -> x 1.405
# canvas DIMENSIONS are unchanged; only declared text sizes move. geom_point size / stroke /
# legend key unit() dims are deliberately not scaled.
# ============================================================================
FS_IND  <- 0.756   # individually-placed covariate UMAPs (3.40 x 2.80 in)
FS_COMP <- 1.405   # 2x2 covariate composite             (8.60 x 5.20 in)

# ============================================================================
# axis_arrows() — the house UMAP axis key, an L of two closed arrows in the lower
# left labelled "UMAP 1" / "UMAP 2", in place of axes on an embedding. Copied from
# NHD_frontal_hippo_rebuild/scripts/22_publication_theme_FH.R (which took it
# verbatim from Figure 5's 64_fig4_neuron_umap_FH.R) so this figure cannot drift
# from the main deck; this tree sources a different theme, hence the local copy.
#
# House RULE: one key per composite, at the reading entry point -- the leftmost
# panel side-by-side, the lower panel stacked, the lower-left tile of a grid.
# Never repeat it on every sub-panel. A panel that stands alone carries its own.
# The panel needs clip = "off"; the labels sit just outside the data range.
# ============================================================================
axis_arrows <- function(xr, yr, frac = 0.16, lab = c("UMAP 1", "UMAP 2"),
                        lab_size = 2.5, head_in = 0.055, linewidth = 0.4) {
  rng <- max(diff(xr), diff(yr))
  x0  <- xr[1] - 0.01 * diff(xr); y0 <- yr[1] - 0.01 * diff(yr)
  a   <- grid::arrow(length = grid::unit(head_in, "in"), type = "closed")
  list(
    ggplot2::annotate("segment", x = x0, xend = x0 + frac * rng, y = y0, yend = y0,
                      arrow = a, linewidth = linewidth),
    ggplot2::annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + frac * rng,
                      arrow = a, linewidth = linewidth),
    ggplot2::annotate("text", x = x0, y = y0 - 0.03 * rng, hjust = 0, vjust = 1,
                      label = lab[1], size = lab_size, colour = "black"),
    ggplot2::annotate("text", x = x0 - 0.03 * rng, y = y0, angle = 90, hjust = 0,
                      vjust = 0, label = lab[2], size = lab_size, colour = "black"))
}

# Same lower-left breathing room on every covariate UMAP, so the key sits in clear
# space and all tiles of the 2x2 stay at one zoom (expanding only the tile that
# carries the key makes the grid stop reading as a single embedding).
UMAP_EXPAND <- ggplot2::expansion(mult = c(0.13, 0.02))

# Thin dark outline on every dot: shape 21 + grey20 border;
# each ggplot() passed in maps the data covariate to `fill` (not colour).
umap_base <- function(p, ptsize = 0.42, fs = 1) p +
  geom_point(size = ptsize, alpha = 0.85, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_x_continuous(expand = UMAP_EXPAND) +
  scale_y_continuous(expand = UMAP_EXPAND) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_cartesian(clip = "off") +
  theme_pub(base_size = 9 * fs) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_text(size = 8 * fs))

# Swap the written axis titles for the house arrow key. theme_pub sets
# axis.title.x/.y as explicit element_text children, so blanking the parent
# axis.title is silently ignored -- blank the children.
no_axis_title <- function(p) p + theme(axis.title.x = element_blank(),
                                       axis.title.y = element_blank())
add_axis_key <- function(p, df, fs) {
  no_axis_title(p) +
    axis_arrows(xr = range(df$UMAP_1), yr = range(df$UMAP_2),
                lab_size  = 9 * fs / 2.845,
                head_in   = 0.055 * 9 * fs / 8,
                linewidth = 0.4 * 9 * fs / 8)
}

# ============================================================================
# Covariate UMAP factory — one builder, two font scales (see FS_* above).
# Returns the six covariate UMAPs (cell type, diagnosis, donor, Sex, Age, Batch)
# with every text size multiplied by `fs`. Base declared sizes (fs = 1) are the
# ============================================================================
mk_cov_umaps <- function(fs) {
  ttl  <- function() element_text(size = 11 * fs, hjust = 0.5, face = "plain")
  bigdot <- list(shape = 21, colour = "grey20", stroke = 0.3, size = 3, alpha = 1)

  p_ct <- umap_base(ggplot(umap_sub, aes(UMAP_1, UMAP_2,
                                         fill = factor(cell_type, levels = ct_present))),
                    fs = fs) +
    scale_fill_manual(values = PAL_CT, name = NULL, drop = FALSE, labels = CT_DISPLAY) +
    guides(fill = guide_legend(ncol = 1, override.aes = bigdot)) +
    labs(subtitle = NULL) + ggtitle("Cell type") +
    theme(plot.title = ttl(),
          legend.text = element_text(size = 8.5 * fs),
          legend.key.size = unit(0.38, "cm"))

  p_diag <- umap_base(ggplot(umap_sub, aes(UMAP_1, UMAP_2, fill = Diagnosis)), fs = fs) +
    scale_fill_manual(values = PAL_COND, name = NULL) +
    guides(fill = guide_legend(override.aes = bigdot)) +
    ggtitle("Diagnosis") +
    theme(plot.title = ttl(),
          legend.text = element_text(size = 8.5 * fs),
          legend.key.size = unit(0.38, "cm"))

  # Donor has 14 levels -> 2-column key (CON blue / NHD red families read as two
  # columns) so the legend stays legible without overrunning the panel.
  p_donor <- umap_base(ggplot(umap_sub, aes(UMAP_1, UMAP_2, fill = Donor)), fs = fs) +
    scale_fill_manual(values = PAL_DONOR, name = NULL) +
    guides(fill = guide_legend(ncol = 2,
                               override.aes = modifyList(bigdot, list(size = 2.4)))) +
    ggtitle("Donor") +
    theme(plot.title = ttl(),
          legend.text = element_text(size = 7.2 * fs),
          legend.key.size = unit(0.30, "cm"))

  p_sex <- umap_base(ggplot(umap_sub, aes(UMAP_1, UMAP_2, fill = Sex)), fs = fs) +
    scale_fill_manual(values = PAL_SEX, name = NULL, na.value = "grey80") +
    guides(fill = guide_legend(override.aes = bigdot)) +
    ggtitle("Sex") +
    theme(plot.title = ttl(),
          legend.text = element_text(size = 8.5 * fs),
          legend.key.size = unit(0.38, "cm"))

  p_age <- umap_base(ggplot(umap_sub, aes(UMAP_1, UMAP_2, fill = Age)), fs = fs) +
    scale_fill_viridis_c(option = "C", name = "Age (y)", na.value = "grey85") +
    ggtitle("Age") +
    theme(plot.title = ttl(),
          legend.text  = element_text(size = 8   * fs),
          legend.title = element_text(size = 8.5 * fs),
          legend.key.size = unit(0.38, "cm"))

  p_batch <- umap_base(ggplot(umap_sub, aes(UMAP_1, UMAP_2, fill = factor(Batch))), fs = fs) +
    scale_fill_brewer(palette = "Dark2", name = NULL, na.value = "grey80") +
    guides(fill = guide_legend(override.aes = bigdot)) +
    ggtitle("Batch") +
    theme(plot.title = ttl(),
          legend.text = element_text(size = 8.5 * fs),
          legend.key.size = unit(0.38, "cm"))

  list(ct = p_ct, diag = p_diag, donor = p_donor,
       sex = p_sex, age = p_age, batch = p_batch)
}

IND  <- mk_cov_umaps(FS_IND)    # for the 3.40 x 2.80 individual panel saves
COMP <- mk_cov_umaps(FS_COMP)   # for the 8.60 x 5.20 2x2 composite

pB_ct    <- IND$ct;    pB_diag  <- IND$diag;  pB_donor <- IND$donor
pB_sex   <- IND$sex;   pB_age   <- IND$age;   pB_batch <- IND$batch

# The standalone annotated cell-type "hero" UMAP (former panel c) was
# removed — it is redundant with the "Cell type" UMAP already in the panel-b
# covariate composite. Remaining panels re-lettered a (QC) / b (covariate
# composite) / c (marker dotplot) / d (composition).

# ============================================================================
# Panel c — canonical marker dot plot
# ============================================================================
dot$cell_type <- factor(dot$cell_type, levels = rev(CT_ORDER))
dot$gene      <- factor(dot$gene, levels = intersect(unlist(MARKERS),
                                                     levels(dot$gene)))
blk_sizes <- vapply(MARKERS, function(g) sum(g %in% levels(dot$gene)), integer(1))
sep_x     <- head(cumsum(blk_sizes), -1) + 0.5

pD_dot <- ggplot(dot, aes(gene, cell_type, size = pct_exp, fill = avg_exp)) +
  geom_vline(xintercept = sep_x, colour = "grey88", linewidth = 0.25) +
  geom_point(shape = 21, stroke = 0.25, colour = "grey20", alpha = 0.95) +
  scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#C0392B",
                       midpoint = 0, name = "scaled\nmean expr",
                       limits = c(-1.5, 1.5), oob = scales::squish) +
  scale_size_continuous(range = c(0.4, 5), name = "% cells\nexpressing",
                        breaks = c(10, 30, 60, 90), limits = c(0, 100)) +
  scale_y_discrete(labels = CT_DISPLAY) +      # Neuron_Ex/Inh -> full names
  labs(x = NULL, y = NULL) +
  # apparent-type parity: saved 5.60 in, placed 3.889 in in
  # Fig2_suppl.pdf => scale 0.694. Axis ticks were 8.5 / 9.5 pt = 5.90 / 6.60 pt
  # on the page (mean 6.25) - the largest type in the figure. House target
  # 5.1 pt, so all text x 0.816. Canvas unchanged; dot size range and key
  # unit() dims not scaled.
  theme_pub(base_size = 7.35) +
  theme(axis.text.x = element_text(size = 6.9, angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(size = 7.8),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.25),
        legend.key.size = unit(0.38, "cm"),
        legend.title = element_text(size = 7.4),
        legend.text  = element_text(size = 7.0))

# Final SuppFig2 = three panels: a (QC violins) / b (covariate UMAP
# composite) / c (canonical marker dot plot).

# ============================================================================
# Save individual panels
# ============================================================================
save_panel <- function(p, name, w, h) {
  ggsave(file.path(PANELS, paste0(name, ".pdf")), p, width = w, height = h,
         useDingbats = FALSE, limitsize = FALSE)
  ggsave(file.path(PANELS, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, limitsize = FALSE, device = ragg::agg_png)
}
save_panel(pA, "SuppFig2_a_QC_violins", 4.6, 3.6)
# Panel-b covariate UMAPs saved individually (enlarged for legibility) + as the
# assembled 6-UMAP composite below.
# stands alone in the assembly, so it carries its own key rather than borrowing one
save_panel(add_axis_key(pB_ct, umap_sub, FS_IND), "SuppFig2_b_celltype", 3.4, 2.8)
save_panel(pB_diag,  "SuppFig2_b_diagnosis", 3.4, 2.8)
save_panel(pB_donor, "SuppFig2_b_donor", 3.4, 2.8)
save_panel(pB_sex,   "SuppFig2_b_sex", 3.4, 2.8)
save_panel(pB_age,   "SuppFig2_b_age", 3.4, 2.8)
save_panel(pB_batch, "SuppFig2_b_batch", 3.4, 2.8)
# The 3360x1920 px raster still linked in Fig2_suppl.pdf is the PRE-request render;
# this panel is the one place in SuppFig2 that does need re-placing in Illustrator.
save_panel(pD_dot,   "SuppFig2_c_marker_dotplot", 5.54, 2.88)

# 2x2 COVARIATE composite, in the SuppFig1_b_covariateUMAP grammar:
# plain title above each UMAP, legend to its right, no axes. Age / Diagnosis / Sex /
# Donor -- the four donor-level covariates that actually vary in the Zhou cohort
# (11 controls + 3 NHD donors, ages 50-80, both sexes), which is why this panel belongs
# to the Zhou supplement and not to ours: in our own atlas there is one donor per arm, so
# age and donor are collinear with diagnosis and sex is constant.
# Composite members use the LARGER font scale (FS_COMP) because this canvas is
# placed at scale 0.427 - see the apparent-type parity block above.
# And that
# reference panel carries no axis titles: a plain title above each UMAP, the key to its
# right, and nothing else. Four repetitions of "UMAP 1 / UMAP 2" say nothing a reader
# needs and eat the space the keys want, so they come off here (the individually placed
# SuppFig2_b_* UMAPs keep theirs -- they stand alone).
# ncol = 2 lays these out  age | diagnosis  /  sex | donor,  so the lower-left tile
# -- the reading entry point -- is Sex. It carries the axis key; the other three
# show no axis furniture at all, exactly as in SuppFig1_b.
pB_cov <- patchwork::wrap_plots(
            list(no_axis_title(COMP$age),
                 no_axis_title(COMP$diag),
                 add_axis_key(COMP$sex, umap_sub, FS_COMP),
                 no_axis_title(COMP$donor)),
            ncol = 2)
save_panel(pB_cov, "SuppFig2_a_covariate_UMAP", 8.6, 5.2)

# ============================================================================
# Assemble full figure — no caption banner (caption lives in the legend .md),
# no bold, plain tags. three panels, each with generous real estate:
#   a  QC violins (full-width row)
#   b  covariate UMAP composite (6 UMAPs, 2x3, full-width row — the main block)
#   c  canonical marker dot plot (full-width row)
# The 6-UMAP grid is wrapped as one taggable element so plot_annotation does
# not recurse into its mini-UMAPs and steal the b/c tags.
# ============================================================================
pB_grid <- (pB_ct | pB_diag | pB_donor) / (pB_sex | pB_age | pB_batch)
pB1  <- patchwork::wrap_elements(full = pB_grid)
row1 <- pA                 # a  QC violins, full width
row2 <- pB1                # b  covariate UMAP composite, full width (main block)
row3 <- pD_dot             # c  canonical marker dot plot, full width

figS2 <- row1 / row2 / row3 +
  plot_layout(heights = c(0.80, 1.45, 0.62)) +
  plot_annotation(tag_levels = list(c("a", "b", "c"))) &
  # Panel LETTERS harmonized to the manually-assembled main figures (Fig4.pdf):
  # Big, bold, black, top-left. This is the one allowed bold element (panel tags);
  # all other text (titles/axes/legends/strips) stays face="plain".
  theme(plot.tag = element_text(face = "bold", size = 18, colour = "black"))

# Illustrator (SuppFig2_Zhou_QC.{pdf,png}) with his own panel order. This script
# now writes only an *_AUTOPREVIEW copy so re-rendering panels never clobbers his
# manual composite. Individual panels/ PNGs (which he places) still refresh above.
W <- 10.0; H <- 12.6
ggsave(file.path(FIGDIR, "SuppFig2_Zhou_QC_AUTOPREVIEW.pdf"), figS2,
       width = W, height = H, useDingbats = FALSE, limitsize = FALSE)
ggsave(file.path(FIGDIR, "SuppFig2_Zhou_QC_AUTOPREVIEW.png"), figS2,
       width = W, height = H, dpi = 600, limitsize = FALSE, device = ragg::agg_png)

# ============================================================================
# Manifest
# ============================================================================
geno_lines <- apply(qc$donor_geno, 1, function(r)
  sprintf("  - %-16s %-5s %s", r["Donor"], r["Diagnosis"],
          ifelse(is.na(r["Genotype"]), "", r["Genotype"])))
manifest <- c(
  "# Supplementary Figure 2 — Zhou 2023 external cohort QC + annotation",
  "",
  "**Cohort:** Zhou et al. 2023 *Nat Immunol* — human OCCIPITAL cortex snRNA-seq.",
  sprintf("**Design:** %d NHD (DAP12-null) vs %d CON donors; %s nuclei; %d annotated cell types.",
          qc$n_nhd, qc$n_con, format(qc$n_total, big.mark = ","), length(qc$ct_counts)),
  sprintf("**Diagnosis cells:** CON = %s, NHD = %s.",
          format(qc$diag_cells["CON"], big.mark = ","),
          format(qc$diag_cells["NHD"], big.mark = ",")),
  sprintf("**Mito field:** %s (Zhou-gated, 0-10%%).", qc$mito_col),
  "",
  "## Donors (diagnosis / genotype)",
  geno_lines,
  "",
  "## Legend (caption — belongs here, NOT embedded in the figure)",
  paste("**Supplementary Figure 2. External-cohort QC and cell-type",
        "annotation (Zhou et al. 2023, occipital cortex; 3 NHD DAP12-null vs",
        "11 CON).** (a) Per-donor QC violins (UMIs, Genes, % mitochondrial),",
        "one violin per donor, CON vs NHD grouped (dashed divider) with white",
        "boxplot inset; values capped at the per-metric 99th percentile for",
        "display; fill = diagnosis (CON blue / NHD red). (b) Covariate UMAPs",
        "(cell type, diagnosis, donor, Sex, Age, Batch) on the Harmony",
        "embedding, arranged as a 2x3 composite. (c)",
        "Canonical marker dot plot (size = % expressing, fill = scaled mean",
        "expression)."),
  "",
  "## Panels",
  "- **a** — per-donor QC violins: UMIs, Genes, % mito; 3 stacked rows over one shared donor x-axis; CON vs NHD grouped (dashed divider), white boxplot inset; fill = pale module-violin PAL_COND (CON #82B2D6 / NHD #D6837A, consistent with the condition UMAPs).",
  "- **b** — covariate UMAP composite (2x3): cell type, diagnosis, donor, Sex, Age, Batch (Harmony-integrated embedding); enlarged high-contrast points + legible legends (SuppFig3 Visium UMAP style).",
  "- **c** — canonical marker dot plot (cell type x marker; size = % expressing, fill = scaled mean expr).",
  "",
  "## Cell-type counts",
  paste0("  - ", names(qc$ct_counts), ": ", format(as.integer(qc$ct_counts), big.mark = ",")),
  "",
  if (length(missing_mk)) paste0("**Markers absent from object:** ",
                                 paste(missing_mk, collapse = ", ")) else
    "All canonical markers present in the Zhou object.",
  "",
  "## Files",
  "- `SuppFig2_Zhou_QC.{pdf,png}` — assembled figure.",
  "- `panels/` — individual panel PDFs + PNGs (a QC violins, b covariate UMAPs [6], c marker dotplot).",
  "",
  "## Build",
  "- Script: `manuscript_7fig/scripts/70_suppfig2_zhou_QC.R`",
  "- Source object: `NHD_repo/NHD_complete_harmony+PMI.rds`",
  "- Cache: `manuscript_7fig/data/_cache_SuppFig2_zhou_QC.rds` (mtime-guarded).",
  "- Style: `scripts/22_publication_theme.R` (theme_pub, PAL_CELLTYPE, PAL_COND);",
  "  visual sibling of Supp Fig 1 (`scripts/70_suppfig1_NHD_QC.R`).",
  sprintf("- Built: %s", format(Sys.time(), "%Y-%m-%d %H:%M")))
writeLines(manifest, file.path(FIGDIR, "SuppFig2_panels.md"))

cat("\n=== DONE ===\n")
cat("wrote:\n  ", file.path(FIGDIR, "SuppFig2_Zhou_QC_AUTOPREVIEW.pdf"),
    " (preview; the assembled figure is laid out by hand)",
    "\n  ", file.path(FIGDIR, "SuppFig2_Zhou_QC_AUTOPREVIEW.png"),
    "\n  ", file.path(FIGDIR, "SuppFig2_panels.md"),
    "\n  panels/ -> ", PANELS, "\n", sep = "")
cat("\n-- provenance --\n"); print(utils::sessionInfo())
