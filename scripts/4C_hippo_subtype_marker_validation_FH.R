#!/usr/bin/env Rscript
# =============================================================================
# 4C_hippo_subtype_marker_validation_FH.R — Figure 4 panel c:
# hippocampal-neuron annotation-validation dot-plot, frontal + hippocampus rebuild.
# Faithful port of manuscript_7fig/scripts/65_fig4_hippo_subtype_marker_dotplot.R.
# -----------------------------------------------------------------------------
# Why: no external NHD-hippocampus DE cohort exists, so this is the honest
# substitute for a hippo cross-cohort validator: it shows our 6 de-novo
# hippocampal subtypes are each defined by canonical human-hippocampus private
# markers (a clean block-diagonal) — the annotation is marker-grounded, not circular.
#
# What: dot-plot. y = canonical private markers, row-blocked by the subtype they
# mark (grey left strips). x = the 6 final FH hippo subtypes (from 60). dot fill =
# scaled avg expression (z across subtypes, DotPlot semantics, clipped +/-2.5). dot
# SIZE = % of nuclei expressing.  pan-gate block (SLC17A7/GAD1/GAD2) = Ex vs Inh.
#
# What changed vs the 3-region reference (6 subtypes, not 7): EC-like did not
# resolve on the FH atlas (60_neuron_reannotation_FH), so its subtype column and its
# EC-like marker column is DROPPED (reported).  All other blocks retained.
# Subiculum block uses SLC17A6/ROBO1/TLE4.
#
# DATA: subtype labels = per-barcode `neuron_subtype` in
# data/_cache_hippo_neuron_umap_FH.rds.  Expression = atlas hippo-neuron subset,
# read once if the small marker cache is stale, then cached
# (data/_cache_hippo_marker_dotplot_FH.rds) with atlas-mtime + marker-set + label-set
# guard so re-runs are atlas-light and reproducible.  seed 42.
#
# Output: figures/Figure_5/panels/F5f2_hippo_subtype_markers.{png,pdf}
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")
# =============================================================================
# apparent-type parity PASS  ("adjust the font/panel size ratio
#   so everything looks more or less uniform, favouring the readability of the figure")
# Measured from the assembled Fig5.pdf: each panel is placed as a linked image with a
# known scale (parsed from the page's placement matrices). On-page type was 3.2-3.9 pt
# for most panels but 6.4-6.9 pt for the chord -- a 2x spread, and most of it below the
# print-legibility floor. Declared sizes here are therefore set to
#     declared = apparent_target / placement_scale
# with one set of apparent targets shared by every Figure-5 panel:
#     tick / legend text 5.2 pt | legend title 5.4 | axis title & facet strip 6.0
#     panel title 6.4 | in-panel gene labels & annotations 5.4
# canvas size and PIXEL DIMENSIONS are unchanged, so re-linking in Illustrator is a
# no-op: the frame keeps its box and only the type grows.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggh4x); library(stringr)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("MISSING project root" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))

DDIR  <- file.path(PROJ, "data")
PANEL <- file.path(PROJ, "figures", "Figure_5", "panels")
LOGD  <- file.path(PROJ, "logs")
for (d in c(DDIR, PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

HIP_CACHE  <- file.path(DDIR, "_cache_hippo_neuron_umap_FH.rds")     # frozen labels (from 60)
ATLAS      <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
MARK_CACHE <- file.path(DDIR, "_cache_hippo_marker_dotplot_FH.rds")
stopifnot("MISSING hippo label cache (run 60_neuron_reannotation_FH.R first)" = file.exists(HIP_CACHE),
          "MISSING atlas NHD_FH_harmony.rds" = file.exists(ATLAS))

# =============================================================================
# Canonical private-marker blocks (6 FH subtypes; EC-like dropped — not resolved).
# Subtype (x) order == block (y) order so on-diagonal dots line up. pan-gate last.
# =============================================================================
# The cluster formerly labelled "Cajal-Retzius" failed the private-marker
# gate in 60_neuron_reannotation_FH.R and is now an unresolved low-complexity cluster.
# Its column stays on this panel and the canonical Cajal-Retzius marker block stays with
# it, because together they are the evidence: the column carries RELN (promiscuous) but
# not TP73 (private), which is exactly why the label was withdrawn.
UNRESOLVED_LAB <- "Unresolved"
SUBTYPE_ORDER <- c("CA3", "Subiculum", UNRESOLVED_LAB,
                   "Inh-CGE (VIP)", "Inh-CGE (LAMP5)", "Inh-MGE")
SUBTYPE_DISPLAY <- c("CA3" = "CA3", "Subiculum" = "Subiculum",
                     "Inh-CGE (VIP)" = "Inh-CGE (VIP)",
                     "Inh-CGE (LAMP5)" = "Inh-CGE (LAMP5)", "Inh-MGE" = "Inh-MGE")
# x labels are angled on a 2.8-in panel: a two-line label overprints its neighbour
# ("(low complexity)" collided straight through "Inh-CGE (VIP)"). One word on the panel;
# "low complexity" is defined in the figure legend.
SUBTYPE_DISPLAY[UNRESOLVED_LAB] <- "Unresolved"
MARKER_BLOCKS <- list(
  "CA3"             = c("HS3ST4", "GRIK4", "NECAB1"),
  # Subiculum: SLC17A6 (VGLUT2) is a cleaner deep-projection/subiculum marker than
  # FN1. TLE4/ROBO1
  # already carry the block. SLC17A6 also frees the marker cache guard to rebuild.
  "Subiculum"       = c("SLC17A6", "ROBO1", "TLE4"),
  "Cajal-Retzius"   = c("TP73", "RELN", "LHX1"),
  "Inh-CGE (VIP)"   = c("ADARB2", "VIP", "CCK", "CNR1"),
  "Inh-CGE (LAMP5)" = c("LAMP5", "ID2", "NDNF"),
  "Inh-MGE"         = c("LHX6", "SST", "PVALB", "GAD1", "GAD2"),
  "pan-gate"        = c("SLC17A7", "GAD1", "GAD2"))
BLOCK_ORDER   <- names(MARKER_BLOCKS)
BLOCK_DISPLAY <- c(
  "CA3" = "CA3", "Subiculum" = "Subiculum", "Cajal-Retzius" = "Cajal-\nRetzius",
  "Inh-CGE (VIP)"   = "Inh-CGE\n(VIP)",  "Inh-CGE (LAMP5)" = "Inh-CGE\n(LAMP5)",
  "Inh-MGE"         = "Inh-MGE",         "pan-gate"        = "Ex / Inh")
row_tbl <- do.call(rbind, lapply(BLOCK_ORDER, function(b)
  data.frame(block = b, gene = MARKER_BLOCKS[[b]], stringsAsFactors = FALSE)))
row_tbl$key <- paste(row_tbl$block, row_tbl$gene, sep = "__")
UNIQUE_MARKERS <- unique(row_tbl$gene)

# =============================================================================
# Load frozen labels; keep only the 6 hippo subtypes.
# =============================================================================
hip <- readRDS(HIP_CACHE)
stopifnot("hippo cache missing meta$neuron_subtype" = "neuron_subtype" %in% colnames(hip$meta))
lab_vec <- setNames(hip$meta$neuron_subtype, rownames(hip$meta))
stopifnot("expected exactly 6 hippo subtype labels" =
            setequal(unique(lab_vec), SUBTYPE_ORDER))
cat("== hippo subtype labels (frozen, from 60) ==\n")
print(table(factor(lab_vec, levels = SUBTYPE_ORDER)))
cat("NOTE: EC-like NOT resolved on FH atlas -> 6 subtypes (its EC-like marker column dropped)\n")

# =============================================================================
# Marker dot table: reload from small cache if fresh, else compute + cache.
# =============================================================================
atlas_mtime <- file.mtime(ATLAS)
fresh_cache <- FALSE
if (file.exists(MARK_CACHE)) {
  mc <- readRDS(MARK_CACHE)
  fresh_cache <- !is.null(mc$atlas_mtime) && mc$atlas_mtime >= atlas_mtime &&
    all(UNIQUE_MARKERS %in% (mc$markers_requested %||% character(0))) &&
    setequal(mc$labels %||% character(0), SUBTYPE_ORDER)
  if (!fresh_cache) cat("== marker cache STALE (atlas newer / marker or label set changed) -> rebuild ==\n")
}
if (fresh_cache) {
  cat("== marker cache fresh -> ATLAS-LIGHT reload ==\n")
  dp <- mc$dp; absent <- mc$markers_absent
} else {
  cat("== loading atlas (marker cache stale) ==\n")
  suppressPackageStartupMessages(library(Seurat))
  atl <- readRDS(ATLAS)
  bc  <- intersect(names(lab_vec), colnames(atl))
  cat(sprintf("   hippo-neuron barcodes: %d of %d present in atlas\n", length(bc), length(lab_vec)))
  stopifnot("no hippo barcodes found in atlas — barcode mismatch" = length(bc) > 0)
  o <- subset(atl, cells = bc); rm(atl); gc()
  DefaultAssay(o) <- "RNA"
  if (length(SeuratObject::Layers(o, assay = "RNA")) > 1) o <- SeuratObject::JoinLayers(o)
  o <- NormalizeData(o, assay = "RNA", verbose = FALSE)
  o$hip_label <- factor(lab_vec[colnames(o)], levels = SUBTYPE_ORDER)
  stopifnot("NA labels after subset" = !anyNA(o$hip_label))

  present <- intersect(UNIQUE_MARKERS, rownames(o))
  absent  <- setdiff(UNIQUE_MARKERS, rownames(o))
  if (length(absent)) cat("   markers ABSENT from atlas (skipped):", paste(absent, collapse = ", "), "\n")
  stopifnot("no requested markers present in atlas" = length(present) > 0)

  dpo <- Seurat::DotPlot(o, features = present, group.by = "hip_label", col.min = -2.5, col.max = 2.5)
  dp  <- dpo$data
  dp  <- data.frame(gene = as.character(dp$features.plot), subtype = as.character(dp$id),
                    avg.exp = dp$avg.exp, pct.exp = dp$pct.exp, scaled = dp$avg.exp.scaled,
                    stringsAsFactors = FALSE)
  n_nan <- sum(is.nan(dp$scaled))
  if (n_nan > 0) {
    inv <- unique(dp$gene[is.nan(dp$scaled)])
    cat(sprintf("   %d (gene x subtype) scaled values NaN (invariant gene[s]: %s) -> set to 0 (neutral)\n",
                n_nan, paste(inv, collapse = ", ")))
    dp$scaled[is.nan(dp$scaled)] <- 0
  }
  saveRDS(list(dp = dp, markers_requested = UNIQUE_MARKERS, markers_absent = absent,
               labels = SUBTYPE_ORDER, atlas_mtime = atlas_mtime, built = Sys.time()),
          paste0(MARK_CACHE, ".tmp"))
  file.rename(paste0(MARK_CACHE, ".tmp"), MARK_CACHE)
  cat("   wrote", basename(MARK_CACHE), "\n")
  rm(o); gc()
}

# =============================================================================
# Assemble plotting frame; join per-gene dot stats onto the block table.
# =============================================================================
plot_df <- row_tbl %>%
  left_join(dp, by = "gene", relationship = "many-to-many") %>%
  filter(!is.na(scaled))
stopifnot("empty plot frame" = nrow(plot_df) > 0)
plot_df$subtype <- factor(plot_df$subtype, levels = SUBTYPE_ORDER)
plot_df$block   <- factor(plot_df$block,   levels = BLOCK_ORDER)
present_keys    <- row_tbl$key[row_tbl$gene %in% dp$gene]
plot_df$key     <- factor(plot_df$key, levels = rev(present_keys))
gene_of_key     <- setNames(row_tbl$gene, row_tbl$key)

cat(sprintf("plot frame: %d tiles (%d gene rows x <= %d subtypes)\n",
            nrow(plot_df), nlevels(droplevels(plot_df$key)), length(SUBTYPE_ORDER)))
cat(sprintf("scaled-expr range: [%.2f, %.2f]  pct.exp range: [%.1f, %.1f]\n",
            min(plot_df$scaled), max(plot_df$scaled), min(plot_df$pct.exp), max(plot_df$pct.exp)))

# =============================================================================
# House marker-dotplot: grey left-strip pills, italic gene y-labels.
# =============================================================================
GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
strip_spec <- strip_nested(
  background_y = elem_list_rect(fill = "grey90", colour = NA),
  text_y       = elem_list_text(colour = "black", face = "plain", size = 9.4, angle = 0, lineheight = 0.82))

# Detection FLOOR. `scaled` is a z across only six columns, so a gene detected in ~1%
# of nuclei everywhere still yields a confident-looking red in whichever column happens
# to be highest -- TP73 rendered red in the Cajal-Retzius column at 1.1% detected. Dots
# below DET_FLOOR are therefore drawn hollow and unfilled: the colour channel is reserved
# for values that clear detection. Size still encodes %detected in both layers.
DET_FLOOR <- 10
plot_lo <- plot_df[plot_df$pct.exp <  DET_FLOOR, , drop = FALSE]
plot_hi <- plot_df[plot_df$pct.exp >= DET_FLOOR, , drop = FALSE]
cat(sprintf("detection floor %d%%: %d of %d dots drawn hollow (below floor)\n",
            DET_FLOOR, nrow(plot_lo), nrow(plot_df)))

p <- ggplot(plot_df, aes(x = subtype, y = key)) +
  geom_point(data = plot_lo, aes(size = pct.exp), shape = 21, fill = "white",
             colour = "grey70", stroke = 0.25) +
  geom_point(data = plot_hi, aes(fill = scaled, size = pct.exp), shape = 21,
             colour = "grey35", stroke = 0.25) +
  facet_nested(rows = vars(block), scales = "free_y", space = "free_y", switch = "y",
               strip = strip_spec, labeller = labeller(block = as_labeller(BLOCK_DISPLAY))) +
  scale_fill_gradient2(low = GSEA_STEEL, mid = "grey95", high = GSEA_ROSE, midpoint = 0,
                       limits = c(-2.5, 2.5), breaks = c(-2, 0, 2), name = "Scaled expr") +
  scale_size_continuous(range = c(0.3, 4.6), limits = c(0, 100), breaks = c(25, 75), name = "% expr") +
  scale_y_discrete(labels = function(k) gene_of_key[k]) +
  # Subtype labels moved to the bottom. Rotated bottom labels anchor
  # at their right end (hjust = 1, vjust = 1) -- the mirror of the top-axis settings.
  scale_x_discrete(labels = function(x) unname(SUBTYPE_DISPLAY[x])) +
  guides(fill = guide_colourbar(order = 1, title.position = "left", title.vjust = 0.9,
                                barwidth = unit(1.15, "cm"), barheight = unit(0.26, "cm")),
         size = guide_legend(order = 2, title.position = "left", nrow = 1)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 14.2) +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, vjust = 1, size = 9.3, colour = "black"),
    axis.text.y     = element_text(size = 8.0, face = "italic", colour = "black"),
    axis.line       = element_blank(), axis.ticks = element_line(linewidth = 0.25),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.spacing.y = unit(1.5, "pt"), strip.placement = "outside", strip.clip = "off",
    legend.position = "bottom", legend.box = "vertical", legend.key.size = unit(0.30, "cm"),
    legend.box.just = "center", legend.justification = "center",
    legend.text = element_text(size = 8.2), legend.title = element_text(size = 8.2, lineheight = 0.85),
    legend.margin = margin(t = 0, b = 0), legend.box.spacing = unit(0.10, "cm"),
    legend.box.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.02, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    plot.margin = margin(4, 4, 4, 4))

W <- 2.79; H <- 4.2   # 20% then a further 15% slimmer
BN <- "F5f2_hippo_subtype_markers"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600, device = ragg::agg_png)
cat(sprintf("Saved %s.{pdf,png}  (%.1f x %.1f in)\n", BN, W, H))

# =============================================================================
# On-diagonal sanity check (printed).
# =============================================================================
cat("\n== on-diagonal check: mean scaled expr of each subtype's OWN markers ==\n")
for (s in SUBTYPE_ORDER) {
  own <- plot_df %>% filter(block == s)
  if (!nrow(own)) next
  by_sub <- tapply(own$scaled, own$subtype, mean)
  top <- names(which.max(by_sub))
  cat(sprintf("   %-18s own-marker argmax = %-18s (self=%.2f) %s\n",
              s, top, by_sub[[s]], ifelse(top == s, "OK", "<-- check")))
}

# =============================================================================
# provenance
# =============================================================================
prov <- file.path(LOGD, "4C_hippo_subtype_marker_validation_FH_provenance.txt")
sink(prov)
cat("4C_hippo_subtype_marker_validation_FH.R provenance\n")
cat("run time:", format(Sys.time()), "\n")
cat("PROJ:", PROJ, "\n")
cat("atlas:", ATLAS, " mtime:", format(atlas_mtime), "\n")
cat("label cache:", HIP_CACHE, " mtime:", format(file.mtime(HIP_CACHE)), "\n")
cat("marker cache:", MARK_CACHE, " mtime:", format(file.mtime(MARK_CACHE)), "\n")
cat("cache reused (atlas-light):", fresh_cache, "\n\n")
cat("subtypes (x, n; 6 — EC-like NOT resolved on FH):\n"); print(table(factor(lab_vec, levels = SUBTYPE_ORDER)))
cat("\nmarker blocks (y):\n")
for (b in BLOCK_ORDER) cat(sprintf("  %-18s %s\n", b, paste(MARKER_BLOCKS[[b]], collapse = ", ")))
cat("\nmarkers requested (unique):", length(UNIQUE_MARKERS), "\n")
if (length(absent)) cat("markers ABSENT (dropped):", paste(absent, collapse = ", "), "\n") else
  cat("markers ABSENT: none — all present\n")
cat("\n--- sessionInfo() ---\n"); print(sessionInfo())
sink()
cat("wrote provenance ->", prov, "\n")
cat("=== DONE ===\n", file = stderr())
