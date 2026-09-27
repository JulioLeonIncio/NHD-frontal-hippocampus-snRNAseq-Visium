#!/usr/bin/env Rscript
# =============================================================================
# 3G_oligo_gene_panel_FH.R — Figure 3g: oligodendrocyte per-cell gene panel, frontal + hippocampus. Port of 09l_oligo_gene_heatmap_horizontal.R (atlas-light; reads the FH oligo MAST tables — no Seurat object needed).
# -----------------------------------------------------------------------------
# Conceptual REFRAME (load-bearing): the NHD myelin loss is compositional, not a
# per-cell collapse of myelin genes. This panel makes that explicit:
#   * the "Structural myelin" module reads flat-to-mildly-UP per cell in NHD
#     (PLP1/MOG/MAG/CNP/MOBP/OPALIN avg_log2FC >= 0 in both regions) — so the
#     panel is not allowed to imply a per-cell myelin collapse;
#   * the real per-cell program loss is in "Cholesterol / SREBP2" (all down);
#   * the NHD stress module (CRYAB/FTL/APOD/CLU) is up.
# The compositional loss (fewer high-myelin oligodendrocytes) is carried by the
# UMAP/density/state panels (3h/3i), not here. The module order + the legend note
# frame myelin as the per-cell-flat readout.
#
# rows (tracks) = Oligo:Frontal, Oligo:Hippocampus   (2 regions; oligo half only)
# cols (genes)  = curated modules (myelin / maturation / cholesterol / stress)
# fill          = avg_log2FC (NHD/CON) from per-nucleus MAST
# grey "q*"     = house effect-gated mark (harmonised with the sibling cholesterol
#                 dotplot / lineage dotmatrix): q(per-nucleus MAST padj) < 0.05 and
#                 |Cliff's delta| < 0.15 — flags a robust-but-small per-cell shift.
#
# reviewer-facing edge cases handled + logged:
#   * genes absent from a region table -> NA -> grey (not silently 0)
#   * padj==NA / padj==1 (not significant) -> no q* mark, value still shown
#   * padj==0 underflow -> floored to xmin for the gate only (never reported as 0)
#   * detection-floor greying (an avg_log2FC on a near-absent gene is noise)
#   * colour scale clipped to a stated robust quantile (printed to log)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ComplexHeatmap); library(circlize); library(grid)
})

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # PAL_REGION_PALE, REGION_FULL, REGION_ORDER
PANEL <- file.path(PROJ, "figures", "Figure_4", "panels")   # oligo -> Figure 4
LOGD  <- file.path(PROJ, "logs")
for (d in c(PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- read FH oligo MAST tables (one per region) -----------------------------
REGIONS <- REGION_ORDER                        # c("Frontal","Hippo")
# OPC added as two further tracks. The oligodendrocyte rows carry
# the myelin/cholesterol story; the OPC rows are the lineage-internal control showing
# which of those shifts belong to the myelinating compartment specifically.
LINEAGES <- c("Oligo", "OPC")
mast <- bind_rows(lapply(LINEAGES, function(ct) {
  bind_rows(lapply(REGIONS, function(r) {
    f <- file.path(PROJ, "tables", "mast_dual", sprintf("MAST_%s_%s.csv", ct, r))
    stopifnot("MISSING MAST table — prep step failed" = file.exists(f))
    d <- read.csv(f); d$Region <- r; d$celltype <- ct; d
  }))
}))
cat(sprintf("MAST rows: %d (Oligo %d, OPC %d)\n", nrow(mast),
            sum(mast$celltype == "Oligo"), sum(mast$celltype == "OPC")))

# ---- curated gene modules --------
# Left-to-right developmental / functional order. Myelin first so the panel opens
# on the "per-cell myelin is flat/up" readout (the reframe). Only oligo-identity
# genes above the detection floor in >=1 region; no neuronal-soup / lncRNA.
# "Cholesterol/SREBP2", "Stress / reactive") overprinted into unreadable garble
# above the narrow column blocks. Shortened to one word each so each strip title
# fits its block ("Myelin | Maturation | Cholesterol | Stress"); the full module
# meaning (SREBP2 / reactive) moves to the legend, per house no-in-panel-caption.
# Blocks and their membership are kept in lockstep with
# scripts/_oligo_programmes_FH.R (REBUILD_PLAN entries 33-34).
#
# Changes from the previous six-block version:
#   * OPALIN kept in Myelin; TF and ANLN split OUT into their own `Support` block —
#     they are the two genes that fall while the membrane proteins rise (TF -0.44/-0.55,
#     ANLN -0.37/-0.46), so averaging them into "myelin" cancelled the claim. ANLN also
#     replicates Zhou 2023's ANLN-down in all three of their NHD donors.
#   * the old single `Lipid` block is split into its ARMS: Sphingolipid (acyl-chain),
#     Galactolipid (head-group, the selectivity control), Lipid import (can the cell
#     scavenge what it cannot make?), Sterol export (the LXR discriminator against the
#     sterol-repletion alternative).
#     They now say what the genes do:
#     TF is iron delivery and ANLN the paranodal septin scaffold, so `Iron / scaffold`;
#     SLC44A1/NPC1/LDLR bring lipid precursors in, so `Lipid import`; ABCG1/NPC2 move
#     sterol OUT under LXR, so `Sterol export`. Two lines each so the wider words do not
#     overprint their narrow blocks.
#   * `Myelin` -> `Myelin proteins` and `Sterol` -> `Sterol synthesis`: the panel's whole
#     argument is protein-UP against lipid-DOWN, so the first block must say PROTEINS;
#     and once import/export blocks exist, "Sterol" alone does not say which of the three
#     it is. The lipid blocks now read synthesis | import | export in parallel.
#   * `Maturation` retired to the programme panel (differentiation is a row there now,
#     with BCAS1 added) to keep this panel placeable at ~7.2 in.
# One-word block titles: multi-word titles overprint above narrow blocks;
# the full meaning lives in the legend, per the no-in-panel-caption rule.
gene_groups <- list(
  `Myelin\nproteins`  = c("PLP1","MOBP","MOG","MAG","CNP","CLDN11","OPALIN"),
  `Iron /\nscaffold` = c("TF","ANLN"),
  `Sterol\nsynthesis`  = c("SREBF2","HMGCS1","FDFT1","SQLE","MSMO1","DHCR24"),
  `Sphingo-\nlipid` = c("SCD","SGMS1","ACER3","ELOVL5"),
  `Galacto-\nlipid` = c("UGT8","FA2H","GALC"),
  `Lipid\nimport`  = c("SLC44A1","NPC1","LDLR"),
  `Sterol\nexport`  = c("ABCG1","SREBF1","NPC2"),
  `Stress`  = c("CRYAB","HSPA1A","FTL","APOD","CLU"),
  `OPC`     = c("PDGFRA","LHFPL3","SOX6","BCAN","MEGF11","GPR17")
)
gene_order <- unlist(gene_groups, use.names = FALSE)
gene_split <- factor(rep(names(gene_groups), lengths(gene_groups)),
                     levels = names(gene_groups))

# ---- build matrices (rows = region tracks, cols = genes) --------------------
track_levels <- as.vector(t(outer(LINEAGES, REGIONS, paste, sep = "|")))  # 4 tracks
track_ct <- sub("\\|.*$", "", track_levels)
track_rg <- sub("^.*\\|", "", track_levels)
lfc_long <- mast %>%
  filter(gene %in% gene_order) %>%
  transmute(gene, track = paste(celltype, Region, sep = "|"),
            log2FoldChange = avg_log2FC, padj = p_val_adj,
            delta = cliffs_delta,                       # Cliff's delta for house q* gate
            det = pmax(pct.1, pct.2, na.rm = TRUE))

missing_genes <- setdiff(gene_order, unique(lfc_long$gene))
if (length(missing_genes))
  cat("Genes absent from ALL region tables (dropped):",
      paste(missing_genes, collapse = ", "), "\n")
gene_present <- intersect(gene_order, unique(lfc_long$gene))
gene_split_p <- gene_split[gene_order %in% gene_present]

to_mat <- function(col) {
  m <- lfc_long %>% select(gene, track, val = all_of(col)) %>%
    pivot_wider(names_from = gene, values_from = val, values_fill = NA_real_) %>%
    column_to_rownames("track") %>% as.matrix()
  # Ensure all tracks x genes present (a gene not tested in a region -> NA row cell)
  full <- matrix(NA_real_, length(track_levels), length(gene_present),
                 dimnames = list(track_levels, gene_present))
  full[rownames(m), colnames(m)] <- m[, , drop = FALSE]
  full[track_levels, gene_present, drop = FALSE]
}
mat_lfc   <- to_mat("log2FoldChange")
mat_padj  <- to_mat("padj")
mat_delta <- to_mat("delta")
mat_det   <- to_mat("det")

# ---- detection floor, measured on raw COUNTS -------------------
# Until now `mat_det` came from the MAST tables and a cell was greyed when the gene was
# absent from them. `40_mast_percell_dual_FH.R` runs FindMarkers(min.pct = 0.1,
# logfc.threshold = 0.1), so a gene is dropped for failing the detection floor or for being
# stable — and 11 of the 61 greys were genes detected in >= 10 % of nuclei, MOBP at 72 % among
# them. Greying them erased the galactolipid block of the frontal oligodendrocyte track, which
# is the preserved arm that makes the sterol/sphingolipid failure specific rather than global.
# Detection is therefore measured directly, and the two states are now drawn differently:
#   grey       = below the 10 % detection floor (nothing can be said)
#   open tile  = detected and stable (|log2FC| < 0.1, so MAST returns no row) -> a real near-zero
DET_FLOOR  <- 0.10
source(file.path(PROJ, "scripts", "_receptor_detection_FH.R"))
.det_raw <- receptor_detection(PROJ, file.path(PROJ, "atlas", "NHD_FH_harmony.rds"),
                               genes = gene_present, classes = LINEAGES, regions = REGIONS)
mat_det_raw <- matrix(NA_real_, length(track_levels), length(gene_present),
                      dimnames = list(track_levels, gene_present))
for (tr in track_levels) {
  ct <- sub("\\|.*$", "", tr); rg <- sub("^.*\\|", "", tr)
  i <- match(paste(gene_present, ct, rg), paste(.det_raw$gene, .det_raw$class, .det_raw$region))
  mat_det_raw[tr, ] <- .det_raw$pct_max[i]
}
mat_det <- ifelse(is.na(mat_det), mat_det_raw, pmax(mat_det, mat_det_raw, na.rm = TRUE))
floor_mask  <- is.na(mat_det) | (mat_det < DET_FLOOR)                 # nothing measurable
stable_mask <- !floor_mask & is.na(mat_lfc)                           # detected, |log2FC| < 0.1
n_nottested <- sum(is.na(mat_lfc))
n_floored   <- sum(floor_mask & !is.na(mat_lfc))
cat(sprintf("detection on raw counts: %d cell(s) below the %.0f %% floor -> grey; %d detected but stable -> open tile\n",
            sum(floor_mask), 100 * DET_FLOOR, sum(stable_mask)))
if (sum(stable_mask)) cat("  detected-but-stable:", paste(sprintf("%s/%s %.0f%%",
    col(stable_mask)[stable_mask] |> (\(j) colnames(mat_lfc)[j])(),
    row(stable_mask)[stable_mask] |> (\(i) rownames(mat_lfc)[i])(),
    100 * mat_det[stable_mask]), collapse = ", "), "\n")
# a detected-but-stable cell is a bounded near-zero, not missing data: MAST withheld it because
# |log2FC| < 0.1, so it is drawn at zero (white) inside a dashed outline that says "bounded, not
# estimated". Only cells below the detection floor stay grey.
mat_lfc[stable_mask]  <- 0
mat_lfc[floor_mask]   <- NA
mat_padj[floor_mask]  <- NA
mat_delta[floor_mask] <- NA

# ---- house effect-gated "q*" mark (harmonised with sibling oligo dotplots) ----
# q* is drawn where q(MAST padj) < Q_CUT and |Cliff's delta| < CD_CUT: a robust but
# small per-cell shift. Matches 3_oligo_cholesterol_FH.R exactly (Q_CUT/CD_CUT/XMIN).
Q_CUT  <- 0.05
CD_CUT <- 0.15
XMIN   <- .Machine$double.xmin          # underflow floor for padj==0 (see below)
# reviewer-facing p=0 underflow: a MAST padj printed as exactly 0 is floating-point
# underflow, not a true zero. Floor to XMIN for the gate only (never reported as 0),
# and LOG how many tiles were affected so the handling is auditable.
n_uf <- sum(!is.na(mat_padj) & mat_padj > 0 & mat_padj < XMIN)
n_z  <- sum(!is.na(mat_padj) & mat_padj == 0)
if (n_uf + n_z > 0)
  cat(sprintf("NOTE: %d padj==0 + %d in (0,xmin) underflowed; floored to xmin=%.2e for the q* gate only.\n",
              n_z, n_uf, XMIN))
mat_qgate <- (!is.na(mat_padj) & pmax(mat_padj, XMIN) < Q_CUT) &
             (!is.na(mat_delta) & abs(mat_delta) < CD_CUT)
cat(sprintf("q* mark: %d tiles pass q<%.2f & |delta|<%.2f (effect-gated, NOT per-nucleus stars)\n",
            sum(mat_qgate), Q_CUT, CD_CUT))

# ---- per-cell myelin flatness assertion (reframe guard, logged) -------------
# The panel must not imply per-cell myelin collapse. Assert + log that the
# structural-myelin module is flat-to-UP per cell (median avg_log2FC >= 0) in
# both regions; if this ever flips negative the reframe/caption must be revisited.
my_g <- intersect(gene_groups[["Myelin"]], colnames(mat_lfc))
for (tr in track_levels[track_ct == "Oligo"]) {   # the reframe guard is about oligodendrocytes
  v <- mat_lfc[tr, my_g]; med <- median(v, na.rm = TRUE)
  cat(sprintf("  per-cell myelin module %s: median avg_log2FC = %+.2f (n=%d genes) -> %s\n",
              tr, med, sum(!is.na(v)), ifelse(med >= 0, "FLAT/UP (compositional-loss framing OK)",
                                              "DOWN (REVISIT reframe!)")))
}

# ---- diverging colour scale (clipped to stated robust quantile) -------------
lim <- quantile(abs(mat_lfc), 0.90, na.rm = TRUE); lim <- max(lim, 1)
cat(sprintf("fill scale: symmetric +/- %.2f (90th pct |log2FC|; values clipped for display)\n", lim))
col_fun <- colorRamp2(c(-lim, 0, lim), c("#3E7CB1","white","#D1495B"))

# ---- row annotation (region banners, pale + full names) ---------------------
row_labs <- paste0(track_ct, " : ", unname(REGION_FULL[track_rg]))
RG_COL   <- setNames(unname(PAL_REGION_PALE[REGIONS]), REGIONS)
row_anno <- rowAnnotation(
  Region = track_rg,
  col = list(Region = RG_COL),
  # Hide the "Region" annotation name: ComplexHeatmap draws it at the bottom of the
  # colour track, aligned with the x gene labels, where it overprinted the first
  # gene name (PLP1) in the bottom-left corner. The colour track's meaning is carried
  # by its own legend, so the axis-level name is redundant.
  show_annotation_name = FALSE,
  simple_anno_size = unit(3.2, "mm"),
  # Suppress the redundant "Region" colour legend: region is already unambiguous
  # from the row labels ("Oligo : Frontal" / "Oligo : Hippocampus") and the row
  # colour bar, so the separate swatch key was a third region cue. Keeps the panel
  # key to just the log2FC scale.
  show_legend = FALSE)

ht <- Heatmap(
  mat_lfc, name = "log2 FC\n(NHD/CON)", col = col_fun, na_col = "grey85",
  cluster_rows = FALSE, cluster_columns = FALSE,
  column_split = gene_split_p,
  # visual gap between the oligodendrocyte block and the OPC control block
  row_split = factor(track_ct, levels = LINEAGES), row_title = NULL,
  row_gap = unit(1.4, "mm"),
  column_title_gp = gpar(fontsize = 8, fontface = "plain"), column_title_side = "top",
  row_names_side = "left", row_labels = row_labs,
  row_names_gp = gpar(fontsize = 8.3),
  column_names_gp = gpar(fontsize = 8.3, fontface = "italic"), column_names_rot = 60,
  rect_gp = gpar(col = "grey92", lwd = 0.35), left_annotation = row_anno, border = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 7, fontface = "plain"), labels_gp = gpar(fontsize = 6.5),
    legend_height = unit(2.2, "cm"), grid_width = unit(3, "mm")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    # House effect-gated grey "q*" (q<0.05 & |Cliff's delta|<0.15), identical to the
    # sibling cholesterol dotplot / lineage dotmatrix, so this panel reads as the same
    # effect-size idiom — not a per-nucleus significance-stars panel. Small grey text, drawn on top.
    if (isTRUE(mat_qgate[i, j]))
      grid.text("q*", x, y, gp = gpar(fontsize = 6.5, col = "grey45"), vjust = 0.72)
    # detected and stable: an outlined empty tile, never the grey of "not measurable"
    if (isTRUE(stable_mask[i, j]))
      grid.rect(x, y, w * 0.86, h * 0.86, gp = gpar(fill = NA, col = "grey55", lty = "22", lwd = 0.6))
  })

n_gene <- ncol(mat_lfc)
# 10% wider and 25% shorter than the 7.42 x 2.90 version.
# Type sizes are not scaled — rendering at the placed size is what keeps a pt a pt.
W <- max(6.0, 0.155 * n_gene + 2.2) * 0.90 * 1.10
H <- 2.9 * 0.75                                # 4 tracks: (Oligo + OPC) x 2 regions

# House RULE: no in-panel caption. The grey "q*" meaning (q<0.05 & |Cliff's delta|
# <0.15, effect-gated) and the single-NHD-Hippocampus-lane caveat live in the figure
# legend, not on the plot. Only the grey "q*" cell marks remain.
BN <- "F4c_oligo_gene_panel"
pdf(file.path(PANEL, paste0(BN, ".pdf")), width = W, height = H, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE, padding = unit(c(3,3,3,3), "mm"))
invisible(dev.off())
png(file.path(PANEL, paste0(BN, ".png")), width = W, height = H, units = "in", res = 600)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     merge_legend = TRUE, padding = unit(c(3,3,3,3), "mm"))
invisible(dev.off())
cat(sprintf("Saved %s.{pdf,png}  (%d genes, %d tracks, %.2f x %.2f in)\n",
            BN, n_gene, nrow(mat_lfc), W, H))

writeLines(capture.output(sessionInfo()),
           file.path(LOGD, "3G_oligo_gene_panel_FH_sessionInfo.txt"))
cat("=== DONE ===\n", file = stderr())
