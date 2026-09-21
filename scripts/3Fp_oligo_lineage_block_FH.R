#!/usr/bin/env Rscript
# =============================================================================
# 3Fp_oligo_lineage_block_FH.R — Figure 3f: oligodendrocyte maturation- continuum signature block (NHD-vs-CON), frontal + hippocampus.
# Port of 09b_v2_lineage_block_final.R + 62_fig3_oligo_lineage_dotmatrix.R,
# REFRAMED for the continuum (oligodendrocytes are one continuum, not subtypes).
# Reads the shared oligo cache built by 3F_oligo_embed_states_FH.R (atlas-light)
# and the OPC cache built by 3M_prep_opc_and_scales_FH.R.
#
# OPC added as a second lineage column. OPC is a separate
# population from the oligodendrocyte continuum -- that is exactly why it belongs
# here: it is the internal control that says which programme shifts are specific
# to the myelinating compartment and which run through the whole lineage. Both
# views now show Oligodendrocyte and OPC side by side within each region.
# -----------------------------------------------------------------------------
# Renders two co-registered views of the same per-nucleus signature scores:
#   3f_oligo_lineage_block      : CON-vs-NHD violin grid (signature x region)
#   F4d_oligo_lineage_dotmatrix  : compact Cliff's-delta dot-matrix (same numbers)
#
# Signatures (oligo maturation continuum + the compositional myelin story):
#   Myelination           — canonical structural myelin (per-cell ~flat/UP: reframe)
#   Oligo differentiation — OPC->OL fate TFs
#   Cholesterol synthesis — SREBP2 sterol cascade (the real per-cell down program)
#   Myelin lipid synthesis— ceramide/galactolipid myelin-lipid enzymes
#   Stress / reactive     — heat-shock/ferritin NHD stress module (up)
#
# House rules honored:
#   * Effect size (Cliff's delta, NHD-CON) + gated house "q*" mark — no per-nucleus
#     significance stars (single-donor pseudoreplication). The gated rule: */**/***
#     only when q<alpha and |delta|>=0.15; grey "q*" when q<alpha and |delta|<0.15.
#   * Region strips = pale banners + black text + full names; single-lane-Hippo
#     honesty -> grey Hippo strip + "(low n)" (Hippo NHD = one lane, 1345 nuclei).
#   * no bold, no captions in panels.
#   * p/q=0 underflow floored + -log10 capped, both logged + stated in legend.
# -----------------------------------------------------------------------------
# apparent-type lift 3Fp_oligo_lineage_block_FH.R: every text size in this script scaled
# by 1.12 so the panel reads at ~5.0 pt on the assembled page, matching the
# Figure-5 dot panels. Point sizes, line widths and unit() dimensions are not
# touched. Canvas size unchanged, so re-linking is a no-op.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(scales); library(ggh4x); library(stringr)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # PAL_REGION, REGION_FULL, REGION_ORDER, strip_region_x, region_lowconf_labeller
DDIR  <- file.path(PROJ, "data")
PANEL <- file.path(PROJ, "figures", "Figure_4", "panels")   # oligo → Fig4
LOGD  <- file.path(PROJ, "logs")
for (d in c(PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
CACHE <- file.path(DDIR, "_cache_oligo_umap.rds")
stopifnot("MISSING atlas" = file.exists(ATLAS),
          "MISSING oligo cache — run 3F_oligo_embed_states_FH.R first" = file.exists(CACHE))
if (file.mtime(CACHE) < file.mtime(ATLAS))
  stop("STALE oligo cache: ", basename(CACHE), " predates the atlas. Re-run 3F_oligo_embed_states_FH.R.")

REGIONS <- REGION_ORDER
COND_COL <- c(CON = "#82B2D6", NHD = "#D6837A")   # Fig 2c soft blue/red

# ── single-lane-Hippocampus honesty ──────────────────────────────────────────
# Hippo NHD rests on a single lane (TFHS000535, 1345 oligo nuclei). Flag the Hippo
# facet low-confidence (grey banner + "(low n)") — the per-nucleus n is large but
# it is one donor lane, so region-level inference is weaker there.
HIPPO_NHD_LANES <- 1L
low_conf <- if (HIPPO_NHD_LANES < 2) "Hippo" else character(0)

# ── signature gene sets (from curated_signatures.rds + literature) ───────────
sigs <- readRDS(file.path(DDIR, "curated_signatures.rds"))
# Gene sets live in one file so the Visium reconciliation panel scores exactly the same
# programmes.
source(file.path(PROJ, "scripts", "_oligo_programmes_FH.R"))
gene_sets <- oligo_programme_sets(sigs)[oligo_programme_order()]

# ── load both lineage caches, score signatures per nucleus ───────────────────
# Identical scoring on both objects (RNA LogNormalize -> mean expression across
# the detected genes of each set) so Oligodendrocyte and OPC are directly
# comparable within a facet.
OPC_CACHE <- file.path(DDIR, "_cache_opc_FH.rds")
stopifnot("MISSING OPC cache — run 3M_prep_opc_and_scales_FH.R first" = file.exists(OPC_CACHE))
if (file.mtime(OPC_CACHE) < file.mtime(ATLAS))
  stop("STALE OPC cache: ", basename(OPC_CACHE), " predates the atlas. Re-run 3M_prep_opc_and_scales_FH.R.")

LINEAGE_LEVELS <- c("Oligodendrocyte", "OPC")
load_norm <- function(path) {
  x <- readRDS(path)
  DefaultAssay(x) <- "RNA"
  if (length(SeuratObject::Layers(x, assay = "RNA")) > 1) x <- SeuratObject::JoinLayers(x)
  NormalizeData(x, verbose = FALSE)
}
o   <- load_norm(CACHE)
opc <- load_norm(OPC_CACHE)
data_mat <- GetAssayData(o, assay = "RNA", layer = "data")

feats <- lapply(gene_sets, function(g) intersect(g, rownames(data_mat)))
for (nm in names(feats))
  cat(sprintf("signature %-24s %d/%d genes: %s\n", nm, length(feats[[nm]]),
              length(gene_sets[[nm]]), paste(feats[[nm]], collapse=",")))
feats <- feats[sapply(feats, length) >= 3]
stopifnot("no signatures with >=3 genes" = length(feats) > 0)

score_one <- function(x, lineage_label) {
  dm <- GetAssayData(x, assay = "RNA", layer = "data")
  miss <- setdiff(unlist(feats), rownames(dm))
  if (length(miss))
    cat(sprintf("  [%s] %d signature genes absent from this object: %s\n",
                lineage_label, length(miss), paste(miss, collapse = ",")))
  sm <- sapply(feats, function(g) Matrix::colMeans(dm[intersect(g, rownames(dm)), , drop = FALSE]))
  cbind(x@meta.data[, c("Region","Condition")], lineage = lineage_label,
        as.data.frame(sm))
}
score_df <- bind_rows(score_one(o, "Oligodendrocyte"), score_one(opc, "OPC")) %>%
  pivot_longer(cols = all_of(names(feats)), names_to = "signature", values_to = "score") %>%
  mutate(Condition = factor(Condition, levels = c("CON","NHD")),
         Region    = factor(Region,    levels = REGIONS),
         lineage   = factor(lineage,   levels = LINEAGE_LEVELS),
         signature = factor(signature, levels = names(feats)))
cat(sprintf("scored %d nuclei (Oligodendrocyte %d, OPC %d)\n",
            ncol(o) + ncol(opc), ncol(o), ncol(opc)))

# ── Cliff's delta (NHD - CON) per signature x region + BH-q gated label ──────
cliff_delta <- function(x, y) {   # delta = 2*AUC-1 ; +delta = higher in x (NHD)
  x <- x[!is.na(x)]; y <- y[!is.na(y)]; nx <- length(x); ny <- length(y)
  if (nx < 3 || ny < 3) return(NA_real_)
  r <- rank(c(x, y)); U <- sum(r[seq_len(nx)]) - nx*(nx+1)/2
  (2*U/(nx*ny)) - 1
}
stat <- score_df %>%
  group_by(signature, Region, lineage) %>%
  summarise(
    p       = tryCatch(wilcox.test(score[Condition=="NHD"], score[Condition=="CON"],
                                   exact = FALSE)$p.value, error = function(e) NA_real_),
    cliff_d = cliff_delta(score[Condition=="NHD"], score[Condition=="CON"]),
    y_min = min(score, na.rm=TRUE), y_max = max(score, na.rm=TRUE),
    n_CON = sum(Condition=="CON"), n_NHD = sum(Condition=="NHD"),
    .groups = "drop") %>%
  mutate(q_BH   = p.adjust(p, method = "BH"),   # BH across all signature x region x lineage cells
         passes = !is.na(cliff_d) & abs(cliff_d) >= 0.15,
         # Drop the per-nucleus
         # ***/**/* ladder — single-donor pseudoreplication cannot support a
         # donor-level significance claim. Canonical convention (matches
         # 2B_micro_program_violins_FH.R): Cliff's delta label + only the grey
         # effect-gated "q*" mark (q<0.05 and |delta|<0.15). q_BH stays in the CSV.
         label  = case_when(q_BH < 0.05 & !passes ~ "q*", TRUE ~ ""),
         delta_lab = sprintf("italic(delta) == '%+.2f'", cliff_d))

# ── degenerate-value guards (reviewer-facing) — p/q underflow + -log10 cap ────
XMIN <- .Machine$double.xmin
CAPL <- 50                          # ceiling on -log10(q): q = 1e-50 (stated in legend)
n_zero  <- sum(!is.na(stat$q_BH) & stat$q_BH == 0)
n_under <- sum(!is.na(stat$q_BH) & stat$q_BH > 0 & stat$q_BH < XMIN)
raw_mlp <- ifelse(is.na(stat$q_BH), NA_real_, -log10(pmax(stat$q_BH, XMIN)))
n_cap   <- sum(!is.na(raw_mlp) & raw_mlp > CAPL)
cat(sprintf("size guard: %d q==0 + %d in (0,xmin) floored to xmin (-log10=%.1f); %d capped at %d (q=1e-%d)\n",
            n_zero, n_under, -log10(XMIN), n_cap, CAPL, CAPL))
stat$dotsize <- pmin(raw_mlp, CAPL)

# ── depth-matched control delta (not plotted; written to the CSV) ────────────
# Every delta above is computed on LogNormalize scores, i.e. on each gene's share of its
# nucleus. NHD frontal oligodendrocytes carry ~half the transcriptome of a control
# oligodendrocyte (median 1520 -> 799 UMI), so a share can move without the quantity
# moving (this is the whole subject of 4_myelin_scales_bridge). The control below thins
# every nucleus to an identical UMI depth and recomputes the same delta on raw counts;
# a signature whose two deltas agree is not reading the contraction. Both columns ship in
# the CSV so a legend or a reviewer reply can quote either, explicitly labelled.
DEPTH_T <- 750L
depth_matched_delta <- function(obj, lineage_label) {
  cts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  m   <- obj@meta.data
  tot <- Matrix::colSums(cts)
  out <- list()
  for (rg in REGIONS) {
    keep <- m$Region == rg & tot >= DEPTH_T
    if (sum(keep) < 50) next
    ds <- Seurat::SampleUMI(cts[, keep, drop = FALSE], max.umi = DEPTH_T,
                            upsample = FALSE, verbose = FALSE)
    rownames(ds) <- rownames(cts)
    mm <- m[keep, , drop = FALSE]
    cat(sprintf("  depth-match %s %s: kept %d/%d nuclei (CON %d / NHD %d)\n",
                lineage_label, rg, sum(keep), sum(m$Region == rg),
                sum(mm$Condition == "CON"), sum(mm$Condition == "NHD")))
    for (nm in names(feats)) {
      g <- intersect(feats[[nm]], rownames(ds))
      v <- Matrix::colSums(ds[g, , drop = FALSE])
      out[[length(out) + 1]] <- data.frame(
        signature = nm, Region = rg, lineage = lineage_label,
        cliff_d_depth_matched = cliff_delta(v[mm$Condition == "NHD"],
                                            v[mm$Condition == "CON"]),
        ratio_depth_matched   = mean(v[mm$Condition == "NHD"]) /
                                mean(v[mm$Condition == "CON"]),
        n_depth_matched       = sum(keep))
    }
  }
  bind_rows(out)
}
cat(sprintf("\ndepth-matched control at T = %d UMI:\n", DEPTH_T))
dmc <- bind_rows(depth_matched_delta(o, "Oligodendrocyte"),
                 depth_matched_delta(opc, "OPC"))
stat <- stat %>% left_join(dmc, by = c("signature","Region","lineage"))
cat("\nplotted (share-based) vs depth-matched delta:\n")
print(as.data.frame(stat %>%
  mutate(cliff_d = round(cliff_d, 3),
         cliff_d_depth_matched = round(cliff_d_depth_matched, 3),
         ratio_depth_matched = round(ratio_depth_matched, 2)) %>%
  select(signature, Region, lineage, cliff_d, cliff_d_depth_matched, ratio_depth_matched)),
  row.names = FALSE)

# ── write the delta/q/n table (legend numbers key to this csv) ───────────────
TDIR <- file.path(PROJ, "tables"); dir.create(TDIR, showWarnings = FALSE, recursive = TRUE)
stat_export <- stat %>% select(signature, Region, lineage, p, q_BH, cliff_d,
                               cliff_d_depth_matched, ratio_depth_matched, n_depth_matched,
                               passes, label, n_CON, n_NHD)
write.csv(stat_export, file.path(TDIR, "oligo_lineage_deltas_FH.csv"), row.names = FALSE)
cat("wrote tables/oligo_lineage_deltas_FH.csv (", nrow(stat_export), " rows)\n", sep = "")
cat("\nper (signature x region) delta / q / label:\n")
print(as.data.frame(stat %>% mutate(cliff_d = round(cliff_d,3), q_BH = signif(q_BH,3)) %>%
                    select(signature, Region, lineage, cliff_d, q_BH, label, n_CON, n_NHD)), row.names = FALSE)

# =============================================================================
# View 1 — CON-vs-NHD violin grid (signature rows x region cols)
# =============================================================================
sig_ann <- stat %>%
  mutate(y_top   = y_max + 0.06*(y_max - y_min),
         y_delta = y_max + 0.19*(y_max - y_min))
strip_lo <- strip_region_x(REGIONS, fontsize = 9.0, clip = "off")  # Hippo = normal pale tone; caveat in legend
DODGE <- position_dodge(width = 0.78)
p_viol <- ggplot(score_df, aes(x = lineage, y = score, fill = Condition)) +
  geom_violin(scale = "width", linewidth = 0.25, alpha = 0.85, trim = FALSE,
              width = 0.68, position = DODGE) +
  # Fill = "white" is set outside aes(), which removes the inherited
  # fill = Condition mapping -- without an explicit group the layer would dodge by
  # `lineage` alone and draw one box per lineage over the pooled CON+NHD nuclei
  # (a single box sat between the two violins). The
  # explicit group restores one box per violin.
  geom_boxplot(aes(group = interaction(lineage, Condition)),
               width = 0.13, outlier.shape = NA, fill = "white",
               colour = "grey25", linewidth = 0.3, fatten = 1.8, position = DODGE) +
  # Only the grey effect-gated "q*" mark (no per-nucleus */**/***, no n.s.).
  geom_text(data = subset(sig_ann, label == "q*"),
            aes(x = lineage, y = y_top, label = label), inherit.aes = FALSE,
            size = 2.2, colour = "grey50", vjust = 0.4) +
  geom_text(data = sig_ann, aes(x = lineage, y = y_delta, label = delta_lab),
            inherit.aes = FALSE, parse = TRUE, hjust = 0.5, vjust = 0, size = 2.0, colour = "black") +
  facet_grid2(signature ~ Region, scales = "free_y", switch = "y",
              labeller = labeller(Region = region_full_lowconf_labeller(low_conf)),
              strip = strip_lo) +
  scale_fill_manual(values = COND_COL, name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.42)) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), expand = expansion(mult = c(0.04, 0.30))) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8.4) +
  theme(strip.clip = "off",
        strip.text.y.left = element_text(size = 9.0, face = "plain", colour = "black", angle = 0, hjust = 1),
        strip.background.y = element_rect(fill = NA, colour = NA), strip.placement = "outside",
        axis.text.x = element_text(size = 7.6, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7.3), axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x = unit(0.3, "lines"), panel.spacing.y = unit(0.35, "lines"),
        legend.position = "bottom", legend.key.size = unit(0.32, "cm"),
        legend.text = element_text(size = 9.0), legend.margin = margin(0,0,0,0),
        plot.margin = margin(6, 6, 4, 4))
# The violin grid is the supplement sibling of the dot-matrix. Write it straight into
# panels/_supp/ so a re-render cannot leave a duplicate in the main-panel folder.
SUPP <- file.path(PANEL, "_supp"); dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)
BN1 <- "3f_oligo_lineage_block"
ggsave(file.path(SUPP, paste0(BN1, ".pdf")), p_viol, width = 4.9, height = 6.4, useDingbats = FALSE)
ggsave(file.path(SUPP, paste0(BN1, ".png")), p_viol, width = 4.9, height = 6.4, dpi = 600, device = ragg::agg_png)
cat("Saved _supp/", BN1, ".{pdf,png}\n", sep = "")

# =============================================================================
# View 2 — compact Cliff's-delta dot-matrix (same numbers as the violins)
# =============================================================================
L <- max(ceiling(max(abs(stat$cliff_d), na.rm = TRUE) * 20) / 20, 0.15)
cat(sprintf("dotmatrix fill: symmetric +/- %.2f (max |delta|=%.3f)\n", L, max(abs(stat$cliff_d), na.rm=TRUE)))
dm <- stat %>% mutate(delta_c = pmax(pmin(cliff_d, L), -L),
                      Region = factor(Region, levels = REGIONS),
                      lineage = factor(lineage, levels = LINEAGE_LEVELS),
                      signature = factor(signature, levels = names(feats)))
GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
# Only the grey effect-gated "q*" mark on the dot-matrix (no */**/***).
stars_qgrey <- dm %>% filter(label == "q*")
# Narrow 2-column dot-matrix: shrink strip fontsize + clip="off" so the full
# "Hippocampus" banner is not truncated.
strip_lo2 <- strip_region_x(REGIONS, fontsize = 7.3, clip = "off")
ylab_map <- setNames(str_wrap(names(feats), width = 18), names(feats))

# Full region name only (no suffix); Frontal / Hippocampus both single-line.
dm_region_labeller <- ggplot2::as_labeller(function(x)
  ifelse(x %in% names(REGION_FULL), unname(REGION_FULL[x]), x))

p_dm <- ggplot(dm, aes(x = lineage, y = signature)) +
  geom_point(aes(fill = delta_c, size = dotsize), shape = 21, colour = "grey35", stroke = 0.3) +
  geom_text(data = stars_qgrey, aes(x = lineage, y = signature, label = label),
            inherit.aes = FALSE, nudge_y = 0.30, size = 2.6, colour = "grey50", vjust = 0.5) +
  facet_wrap2(~ Region, nrow = 1,
              labeller = dm_region_labeller, strip = strip_lo2) +
  scale_fill_gradient2(low = GSEA_STEEL, mid = "white", high = GSEA_ROSE, midpoint = 0,
                       limits = c(-L, L), breaks = c(-L, 0, L),
                       labels = c(sprintf("%.2f\n(CON-up)", -L), "0", sprintf("+%.2f\n(NHD-up)", L)),
                       name = expression(atop("Cliff's " * delta, "(NHD - CON)"))) +
  scale_size_continuous(range = c(0.8, 5.2), limits = c(0, CAPL), breaks = c(1.3, 10, 50),
                        labels = c("0.05", "1e-10", "<=1e-50"), name = "BH q") +
  scale_x_discrete(labels = c("Oligodendrocyte" = "Oligo", "OPC" = "OPC"),
                   expand = expansion(add = 0.55)) +
  scale_y_discrete(limits = rev(names(feats)), labels = function(x) ylab_map[x]) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(0.30,"cm"), barheight = unit(1.7,"cm")),
         size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 10.1) +
  # Pin strip.text so both region banners render at the intended 6.5 pt. ggh4x
  # elem_list_text() styles only the first strip; without this override the second
  # facet ("Hippocampus") inherits the theme_pub base_size (9 pt) fallback, making
  # the sibling banners mismatched in size.
  # Axis.text 6.3, panel.spacing.x 1.5 pt, legend text and
  # title 6.2, legend.box.just left. The y labels were carrying 8.5 pt against siblings at
  # 6.3, which is what made this panel read loose next to Fig 3's.
  theme(strip.text = element_text(size = 7.3, face = "plain"),
        axis.text.x = element_text(size = 7.1, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 7.1, colour = "black", lineheight = 0.82),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x = unit(1.5, "pt"), strip.clip = "off",
        legend.position = "right", legend.box = "vertical", legend.box.just = "left",
        legend.margin = margin(t = -2, r = 0, b = 0, l = 0),
        legend.key.size = unit(0.4, "lines"), legend.text = element_text(size = 6.9),
        legend.title = element_text(size = 6.9, lineheight = 0.9), plot.margin = margin(5, 4, 3, 3))
BN2 <- "F4d_oligo_lineage_dotmatrix"
# Width 3.4 -> 4.0 in so the two full-name banners (esp. "Hippocampus
# (low n)") are not clipped even with the strip fontsize reduced.
ggsave(file.path(PANEL, paste0(BN2, ".pdf")), p_dm, width = 2.91, height = 2.21, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN2, ".png")), p_dm, width = 2.91, height = 2.21, dpi = 600, device = ragg::agg_png)
cat("Saved ", BN2, ".{pdf,png}\n", sep = "")

writeLines(capture.output(sessionInfo()), file.path(LOGD, "3Fp_oligo_lineage_block_FH_sessionInfo.txt"))
cat("=== DONE ===\n", file = stderr())
