#!/usr/bin/env Rscript
# =============================================================================
# 3I_oligo_state_dotplot_FH.R — Figure 3i: marker dotplot for the oligo maturation-continuum positions + the 2 NHD states. frontal + hippocampus.
# Port of 09h_oligo_final_labeled.R's dotplot (block-grouped, house grammar).
# Reads the shared cache built by 3F_oligo_embed_states_FH.R (atlas-light).
# -----------------------------------------------------------------------------
# Conceptual reframe: the columns are POSITIONS along one continuum + the 2 NHD
# states, not invented discrete subtypes. Marker blocks:
#   OPC-proximal (early)  : early-maturation flux markers (BCAS1/TCF7L2/...)
#   Mature (continuum)    : mature-oligo identity (PLP1/MOBP/...) — the continuum bulk
#   High-myelin (NHD)     : the high-structural-myelin genes that peak in this state
#   Stress (NHD)          : heat-shock/ferritin stress genes (CRYAB/FTL/APOD/...)
# no per-nucleus significance stars (single-donor pseudoreplication). This is a
# descriptive marker dotplot (avg expression + % cells), not a test.
#
# Every marker is oligo-credible and passes a min-%-cells floor in >=1 state
# (dropout/near-binary genes are excluded up front — the Fig1 oligo-callout trap).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(ggh4x)   # strip_themed(): pale per-state row-strip fills (house banner idiom)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
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
o <- readRDS(CACHE)
stopifnot("cache missing oligo_state" = "oligo_state" %in% colnames(o@meta.data))

STATE_ORDER <- levels(o$oligo_state)
cat("oligo_state levels:", paste(STATE_ORDER, collapse=" | "), "\n")
cat("counts:\n"); print(table(o$oligo_state))

DefaultAssay(o) <- "RNA"
if (length(SeuratObject::Layers(o, assay = "RNA")) > 1) o <- SeuratObject::JoinLayers(o)
Idents(o) <- "oligo_state"

# ---- marker BLOCKS (block label -> genes). Block order == state column order --
# so the marker->state diagonal reads top-to-bottom. Only oligo-credible markers
# that discriminate their column; dropout/near-binary genes filtered below.
# Markers chosen from the per-state %-detection table (all >MIN_PCT in >=1 state,
# and discriminating their column): early = ENPP6/SLC5A11/GPC5 (decline along
# maturation); mature continuum = the mature-oligo identity bulk; high-myelin =
# OPALIN (12->55% along the axis, the clearest single discriminator) + MOG/CNP;
# stress = the ferritin/heat-shock axis CRYAB/FTL/APOD/B2M (module-defined state).
GENE_BLOCKS <- list(
  "OPC-proximal\n(early)" = c("ENPP6","SLC5A11","GPC5","BCAS1"),
  "Mature\n(continuum)"   = c("PLP1","MOBP","QDPR","OLIG1"),
  "High-myelin\n(NHD)"    = c("OPALIN","MOG","CNP"),
  "Stress\n(NHD)"         = c("CRYAB","FTL","APOD","B2M")
)
gene_lut  <- utils::stack(GENE_BLOCKS)
req_genes <- unlist(GENE_BLOCKS, use.names = FALSE)

# ---- min-%-cells detectability floor (Fig1 dropout-callout trap) -------------
# Drop any marker whose MAX % cells across the 4 states is below MIN_PCT — those
# are near-binary dropout genes that make a meaningless tiny dot.
MIN_PCT <- 10   # percent
present  <- req_genes[req_genes %in% rownames(o)]
absent   <- setdiff(req_genes, present)
if (length(absent)) cat("markers absent from assay (skipped):", paste(absent, collapse=", "), "\n")
stopifnot("no markers present" = length(present) > 0)

pd0 <- DotPlot(o, features = present, group.by = "oligo_state", dot.scale = 4)$data
maxpct <- pd0 %>% group_by(features.plot) %>% summarise(mx = max(pct.exp), .groups="drop")
low_det <- maxpct$features.plot[maxpct$mx < MIN_PCT]
if (length(low_det))
  cat(sprintf("markers below %d%% detection floor in ALL states (dropped as dropout): %s\n",
              MIN_PCT, paste(low_det, collapse=", ")))
dot_genes <- setdiff(present, as.character(low_det))
stopifnot("all markers filtered out" = length(dot_genes) > 0)

pd <- DotPlot(o, features = dot_genes, group.by = "oligo_state", dot.scale = 4)$data
# reviewer-facing: zero-variance gene -> NaN scaled expr -> neutral 0 (logged).
n_bad <- sum(!is.finite(pd$avg.exp.scaled))
if (n_bad) { cat(sprintf("%d/%d scaled-expr non-finite (zero-variance) -> set to 0 (neutral)\n",
                         n_bad, nrow(pd)))
  pd$avg.exp.scaled[!is.finite(pd$avg.exp.scaled)] <- 0 }

pd$gene <- as.character(pd$features.plot)
pd$gene_group <- factor(gene_lut$ind[match(pd$gene, as.character(gene_lut$values))],
                        levels = names(GENE_BLOCKS))
pd$gene <- factor(pd$gene, levels = rev(dot_genes))   # first-listed marker at TOP of its block
pd$id   <- factor(pd$id,   levels = STATE_ORDER)

# Pale per-state row-strip fills — the 4 marker blocks map onto the 4 continuum
# states, so tint each block's strip with its own state colour blended 0.55 -> white
# (the same lighten used for the region banners), replacing the default grey92. Text
# stays plain black. Brings the 3i facet strips into the Fig4 pale-banner family
# Strip order == GENE_BLOCKS order.
PAL_STATE_ROW <- c("#CCBB44","#4477AA","#228833","#AA3377")   # early/mature/high-myelin/stress
strip_fills   <- .region_lighten(PAL_STATE_ROW, 0.55)
strip_state   <- ggh4x::strip_themed(
  background_y = lapply(strip_fills, function(f)
    element_rect(fill = f, colour = "grey80", linewidth = 0.25)))

p_dot <- ggplot(pd, aes(x = id, y = gene, size = pct.exp, fill = avg.exp.scaled)) +
  geom_point(shape = 21, colour = "grey25", stroke = 0.25) +
  ggh4x::facet_grid2(rows = vars(gene_group), scales = "free_y", space = "free",
                     switch = "y", strip = strip_state) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "grey96", high = "#B40426", midpoint = 0,
                       name = "Avg\nexpression",
                       guide = guide_colorbar(order = 1, barheight = unit(1.5,"cm"), barwidth = unit(0.28,"cm"))) +
  scale_size_continuous(name = "% cells", range = c(0.3, 4), limits = c(0, 100),
                        breaks = c(0, 25, 50, 75, 100),
                        guide = guide_legend(order = 2, override.aes = list(fill = "grey50"))) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 7.2, colour = "black"),
    axis.text.y = element_text(size = 7.5, face = "italic", colour = "black"),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.spacing.y = unit(0.12, "lines"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 6.6, face = "plain",
                                     colour = "black", lineheight = 0.85,
                                     margin = margin(r = 2, l = 2)),
    # Strip.background.y intentionally not set here — the pale per-state fills are
    # supplied by strip_themed(background_y=...) above (would be overridden otherwise).
    legend.title = element_text(size = 7), legend.text = element_text(size = 6.5),
    legend.key.size = unit(0.3, "cm"), plot.margin = margin(4, 6, 4, 4))

BN <- "F4f_oligo_state_dotplot"
ggsave(file.path(PANEL, paste0(BN, ".png")), p_dot, width = 3.0, height = 3.15, dpi = 600, device = ragg::agg_png)
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p_dot, width = 3.0, height = 3.15, useDingbats = FALSE)
cat(sprintf("Saved %s.{png,pdf}  (%d markers x %d states)\n", BN, length(dot_genes), nlevels(pd$id)))

writeLines(capture.output(sessionInfo()), file.path(LOGD, "3I_oligo_state_dotplot_FH_sessionInfo.txt"))
cat("=== DONE ===\n", file = stderr())
