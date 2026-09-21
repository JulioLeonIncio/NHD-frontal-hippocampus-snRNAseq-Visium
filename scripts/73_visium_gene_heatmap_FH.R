#!/usr/bin/env Rscript
# =============================================================================
# 73_visium_gene_heatmap_FH.R — Figure 6: "F6d_supplychain_genes".
# The gene-level companion that sits beside 6_supplychain_spatial.
# -----------------------------------------------------------------------------
# Why IT EXISTS ("next to it I want a gene heatmap, one big one").
# The spatial panel shows seven programme SCORES in space; it cannot show which genes
# carry them. This panel puts every member gene of the same seven programmes on one
# grid, in the same row order, so the two panels read as a pair: the map says where,
# the heatmap says which.
#
# IT is built from the same object as the maps. Same spots, same 3,000-UMI thinning,
# same programme definitions (_oligo_programmes_FH.R via the prep). Nothing here is
# recomputed from a second source, so the two panels cannot drift apart.
#
# What the colour is. log2(NHD / CON) of the per-spot mean, per gene, per section pair.
# Not a z-score across genes: a z would make a gene that barely moves look identical to
# one that halves, because z is blind to magnitude. The diverging scale is symmetric and
# capped so one extreme gene cannot flatten the rest, and capped values are marked.
#
# Columns are the figure's own BANKSY layer DOMAINS (decided after
# the caveat below was raised). Same annotation as panels b and c: seurat_clusters mapped
# through the same curated map as Visium/figure_Visium/assemble_figure_Visium.R — the
# cluster -> domain map is read from integrated_harmony/integrated_cluster_identity.csv
# (never hard-coded here), ordered superficial -> deep exactly as panel c orders it.
# One laminar axis for the whole figure, which is the point.
#
# The domains are
# strongly confounded with condition, because they are an outcome of the pathology rather
# than a neutral coordinate -- panels b and c show laminar architecture dissolving in NHD,
# so NHD spots are assigned to different domains.
# QC-passed spots per domain, CON vs NHD, differ by up to an order of magnitude in either
# direction (spot counts are read from the live CSVs, printed below as "spots per Banksy
# layer domain" and shipped in the source CSV). So a column is not a like-for-like
# contrast: a sparsely populated domain in one condition is a residue of a layer, not the
# layer. Read this panel as where each programme sits within each condition's own laminar
# organisation -- not as a within-layer NHD-versus-CON test. A condition x domain cell
# with no spots is drawn empty rather than filled or dropped, so the imbalance is visible
# instead of hidden. A condition-neutral geometric-depth axis would be the alternative if
# a quantitative laminar contrast were ever needed.
#
# n = 2 sections per condition, one donor per condition. Spots are pooled within a
# condition x band, so this panel trades the between-section spread (which the sibling
# 6_depth_matched_markers panel still shows per section) for laminar resolution.
# No p-values: condition is perfectly aliased with donor.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggh4x)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_oligo_programmes_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_6", "panels")

# Every spot at native depth (white matter must be shown, and control white
# matter never reaches the depth-matching floor): per spot, gene UMIs per 10,000 UMIs — a
# depth-free ratio — for the union of the programme genes below, from the integrated object
# (cached on the object's mtime). Domain means over hundreds of spots are stable even where
# single spots are shallow (control WM ~200 UMI).
IH_RDS <- file.path(dirname(PROJ), "Visium", "integrated_harmony", "NHD_frontal_integrated_harmony.rds")
CACHE_D <- file.path(PROJ, "data", "_cache_73_gene_umi_per10k_FH.rds")
.oli <- oligo_programme_sets(readRDS(file.path(PROJ, "data", "curated_signatures.rds")))
SETS <- c(list("Complement / MHC-II" = c("C1QA","C1QB","C1QC","CD74","HLA-DRA","HLA-DRB1",
                                         "HLA-DPA1","HLA-DMB","C3","TYROBP","AIF1"),
               "Reactive astrocyte"  = c("GFAP","SERPINA3","VIM","CD44","CLU","AQP4")),
          .oli[c("Structural myelin","Cholesterol (sterol arm)",
                 "Fatty-acid / sphingomyelin","Galactolipid","Lipid uptake / salvage")])
ALLG <- unique(unlist(SETS))
KEY_D <- paste(file.info(IH_RDS)$mtime, paste(sort(ALLG), collapse = ","))
if (file.exists(CACHE_D) && identical(readRDS(CACHE_D)$key, KEY_D)) { d <- readRDS(CACHE_D)$d } else {
  suppressPackageStartupMessages(library(Seurat)); vv <- readRDS(IH_RDS); DefaultAssay(vv) <- "Spatial"; vv <- JoinLayers(vv)
  cnt <- GetAssayData(vv, layer = "counts"); tot <- Matrix::colSums(cnt); g <- intersect(ALLG, rownames(cnt))
  d <- data.frame(spot_id = colnames(vv), Condition = vv$condition, stringsAsFactors = FALSE, check.names = FALSE)
  for (gg in g) d[[paste0("gene_", gg)]] <- as.numeric(cnt[gg, ] / tot * 1e4)
  saveRDS(list(key = KEY_D, d = d), CACHE_D); rm(vv, cnt); invisible(gc()) }
gcols <- grep("^gene_", names(d), value = TRUE)
stopifnot("programme genes missing from the object" = length(gcols) > 40)

# the same seven blocks, in the same order, as 6_supplychain_spatial
have <- sub("^gene_", "", gcols)
SETS <- lapply(SETS, function(g) intersect(g, have))
cat("genes available per programme:\n")
for (k in names(SETS)) cat(sprintf("  %-28s %2d\n", k, length(SETS[[k]])))
SETS <- SETS[lengths(SETS) > 0]

source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))   # cluster -> domain map + level order from integrated_harmony/*.csv
BANDS   <- DOM_BANDS                                 # grey-matter domains in depth_rank order, then white matter (as panel c)
C2L <- file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN", "proportions_q05.csv")
stopifnot("MISSING the c2l table that carries seurat_clusters" = file.exists(C2L))
dom <- read.csv(C2L, stringsAsFactors = FALSE) %>%
  transmute(spot_id, band = unname(DOM_MAP[as.character(seurat_clusters)]))
d <- d %>% left_join(dom, by = "spot_id") %>%
  filter(band %in% BANDS) %>% mutate(band = factor(band, levels = BANDS))
cat("spots per Banksy layer domain (NOT balanced — see the header):\n")
print(table(d$band, d$Condition))

long <- bind_rows(lapply(names(SETS), function(k)
  d %>% filter(!is.na(band)) %>%
    select(band, Condition, all_of(paste0("gene_", SETS[[k]]))) %>%
    pivot_longer(-c(band, Condition), names_to = "gene", values_to = "e") %>%
    mutate(gene = sub("^gene_", "", gene), programme = k)))

# The previous version collapsed the two conditions into one log2 ratio, which hid the
# levels, and on a single shared scale the astrocyte block saturated at +3.6 while the
# myelin block sat near -0.3 and rendered almost white.
# The scale unit is the pathway BLOCK ("do it per pathway block or
# program"), which is what makes genes comparable within a block: a gene that moves twice
# as much as its neighbour now looks twice as strong, which a per-gene z cannot show.
#
# It is done in two steps, and the second step is not optional. Standardising raw
# expression across a block would let the abundant members set both the location and the
# scale -- the per-gene mean expression inside a block spans 24.7-fold for structural
# myelin (mal 0.090 vs PLP1 2.228) and 17.2-fold for the astrocyte pair (SERPINA3 0.112
# vs GFAP 1.926), so PLP1 and GFAP would define their blocks and every other member would
# render as a flat offset regardless of condition. So:
#   1. CENTRE each gene on its own mean, removing abundance as a source of location;
#   2. scale by the block's pooled SD of those centred values -- one shared unit per
#      pathway, which is the comparability the block-level choice is for.
# Magnitude is still not comparable between programmes. The absolute log2(NHD/CON) per
# programme per band is printed below and shipped in the source CSV; quote that.
per_band <- long %>%
  group_by(programme, gene, band, Condition) %>%
  summarise(m = mean(e), .groups = "drop")

raw_lfc <- per_band %>% pivot_wider(names_from = Condition, values_from = m) %>%
  filter(CON > 0) %>% mutate(lfc = log2((NHD + 1e-4)/(CON + 1e-4)))
cat("\n== absolute log2(NHD/CON) per programme per band (NOT on the panel; for the legend) ==\n")
print(as.data.frame(raw_lfc %>% group_by(programme, band) %>%
      summarise(median_lfc = round(median(lfc), 2), .groups = "drop") %>%
      pivot_wider(names_from = band, values_from = median_lfc)), row.names = FALSE)

keep <- raw_lfc %>% distinct(gene)
mat <- per_band %>% semi_join(keep, by = "gene") %>%
  group_by(gene) %>% mutate(dm = m - mean(m)) %>%                   # 1. centre per gene
  group_by(programme) %>%                                           # 2. scale per block
  mutate(.sd = stats::sd(dm),
         z = if (is.na(.sd[1]) || .sd[1] == 0) 0 else dm / .sd[1]) %>%
  ungroup() %>%
  mutate(programme = factor(programme, levels = names(SETS)),
         Condition = factor(Condition, levels = c("CON","NHD")))

# Nothing is CLIPPED ("I don't like these x inside the heatmap, so
# scale it differently"). The panel previously capped at |z| = 2 and marked the 25 clipped
# cells with a cross. Extending the limits linearly to the true maximum instead would have
# emptied the panel: |z| runs to 3.72 but its median is 0.48 and its 90th percentile 1.79,
# so a linear ramp to 3.72 puts most of the data in the palest tenth of the scale.
# The colour position is therefore a sign-preserving compressive function of z,
#     pos = sign(z) * (|z| / max|z|) ^ POW ,
# which keeps the mid-range separable while still resolving the tails. The legend is
# labelled in z units at the positions those z values actually occupy, so the mapping is
# non-linear but readable and nothing is hidden.
POW <- 0.60
ZMAX <- max(abs(mat$z))
compress <- function(v) sign(v) * (abs(v) / ZMAX)^POW
mat <- mat %>% mutate(lfc_c = compress(z))
BR_Z  <- c(-3, -1, 0, 1, 3)
BR_Z  <- BR_Z[abs(BR_Z) <= ZMAX]
cat(sprintf("\ngenes drawn: %d | |z| max %.2f, median %.2f | compressive colour (power %.2f), nothing clipped\n",
            dplyr::n_distinct(mat$gene), ZMAX, median(abs(mat$z)), POW))

# order genes within a block by their NHD-minus-CON z, so the block structure reads
ord <- mat %>% group_by(programme, gene) %>%
  summarise(v = mean(z[Condition=="NHD"]) - mean(z[Condition=="CON"]), .groups = "drop") %>%
  arrange(programme, v) %>% pull(gene)
mat$gene <- factor(mat$gene, levels = unique(ord))

GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
CELL <- 0.091   # cell width, inches: 51 gene columns + strips + left labels = ~4.8 in, placed at 100 % beside e
CELL_H <- 0.125 # cell height, inches: rows taller than wide so the five domain rows read at page scale ("d looks squeezed")
# Horizontal layout: genes along x (ordered within programme blocks), domains as
# rows (superficial at the top, white matter at the bottom), CON over NHD as row facets, programme
# blocks as column facets, one colour bar underneath. Canvas measured from the gtable.
mat <- mat %>% mutate(band = factor(band, levels = rev(levels(factor(band, levels = BANDS)))))
p <- ggplot(mat, aes(x = gene, y = band, fill = lfc_c)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  scale_fill_gradient2(low = GSEA_STEEL, mid = "white", high = GSEA_ROSE, midpoint = 0,
                       limits = c(-1, 1), breaks = compress(BR_Z),
                       labels = ifelse(BR_Z == 0, "0", sprintf("%+.0f", BR_Z)),
                       name = "centred expression, scaled within pathway (SD)",
                       guide = guide_colourbar(title.position = "left", title.vjust = 0.85, direction = "horizontal")) +
  ggh4x::facet_grid2(Condition ~ programme, scales = "free_x", space = "free_x", switch = "y",   # CON / NHD strips on the left
                     labeller = labeller(programme = function(x) { x <- sub(" \\(", "\n(", x); x <- sub("^Reactive astrocyte$", "Reactive\nastrocyte", x)
                       x <- sub("^Structural myelin$", "Structural\nmyelin", x); x <- sub("^Fatty-acid / sphingomyelin$", "Fatty acid/\nsphingo-\nmyelin", x)
                       x <- sub("^Galactolipid$", "Galacto-\nlipid", x); x <- sub("^Lipid uptake / salvage$", "Lipid\nuptake /\nsalvage", x); x }),
                     strip = ggh4x::strip_themed(
                       background_x = ggh4x::elem_list_rect(fill = "grey94", colour = NA),
                       text_x = ggh4x::elem_list_text(size = 5.4, colour = "black", lineheight = 0.85),
                       background_y = ggh4x::elem_list_rect(fill = NA, colour = NA),
                       text_y = ggh4x::elem_list_text(size = 6.4, colour = "black", angle = 0, hjust = 1))) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8) +
  theme(axis.text.x = element_text(size = 5.4, colour = "black", face = "italic", angle = 90, hjust = 1, vjust = 0.5),   # 5.4 pt at 100 % (0.083-in columns)
        axis.text.y = element_text(size = 6.0, colour = "black"),
        axis.ticks = element_blank(), axis.line = element_blank(),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.2),
        panel.spacing.x = unit(0.08, "cm"), panel.spacing.y = unit(0.08, "cm"),
        strip.placement = "outside", strip.clip = "off",
        legend.position = "bottom", legend.key.width = unit(0.7, "cm"),
        legend.key.height = unit(0.18, "cm"), legend.title = element_text(size = 6.0, colour = "black"),
        legend.text = element_text(size = 6.0, colour = "black"), legend.margin = margin(t = -1),
        legend.box.spacing = unit(0.06, "cm"),
        plot.margin = margin(1, 2, 1, 1)) +
  ggh4x::force_panelsizes(rows = unit(length(BANDS) * CELL_H, "in"),
                          cols = unit(CELL * (mat %>% distinct(programme, gene) %>% count(programme) %>%
                                                 arrange(factor(programme, levels = levels(mat$programme))) %>% pull(n)), "in"))   # square cells, CELL in
g_d <- ggplotGrob(p)
W <- sum(grid::convertWidth(g_d$widths, "in", valueOnly = TRUE)) + 0.06
H <- sum(grid::convertHeight(g_d$heights, "in", valueOnly = TRUE)) + 0.06
BN <- "F6d_supplychain_genes"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nSaved %s.{png,pdf} (%.2f x %.2f in)\n", BN, W, H))
write.csv(raw_lfc %>% select(programme, gene, band, CON, NHD, lfc),
          file.path(TDIR, "figF6d_supplychain_genes_FH.csv"), row.names = FALSE)
cat("=== DONE ===\n")
