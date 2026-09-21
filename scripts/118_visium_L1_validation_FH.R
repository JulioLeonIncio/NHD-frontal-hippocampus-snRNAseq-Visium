#!/usr/bin/env Rscript
# =============================================================================
# 118_visium_L1_validation_FH.R — Does the recovered rim behave as layer I / pia, and does an L1 domain separate from white matter?
# Groups: every cluster of the current integrated object, split by tissue_call_source (spaceranger /
#   he_mask / expression), plus — when integrated_cluster_identity.csv covers every cluster — the curated
#   domains. Per group: UMI per 10,000 for six gene sets: pia / leptomeningeal (PTGDS, DCN, COL1A2,
#   SLC6A13, CEMIP), glia limitans (GFAP, AQP4, ID3), L1 neuronal (RELN, NDNF, CXCL14, CPLX3; expected
#   uninformative at Visium depth), myelin (MBP, PLP1, MOBP, MOG), pan-neuronal (SNAP25, RBFOX3, SLC17A7),
#   erythroid (HBB); median distance to the tissue edge (um, from propose_identity_20260918.R); n per section.
# Acceptance for an "L1" cluster (recorded in IDENTITY_CURATION_20260918.md): myelin <= 15 % of interior WM,
#   GFAP 2-3x L2/3, PTGDS ~2x L2/3, >= 70 % of spots within 200 um of the edge, present in all four sections,
#   holds >= 50 % of the recovered spots.
# Outputs: tables/visium_L1_validation_FH.csv; figures/_diagnostics/L1_validation_<sec>.png (spot maps:
#   cluster/domain, tissue-call source, glia-limitans score, myelin score in one frame per section).
# Run: Rscript scripts/118_visium_L1_validation_FH.R
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr); library(ggplot2); library(patchwork) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R")); source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))
IH <- file.path(dirname(PROJ), "Visium", "integrated_harmony")
OUTD <- file.path(PROJ, "figures", "_diagnostics"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(file.path(IH, "NHD_frontal_integrated_harmony.rds"))
es <- readRDS(file.path(IH, "_spot_edge_src_20260918.rds"))
stopifnot("edge/source table older than the integrated object — re-run propose_identity_20260918.R" =
          file.mtime(file.path(IH, "_spot_edge_src_20260918.rds")) > file.mtime(file.path(IH, "NHD_frontal_integrated_harmony.rds")),
          "edge/source table does not cover every spot" = all(colnames(obj) %in% names(es$edge_um)))
obj$edge_um <- unname(es$edge_um[colnames(obj)]); obj$src <- unname(setNames(es$src, names(es$edge_um))[colnames(obj)])
idf <- read.csv(file.path(IH, "integrated_cluster_identity.csv"))
map <- setNames(idf$identity, sub("^g", "", idf$cluster))
has_dom <- all(levels(obj$seurat_clusters) %in% names(map))
obj$domain <- if (has_dom) unname(map[as.character(obj$seurat_clusters)]) else NA_character_
cat("object:", ncol(obj), "spots; clusters:", nlevels(obj$seurat_clusters), "; curated domains cover all clusters:", has_dom, "\n")

SETS <- list(pia = c("PTGDS","DCN","COL1A2","SLC6A13","CEMIP"), glia_limitans = c("GFAP","AQP4","ID3"),
             L1_neuronal = c("RELN","NDNF","CXCL14","CPLX3"), myelin = c("MBP","PLP1","MOBP","MOG"),
             pan_neuronal = c("SNAP25","RBFOX3","SLC17A7"), erythroid = c("HBB"))
cnt <- GetAssayData(obj, assay = "Spatial", layer = "counts"); tot <- Matrix::colSums(cnt)
per10k <- sapply(SETS, function(g) { g <- intersect(g, rownames(cnt)); 1e4 * Matrix::colSums(cnt[g, , drop = FALSE]) / tot })
md <- cbind(obj@meta.data, per10k)
grp <- function(d, by) d %>% group_by(across(all_of(by))) %>% summarise(n = n(), n_sections = n_distinct(sample_id),
  CON1 = sum(sample_id == "CON_Frontal1"), CON2 = sum(sample_id == "CON_Frontal2"), NHD1 = sum(sample_id == "NHD_Frontal1"), NHD2 = sum(sample_id == "NHD_Frontal2"),
  across(all_of(names(SETS)), ~ round(mean(.x), 2)), median_edge_um = round(median(edge_um, na.rm = TRUE)), frac_within_200um = round(mean(edge_um <= 200, na.rm = TRUE), 2),
  median_UMI = round(median(nCount_Spatial)), .groups = "drop")
t1 <- grp(md, c("seurat_clusters", "src")) %>% mutate(level = "cluster_x_source")
t2 <- grp(md, "seurat_clusters") %>% mutate(level = "cluster")
out <- bind_rows(t2, t1)
if (has_dom) out <- bind_rows(out, grp(md, c("domain", "condition")) %>% mutate(level = "domain_x_condition"), grp(md, "domain") %>% mutate(level = "domain"))
write.csv(out, file.path(PROJ, "tables", "visium_L1_validation_FH.csv"), row.names = FALSE)
cat("\n== per cluster (UMI per 10,000; edge distance) ==\n"); print(as.data.frame(t2), digits = 3)
if (has_dom) { cat("\n== per domain x condition ==\n"); print(as.data.frame(out %>% filter(level == "domain_x_condition")), digits = 3) }
cat("\n== recovered spots by cluster ==\n"); print(table(md$src, md$seurat_clusters))

# spot maps per section, one frame (image coordinates from the Seurat images slot)
for (nm in names(obj@images)) {
  co <- GetTissueCoordinates(obj, image = nm); co <- co[intersect(rownames(co), colnames(obj)), ]
  d <- cbind(co, md[rownames(co), ]); names(d)[1:2] <- c("y", "x")
  d$grp <- if (has_dom) d$domain else as.character(d$seurat_clusters)
  base <- function(p) p + coord_fixed() + scale_y_reverse() + theme_void(base_size = 8) + theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))
  p1 <- base(ggplot(d, aes(x, y, colour = grp)) + geom_point(size = 1.6, stroke = 0) + labs(title = paste(nm, if (has_dom) "domain" else "cluster"), colour = NULL) + guides(colour = guide_legend(override.aes = list(size = 2), nrow = 2)))
  if (has_dom) p1 <- p1 + scale_colour_manual(values = DOM_PAL_SAFE <- { source(file.path(PROJ, "scripts", "_visium_domains_FH.R")); DOM_PAL })
  p2 <- base(ggplot(d, aes(x, y, colour = src)) + geom_point(size = 1.6, stroke = 0) + scale_colour_manual(values = c(spaceranger = "grey80", he_mask = "#2CA02C", expression = "#D62728")) + labs(title = "tissue-call source", colour = NULL) + guides(colour = guide_legend(override.aes = list(size = 2))))
  p3 <- base(ggplot(d, aes(x, y, colour = pmin(glia_limitans, quantile(glia_limitans, .98)))) + geom_point(size = 1.6, stroke = 0) + scale_colour_viridis_c(option = "B") + labs(title = "glia limitans (GFAP AQP4 ID3) per 10k", colour = NULL))
  p4 <- base(ggplot(d, aes(x, y, colour = pmin(myelin, quantile(myelin, .98)))) + geom_point(size = 1.6, stroke = 0) + scale_colour_viridis_c(option = "D") + labs(title = "myelin per 10k", colour = NULL))
  p5 <- base(ggplot(d, aes(x, y, colour = pmin(pia, quantile(pia, .98)))) + geom_point(size = 1.6, stroke = 0) + scale_colour_viridis_c(option = "A") + labs(title = "pia (PTGDS DCN COL1A2 SLC6A13 CEMIP) per 10k", colour = NULL))
  ggsave(file.path(OUTD, paste0("L1_validation_", nm, ".png")), (p1 | p2 | p3) / (p4 | p5 | plot_spacer()), width = 15, height = 10.5, dpi = 150, bg = "white")
}
cat("=== DONE ===\n")
