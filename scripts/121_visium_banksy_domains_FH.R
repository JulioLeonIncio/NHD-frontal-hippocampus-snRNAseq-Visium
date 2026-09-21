#!/usr/bin/env Rscript
# =============================================================================
# 121_visium_banksy_domains_FH.R — BANKSY + Harmony domain segmentation of the NHD Visium sections with the SCZ-project recipe (Visium/_banksy_reference_SCZ/Bansky_alone.R; banksy_NHD_frontal/banksy_NHD_frontal.R):
#   k_geom = 6, agf = FALSE, lambda = 0.2 (domain mode), 20 PCs, Harmony over section, Leiden on the Harmony
#   embedding at several resolutions; layer annotation with the SCZ dictionary, whose L1 entry is astrocytic /
#   interlaminar (GFAP, AQP4, HOPX, FABP7) plus NDNF/RELN/CXCL14 and whose WM entry is myelin.
# Evidence run to decide whether the Fig 6 domain map should come from BANKSY (as in SCZ) instead of
#   the expression-only Seurat/Harmony clustering that merges the gliotic layer I of the NHD sulci with white matter.
#   Reads the current outs (VISIUM_OUTS_MANIFEST.json; two-tier tissue call) and the house QC floor (50 / 25, no mito cut).
# OUT (Visium/banksy_NHD_frontal_20260918/): spe rds; per-spot csv (spot_id, section, cluster at each resolution,
#   argmax layer per cluster at each resolution); layer-score tables; domain maps per resolution; overlay of the argmax
#   layers against the pathologists' L1/WM lines (tables/annotation_dashes_20260918_mapunits.csv) with the same zone
#   statistics as 119 (composition of the L1 band, WM band).
# Run: OMP_NUM_THREADS=4 Rscript scripts/121_visium_banksy_domains_FH.R
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(SpatialExperiment); library(SingleCellExperiment); library(Banksy); library(harmony)
  library(Matrix); library(dplyr); library(ggplot2); library(patchwork); library(jsonlite) })
SEED <- 1000; set.seed(SEED)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))
OUT <- file.path(VIS_BASE, "banksy_NHD_frontal_20260918"); dir.create(OUT, showWarnings = FALSE)
TAB <- file.path(PROJ, "tables"); C2L <- file.path(VIS_BASE, "cell2location", "c2l_MAIN")
QC_MIN_COUNTS <- 50; QC_MIN_FEATURES <- 25          # house floor (integrated_harmony.R); no mito ceiling
K_GEOM <- 6; USE_AGF <- FALSE; LAMBDA <- 0.2; NPCS <- 20; NORM_SCALE <- 5000; N_HVG <- 2000
RES_SET <- c(0.5, 0.65, 0.8, 1.0)
say <- function(...) cat(sprintf(...), "\n")

# ---- load (manifest outs), QC, SPE ----------------------------------------------------------
to_spe <- function(seu, name, dir) {
  pos <- read.csv(file.path(dir, "spatial", "tissue_positions_list.csv"), header = FALSE,
                  col.names = c("barcode","in_tissue","array_row","array_col","pxl_row_in_fullres","pxl_col_in_fullres"))
  pos <- pos[match(colnames(seu), pos$barcode), ]
  spe <- SpatialExperiment(assays = list(counts = GetAssayData(seu, assay = "Spatial", layer = "counts")),
                           colData = data.frame(sample_id = name, barcode = colnames(seu), nCount = seu$nCount_Spatial, nFeature = seu$nFeature_Spatial,
                                                array_row = pos$array_row, array_col = pos$array_col, pxl_row = pos$pxl_row_in_fullres, pxl_col = pos$pxl_col_in_fullres),
                           spatialCoords = as.matrix(pos[, c("array_col","array_row")]))
  colnames(spe) <- paste0(name, "_", colnames(seu)); spe
}
spe_list <- lapply(VIS_SECTIONS, function(nm) {
  d <- vis_outs(nm); seu <- Load10X_Spatial(d, filename = "filtered_feature_bc_matrix.h5", slice = nm)
  n0 <- ncol(seu); seu <- subset(seu, subset = nCount_Spatial >= QC_MIN_COUNTS & nFeature_Spatial >= QC_MIN_FEATURES)
  say("  %s: %d -> %d spots", nm, n0, ncol(seu)); to_spe(seu, nm, d) })
names(spe_list) <- VIS_SECTIONS
excl <- function(g) unique(c(grep("^MT-", g, value = TRUE), grep("^RP[SL][0-9]", g, value = TRUE), grep("^MTRNR", g, value = TRUE)))
spe_list <- lapply(spe_list, function(spe) spe[setdiff(rownames(spe), excl(rownames(spe))), ])
hv <- list()
spe_list <- lapply(names(spe_list), function(nm) { spe <- spe_list[[nm]]; cts <- counts(spe)
  assay(spe, "normcounts") <- as(Matrix::t(Matrix::t(cts) / pmax(Matrix::colSums(cts), 1)) * NORM_SCALE, "dgCMatrix")
  seu <- CreateSeuratObject(counts = cts, assay = "Spatial") |> NormalizeData(verbose = FALSE) |> FindVariableFeatures(nfeatures = N_HVG, verbose = FALSE)
  hv[[nm]] <<- VariableFeatures(seu); spe }); names(spe_list) <- VIS_SECTIONS
hvgs <- intersect(Reduce(union, hv), Reduce(intersect, lapply(spe_list, rownames))); say("union HVGs: %d", length(hvgs))
spe_list <- lapply(spe_list, function(spe) computeBanksy(spe[hvgs, ], assay_name = "normcounts", compute_agf = FALSE, k_geom = K_GEOM))
spe <- do.call(cbind, spe_list); say("joint: %d spots", ncol(spe))
spe <- runBanksyPCA(spe, use_agf = USE_AGF, lambda = LAMBDA, npcs = NPCS, seed = SEED)
pca_name <- sprintf("PCA_M%s_lam%s", as.numeric(USE_AGF), LAMBDA)
H <- harmony::RunHarmony(data_mat = reducedDim(spe, pca_name), meta_data = as.data.frame(colData(spe)), vars_use = "sample_id", max_iter = 20, verbose = FALSE)
reducedDim(spe, "Harmony_BANKSY") <- H
for (r in RES_SET) spe <- clusterBanksy(spe, dimred = "Harmony_BANKSY", resolution = r, seed = SEED)
cn <- clusterNames(spe); say("cluster columns: %s", paste(cn, collapse = ", "))

# ---- SCZ layer dictionary (Bansky_alone.R step 16) + argmax per cluster ------------------------------
layer_markers <- list(L1 = c("NDNF","RELN","CXCL14","GFAP","AQP4","HOPX","FABP7"), L2_3 = c("CUX2","RASGRF2","MDGA1","LAMP5","CARTPT"),
  L4 = c("RORB","PLCH1","RSPO1"), L5 = c("BCL11B","FEZF2","PCP4","CRYM","TOX","ETV1"), L6 = c("TLE4","FOXP2","SEMA3E","NR4A2","OPRK1"),
  WM = c("MBP","MOBP","PLP1","MOG","MAG","OPALIN"), Vasc = c("HBA1","HBA2","HBB","CLDN5","VWF","PECAM1"))
# scores need genes outside the HVG set: recompute from the full normalised matrices
full <- do.call(cbind, lapply(VIS_SECTIONS, function(nm) { d <- vis_outs(nm); seu <- Load10X_Spatial(d, filename = "filtered_feature_bc_matrix.h5", slice = nm)
  m <- GetAssayData(seu, assay = "Spatial", layer = "counts"); colnames(m) <- paste0(nm, "_", colnames(m)); m }))
full <- full[, colnames(spe)]; norm_full <- Matrix::t(Matrix::t(full) / pmax(Matrix::colSums(full), 1)) * NORM_SCALE
z_scores <- function(cl) { cls <- factor(colData(spe)[[cl]])
  sc <- sapply(names(layer_markers), function(L) { g <- intersect(layer_markers[[L]], rownames(norm_full))
    v <- Matrix::colMeans(norm_full[g, , drop = FALSE]); tapply(v, cls, mean) })
  sc <- as.data.frame(sc); rownames(sc) <- paste0("c", levels(cls)); sc }
per_spot <- data.frame(spot_id = colnames(spe), section = colData(spe)$sample_id, array_row = colData(spe)$array_row, array_col = colData(spe)$array_col)
for (cl in cn) {
  sc <- z_scores(cl); stopifnot("a resolution collapsed to a single cluster" = nrow(sc) >= 2); scz <- scale(as.matrix(sc))
  bad <- colnames(scz)[colSums(is.na(scz)) == nrow(scz)]; if (length(bad)) say("  WARNING: zero-variance layer column(s) %s cannot win any argmax at %s", paste(bad, collapse = ","), cl)
  best <- colnames(scz)[apply(scz, 1, which.max)]
  write.csv(cbind(cluster = rownames(sc), round(sc, 3), argmax = best), file.path(OUT, paste0("layer_scores_", cl, ".csv")), row.names = FALSE)
  per_spot[[cl]] <- as.character(colData(spe)[[cl]]); per_spot[[paste0("layer_", cl)]] <- best[match(paste0("c", per_spot[[cl]]), rownames(sc))]
  say("\n%s: %d clusters -> %s", cl, nrow(sc), paste(rownames(sc), best, sep = "=", collapse = " "))
}
# spot coordinates in the Fig-6 map frame (01f export) for overlays
xy <- read.csv(file.path(C2L, "spot_coords.csv")); per_spot <- per_spot %>% left_join(xy %>% transmute(spot_id = cell, x, y), by = "spot_id")
write.csv(per_spot, file.path(OUT, "banksy_per_spot_layers_20260918.csv"), row.names = FALSE)
saveRDS(spe, file.path(OUT, "BANKSY_Harmony_NHD_frontal_spe_20260918.rds"))

# ---- zone statistics against the pathologists' lines (same components/labels as 119) --------------------
dash <- read.csv(file.path(TAB, "annotation_dashes_20260918_mapunits.csv")); UM <- 1000 / 81.5
source(file.path(PROJ, "scripts", "_visium_L1_annotation_rule_FH.R"))   # is_L1(), H_LINK — single source
LPAL <- c(L1 = "#64B5CD", L2_3 = "#4C72B0", L4 = "#DD8452", L5 = "#55A868", L6 = "#C44E52", WM = "#8172B3", Vasc = "#DA8BC3")
zone_all <- list()
for (sec in VIS_SECTIONS) {
  d <- dash[dash$section == sec, ]; d <- d[!duplicated(round(d[, c("x","y")] / 2)), ]
  hc <- hclust(dist(d[, c("x","y")]), method = "single"); d$comp <- cutree(hc, h = H_LINK)
  s <- per_spot[per_spot$section == sec & !is.na(per_spot$x), ]; xr <- range(s$x); yr <- range(s$y)
  ce <- d %>% group_by(comp) %>% summarise(cx = mean(x), cy = mean(y), .groups = "drop") %>% mutate(fx = (cx - xr[1]) / diff(xr), fy = (cy - yr[1]) / diff(yr), label = ifelse(is_L1(sec, fx, fy), "L1", "WM"))
  d$label <- ce$label[match(d$comp, ce$comp)]
  nn <- function(px, py, q) if (!nrow(q)) Inf else apply(sqrt(outer(px, q$x, "-")^2 + outer(py, q$y, "-")^2), 1, min) * UM
  s$d_L1 <- nn(s$x, s$y, d[d$label == "L1", ]); s$d_WM <- nn(s$x, s$y, d[d$label == "WM", ])
  for (cl in cn) {
    s$lay <- s[[paste0("layer_", cl)]]
    zone_all[[paste(sec, cl)]] <- bind_rows(
      s %>% filter(d_L1 <= 100) %>% count(layer = lay) %>% mutate(zone = "L1 band (<=100 um)", section = sec, res = cl),
      s %>% filter(d_WM <= 100) %>% count(layer = lay) %>% mutate(zone = "WM band (<=100 um)", section = sec, res = cl))
  }
  cl <- cn[2]  # res 0.65 (the SCZ setting) for the overlay
  s$layer <- factor(s[[paste0("layer_", cl)]], levels = names(LPAL))
  p <- ggplot(s, aes(x, y, colour = layer)) + geom_point(size = 1.2, stroke = 0) + geom_point(data = d, aes(x, y), colour = ifelse(d$label == "L1", "black", "grey40"), size = 0.25, inherit.aes = FALSE) +
    scale_colour_manual(values = LPAL, drop = FALSE) + coord_equal() + theme_void(base_size = 9) + theme(legend.position = "bottom") + labs(title = paste(sec, "— BANKSY layers (", cl, ") vs the pathologists' lines"))
  ggsave(file.path(OUT, paste0("banksy_vs_annotation_", sec, ".png")), p, width = 8, height = 8.5, dpi = 150, bg = "white")
}
zone <- bind_rows(zone_all) %>% group_by(section, res, zone) %>% mutate(frac = round(n / sum(n), 3)) %>% ungroup()
write.csv(zone, file.path(OUT, "banksy_zone_composition_vs_annotation.csv"), row.names = FALSE)
cat("\n== composition of the annotated L1 band by BANKSY layer (res 0.65) ==\n"); print(as.data.frame(zone %>% filter(res == cn[2], grepl("L1", zone)) %>% arrange(section, -n)))
cat("\n== composition of the annotated WM band by BANKSY layer (res 0.65) ==\n"); print(as.data.frame(zone %>% filter(res == cn[2], grepl("WM", zone)) %>% arrange(section, -n)))
cat("=== DONE ===\n")
