#!/usr/bin/env Rscript
# =============================================================================
# integrated_harmony.R — NHD Visium frontal — standard INTEGRATED pipeline (Seurat + Harmony, all 4 CON+NHD samples). Overview/atlas view (mirrors SCZ Fig1_supp b/c):
#   - clusters UMAP (by cluster / condition / sample)
#   - spatial clusters on tissue (per sample)
# this integrated view is for overview; the layer-LOSS claim uses the
# per-condition analysis (per_condition_layers.R) to avoid integration bias.
# Sections are read from the outs folders named in Visium/VISIUM_OUTS_MANIFEST.json (helper
# NHD_frontal_hippo_rebuild/scripts/_visium_outs_FH.R).
# <sec>_masked/outs folders of scripts/117b: Space Ranger's automatic tissue call had excluded the pale
# layer-I / pial rim (CON 379 + 1,002 barcodes at median 4,334 / 3,297 UMI), so the call is re-made in
# three tiers (spaceranger / he_mask / expression; per-barcode source in <outs>/tissue_call.csv).
# NHD_Frontal2 inherits the 108d reflection + registered tissue call through NHD_Frontal2_registered.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(ggplot2); library(patchwork); library(dplyr)
})
set.seed(1000)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
BASE <- file.path(dirname(Sys.getenv("NHD_PROJ")), "Visium")
OUT  <- file.path(BASE, "integrated_harmony"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

source(file.path(dirname(BASE), "NHD_frontal_hippo_rebuild", "scripts", "_visium_outs_FH.R"))   # VIS_SECTIONS, vis_outs()
SAMPLES <- list(
  CON_Frontal1 = list(dir=vis_outs("CON_Frontal1"), cond="CON"),
  CON_Frontal2 = list(dir=vis_outs("CON_Frontal2"), cond="CON"),
  NHD_Frontal1 = list(dir=vis_outs("NHD_Frontal1"), cond="NHD"),
  NHD_Frontal2 = list(dir=vis_outs("NHD_Frontal2"), cond="NHD")
)
for (nm in names(SAMPLES)) {                                   # log the tissue-call tiers (115 parses this log)
  tc <- file.path(SAMPLES[[nm]]$dir, "tissue_call.csv")
  if (file.exists(tc)) { t <- read.csv(tc); t <- t[t$in_tissue_analysis == 1, ]
    cat(sprintf("  tissue call %s: %s\n", nm, paste(names(table(t$tissue_call_source)), table(t$tissue_call_source), sep="=", collapse=" "))) }
}
QC_MIN_COUNTS<-50; QC_MIN_FEAT<-25; NDIMS<-30; RES<-0.5   # permissive floor to retain low-count white matter
# The floor is 50/25 (permissive on purpose, to retain low-count white matter); this is the
# value applied on the line above and the value quoted in the Methods.
# No mito ceiling (condition-confounded; high-mito in NHD is decoupled from library size
# = candidate biology). percent.mt is REGRESSED OUT for clustering (not used to drop spots).

CKPT <- file.path(OUT, "_presub_integrated.rds")
RESUME <- Sys.getenv("RESUME") == "1" && file.exists(CKPT) &&
  all(file.mtime(CKPT) > sapply(SAMPLES, function(x) file.mtime(file.path(x$dir, "filtered_feature_bc_matrix.h5"))))
if (RESUME) cat("== RESUME=1: reading the checkpoint", basename(CKPT), "(newer than every filtered h5) ==\n")
if (!RESUME) {
cat("== load + QC (consensus: no mito cut; regress percent.mt) ==\n")
seu_list <- lapply(names(SAMPLES), function(nm){
  d <- SAMPLES[[nm]]$dir
  s <- Load10X_Spatial(d, filename="filtered_feature_bc_matrix.h5", slice=nm)
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern="^MT-")
  n0 <- ncol(s)
  s <- subset(s, subset = nCount_Spatial>=QC_MIN_COUNTS & nFeature_Spatial>=QC_MIN_FEAT)  # no mito cut
  s$sample_id <- nm; s$condition <- SAMPLES[[nm]]$cond
  cat(sprintf("  %s: %d -> %d spots\n", nm, n0, ncol(s)))
  s
})
names(seu_list) <- names(SAMPLES)

cat("== merge + normalize + PCA ==\n")
obj <- merge(seu_list[[1]], y=seu_list[-1], add.cell.ids=names(SAMPLES))
obj <- JoinLayers(obj)
obj <- NormalizeData(obj, verbose=FALSE)
obj <- FindVariableFeatures(obj, nfeatures=2000, verbose=FALSE)
obj <- ScaleData(obj, vars.to.regress="percent.mt", verbose=FALSE)  # consensus: regress mito, don't filter
obj <- RunPCA(obj, npcs=NDIMS, verbose=FALSE)

cat("== Harmony integration (group.by sample_id) ==\n")
obj <- RunHarmony(obj, group.by.vars="sample_id", reduction.use="pca",
                  reduction.save="harmony", verbose=FALSE)

cat("== cluster + UMAP on harmony ==\n")
obj <- FindNeighbors(obj, reduction="harmony", dims=1:NDIMS, verbose=FALSE)
obj <- FindClusters(obj, resolution=RES, verbose=FALSE)
cat("  clusters at resolution", RES, ":", nlevels(obj$seurat_clusters), " | ", paste(table(obj$seurat_clusters),collapse="/"), "\n")
saveRDS(obj, CKPT)   # checkpoint before the seeded sub-clustering (inspection / RESUME without a 25-min re-run)
} else { obj <- readRDS(CKPT); cat("  clusters at resolution", RES, ":", nlevels(obj$seurat_clusters), " | ", paste(table(obj$seurat_clusters),collapse="/"), "\n") }
# ---- upper-cortex sub-clustering --------------------------------------------------------
# At resolution 0.5 the control upper cortex (layers II/III and IV) forms one cluster: the pial
# band with higher CUX2/LAMP5 and the band beneath it with higher RORB/PCP4 are merged. That
# cluster is sub-clustered on the same graph (resolution 0.3), which separates the two bands
# (superficial L2/3, deeper L4) in both control sections. The cluster is identified by its
# marker profile, not by its number: the largest cluster in which both the L2/3 markers and
# RORB sit in the top quartile of cluster means. Sub-cluster ids: "<k>" keeps the superficial
# band, "<k>_L4" is the RORB-high band; both are re-coded to plain integers below.
SUB_RES <- 0.3
mk <- c("CUX2","LAMP5","RORB"); ex <- GetAssayData(obj, layer="data")[mk, ]
cm <- t(sapply(levels(obj$seurat_clusters), function(k) Matrix::rowMeans(ex[, obj$seurat_clusters == k, drop=FALSE])))
q75 <- apply(cm, 2, quantile, 0.75)
cand <- rownames(cm)[cm[, "RORB"] >= q75["RORB"] & (cm[, "CUX2"] >= q75["CUX2"] | cm[, "LAMP5"] >= q75["LAMP5"])]
cat("  cluster means (CUX2 / LAMP5 / RORB):\n"); print(round(cm, 3))
cand <- cand[which.max(table(obj$seurat_clusters)[cand])]
stopifnot("no upper-cortex cluster with both L2/3 and L4 markers in the top quartile" = length(cand) == 1)
cat("  upper-cortex cluster to sub-cluster:", cand, "(n =", sum(obj$seurat_clusters == cand), ")\n")
set.seed(1000)
obj <- FindSubCluster(obj, cluster=cand, graph.name="Spatial_snn", subcluster.name="sub_upper", resolution=SUB_RES)
subs <- unique(obj$sub_upper[obj$seurat_clusters == cand])
cat("  sub-clusters of", cand, "at resolution", SUB_RES, ":", paste(subs, table(obj$sub_upper[obj$seurat_clusters == cand])[subs], sep="=", collapse=" "), "\n")
mk2 <- intersect(c("CUX2","LAMP5","RORB","PCP4","TLE4","PTGDS","DCN","GFAP","AQP4","ID3","MBP","PLP1","SNAP25","HBB"), rownames(obj))
ex2 <- GetAssayData(obj, layer="data")[mk2, ]
cat("  sub-cluster marker means (log-normalised) and composition:\n")
print(round(t(sapply(subs, function(s2) Matrix::rowMeans(ex2[, obj$sub_upper == s2, drop=FALSE]))), 3))
print(table(obj$sub_upper[obj$seurat_clusters == cand], obj$sample_id[obj$seurat_clusters == cand]))
if (file.exists(file.path(SAMPLES[[1]]$dir, "tissue_call.csv"))) {                  # share of recovered (he_mask / expression) spots per sub-cluster
  tcs <- do.call(rbind, lapply(names(SAMPLES), function(nm) { t <- read.csv(file.path(SAMPLES[[nm]]$dir, "tissue_call.csv")); t$cell <- paste0(nm, "_", t$barcode); t[, c("cell","tissue_call_source")] }))
  src <- tcs$tissue_call_source[match(colnames(obj), tcs$cell)]
  print(table(obj$sub_upper[obj$seurat_clusters == cand], src[obj$seurat_clusters == cand]))
}
# With the layer-I rim recovered the upper-cortex cluster splits into more than two sub-clusters
# (1,502 / 1,388 / 499 / 107 at resolution 0.3). The two bands are now defined by the largest gap in the
# sub-cluster RORB means (sorted): sub-clusters above the gap = the RORB-high band (one new id), below = the
# superficial band (keeps the parent id).
stopifnot("upper-cortex cluster did not split" = length(subs) >= 2)
rb <- sort(sapply(subs, function(s) mean(ex["RORB", obj$sub_upper == s])))
gap <- which.max(diff(rb)); l4 <- names(rb)[(gap + 1):length(rb)]
cat("  RORB means:", paste(names(rb), round(rb, 3), sep = "=", collapse = " "), "| largest gap after", names(rb)[gap], "-> RORB-high band =", paste(l4, collapse = "+"), "\n")
newid <- as.character(obj$seurat_clusters); k_new <- as.character(nlevels(obj$seurat_clusters))   # next free integer
newid[obj$sub_upper %in% l4] <- k_new
obj$seurat_clusters <- factor(newid, levels = as.character(sort(as.integer(unique(newid)))))
Idents(obj) <- "seurat_clusters"
cat("  sub-clustered", cand, "->", cand, "(superficial, CUX2/LAMP5) and", k_new, "(RORB/PCP4 band):",
    sum(newid == cand), "/", sum(newid == k_new), "spots\n")

# ---- layer I / pia vs parenchymal vascular-immune sub-clustering ----------------------
# The recovered pial rim (117b) and the NHD sulcal leptomeninges fall into one cluster together with the
# parenchymal vascular/immune spots (DCN/ID3/PTGDS and CLDN5/PECAM1 both in the top quartile). That cluster
# is sub-clustered on the same graph (resolution 0.3, seeded) and each sub-cluster is classed by its
# position: median distance to the tissue edge (nearest non-tissue barcode or capture-area border, from
# <outs>/tissue_call.csv) below 0.75 x the tissue-wide median = PIAL (own new id, curated as L1/pia);
# otherwise it keeps the parent id (Vasc/immune). Pial sub-clusters are kept as separate ids because the
# control glia limitans (AQP4/GFAP, pial vessels) and the NHD leptomeninges (DCN/COL1A2, CD74/C1QB) differ.
edge_um <- local({
  tc <- do.call(rbind, lapply(names(SAMPLES), function(nm) { t <- read.csv(file.path(SAMPLES[[nm]]$dir, "tissue_call.csv")); t$cell <- paste0(nm, "_", t$barcode); t$section <- nm; t }))
  e <- setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  for (nm in names(SAMPLES)) {
    sfj <- jsonlite::fromJSON(file.path(SAMPLES[[nm]]$dir, "spatial", "scalefactors_json.json")); ppu <- sfj$spot_diameter_fullres / 55
    t <- tc[tc$section == nm, ]; inn <- t[t$in_tissue_analysis == 1, ]; outp <- t[t$in_tissue_analysis == 0, ]
    cells <- intersect(inn$cell, colnames(obj)); inn <- inn[match(cells, inn$cell), ]
    d_out <- if (nrow(outp)) apply(sqrt(outer(inn$pxl_row_in_fullres, outp$pxl_row_in_fullres, "-")^2 + outer(inn$pxl_col_in_fullres, outp$pxl_col_in_fullres, "-")^2), 1, min) / ppu else Inf
    d_edge <- pmin(pmin(inn$array_row, 77 - inn$array_row) * 86.6, pmin(inn$array_col, 127 - inn$array_col) * 50)   # hex lattice: rows 86.6 um, staggered cols 50 um
    e[cells] <- pmin(d_out, d_edge)
  }
  e
})
obj$edge_um <- edge_um
mkp <- c("DCN","ID3","PTGDS","CLDN5","PECAM1","PDGFRB"); exp_ <- GetAssayData(obj, layer="data")[mkp, ]
cmp <- t(sapply(levels(obj$seurat_clusters), function(k) Matrix::rowMeans(exp_[, obj$seurat_clusters == k, drop=FALSE])))
pia_score <- rowMeans(scale(cmp[, c("DCN","ID3","PTGDS")])); vasc_score <- rowMeans(scale(cmp[, c("CLDN5","PECAM1","PDGFRB")]))
candp <- rownames(cmp)[pia_score >= quantile(pia_score, 0.75) & vasc_score >= quantile(vasc_score, 0.75)]
candp <- candp[which.max(table(obj$seurat_clusters)[candp])]
stopifnot("no pia/vascular cluster (top-quartile DCN/ID3/PTGDS and CLDN5/PECAM1/PDGFRB)" = length(candp) == 1)
cat("  pia/vascular cluster to sub-cluster:", candp, "(n =", sum(obj$seurat_clusters == candp), ")\n")
set.seed(1000)
obj <- FindSubCluster(obj, cluster=candp, graph.name="Spatial_snn", subcluster.name="sub_pia", resolution=SUB_RES)
subp <- sort(unique(obj$sub_pia[obj$seurat_clusters == candp]))
med_all <- median(obj$edge_um, na.rm=TRUE)
edp <- sapply(subp, function(s2) median(obj$edge_um[obj$sub_pia == s2], na.rm=TRUE))
pial <- names(edp)[edp < 0.75 * med_all]
cat("  sub-clusters:", paste(subp, table(obj$sub_pia[obj$seurat_clusters == candp])[subp], sep="=", collapse=" "),
    "| median edge distance (um):", paste(names(edp), round(edp), sep="=", collapse=" "), "| tissue-wide median", round(med_all), "-> pial:", paste(pial, collapse=" "), "\n")
stopifnot("no pial sub-cluster found" = length(pial) >= 1)
newid <- as.character(obj$seurat_clusters)
if (length(pial) == length(subp)) {
  cat("  every sub-cluster of", candp, "is pial: the cluster is already a pure pial cluster and is kept whole (curate as L1/pia)\n")
} else {
  for (s2 in pial) { k_new <- as.character(max(as.integer(newid)) + 1); newid[obj$sub_pia == s2] <- k_new; cat("   ", s2, "->", k_new, "(", sum(obj$sub_pia == s2), "spots )\n") }
}
obj$seurat_clusters <- factor(newid, levels = as.character(sort(as.integer(unique(newid)))))
Idents(obj) <- "seurat_clusters"
obj <- RunUMAP(obj, reduction="harmony", dims=1:NDIMS, verbose=FALSE)
cat("  clusters:", nlevels(obj$seurat_clusters), " | ", paste(table(obj$seurat_clusters),collapse="/"), "\n")

saveRDS(obj, file.path(OUT, "NHD_frontal_integrated_harmony.rds"))

## ---- UMAP panels (clusters / condition / sample) ----
pal <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2",
         "#7f7f7f","#bcbd22","#17becf","#aec7e8","#ffbb78","#98df8a","#ff9896")
stopifnot("more clusters than palette entries — extend pal" = nlevels(obj$seurat_clusters) <= length(pal))
u1 <- DimPlot(obj, reduction="umap", group.by="seurat_clusters", label=TRUE, cols=pal) + ggtitle("Clusters") + coord_equal()
u2 <- DimPlot(obj, reduction="umap", group.by="condition", cols=c(CON="#B0BEC5", NHD="#FFB27A")) + ggtitle("Condition") + coord_equal()
u3 <- DimPlot(obj, reduction="umap", group.by="sample_id") + ggtitle("Sample") + coord_equal()
ggsave(file.path(OUT,"UMAP_integrated.png"), (u1|u2|u3), width=16, height=5, dpi=300)
ggsave(file.path(OUT,"UMAP_clusters.png"), u1, width=6.5, height=5.5, dpi=300)

## ---- spatial clusters on tissue (per sample) ----
sp <- SpatialDimPlot(obj, group.by="seurat_clusters", combine=FALSE,
                     pt.size.factor=1.6, image.alpha=0.6, stroke=0)
sp <- lapply(seq_along(sp), function(i) sp[[i]] + scale_fill_manual(values=pal) +
               ggtitle(names(SAMPLES)[i]) + theme(legend.position="none",
                                                  plot.title=element_text(hjust=.5,size=11)))
ggsave(file.path(OUT,"spatial_clusters.png"), wrap_plots(sp, nrow=2) , width=11, height=11, dpi=300)

## cluster composition by condition/sample (like Fig1_supp)
comp <- as.data.frame(prop.table(table(obj$sample_id, obj$seurat_clusters), 1)*100)
colnames(comp) <- c("sample","cluster","pct")
pc <- ggplot(comp, aes(sample, pct, fill=cluster)) + geom_col() + scale_fill_manual(values=pal) +
  theme_bw() + theme(axis.text.x=element_text(angle=45,hjust=1)) + labs(y="% spots", title="Cluster composition")
ggsave(file.path(OUT,"cluster_composition.png"), pc, width=7, height=5, dpi=300)
write.csv(comp, file.path(OUT,"cluster_composition.csv"), row.names=FALSE)

cat("\n=== DONE ===\n", file=stderr())
cat("== integrated outputs in:", OUT, "\n")
