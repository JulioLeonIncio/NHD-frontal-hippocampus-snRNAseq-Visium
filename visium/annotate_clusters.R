#!/usr/bin/env Rscript
# annotate_clusters.R — (b) Annotate the 9 integrated Harmony clusters: DGE per cluster + canonical marker dotplot (SCZ Fig1_supp-d style) + identity assignment.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
OUT <- file.path(dirname(Sys.getenv("NHD_PROJ")), "Visium", "integrated_harmony")
obj <- readRDS(file.path(OUT,"NHD_frontal_integrated_harmony.rds"))
Idents(obj) <- "seurat_clusters"

cat("== FindAllMarkers ==\n")
mk <- FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.5, verbose=FALSE)
write.csv(mk, file.path(OUT,"integrated_cluster_DGE.csv"), row.names=FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n=8) %>% ungroup()
write.csv(top, file.path(OUT,"integrated_cluster_top8.csv"), row.names=FALSE)

MARK <- list(
  "L2/3"=c("CUX2","LAMP5"), "L4"=c("RORB","NEFH"), "L5"=c("PCP4","FEZF2","BCL11B"),
  "L6"=c("TLE4","FOXP2"), "ExN"=c("RBFOX3","SLC17A7","SNAP25"), "InN"=c("GAD1","GAD2"),
  "Oligo/WM"=c("MBP","PLP1","MOBP","MOG"), "OPC"=c("PDGFRA","OLIG1"),
  "Astro"=c("GFAP","AQP4","SLC1A2"), "Micro"=c("TYROBP","C1QB","CSF1R","P2RY12"),
  "Vasc"=c("CLDN5","PECAM1","PDGFRB"))
feat <- intersect(unlist(MARK), rownames(obj))

dp <- DotPlot(obj, features=feat) + RotatedAxis() +
  theme(axis.text.x=element_text(size=7,face="italic")) + ggtitle("Integrated clusters x canonical markers")
ggsave(file.path(OUT,"integrated_cluster_dotplot.png"), dp, width=14, height=5, dpi=300)
dp2 <- DotPlot(obj, features=unique(top$gene)) + RotatedAxis() +
  theme(axis.text.x=element_text(size=6,face="italic")) + ggtitle("Top DE markers per integrated cluster")
ggsave(file.path(OUT,"integrated_topDE_dotplot.png"), dp2, width=18, height=5, dpi=300)

## identity assignment: z-score avg expr across clusters, mean per category, argmax
avg <- AverageExpression(obj, assays="Spatial", features=feat, group.by="seurat_clusters")[[1]]
z <- t(scale(t(as.matrix(avg))))           # gene x cluster
catz <- t(sapply(names(MARK), function(k){ g<-intersect(MARK[[k]],rownames(z)); if(length(g)) colMeans(z[g,,drop=FALSE]) else rep(NA_real_,ncol(z)) }))
ident <- rownames(catz)[apply(catz,2,function(x) which.max(replace(x, is.na(x), -Inf)))]
ann <- data.frame(cluster=colnames(catz), identity=ident, t(round(catz,2)), check.names=FALSE)
write.csv(ann, file.path(OUT,"integrated_cluster_identity.csv"), row.names=FALSE)
cat("== cluster -> identity ==\n"); print(ann[,c("cluster","identity")])

## composition of identities by condition
comp <- as.data.frame(table(obj$seurat_clusters, obj$condition)); colnames(comp)<-c("cluster","condition","n")
comp$identity <- ann$identity[match(comp$cluster, ann$cluster)]
write.csv(comp, file.path(OUT,"integrated_cluster_condition_counts.csv"), row.names=FALSE)
cat("\n=== DONE ===\n", file=stderr()); cat("== annotation outputs in:", OUT, "\n")
