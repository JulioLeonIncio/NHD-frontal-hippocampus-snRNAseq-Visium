#!/usr/bin/env Rscript
# refresh_identity.R — Lightweight refresh of cluster identity from the current integrated object (identity CSV was Jun 10).
# Mirrors annotate_clusters.R identity logic (module-score argmax) without the
# CPU-heavy FindAllMarkers, so it won't disturb the running cell2location job.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
suppressPackageStartupMessages({ library(Seurat) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
OUT <- file.path(dirname(Sys.getenv("NHD_PROJ")), "Visium", "integrated_harmony")
obj <- readRDS(file.path(OUT,"NHD_frontal_integrated_harmony.rds"))
Idents(obj) <- "seurat_clusters"
MARK <- list(
  "L2/3"=c("CUX2","LAMP5"), "L4"=c("RORB","NEFH"), "L5"=c("PCP4","FEZF2","BCL11B"),
  "L6"=c("TLE4","FOXP2"), "ExN"=c("RBFOX3","SLC17A7","SNAP25"), "InN"=c("GAD1","GAD2"),
  "Oligo/WM"=c("MBP","PLP1","MOBP","MOG"), "OPC"=c("PDGFRA","OLIG1"),
  "Astro"=c("GFAP","AQP4","SLC1A2"), "Micro"=c("TYROBP","C1QB","CSF1R","P2RY12"),
  "Vasc"=c("CLDN5","PECAM1","PDGFRB"))
feat <- intersect(unlist(MARK), rownames(obj))
avg <- AverageExpression(obj, assays="Spatial", features=feat, group.by="seurat_clusters")[[1]]
z <- t(scale(t(as.matrix(avg))))
catz <- t(sapply(names(MARK), function(k){ g<-intersect(MARK[[k]],rownames(z))
  if(length(g)) colMeans(z[g,,drop=FALSE]) else rep(NA_real_,ncol(z)) }))
ident <- rownames(catz)[apply(catz,2,function(x) which.max(replace(x, is.na(x), -Inf)))]
clk <- colnames(catz)                                  # already "g0".."g8" (AverageExpression prefix)
ann <- data.frame(cluster=clk, identity=ident, t(round(catz,2)), check.names=FALSE)
write.csv(ann, file.path(OUT,"integrated_cluster_identity.csv"), row.names=FALSE)
comp <- as.data.frame(table(obj$seurat_clusters, obj$condition)); colnames(comp)<-c("cluster","condition","n")
comp$identity <- ann$identity[match(paste0("g",comp$cluster), ann$cluster)]
write.csv(comp, file.path(OUT,"integrated_cluster_condition_counts.csv"), row.names=FALSE)
cat("== cluster -> identity (refreshed from current object) ==\n"); print(ann[,c("cluster","identity")])
cat("\n== full module-score matrix (category x cluster) ==\n"); print(round(catz,2))
cat("\n== top-2 categories per cluster (margin check) ==\n")
for(j in seq_len(ncol(catz))){ x<-catz[,j]; o<-order(x,decreasing=TRUE)
  cat(sprintf("  %s: %s(%.2f) > %s(%.2f)\n", colnames(catz)[j], rownames(catz)[o[1]], x[o[1]], rownames(catz)[o[2]], x[o[2]])) }
cat("\n=== DONE ===\n", file=stderr())
