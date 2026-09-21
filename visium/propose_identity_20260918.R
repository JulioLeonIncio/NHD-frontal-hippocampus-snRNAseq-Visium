#!/usr/bin/env Rscript
# propose_identity_20260918.R — Propose cluster identities for the L1-recovered integration.
# Writes only integrated_cluster_identity_PROPOSED_20260918.csv + crosswalk_L1_20260918.csv +
# cluster_evidence_20260918.csv; the live integrated_cluster_identity.csv is curated by hand from these
# (see IDENTITY_CURATION_20260918.md). refresh_identity.R must still not be run.
#
# Evidence per cluster: (1) module-score argmax with the dictionary extended by an "L1/pia" category
# (leptomeningeal / VLMC genes + AQP4; PTGDS/GFAP/ID3 excluded as WM-shared; the L1 neuronal markers RELN/NDNF/
# CPLX3 are below Visium depth in this run) — the old 11 categories are unchanged; (2) myelin z (so L1/pia is read against
# WM); (4) share of
# recovered spots (tissue_call_source he_mask / expression) per cluster and per section; (5) median distance
# of the cluster's spots to the nearest non-tissue barcode or capture-area edge (um) — L1 must hug the edge;
# (6) per-section composition and NHD fraction.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
OUT <- file.path(dirname(Sys.getenv("NHD_PROJ")), "Visium", "integrated_harmony")
VIS <- dirname(OUT); PROJ <- file.path(dirname(VIS), "NHD_frontal_hippo_rebuild")
source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))
stopifnot(!file.exists(file.path(OUT, "integrated_cluster_identity_PROPOSED_20260918.csv")) ||
          Sys.getenv("FORCE") == "1")
obj <- readRDS(file.path(OUT, "NHD_frontal_integrated_harmony.rds")); Idents(obj) <- "seurat_clusters"
cat("object:", ncol(obj), "spots,", nlevels(obj$seurat_clusters), "clusters\n")

MARK <- list(
  # L1/pia = leptomeningeal / VLMC genes (Yang 2022 PMID 35165441; Kearns 2023 PMID 37923721) + AQP4 (L1-enriched in
  # DLPFC Visium, Maynard 2021 PMID 33558695). PTGDS / GFAP / ID3 are deliberately not here: they are shared with
  # (gliotic, demyelinated) white matter and would separate L1 from L2/3 but not from WM.
  "L1/pia"=c("COL1A2","COL1A1","COL3A1","SLC6A13","CEMIP","DCN","LUM","LAMA2","SPARC","HTRA1","AQP4"),
  "L2/3"=c("CUX2","LAMP5"), "L4"=c("RORB","NEFH"), "L5"=c("PCP4","FEZF2","BCL11B"),
  "L6"=c("TLE4","FOXP2"), "ExN"=c("RBFOX3","SLC17A7","SNAP25"), "InN"=c("GAD1","GAD2"),
  "Oligo/WM"=c("MBP","PLP1","MOBP","MOG"), "OPC"=c("PDGFRA","OLIG1"),
  "Astro"=c("GFAP","AQP4","SLC1A2"), "Micro"=c("TYROBP","C1QB","CSF1R","P2RY12"),
  "Vasc"=c("CLDN5","PECAM1","PDGFRB"))
feat <- intersect(unlist(MARK), rownames(obj))
avg <- AverageExpression(obj, assays="Spatial", features=feat, group.by="seurat_clusters")[[1]]
z <- t(scale(t(as.matrix(avg))))
catz <- t(sapply(names(MARK), function(k){ g<-intersect(MARK[[k]],rownames(z)); colMeans(z[g,,drop=FALSE]) }))
ident <- rownames(catz)[apply(catz,2,which.max)]
clk <- colnames(catz)
ann <- data.frame(cluster=clk, argmax=ident, t(round(catz,2)), check.names=FALSE)

old <- readRDS(file.path(OUT, "_pre_L1_20260918", "NHD_frontal_integrated_harmony.rds"))
oid <- read.csv(file.path(OUT, "_pre_L1_20260918", "integrated_cluster_identity.csv"))
olab <- setNames(oid$identity, sub("^g","",oid$cluster))[as.character(old$seurat_clusters)]; names(olab) <- colnames(old)
shared <- intersect(colnames(obj), colnames(old)); cat("shared spots with the 2026-09-16 object:", length(shared), "of", ncol(obj), "\n")
cw <- table(new = as.character(obj$seurat_clusters[shared]), old = olab[shared])
write.csv(as.data.frame.matrix(cw), file.path(OUT, "crosswalk_L1_20260918.csv"))
top_old <- apply(cw, 1, function(x) names(x)[which.max(x)]); purity <- round(apply(cw, 1, max) / rowSums(cw), 2)

# ---- recovered-spot share, edge distance, composition ----------------------------------------
tcs <- do.call(rbind, lapply(VIS_SECTIONS, function(nm) { t <- read.csv(file.path(vis_outs(nm), "tissue_call.csv")); t$cell <- paste0(nm, "_", t$barcode); t$section <- nm; t }))
rownames(tcs) <- tcs$cell
src <- tcs[colnames(obj), "tissue_call_source"]
# Distance (um) to the nearest barcode that is not in tissue, or to the capture-area edge (array rows 0/77, cols 0/127)
sf <- lapply(VIS_SECTIONS, function(nm) jsonlite::fromJSON(file.path(vis_outs(nm), "spatial", "scalefactors_json.json"))); names(sf) <- VIS_SECTIONS
edge_um <- rep(NA_real_, ncol(obj)); names(edge_um) <- colnames(obj)
for (nm in VIS_SECTIONS) {
  t <- tcs[tcs$section == nm, ]; px_per_um <- sf[[nm]]$spot_diameter_fullres / 55
  inn <- t[t$in_tissue_analysis == 1, ]; outp <- t[t$in_tissue_analysis == 0, ]
  cells <- intersect(inn$cell, colnames(obj)); inn <- inn[cells, ]
  d_out <- if (nrow(outp)) { D <- sqrt(outer(inn$pxl_row_in_fullres, outp$pxl_row_in_fullres, "-")^2 + outer(inn$pxl_col_in_fullres, outp$pxl_col_in_fullres, "-")^2); apply(D, 1, min) / px_per_um } else rep(Inf, nrow(inn))
  # capture-area edge: array_row 0..77 (100 um pitch along rows = 2 array cols), array_col 0..127
  d_edge <- pmin(pmin(inn$array_row, 77 - inn$array_row) * 86.6,      # hex lattice: rows 86.6 um apart,
                 pmin(inn$array_col, 127 - inn$array_col) * 50)        # staggered cols 50 um apart (100 um centre-to-centre)
  edge_um[cells] <- pmin(d_out, d_edge)
}
obj$edge_um <- edge_um; obj$src <- src
md <- obj@meta.data
ev <- do.call(rbind, lapply(levels(obj$seurat_clusters), function(k) {
  m <- md[md$seurat_clusters == k, ]
  data.frame(cluster = paste0("g", k), n = nrow(m), n_CON = sum(m$condition == "CON"), n_NHD = sum(m$condition == "NHD"),
             frac_NHD = round(mean(m$condition == "NHD"), 3),
             CON1 = sum(m$sample_id == "CON_Frontal1"), CON2 = sum(m$sample_id == "CON_Frontal2"), NHD1 = sum(m$sample_id == "NHD_Frontal1"), NHD2 = sum(m$sample_id == "NHD_Frontal2"),
             n_recovered = sum(m$src %in% c("he_mask","expression")), frac_recovered = round(mean(m$src %in% c("he_mask","expression")), 3),
             share_of_all_recovered = round(sum(m$src %in% c("he_mask","expression")) / sum(md$src %in% c("he_mask","expression")), 3),
             median_edge_um = round(median(m$edge_um, na.rm = TRUE)), frac_within_200um = round(mean(m$edge_um <= 200, na.rm = TRUE), 2),
             median_UMI = round(median(m$nCount_Spatial)), median_pct_mt = round(median(m$percent.mt), 1),
             top_old_label = top_old[k], purity_vs_old = purity[k], n_shared = rowSums(cw)[k])
}))
ann <- merge(ann, ev, by = "cluster"); ann <- ann[order(as.integer(sub("^g","",ann$cluster))), ]
ann$myelin_z <- ann[["Oligo/WM"]]
# proposal: argmax, with the L1/pia call accepted only where the cluster is myelin-low and edge-hugging
ann$proposed <- ann$argmax
ann$proposed[ann$argmax == "L1/pia" & !(ann$myelin_z < 0 & ann$frac_within_200um >= 0.5)] <- NA
ann$proposed[ann$argmax == "Oligo/WM"] <- "WM"; ann$proposed[ann$argmax == "Vasc"] <- "Vasc/immune"; ann$proposed[ann$argmax %in% c("Astro","Micro","OPC","ExN")] <- NA
ann$note <- ifelse(is.na(ann$proposed), "CURATE BY HAND", "")
write.csv(ann, file.path(OUT, "cluster_evidence_20260918.csv"), row.names = FALSE)
write.csv(ann[, c("cluster","proposed","argmax","top_old_label","purity_vs_old","n","frac_NHD","frac_recovered","share_of_all_recovered","median_edge_um","frac_within_200um","myelin_z","note")],
          file.path(OUT, "integrated_cluster_identity_PROPOSED_20260918.csv"), row.names = FALSE)
cat("\n== proposal ==\n"); print(ann[, c("cluster","proposed","argmax","top_old_label","purity_vs_old","n","n_CON","n_NHD","frac_recovered","share_of_all_recovered","median_edge_um","frac_within_200um","myelin_z","median_UMI")], row.names = FALSE)
cat("\n== module scores (category x cluster) ==\n"); print(round(catz, 2))
cat("\n== top-2 categories per cluster ==\n")
for (j in seq_len(ncol(catz))) { x <- catz[, j]; o <- order(x, decreasing = TRUE); cat(sprintf("  %s: %s(%.2f) > %s(%.2f)\n", colnames(catz)[j], rownames(catz)[o[1]], x[o[1]], rownames(catz)[o[2]], x[o[2]])) }
cat("\n== where did the recovered spots go (source x cluster) ==\n"); print(table(md$src, md$seurat_clusters))
cat("\n== per-section spot counts ==\n"); print(table(md$sample_id))
saveRDS(list(edge_um = edge_um, src = src), file.path(OUT, "_spot_edge_src_20260918.rds"))
cat("\n=== DONE ===\n")
