#!/usr/bin/env Rscript
# =============================================================================
# 120_visium_L1_histology_refinement_FH.R — Histology-guided refinement of the Visium domain map.
# Why: the pathologists' layer-I annotation (H&E, drawn blind to the transcriptome) marks the
#   sulcal pial bands of NHD_Frontal1/2. Inside those bands the clustering left a collar of spots labelled WM or L6,
#   labels that are anatomically impossible at a pial surface. Their transcriptome is that of the adjacent L1/pia
#   spots and not of the section's deep white matter (119 / 119b: NHD1 myelin 12 vs 28 per 10k, neuropil 22 vs 34;
#   NHD2 myelin 23 vs 54, neuropil 17 vs 23, GFAP/AQP4 52 vs 29, immune 11 vs 4): sub-pial gliotic layer I that, in
#   this donor, shares an expression state with the gliotic, demyelinated white matter, so the marker-argmax domain
#   call cannot separate them. The correction is therefore made from the histology, by one rule applied identically
#   to all four sections, and is fully reversible (the cluster ids are kept; the moved spots get a NEW cluster id).
#
# RULE (RULE env, default "banksy"):
#   banksy   : PRIMARY. L1/pia = the pial clusters of the neighbourhood-aware BANKSY segmentation run with the SCZ recipe
#              (121; lambda 0.2, k_geom 6, Harmony, res 0.65; pial = >= 30 % of the cluster's spots within 100 um of the pia
#              [tissue edge or L1 line]) restricted to spots inside the pathologists' closed layer-I band or within 300 um of
#              an L1 line (v4). Any current label moves (WM, L6, L2/3, Vasc/immune, ...); validated against the blind
#              annotation: 88 % / 89 % of the NHD 1 / NHD 2 layer-I bands, WM bands unchanged.
#   interior : sensitivity. Only WM/L6-labelled spots inside the closed band between the pathologists' two L1 lines (119b).
#   100um    : as interior or within 100 um of an L1 line.
#   The per-spot table records both the banksy and the interior membership so every downstream number can be recomputed
#   under the sensitivity rule without touching the object.
#
# What it does: reads the integrated object, moves the qualifying spots to a new cluster id (= max id + 1),
#   appends "g<id>" = "L1/pia" to integrated_cluster_identity.csv (module scores NA; note column), writes the object
#   back, and writes tables/visium_L1_refinement_20260918.csv
#   (spot_id, section, previous cluster, previous domain, rule). Re-run 01f afterwards (visium_meta.csv must carry the new id).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
IH <- file.path(dirname(PROJ), "Visium", "integrated_harmony"); TAB <- file.path(PROJ, "tables")
RULE <- Sys.getenv("RULE", "banksy"); stopifnot(RULE %in% c("banksy", "interior", "100um"))
BK <- file.path(dirname(PROJ), "Visium", "banksy_NHD_frontal_20260918", "banksy_per_spot_layers_20260918.csv"); MAX_UM <- 300
BK_COL <- "clust_Harmony_BANKSY_k50_res0.65"; PIAL_FRAC <- 0.30   # a BANKSY cluster is pial when >= 30 % of its spots lie within 100 um of the pia
OBJ <- file.path(IH, "NHD_frontal_integrated_harmony.rds"); IDF <- file.path(IH, "integrated_cluster_identity.csv")
mem <- read.csv(file.path(TAB, "visium_L1band_membership_20260918.csv")); ps <- read.csv(file.path(TAB, "visium_sulcus_check_perspot_FH.csv"))
stopifnot(file.mtime(file.path(TAB, "visium_L1band_membership_20260918.csv")) > file.mtime(OBJ), file.mtime(file.path(TAB, "visium_sulcus_check_perspot_FH.csv")) > file.mtime(OBJ))
obj <- readRDS(OBJ); idf <- read.csv(IDF, check.names = FALSE)
if (!"note" %in% names(idf)) idf$note <- ""
stopifnot(all(paste0("g", levels(obj$seurat_clusters)) %in% idf$cluster), "identity table already carries a refinement row — restore the pre-refinement object and csv first" = !any(grepl("refin", idf$note)))
dom <- unname(setNames(idf$identity, idf$cluster)[paste0("g", as.character(obj$seurat_clusters))])
rownames(mem) <- mem$spot_id; rownames(ps) <- ps$spot_id
inside <- mem[colnames(obj), "inside_L1_band"] == 1; near <- ps[colnames(obj), "d_L1"] <= 100
bk <- read.csv(BK); rownames(bk) <- bk$spot_id; stopifnot(all(colnames(obj) %in% bk$spot_id))
es <- readRDS(file.path(IH, "_spot_edge_src_20260918.rds")); edge <- unname(es$edge_um[colnames(obj)])
d_pia <- pmin(ps[colnames(obj), "d_L1"], edge, na.rm = TRUE)          # for the pial-cluster identification (position of the cluster)
# spot-level position from the pathologists' annotation only: inside the closed layer-I band (119b) or within
# 300 um of an L1 line. The capture-area border is not used (CON1 cut faces / WM-border spots would enter). No transcript gate.
pos_ok <- inside | ps[colnames(obj), "d_L1"] <= MAX_UM
stopifnot("BANKSY table is older than the tissue-call manifest — re-run 121" = file.mtime(BK) > file.mtime(file.path(dirname(PROJ), "Visium", "VISIUM_OUTS_MANIFEST.json")))   # BANKSY must postdate the tissue call
bkc <- as.character(bk[colnames(obj), BK_COL]); stopifnot(!anyNA(bkc))
pial_frac <- tapply(d_pia <= 100, bkc, mean); PIAL_CL <- names(pial_frac)[pial_frac >= PIAL_FRAC]
cat("BANKSY clusters: fraction of spots within 100 um of the pia:\n"); print(round(sort(pial_frac, decreasing = TRUE), 2)); cat("-> pial clusters:", paste(PIAL_CL, collapse = ", "), "\n")
stopifnot("no pial BANKSY cluster found" = length(PIAL_CL) >= 1, "too many pial clusters" = length(PIAL_CL) <= 3)
sel_banksy   <- dom != "L1/pia" & bkc %in% PIAL_CL & pos_ok
sel_interior <- dom %in% c("WM", "L6") & inside
sel_100um    <- dom %in% c("WM", "L6") & (inside | near)
sel <- switch(RULE, banksy = sel_banksy, interior = sel_interior, `100um` = sel_100um)
cat("rule:", RULE, "| spots re-assigned to L1/pia:", sum(sel), "\n"); print(table(section = obj$sample_id[sel], previous = dom[sel]))
k_new <- as.character(max(as.integer(levels(obj$seurat_clusters))) + 1)
prev <- data.frame(spot_id = colnames(obj), section = obj$sample_id, previous_cluster = as.character(obj$seurat_clusters), previous_domain = dom,
                   moved_primary = sel, rule_primary = RULE, banksy_pial = sel_banksy, interior_WM_L6 = sel_interior, within100um_WM_L6 = sel_100um, d_pia_um = round(d_pia))
write.csv(prev, file.path(TAB, "visium_L1_refinement_20260918.csv"), row.names = FALSE)
ts <- format(Sys.time(), "%Y%m%d_%H%M%S"); file.copy(OBJ, file.path(IH, paste0("_pre_refinement_", ts, ".rds")))
newid <- as.character(obj$seurat_clusters); newid[sel] <- k_new
obj$seurat_clusters <- factor(newid, levels = as.character(sort(as.integer(unique(newid))))); Idents(obj) <- "seurat_clusters"
obj$histology_refined_L1 <- unname(sel)
saveRDS(obj, paste0(OBJ, ".tmp")); file.rename(paste0(OBJ, ".tmp"), OBJ)
row <- idf[1, ]; row[] <- NA; row$cluster <- paste0("g", k_new); row$identity <- "L1/pia"
idf <- rbind(idf, row)
idf$note[idf$cluster == paste0("g", k_new)] <- sprintf("L1/pia refinement %s (script 120, rule %s): %s", ts, RULE,
  if (RULE == "banksy") "BANKSY pial cluster (121; SCZ recipe) within 300 um of the pial surface; sensitivity = histology band interior" else "WM/L6-labelled spots inside the pathologists' layer-I band")
write.csv(idf, IDF, row.names = FALSE)
cat("object rewritten:", ncol(obj), "spots,", nlevels(obj$seurat_clusters), "clusters; new cluster", k_new, "= L1/pia (", sum(sel), "spots ). Previous object:", paste0("_pre_refinement_", ts, ".rds"), "\n=== DONE ===\n")
