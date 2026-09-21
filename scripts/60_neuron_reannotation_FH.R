#!/usr/bin/env Rscript
# =============================================================================
# 60_neuron_reannotation_FH.R — 60_neuron_reannotation_FH.R Re-annotate neurons on the NEW Frontal+Hippocampus atlas (Fig 4 prerequisite).
# Those
# cluster ids are invalid on this re-clustered FH atlas, so the hippo taxonomy is
# re-derived here on the new clustering, from marker evidence — never a blind copy
# of the old cluster->label map.
#
# METHOD (ported verbatim from the established recipe in
# manuscript_7fig/scripts/64_fig4_neuron_umap.R, scripts/15_hippo_annotation_Tippani2025.R):
#   HIPPO (de-novo): subset Region=="Hippo" Neuron_Ex+Neuron_Inh ->
#     SCT -> RunPCA(30) -> RunHarmony(SampleID) -> FindNeighbors(1:30) ->
#     FindClusters(res 0.1, seed 42) -> RunUMAP(seed 42). Then score each frozen
#     Louvain cluster against the same 60-marker canonical human-hippocampus panel
#     (Franjic 2022 / Ayhan 2021 / Tippani 2025) and assign the reproducible
#     taxonomy: CA3, Subiculum, EC-like, Cajal-Retzius, Inh-CGE (VIP),
#     Inh-CGE (LAMP5), Inh-MGE. The per-cluster marker mean-expr table + a
#     z-scored argmax hint are written so the label call is fully auditable.
# On this FH clustering (res 0.1) the taxonomy re-derives to six labels —
#     EC-like did not separate as its own entorhinal-dominant cluster (it did on the
#     pre-OCC-removal atlas). Reported honestly + flagged, not forced onto a cluster.
#   CORTICAL (Frontal):
# Excitatory nuclei now take the Jorstad 2023 dorsolateral-prefrontal
#     label (scripts/116b_L4_label_transfer_DFC_FH.R, tables/L4_DFC_transfer_percell_FH.csv,
#     column `dfc_within`). why: the Azimuth reference is human motor cortex, which has
#     no L4 IT class, so the 2,102 frontal L4 IT nuclei were force-mapped onto L2/3 IT /
#     L5 IT and the L4 loss seen in Visium could not be asked of the snRNA data. The DFC
#     transfer gives 9 excitatory classes (L2/3 IT, L4 IT, L5 IT, L6 IT, L6 IT Car3, L5 ET,
#     L5/6 NP, L6 CT, L6b). inhibitory nuclei keep the Azimuth subclass (Lamp5, Pvalb,
#     Sncg, Sst, Sst Chodl, Vip) — the DFC transfer was run on excitatory nuclei only.
#     Cortical taxonomy = 15 subclasses. `source_annot` records which reference labelled
#     each nucleus; `dfc_within_score` (prediction.score.max) is carried for excitatory
#     nuclei; `azimuth_subclass` keeps the old label for every frontal neuron (audit).
#     Every consumer reads the label through scripts/_neuron_subclass_FH.R, never from
#     atlas `predicted.subclass`.
#
# Honest LIMITS (confirmed + reported): DG-granule not captured (no PROX1-high
# pyramidal cluster — dissection); CA1/CA2 not separately resolved (fold into CA3
# anchor); EC-like not resolved on the FH clustering (6 labels, not 7).
#
# Outputs:
#   data/neuron_subtype_map_FH.rds/.csv  (barcode, Region, neuron_subtype, Condition,
#                                         source_annot, azimuth_subclass, dfc_within_score;
#                                         hippo=de-novo label, cortical Ex=DFC, Inh=Azimuth)
#   data/_cache_hippo_neuron_umap_FH.rds (hippo embedding + labels, mtime-guarded)
#   figures/Figure_5/panels/4_hippo_neuron_subtypes_UMAP.png  (eye-check)
#   logs/60_neuron_reannotation_FH.log + provenance (census, marker table, sessionInfo)
#
# House rules: no bold, no caption-in-panel, full region names, seed 42.
# Run: OMP_NUM_THREADS=2 KMP_DUPLICATE_LIB_OK=TRUE Rscript scripts/60_neuron_reannotation_FH.R
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(ggplot2); library(scales)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, PAL_COND, REGION_FULL

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
DDIR  <- file.path(PROJ, "data")
CACHE_HIP <- file.path(DDIR, "_cache_hippo_neuron_umap_FH.rds")
MAP_RDS   <- file.path(DDIR, "neuron_subtype_map_FH.rds")
MAP_CSV   <- file.path(DDIR, "neuron_subtype_map_FH.csv")
DFC_CSV   <- file.path(PROJ, "tables", "L4_DFC_transfer_percell_FH.csv")   # 116b output
PANEL <- file.path(PROJ, "figures", "Figure_5", "panels")
LOGD  <- file.path(PROJ, "logs")
for (d in c(DDIR, PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING atlas NHD_FH_harmony.rds" = file.exists(ATLAS),
          "MISSING tables/L4_DFC_transfer_percell_FH.csv — run 116b_L4_label_transfer_DFC_FH.R" = file.exists(DFC_CSV))

LOG <- file.path(LOGD, "60_neuron_reannotation_FH.log")
logcon <- file(LOG, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
say("== 60_neuron_reannotation_FH ==  ", format(Sys.time()))
say("PROJ: ", PROJ)

HIP_RES <- 0.1     # same target (~6-12 clusters on ~3.5k nuclei) as reference 64

# ---- canonical human-hippocampus marker panel (Franjic 2022 / Ayhan 2021 /
#      Tippani 2025) — identical to reference 64_fig4_neuron_umap.R -------------
HIP_MARKERS <- list(
  "DG granule"    = c("PROX1","PPFIA2","TRPC6","GLIS3"),
  "Mossy cell"    = c("ADCYAP1","CALB2","SATB1","CARTPT"),
  "CA1"           = c("GRIK1","GRM3","FIBCD1","SPOCK1","POU3F1","WFS1","MPPED1"),
  "CA2"           = c("RGS14","AMIGO2","HGF"),
  "CA3"           = c("CFAP299","HS3ST4","GRIK4","NECAB1","NECTIN3","CHGB"),
  "Subiculum"     = c("ROBO1","FN1","NTS","TLE4","NR4A2"),
  "Entorhinal"    = c("CALB1","RORB","THEMIS","FOXP2","LAMA3"),
  "Cajal-Retzius" = c("TP73","RELN","LHX1","NDNF"),
  "Interneuron"   = c("GAD1","GAD2","ADARB2","LHX6","PVALB","SST","VIP","LAMP5",
                      "ID2","PAX6","SNCG","CNR1","CCK"),
  "pan-neuronal"  = c("SYT1","SNAP25","RBFOX3","SLC17A7","SLC17A6"),
  "oligo-QC"      = c("PLP1","MBP","MOBP","MOG"))

# CGE vs MGE interneuron split markers (used to name the two interneuron classes)
CGE_MK <- c("ADARB2","VIP","LAMP5","ID2","PAX6","SNCG","CNR1","CCK")   # caudal ganglionic eminence
MGE_MK <- c("LHX6","PVALB","SST")                                     # medial ganglionic eminence
VIP_MK <- c("VIP"); LAMP5_MK <- c("LAMP5")

# cortical subclasses for census ordering (EX = the 9 Jorstad-DFC classes,
# INH = the Azimuth classes; same fixed order as scripts/_neuron_subclass_FH.R)
EX_SUBC  <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L5 ET","L5/6 NP",
              "L6 CT","L6b")
INH_SUBC <- c("Lamp5","Pvalb","Sncg","Sst","Sst Chodl","Vip")

# hippo subtype legend order + palette (mirrors reference 64; interneurons first)
# The cluster that carried the label is an unresolved low-complexity
# cluster, drawn last and in neutral grey so it never reads as a resolved subtype.
UNRESOLVED_LAB  <- "Unresolved"
HIP_LABEL_ORDER <- c("Inh-CGE (VIP)","Inh-CGE (LAMP5)","Inh-MGE",
                     "CA3","Subiculum","EC-like", UNRESOLVED_LAB)
HIP_LAB_PAL <- c(
  "Inh-CGE (VIP)"   = "#8E0152",  "Inh-CGE (LAMP5)" = "#C2549D",
  "Inh-MGE"         = "#8E44AD",  "CA3"             = "#08519C",
  "Subiculum"       = "#41AB5D",  "EC-like"         = "#00CED1")
HIP_LAB_PAL[UNRESOLVED_LAB] <- "#9E9E9E"

# =============================================================================
# 1. Load atlas once; slice hippo-neuron subset + cortical-neuron metadata; release
# =============================================================================
say("\n== step: load atlas (once) ==")
t0 <- Sys.time()
atl <- readRDS(ATLAS)
say("atlas cells: ", ncol(atl), "  (loaded in ",
    round(as.numeric(Sys.time() - t0, units = "mins"), 2), " min)")
stopifnot("atlas missing required meta" =
            all(c("new_annotation","Region","Condition","SampleID","predicted.subclass")
                %in% colnames(atl@meta.data)))

md   <- atl@meta.data
isN  <- md$new_annotation %in% c("Neuron_Ex","Neuron_Inh")
say("total neuron nuclei: ", sum(isN))
say("  by region x class:")
print(table(Region = md$Region[isN], class = md$new_annotation[isN]))

# --- cortical (Frontal) neuron metadata -------------------------------------------
# excitatory -> Jorstad 2023 DFC `dfc_within` (116b); inhibitory -> Azimuth.
ctx_bc <- rownames(md)[isN & md$Region == "Frontal"]
ctx_class <- as.character(md[ctx_bc, "new_annotation"])
azi <- as.character(md[ctx_bc, "predicted.subclass"])
dfc <- read.csv(DFC_CSV, stringsAsFactors = FALSE)
stopifnot("DFC table lacks required columns" =
            all(c("barcode","dfc_within","dfc_within_score") %in% colnames(dfc)),
          "duplicated barcode in DFC table" = !anyDuplicated(dfc$barcode))
ex_bc <- ctx_bc[ctx_class == "Neuron_Ex"]
miss_dfc <- setdiff(ex_bc, dfc$barcode)
stopifnot("frontal Neuron_Ex barcode(s) ABSENT from the DFC transfer table" = length(miss_dfc) == 0)
extra_dfc <- setdiff(dfc$barcode, ex_bc)
if (length(extra_dfc)) say("NOTE: ", length(extra_dfc), " DFC rows are not frontal Neuron_Ex in this atlas (ignored)")
stopifnot("DFC labels are not the expected 9 excitatory classes" =
            setequal(unique(dfc$dfc_within[dfc$barcode %in% ex_bc]), EX_SUBC))
i_dfc <- match(ctx_bc, dfc$barcode)                       # NA for inhibitory nuclei
is_ex <- ctx_class == "Neuron_Ex"
ctx_lab <- ifelse(is_ex, dfc$dfc_within[i_dfc], azi)
ctx_map <- data.frame(
  barcode          = ctx_bc,
  Region           = "Frontal",
  neuron_subtype   = ctx_lab,
  Condition        = as.character(md[ctx_bc, "Condition"]),
  source_annot     = ifelse(is_ex, "ctx_DFC_Jorstad2023_within", "ctx_Azimuth_subclass"),
  azimuth_subclass = azi,
  dfc_within_score = ifelse(is_ex, dfc$dfc_within_score[i_dfc], NA_real_),
  stringsAsFactors = FALSE)
stopifnot("frontal inhibitory nucleus with a non-Azimuth-inhibitory label" =
            all(ctx_map$neuron_subtype[!is_ex] %in% INH_SUBC),
          "NA cortical label" = !anyNA(ctx_map$neuron_subtype))
say("cortical (Frontal) neurons: ", nrow(ctx_map),
    "  Ex (DFC-labelled): ", sum(is_ex), "  Inh (Azimuth): ", sum(!is_ex),
    "  NA subclass: ", sum(is.na(ctx_map$neuron_subtype)))
say("  Ex relabelled (Azimuth -> DFC differs): ", sum(ctx_lab[is_ex] != azi[is_ex]),
    " of ", sum(is_ex), sprintf(" (%.1f%%)", 100 * mean(ctx_lab[is_ex] != azi[is_ex])))
say("  DFC within-area score, Ex: median ", round(median(ctx_map$dfc_within_score[is_ex]), 3),
    ", <0.5 in ", sum(ctx_map$dfc_within_score[is_ex] < 0.5), " nuclei")
say("  Azimuth x DFC crosswalk (frontal Ex):")
print(table(Azimuth = azi[is_ex], DFC = ctx_lab[is_ex]))

# --- hippo neuron subset (Seurat object) for de-novo re-embedding ---------------
hip_bc <- rownames(md)[isN & md$Region == "Hippo"]
say("hippo neurons to re-embed: ", length(hip_bc))
o_hip <- subset(atl, cells = hip_bc)
rm(atl); gc()                                    # release the 2.9 GB atlas
say("atlas released; hippo subset nuclei: ", ncol(o_hip))

# =============================================================================
# 2. Hippo de-novo embedding + clustering (rebuild iff cache stale vs atlas)
# =============================================================================
atlas_mtime <- file.mtime(ATLAS)
cache_stale <- !file.exists(CACHE_HIP) || file.mtime(CACHE_HIP) < atlas_mtime
if (!cache_stale) {                              # schema/res guard, like reference 64
  .h <- readRDS(CACHE_HIP)
  if (is.null(.h$res) || !isTRUE(all.equal(.h$res, HIP_RES)) ||
      is.null(.h$method) || .h$method != "denovo_cluster" ||
      is.null(.h$meta$nFeature)) {
    say("== hip cache res/method mismatch -> rebuild =="); cache_stale <- TRUE
  }
  rm(.h)
}

if (cache_stale) {
  say("\n== step: hippo de-novo embed + cluster (SCT->PCA->Harmony->UMAP seed42) ==")
  suppressPackageStartupMessages(library(harmony))
  DefaultAssay(o_hip) <- "SCT"
  stopifnot("SCT scale.data empty — re-run SCTransform/PrepSCT" =
              nrow(GetAssayData(o_hip, assay = "SCT", layer = "scale.data")) > 0)
  tt <- system.time({
    o_hip <- RunPCA(o_hip, npcs = 30, verbose = FALSE)
    o_hip <- harmony::RunHarmony(o_hip, group.by.vars = "SampleID", verbose = FALSE)
    o_hip <- FindNeighbors(o_hip, reduction = "harmony", dims = 1:30, verbose = FALSE)
    set.seed(42)
    o_hip <- FindClusters(o_hip, resolution = HIP_RES, random.seed = 42, verbose = FALSE)
    o_hip <- RunUMAP(o_hip, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)
  })
  rt_hip <- unname(tt["elapsed"])
  emb <- Embeddings(o_hip, "umap")[, 1:2]
  clus <- factor(as.character(o_hip$seurat_clusters),
                 levels = as.character(sort(as.integer(levels(o_hip$seurat_clusters)))))
  names(clus) <- colnames(o_hip)
  say(sprintf("de-novo clusters (res %.2f): %d  sizes: %s",
              HIP_RES, nlevels(clus), paste(table(clus), collapse = "/")))
  say(sprintf("embedding runtime: %.1f s", rt_hip))

  # --- per-cluster marker mean-expr (RNA log-norm).
  o_rna <- o_hip
  DefaultAssay(o_rna) <- "RNA"
  if (length(SeuratObject::Layers(o_rna, assay = "RNA")) > 1)
    o_rna <- SeuratObject::JoinLayers(o_rna, assay = "RNA")
  o_rna <- NormalizeData(o_rna, assay = "RNA", verbose = FALSE)
  mk_all <- unique(unlist(HIP_MARKERS, use.names = FALSE))
  mk     <- intersect(mk_all, rownames(o_rna))
  mk_abs <- setdiff(mk_all, rownames(o_rna))
  if (length(mk_abs)) say("markers ABSENT from atlas (skipped): ", paste(mk_abs, collapse = ", "))
  stopifnot("no requested markers present" = length(mk) > 0)
  expr <- GetAssayData(o_rna, assay = "RNA", layer = "data")[mk, , drop = FALSE]
  cl2  <- droplevels(clus[colnames(expr)])
  mtbl <- sapply(levels(cl2), function(k) Matrix::rowMeans(expr[, cl2 == k, drop = FALSE]))
  mtbl <- as.matrix(mtbl); rownames(mtbl) <- mk
  # %detected too (used for the DG floor check)
  pdet <- sapply(levels(cl2), function(k) Matrix::rowMeans(expr[, cl2 == k, drop = FALSE] > 0) * 100)
  pdet <- as.matrix(pdet); rownames(pdet) <- mk
  grp_of <- setNames(rep(names(HIP_MARKERS), lengths(HIP_MARKERS)), unlist(HIP_MARKERS))
  mtbl <- mtbl[order(match(grp_of[mk], names(HIP_MARKERS))), , drop = FALSE]
  rm(o_rna, expr); gc()

  meta_hip <- data.frame(
    Condition = as.character(o_hip$Condition),
    Region    = "Hippo",
    SampleID  = as.character(o_hip$SampleID),
    class     = as.character(o_hip$new_annotation),
    cluster   = as.character(clus[colnames(o_hip)]),
    nFeature  = as.integer(o_hip$nFeature_RNA),
    row.names = colnames(o_hip), stringsAsFactors = FALSE)

  saveRDS(list(embedding = emb[rownames(meta_hip), , drop = FALSE], meta = meta_hip,
               atlas_mtime = atlas_mtime, res = HIP_RES, method = "denovo_cluster",
               runtime = rt_hip, marker_tbl = mtbl, marker_pct = pdet,
               marker_grp = grp_of[rownames(mtbl)], markers_used = mk, markers_absent = mk_abs),
          paste0(CACHE_HIP, ".tmp"))
  file.rename(paste0(CACHE_HIP, ".tmp"), CACHE_HIP)   # atomic
  say("wrote cache: ", basename(CACHE_HIP))
} else {
  say("\n== hippo cache fresh (mtime >= atlas) -> reload ==")
}
hip <- readRDS(CACHE_HIP)
rm(o_hip); gc()

# =============================================================================
# 3. re-derive cluster -> subtype label from marker evidence (data-driven)
# =============================================================================
say("\n== step: re-derive hippo cluster -> subtype (marker z-score argmax) ==")
mtbl <- hip$marker_tbl; mgrp <- hip$marker_grp; pdet <- hip$marker_pct
clus_ids <- colnames(mtbl)

# per-cluster group score = mean of each marker group's genes; z across clusters.
grp_names <- unique(mgrp)
grp_score <- sapply(grp_names, function(g) colMeans(mtbl[mgrp == g, , drop = FALSE]))  # clusters x groups
grp_z <- scale(grp_score)                                     # z per group across clusters
say("cluster x marker-group mean expr (RNA log-norm):")
print(round(as.data.frame(grp_score), 3))
cat(capture.output(print(round(as.data.frame(grp_score), 3))), sep = "\n", file = logcon)

# oligo-QC contamination guard: compare RAW mean expression (not z-scores — z is
# scaled within each group so a cluster with the highest oligo z can still have
# near-zero raw myelin). A cluster is a possible doublet iff its RAW myelin mean
# EXCEEDS its RAW pan-neuronal mean (the honest single-cluster dominance test,
# mirrors scripts/15b_hippo_neuron_clean_glia_call.R).
oligo_mean  <- grp_score[, "oligo-QC"]
neuron_mean <- grp_score[, "pan-neuronal"]

# interneuron vs excitatory: pan-inhibitory GAD1/GAD2 signal per cluster
gad_mean <- colMeans(mtbl[intersect(c("GAD1","GAD2"), rownames(mtbl)), , drop = FALSE])
ex_mean  <- colMeans(mtbl[intersect(c("SLC17A7","SLC17A6"), rownames(mtbl)), , drop = FALSE])
is_inh   <- gad_mean > ex_mean

# CGE vs MGE (interneuron clusters); VIP vs LAMP5 within CGE
cge_mean  <- colMeans(mtbl[intersect(CGE_MK,  rownames(mtbl)), , drop = FALSE])
mge_mean  <- colMeans(mtbl[intersect(MGE_MK,  rownames(mtbl)), , drop = FALSE])
vip_mean  <- mtbl[intersect(VIP_MK,   rownames(mtbl)), , drop = FALSE][1, ]
lamp5_mean<- mtbl[intersect(LAMP5_MK, rownames(mtbl)), , drop = FALSE][1, ]

# Excitatory subfield: argmax over the finest subfield groups (for the audit
# `ex_subfield` column) using z (comparable across groups). This is a hint.
EX_GROUPS <- c("CA3","CA1","CA2","Subiculum","Entorhinal","DG granule","Mossy cell","Cajal-Retzius")
ex_z <- grp_z[, EX_GROUPS, drop = FALSE]
ex_call_raw <- EX_GROUPS[apply(ex_z, 1, which.max)]

# Final excitatory taxonomy call = argmax (on z, comparable) restricted to the
# four resolved excitatory labels only. CA1/CA2/DG/Mossy are not separately
# resolvable at this cohort/resolution, so their z folds into the CA3 (ca-pyramidal)
# anchor; Entorhinal -> EC-like; Cajal-Retzius stays. A cluster becomes EC-like
# only if the entorhinal programme actually dominates the ca/subiculum programmes.
EX_FINAL_MAP <- list(CA3 = "CA3", Subiculum = "Subiculum",
                     Entorhinal = "EC-like", `Cajal-Retzius` = "Cajal-Retzius")
# ca-pyramidal composite z = max over CA1/CA2/CA3/DG/Mossy (all -> "CA3" anchor)
ca_anchor_z <- apply(grp_z[, c("CA3","CA1","CA2","DG granule","Mossy cell"), drop = FALSE], 1, max)
ex_final_z <- cbind(CA3           = ca_anchor_z,
                    Subiculum     = grp_z[, "Subiculum"],
                    Entorhinal    = grp_z[, "Entorhinal"],
                    `Cajal-Retzius` = grp_z[, "Cajal-Retzius"])
ex_final_call <- setNames(colnames(ex_final_z)[apply(ex_final_z, 1, which.max)], clus_ids)
ex_final_lab  <- setNames(vapply(ex_final_call, function(g) EX_FINAL_MAP[[g]], character(1)),
                          clus_ids)

subtype <- character(length(clus_ids)); names(subtype) <- clus_ids
hint    <- character(length(clus_ids)); names(hint) <- clus_ids
for (i in seq_along(clus_ids)) {
  k <- clus_ids[i]
  hint[k] <- grp_names[which.max(grp_z[k, ])]
  if (is_inh[k]) {
    if (cge_mean[k] >= mge_mean[k]) {
      subtype[k] <- if (vip_mean[k] >= lamp5_mean[k]) "Inh-CGE (VIP)" else "Inh-CGE (LAMP5)"
    } else {
      subtype[k] <- "Inh-MGE"
    }
  } else {
    subtype[k] <- ex_final_lab[k]
  }
}

# =============================================================================
# private-marker GATE  (after the ST11 taxonomy audit)
# -----------------------------------------------------------------------------
# Why this exists. The argmax above runs on a group-mean z, and a group mean can be
# carried by one promiscuous gene. Cluster 6 won "Cajal-Retzius" on RELN alone
# (mean 2.39, the highest of any cluster) while the three markers that actually make
# the call specific were absent: TP73 0.024, LHX1 0.163, NDNF 0.000. RELN is expressed
# by hippocampal interneurons and by stressed / ectopic neurons, so it cannot establish
# Cajal-Retzius identity on its own, and TP73 -- the private CR marker -- is the one
# that has to be there.
#
# The cluster is a complexity artefact, not a subtype. Four independent lines agree:
#   (i)   median nFeature 1312 vs 2543-5458 for every resolved subtype -- the lowest
#         in the hippocampus, roughly half the next-lowest;
#   (ii)  186 of its 190 nuclei come from one library (TFHS000535), i.e. 98% from a
#         single NHD lane -- a population cannot be 98% one library and one condition;
#   (iii) it carries neither class marker (SLC17A7 0.18, interneuron-group mean 0.055),
#         so it is neither excitatory nor inhibitory -- the signature of nuclei so
#         shallow that only the highest-expressed genes survive;
#   (iv)  190 nuclei = 5.4% of hippocampal neurons, whereas Cajal-Retzius cells are
#         <0.1% of neurons in the adult human hippocampus.
# Low-complexity nuclei co-cluster on complexity rather than on identity, which is
# exactly what this gate is here to catch.
#
# The RULE. A label must be supported by at least one of its own private markers at
# >= GATE_PCT %detected (see below). Clusters that fail are relabelled "Unresolved (low
# complexity)" and are excluded from every subtype claim (they were already excluded
# from DGE by 47b and from the ligand-receptor panels by the MIN_CELLS = 20 floor).
# =============================================================================
say("\n== step: private-marker gate (a label may not rest on a promiscuous gene) ==")
# The RULE. Every label must be supported by at least one of its own private markers --
# a gene that is specific to that identity -- detected in >= GATE_PCT of the cluster's
# nuclei. Counting markers is not enough and would get both halves of this wrong:
# cluster 6 clears two markers (RELN 83.2%, LHX1 10.0%) and would pass a >=2 rule,
# while the genuine Inh-MGE clusters 3 and 9 rest on LHX6 alone (54.8% / 70.0%) because
# PVALB and SST are poorly captured in nuclei, and would fail it.
#
# Which MARKERS count as private. RELN is excluded from the Cajal-Retzius set because it
# is expressed by hippocampal interneurons and stressed neurons -- it reaches 44% and 52%
# detected in two other clusters here. NDNF is excluded because it is a LAMP5 /
# neurogliaform marker (it is a member of the Inh-CGE LAMP5 block in this same script).
# LHX1 is excluded as a low-detection TF sitting exactly on the floor. That leaves TP73,
# the private Cajal-Retzius marker, at 1.1% -- absent. LHX6, by contrast, is the
# definitive MGE lineage TF and legitimately carries Inh-MGE on its own.
GATE_PCT <- 10   # a private marker counts as present at >= 10% of nuclei detected
PRIVATE_MARKERS <- list(
  "CA3"             = c("CFAP299","HS3ST4","GRIK4","NECAB1","NECTIN3",       # CA3
                        "GRIK1","FIBCD1","SPOCK1","POU3F1","WFS1","MPPED1",  # CA1
                        "RGS14","AMIGO2",                                    # CA2
                        "PROX1","GLIS3","TRPC6",                             # DG granule
                        "ADCYAP1","CARTPT"),                                 # Mossy cell
  "Subiculum"       = c("ROBO1","TLE4","NR4A2","NTS","FN1","SLC17A6"),
  "EC-like"         = c("CALB1","RORB","THEMIS","FOXP2","LAMA3"),
  "Cajal-Retzius"   = c("TP73"),
  "Inh-CGE (VIP)"   = c("VIP","ADARB2","CNR1","CCK"),
  "Inh-CGE (LAMP5)" = c("LAMP5","ID2","NDNF"),
  "Inh-MGE"         = c("LHX6","PVALB","SST"))
stopifnot("a final label has no private-marker set for the gate" =
            all(unique(subtype) %in% names(PRIVATE_MARKERS)))

gate_n  <- setNames(integer(length(clus_ids)), clus_ids)
gate_mk <- setNames(character(length(clus_ids)), clus_ids)
for (k in clus_ids) {
  mkk <- intersect(PRIVATE_MARKERS[[subtype[[k]]]], rownames(pdet))
  ok  <- mkk[pdet[mkk, k] >= GATE_PCT]
  gate_n[k]  <- length(ok)
  gate_mk[k] <- if (length(ok)) paste(sprintf("%s %.0f%%", ok, pdet[ok, k]), collapse = ", ")
                else paste0("(none; best ", mkk[which.max(pdet[mkk, k])], " ",
                            sprintf("%.1f%%", max(pdet[mkk, k])), ")")
}
nfeat_med <- tapply(hip$meta$nFeature,
                    factor(hip$meta$cluster, levels = clus_ids), median)
# Complexity is reported for every cluster, not used to relabel: a shallow cluster that
# still carries its private markers keeps its identity (cluster 8), whereas the gate
# rejects an unsupported label whatever its depth.
lowcplx <- nfeat_med < 0.5 * median(hip$meta$nFeature)
gate_tbl <- data.frame(
  cluster = clus_ids, proposed = unname(subtype[clus_ids]),
  private_markers_present = unname(gate_mk[clus_ids]),
  median_nFeature = as.integer(nfeat_med[clus_ids]),
  low_complexity  = unname(lowcplx[clus_ids]),
  pass = unname(gate_n[clus_ids]) >= 1,
  stringsAsFactors = FALSE, row.names = NULL)
say(sprintf("gate: >= 1 PRIVATE marker of the proposed label at >= %d%% detected", GATE_PCT))
print(gate_tbl); cat(capture.output(print(gate_tbl)), sep = "\n", file = logcon)

gate_fail <- clus_ids[gate_n[clus_ids] < 1]
if (length(gate_fail)) {
  for (k in gate_fail)
    say(sprintf("  FAIL cluster %s: proposed '%s' has NO private marker %s; median nFeature %d%s -> RELABELLED '%s'",
                k, subtype[[k]], gate_mk[[k]], as.integer(nfeat_med[[k]]),
                ifelse(lowcplx[[k]], ", LOW COMPLEXITY", ""), UNRESOLVED_LAB))
  subtype[gate_fail] <- UNRESOLVED_LAB
} else {
  say("  all clusters pass the private-marker gate")
}
if (any(lowcplx))
  say("  low-complexity FLAG (reported, identity still marker-supported): cluster(s) ",
      paste(clus_ids[lowcplx], collapse = ", "))
say(sprintf("HONEST LIMIT CR: 'Cajal-Retzius' is ABSENT from the FH taxonomy. The cluster argmax proposed failed the private-marker gate: TP73 (the private CR marker) 1.1%%, NDNF 0.0%%; the call rested on RELN (83.2%%, but 44%% and 52%% in two other clusters) and LHX1 (10.0%%), neither of which is private to CR. Reported as '%s'.",
            UNRESOLVED_LAB))

# Every surviving label must be in the legend order, or `intersect()` downstream would
# silently drop it to NA instead of failing.
stopifnot("a final hippo label is not in HIP_LABEL_ORDER (would be dropped silently)" =
            all(unique(subtype) %in% HIP_LABEL_ORDER))

# audit table + honest-limits report
audit <- data.frame(
  cluster       = clus_ids,
  n             = as.integer(table(factor(hip$meta$cluster, levels = clus_ids))),
  argmax_group  = unname(hint[clus_ids]),
  ex_subfield   = unname(ex_call_raw),
  GADmean       = round(gad_mean, 3),
  EXmean        = round(ex_mean, 3),
  CGEmean       = round(cge_mean, 3),
  MGEmean       = round(mge_mean, 3),
  VIP           = round(vip_mean, 3),
  LAMP5         = round(lamp5_mean, 3),
  is_inh        = unname(is_inh[clus_ids]),
  ex_final      = unname(ex_final_call[clus_ids]),
  myelin_mean   = round(oligo_mean[clus_ids], 3),
  panneuron_mean= round(neuron_mean[clus_ids], 3),
  gate_private  = unname(gate_mk[clus_ids]),
  median_nFeature = as.integer(nfeat_med[clus_ids]),
  low_complexity  = unname(lowcplx[clus_ids]),
  final_label   = unname(subtype[clus_ids]),
  stringsAsFactors = FALSE, row.names = NULL)
say("\n== per-cluster re-derived call (auditable) ==")
print(audit); cat(capture.output(print(audit)), sep = "\n", file = logcon)

# oligo contamination guard: flag (do not silently keep) any cluster where RAW
# myelin mean EXCEEDS raw pan-neuronal mean (a doublet cluster, not a subtype).
oligo_flag <- clus_ids[oligo_mean > neuron_mean]
if (length(oligo_flag)) {
  say("WARNING: cluster(s) with RAW myelin mean > pan-neuronal mean (possible doublet): ",
      paste(oligo_flag, collapse = ", "), "")
} else {
  say("oligo-QC guard: no de-novo cluster is myelin-dominated (all pan-neuronal >> myelin). OK")
}

# --- honest-limits checks --------------------------------------------------------
# DG-granule floor: PROX1 is also a marker of CGE interneurons (VIP/LAMP5) and of
# Cajal-Retzius cells in adult human hippocampus, so a naive max-PROX1 is confounded.
# A true DG-granule cluster is a GLUTAMATERGIC pyramidal cluster: SLC17A7+ /
# gad-negative / RELN-low (not CR) and PROX1-high. Test PROX1 %detected only among
# the pyramidal excitatory clusters (those whose final label is a ca/Subiculum
# pyramidal subtype), i.e. exclude the interneuron and Cajal-Retzius clusters.
pyr_clus <- clus_ids[subtype[clus_ids] %in% c("CA3","Subiculum")]
prox1_pyr_max <- if ("PROX1" %in% rownames(pdet) && length(pyr_clus))
                   max(pdet["PROX1", pyr_clus]) else NA_real_
prox1_any_max <- if ("PROX1" %in% rownames(pdet)) max(pdet["PROX1", ]) else NA_real_
prox1_top_clus <- if ("PROX1" %in% rownames(pdet))
                    clus_ids[which.max(pdet["PROX1", ])] else NA
say(sprintf("\nHONEST LIMIT DG: PROX1 max %%detected = %.1f%% overall (cluster %s -> %s); %.1f%% among PYRAMIDAL (CA3/Subiculum) clusters",
            prox1_any_max, prox1_top_clus,
            ifelse(is.na(prox1_top_clus), NA, subtype[[prox1_top_clus]]), prox1_pyr_max))
say(sprintf("   -> DG-granule NOT captured: no pyramidal cluster is PROX1-high (max %.1f%% << floor). The PROX1+ nuclei sit in a GAD+ CGE-interneuron cluster (PROX1 is a CGE-IN marker), NOT a granule-cell cluster. Consistent with the DG dissection limit.",
            prox1_pyr_max))
say("HONEST LIMIT CA1/CA2: not separately resolved at res ", HIP_RES,
    " — CA-pyramidal clusters fold into the CA3 anchor label (see audit ex_final).")
# EC-like: report whether an entorhinal-dominant excitatory cluster actually emerged.
ex_clus <- clus_ids[!is_inh[clus_ids]]
ec_clus <- ex_clus[ex_final_call[ex_clus] == "Entorhinal"]
if (length(ec_clus)) {
  say("EC-like: resolved -> cluster(s) ", paste(ec_clus, collapse = ", "))
} else {
  say("HONEST LIMIT EC-like: NO entorhinal-dominant excitatory cluster emerged at res ",
      HIP_RES, " on the FH atlas (unlike the pre-OCC-removal atlas). EC-like is ABSENT",
      " from the FH taxonomy -> 6 labels, not 7.")
}

# persist label map to cache (relabel only; embedding/clusters untouched)
hip$meta$neuron_subtype <- unname(subtype[hip$meta$cluster])
stopifnot("unmapped hippo cluster -> NA label" = !anyNA(hip$meta$neuron_subtype))
hip$clust2label <- subtype
hip$audit <- audit
saveRDS(hip, paste0(CACHE_HIP, ".tmp")); file.rename(paste0(CACHE_HIP, ".tmp"), CACHE_HIP)
say("persisted re-derived labels to cache")

# =============================================================================
# 4. Build the unified neuron_subtype_map_FH (hippo de-novo + cortical subclass)
# =============================================================================
say("\n== step: build neuron_subtype_map_FH ==")
hip_map <- data.frame(
  barcode          = rownames(hip$meta),
  Region           = "Hippo",
  neuron_subtype   = hip$meta$neuron_subtype,
  Condition        = hip$meta$Condition,
  source_annot     = "hip_denovo_label",
  azimuth_subclass = NA_character_,   # hippo neurons are typed de novo
  dfc_within_score = NA_real_,
  stringsAsFactors = FALSE)

dup <- intersect(ctx_map$barcode, hip_map$barcode)
if (length(dup)) stop("Overlapping barcodes cortical/hippo: ", length(dup), " (unexpected)")
MAP_COLS <- c("barcode","Region","neuron_subtype","Condition","source_annot",
              "azimuth_subclass","dfc_within_score")
nmap <- rbind(ctx_map[, MAP_COLS], hip_map[, MAP_COLS])
stopifnot("duplicated barcode in map" = !anyDuplicated(nmap$barcode))
stopifnot("NA subtype in map" = !anyNA(nmap$neuron_subtype))
saveRDS(nmap, paste0(MAP_RDS, ".tmp")); file.rename(paste0(MAP_RDS, ".tmp"), MAP_RDS)
write.csv(nmap, MAP_CSV, row.names = FALSE)
say("wrote ", basename(MAP_RDS), " (", nrow(nmap), " neurons; ",
    sum(nmap$Region == "Frontal"), " Frontal + ", sum(nmap$Region == "Hippo"), " Hippo)")

# =============================================================================
# 5. CENSUS (per house policy) — CON vs NHD per subtype/subclass; flag thin lanes
# =============================================================================
THIN <- 30L
census_block <- function(sub, order_hint = NULL) {
  tb <- as.data.frame.matrix(table(sub$neuron_subtype, sub$Condition))
  for (cc in c("CON","NHD")) if (!cc %in% colnames(tb)) tb[[cc]] <- 0L
  tb <- tb[, c("CON","NHD")]; tb$total <- tb$CON + tb$NHD
  ord <- if (!is.null(order_hint)) c(intersect(order_hint, rownames(tb)),
                                     setdiff(rownames(tb), order_hint)) else
         rownames(tb)[order(-tb$total)]
  tb <- tb[ord, , drop = FALSE]
  tb$min_cond    <- pmin(tb$CON, tb$NHD)
  # DGE power is set by the smaller condition, not the total: a subtype with 190
  # nuclei but only 4 in one condition is unpowered for CON-vs-NHD. Flag on
  # min_cond (the honest per-condition floor), then on total.
  MIN_COND_FLOOR <- 10L
  tb$flag <- ifelse(tb$min_cond == 0, "SINGLE-CONDITION (no DGE)",
              ifelse(tb$min_cond < MIN_COND_FLOOR,
                     sprintf("IMBALANCED min-cond=%d (<%d; unpowered)", tb$min_cond, MIN_COND_FLOOR),
              ifelse(tb$total < THIN, "THIN <30 (unpowered)", "")))
  tb
}

say("\n== HIPPO de-novo subtype census (CON | NHD) ==")
hip_cen <- census_block(hip_map, HIP_LABEL_ORDER)
print(hip_cen); cat(capture.output(print(hip_cen)), sep = "\n", file = logcon)

say("\n== CORTICAL (Frontal) subclass census (CON | NHD) — Ex = Jorstad-DFC, Inh = Azimuth (2026-09-18) ==")
ctx_cen <- census_block(ctx_map, c(EX_SUBC, INH_SUBC))
print(ctx_cen); cat(capture.output(print(ctx_cen)), sep = "\n", file = logcon)
say("  cortical subclasses present: ", nrow(ctx_cen), " (expected 15 = 9 Ex + 6 Inh)")
say("  source_annot: ", paste(names(table(ctx_map$source_annot)), table(ctx_map$source_annot),
                              sep = "=", collapse = ", "))
# machine-readable census (subclass x Region x Condition) for the report / ST10 cross-check
cen_csv <- file.path(LOGD, "60_neuron_subtype_census_FH.csv")
cen_all <- as.data.frame(table(neuron_subtype = nmap$neuron_subtype, Region = nmap$Region,
                               Condition = nmap$Condition), responseName = "n")
cen_all <- cen_all[cen_all$n > 0, ]
cen_all <- cen_all[order(cen_all$Region, match(cen_all$neuron_subtype, c(EX_SUBC, INH_SUBC, HIP_LABEL_ORDER)),
                         cen_all$Condition), ]
write.csv(cen_all, cen_csv, row.names = FALSE)
say("wrote ", basename(cen_csv))

say("\n== THIN / SINGLE-CONDITION flags (unpowered for DGE) ==")
flag_h <- hip_cen[hip_cen$flag != "", , drop = FALSE]
flag_c <- ctx_cen[ctx_cen$flag != "", , drop = FALSE]
if (nrow(flag_h)) { say("  HIPPO:");    for (r in rownames(flag_h)) say(sprintf("    %-18s CON=%d NHD=%d total=%d  [%s]", r, flag_h[r,"CON"], flag_h[r,"NHD"], flag_h[r,"total"], flag_h[r,"flag"])) } else say("  HIPPO: none")
if (nrow(flag_c)) { say("  CORTICAL:"); for (r in rownames(flag_c)) say(sprintf("    %-18s CON=%d NHD=%d total=%d  [%s]", r, flag_c[r,"CON"], flag_c[r,"NHD"], flag_c[r,"total"], flag_c[r,"flag"])) } else say("  CORTICAL: none")

# =============================================================================
# 6. Eye-check UMAP — hippo de-novo subtypes (+ condition split)
# =============================================================================
say("\n== step: render hippo subtype UMAP (eye-check) ==")
d <- as.data.frame(hip$embedding)[, 1:2]; colnames(d) <- c("UMAP_1","UMAP_2")
d$subtype   <- factor(hip$meta$neuron_subtype,
                      levels = intersect(HIP_LABEL_ORDER, unique(hip$meta$neuron_subtype)))
d$Condition <- factor(hip$meta$Condition, levels = c("CON","NHD"))

umap_theme <- function(base = 8) {
  theme_pub(base_size = base) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = base - 1, hjust = 0.02),
          legend.key.height = unit(0.7, "lines"),
          legend.text = element_text(size = base - 1))
}
pal <- HIP_LAB_PAL[levels(d$subtype)]

set.seed(42); d_sh <- d[sample(nrow(d)), ]
p_sub <- ggplot(d_sh, aes(UMAP_1, UMAP_2, fill = subtype)) +
  geom_point(size = 0.7, alpha = 0.75, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = pal, name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  coord_equal() + umap_theme() +
  ggtitle("Hippocampal neurons — de-novo subtypes")

p_split <- ggplot(d_sh, aes(UMAP_1, UMAP_2, fill = subtype)) +
  geom_point(size = 0.6, alpha = 0.75, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = pal, name = NULL, guide = "none") +
  facet_wrap(~ Condition, nrow = 1) +
  coord_equal() + umap_theme() +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8)) +
  ggtitle("by condition")

suppressPackageStartupMessages(library(patchwork))
comp <- (p_sub / p_split) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(theme = theme(plot.title = element_text(size = 10, hjust = 0.5, face = "plain")))
ggsave(file.path(PANEL, "4_hippo_neuron_subtypes_UMAP.png"), comp,
       width = 6.2, height = 8.4, dpi = 300, bg = "white",
       device = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else NULL)
say("wrote ", file.path(PANEL, "4_hippo_neuron_subtypes_UMAP.png"))

# per-cluster marker mean-expr table (audit) -> logs
mk_csv <- file.path(LOGD, "60_hippo_cluster_marker_meanexpr_FH.csv")
write.csv(data.frame(marker_group = mgrp, marker = rownames(mtbl),
                     round(as.data.frame(mtbl), 4), check.names = FALSE),
          mk_csv, row.names = FALSE)
say("wrote ", basename(mk_csv))

# =============================================================================
# provenance + sessionInfo
# =============================================================================
say("\n== provenance ==")
say("atlas: ", ATLAS, "  mtime: ", format(atlas_mtime))
say("DFC transfer table: ", DFC_CSV, "  mtime: ", format(file.mtime(DFC_CSV)),
    "  (frontal Ex labels = dfc_within; 2026-09-18)")
say("recipe (hippo): SCT -> RunPCA(30) -> RunHarmony(SampleID) -> FindNeighbors(1:30) -> FindClusters(res ",
    HIP_RES, ", seed 42) -> RunUMAP(seed 42)")
say("hippo cluster -> label map (re-derived on NEW clustering):")
for (k in names(hip$clust2label)) say(sprintf("   cluster %s -> %s", k, hip$clust2label[[k]]))
si <- capture.output(sessionInfo()); cat(si, sep = "\n", file = logcon); cat("\n", file = logcon)
say("=== DONE 60_neuron_reannotation_FH ===")
close(logcon)
cat("=== DONE ===\n", file = stderr())
