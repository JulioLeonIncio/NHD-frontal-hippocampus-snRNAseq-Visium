#!/usr/bin/env Rscript
# =============================================================================
# 3M_prep_opc_and_scales_FH.R — Shared prep for the Figure-4 additions (Frontal + Hippocampus rebuild; no OCC anywhere).
# -----------------------------------------------------------------------------
# Loads the FH atlas once and emits three small, reusable artifacts:
#
#   1. data/_cache_opc_FH.rds
#        OPC nuclei (RNA counts + LogNormalize data), so the module-score panels
#        (3Fp lineage block / dot-matrix) can show OPC beside Oligo without any
#        further atlas load.  Mirrors how _cache_oligo_umap.rds is consumed.
#
#   2. tables/myelin_scales_per_nucleus_FH.csv
#        Per nucleus: total UMI, myelin UMI, cell type, region, condition, lane.
#        This is the source for the new "three denominators" panel, which
#        resolves the apparent contradiction between Figure 4 (per-cell myelin
#        transcripts read up in NHD) and Figure 6/Visium (tissue myelin COLLAPSES).
#
#   3. tables/myelin_scales_per_spot_visium_FH.csv
#        The same quantity measured on Visium spots (myelin UMI / total UMI per
#        spot), i.e. myelin's share of the tissue transcriptome.  Optional: the
#        Visium object lives outside the rebuild tree, so its absence is a
#        warning, not an error.
#
# Why this exists.  Three different denominators were being compared
# as if they were one measurement:
#   tissue  : myelin per unit tissue            -> DOWN   (Visium; oligo loss)
#   nucleus : myelin UMIs per oligo nucleus     -> ~FLAT  (20 -> 22 median)
#   share   : myelin's share of that nucleus    -> up ~1.9x (the Fig-4 log2FC)
# The bridge is that NHD frontal oligodendrocytes carry roughly half the
# transcriptome of a control oligodendrocyte (median 1520 -> 799 UMI), and are
# ~4x rarer.  All three numbers are correct and the framing follows from them.
# Verified not to be a technical artifact:
#   * depth-matched binomial downsampling makes the share effect LARGER, not
#     smaller (Cliff's d +0.78 at T=750, retaining only well-covered nuclei);
#   * the complexity drop is cell-type-specific within LANE (oligo/other nFeature
#     ratio 0.22-0.26 in both NHD frontal lanes vs 0.43-0.60 in CON lanes), so it
#     cannot be capture chemistry or sequencing depth;
#   * ambient/debris is ruled out (myelin is 0.009% of the astrocyte and 0.048%
#     of the myeloid transcriptome at matched depth, vs 2.5% in oligodendrocytes).
#
# Atlas is the single source; every cache carries an mtime staleness guard.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))

DDIR <- file.path(PROJ, "data")
TDIR <- file.path(PROJ, "tables")
LOGD <- file.path(PROJ, "logs")
for (d in c(DDIR, TDIR, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ATLAS   <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
OPC_RDS <- file.path(DDIR, "_cache_opc_FH.rds")
NUC_CSV <- file.path(TDIR, "myelin_scales_per_nucleus_FH.csv")
SPOT_CSV<- file.path(TDIR, "myelin_scales_per_spot_visium_FH.csv")
stopifnot("MISSING atlas NHD_FH_harmony.rds" = file.exists(ATLAS))

# The myelin gene set is the same curated list used by the Fig-4 lineage panels and
# by the Visium myelin heatmap, so the three scales measure one quantity.
# Use the AUDITED structural-myelin set, not sigs$Myelination_genes.
# The curated set carried TF (delta -0.44/-0.55) and ANLN (-0.37/-0.46) — both strongly
# down — and omitted OPALIN (+0.46/+0.29), the one gene genuinely induced. Averaging an
# up arm against a down arm drove the module to ~0. Every quantity in this file feeds
# panels that assert "the cell still orders myelin", so it must measure the membrane
# proteins only.
sigs <- readRDS(file.path(DDIR, "curated_signatures.rds"))
source(file.path(PROJ, "scripts", "_oligo_programmes_FH.R"))
MYE  <- oligo_programme_sets(sigs)[["Structural myelin"]]
cat("structural myelin gene set (", length(MYE), "):", paste(MYE, collapse = ", "), "\n")
cat("  (excluded by the direction audit and reported separately:",
    paste(oligo_myelin_support_genes, collapse = ", "), ")\n")

# ---------------------------------------------------------------------------
# 1 + 2.  snRNA atlas
# ---------------------------------------------------------------------------
# Version the gene set and
# force a rebuild when it changes. (Dropbox also rewrites mtimes on CloudStorage, so
# mtime alone is not trustworthy in this project.)
MYE_KEY  <- paste0("structural_myelin:", paste(sort(MYE), collapse = ","))
KEY_FILE <- file.path(DDIR, "_myelin_scales_geneset.key")
key_stale <- !file.exists(KEY_FILE) || !identical(readLines(KEY_FILE, warn = FALSE)[1], MYE_KEY)
if (key_stale) cat("gene set changed since the last run — forcing a rebuild\n")

need_opc <- !file.exists(OPC_RDS) || file.mtime(OPC_RDS) < file.mtime(ATLAS)
need_nuc <- !file.exists(NUC_CSV) || file.mtime(NUC_CSV) < file.mtime(ATLAS) || key_stale

if (need_opc || need_nuc) {
  cat("Loading atlas (this is the only atlas load in the Fig-4 additions)...\n")
  obj <- readRDS(ATLAS)
  DefaultAssay(obj) <- "RNA"
  if (length(SeuratObject::Layers(obj, assay = "RNA")) > 1)
    obj <- SeuratObject::JoinLayers(obj)
  stopifnot("atlas is not OCC-free — wrong object" = !any(grepl("OCC", obj$Region)))
  cat(sprintf("atlas: %d nuclei, regions %s\n", ncol(obj),
              paste(sort(unique(as.character(obj$Region))), collapse = "/")))

  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
  mye <- intersect(MYE, rownames(cnt))
  cat(sprintf("myelin genes present in atlas: %d/%d\n", length(mye), length(MYE)))

  if (need_nuc) {
    md <- obj@meta.data
    out <- data.frame(
      cell       = colnames(cnt),
      cell_type  = as.character(md$new_annotation),
      Region     = as.character(md$Region),
      Condition  = as.character(md$Condition),
      SampleID   = as.character(md$SampleID),
      total_umi  = as.numeric(Matrix::colSums(cnt)),
      myelin_umi = as.numeric(Matrix::colSums(cnt[mye, , drop = FALSE])),
      nFeature   = as.numeric(md$nFeature_RNA),
      stringsAsFactors = FALSE)
    tmp <- paste0(NUC_CSV, ".tmp"); write.csv(out, tmp, row.names = FALSE)
    file.rename(tmp, NUC_CSV)
    writeLines(MYE_KEY, KEY_FILE)
    cat(sprintf("Wrote %s (%d nuclei)\n", basename(NUC_CSV), nrow(out)))
  }

  if (need_opc) {
    o <- subset(obj, subset = new_annotation == "OPC")
    o <- DietSeurat(o, assays = "RNA", layers = "counts")
    DefaultAssay(o) <- "RNA"
    o <- NormalizeData(o, verbose = FALSE)
    cat(sprintf("OPC cache: %d nuclei (Frontal %d, Hippo %d; CON %d, NHD %d)\n",
                ncol(o), sum(o$Region == "Frontal"), sum(o$Region == "Hippo"),
                sum(o$Condition == "CON"), sum(o$Condition == "NHD")))
    tmp <- paste0(OPC_RDS, ".tmp"); saveRDS(o, tmp); file.rename(tmp, OPC_RDS)
    cat("Wrote", basename(OPC_RDS), "\n")
    rm(o)
  }
  rm(obj, cnt); gc(verbose = FALSE)
} else {
  cat("snRNA artifacts are current (newer than the atlas) — skipped.\n")
}

# ---------------------------------------------------------------------------
# 3.  Visium spots — the same quantity at tissue scale
# ---------------------------------------------------------------------------
VIS <- file.path(dirname(PROJ), "Visium", "integrated_harmony",
                 "NHD_frontal_integrated_harmony.rds")
if (!file.exists(VIS)) {
  warning("Visium object not found at ", VIS,
          " — the tissue-scale row of the panel will be omitted.")
} else if (file.exists(SPOT_CSV) && file.mtime(SPOT_CSV) >= file.mtime(VIS) && !key_stale) {
  cat("Visium spot table is current — skipped.\n")
} else {
  cat("Loading Visium object...\n")
  v <- readRDS(VIS)
  DefaultAssay(v) <- "Spatial"
  vc  <- GetAssayData(v, assay = "Spatial", layer = "counts")
  myev <- intersect(MYE, rownames(vc))
  cat(sprintf("myelin genes present in Visium: %d/%d\n", length(myev), length(MYE)))
  vm <- v@meta.data
  sid <- if ("sample_id" %in% colnames(vm)) as.character(vm$sample_id) else as.character(vm$orig.ident)
  spot <- data.frame(
    spot       = colnames(vc),
    sample_id  = sid,
    Condition  = ifelse(grepl("^NHD", sid), "NHD", "CON"),
    total_umi  = as.numeric(Matrix::colSums(vc)),
    myelin_umi = as.numeric(Matrix::colSums(vc[myev, , drop = FALSE])),
    stringsAsFactors = FALSE)
  tmp <- paste0(SPOT_CSV, ".tmp"); write.csv(spot, tmp, row.names = FALSE)
  file.rename(tmp, SPOT_CSV)
  cat(sprintf("Wrote %s (%d spots, %d sections)\n", basename(SPOT_CSV),
              nrow(spot), length(unique(spot$sample_id))))
  rm(v, vc); gc(verbose = FALSE)
}

writeLines(capture.output(sessionInfo()),
           file.path(LOGD, "3M_prep_opc_and_scales_FH_sessionInfo.txt"))
cat("=== DONE ===\n")
