#!/usr/bin/env Rscript
# =============================================================================
# 102_export_deposition_bundle_FH.R — Processed-data bundle for repository deposition (GEO / DDBJ-GEA style), so the Data Availability statement can be filled the day the
# repository is chosen.
# -----------------------------------------------------------------------------
# What is EXPORTED (manuscript/Final_fig_and_tables/Data_deposition/):
#   snRNA/   raw UMI counts of the 55,354 post-QC nuclei as 10x-style MTX
#            (matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz), the per-nucleus
#            metadata, and the
#            per-library sample sheet.
#   Visium/  per section: Space Ranger filtered_feature_bc_matrix.h5 + spatial/ folder
#            copied verbatim from the outs named in Visium/VISIUM_OUTS_MANIFEST.json (
#            <sec>_masked/outs, three-tier tissue call; tissue_call.csv + TISSUE_CALL_provenance.json
#            ship per section), plus spot metadata (condition, section, spatial domain,
#            depth, cell2location abundances where present) from the analysed object.
#   GEO_samples_metadata.csv   one row per library/section in GEO "samples" layout with
#            characteristics filled from Table 1 / ST12; raw-file columns left blank —
#            FASTQs are at RIKEN and are listed by the submitter.
#   README.md + md5sums.txt
#
# The counts come from the same objects every figure reads (atlas/NHD_FH_harmony.rds for
# metadata + embeddings; the RNA assay counts layer for UMIs), so the deposit is the
# analysed data, not a re-processed copy. Nothing here is filtered or normalised.
#
# Run:  Rscript scripts/102_export_deposition_bundle_FH.R   (needs ~8 GB RAM, ~10 min)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table); library(tools) })

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
say <- function(...) cat(sprintf(...), "\n")
ROOT <- dirname(PROJ)
source(file.path(PROJ, "scripts", "_visium_outs_FH.R"))   # VIS_SECTIONS + vis_outs(sec) from VISIUM_OUTS_MANIFEST.json
OUT  <- file.path(PROJ, "manuscript", "Final_fig_and_tables", "Data_deposition")
SN   <- file.path(OUT, "snRNA"); VI <- file.path(OUT, "Visium")
for (d in c(OUT, SN, VI)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- snRNA ------------------------------------------------------------------
ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
stopifnot(file.exists(ATLAS))
META_CSV <- file.path(SN, "nucleus_metadata.csv.gz")
sn_current <- all(file.exists(file.path(SN, c("matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz")), META_CSV)) &&
  min(file.mtime(file.path(SN, c("matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz", "nucleus_metadata.csv.gz")))) > file.mtime(ATLAS)
if (sn_current) {
  say("snRNA export already newer than the atlas — reusing it (delete snRNA/ to force a rebuild)")
  meta <- fread(META_CSV)
  if ("percent_mito" %in% names(meta)) { meta[, percent_mito := NULL]; fwrite(meta, META_CSV) }
  if (!"matter" %in% names(meta)) {
    lm <- fread(file.path(ROOT, "manuscript_7fig", "tables", "sample_metadata_map.csv"))
    lm[, matter := fifelse(grepl("_g$", Sample), "grey", fifelse(grepl("_w$", Sample), "white", "whole"))]
    meta[, matter := lm$matter[match(SampleID, lm$SampleID)]]
    setcolorder(meta, c("barcode","SampleID","Condition","Region","matter")); fwrite(meta, META_CSV)
  }
} else {
say("loading atlas %s (%.1f GB) ...", basename(ATLAS), file.size(ATLAS) / 1e9)
a <- readRDS(ATLAS)
stopifnot(ncol(a) == 55354, "RNA" %in% names(a@assays))
DefaultAssay(a) <- "RNA"
# Seurat v5: the RNA assay is split into one counts layer per library — join before reading
a[["RNA"]] <- JoinLayers(a[["RNA"]])
cnt <- LayerData(a, assay = "RNA", layer = "counts")
stopifnot(inherits(cnt, "dgCMatrix"), ncol(cnt) == 55354, max(cnt@x %% 1) == 0)   # raw integer UMIs
say("counts: %d genes x %d nuclei, %s UMIs", nrow(cnt), ncol(cnt), format(sum(cnt), big.mark = ","))

# neuron_subtype is attached for all neurons after the branch (see below) from the
# single-source map; here only the atlas columns are exported.
md <- a@meta.data
um <- Embeddings(a, "umap")
keep_cols <- intersect(c("SampleID","Condition","Region","matter","hippo_area","new_annotation",
                         "predicted.subclass","predicted.subclass.score","seurat_clusters",
                         "nCount_RNA","nFeature_RNA","percent.mt"), names(md))   # percent_mito is a duplicate of percent.mt
meta <- data.table(barcode = colnames(a), md[, keep_cols, drop = FALSE],
                   UMAP_1 = um[, 1], UMAP_2 = um[, 2])
stopifnot(!anyDuplicated(meta$barcode), nrow(meta) == ncol(cnt))
say("libraries: %s", paste(names(table(meta$SampleID)), table(meta$SampleID), collapse = "; "))

writeMM(cnt, file.path(SN, "matrix.mtx"))
fwrite(data.table(gene = rownames(cnt), gene2 = rownames(cnt), type = "Gene Expression"),
       file.path(SN, "features.tsv"), sep = "\t", col.names = FALSE)
fwrite(data.table(colnames(cnt)), file.path(SN, "barcodes.tsv"), sep = "\t", col.names = FALSE)
for (f in c("matrix.mtx","features.tsv","barcodes.tsv")) { p <- file.path(SN, f); if (file.exists(paste0(p, ".gz"))) unlink(paste0(p, ".gz"))
  stopifnot(system2("gzip", shQuote(p)) == 0, file.exists(paste0(p, ".gz"))) }
fwrite(meta, META_CSV)
rm(a, cnt); invisible(gc())
}

# Applied in both branches (fresh export or reuse) so a reused metadata file is upgraded in
# place. `predicted.subclass` (Azimuth, motor-cortex reference) is renamed `azimuth_subclass`
# and kept as a separate column; `neuron_subtype` carries the DFC/Azimuth/de-novo label used
# by every figure, `neuron_subtype_source` says which reference labelled each nucleus.
source(file.path(PROJ, "scripts", "_neuron_subclass_FH.R"))
if ("predicted.subclass" %in% names(meta)) setnames(meta, "predicted.subclass", "azimuth_subclass")
if ("predicted.subclass.score" %in% names(meta)) setnames(meta, "predicted.subclass.score", "azimuth_subclass_score")
meta[, neuron_subtype := unname(SUBCLASS_MAP[barcode])]
meta[, neuron_subtype_source := unname(SUBCLASS_SOURCE_MAP[barcode])]
.nmap <- readRDS(NEURON_MAP_PATH)
meta[, dfc_within_score := .nmap$dfc_within_score[match(barcode, .nmap$barcode)]]
stopifnot("dfc_within_score must be present for every DFC-labelled nucleus" =
            !anyNA(meta$dfc_within_score[meta$neuron_subtype_source %in% "ctx_DFC_Jorstad2023_within"]))
# column placement: the three Figure-5 label columns sit right after the Azimuth pair.
.lab_cols <- c("neuron_subtype","neuron_subtype_source","dfc_within_score")
.anchor   <- intersect(c("azimuth_subclass_score","azimuth_subclass"), names(meta))[1]
if (!is.na(.anchor)) {
  .rest <- setdiff(names(meta), .lab_cols); .i <- match(.anchor, .rest)
  setcolorder(meta, append(.rest, .lab_cols, after = .i))
}
stopifnot("every Neuron_Ex/Neuron_Inh nucleus must carry a neuron_subtype" =
            !anyNA(meta$neuron_subtype[meta$new_annotation %in% c("Neuron_Ex","Neuron_Inh")]),
          "neuron_subtype on a non-neuron nucleus" =
            all(is.na(meta$neuron_subtype[!meta$new_annotation %in% c("Neuron_Ex","Neuron_Inh")])))
say("neuron_subtype attached for %d nuclei (%s)", sum(!is.na(meta$neuron_subtype)),
    paste(names(table(meta$neuron_subtype_source)), table(meta$neuron_subtype_source), sep = "=", collapse = "; "))
fwrite(meta, META_CSV)

# matter (grey / white / whole hippocampus) is not a column of the atlas; every script derives
# it from the lane map, so the deposit does the same and from the same file.
lm <- fread(file.path(ROOT, "manuscript_7fig", "tables", "sample_metadata_map.csv"))
lm[, matter := fifelse(grepl("_g$", Sample), "grey", fifelse(grepl("_w$", Sample), "white", "whole"))]
meta[, matter := lm$matter[match(SampleID, lm$SampleID)]]
stopifnot(!anyNA(meta$matter))
setcolorder(meta, c("barcode","SampleID","Condition","Region","matter"))
# per-library sample sheet from the shipped QC table (ST12) — never re-derived
st12 <- fread(file.path(PROJ, "supplementary_tables", list.files(file.path(PROJ, "supplementary_tables"), pattern = "^ST12_")[1]))
fwrite(st12, file.path(SN, "library_sample_sheet.csv"))

# ---- Visium -------------------------------------------------------------------
VIS <- file.path(ROOT, "Visium", "integrated_harmony", "NHD_frontal_integrated_harmony.rds")
stopifnot(file.exists(VIS))
stopifnot("Visium object is older than the tissue-call manifest (117b) — re-run integrated_harmony.R" =
          file.mtime(VIS) > file.mtime(file.path(ROOT, "Visium", "VISIUM_OUTS_MANIFEST.json")))
v <- readRDS(VIS)
vm <- v@meta.data
sec_col <- intersect(c("sample_id","orig.ident","Sample","section"), names(vm))[1]
say("Visium object: %d spots, %d sections (%s)", ncol(v), length(unique(vm[[sec_col]])), sec_col)
c2l <- file.path(ROOT, "Visium", "cell2location", "c2l_MAIN", "abundance_q05.csv")
spot <- data.table(spot_id = colnames(v), vm)
CT8 <- c("Astro","Endo","Micro-PVM","Neuron_Ex","Neuron_Inh","OPC","Oligo","Pericytes")
if (file.exists(c2l)) {
  ab <- fread(c2l)[, c("spot_id", CT8), with = FALSE]
  setnames(ab, CT8, paste0("c2l_q05_", CT8))
  spot <- merge(spot, ab, by = "spot_id", all.x = TRUE, sort = FALSE)
}
# spatial-domain identity per integrated cluster (the labels used in Fig 6a/b)
idf <- file.path(ROOT, "Visium", "integrated_harmony", "integrated_cluster_identity.csv")
if (file.exists(idf) && "seurat_clusters" %in% names(spot)) {
  id <- fread(idf)[, .(cluster = sub("^g", "", cluster), identity)]
  spot[, spatial_domain := id$identity[match(as.character(seurat_clusters), id$cluster)]]
  stopifnot("spatial_domain mapping matched nothing (cluster prefix changed?)" = !all(is.na(spot$spatial_domain)))
}
fwrite(spot, file.path(VI, "spot_metadata.csv.gz"))
sections <- VIS_SECTIONS
for (s in sections) {
  src <- vis_outs(s)   # outs per section from the manifest; vis_outs() stops if the folder is missing
  if (!dir.exists(src)) { warning("no Space Ranger outs for ", s); next }
  dst <- file.path(VI, s); dir.create(dst, showWarnings = FALSE)
  for (f in c("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5", "metrics_summary.csv"))
    if (file.exists(file.path(src, f))) file.copy(file.path(src, f), file.path(dst, f), overwrite = TRUE)
  if (dir.exists(file.path(src, "spatial"))) file.copy(file.path(src, "spatial"), dst, recursive = TRUE, overwrite = TRUE)
  # Every section: the per-barcode three-tier tissue call (tissue_call_source = spaceranger | he_mask |
  # expression | excluded) and its provenance (117b) — so a reader can rebuild the filtered matrix from the raw one
  for (f in c("tissue_call.csv", "TISSUE_CALL_provenance.json")) {
    stopifnot("tissue call file missing in outs" = file.exists(file.path(src, f)))
    file.copy(file.path(src, f), file.path(dst, f), overwrite = TRUE) }
  # NHD_Frontal2 only: the delivered-vs-reflected barcode pairing and the registration provenance (108d), kept as before
  for (f in c("in_tissue_registered.csv", "REGISTRATION_provenance.json"))
    if (file.exists(file.path(src, f))) file.copy(file.path(src, f), file.path(dst, f), overwrite = TRUE)
  say("Visium %s: %s", s, paste(list.files(dst, recursive = TRUE), collapse = ", "))
}
rm(v); invisible(gc())

# ---- GEO-style samples sheet ------------------------------------------------------
t1 <- fread(file.path(PROJ, "manuscript", "Final_fig_and_tables", "Tables", "Table1_donor_sample_FH.csv"))
g <- function(ch, col) t1[Characteristic == ch][[col]][1]
# libraries from the atlas metadata itself (SampleID = TFHS id), so the sheet cannot disagree
# with the deposited barcodes; ST12 (Donor x Region x Matter) ships beside it as the QC sheet.
mcol <- intersect(c("matter","Matter"), names(meta))[1]
lib <- meta[, .(n_nuclei = .N), by = c("SampleID", "Condition", "Region", if (!is.na(mcol)) mcol)][order(SampleID)]
tissue_lab <- function(L) {   # curator-readable, not the analysis tokens
  m <- if (!is.na(mcol)) L[[mcol]] else "whole"
  if (L$Region == "Frontal") paste0("prefrontal cortex, ", m, " matter") else "hippocampus"
}
rows <- list()
for (i in seq_len(nrow(lib))) {
  L <- lib[i]; cond <- L$Condition
  rows[[length(rows) + 1]] <- data.table(
    `Sample name` = L$SampleID, title = sprintf("%s %s snRNA-seq (%s)", cond, tissue_lab(L), L$SampleID),   # SampleID keeps titles unique (two control hippocampal captures)
    `source name` = "human brain, post-mortem",
    organism = "Homo sapiens", `characteristics: diagnosis` = if (cond == "NHD") "Nasu-Hakola disease" else "control",
    `characteristics: genotype` = if (cond == "NHD") g("TYROBP / DAP12 genotype", "NHD donor") else "TYROBP wild-type",
    `characteristics: sex` = "female", `characteristics: age` = if (cond == "NHD") g("Age at death (years)", "NHD donor") else g("Age at death (years)", "Control donor"),
    `characteristics: post-mortem interval` = if (cond == "NHD") g("Post-mortem interval", "NHD donor") else g("Post-mortem interval", "Control donor"),
    `characteristics: tissue` = tissue_lab(L), `post-QC nuclei` = L$n_nuclei,
    molecule = "nuclear RNA", `library strategy` = "snRNA-seq (10x Chromium 5' PE)",
    `processed data file` = "snRNA/matrix.mtx.gz; snRNA/nucleus_metadata.csv.gz", `raw file` = "[FASTQ at RIKEN — to be listed]")
}
for (s in sections) rows[[length(rows) + 1]] <- data.table(
  `Sample name` = s, title = paste("Visium prefrontal cortex", s), `source name` = "human brain, post-mortem",
  organism = "Homo sapiens", `characteristics: diagnosis` = if (grepl("NHD", s)) "Nasu-Hakola disease" else "control",
  `characteristics: genotype` = if (grepl("NHD", s)) g("TYROBP / DAP12 genotype", "NHD donor") else "TYROBP wild-type",
  `characteristics: sex` = "female", `characteristics: age` = if (grepl("NHD", s)) g("Age at death (years)", "NHD donor") else g("Age at death (years)", "Control donor"),
  `characteristics: post-mortem interval` = if (grepl("NHD", s)) g("Post-mortem interval", "NHD donor") else g("Post-mortem interval", "Control donor"),
  `characteristics: tissue` = "prefrontal cortex", molecule = "polyA RNA", `library strategy` = "Visium spatial (Space Ranger 1.3.1)",
  `processed data file` = sprintf("Visium/%s/filtered_feature_bc_matrix.h5; Visium/%s/spatial/; Visium/%s/tissue_call.csv", s, s, s), `raw file` = "[FASTQ at RIKEN — to be listed]")   # tissue_call.csv ships per section
fwrite(rbindlist(rows, fill = TRUE), file.path(OUT, "GEO_samples_metadata.csv"))

# ---- README + checksums ---------------------------------------------------------
files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
files <- files[!grepl("md5sums.txt$|README.md$", files)]
md5 <- md5sum(files)
writeLines(sprintf("%s  %s", md5, sub(paste0(OUT, "/"), "", names(md5), fixed = TRUE)), file.path(OUT, "md5sums.txt"))
writeLines(c("# NHD frontal + hippocampus — processed data for deposition",
             sprintf("Built %s by scripts/102_export_deposition_bundle_FH.R.", format(Sys.time(), "%Y-%m-%d %H:%M")),
             "", "## snRNA/",
             "matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz — raw UMI counts, 10x MTX layout, 55,354 post-QC nuclei (the analysed atlas; QC = nFeature_RNA > 200 & < 6,500, percent mitochondrial < 10, scDblFinder singlets).",
             "nucleus_metadata.csv.gz — one row per nucleus: library (SampleID = TFHS id), condition, region, matter, cell type, neuron_subtype (the neuronal subclass/subtype used in the figures: frontal excitatory = Jorstad 2023 dorsolateral prefrontal cortex label transfer, frontal inhibitory = Azimuth human motor-cortex subclass, hippocampal = de-novo subtype) with neuron_subtype_source, azimuth_subclass (+ azimuth_subclass_score; the raw Azimuth motor-cortex call for every nucleus), dfc_within_score (Jorstad DLPFC transfer prediction.score.max, frontal excitatory nuclei only), QC metrics, Harmony-UMAP coordinates, cluster.",
             sprintf("library_sample_sheet.csv — per-library sequencing and QC metrics (identical to %s).", { .e <- new.env(); sys.source(file.path(PROJ, "scripts", "_sd_map_FH.R"), envir = .e); get("sd_ref", .e)("ST12") }),   # number rendered from the fixed map, never typed
             "", "## Visium/",
             # All four sections carry the three-tier tissue call
             "Per section: the raw (all-barcode) Space Ranger count matrix, the filtered in-tissue matrix, metrics_summary.csv and the spatial/ folder (H&E images, scale factors, barcode positions with the in-tissue call used in the analysis). All four sections carry a two-tier tissue call: Space Ranger's image-based call, plus barcodes whose 55-um footprint lies >= 50 % on a relaxed H&E tissue mask and that contain cortical neuropil (SNAP25/RBFOX3/SLC17A7 >= 12 per 10,000; this recovers the pale layer-I rim that automatic tissue detection had excluded while rejecting the detached leptomeninges outside the cortex); the per-barcode tier is in tissue_call.csv (column tissue_call_source = spaceranger | he_mask | excluded, with the delivered call in_tissue_spaceranger, the analysis call in_tissue_analysis and neuropil_per10k) and the thresholds and counts in TISSUE_CALL_provenance.json. spot_metadata.csv.gz additionally flags the spots re-assigned to the L1/pia domain by the BANKSY-guided refinement (histology_refined_L1) and the 14 NHD_Frontal2 spots under an exogenous dark particulate contaminant of ~350 um at the section's lower-right corner (contaminant_particle). The raw matrix is the unmodified Space Ranger output; the filtered matrix is built from it with that tissue call. For NHD_Frontal2 the barcode map is additionally reflected about the array centre (Methods); it carries in_tissue_registered.csv (per barcode: the delivered and the H&E-derived tissue call) and REGISTRATION_provenance.json. spot_metadata.csv.gz carries the analysed object's spot annotations and cell2location abundances.",
             "", "## GEO_samples_metadata.csv",
             "One row per library/section in GEO 'samples' layout; raw-file columns are placeholders until the FASTQs are listed from RIKEN.",
             "", "Checksums: md5sums.txt"),
           file.path(OUT, "README.md"))
say("bundle size: %.2f GB in %d files", sum(file.size(files)) / 1e9, length(files))
cat("\n== sessionInfo ==\n"); print(sessionInfo())
cat("\n=== DONE ===\n", file = stderr())
