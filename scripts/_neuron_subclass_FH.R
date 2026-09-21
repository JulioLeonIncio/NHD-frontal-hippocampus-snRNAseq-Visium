# =============================================================================
# _neuron_subclass_FH.R — Single source for the per-nucleus neuron subclass / subtype label used by every Figure 5 / supplementary script.
# -----------------------------------------------------------------------------
# The cortical (Frontal) excitatory labels were switched from the Azimuth
# human motor-cortex reference (which has no L4 IT class) to the labels transferred
# from the Jorstad 2023 dorsolateral prefrontal cortex reference (scripts/116b,
# `dfc_within`). Frontal inhibitory nuclei keep their Azimuth subclass; hippocampal
# neurons keep the de-novo taxonomy of 60_neuron_reannotation_FH.R. All three live in
# one file, data/neuron_subtype_map_FH.rds (written by 60), and every consumer must
# read the label from that map — never from atlas `predicted.subclass`, which is now
# only the Azimuth QC column (kept in the deposition export as `azimuth_subclass`).
#
# Sourced by: 45_mast_subclass_FH.R, 47b_full_expressed_subclass_FH.R,
#   64_fig4_neuron_umap_FH.R, 70_suppfig1_NHD_QC_FH.R, 90_build_supplementary_tables_FH.R
#   (ST10), 102_export_deposition_bundle_FH.R, 4A / twin / 46 (subclass lists).
# Requires nothing but the map; resolves the rebuild root from `PROJ` / `DIR` if the
# caller defined them, else from the two known machine paths.
#
# Exports:
#   SUBCLASS_MAP         named chr, barcode -> neuron_subtype, all neurons in the map
#   SUBCLASS_SOURCE_MAP  named chr, barcode -> source_annot (provenance per nucleus)
#   CTX_SUBCLASSES       the 15 cortical names in fixed depth order (9 Ex + 6 Inh)
#   CTX_EX_SUBCLASSES / CTX_INH_SUBCLASSES  the two halves of CTX_SUBCLASSES
#   CTX_SUBCLASS_PAL     named hex palette for the 15 cortical names (64 + 70 share it;
#                        L4 IT = #8DD3C7, Sst Chodl = #B3B3B3 — the former 70 collision fixed)
#   CTX_SUBCLASS_SOURCE  one-line provenance string for logs / table notes
#   NEURON_MAP_PATH / NEURON_MAP_MTIME  for cache-invalidation guards
#   subclass_of(barcodes)  chr vector of labels (NA where a barcode is not a neuron
#                          in the map); names = barcodes
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
.nsc_root <- local({
  cand <- character(0)
  if (exists("PROJ", inherits = TRUE) && is.character(get("PROJ", inherits = TRUE)))
    cand <- c(cand, get("PROJ", inherits = TRUE), file.path(get("PROJ", inherits = TRUE), "NHD_frontal_hippo_rebuild"))
  if (exists("DIR", inherits = TRUE) && is.character(get("DIR", inherits = TRUE)))
    cand <- c(cand, get("DIR", inherits = TRUE))
  cand <- c(cand,
    Sys.getenv("NHD_PROJ"),
    Sys.getenv("NHD_PROJ"))
  ok <- cand[file.exists(file.path(cand, "data", "neuron_subtype_map_FH.rds"))]
  if (!length(ok)) stop("_neuron_subclass_FH.R: MISSING data/neuron_subtype_map_FH.rds under any candidate root ",
                        "(run scripts/60_neuron_reannotation_FH.R first). Tried: ", paste(cand, collapse = " | "))
  ok[1]
})
NEURON_MAP_PATH  <- file.path(.nsc_root, "data", "neuron_subtype_map_FH.rds")
NEURON_MAP_MTIME <- file.mtime(NEURON_MAP_PATH)

.nsc_map <- readRDS(NEURON_MAP_PATH)
stopifnot("neuron_subtype_map_FH.rds lacks required columns" =
            all(c("barcode","Region","neuron_subtype","Condition","source_annot") %in% colnames(.nsc_map)),
          "duplicated barcode in neuron_subtype_map_FH.rds" = !anyDuplicated(.nsc_map$barcode),
          "NA neuron_subtype in neuron_subtype_map_FH.rds" = !anyNA(.nsc_map$neuron_subtype))

SUBCLASS_MAP        <- setNames(as.character(.nsc_map$neuron_subtype), .nsc_map$barcode)
SUBCLASS_SOURCE_MAP <- setNames(as.character(.nsc_map$source_annot),   .nsc_map$barcode)

# Fixed cortical order: excitatory by cortical depth (IT first, then the L5/L6
# projection classes), then the six GABAergic subclasses in the shipped order
# . This is the legend / axis /
# table order everywhere; scripts subset it (e.g. drop L5 ET, Sst Chodl for n) but
# never re-order it.
CTX_EX_SUBCLASSES  <- c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L5 ET","L5/6 NP","L6 CT","L6b")
CTX_INH_SUBCLASSES <- c("Pvalb","Sst","Sst Chodl","Vip","Lamp5","Sncg")
CTX_SUBCLASSES     <- c(CTX_EX_SUBCLASSES, CTX_INH_SUBCLASSES)

# Cortical subclass palette (deterministic; from 64_fig4_neuron_umap_FH.R, restricted to
# the 15 names). 70_suppfig1 used to give #8DD3C7 to Sst Chodl, colliding with L4 IT —
# both scripts now read this palette.
CTX_SUBCLASS_PAL <- c(
  "L2/3 IT" = "#1F78C8", "L4 IT" = "#8DD3C7", "L5 IT"     = "#E31A1C", "L6 IT" = "#FF7F00",
  "L6 IT Car3" = "#00CED1", "L5 ET" = "#FFD700", "L5/6 NP" = "#33A02C", "L6 CT" = "#6A33C2",
  "L6b" = "#B15928",
  "Lamp5" = "#A6CEE3", "Pvalb" = "#C8308C", "Sncg" = "#CAB2D6", "Sst" = "#565656",
  "Sst Chodl" = "#B3B3B3", "Vip" = "#F781BF")
stopifnot("CTX_SUBCLASS_PAL must cover CTX_SUBCLASSES exactly" = setequal(names(CTX_SUBCLASS_PAL), CTX_SUBCLASSES),
          "CTX_SUBCLASS_PAL has a duplicated hue" = !anyDuplicated(toupper(CTX_SUBCLASS_PAL)))
CTX_SUBCLASS_PAL <- CTX_SUBCLASS_PAL[CTX_SUBCLASSES]

CTX_SUBCLASS_SOURCE <- paste0(
  "Frontal excitatory subclass = Jorstad 2023 DLPFC label transfer (Seurat TransferData, ",
  "WithinArea_subclass; scripts/116b, tables/L4_DFC_transfer_percell_FH.csv); ",
  "Frontal inhibitory subclass = Azimuth human motor-cortex subclass; ",
  "hippocampal subtype = de-novo (60_neuron_reannotation_FH.R). Map: data/neuron_subtype_map_FH.rds (",
  format(NEURON_MAP_MTIME, "%Y-%m-%d"), ").")

# Guard: the map's Frontal labels must be exactly the 15 cortical names (fail loud if
# 60 was re-run with a different label set, rather than silently greying out a class).
.nsc_ctx_seen <- unique(.nsc_map$neuron_subtype[.nsc_map$Region == "Frontal"])
if (length(setdiff(.nsc_ctx_seen, CTX_SUBCLASSES)))
  stop("_neuron_subclass_FH.R: Frontal labels outside CTX_SUBCLASSES: ",
       paste(setdiff(.nsc_ctx_seen, CTX_SUBCLASSES), collapse = ", "))
if (!"L4 IT" %in% .nsc_ctx_seen)
  stop("_neuron_subclass_FH.R: map carries NO 'L4 IT' Frontal nuclei — 60_neuron_reannotation_FH.R ",
       "was not re-run with the DFC transfer (expected the 2026-09-18 map).")

subclass_of <- function(barcodes) {
  out <- unname(SUBCLASS_MAP[as.character(barcodes)])
  names(out) <- as.character(barcodes)
  out
}
rm(.nsc_ctx_seen)
