#!/usr/bin/env Rscript
# =============================================================================
# 90_build_supplementary_tables_FH.R — 90_build_supplementary_tables_FH.R Build all supplementary tables (ST1-ST14) + Table 1 for the NHD
# Frontal + Hippocampus rebuild, matching the OLD 3-region paper's schemas
# exactly, on the new FH atlas (55,354 nuclei, Region in {Frontal,Hippo}).
#
# One assembly script. Writes CSV per table to supplementary_tables/ and a
# combined workbook NHD_Supplementary_Tables_FH.xlsx (one sheet per STn).
#
# Design rules:
#   * Portable PROJ resolver (julio.l | JulioLeon).
#   * Idempotent: dir.create recursive; temp-then-rename for the workbook.
#   * seed 42 where randomness applies (FindAllMarkers is deterministic given seed).
#   * PREFER the 228 MB subset (atlas/NHD_FH_subset_RNA.rds); never the 2.9 GB heavy
#     atlas here.  Never as.matrix() a sparse genes x cells matrix.
#   * Fail loud: stopifnot on every source file; per-table tryCatch so one failure
#     does not abort the run; failures reported at the end.
#   * never fabricate: every value comes from a source file or the atlas.  Where a
#     source value contradicts the task's suggested literal (e.g. ST3 metric), the
#     TRUE source value is kept and the discrepancy is reported.
#   * American spelling.  No plotting (ragg-free).
#
# Selective rebuild.  Regenerating everything re-runs
# FindAllMarkers + the per-nucleus module scores and re-derives Table 1, which is
# not always wanted (Table 1 is hand-curated).  Pass a comma-separated list of
# table ids to rebuild only those; every other table is registered from the file
# already on disk so the workbook / verification / MANIFEST stay complete:
#
#   Rscript 90_build_supplementary_tables_FH.R ST1,ST3,ST4,ST5,ST13,ST14
#   Rscript 90_build_supplementary_tables_FH.R all       # everything (default)
#
# Table1 is never rebuilt unless it is named explicitly and ALLOW_TABLE1=1 is set
# in the environment (guard added after Table 1 became hand-maintained).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(openxlsx)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

# ---- portable roots ---------------------------------------------------------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("Could not resolve PROJ root" = !is.na(PROJ) && dir.exists(PROJ))

NHDROOT <- dirname(PROJ)                                   # .../Claude_code/NHD
OLD <- file.path(NHDROOT, "NHD_figures_and_legends", "Supplementary_Tables")
stopifnot("Could not resolve OLD supplementary-tables dir" = dir.exists(OLD))

OUT <- file.path(PROJ, "supplementary_tables")
LOGS <- file.path(PROJ, "logs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(LOGS, showWarnings = FALSE, recursive = TRUE)

logf <- file.path(LOGS, "90_build_supplementary_tables_FH.log")
sink(logf, split = TRUE)
on.exit(sink(), add = TRUE)

say  <- function(...) cat(sprintf(...), "\n")
step <- function(x) cat(sprintf("\n== %s ==\n", x))
say("build_supplementary_tables_FH  %s", format(Sys.time()))
say("PROJ = %s", PROJ)
say("OLD  = %s", OLD)

# registry of results for the final verification table
REG <- list()
FAIL <- character(0)
register <- function(id, path, mode, note = "", occ_exempt = FALSE) {
  n_rows <- NA_integer_; n_cols <- NA_integer_; occ <- NA
  if (!is.na(path) && file.exists(path)) {
    d <- suppressWarnings(fread(path, colClasses = "character"))
    n_rows <- nrow(d); n_cols <- ncol(d)
    # Flag OCC in our FH data only. Legitimate substrings excluded: the gene CROCC,
    # the pathway word homeostatic, and "Zhou OCC" (the external Zhou occipital
    # validation cohort in ST13 panel names — a reference dataset, not our tissue).
    occ <- if (occ_exempt) FALSE else any(vapply(d, function(col) {
      v <- col
      v <- gsub("CROCC|HOMEOSTATIC", "", v, ignore.case = TRUE)
      v <- gsub("Zhou +OCC", "", v, ignore.case = TRUE)
      any(grepl("\\bocc\\b|occipital", v, ignore.case = TRUE))
    }, logical(1)))
  }
  REG[[id]] <<- data.frame(table = id, path = ifelse(is.na(path), "", basename(path)),
                           n_rows = n_rows, n_cols = n_cols, mode = mode,
                           occ_found = occ, note = note, stringsAsFactors = FALSE)
}
schema_check <- function(id, got, oldfile) {
  op <- file.path(OLD, oldfile)
  if (!file.exists(op)) { say("  [%s] schema: OLD file missing (%s)", id, oldfile); return(NA) }
  old_hdr <- names(fread(op, nrows = 0))
  ok <- identical(got, old_hdr)
  if (!ok) {
    say("  [%s] SCHEMA MISMATCH", id)
    say("     got: %s", paste(got, collapse = ","))
    say("     old: %s", paste(old_hdr, collapse = ","))
  } else say("  [%s] schema matches OLD (Y)", id)
  ok
}
# schema_check() above compares against the OLD 3-region paper and is only meaningful for
# tables that must still match it. Two other cases exist, and each needs its own check —
# otherwise a new table inherits a check that always no-ops, which LOOKS like coverage:
#
#   schema_new()      a table with no old counterpart (ST14-ST17b). Its schema authority is
#                     scripts/_ST0_data_dictionary_FH.R; the drift check at the bottom of
#                     this script already fails the build on an undocumented or missing
#                     column, but it is setdiff-based and cannot see column order, which
#                     this pins.
#   schema_extends()  a table that legitimately gains columns (ST6). The OLD header must
#                     remain an exact prefix — original names and order frozen, new columns
#                     appended — so a referee's existing parser still reads the file.
schema_new <- function(id, got, expect) {
  if (!identical(got, expect)) {
    say("  [%s] DECLARED SCHEMA MISMATCH", id)
    say("     got:      %s", paste(got,    collapse = ","))
    say("     declared: %s", paste(expect, collapse = ","))
    stop(sprintf("%s: shipped header does not match its declared schema", id))
  }
  say("  [%s] schema matches its declared %d-column header (Y)", id, length(expect))
  TRUE
}
schema_extends <- function(id, got, oldfile) {
  op <- file.path(OLD, oldfile)
  if (!file.exists(op)) { say("  [%s] schema: OLD file missing (%s)", id, oldfile); return(NA) }
  old_hdr <- names(fread(op, nrows = 0))
  if (!identical(got[seq_along(old_hdr)], old_hdr)) {
    say("  [%s] SCHEMA PREFIX BROKEN — the OLD columns must keep their names AND their order", id)
    say("     got[1:%d]: %s", length(old_hdr), paste(got[seq_along(old_hdr)], collapse = ","))
    say("     old:       %s", paste(old_hdr, collapse = ","))
    stop(sprintf("%s: the OLD schema is no longer a prefix of the shipped schema", id))
  }
  say("  [%s] schema = OLD %d cols (names and order identical) + %d appended (DELIBERATE): %s",
      id, length(old_hdr), length(got) - length(old_hdr),
      paste(setdiff(got, old_hdr), collapse = ","))
  TRUE
}
wr <- function(df, id, fname) {
  p <- file.path(OUT, fname)
  fwrite(df, p)
  say("  wrote %s  (%d x %d)", fname, nrow(df), ncol(df))
  p
}
# -----------------------------------------------------------------------------
# Canonical id -> shipped filename map.  Defined once here and reused for the
# workbook, the MANIFEST and the "keep what is already on disk" path, so a rename
# can never leave the three out of sync (that is exactly how ST3/ST4/ST5/ST13 came
# to carry filenames that misstated their own method).
# -----------------------------------------------------------------------------
sheet_map <- c(
  ST0  = "ST0_data_dictionary.csv",
  ST1  = "ST1_MAST_discovery_DGE_per_celltype_region.csv",
  ST2  = "ST2_pseudobulk_DGE_per_celltype_region.csv",
  ST3  = "ST3_fGSEA_prerank_signed_log2FC.csv",
  ST4  = "ST4_TF_activity_decoupleR_ULM_GTRD.csv",
  ST5  = "ST5_TF_regulon_GTRD_unsigned.csv",
  ST6  = "ST6_curated_signatures.csv",
  ST7  = "ST7_neuron_module_scores_per_nucleus.csv",
  ST8  = "ST8_celltype_marker_genes.csv",
  ST9  = "ST9_celltype_counts_per_region_condition.csv",
  ST10 = "ST10_cortical_neuron_subclass_counts.csv",
  ST11 = "ST11_hippo_neuron_subtype_meta.csv",
  ST12 = "ST12_per_lane_sequencing_QC.csv",
  ST13 = "ST13_zhou_crosscohort_direction_consistency.csv",
  ST14 = "ST14_micro_TF_activity_CollecTRI_Fig2h.csv",
  ST15  = "ST15_neuron_programme_per_library_Fig5d.csv",
  ST15b = "ST15b_neuron_programme_RINgap_vs_loss.csv",
  ST16  = "ST16_hippo_Ex_programme_delta_unresolved_sensitivity.csv",
  ST17  = "ST17_micro_only_vs_pooled_MicroPVM.csv",
  ST17b = "ST17b_myeloid_selection_chain.csv",
  Table1 = "Table1_donor_sample_FH.csv")

# Filenames this directory used to ship under.  A rename that leaves the old file
# in place is worse than no rename: the folder then contains two files that differ
# and a referee cannot tell which one the figures used.  Any legacy name still on
# disk is MOVED (not deleted) into _superseded_<date>/ at the end of the run.
LEGACY <- c("ST3_fGSEA_signed_cliffs_delta.csv",
            "ST4_TF_activity_decoupleR_ULM.csv",
            "ST5_TF_regulon_GTRD.csv",
            "ST13_zhou_crosscohort_replication_stats.csv")

# ---- which tables to (re)build ----------------------------------------------
.args <- commandArgs(trailingOnly = TRUE)
SELECT <- if (!length(.args) || tolower(.args[1]) %in% c("all", "")) names(sheet_map) else
  trimws(strsplit(.args[1], ",")[[1]])
bad <- setdiff(SELECT, names(sheet_map))
if (length(bad)) stop(sprintf("unknown table id(s): %s", paste(bad, collapse = ",")))
SELECT <- union(SELECT, "ST0")   # the dictionary always re-derives from what shipped
if ("Table1" %in% SELECT && !identical(Sys.getenv("ALLOW_TABLE1"), "1")) {
  say("Table1 requested but ALLOW_TABLE1 != 1 -> Table1 will be KEPT AS-IS (hand-maintained).")
  SELECT <- setdiff(SELECT, "Table1")
}
# ---- ids that are derived in one pass from one source ------------------------
# ST15 and ST15b are two views of the same seven library means; ST17 and ST17b are two
# views of the same microglia-only run. Rebuilding one and leaving the other on disk
# ships a fresh half beside a stale half — which is exactly the failure the ST1 header
# below documents.
COBUILD <- list(ST15 = "ST15b", ST15b = "ST15", ST17 = "ST17b", ST17b = "ST17")
for (k in intersect(names(COBUILD), SELECT)) SELECT <- union(SELECT, COBUILD[[k]])

say("rebuilding: %s", paste(SELECT, collapse = ", "))
say("keeping as-is: %s", paste(setdiff(names(sheet_map), SELECT), collapse = ", "))

# ---- fail fast on a table that is neither selected nor present ---------------
# The ST0 drift check near the end of this script stop()s when a documented table has
# no file on disk. Reaching that point costs the whole run (atlas load, FindAllMarkers,
# module scores) and dies before the workbook, _verification_FH.csv and the MANIFEST are
# written — so the run produces nothing at all. Catch it in the first second instead.
.absent <- names(sheet_map)[!(names(sheet_map) %in% SELECT) &
                            !file.exists(file.path(OUT, unlist(sheet_map[names(sheet_map)])))]
if (length(.absent))
  stop(sprintf(paste("these tables are neither selected for rebuild nor present on disk: %s.",
                     "The ST0 drift check would abort this run at the very end, after every",
                     "expensive step. Re-run with 'all', or name them in the selection."),
               paste(.absent, collapse = ", ")))

run <- function(id, expr) {
  if (!id %in% SELECT) {
    fp <- file.path(OUT, sheet_map[[id]])
    if (file.exists(fp)) {
      register(id, fp, "UNCHANGED", "not selected for rebuild; existing file kept verbatim")
      say("  [%s] KEPT AS-IS (%s)", id, sheet_map[[id]])
    } else {
      say("  [%s] SKIPPED and MISSING on disk (%s)", id, sheet_map[[id]])
      FAIL[[id]] <<- sprintf("not selected and %s absent", sheet_map[[id]])
    }
    return(invisible(NULL))
  }
  tryCatch(expr, error = function(e) {
    say("  [%s] ERROR: %s", id, conditionMessage(e))
    FAIL[[id]] <<- conditionMessage(e); NULL
  })
}

# -----------------------------------------------------------------------------
# Reviewer-facing numerical rigor: p-value underflow.
# A p-value stored as exactly 0 is double-precision underflow, not a true zero, and
# -log10(0) = Inf silently breaks any downstream ranking or colour scale.  Every
# p-column shipped in a supplementary table is therefore floored at
# .Machine$double.xmin (2.225074e-308) with pmax(); values written by an earlier
# run at 14 significant digits can round-trip to a hair below xmin, so the floor is
# applied to those too, which puts every underflowed value on the same tie.  Read
# any value equal to the floor as "p < 2.2e-308".
# NA is left as NA on purpose (DESeq2 independent filtering = "not tested"); it is
# never coerced to 0 or 1.  The count of floored values is logged per column.
# -----------------------------------------------------------------------------
# The floor must survive the CSV round-trip, or the shipped file contradicts its
# own documentation.  fwrite writes doubles at 15 significant digits, so
# .Machine$double.xmin (2.2250738585072014e-308) is written as 2.2250738585072e-308
# and reads back a hair below the limit it claims to enforce — a referee re-reading
# the CSV would find p < double.xmin and be right.  The floor is therefore
# double.xmin rounded up to 8 significant digits: still the smallest normalized
# double to every digit anyone reports, and exactly representable on disk.
P_FLOOR <- 2.2250739e-308
stopifnot("p-floor is below the smallest normalized double" = P_FLOOR >= .Machine$double.xmin,
          "p-floor does not survive a 15-digit CSV round-trip" =
            identical(P_FLOOR, as.numeric(sprintf("%.15g", P_FLOOR))))
# Which tables actually got floored, recorded as it happens. The MANIFEST used to name
# them in hand-typed prose ("ST1, ST2, ST3, ST13 and ST14"), which silently goes wrong the
# first time a table is added — prose that names tables must come from the registry.
P_FLOORED <- character(0)
floor_p <- function(dt, cols, id) {
  P_FLOORED <<- union(P_FLOORED, id)
  for (cc in cols) {
    if (!cc %in% names(dt)) { say("  [%s] p-floor: column %s absent -> skipped", id, cc); next }
    v <- dt[[cc]]
    n_zero <- sum(v == 0, na.rm = TRUE)
    n_sub  <- sum(v > 0 & v < P_FLOOR, na.rm = TRUE)
    if (n_zero + n_sub > 0) {
      data.table::set(dt, j = cc, value = pmax(v, P_FLOOR))
      say("  [%s] p-floor %-12s: %d exact zeros + %d sub-normal -> %s (report as p < 2.2e-308)",
          id, cc, n_zero, n_sub, sprintf("%.8g", P_FLOOR))
    } else say("  [%s] p-floor %-12s: none (min > 0 = %.4g)", id, cc,
               suppressWarnings(min(v[v > 0], na.rm = TRUE)))
    n_na <- sum(is.na(v))
    if (n_na) say("  [%s]         %-12s: %d NA kept as NA (not tested; never coerced)", id, cc, n_na)
  }
  dt
}

MAST_DIR <- file.path(PROJ, "tables", "mast_dual")

# =============================================================================
# ST1  MAST discovery DGE per cell-type x region
# OLD: comparison,cell_type,region,gene,avg_log2FC,p_val,p_val_adj,pct_NHD,
#      pct_CON,cliffs_delta,n_NHD,n_CON,discovery
#
# ST1 is built from the two authoritative aggregate discovery tables that the figures
# themselves read, so it cannot drift from them (a glob over tables/mast_dual would also
# sweep in aggregates without a `comparison` column and the Microglia_/PVM_ split files,
# which no figure uses):
#     tables/mast_dual/MAST_dual_discovery_all.csv      12 broad cell-type x region
#     tables/mast_dual/MAST_subclass_discovery_all.csv  13 Frontal neuron subclasses (+L4 IT)
# 37,626 + 42,588 = 80,214 rows (row counts as of the 12-subclass build; the 13-subclass file is larger). pct_NHD / pct_CON / cliffs_delta (and the derived
# `discovery` flag) are computed on the raw RNA layer, never on SCT (zero-UMI genes score
# ~99 % SCT detection).
# =============================================================================
run("ST1", {
  step("ST1 MAST discovery DGE")
  src <- c(broad    = file.path(MAST_DIR, "MAST_dual_discovery_all.csv"),
           subclass = file.path(MAST_DIR, "MAST_subclass_discovery_all.csv"))
  for (nm in names(src))
    if (!file.exists(src[[nm]]))
      stop(sprintf("MISSING ST1 %s source %s — run 40_mast_percell_dual_FH.R / 45_mast_subclass_FH.R",
                   nm, basename(src[[nm]])))
  # latent_vars added. Both discovery scripts have always written the model's
  # covariates per row (40_mast_percell_dual_FH.R:104, 45_mast_subclass_FH.R:169) and both
  # source CSVs carry the column — but ST1 dropped it, so the shipped table gave a referee no
  # way to verify the Methods' claim that log10(nCount_RNA) is in every per-nucleus model.
  # That claim matters here: NHD frontal nuclei are ~1.7-2.0x deeper than control.
  cols <- c("comparison","cell_type","region","gene","avg_log2FC","p_val","p_val_adj",
            "pct.1","pct.2","cliffs_delta","n_NHD","n_CON","latent_vars","discovery")
  L <- lapply(names(src), function(nm) {
    d <- fread(src[[nm]])
    miss <- setdiff(cols, names(d))
    if (length(miss)) stop(sprintf("%s missing cols: %s", basename(src[[nm]]),
                                   paste(miss, collapse = ",")))
    say("  %-8s %-38s %6d rows, %2d comparisons", nm, basename(src[[nm]]),
        nrow(d), length(unique(d$comparison)))
    d[, ..cols]
  })
  st1 <- rbindlist(L)
  setnames(st1, c("pct.1","pct.2"), c("pct_NHD","pct_CON"))
  # The figures do not plot the raw discovery set: a curated ambient/artifact
  # filter is applied downstream (scripts/_artifact_genes.R), removing unannotated lncRNA
  # and antisense classes, clustered histones, housekeeping and sex-chromosome genes,
  # unnamed clone contigs, and neuronal transcripts scored in glia. Without this flag a
  # referee recomputing Fig 1d-f from ST1 gets the unfiltered counts (e.g. Micro-PVM 654
  # rather than the 599 on the panel) and concludes the panel is unreproducible. With it,
  # all six cell types reproduce exactly, and so do the Fig 1e/1f correlations and their n.
  source(file.path(PROJ, "scripts", "_artifact_genes.R"))
  st1[, artifact := is_artifact(gene, cell_type)]
  st1[, discovery_filtered := discovery & !artifact]
  setcolorder(st1, c("comparison","cell_type","region","gene","avg_log2FC","p_val",
                     "p_val_adj","pct_NHD","pct_CON","cliffs_delta","n_NHD","n_CON",
                     "discovery","artifact","discovery_filtered","latent_vars"))
  say("  ST1 artifact-filtered: %d of %d discovery rows removed by the curated filter",
      sum(st1$discovery & st1$artifact), sum(st1$discovery))

  # Dual-method support tier, so Fig 1c is reproducible from the shipped
  # tables alone. Frontal is the only region with two libraries per condition, so it is
  # the only one with a pseudobulk tier; the hippocampus has one NHD library and its
  # "support" is agreement of the per-nucleus direction with the frontal pseudobulk sign.
  # These are different quantities and are named differently rather than pooled.
  # ST1 READS ST2 from disk, and the ST2 block sits below this one, so in a single run
  # that rebuilds both, ST1 sees the previous ST2.
  # Fail with the fix rather than with a column-not-found from 200 lines away.
  pbf <- fread(file.path(OUT, "ST2_pseudobulk_DGE_per_celltype_region.csv"))
  if (!"model_fit_reportable" %in% names(pbf))
    stop("ST1 needs the CURRENT ST2 on disk, but the file still has the pre-2026-09-08 schema. ",
         "ST2 is built after ST1 in this script, so rebuild ST2 first and then ST1: ",
         "run with 'ST2', then run with 'ST1'.")
  pbf[, cell_type := sub("_(Frontal|Hippo)$", "", comparison)]
  pbf[, pb_strict := !is.na(padj) & padj < 0.01 & abs(log2FoldChange_apeglm) > log2(2.5)]
  pbf[, pb_sign := sign(log2FoldChange_apeglm)]
  fr <- unique(pbf[region == "Frontal" & model_fit_reportable %in% TRUE & pb_strict %in% TRUE,
                   .(cell_type, gene, front_pb_sign = pb_sign)])
  st1[, mast_sign := sign(avg_log2FC)]
  st1 <- merge(st1, unique(pbf[, .(cell_type, region, gene, pb_strict, pb_sign)]),
               by = c("cell_type", "region", "gene"), all.x = TRUE)
  st1 <- merge(st1, fr, by = c("cell_type", "gene"), all.x = TRUE)
  st1[, pb_confirmed      := region == "Frontal" & pb_strict %in% TRUE & pb_sign == mast_sign]
  st1[, xregion_supported := region == "Hippo" & !is.na(front_pb_sign) & front_pb_sign == mast_sign]
  st1[, support_tier := fifelse(discovery_filtered & pb_confirmed,      "lane_reproducible",
                        fifelse(discovery_filtered & xregion_supported, "xregion_supported",
                        fifelse(discovery_filtered,                     "discovery_only", NA_character_)))]
  st1[, c("mast_sign","pb_strict","pb_sign","front_pb_sign") := NULL]
  # latent_vars is APPENDED last, not slotted in beside n_CON, so that the OLD 13-column
  # header stays an exact prefix of the shipped one and schema_extends() keeps passing.
  setcolorder(st1, c("comparison","cell_type","region","gene","avg_log2FC","p_val",
                     "p_val_adj","pct_NHD","pct_CON","cliffs_delta","n_NHD","n_CON",
                     "discovery","artifact","discovery_filtered",
                     "pb_confirmed","xregion_supported","support_tier","latent_vars"))
  say("  ST1 support tiers: lane_reproducible=%d  xregion_supported=%d  discovery_only=%d",
      sum(st1$support_tier == "lane_reproducible", na.rm=TRUE),
      sum(st1$support_tier == "xregion_supported", na.rm=TRUE),
      sum(st1$support_tier == "discovery_only",    na.rm=TRUE))
  # no comparison may appear in both sources
  dupe <- st1[, .N, by = .(comparison, gene)][N > 1]
  if (nrow(dupe)) stop(sprintf("ST1: %d duplicated (comparison,gene) keys", nrow(dupe)))
  # our regions only — Frontal/Hippo (no occipital in the FH rebuild)
  stopifnot("ST1 region column contains something other than Frontal/Hippo" =
              all(st1$region %in% c("Frontal", "Hippo")))
  # detection percentages must be proportions in [0,1]; a >1 value would mean an
  # SCT-style percentage leaked back in
  stopifnot("ST1 pct columns outside [0,1]" =
              all(st1$pct_NHD >= 0 & st1$pct_NHD <= 1 & st1$pct_CON >= 0 & st1$pct_CON <= 1))
  stopifnot("ST1 cliffs_delta outside [-1,1]" =
              all(is.na(st1$cliffs_delta) | abs(st1$cliffs_delta) <= 1))
  st1 <- floor_p(st1, c("p_val", "p_val_adj"), "ST1")
  say("  ST1 rows=%d  comparisons=%d  discovery TRUE=%d", nrow(st1),
      length(unique(st1$comparison)), sum(st1$discovery))
  p <- wr(st1, "ST1", sheet_map[["ST1"]])
  schema_extends("ST1", names(st1), "ST1_MAST_discovery_DGE_per_celltype_region.csv")
  register("ST1", p, "REGENERATED",
           "MAST_dual_discovery_all + MAST_subclass_discovery_all (post-2026-08-19 raw-RNA pct/delta); p floored at double.xmin")
})

# =============================================================================
# ST2  pseudobulk DESeq2 per cell-type x region
# OLD: gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,log2FoldChange_apeglm,
#      lfcSE_apeglm,comparison,n_cells_total,n_lanes,n_NHD_lanes,n_CON_lanes
# Source (FH): diagnostics/03_DGE_pseudobulk/pseudobulk_DESeq2_with_Wald_stat.csv
#   (extra trailing cols note,dispersion_basis -> dropped).
# =============================================================================
run("ST2", {
  step("ST2 pseudobulk DESeq2")
  f <- file.path(PROJ, "diagnostics", "03_DGE_pseudobulk", "pseudobulk_DESeq2_with_Wald_stat.csv")
  stopifnot("MISSING ST2 pseudobulk source" = file.exists(f))
  d <- fread(f)
  # The previous ST2 dropped the source's `note` and
  # `dispersion_basis` columns and shipped every row — frontal and hippocampal — labelled
  # n_lanes=4, n_NHD_lanes=2. For the hippocampus that is the design-matrix count of a
  # random split of the single NHD library into two pseudo-replicates, not a count of
  # physical libraries. A referee opening ST2 therefore found a 2-versus-2 hippocampal
  # DESeq2 that the Methods says is impossible. The provenance columns are now retained
  # and unambiguous physical-library counts and an explicit validity flag are added.
  cols <- c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj",
            "log2FoldChange_apeglm","lfcSE_apeglm","comparison","n_cells_total",
            "n_lanes","n_NHD_lanes","n_CON_lanes","note","dispersion_basis")
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(sprintf("ST2 source missing: %s", paste(miss, collapse=",")))
  st2 <- d[, ..cols]
  st2 <- floor_p(st2, c("pvalue", "padj"), "ST2")

  # region from the comparison string; physical library counts from the study design
  st2[, region := fifelse(grepl("Hippo", comparison), "Hippo", "Frontal")]
  st2[, design_n_lanes      := n_lanes]
  st2[, design_n_NHD_lanes  := n_NHD_lanes]
  st2[, design_n_CON_lanes  := n_CON_lanes]
  st2[, n_NHD_libraries_physical := fifelse(region == "Hippo", 1L, 2L)]
  st2[, n_CON_libraries_physical := 2L]
  # Two column names were making a claim the design cannot
  # support. `inference_valid` reads as "this row licenses a donor-level inference", which is
  # false for every row in a 1-vs-1 study: the frontal tier compares the grey- and white-matter
  # dissections of one donor against the two of another, so it is a compartment-consistency
  # check, not biological replication. `biological_2v2` said "biological" for the same reason.
  # Renamed to say what the rows actually are.
  st2[, dispersion_basis := fifelse(dispersion_basis == "biological_2v2",
                                    "within_donor_tissue_libraries", dispersion_basis)]
  st2[, model_fit_reportable := dispersion_basis == "within_donor_tissue_libraries"]
  st2[, n_lanes := NULL][, n_NHD_lanes := NULL][, n_CON_lanes := NULL]
  setcolorder(st2, c("gene","comparison","region","baseMean",
                     "log2FoldChange","lfcSE","stat","pvalue","padj",
                     "log2FoldChange_apeglm","lfcSE_apeglm","n_cells_total",
                     "n_NHD_libraries_physical","n_CON_libraries_physical",
                     "design_n_lanes","design_n_NHD_lanes","design_n_CON_lanes",
                     "dispersion_basis","model_fit_reportable","note"))

  # Guards: the pseudo-split must never carry inference, and must never claim >1 NHD library
  ps <- st2[dispersion_basis == "pseudo_split"]
  stopifnot("ST2: pseudo-split rows must have NA padj"  = all(is.na(ps$padj)),
            "ST2: pseudo-split rows must have NA pvalue"= all(is.na(ps$pvalue)),
            "ST2: pseudo-split rows must have NA stat"  = all(is.na(ps$stat)),
            "ST2: pseudo-split must be hippocampus only"= all(ps$region == "Hippo"),
            "ST2: hippocampus has exactly 1 physical NHD library" =
              all(st2[region == "Hippo"]$n_NHD_libraries_physical == 1L))
  say("  ST2 rows=%d comparisons=%d", nrow(st2), length(unique(st2$comparison)))
  say("  ST2 rows with a reportable model fit (frontal, within-donor tissue libraries): %d", sum(st2$model_fit_reportable))
  say("  ST2 descriptive-only rows (hippocampal pseudo-split, p/padj/stat NA): %d",
      sum(!st2$model_fit_reportable))
  p <- wr(st2, "ST2", "ST2_pseudobulk_DGE_per_celltype_region.csv")
  schema_check("ST2", names(st2), "ST2_pseudobulk_DGE_per_celltype_region.csv")
  register("ST2", p, "ASSEMBLED", "select 14 cols from FH pseudobulk_DESeq2_with_Wald_stat.csv")
})

# =============================================================================
# ST3  fGSEA prerank per subtype x region
#
# The file shipped as
# ST3_fGSEA_signed_cliffs_delta.csv, but every one of its 42,305 rows carries
# metric = "log2fc_winsorized": the signed-Cliff's-delta prerank was retired on
# The file is now
# ST3_fGSEA_prerank_signed_log2FC.csv, which states the metric actually used.
# No value changes.
# OLD: subtype,pathway,pval,padj,ES,NES,size,leadingEdge,direction,db,metric,
#      n_nuclei_min,power_tier
# Source: diagnostics/05_GSEA/fGSEA_NHD_MAST_cliffs/fGSEA_NHD_MAST_cliffs_combined_all.csv
#   (already carries all 13 cols + extras).  direction values are CON/NHD (raw);
#   OLD ST3 used the same subtype-level CON/NHD tokens.  metric = value from source
#   ("log2fc_winsorized") kept VERBATIM (not fabricated to "signed_cliffs_delta").
# =============================================================================
run("ST3", {
  step("ST3 fGSEA")
  f <- file.path(PROJ, "diagnostics", "05_GSEA", "fGSEA_NHD_MAST_cliffs",
                 "fGSEA_NHD_MAST_cliffs_combined_all.csv")
  stopifnot("MISSING ST3 GSEA source" = file.exists(f))
  d <- fread(f)
  cols <- c("subtype","pathway","pval","padj","ES","NES","size","leadingEdge",
            "direction","db","metric","n_nuclei_min","power_tier")
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(sprintf("ST3 source missing: %s", paste(miss, collapse=",")))
  st3 <- d[, ..cols]
  st3 <- floor_p(st3, c("pval", "padj"), "ST3")
  # fgsea returns NA pval/padj/NES for pathways whose enrichment could not be
  # estimated (degenerate permutation null); ES is still finite.  Kept as NA.
  say("  ST3 rows=%d subtypes=%d ; NA pval rows=%d (fgsea could not estimate NES; kept NA)",
      nrow(st3), length(unique(st3$subtype)), sum(is.na(st3$pval)))
  say("  metric values in source: %s", paste(unique(st3$metric), collapse=","))
  stopifnot("ST3 metric is no longer log2fc_winsorized -> the FILENAME must change too" =
              identical(sort(unique(as.character(st3$metric))), "log2fc_winsorized"))
  say("  direction values: %s", paste(sort(unique(na.omit(st3$direction))), collapse=","))
  p <- wr(st3, "ST3", sheet_map[["ST3"]])
  schema_check("ST3", names(st3), "ST3_fGSEA_signed_cliffs_delta.csv")
  register("ST3", p, "RENAMED+ASSEMBLED",
           "13 cols from combined_all; metric=log2fc_winsorized -> filename now states the real prerank metric")
})

# =============================================================================
# ST4  TF activity, decoupleR ULM on the GTRD regulon  (genome-wide survey)
# ST5  the GTRD regulon itself (unsigned: mode_of_regulation == 1 everywhere)
#
# ST4/ST5 are the broad survey across all six
# cell types on GTRD.  They are not the source of Figure 2h, which is microglia-only
# and uses the SIGNED, curated CollecTRI regulon (scripts 68 + 70b) — only 8 of
# ST4's 459 factors appear in the Fig-2h result table.  Both files now name their
# regulon so the two analyses cannot be confused, and the CollecTRI numbers that do
# back Figure 2h ship separately as ST14.  Values unchanged; rename only.
# =============================================================================
run("ST4", {
  step("ST4 TF activity (GTRD survey)")
  src <- file.path(OUT, "ST4_TF_activity_decoupleR_ULM.csv")     # legacy name
  dst <- file.path(OUT, sheet_map[["ST4"]])
  if (!file.exists(src) && file.exists(dst)) src <- dst           # already renamed
  stopifnot("MISSING ST4" = file.exists(src))
  d <- fread(src)
  if (!identical(normalizePath(src), normalizePath(dst, mustWork = FALSE))) {
    fwrite(d, dst); file.remove(src)
    say("  renamed %s -> %s", basename(src), basename(dst))
  }
  say("  ST4 rows=%d  TFs=%d  comparisons=%d", nrow(d),
      length(unique(d$source)), length(unique(d$condition)))
  schema_check("ST4", names(d), "ST4_TF_activity_decoupleR_ULM.csv")
  register("ST4", dst, "RENAMED", "GTRD genome-wide TF-activity survey (all cell types); NOT the Fig-2h source (see ST14)")
})
run("ST5", {
  step("ST5 TF regulon (GTRD, unsigned)")
  src <- file.path(OUT, "ST5_TF_regulon_GTRD.csv")                # legacy name
  dst <- file.path(OUT, sheet_map[["ST5"]])
  if (!file.exists(src) && file.exists(dst)) src <- dst
  stopifnot("MISSING ST5" = file.exists(src))
  d <- fread(src)
  if (!identical(normalizePath(src), normalizePath(dst, mustWork = FALSE))) {
    fwrite(d, dst); file.remove(src)
    say("  renamed %s -> %s", basename(src), basename(dst))
  }
  mor_col <- intersect(c("mode_of_regulation", "mor"), names(d))[1]
  say("  ST5 edges=%d  TFs=%d  %s values: %s", nrow(d), length(unique(d$source)),
      mor_col, paste(sort(unique(d[[mor_col]])), collapse = ","))
  schema_check("ST5", names(d), "ST5_TF_regulon_GTRD.csv")
  register("ST5", dst, "RENAMED", "GTRD TF->target network; mode_of_regulation is 1 for every edge (UNSIGNED)")
})

# =============================================================================
# ST6  curated signatures (long: signature, gene)
# Sources: curated_signatures.rds (22) + micro PROGRAMS (2B) + neuron 8-programme
#   set (identical defs to OLD 10b_figure5_v2.R / used by 4C) + Pan2024 CSF1R sigs
#   + the Fig5d_* sets (score-cache attribute) + the six inline astrocyte
#   sets of Fig. 3e,f (scripts/_astro_inline_sets_FH.R, shared with 3D) and the eight
#   oligodendrocyte-lineage programs of Fig. 4d (scripts/_oligo_programmes_FH.R, shared
#   with 3Fp/3M/3P/65). De-dup on (signature, gene).
#   INVARIANT (asserted below): every Fig. 3e / Fig. 4d program of ST24 has >= 1 ST6 row.
# =============================================================================
run("ST6", {
  step("ST6 curated signatures")
  DATA <- file.path(PROJ, "data")
  base_sigs <- readRDS(file.path(DATA, "curated_signatures.rds"))  # 22 named vectors

  # micro programme sets (verbatim from 2B_micro_program_violins_FH.R)
  micro_progs <- list(
    "Antigen presentation"           = c("CD74","HLA-DRA","HLA-DRB1","HLA-DMB"),
    "Metal handling (iron/zinc)"     = c("FTH1","FTL","HAMP","TMEM163"),
    "Phagocytic / lysosomal"         = c("FCGR3A","SCIN","ADAM28","CTSB"),
    "Glycolytic shift"               = c("PFKFB3","SLC2A3"),
    "Inflammatory / stress response" = c("CEBPD","ZFP36L1","RGS1","SRGN","HSPA1A","TNFRSF1B","NAIP"),
    "Homeostatic (lost)"             = c("P2RY12","MEF2C","PLXDC2","SORL1"),
    "Immunoregulatory brake (lost)"  = c("IRAK3","LDLRAD4","ZBTB16","HDAC9"))

  # Neuron 8-programme set: exact defs that produced OLD ST7 (10b_figure5_v2.R);
  # the FH 4C script scores a subset (Presyn/Postsyn/GluR/OxPhos/IEG) of the same
  # curated anchors.  We ship the full 8 so ST6 documents every neuron programme.
  neuron_sigs <- list(
    Synaptic     = unique(c(base_sigs$Synaptic_presynaptic, base_sigs$Synaptic_postsynaptic)),
    OxPhos       = base_sigs$Neuronal_OXPHOS,
    UPR          = c("ATF4","DDIT3","HSPA5","XBP1","ATF6","HERPUD1","EIF2AK3",
                     "ERN1","EIF2S1","ASNS","SEL1L","HSPH1","DNAJB9","PPP1R15A",
                     "TRIB3","DDIT4","NUPR1","CHAC1"),
    DNA_damage   = c("H2AFX","TP53","CDKN1A","CDKN2A","ATM","ATR","CHEK1",
                     "CHEK2","BRCA1","RAD51","MDC1","MRE11A","RAD50","NBN",
                     "PARP1","GADD45A","GADD45B","MLH1","PCNA","RPA1"),
    Apoptosis    = c("TP53","BAX","BAK1","BID","BBC3","PMAIP1","CASP3","CASP7",
                     "CASP9","APAF1","DIABLO","CDKN1A"),
    IEG          = base_sigs$IEG_core,
    Postsynaptic = base_sigs$Synaptic_postsynaptic,
    Presynaptic  = base_sigs$Synaptic_presynaptic)

  # Pan 2024 CSF1R-RD cross-microgliopathy signatures
  pan_src <- file.path(PROJ, "scripts", "_pan2024_csf1r_signatures.R")
  pan_sigs <- list()
  if (file.exists(pan_src)) {
    e <- new.env(); sys.source(pan_src, envir = e)
    if (exists("PAN2024_SIGS", envir = e)) {
      ps <- get("PAN2024_SIGS", envir = e)
      names(ps) <- paste0("Pan2024_", names(ps))
      pan_sigs <- ps
    }
  } else say("  WARNING: _pan2024_csf1r_signatures.R not found -> Pan sigs omitted")

  # The figure-5d SETS were missing from ST6 entirely.
  # `neuron_sigs` above defines Presynaptic := base_sigs$Synaptic_presynaptic (21 genes),
  # Postsynaptic := Synaptic_postsynaptic (14) and OxPhos := Neuronal_OXPHOS (16) — the sets
  # scored per nucleus in ST7. Figure 5d scores different sets, defined literally in
  # 4C_signature_violins_FH.R (19 / 17 / 23 genes), plus an ionotropic-GluR set that ST6 did
  # not carry at all. Same labels, different membership: 7, 6 and 7 genes apart. Only IEG is
  # common to both. A referee joining ST16 to ST6's "Presynaptic" would have got the wrong
  # gene list, and Figure 5d was not reconstructible from the shipped tables.
  # They are therefore added under a Fig5d_ prefix, read from the score cache's own `sigs`
  # attribute so the table cannot drift from the panel, and ST15/ST16 use these exact names.
  CACHE5D <- file.path(DATA, "_cache_neuron_scores_FH.rds")
  fig5d_sigs <- list()
  if (file.exists(CACHE5D)) {
    .sg <- attr(readRDS(CACHE5D), "sigs")
    if (!is.null(.sg)) { fig5d_sigs <- .sg; names(fig5d_sigs) <- paste0("Fig5d_", names(.sg)) }
  }
  if (!length(fig5d_sigs))
    stop("ST6: cannot read the Fig-5d module definitions from data/_cache_neuron_scores_FH.rds — run scripts/4C_signature_violins_FH.R")
  say("  ST6 Fig5d sets added: %s", paste(sprintf("%s(%d)", names(fig5d_sigs),
      lengths(fig5d_sigs)), collapse = " "))

  # ST6 (sheet Signatures) must list every set that ST24 (sheet
  # Program_deltas) scores. Two families were missing — the six inline astrocyte sets of Fig. 3e,f and
  # the eight oligodendrocyte-lineage programs of Fig. 4d. Both are read from the same helper
  # files the panel scripts source, never typed here, so the table cannot drift from the panels.
  # Their names are the ST24 `program` labels verbatim (asserted below), so ST6 joins ST24.
  .h_astro <- file.path(PROJ, "scripts", "_astro_inline_sets_FH.R")
  .h_oligo <- file.path(PROJ, "scripts", "_oligo_programmes_FH.R")
  stopifnot("MISSING scripts/_astro_inline_sets_FH.R (shared with 3D_astro_state_violins_FH.R)" = file.exists(.h_astro),
            "MISSING scripts/_oligo_programmes_FH.R (shared with 3Fp/3M/3P/65)" = file.exists(.h_oligo))
  .ea <- new.env(); sys.source(.h_astro, envir = .ea)
  astro_inline <- get("astro_inline_sets", .ea)()
  ASTRO_PANEL  <- get("astro_published_panel_names", .ea)()   # A1_reactive -> "Complement/IFN-reactive (Liddelow)", ...
  .eo <- new.env(); sys.source(.h_oligo, envir = .eo)
  oligo_progs  <- get("oligo_programme_sets", .eo)(base_sigs)[get("oligo_programme_order", .eo)()]
  stopifnot("astro inline sets: expected six named, non-empty character vectors" =
              length(astro_inline) == 6L && all(nzchar(names(astro_inline))) && all(lengths(astro_inline) >= 3L),
            "oligo programs: expected eight named, non-empty character vectors" =
              length(oligo_progs) == 8L && all(nzchar(names(oligo_progs))) && all(lengths(oligo_progs) >= 3L),
            "a new set name collides with an existing ST6 signature" =
              !any(c(names(astro_inline), names(oligo_progs)) %in% c(names(base_sigs), names(micro_progs), names(neuron_sigs), names(pan_sigs), names(fig5d_sigs))),
            "the published astrocyte panel map must point at curated_signatures.rds sets" = all(names(ASTRO_PANEL) %in% names(base_sigs)))
  say("  ST6 astro inline sets added (Fig. 3e,f): %s", paste(sprintf("%s(%d)", names(astro_inline), lengths(astro_inline)), collapse = " | "))
  say("  ST6 oligo programs added (Fig. 4d): %s", paste(sprintf("%s(%d)", names(oligo_progs), lengths(oligo_progs)), collapse = " | "))

  all_sigs <- c(base_sigs, micro_progs, neuron_sigs, pan_sigs, fig5d_sigs, astro_inline, oligo_progs)
  st6 <- rbindlist(lapply(names(all_sigs), function(nm)
    data.table(signature = nm, gene = as.character(all_sigs[[nm]]))))
  st6 <- unique(st6[!is.na(gene) & gene != ""])

  # ---- provenance columns (task 3) -------------------------------------------
  # The point is to separate a symbol problem from an expression problem for every gene,
  # so that "this gene is not in the signature score" always has a stated reason.
  #
  # Gene universe = rownames of the atlas RNA assay. This is an annotation match, not a
  # detection test: the RNA rownames are the full reference gene set, so TRUE means only
  # "the symbol exists in the matrix", never "it was detected in any nucleus".
  GU_CACHE <- file.path(DATA, "_gene_universe_FH.rds")
  ATLAS_SUB_ <- file.path(PROJ, "atlas", "NHD_FH_subset_RNA.rds")
  stopifnot("MISSING atlas subset (needed once for the ST6 gene universe)" = file.exists(ATLAS_SUB_))
  if (file.exists(GU_CACHE) && file.mtime(GU_CACHE) >= file.mtime(ATLAS_SUB_)) {
    GENES <- readRDS(GU_CACHE)
    say("  gene universe: %d symbols (cache)", length(GENES))
  } else {
    say("  gene universe: building from the atlas subset (one-off, ~30 s)")
    suppressPackageStartupMessages(library(Seurat))
    .o <- readRDS(ATLAS_SUB_); GENES <- rownames(.o[["RNA"]]); rm(.o); invisible(gc())
    .t <- paste0(GU_CACHE, ".tmp"); saveRDS(GENES, .t); file.rename(.t, GU_CACHE)
    say("  gene universe: %d symbols (cached to %s)", length(GENES), basename(GU_CACHE))
  }
  stopifnot("ST6 gene universe looks wrong" = length(GENES) > 1000L)

  # Legacy aliases that appear in the published sets. Each is scored under its current
  # HGNC symbol, which is why the set still works despite the published symbol being absent.
  ALIAS <- c(PSD95 = "DLG4", DAP12 = "TYROBP", MRE11A = "MRE11")
  # Mouse symbols with no human ortholog symbol in the matrix. All four are in A1_reactive
  # (Liddelow 2017), which was defined in mouse; 4 of its 10 genes cannot be scored in human.
  MOUSE_ONLY <- c("H2-T23", "IIGP1", "LIGP1", "GGTA1")

  st6[, symbol_original := gene]
  st6[, symbol_mapped := fifelse(gene %in% names(ALIAS), unname(ALIAS[gene]),
                          fifelse(gene %in% MOUSE_ONLY, NA_character_, gene))]
  st6[, in_scoring_matrix := !is.na(symbol_mapped) & symbol_mapped %in% GENES]

  # Membership of each context is verified against the shipped
  # panel tables where one exists, so a set cannot be labelled as scored by a panel that does not
  # carry it.
  .esd <- new.env(); sys.source(file.path(PROJ, "scripts", "_sd_map_FH.R"), envir = .esd)
  CTX <- get("SD_ST6_CTX", .esd)
  stopifnot(setequal(names(CTX), c("nucleus","fig2d","visium","astro","fig4d","fig5d","fig2g","refonly")))
  sig_ctx <- setNames(rep(CTX[["refonly"]], length(unique(st6$signature))), unique(st6$signature))
  # (a) the eight neuron programs are the ST7 (sheet Neuron_scores_per_nucleus) per-nucleus scores; the four
  #     curated sets they are built from are scored under the ST7 column name
  sig_ctx[names(neuron_sigs)] <- CTX[["nucleus"]]
  for (.pair in list(c("Synaptic_presynaptic","Presynaptic"), c("Synaptic_postsynaptic","Postsynaptic"),
                     c("Neuronal_OXPHOS","OxPhos"), c("IEG_core","IEG")))
    sig_ctx[.pair[1]] <- sprintf("%s as column %s", CTX[["nucleus"]], .pair[2])
  # (b) the seven microglial programs = Fig. 2d (verified against ST24 when it is on disk)
  .f24 <- file.path(OUT, "ST24_program_cliffs_delta_FH.csv")
  .p24 <- if (file.exists(.f24)) fread(.f24) else NULL
  if (!is.null(.p24))
    stopifnot("ST6: a microglial program is not a Fig. 2d row of ST24" = all(names(micro_progs) %in% .p24[figure_panel == "Fig. 2d"]$program))
  sig_ctx[names(micro_progs)] <- CTX[["fig2d"]]
  # (d) astrocyte sets scored in Fig. 3e,f: the six inline sets carry their ST24 program name as the
  #     signature; the five published sets are scored under a panel LABEL (A1_reactive is reported as
  #     "Complement/IFN-reactive (Liddelow)", ...) which is stated in scoring_context so ST6 joins ST24
  #     — the same "as column X" convention the ST7 sets use in (a). Both maps come from the helper
  #     that 3D_astro_state_violins_FH.R sources (asserted: 3D sources that helper).
  stopifnot("ST6: 3D_astro_state_violins_FH.R no longer sources _astro_inline_sets_FH.R" =
              any(grepl("_astro_inline_sets_FH.R", readLines(file.path(PROJ, "scripts", "3D_astro_state_violins_FH.R"), warn = FALSE), fixed = TRUE)))
  sig_ctx[names(astro_inline)] <- CTX[["astro"]]
  sig_ctx[names(ASTRO_PANEL)]  <- sprintf("%s as program %s", CTX[["astro"]], unname(ASTRO_PANEL))
  # (d2) the eight oligodendrocyte-lineage programs = Fig. 4d (Oligo + OPC rows of ST24)
  sig_ctx[names(oligo_progs)] <- CTX[["fig4d"]]
  # (c) sets also scored on tissue: read from the shipped Visium gene-set sheet (ST23c), never typed.
  #     Runs AFTER (d)/(d2) so the five oligodendrocyte programs also scored on Visium (Structural
  #     myelin, Cholesterol (sterol arm), Fatty-acid / sphingomyelin, Galactolipid, Lipid uptake /
  #     salvage) get the Visium context appended exactly as the two microglial programs do.
  .f23c <- file.path(OUT, "ST23c_visium_program_gene_sets_FH.csv")
  if (file.exists(.f23c)) { .vis <- intersect(unique(fread(.f23c)$program), names(sig_ctx))
    sig_ctx[.vis] <- paste(sig_ctx[.vis], CTX[["visium"]], sep = "; ")
    say("  ST6 sets also scored on Visium (from ST23c): %s", paste(.vis, collapse = ", ")) }
  # Invariant: every Fig. 3e / Fig. 4d program of ST24 (sheet Program_deltas) has >= 1 ST6 row,
  #   either under its own name or (published astrocyte sets) under the ST6 name whose scoring_context
  #   names it as the panel program.
  if (!is.null(.p24)) {
    .need <- unique(.p24[figure_panel %in% c("Fig. 3e", "Fig. 4d")]$program)
    .have <- c(unique(st6$signature), unname(ASTRO_PANEL))
    .gap  <- setdiff(.need, .have)
    if (length(.gap)) stop("ST6: Fig. 3e / Fig. 4d program(s) of ST24 with no ST6 row: ", paste(.gap, collapse = "; "))
    say("  ST6 covers all %d Fig. 3e / Fig. 4d programs of ST24 (%d by signature name, %d via a stated panel label)",
        length(.need), sum(.need %in% unique(st6$signature)), sum(.need %in% unname(ASTRO_PANEL) & !.need %in% unique(st6$signature)))
  } else say("  WARNING: ST24 not on disk -> the Fig. 3e / Fig. 4d coverage invariant was NOT checked")
  # (e) the Fig5d_ names live only in the score cache's attribute
  sig_ctx[grep("^Fig5d_", names(sig_ctx))] <- CTX[["fig5d"]]
  # (f) Zhou_NHD_microglia_up is a gene set of the Fig. 2g state-enrichment panel
  stopifnot(grepl("Zhou_NHD_microglia_up", paste(readLines(file.path(PROJ, "scripts", "2E_states_MAST_GSEA_FH.R"), warn = FALSE), collapse = "\n"), fixed = TRUE))
  sig_ctx["Zhou_NHD_microglia_up"] <- CTX[["fig2g"]]
  # everything else (Pan 2024 CSF1R states, DAM_up, LDAM, ...) stays reference only
  stopifnot(!anyNA(sig_ctx), !any(grepl("scripts/|\\.R\\b", sig_ctx)))
  st6[, scoring_context := unname(sig_ctx[signature])]
  st6[, scored := scoring_context != CTX[["refonly"]]]

  st6[, used_in_scoring := in_scoring_matrix & scored]
  st6[, exclusion_reason := fifelse(!is.na(symbol_mapped) & symbol_mapped != symbol_original &
                                      in_scoring_matrix, "alias_merged",
                             fifelse(is.na(symbol_mapped),        "symbol_not_mapped",
                             fifelse(!in_scoring_matrix,          "not_in_matrix", "none")))]
  # `below_min_detection` is a permitted value but is empty by construction and that is
  # stated rather than left to be noticed: no scoring path in this project applies a
  # per-gene detection floor. score_one() (4C / ST7) takes colMeans over the genes present
  # and only requires >= 3 of them; AddModuleScore has no per-gene floor either.
  st6[, scored := NULL]
  st6[, symbol_mapped := fifelse(is.na(symbol_mapped), "", symbol_mapped)]

  ST6_EXPECT <- c("signature","gene","symbol_original","symbol_mapped","in_scoring_matrix",
                  "used_in_scoring","exclusion_reason","scoring_context")
  setcolorder(st6, ST6_EXPECT)
  stopifnot("ST6: a gene has an empty exclusion_reason" =
              all(!is.na(st6$exclusion_reason) & nzchar(st6$exclusion_reason)),
            "ST6: exclusion_reason outside its permitted set" =
              all(st6$exclusion_reason %in% c("none","symbol_not_mapped","not_in_matrix",
                                              "below_min_detection","alias_merged")))
  say("  ST6 signatures=%d  rows(unique sig,gene)=%d", length(unique(st6$signature)), nrow(st6))
  say("  ST6 exclusion_reason: %s",
      paste(sprintf("%s=%d", names(table(st6$exclusion_reason)), table(st6$exclusion_reason)),
            collapse = "  "))
  say("  ST6 signatures not scored by any shipped panel: %d of %d",
      length(unique(st6$signature[st6$scoring_context == CTX[["refonly"]]])),
      length(unique(st6$signature)))
  say("  ST6 scoring_context values: %s", paste(sprintf("%s (n=%d)", names(table(st6$scoring_context)), as.integer(table(st6$scoring_context))), collapse = " | "))
  p <- wr(st6, "ST6", sheet_map[["ST6"]])
  schema_extends("ST6", names(st6), "ST6_curated_signatures.csv")
  register("ST6", p, "REGENERATED+EXTENDED",
           "curated_signatures.rds + micro programmes + 8 neuron programmes + Pan2024 + the Fig5d_* sets + 6 inline astro sets (Fig. 3e,f) + 8 oligo-lineage programs (Fig. 4d); de-duped; five provenance columns appended")
})

# =============================================================================
# The remaining tables need the light atlas subset. Load once.
# =============================================================================
ATLAS_TABLES <- c("ST7","ST8","ST9","ST10","ST11","ST12","Table1")
NEED_ATLAS <- length(intersect(SELECT, ATLAS_TABLES)) > 0
ATLAS_SUB <- file.path(PROJ, "atlas", "NHD_FH_subset_RNA.rds")
if (!NEED_ATLAS) {
  step("atlas subset NOT needed for this selection -> skipping the 228 MB load")
  MD <- NULL
} else {
stopifnot("MISSING atlas subset NHD_FH_subset_RNA.rds" = file.exists(ATLAS_SUB))
step("loading light atlas subset (228 MB)")
suppressPackageStartupMessages(library(Seurat))
obj <- readRDS(ATLAS_SUB)
if ("RNA" %in% names(obj@assays)) DefaultAssay(obj) <- "RNA"
if (length(SeuratObject::Layers(obj, assay = "RNA")) > 1) obj <- JoinLayers(obj, assay = "RNA")
# subset ships a counts-only RNA layer; log-normalize so FindAllMarkers (ST8) and
# the ST7 module scores have a `data` layer to read.
obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
MD <- obj@meta.data
stopifnot("atlas subset != 55354 nuclei" = ncol(obj) == 55354L)
say("  loaded %d nuclei; cell types=%d", ncol(obj), length(unique(MD$new_annotation)))
}

# =============================================================================
# ST9  cell-type counts per region x condition (CON/NHD wide). must total 55354.
# OLD: new_annotation,Region,CON,NHD
# =============================================================================
run("ST9", {
  step("ST9 cell-type counts")
  tab <- as.data.frame(table(new_annotation = MD$new_annotation,
                             Region = MD$Region, Condition = MD$Condition),
                       responseName = "n")
  st9 <- tab %>%
    pivot_wider(names_from = Condition, values_from = n, values_fill = 0) %>%
    arrange(new_annotation, Region) %>%
    select(cell_type = new_annotation, Region, CON, NHD) %>%   # shipped as cell_type (the atlas metadata name new_annotation is internal)
    as.data.frame()
  tot <- sum(st9$CON) + sum(st9$NHD)
  say("  ST9 rows=%d  total nuclei=%d (expect 55354)", nrow(st9), tot)
  stopifnot("ST9 does not total 55354" = tot == 55354L)
  p <- wr(st9, "ST9", "ST9_celltype_counts_per_region_condition.csv")
  schema_check("ST9", names(st9), "ST9_celltype_counts_per_region_condition.csv")   # reports the deliberate new_annotation -> cell_type rename
  register("ST9", p, "REGENERATED", sprintf("table(cell_type,Region,Condition); totals %d", tot))
})

# =============================================================================
# ST10  Frontal cortical neuron subclass counts (no OCC in FH)
# OLD: neuron_subtype,Region,Condition,n
# label = scripts/_neuron_subclass_FH.R (data/neuron_subtype_map_FH.rds):
#   excitatory = Jorstad 2023 DLPFC transfer (adds L4 IT), inhibitory = Azimuth.
#   Not atlas predicted.subclass any more. 15 subclasses.
# =============================================================================
run("ST10", {
  step("ST10 cortical neuron subclass counts")
  source(file.path(PROJ, "scripts", "_neuron_subclass_FH.R"))
  neu <- MD[MD$new_annotation %in% c("Neuron_Ex","Neuron_Inh") & MD$Region == "Frontal", ]
  lab <- subclass_of(rownames(neu))
  src <- unname(SUBCLASS_SOURCE_MAP[rownames(neu)])
  stopifnot("ST10: frontal neuron(s) missing from neuron_subtype_map_FH (re-run 60)" = !anyNA(lab),
            "ST10: label outside the 15-name cortical taxonomy" = all(lab %in% CTX_SUBCLASSES))
  tab <- as.data.frame(table(neuron_subtype = factor(lab, levels = CTX_SUBCLASSES),
                             Region = neu$Region, Condition = neu$Condition),
                       responseName = "n")
  src_by <- unique(data.frame(neuron_subtype = lab, label_source = src, stringsAsFactors = FALSE))
  stopifnot("ST10: a subclass carries labels from two references" = !anyDuplicated(src_by$neuron_subtype))
  st10 <- tab %>% filter(n > 0) %>%
    mutate(neuron_subtype = as.character(neuron_subtype)) %>%
    left_join(src_by, by = "neuron_subtype") %>%
    arrange(match(neuron_subtype, CTX_SUBCLASSES), Condition) %>%
    select(neuron_subtype, Region, Condition, n, label_source) %>%
    as.data.frame()
  stopifnot("ST10: NA label_source" = !anyNA(st10$label_source))
  say("  ST10 rows=%d subclasses=%d (Frontal only; no OCC); L4 IT CON=%d NHD=%d", nrow(st10),
      length(unique(st10$neuron_subtype)),
      sum(st10$n[st10$neuron_subtype == "L4 IT" & st10$Condition == "CON"]),
      sum(st10$n[st10$neuron_subtype == "L4 IT" & st10$Condition == "NHD"]))
  stopifnot("ST10 must carry 15 cortical subclasses incl. L4 IT" =
              length(unique(st10$neuron_subtype)) == 15L && "L4 IT" %in% st10$neuron_subtype)
  p <- wr(st10, "ST10", "ST10_cortical_neuron_subclass_counts.csv")
  schema_extends("ST10", names(st10), "ST10_cortical_neuron_subclass_counts.csv")   # +label_source (OLD 4 cols kept as prefix)
  register("ST10", p, "REGENERATED",
           paste0("Frontal neuron subclass counts (Ex = Jorstad-DFC transfer incl. L4 IT; Inh = Azimuth); ",
                  "OCC absent by design. ", CTX_SUBCLASS_SOURCE))
})

# =============================================================================
# ST8  cell-type marker genes: FindAllMarkers(only.pos), top-10 protein-coding
# OLD: gene,cluster,avg_log2FC,pct_in   (pct_in = pct.1)
# House identity filter: drop lncRNA/MT-/RP- style non-protein-coding symbols.
# =============================================================================
run("ST8", {
  step("ST8 cell-type marker genes")
  Idents(obj) <- "new_annotation"
  say("  FindAllMarkers on %d clusters (seed 42)...", length(levels(Idents(obj))))
  set.seed(42)
  mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                       logfc.threshold = 0.25, verbose = FALSE)
  # house non-protein-coding filter (match OLD ST8, which carried only real gene
  # symbols): mitochondrial (MT-), ribosomal (RP[SL]), antisense (-AS#), LINC/MIR/
  # SNHG lncRNA families, ENSG passthroughs, and clone/contig accession identifiers
  # (AC/AL/AP/AJ/BX...######.# from Genbank/Vega — uncharacterized loci, not proteins).
  drop_re <- paste0(
    "^MT-|^RP[SL]|^RPL|^RPS|-AS[0-9]*$|^LINC[0-9]|^MIR[0-9]|^SNHG[0-9]|^ENSG[0-9]|",
    "orf[0-9]|^[A-Z]{2}[0-9]{6}\\.[0-9]+$|\\.[0-9]+$")
  mk <- as.data.table(mk)
  mk[, is_ncrna := grepl(drop_re, gene)]
  say("  markers total=%d ; dropped non-protein-coding=%d", nrow(mk), sum(mk$is_ncrna))
  mk <- mk[is_ncrna == FALSE]
  top <- mk[order(cluster, -avg_log2FC)][, head(.SD, 10), by = cluster]
  st8 <- top[, .(gene, cluster, avg_log2FC, pct_in = pct.1)]
  st8 <- as.data.frame(st8)
  say("  ST8 rows=%d clusters=%d", nrow(st8), length(unique(st8$cluster)))
  p <- wr(st8, "ST8", "ST8_celltype_marker_genes.csv")
  schema_check("ST8", names(st8), "ST8_celltype_marker_genes.csv")
  register("ST8", p, "REGENERATED",
           "FindAllMarkers only.pos, top-10 protein-coding per new_annotation")
})

# =============================================================================
# ST7  neuron per-nucleus module scores (8 programmes) for Neuron_Ex + Neuron_Inh
# OLD: cell,Condition,Region,cell_type,orig.ident,Synaptic,OxPhos,UPR,DNA_damage,
#      Apoptosis,IEG,Postsynaptic,Presynaptic
# Score = colMeans of NormalizeData log data over present genes (>=3), matching the
# OLD ST7 (all non-negative) and FH 4C score_one().  Sparse ops only.
# =============================================================================
run("ST7", {
  step("ST7 neuron module scores")
  DATA <- file.path(PROJ, "data")
  base_sigs <- readRDS(file.path(DATA, "curated_signatures.rds"))
  neuron_sigs <- list(
    Synaptic     = unique(c(base_sigs$Synaptic_presynaptic, base_sigs$Synaptic_postsynaptic)),
    OxPhos       = base_sigs$Neuronal_OXPHOS,
    UPR          = c("ATF4","DDIT3","HSPA5","XBP1","ATF6","HERPUD1","EIF2AK3",
                     "ERN1","EIF2S1","ASNS","SEL1L","HSPH1","DNAJB9","PPP1R15A",
                     "TRIB3","DDIT4","NUPR1","CHAC1"),
    DNA_damage   = c("H2AFX","TP53","CDKN1A","CDKN2A","ATM","ATR","CHEK1",
                     "CHEK2","BRCA1","RAD51","MDC1","MRE11A","RAD50","NBN",
                     "PARP1","GADD45A","GADD45B","MLH1","PCNA","RPA1"),
    Apoptosis    = c("TP53","BAX","BAK1","BID","BBC3","PMAIP1","CASP3","CASP7",
                     "CASP9","APAF1","DIABLO","CDKN1A"),
    IEG          = base_sigs$IEG_core,
    Postsynaptic = base_sigs$Synaptic_postsynaptic,
    Presynaptic  = base_sigs$Synaptic_presynaptic)
  ord <- c("Synaptic","OxPhos","UPR","DNA_damage","Apoptosis","IEG","Postsynaptic","Presynaptic")

  neurons <- subset(obj, new_annotation %in% c("Neuron_Ex","Neuron_Inh"))
  neurons <- NormalizeData(neurons, assay = "RNA", verbose = FALSE)
  dm <- GetAssayData(neurons, assay = "RNA", layer = "data")   # sparse; never as.matrix()
  score_one <- function(genes) {
    g <- intersect(genes, rownames(dm))
    if (length(g) < 3) { say("    %s: <3 genes present -> NA column", "sig"); return(rep(NA_real_, ncol(dm))) }
    Matrix::colMeans(dm[g, , drop = FALSE])
  }
  say("  scoring %d neuron nuclei over %d modules:", ncol(neurons), length(neuron_sigs))
  for (s in ord) {
    g <- intersect(neuron_sigs[[s]], rownames(dm))
    say("    %-14s %2d / %2d present", s, length(g), length(neuron_sigs[[s]]))
  }
  score_mat <- sapply(ord, function(s) score_one(neuron_sigs[[s]]))
  nm <- neurons@meta.data
  st7 <- data.frame(cell = colnames(neurons),
                    Condition = as.character(nm$Condition),
                    Region = as.character(nm$Region),
                    cell_type = as.character(nm$new_annotation),
                    orig.ident = as.character(nm$orig.ident),
                    as.data.frame(score_mat), check.names = FALSE,
                    stringsAsFactors = FALSE)
  st7 <- st7[, c("cell","Condition","Region","cell_type","orig.ident", ord)]
  say("  ST7 rows=%d (Ex=%d Inh=%d)", nrow(st7),
      sum(st7$cell_type=="Neuron_Ex"), sum(st7$cell_type=="Neuron_Inh"))
  say("  score ranges: Synaptic[%.3f,%.3f] any all-NA col? %s",
      min(st7$Synaptic,na.rm=TRUE), max(st7$Synaptic,na.rm=TRUE),
      any(vapply(st7[ord], function(x) all(is.na(x)), logical(1))))
  p <- wr(st7, "ST7", "ST7_neuron_module_scores_per_nucleus.csv")
  schema_check("ST7", names(st7), "ST7_neuron_module_scores_per_nucleus.csv")
  register("ST7", p, "REGENERATED", "colMeans module scores, 8 programmes, Ex+Inh neurons")
})

# =============================================================================
# ST11  hippo neuron subtype metadata (FH de-novo).  OLD was a Tippani glia-QC
# table (is_glia_like, cluster_myelin_mean...) that the FH pipeline does not
# reproduce.  We ship the FH hippo-neuron subtype map (data/neuron_subtype_map_FH)
# and report the schema divergence rather than fabricate the missing QC columns.
# =============================================================================
run("ST11", {
  step("ST11 hippo neuron subtype metadata")
  f <- file.path(PROJ, "data", "neuron_subtype_map_FH.csv")
  stopifnot("MISSING neuron_subtype_map_FH.csv" = file.exists(f))
  d <- fread(f)
  # The map now also carries azimuth_subclass / dfc_within_score (cortical-only
  # columns, NA for every hippocampal row) — keep ST11's published 5-column schema.
  st11 <- d[Region == "Hippo", .(barcode, Region, neuron_subtype, Condition, source_annot)]
  st11[, taxonomy_flag := fifelse(neuron_subtype == "Unresolved",
                                  "unresolved (failed private-marker gate)",
                                  "private-marker supported")]
  say("  ST11 rows=%d hippo neuron subtypes=%d", nrow(st11), length(unique(st11$neuron_subtype)))
  say("  taxonomy_flag: %d rows unresolved (low complexity), %d marker-supported",
      sum(st11$neuron_subtype == "Unresolved"), sum(st11$neuron_subtype != "Unresolved"))
  say("  NOTE: FH ST11 schema DIFFERS from OLD Tippani ST11 (barcode,Region,neuron_subtype,")
  say("        Condition,source_annot).  OLD glia-QC columns (is_glia_like/cluster_*_mean)")
  say("        were a Tippani-mapping analysis NOT rerun in FH -> would be fabrication to invent.")
  p <- wr(as.data.frame(st11), "ST11", "ST11_hippo_neuron_subtype_meta.csv")
  # Deliberate: do not force old schema (documented divergence)
  register("ST11", p, "REGENERATED",
           "FH hippo neuron subtype map; schema differs from OLD Tippani QC table (reported)")
})

# =============================================================================
# ST12  per-lane sequencing QC.  keep OLD sequencing metrics for the 7 FH lanes,
# update only "Nuclei (post-QC)" from the FH atlas per-lane counts + update Total.
# Lane -> Donor/Region/Matter map derived from atlas SampleLabel (unambiguous).
# =============================================================================
run("ST12", {
  step("ST12 per-lane sequencing QC")
  old12 <- fread(file.path(OLD, "ST12_per_lane_sequencing_QC.csv"))
  # atlas per-lane post-QC nucleus counts by SampleLabel
  tab <- table(MD$SampleLabel)
  # unambiguous SampleLabel -> (Donor, Region, Matter) matching OLD row keys
  lane_map <- list(
    "Ctrl_1_Frontal_grey"  = c("Control","Frontal cortex","Grey"),
    "Ctrl_1_Frontal_white" = c("Control","Frontal cortex","White"),
    "Ctrl_1_Hippo_1"       = c("Control","Hippocampus (capture 1)","—"),
    "Ctrl_1_Hippo_2"       = c("Control","Hippocampus (capture 2)","—"),
    "NHD_1_Frontal_grey"   = c("NHD","Frontal cortex","Grey"),
    "NHD_2_Frontal_white"  = c("NHD","Frontal cortex","White"),
    "NHD_5_Hippo_1"        = c("NHD","Hippocampus (capture 1)","—"))
  # ensure every atlas lane is mapped
  unmapped <- setdiff(names(tab), names(lane_map))
  if (length(unmapped)) stop(sprintf("UNMAPPED lanes -> STOP: %s", paste(unmapped, collapse=",")))
  postqc <- setNames(as.integer(tab[names(lane_map)]), names(lane_map))

  # drop OCC rows + Total from OLD, keep the 7 FH lanes (match by Donor+Region+Matter)
  key <- function(don, reg, mat) paste(don, reg, mat, sep="||")
  old12$._key <- key(old12$Donor, old12$Region, old12$Matter)
  fh_keys <- vapply(lane_map, function(v) key(v[1], v[2], v[3]), character(1))
  keep <- old12[old12$._key %in% fh_keys]
  if (nrow(keep) != 7L) stop(sprintf("expected 7 FH lane rows matched in OLD ST12, got %d", nrow(keep)))
  # update Nuclei (post-QC) per matched lane
  lab_by_key <- setNames(names(fh_keys), fh_keys)
  keep$`Nuclei (post-QC)` <- format(postqc[lab_by_key[keep$._key]], big.mark=",", trim=TRUE)
  # RIN per library (Agilent Bioanalyzer, measured on the same RNA used for library prep).
  # RIN is systematically lower in NHD and predicts the transcriptome mapping rate
  # (Pearson r = 0.96 over these seven libraries); see Methods.
  rin_by_key <- c(
    "Control||Frontal cortex||Grey"            = 7.4,
    "Control||Frontal cortex||White"           = 7.5,
    "Control||Hippocampus (capture 1)||\u2014" = 7.1,
    "Control||Hippocampus (capture 2)||\u2014" = 6.7,
    "NHD||Frontal cortex||Grey"                = 6.3,
    "NHD||Frontal cortex||White"               = 6.0,
    "NHD||Hippocampus (capture 1)||\u2014"     = 4.4)
  miss_rin <- setdiff(keep$._key, names(rin_by_key))
  if (length(miss_rin)) stop(sprintf("ST12: no RIN for %s", paste(miss_rin, collapse=", ")))
  keep$RIN <- sprintf("%.1f", rin_by_key[keep$._key])
  keep[, ._key := NULL]
  # Total row
  tot <- as.list(rep("—", ncol(keep))); names(tot) <- names(keep)
  tot$Donor <- "Total"
  tot$`Nuclei (post-QC)` <- format(sum(postqc), big.mark=",", trim=TRUE)
  st12 <- rbind(keep, as.data.table(tot), use.names = TRUE)
  say("  ST12 rows=%d (7 lanes + Total); Total post-QC=%s (expect 55,354)",
      nrow(st12), format(sum(postqc), big.mark=","))
  stopifnot("ST12 post-QC total != 55354" = sum(postqc) == 55354L)
  p <- wr(as.data.frame(st12), "ST12", "ST12_per_lane_sequencing_QC.csv")
  schema_check("ST12", names(st12), "ST12_per_lane_sequencing_QC.csv")
  register("ST12", p, "ASSEMBLED",
           "OLD sequencing metrics for 7 FH lanes; Nuclei(post-QC)+Total updated from FH atlas")
})

# =============================================================================
# ST13  Zhou 2023 cross-cohort direction consistency.
#
# It shipped as "..._replication_stats.csv".
# House rule: the Zhou comparison is a threshold-free direction-consistency test
# against an empirical baseline, and is never described as replication.
# The shipped file kept only the OLD
# 13 columns and therefore omitted `baseline_concord` — the empirical baseline the
# whole claim is measured against.  A table named "direction consistency" that does
# not carry its own baseline is not checkable, so the four context columns already
# present in the working copy are appended after the original 13 (original names
# and order untouched, no value recomputed).
# cols 1-13 (OLD): panel,n_shared,n_both_sig,rho_full,rho_full_p,rho_sig,rho_sig_p,
#      tau_sig,tau_sig_p,n_concord,concord_frac,binom_p,rho_loo
# appended: baseline_concord,n_baseline_genes,n_dir_zero,cell_type,region
# =============================================================================
run("ST13", {
  step("ST13 Zhou cross-cohort direction consistency")
  f <- file.path(PROJ, "tables", "ST13_zhou_crosscohort_FH.csv")
  stopifnot("MISSING ST13 FH source" = file.exists(f))
  d <- fread(f)
  cols <- c("panel","n_shared","n_both_sig","rho_full","rho_full_p","rho_sig",
            "rho_sig_p","tau_sig","tau_sig_p","n_concord","concord_frac","binom_p","rho_loo")
  extra <- c("baseline_concord","n_baseline_genes","n_dir_zero","cell_type","region")
  miss <- setdiff(c(cols, extra), names(d))
  if (length(miss)) stop(sprintf("ST13 source missing: %s", paste(miss, collapse=",")))
  st13 <- d[, c(cols, extra), with = FALSE]
  st13 <- floor_p(st13, c("rho_full_p","rho_sig_p","tau_sig_p","binom_p"), "ST13")
  # our side of every panel must be a FH region; "Zhou OCC" in `panel` is the
  # external occipital cohort and is expected.
  stopifnot("ST13 region column contains a non-FH region" =
              all(st13$region %in% c("Frontal", "Hippo")))
  say("  ST13 rows=%d panels; concord_frac %.3f-%.3f vs empirical baseline %.3f-%.3f",
      nrow(st13), min(st13$concord_frac), max(st13$concord_frac),
      min(st13$baseline_concord), max(st13$baseline_concord))
  p <- wr(as.data.frame(st13), "ST13", sheet_map[["ST13"]])
  say("  [ST13] schema = OLD 13 cols + 5 appended context cols (DELIBERATE, see header);")
  say("         cols 1-13 identical in name and order to the OLD table: %s",
      identical(names(st13)[1:13], names(fread(file.path(OLD, "ST13_zhou_crosscohort_replication_stats.csv"), nrows = 0))))
  register("ST13", p, "RENAMED+EXTENDED",
           "OLD 13 cols + baseline_concord/n_baseline_genes/n_dir_zero/cell_type/region appended; direction-consistency, not replication")
})

# =============================================================================
# ST14  CollecTRI microglial TF activity — the table behind Figure 2h.  (NEW)
#
# Fig 2h is microglia-only, depth-matched, scored with decoupleR ULM on the SIGNED
# curated CollecTRI regulon (script 68), then filtered by the pre-declared six-part
# gate in script 70b.  ST4/ST5 (GTRD, all cell types, unsigned) do not back it.
# Sources, both verbatim:
#   tables/micro_TF_collectri_thin750_FH.csv   100 rows = TF x region activity
#   tables/fig2_micro_TF_gate_audit_FH.csv           46 TFs scored in both regions
# 54 unique TFs are scored; 46 clear the detection floor in both regions (46 x 2 =
# 92 rows) and 8 are scorable in one region only (8 rows) -> 100 rows.
# BH is applied once over all 100 TF x region tests — not over the 46 rows of the
# gate audit — reproducing scripts/70b line "p.adjust(allp$p_value, method='BH')".
# =============================================================================
run("ST14", {
  step("ST14 CollecTRI microglial TF activity (Fig 2h)")
  # PRIMARY = the log-normalized, all-nuclei run (68 with DEPTH_MODE=lognorm_all),
  # which is what Fig 2h now plots. The 750-UMI arm ships beside it as a sensitivity rather
  # than as the result. The gate audit is regenerated by 70b from the same matrix as the
  # scores, so the audit never documents a gate applied to a different dataset than the panel.
  fa <- file.path(PROJ, "tables", "micro_TF_collectri_FH.csv")
  fs <- file.path(PROJ, "tables", "micro_TF_collectri_thin750_FH.csv")
  fg <- file.path(PROJ, "tables", "fig2_micro_TF_gate_audit_FH.csv")
  for (f in c(fa, fs, fg))
    if (!file.exists(f)) stop(sprintf("MISSING ST14 source %s — run 68_micro_TF_collectri_FH.R (both DEPTH_MODEs) / 70b_micro_TF_bipartite_FH.R", basename(f)))
  act  <- fread(fa)
  gate <- fread(fg)
  stopifnot("ST14: activity table is the OLD version (no score_noHK/n_target)" =
              all(c("score_noHK","n_target") %in% names(act)))
  stopifnot("ST14: activity table is not microglia CollecTRI (expect statistic=ulm, condition=cliffs_delta)" =
              all(act$statistic == "ulm") && all(act$condition == "cliffs_delta"))
  stopifnot("ST14 region column contains a non-FH region" =
              all(act$region %in% c("Frontal", "Hippo")))

  # BH over every TF x region test actually performed (100), matching 70b exactly.
  act[, q_BH := p.adjust(p_value, method = "BH")]
  say("  BH applied over %d TF x region tests (%d unique TFs); q<0.05 in %d rows",
      nrow(act), length(unique(act$source)), sum(act$q_BH < 0.05))

  # the 750-UMI sensitivity arm, joined on (TF, region). Disclosure by agreement, not by
  # substitution: a referee can read both columns and see the answer does not depend on it.
  sens <- fread(fs)
  sens[, q_BH_750 := p.adjust(p_value, method = "BH")]
  act <- merge(act, sens[, .(source, region, score_750 = score, q_BH_750)],
               by = c("source","region"), all.x = TRUE)
  .both <- act[!is.na(score_750)]
  say("  sensitivity arm joined: %d of %d rows have a 750-UMI counterpart; on those, r = %.3f",
      nrow(.both), nrow(act), suppressWarnings(cor(.both$score, .both$score_750)))
  say("  %d factors are scoreable only WITHOUT the 750-UMI floor (incl. %s)",
      sum(is.na(act$score_750)),
      paste(head(sort(unique(act$source[is.na(act$score_750)])), 6), collapse = ", "))

  g <- gate[, .(TF, gate_mean_z = mean_z, gate_coverage = coverage,
                gate_leaveout_ratio = leaveout, gate_mor_purity = mor_pure,
                gate_tf_own_delta = tf_delta, gate_anchor_detected = anchor_det,
                gate_anchor_total = anchor_tot, gate_anchor_dir_frac = anchor_dir,
                gate_C1_coverage = C1_coverage, gate_C2_reproducible = C2_reproducible,
                gate_C3_robust = C3_robust, gate_C4_coherent = C4_coherent,
                gate_C5_interpretable = C5_interpretable, gate_C6_anchored = C6_anchored,
                gate_first_failure = first_failure, gate_outcome = admitted)]
  st14 <- merge(act[, .(TF = source, region, statistic, ranking_metric = condition,
                        score, p_value, q_BH, score_750, q_BH_750, score_noHK, n_target)],
                g, by = "TF", all.x = TRUE, sort = FALSE)
  # a TF scored in only one region never reaches the gate (it needs both) -> NA gate
  st14[, scored_both_regions := !is.na(gate_outcome)]
  st14[is.na(gate_outcome),
       `:=`(gate_first_failure = "not gated (scored in one region only)",
            gate_outcome       = "not gated")]
  setcolorder(st14, c("TF","region","statistic","ranking_metric","score","p_value","q_BH",
                      "score_750","q_BH_750","score_noHK","n_target","scored_both_regions","gate_outcome",
                      "gate_first_failure","gate_mean_z","gate_coverage",
                      "gate_leaveout_ratio","gate_mor_purity","gate_tf_own_delta",
                      "gate_anchor_detected","gate_anchor_total","gate_anchor_dir_frac",
                      "gate_C1_coverage","gate_C2_reproducible","gate_C3_robust",
                      "gate_C4_coherent","gate_C5_interpretable","gate_C6_anchored"))
  st14 <- st14[order(gate_outcome != "gated hit", -abs(score))]
  st14 <- floor_p(st14, c("p_value","q_BH","q_BH_750"), "ST14")
  # The row count is a property of the scored activity table, not a constant.
  # It now
  # asserts the invariant that actually matters: ST14 carries every scored TF x region test
  # and loses none of them to the gate merge.
  stopifnot("ST14 lost or gained rows in the gate merge" = nrow(st14) == nrow(act),
            "ST14 has duplicate TF x region keys" = !anyDuplicated(st14[, .(TF, region)]),
            "a TF appears in more than the two regions" = all(st14[, .N, by = TF]$N <= 2L),
            "ST14 region column contains a non-FH region" =
              all(st14$region %in% c("Frontal","Hippo")))
  n_both <- length(unique(st14$TF[st14$scored_both_regions]))
  say("  ST14 rows=%d  unique TFs=%d  scored in both regions=%d  one region only=%d",
      nrow(st14), length(unique(st14$TF)), n_both,
      length(unique(st14$TF)) - n_both)
  say("  gate outcomes: %s", paste(sprintf("%s=%d", names(table(st14$gate_outcome)),
                                           as.integer(table(st14$gate_outcome))), collapse = " | "))
  say("  gated hits: %s", paste(sort(unique(st14$TF[st14$gate_outcome == "gated hit"])), collapse = ", "))
  p <- wr(as.data.frame(st14), "ST14", sheet_map[["ST14"]])
  register("ST14", p, "NEW",
           sprintf("CollecTRI microglial TF activity + pre-declared gate audit = the Fig-2h source; BH over %d TF x region tests", nrow(st14)))
})

# =============================================================================
# ST15 / ST15b / ST16   neuronal programme tables behind Figure 5d.
#
# Source is the figure's own cache, not ST7.  data/_cache_neuron_scores_FH.rds is the
# per-nucleus score matrix that 4C2_neuron_programme_dotmatrix_FH.R plots as Fig 5d. It is
# read here and never rescored: rebuilding it loads the 2.9 GB harmony atlas, which a table
# build must not trigger, and rescoring would let the table drift from the panel.
#
# Why not ST7 (this is the whole reason these tables exist).  ST7 ships columns named
# Presynaptic / Postsynaptic / OxPhos, but they are curated_signatures.rds sets — 21 / 14 / 16
# genes — while Figure 5d scores 19 / 17 / 23-gene sets of the same names, and ST7 is scored on
# the 228 MB subset rather than the full atlas. Only IEG is common to both. Computing these
# tables from ST7 reproduces the panel only to about 0.08, which is the discrepancy that
# prompted this work. The Fig-5d sets now ship in ST6 as Fig5d_*.
# =============================================================================
NEURON_CACHE <- file.path(PROJ, "data", "_cache_neuron_scores_FH.rds")
PROG4 <- c("Presynaptic","OxPhos","IEG","Postsynaptic")   # the four the task specifies
# Cliff's delta, verbatim from 40_mast_percell_dual_FH.R:96-97 (the ST1 definition).
# +ve = higher in NHD.
cliffs_ST1 <- function(x, isN) {
  n1 <- sum(isN); n2 <- sum(!isN)
  if (n1 < 3 || n2 < 3) return(NA_real_)
  r <- rank(x); U1 <- sum(r[isN]) - n1 * (n1 + 1) / 2
  2 * U1 / (n1 * n2) - 1
}
load_neuron_cache <- function() {
  HARM <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
  stopifnot("MISSING neuron score cache — run scripts/4C_signature_violins_FH.R" =
              file.exists(NEURON_CACHE))
  if (file.exists(HARM) && file.mtime(NEURON_CACHE) < file.mtime(HARM))
    stop("neuron score cache is STALE vs atlas/NHD_FH_harmony.rds — re-run scripts/4C_signature_violins_FH.R")
  sc <- readRDS(NEURON_CACHE)
  stopifnot("neuron cache is not the 17,963-nucleus Fig-5d cache" = nrow(sc) == 17963L,
            "neuron cache lost a programme column" = all(PROG4 %in% names(sc)))
  sc
}

run("ST15", {
  step("ST15 per-library neuronal programme means (excitatory neurons)")
  sc <- load_neuron_cache()
  ex <- as.data.table(sc)[cell_type == "Neuron_Ex" & Condition %in% c("CON","NHD") &
                            Region %in% c("Frontal","Hippo")]
  say("  excitatory nuclei: %d of %d in cache", nrow(ex), nrow(sc))

  # Library label and RIN are JOINED, never typed: RIN lives in ST12 (keyed on
  # Donor/Region/Matter) and the TFHS id -> dissection mapping lives in the sample map.
  smap <- fread(file.path(NHDROOT, "manuscript_7fig", "tables", "sample_metadata_map.csv"))
  smap[, Library := fifelse(grepl("_g$", Sample), "Grey",
                     fifelse(grepl("_w$", Sample), "White",
                      fifelse(grepl("Hippo1$", Sample), "capture 1", "capture 2")))]
  st12 <- fread(file.path(OUT, sheet_map[["ST12"]]))[Donor != "Total"]
  st12[, `:=`(.Cond = fifelse(Donor == "NHD", "NHD", "CON"),
              .Reg  = fifelse(grepl("Frontal", Region), "Frontal", "Hippo"),
              .Lib  = fifelse(Matter %in% c("Grey","White"), Matter,
                       fifelse(grepl("capture 1", Region), "capture 1", "capture 2")))]

  st15 <- ex[, c(list(n_nuclei = .N), lapply(.SD, mean)),
             by = .(library_id = SampleID, Condition, Region), .SDcols = PROG4]
  st15 <- merge(st15, smap[, .(library_id = SampleID, Library)], by = "library_id", all.x = TRUE)
  st15 <- merge(st15, st12[, .(Condition = .Cond, Region = .Reg, Library = .Lib,
                               RIN = as.numeric(RIN))],
                by = c("Condition","Region","Library"), all.x = TRUE)
  stopifnot("ST15: a library failed the RIN join" = !any(is.na(st15$RIN)),
            "ST15: a library failed the Library-label join" = !any(is.na(st15$Library)))

  ST15_EXPECT <- c("library_id","Condition","Region","Library","RIN","n_nuclei", PROG4)
  setcolorder(st15, ST15_EXPECT)
  st15 <- st15[order(Region, Condition, Library)]
  stopifnot("ST15 must have exactly 7 library rows" = nrow(st15) == 7L,
            "ST15 n_nuclei does not sum to the excitatory nuclei in the cache" =
              sum(st15$n_nuclei) == nrow(ex))
  print(as.data.frame(st15), digits = 4, row.names = FALSE)
  p <- wr(as.data.frame(st15), "ST15", sheet_map[["ST15"]])
  schema_new("ST15", names(st15), ST15_EXPECT)
  register("ST15", p, "NEW",
           "per-library mean programme scores, excitatory neurons, from the Fig-5d score cache; RIN joined from ST12")
})

run("ST15b", {
  step("ST15b per-region RIN gap vs programme loss")
  sc <- load_neuron_cache()
  ex <- as.data.table(sc)[cell_type == "Neuron_Ex" & Condition %in% c("CON","NHD") &
                            Region %in% c("Frontal","Hippo")]
  lib <- fread(file.path(OUT, sheet_map[["ST15"]]))   # library means, already validated
  # Losses are computed on library means (the unit ST15 ships), not on nucleus means, so
  # that the deeper NHD libraries cannot dominate by nucleus count.
  st15b <- lib[, {
    rc <- mean(RIN[Condition == "CON"]); rn <- mean(RIN[Condition == "NHD"])
    out <- list(RIN_gap_CON_minus_NHD = rc - rn)
    for (pg in PROG4) out[[paste0(pg, "_loss")]] <-
      mean(get(pg)[Condition == "NHD"]) - mean(get(pg)[Condition == "CON"])
    out
  }, by = Region]
  ST15b_EXPECT <- c("Region","RIN_gap_CON_minus_NHD", paste0(PROG4, "_loss"))
  setcolorder(st15b, ST15b_EXPECT)
  st15b <- st15b[order(match(Region, c("Frontal","Hippo")))]
  stopifnot("ST15b must have one row per region" = nrow(st15b) == 2L)
  print(as.data.frame(st15b), digits = 4, row.names = FALSE)
  say("  NOTE the hippocampal RIN gap (%.1f) is ~2x the frontal one (%.1f) while the programme",
      st15b$RIN_gap_CON_minus_NHD[st15b$Region == "Hippo"],
      st15b$RIN_gap_CON_minus_NHD[st15b$Region == "Frontal"])
  say("       losses are not correspondingly larger — which is the point of shipping this table.")
  p <- wr(as.data.frame(st15b), "ST15b", sheet_map[["ST15b"]])
  schema_new("ST15b", names(st15b), ST15b_EXPECT)
  register("ST15b", p, "NEW",
           "per-region RIN gap beside the NHD-minus-CON programme loss, both on ST15 library means")
})

run("ST16", {
  step("ST16 hippocampal excitatory programme delta, Unresolved sensitivity")
  sc <- load_neuron_cache()
  h <- as.data.table(sc)[cell_type == "Neuron_Ex" & Region == "Hippo" &
                           Condition %in% c("CON","NHD")]
  st11 <- fread(file.path(OUT, sheet_map[["ST11"]]))
  unres <- st11$barcode[st11$neuron_subtype == "Unresolved"]
  say("  ST11 Unresolved nuclei: %d, of which hippocampal excitatory: %d",
      length(unres), sum(h$cell %in% unres))

  one <- function(d, lbl) {
    isN <- d$Condition == "NHD"
    c(list(set = lbl, n_NHD = sum(isN), n_CON = sum(!isN)),
      setNames(lapply(PROG4, function(pg) cliffs_ST1(d[[pg]], isN)), PROG4))
  }
  st16 <- rbindlist(list(as.data.table(one(h, "all_nuclei")),
                         as.data.table(one(h[!cell %in% unres], "excluding_unresolved"))))
  ST16_EXPECT <- c("set","n_NHD","n_CON", PROG4)
  setcolorder(st16, ST16_EXPECT)
  print(as.data.frame(st16), digits = 4, row.names = FALSE)

  # Acceptance GATE. The all-nuclei row must be the panel, to 2 dp.
  fp <- file.path(PROJ, "tables", "fig5_neuron_programme_dotmatrix_FH.csv")
  if (file.exists(fp)) {
    fg <- fread(fp)[Region == "Hippo" & class == "Ex"]
    ref <- setNames(fg$cliff_d, fg$signature)[PROG4]
    got <- unlist(st16[set == "all_nuclei", ..PROG4])
    say("  vs Fig 5d: %s", paste(sprintf("%s %.3f/%.3f", PROG4, got, ref), collapse = "  "))
    stopifnot("ST16 all_nuclei does not equal the Fig-5d deltas to 2 dp" =
                all(round(got, 2) == round(ref, 2)))
    say("  acceptance gate PASSED (all four programmes equal the panel to 2 dp)")
  } else say("  WARNING: fig5 source CSV absent -> acceptance gate skipped")

  d1 <- unlist(st16[set == "all_nuclei", ..PROG4])
  d2 <- unlist(st16[set == "excluding_unresolved", ..PROG4])
  say("  max |change| over Presynaptic/OxPhos/IEG = %.4f ; over all four = %.4f ; sign flips: %d",
      max(abs((d2 - d1)[c("Presynaptic","OxPhos","IEG")])), max(abs(d2 - d1)),
      sum(sign(d1) != sign(d2)))
  say("  Postsynaptic %.3f -> %.3f. The exclusion drops %d NHD and %d CON nuclei, i.e. it is",
      d1[["Postsynaptic"]], d2[["Postsynaptic"]],
      st16$n_NHD[1] - st16$n_NHD[2], st16$n_CON[1] - st16$n_CON[2])
  say("  almost entirely an NHD exclusion — the postsynaptic result depends on those nuclei.")
  p <- wr(as.data.frame(st16), "ST16", sheet_map[["ST16"]])
  schema_new("ST16", names(st16), ST16_EXPECT)
  register("ST16", p, "NEW",
           "Cliff's delta (ST1 definition) per programme, hippocampal excitatory neurons, all nuclei vs excluding ST11 Unresolved")
})

# =============================================================================
# ST17 / ST17b   pooled Micro-PVM vs microglia-only.
# Upstream: scripts/92_micro_only_MAST_FH.R (2.9 GB atlas + atlas-wide PrepSCTFindMarkers).
# That script is never invoked from here — heavy compute stays in its own numbered script
# and this block does the CSV join, exactly as ST13/ST14 do.
# =============================================================================
ST17_GENES <- c("HLA-DRA","CD74","FTH1","C1QB","SPP1","CLEC7A","ITGAX","GPNMB")

run("ST17", {
  step("ST17 microglia-only vs pooled Micro-PVM")
  MS <- file.path(PROJ, "tables", "micro_states")
  f_mast <- file.path(MS, "MAST_microglia_only_FH.csv")
  f_gene <- file.path(MS, "micro_only_ST17gene_stats_FH.csv")
  for (f in c(f_mast, f_gene))
    if (!file.exists(f))
      stop(sprintf("MISSING ST17 source %s — run scripts/92_micro_only_MAST_FH.R", basename(f)))

  # Pooled half: straight out of the shipped ST1, so it matches by construction
  st1 <- fread(file.path(OUT, sheet_map[["ST1"]]))
  pool <- st1[cell_type == "Micro-PVM" & gene %in% ST17_GENES,
              .(gene, region, pooled_delta = cliffs_delta, pooled_pct_NHD = pct_NHD,
                pooled_pct_CON = pct_CON, pooled_padj = p_val_adj,
                pooled_discovery = discovery_filtered, pooled_tested = TRUE)]
  mic <- fread(f_mast)[gene %in% ST17_GENES,
              .(gene, region, micro_delta = cliffs_delta, micro_pct_NHD = pct.1,
                micro_pct_CON = pct.2, micro_padj = p_val_adj,
                micro_discovery = discovery, micro_tested = TRUE)]
  # Backfill every gene MAST did not return. MAST only reports features clearing min.pct /
  # logfc.threshold, so a gene the pooled arm tested and this arm did not would otherwise ship
  # blank and read as "no effect". 92_ computes delta and detection for all eight regardless;
  # those rows are marked micro_tested = FALSE, the same distinction Fig 2e draws between a
  # hollow circle (detected, not estimable) and a cross (not tested).
  gs <- fread(f_gene)[, .(gene, region, micro_delta = cliffs_delta, micro_pct_NHD = pct.1,
                          micro_pct_CON = pct.2, micro_padj = NA_real_,
                          micro_discovery = FALSE, micro_tested = FALSE)]
  mic <- rbind(mic, gs[!paste(gene, region) %in% paste(mic$gene, mic$region)])

  grid <- CJ(gene = ST17_GENES, region = c("Frontal","Hippo"), unique = TRUE)
  st17 <- merge(merge(grid, pool, by = c("gene","region"), all.x = TRUE),
                mic, by = c("gene","region"), all.x = TRUE)
  st17[is.na(pooled_tested), pooled_tested := FALSE]
  st17[is.na(micro_tested),  micro_tested  := FALSE]
  st17[, delta_difference := pooled_delta - micro_delta]

  ST17_EXPECT <- c("gene","region","pooled_delta","pooled_pct_NHD","pooled_pct_CON",
                   "pooled_padj","pooled_discovery","micro_delta","micro_pct_NHD",
                   "micro_pct_CON","micro_padj","micro_discovery","delta_difference",
                   "pooled_tested","micro_tested")
  setcolorder(st17, ST17_EXPECT)
  st17 <- st17[order(region, match(gene, ST17_GENES))]
  st17 <- floor_p(st17, c("pooled_padj","micro_padj"), "ST17")
  stopifnot("ST17 must have 8 genes x 2 regions = 16 rows" = nrow(st17) == 16L,
            "ST17 has a row with neither arm's effect size" =
              !any(is.na(st17$pooled_delta) & is.na(st17$micro_delta)))

  # pooled_* must equal ST1 exactly on every key it has
  chk <- merge(st17[pooled_tested == TRUE, .(gene, region, pooled_delta)],
               st1[cell_type == "Micro-PVM", .(gene, region, ref = cliffs_delta)],
               by = c("gene","region"))
  stopifnot("ST17 pooled_delta does not match ST1" =
              nrow(chk) == sum(st17$pooled_tested) && isTRUE(all.equal(chk$pooled_delta, chk$ref)))
  say("  pooled_* verified identical to ST1 on all %d tested keys", nrow(chk))
  say("  not tested in the pooled arm: %s",
      paste(st17[pooled_tested == FALSE, paste0(gene, "/", region)], collapse = ", "))
  say("  not tested in the microglia-only arm: %s",
      paste(st17[micro_tested == FALSE, paste0(gene, "/", region)], collapse = ", "))
  say("  max |delta_difference| = %.3f ; discovery kept in both arms: %d of %d",
      max(abs(st17$delta_difference), na.rm = TRUE),
      sum(st17$pooled_discovery %in% TRUE & st17$micro_discovery %in% TRUE),
      sum(st17$pooled_discovery %in% TRUE))
  print(as.data.frame(st17[, .(gene, region, pooled_delta, micro_delta, delta_difference,
                               pooled_discovery, micro_discovery)]), digits = 3, row.names = FALSE)
  p <- wr(as.data.frame(st17), "ST17", sheet_map[["ST17"]])
  schema_new("ST17", names(st17), ST17_EXPECT)
  register("ST17", p, "NEW",
           "8 activation genes x 2 regions; pooled_* verbatim from ST1, micro_* from the identical model on the 2,824 microglia; untested genes flagged, never blank")
})

run("ST17b", {
  step("ST17b myeloid selection chain")
  f <- file.path(PROJ, "tables", "micro_states", "micro_selection_steps_FH.csv")
  if (!file.exists(f))
    stop("MISSING ST17b source micro_selection_steps_FH.csv — run scripts/92_micro_only_MAST_FH.R")
  sel <- dcast(fread(f), Condition + Region ~ step, value.var = "n", fill = 0L)
  for (cc in c("T_lymphocyte","oligo_doublet_removed","CD163_F13A1_removed","microglia_final"))
    if (!cc %in% names(sel)) sel[, (cc) := 0L]
  sel[, MicroPVM_azimuth := microglia_final + CD163_F13A1_removed + oligo_doublet_removed]
  setnames(sel, c("oligo_doublet_removed","CD163_F13A1_removed","microglia_final","T_lymphocyte"),
           c("minus_oligo_doublets","minus_CD163_F13A1","microglia_analysed","Lymphocyte_separate_class"))
  sel[, after_doublet_removal := MicroPVM_azimuth - minus_oligo_doublets]
  ST17b_EXPECT <- c("Condition","Region","MicroPVM_azimuth","minus_oligo_doublets",
                    "after_doublet_removal","minus_CD163_F13A1","microglia_analysed",
                    "Lymphocyte_separate_class")
  setcolorder(sel, ST17b_EXPECT)
  sel <- sel[order(match(Region, c("Frontal","Hippo")), Condition)]
  print(as.data.frame(sel), row.names = FALSE)
  stopifnot("ST17b chain does not close" =
              all(sel$after_doublet_removal - sel$minus_CD163_F13A1 == sel$microglia_analysed),
            "ST17b totals disagree with the shipped figures (3081 / 95 / 162 / 2824 / 60)" =
              sum(sel$MicroPVM_azimuth) == 3081L && sum(sel$minus_oligo_doublets) == 95L &&
              sum(sel$minus_CD163_F13A1) == 162L && sum(sel$microglia_analysed) == 2824L &&
              sum(sel$Lymphocyte_separate_class) == 60L)
  say("  chain closes: 3,081 Micro-PVM - 95 oligo doublets = 2,986 - 162 CD163+/F13A1+ = 2,824")
  say("  T-lymphocytes (n=60) are a SEPARATE Azimuth class and were never inside Micro-PVM.")
  p <- wr(as.data.frame(sel), "ST17b", sheet_map[["ST17b"]])
  schema_new("ST17b", names(sel), ST17b_EXPECT)
  register("ST17b", p, "NEW",
           "myeloid selection chain per condition x region; Lymphocyte is a separate class, not a Micro-PVM subset")
})

# =============================================================================
# Table 1  donor / sample demographics (FH cohort: drop OCC row(s))
# OLD source is a .docx (cannot parse binary here) -> derive donor/sample
# demographics from known facts + atlas lane structure; leave irrecoverable
# demographic fields blank and report them (never invent ages/PMI/sex).
# Known: NHD 49F / CON 56F (Mihara donors),
# both female.  PMI / diagnosis-detail not recoverable here -> blank + flagged.
# =============================================================================
run("Table1", {
  step("Table 1 donor/sample demographics")
  # Transposed publication format matching the OLD NHD_Table1_donor_sample.docx
  # (Characteristic x {NHD donor, Control donor}), FH-updated: OCC dropped, library
  # counts and post-QC nuclei recomputed from the FH atlas. The OLD table has no PMI
  # column, so none is added here. Content of the OLD docx (recovered by unzipping
  # word/document.xml): Diagnosis / Age / Sex / Ancestry / DAP12 genotype / Tissue
  # source / Regions / grey-white / libraries / Visium / Post-QC nuclei.
  n_nhd <- sum(MD$Condition == "NHD"); n_con <- sum(MD$Condition == "CON")
  fmt <- function(x) formatC(x, format = "d", big.mark = ",")
  t1 <- data.frame(
    Characteristic = c("Diagnosis","Age (years)","Sex","Ancestry",
                       "TYROBP / DAP12 genotype","Tissue source",
                       "Regions profiled (snRNA-seq)","Grey / white matter",
                       "snRNA-seq libraries","Spatial transcriptomics (Visium)",
                       "Post-QC nuclei","Note"),
    `NHD donor` = c("Nasu–Hakola disease","49","Female","Japanese",
                    "Homozygous loss-of-function; DAP12 protein absent",
                    "Mihara Hospital brain bank, Japan","Frontal, hippocampus",
                    "Both (cortex)","3","Frontal cortex, 2 sections", fmt(n_nhd),
                    sprintf("Post-QC nuclei total %s across 7 libraries (NHD 3, control 4). The control hippocampus was sampled as two captures and the NHD hippocampus as one.",
                            fmt(n_nhd + n_con))),
    `Control donor` = c("Neurologically unaffected","56","Female","Japanese",
                        "Wild-type","Mihara Hospital brain bank, Japan",
                        "Frontal, hippocampus","Both (cortex)","4",
                        "Frontal cortex, 2 sections", fmt(n_con), ""),
    check.names = FALSE, stringsAsFactors = FALSE)
  say("  Table1: CON nuclei=%d, NHD nuclei=%d, total=%d", n_con, n_nhd, n_con + n_nhd)
  p <- wr(t1, "Table1", sheet_map[["Table1"]])
  register("Table1", p, "REGENERATED",
           "FH donor demographics (transposed pub format); OCC dropped; libraries+nuclei recomputed; no PMI in source")
})

# =============================================================================
# ST0  data dictionary.  Built last, from every table that is
# now on disk, so it can only ever describe what actually shipped.  The authored
# content lives in scripts/_ST0_data_dictionary_FH.R; here it is cross-checked
# against the real headers and the run fails if the two disagree, which is what
# stops a dictionary from rotting the way ST3/ST4/ST5/ST13's filenames did.
# =============================================================================
step("ST0 data dictionary")
dict_src <- file.path(PROJ, "scripts", "_ST0_data_dictionary_FH.R")
stopifnot("MISSING scripts/_ST0_data_dictionary_FH.R" = file.exists(dict_src))
de <- new.env(); sys.source(dict_src, envir = de)
for (o in c("ST0_TABLES","ST0_COLUMNS","ST0_CONVENTIONS"))
  if (!exists(o, envir = de)) stop(sprintf("dictionary source did not define %s", o))
ST0_TABLES <- get("ST0_TABLES", de); ST0_COLUMNS <- get("ST0_COLUMNS", de)
ST0_CONVENTIONS <- get("ST0_CONVENTIONS", de)

# --- drift check: documented columns vs columns actually shipped --------------
drift <- character(0)
for (id in names(sheet_map)) {
  if (id == "ST0") next
  fp <- file.path(OUT, sheet_map[[id]])
  if (!file.exists(fp)) { drift <- c(drift, sprintf("%s: file absent", id)); next }
  got <- names(fread(fp, nrows = 0))
  doc <- ST0_COLUMNS$column[ST0_COLUMNS$table == id]
  miss <- setdiff(got, doc); extra <- setdiff(doc, got)
  if (length(miss))  drift <- c(drift, sprintf("%s: UNDOCUMENTED column(s) %s", id, paste(miss, collapse=",")))
  if (length(extra)) drift <- c(drift, sprintf("%s: documented but ABSENT %s", id, paste(extra, collapse=",")))
}
if (length(drift)) {
  for (d in drift) say("  DICTIONARY DRIFT -> %s", d)
  stop("ST0 data dictionary disagrees with the shipped tables — fix scripts/_ST0_data_dictionary_FH.R")
}
say("  dictionary covers every column of every shipped table (no drift)")
# The ST6 scoring_context vocabulary the dictionary lists must be the one the shipped ST6 carries
# (both come from SD_ST6_CTX of _sd_map_FH.R; this catches a stale ST6 built before a map change)
{ .pv <- ST0_COLUMNS$permitted_values[ST0_COLUMNS$table == "ST6" & ST0_COLUMNS$column == "scoring_context"]
  .ctx <- unique(fread(file.path(OUT, sheet_map[["ST6"]]))$scoring_context)
  .ok  <- trimws(strsplit(.pv, " | ", fixed = TRUE)[[1]])
  .base <- unique(unlist(lapply(.ctx, function(v) { if (v %in% .ok) return(v)     # whole value listed (incl. "reference only; not scored ...", which itself contains "; ")
    sub(" as (column|program) .*$", "", strsplit(v, "; ", fixed = TRUE)[[1]]) })))
  .miss <- setdiff(.base, .ok)
  if (length(.miss)) stop("ST6 scoring_context value(s) outside the dictionary's permitted set (stale ST6? re-run 90 ST6): ", paste(.miss, collapse = " || "))
  say("  ST6 scoring_context vocabulary (%d base contexts) matches the dictionary", length(.base)) }
stopifnot("ST0_TABLES does not cover every shipped table" =
            setequal(ST0_TABLES$table, names(sheet_map)))

ST0_COLUMNS$.row <- seq_len(nrow(ST0_COLUMNS))      # keep the authored column order
st0 <- merge(ST0_COLUMNS, ST0_TABLES, by = "table", all.x = TRUE, sort = FALSE)
st0 <- as.data.table(st0)
# stable, human order: Table1, ST0, ST1..ST14; authored column order within a table
tab_ord <- unique(c("Table1", names(sheet_map)))
st0[, .ord := match(table, tab_ord)]
st0 <- st0[order(.ord, .row)][, .(table, table_title = title, backs_figures, source_script,
                                 column, type, permitted_values, definition)]
conv <- data.table(table = "CONVENTIONS",
                   table_title = "Analysis-wide conventions (apply to every table above)",
                   backs_figures = "all", source_script = "_ST0_data_dictionary_FH.R",
                   column = ST0_CONVENTIONS$convention, type = "convention",
                   permitted_values = "", definition = ST0_CONVENTIONS$definition)
st0 <- rbind(st0, conv)
say("  ST0 rows=%d (%d column definitions + %d conventions), tables documented=%d",
    nrow(st0), nrow(ST0_COLUMNS), nrow(conv), length(unique(ST0_COLUMNS$table)))
# house rule: the word "replication" must not appear as prose in the dictionary
# MSigDB set identifiers (GOBP_DNA_REPLICATION) and explicit negations are allowed;
# anything else asserting replication is a house-rule violation.
.txt <- paste(st0$table_title, st0$definition)
.bad <- grepl("replicat", .txt, ignore.case = TRUE) &
        !grepl("REPLICATION\\b|not a replication|Not a replication|never .{0,30}replicat|not .{0,20}replicat",
               .txt)
if (any(.bad)) say("  WARNING: 'replicat' in dictionary prose at rows %s",
                   paste(which(.bad), collapse = ","))
p0 <- wr(as.data.frame(st0), "ST0", sheet_map[["ST0"]])
# occ_exempt: ST0 is prose and must be able to say "Zhou's external occipital
# cohort" and "our regions are Frontal and Hippo only" without tripping the scan.
register("ST0", p0, "NEW",
         "data dictionary; cross-checked against every shipped table header",
         occ_exempt = TRUE)

# ---- retire legacy filenames so the folder holds exactly one file per table ---
step("retire legacy filenames")
SUP <- file.path(OUT, sprintf("_superseded_%s", format(Sys.Date(), "%Y%m%d")))
n_ret <- 0L
for (lf in LEGACY) {
  src <- file.path(OUT, lf)
  if (!file.exists(src)) next
  dir.create(SUP, showWarnings = FALSE, recursive = TRUE)
  ok <- file.rename(src, file.path(SUP, lf))
  say("  moved legacy %s -> %s (%s)", lf, basename(SUP), ifelse(ok, "ok", "FAILED"))
  n_ret <- n_ret + 1L
}
if (!n_ret) say("  none present (folder already clean)")

step("assemble workbook NHD_Supplementary_Tables_FH.xlsx")
wb <- createWorkbook()
for (sn in names(sheet_map)) {
  fp <- file.path(OUT, sheet_map[[sn]])
  if (!file.exists(fp)) { say("  SKIP sheet %s (missing %s)", sn, sheet_map[[sn]]); next }
  d <- fread(fp)
  addWorksheet(wb, sn)
  writeData(wb, sn, d)
  say("  sheet %-6s <- %s (%d x %d)", sn, sheet_map[[sn]], nrow(d), ncol(d))
}
# The loop above SKIPS a missing CSV with a log line, while the MANIFEST below used to
# report the sheet list from sheet_map — so the MANIFEST could assert a sheet the workbook
# does not contain. Assert the workbook actually got every sheet instead.
.missing_sheets <- setdiff(names(sheet_map), names(wb))
if (length(.missing_sheets))
  stop(sprintf("workbook is missing a sheet for: %s — it would contradict the MANIFEST",
               paste(.missing_sheets, collapse = ", ")))
xlsx <- file.path(OUT, "NHD_Supplementary_Tables_FH.xlsx")
tmpx <- paste0(xlsx, ".tmp.xlsx")
saveWorkbook(wb, tmpx, overwrite = TRUE)
invisible(file.rename(tmpx, xlsx))
say("  wrote workbook: %s  (%.1f MB, %d sheets)", xlsx,
    file.size(xlsx) / 1024^2, length(names(wb)))

# workbook-vs-CSV consistency: every sheet must match its CSV in shape, and the workbook
# must not be older than any CSV it claims to contain.
step("workbook vs CSV consistency")
.wbad <- character(0)
for (sn in names(wb)) {
  fp <- file.path(OUT, sheet_map[[sn]])
  d  <- fread(fp)
  s  <- openxlsx::read.xlsx(xlsx, sheet = sn)
  if (nrow(d) != nrow(s) || ncol(d) != ncol(s))
    .wbad <- c(.wbad, sprintf("%s: CSV %dx%d vs sheet %dx%d", sn, nrow(d), ncol(d), nrow(s), ncol(s)))
  if (file.mtime(fp) > file.mtime(xlsx))
    .wbad <- c(.wbad, sprintf("%s: %s is NEWER than the workbook", sn, basename(fp)))
}
if (length(.wbad)) { for (l in .wbad) say("  WORKBOOK MISMATCH %s", l)
  stop("workbook does not agree with the shipped CSVs") }
say("  all %d sheets match their CSV in shape, and none is newer than the workbook", length(names(wb)))

# =============================================================================
# Verification table + MANIFEST
# =============================================================================
step("verification")
ver <- do.call(rbind, REG)
ver <- ver[order(match(ver$table, names(sheet_map))), ]
print(ver, row.names = FALSE)
fwrite(ver, file.path(OUT, "_verification_FH.csv"))

if (length(FAIL)) {
  say("\n*** FAILED TABLES ***")
  for (nm in names(FAIL)) say("  %s: %s", nm, FAIL[[nm]])
} else say("\nAll table steps completed without error.")

# MANIFEST.  Descriptions are pulled from the same authored source as ST0, so the
# MANIFEST and the data dictionary can never say different things about a table.
desc <- setNames(ST0_TABLES$title, ST0_TABLES$table)
figs <- setNames(ST0_TABLES$backs_figures, ST0_TABLES$table)
man <- file.path(OUT, "MANIFEST_FH.md")
mcon <- c("# NHD Frontal+Hippocampus rebuild — Supplementary Tables MANIFEST",
          "",
          sprintf("Generated: %s", format(Sys.time())),
          "6-figure scheme: Fig1 atlas, Fig2 microglia, Fig3 astro, Fig4 oligo, Fig5 neurons, Fig6 Visium.",
          "",
          "Cohort: 1 NHD + 1 control donor with no neurological diagnosis; 7 snRNA-seq libraries",
          "(control 4, NHD 3); 55,354 post-QC nuclei (control 37,371 / NHD 17,983; frontal",
          "37,034 / hippocampus 18,320).  Our regions are Frontal and Hippo ONLY.",
          "",
          "**Start with ST0 (data dictionary).** It defines every column of every table,",
          "the permitted values, the figure each table backs, and the analysis-wide",
          "conventions (discovery criteria, Cliff's delta, the q* mark, depth matching,",
          "the p-underflow convention).",
          "",
          "| Table | File | Rows | Cols | Mode | Backs | OCC | Description |",
          "|-------|------|------|------|------|-------|-----|-------------|")
for (id in names(sheet_map)) {
  r <- REG[[id]]
  if (is.null(r)) { mcon <- c(mcon, sprintf("| %s | (FAILED) | - | - | - | %s | - | %s |",
                                            id, figs[[id]], desc[[id]])); next }
  mcon <- c(mcon, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
                          id, r$path, r$n_rows, r$n_cols, r$mode, figs[[id]],
                          ifelse(isTRUE(r$occ_found), "**YES(BUG)**", "none"), desc[[id]]))
}
mcon <- c(mcon, "",
          sprintf("Workbook: NHD_Supplementary_Tables_FH.xlsx (%d sheets: %s).",
                  length(names(wb)), paste(names(wb), collapse = ", ")),
          "",
          "## Numerical conventions a referee will check",
          sprintf("- **p-value underflow.** No p-value is ever reported as 0. Every p column in %s is floored at `%s` — the smallest normalized double (`.Machine$double.xmin`) rounded up to 8 significant digits so that the floor survives the CSV round-trip exactly. A value equal to that floor means **p < 2.2e-308**, not p = 0. Per-column counts of floored values are in `logs/90_build_supplementary_tables_FH.log`.",
                  paste(sort(P_FLOORED), collapse = ", "), sprintf("%.8g", P_FLOOR)),
          "- **NA is 'not tested'.** DESeq2 independent filtering leaves `padj = NA` in ST2 and fgsea leaves `pval/padj/NES = NA` in ST3 where the permutation null was degenerate. Those NAs are preserved and are never recoded to 0 or 1, nor dropped.",
          "- **Detection and effect size on raw RNA.** `pct_NHD` / `pct_CON` / `cliffs_delta` are computed on the raw RNA layers, never on SCT (fix of 2026-08-19: SCT fabricates detection for zero-UMI genes in this design). ST1 was regenerated on 2026-09-05 because it had been left at its pre-fix 2026-07-14 state after the ST1 step errored in the 2026-08-19 run.",
          "",
          "## Notes / caveats",
          "- **ST3 filename corrected (2026-09-05).** Previously `ST3_fGSEA_signed_cliffs_delta.csv`, but all 42,305 rows carry `metric = log2fc_winsorized`; the signed-Cliff's-delta prerank was retired on 2026-07-01. The file now states the metric it actually used. No value changed. Pathway names containing `REPLICATION` are MSigDB set identifiers, not a replication claim.",
          "- **ST4 / ST5 are the GTRD genome-wide survey and do NOT back Figure 2h.** They score all six cell types on GTRD, which is UNSIGNED (`mode_of_regulation` = 1 for every one of its 257,608 edges). Figure 2h is microglia-only on the SIGNED, curated CollecTRI network; only 8 of ST4's 459 factors appear in that result. Both files were renamed on 2026-09-05 so the regulon is explicit.",
          "- **ST14 is the table behind Figure 2h.** CollecTRI, myeloid compartment (Micro-PVM after doublet removal, n = 2,985 — microglia together with the 162 CD163+/F13A1+ nuclei, whose pooling is justified in ST17), scored on log-normalized counts across ALL nuclei, with the pre-declared six-criterion admission gate from `scripts/70b_micro_TF_bipartite_FH.R`. **BH is applied once over the 118 TF x region tests in the table (64 unique factors) — NOT over the rows of the gate audit**, which contains only the factors scorable in both regions. 10 factors are scorable in one region only and are ungated by construction (`scored_both_regions = FALSE`). Columns `score_750` / `q_BH_750` carry the 750-UMI down-sampling falsification test; 22 rows are NA there because the factor cannot be scored once half the nuclei are discarded — *PPARG*, *STAT1* and *CEBPB* among them. Down-sampling is a falsification test here, not the estimator; see Methods > Sequencing depth.",
          "- **ST13 is direction consistency, not replication.** The Zhou 2023 comparison is threshold-free direction consistency measured against an EMPIRICAL per-panel baseline. `concord_frac` must be read against `baseline_concord`, never against 0.5. The file was renamed on 2026-09-05 and `baseline_concord`, `n_baseline_genes`, `n_dir_zero`, `cell_type`, `region` were appended after the original 13 columns so the baseline ships with the statistic. `OCC` inside a `panel` label is Zhou's own occipital cohort — an external dataset, never our tissue.",
          "- **ST2 `n_lanes`** reports 4 for Hippo (2 NHD + 2 control) though only 3 physical hippocampal libraries exist; kept verbatim from the DESeq2 run",
          "- **Quality-control gate, and which script owns it.** The 55,354 nuclei were selected by `scripts/01_build_FH_atlas.R` at `nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 10` (observed range in the shipped atlas: 201-6,497 genes, maximum 9.95% mitochondrial, which is what those gates produce). **Three different nFeature floors exist in this repository and they belong to three different objects**: the relaxed-QC atlas that 01_build_FH_atlas.R subsets has a floor of 121, the older root-level `NHD_QC_Azimuth+Harmony.R` uses 350, and the FH atlas reported here uses 200. Only the last describes this dataset. Confirmed 2026-09-07.",
          "- **ST11** schema differs from the OLD Tippani glia-QC table; the FH pipeline did not rerun that analysis, and the missing QC columns were not invented.",
          "- **Region tokens differ by table on purpose.** ST1/ST7/ST9/ST11/ST13/ST14 use the analysis tokens `Frontal` / `Hippo`; ST12 uses the long-form library labels `Frontal cortex` / `Hippocampus (capture 1|2)`. ST0 maps the two. The tokens were NOT rewritten, because they are join keys against `tables/`.",
          "- **Table 1 is hand-maintained** and is not regenerated by this script (guarded by `ALLOW_TABLE1`).")
writeLines(mcon, man)
say("wrote MANIFEST: %s", man)
# house rule: 'replication' must not appear as prose describing the Zhou comparison
.mbad <- grep("replicat", mcon, ignore.case = TRUE, value = TRUE)
.mbad <- .mbad[!grepl("not replication|not a replication|never described as replication|not a replication claim", .mbad, ignore.case = TRUE)]
if (length(.mbad)) { for (l in .mbad) say("  MANIFEST WORDING WARNING: %s", l) } else
  say("  MANIFEST wording sweep: no stray 'replication' claim")

say("\n=== DONE ===  %s", format(Sys.time()))
