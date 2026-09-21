#!/usr/bin/env Rscript
# =============================================================================
# 100_package_supplementary_data_FH.R — Package the supplementary tables for the journal as "Supplementary Data 1-11": eleven themed multi-sheet workbooks.
# -----------------------------------------------------------------------------
# Why. The journal does not accept large or Excel tables as
# "Supplementary Table": they must be labeled "Supplementary Data N", uploaded as
# separate files, and given a title + legend inside the single Supplementary
# Information PDF.
#
# The map is fixed and lives in scripts/_sd_map_FH.R (SD_WB = workbook number ->
# title; SD_SHEETS = one row per data sheet: former table id -> workbook, sheet name,
# pre-consolidation number). Every citation of another file inside a README or the
# dictionary is rendered by sd_ref() from that map, never typed; the old numbers
# survive only in _build/SD_rename_map.csv (the re-keying aid for the prose).
#
# What it does.
#   1. Reads every shipped table from supplementary_tables/ (the 90_build outputs —
#      never re-derives a number) plus the six depth-diagnostic tables from tables/,
#      running every per-table guard the single-file package had.
#   2. Writes manuscript/Final_fig_and_tables/Supplementary_Data/Supplementary_Data_N.xlsx
#      (N = 1-11), each opening with a README sheet: label, title, sheet list,
#      conventions, then one block per former table (its sheet name(s), what it holds,
#      figures backed, source script, column dictionary) followed by the data sheets.
#   3. Applies sd_rekey_old() to every text cell of every sheet (the packaging-time
#      safety net for a pre-consolidation number typed somewhere upstream) and then
#      asserts that no "Supplementary Data 12-27" survives in any README or in the
#      dictionary.
#   4. Writes to Supplementary_Data/_build/: Supplementary_Data_INDEX.csv (one row per
#      sheet), MANIFEST.md, SD_rename_map.csv (every old label -> new workbook + sheet)
#      and FIND_REPLACE_for_docx.md (a two-pass, collision-free re-keying recipe).
#   5. Removes the superseded Supplementary_Data_12-27.xlsx from the delivery folder
#
# Run:  Rscript scripts/100_package_supplementary_data_FH.R
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(openxlsx); library(tools) })

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
say <- function(...) cat(sprintf(...), "\n")

SRC   <- file.path(PROJ, "supplementary_tables")
TAB   <- file.path(PROJ, "tables")
# Visium domain vocabulary for the ST20/ST21/ST22 dictionaries comes from integrated_domain_levels.csv
# (DOM_LEV / DOM_GREY / DOM_WHITE / DOM_INFO), so a new domain (e.g. an L1 rim) needs no edit here
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))
.list_and <- function(x) if (length(x) <= 1) x else paste0(paste(head(x, -1), collapse = ", "), " and ", tail(x, 1))
DOM_OTHER <- DOM_INFO$domain[DOM_INFO$compartment == "other"]
OUT   <- file.path(PROJ, "manuscript", "Final_fig_and_tables", "Supplementary_Data")
BUILD <- file.path(OUT, "_build")   # every builder aid goes here
dir.create(BUILD, showWarnings = FALSE, recursive = TRUE)
# builder aids written by earlier runs into Supplementary_Data/ itself are moved (not deleted) into _build/
for (.f in c("SD_rename_map.csv", "FIND_REPLACE_for_docx.md", "MANIFEST.md", "Supplementary_Data_INDEX.csv"))
  if (file.exists(file.path(OUT, .f))) { file.rename(file.path(OUT, .f), file.path(BUILD, .f)); say("moved stale %s -> _build/", .f) }

# ---- the fixed map (scripts/_sd_map_FH.R) -----------------------------------
source(file.path(PROJ, "scripts", "_sd_map_FH.R"))   # SD_WB, SD_SHEETS, SD_DROPPED, SD_OLD_ORDER, SD_TITLE_PUBLIC, sd_ref(), sd_rekey_old(), sd_rename_map()
N_WB     <- nrow(SD_WB)
ST_ORDER <- SD_OLD_ORDER                          # the 27 former tables (pre-consolidation order = their old number)
ST_SHIP  <- unique(SD_SHEETS$st)                  # the 25 that are placed in a workbook
stopifnot(N_WB == 11L, length(ST_ORDER) == 27L, length(ST_SHIP) == 25L, setequal(setdiff(ST_ORDER, ST_SHIP), c("ST4", "ST5")),
          setequal(SD_DROPPED$st, c("ST4", "ST5")))
TITLE_PUBLIC <- SD_TITLE_PUBLIC                   # referee-facing title of every former table (the "Holds" line of its README block)
stopifnot(setequal(names(TITLE_PUBLIC), ST_ORDER))
sd_wb_of  <- function(id) sd_workbook_of(id)                                    # former table id -> workbook number
sd_label  <- function(id) sprintf("Supplementary Data %d", sd_wb_of(id))         # ... -> "Supplementary Data N"
sd_file_n <- function(n)  sprintf("Supplementary_Data_%d.xlsx", n)
old_num   <- function(id) match(id, SD_OLD_ORDER)                                # pre-consolidation number (1-27)

# shipped CSVs (one per former table, as 90_build writes them)
st_files <- list.files(SRC, pattern = "^ST[0-9]+[bc]?_.*\\.csv$", full.names = TRUE)
st_files <- st_files[!grepl("^ST19b_", basename(st_files))]   # ST19b ships as sheet Within_Zhou_donors beside ST19
st_files <- st_files[!grepl("^ST22b_", basename(st_files))]   # ST22b ships as sheet Layer_IV_order beside ST22
st_files <- st_files[!grepl("^ST23[bc]_", basename(st_files))]   # ST23b / ST23c ship as Genes_by_band / Gene_sets_Visium beside ST23
st_files <- st_files[!grepl("^ST24b_", basename(st_files))]   # ST24b ships as sheet By_program beside ST24
st_id    <- sub("_.*$", "", basename(st_files))
stopifnot("duplicate ST ids in supplementary_tables/ (stray copy?)" = !anyDuplicated(st_id))
names(st_files) <- st_id
missing <- setdiff(setdiff(ST_SHIP, "ST18"), st_id)   # ST19-ST24 are shipped CSVs like the others (ST20 written by script 111, ST21 by 114, ST22 by 125, ST23 by 126, ST24 by 127)
if (length(missing)) stop("Missing shipped tables in supplementary_tables/: ", paste(missing, collapse = ", "))
# the two tables that left the package must still exist: 103 ships them as reference tables and 90 documents them
stopifnot("ST4 / ST5 csv missing from supplementary_tables/ (103 ships them in the code release)" = all(c("ST4", "ST5") %in% st_id))

# depth diagnostics -> ST18 (six sheets; schemas differ, so they stay separate)
depth_files <- list.files(TAB, pattern = "^depth_[A-F]_.*_FH\\.csv$", full.names = TRUE)
stopifnot("expected six depth_A..F tables in tables/" = length(depth_files) == 6)
depth_sheet <- sub("_FH\\.csv$", "", sub("^depth_", "", basename(depth_files)))  # A_relative_RNA_content ...
stopifnot("depth tables must be exactly one each of A-F" = setequal(substr(depth_sheet, 1, 1), LETTERS[1:6]),
          "depth sheet names must be the ones the fixed map lists for workbook 10" = identical(depth_sheet, sd_sheets_of("ST18")))

# ---- dictionary (ST0) for README sheets -----------------------------------
dict <- fread(st_files[["ST0"]], colClasses = "character")
titles <- unique(dict[table != "CONVENTIONS", .(table, table_title, backs_figures, source_script)])
title_of <- function(id) {
  r <- titles[table == id]
  if (nrow(r)) r$table_title[1] else NA_character_
}
# every shipped table must be documented in ST0, or its README would ship blank
POST_ST0 <- c("ST18","ST19","ST20","ST21","ST22","ST23","ST24")   # built after the dictionary; documented in this script
.undoc <- setdiff(ST_ORDER, POST_ST0)[vapply(setdiff(ST_ORDER, POST_ST0), function(i) is.na(title_of(i)), logical(1))]
if (length(.undoc)) stop("ST0 has no title row for: ", paste(.undoc, collapse = ", "))
# ST18 is not in ST0 (built after the dictionary); document it here.
ST18_TITLE  <- "Sequencing-depth diagnostics: per-compartment NHD-to-control depth ratios, the UMI-floor sweep, the depth-only artifact by normalization, the neuronal-program depth artifact and the myeloid depth comparison (six sheets, A-F)."
ST18_BACKS  <- "Methods > Sequencing depth; Fig 2h sensitivity arm; Fig 5d depth caveat"
ST18_SOURCE <- "96_depth_diagnostics_FH.R (A-E); 95_TF_depth_sensitivity_FH.R (F)"
ST18_DESC   <- data.table(sheet = depth_sheet, description = c(
  "A: cell-type median UMI divided by the median UMI of its own library — isolates compartment RNA content from sequencing depth.",
  "B: what each UMI floor from 400 to 1,500 costs (retention per condition) and buys (residual depth-only artifact).",
  "C: depth-only artifact remaining under raw / log-normalized / thinned treatment, per region.",
  "D: NHD-over-control median depth ratio per compartment and region, with the balance call.",
  "E: for each Fig. 5d program, the loss measured against the loss depth alone can produce (control split at the same depth gap).",
  "F: myeloid depth per condition and region for the Azimuth class and the analysed microglia, with the Kolmogorov-Smirnov comparison."))

# Referee-facing titles: the ST0 titles are working titles and were leaking onto the SI
# contents page and the workbook README. TITLE_PUBLIC (from _sd_map_FH.R) is what the referee reads; the internal
# title stays in the build index as `title_internal`.
# column dictionary for the six depth sheets (they post-date ST0)
ST18_COLS <- rbindlist(list(
  data.table(sheet = "A", column = c("Condition","Region","lane","cell_type","n_nuclei","celltype_median_UMI","lane_median_UMI","ratio_to_lane"),
             definition = c("CON or NHD","Frontal or Hippo","library (TFHS id)","atlas cell type","nuclei in that cell type x library","median UMI per nucleus in the cell type","median UMI per nucleus over the whole library","celltype_median_UMI / lane_median_UMI — compartment RNA content net of sequencing depth")),
  data.table(sheet = "B", column = c("UMI_floor","keep_CON","keep_NHD","n_CON_kept","n_NHD_kept","residual_depth_artifact","pct_artifact_over_0.15","retention_asymmetry_pp"),
             definition = c("minimum UMI per nucleus tested","fraction of control myeloid nuclei retained","fraction of NHD myeloid nuclei retained","control nuclei retained","NHD nuclei retained","median |Cliff's delta| produced by depth alone after thinning to this floor","percent of genes whose depth-only delta exceeds 0.15","keep_NHD minus keep_CON, percentage points")),
  data.table(sheet = "C", column = c("region","treatment","depth_gap_x","n_genes","median_abs_delta","pct_over_0.15","note"),
             definition = c("Frontal or Hippo","raw_counts, log_normalized or thinned","NHD/CON median-depth ratio in the split","genes evaluated","median |Cliff's delta| attributable to depth alone","percent of genes with depth-only |delta| > 0.15","free text")),
  data.table(sheet = "D", column = c("Region","cell_type","med_CON","med_NHD","n_CON","n_NHD","ratio_NHD_over_CON","balance"),
             definition = c("Frontal or Hippo","atlas cell type","median UMI, control nuclei","median UMI, NHD nuclei","control nuclei","NHD nuclei","med_NHD / med_CON","matched, NHD deeper or NHD shallower")),
  data.table(sheet = "E", column = c("region","program","depth_ratio_NHD_over_CON","control_split_x","n_per_arm","real_delta","depth_only_artifact","separable"),
             definition = c("Frontal or Hippo","Fig. 5d program","NHD/CON median-depth ratio in excitatory neurons","depth gap reproduced within control nuclei","nuclei per arm of the control split","observed NHD-minus-control Cliff's delta","delta produced by the depth gap alone","TRUE when |real_delta| exceeds the depth-only artifact")),
  data.table(sheet = "F", column = c("region","cell_set","n_CON","n_NHD","median_UMI_CON","median_UMI_NHD","ratio_NHD_over_CON","KS_D","KS_p"),
             definition = c("Frontal or Hippo","Micro-PVM (Azimuth class) or Microglia (analysis set)","control nuclei","NHD nuclei","median UMI, control","median UMI, NHD","median_UMI_NHD / median_UMI_CON","Kolmogorov-Smirnov D between the two depth distributions","Kolmogorov-Smirnov P value"))))

ST19_TITLE  <- unname(TITLE_PUBLIC[["ST19"]])
ST19_BACKS  <- "Supplementary Fig. 2e"
ST19_SOURCE <- paste0("105_zhou_pseudobulk_cache_FH.R; 106_zhou_crosscohort_arms_FH.R (statistic: _zhou_crosscohort_core_FH.R, shared with ", sd_ref("ST13"), ")")
ST19_COLS <- rbind(
  dict[table == "ST13", .(column, type, permitted_values, definition)],
  data.table(column = c("genomewide_sign_agreement","arm","arm_label","n_NHD_zhou","n_CON_zhou","n_zhou_sig_arm","rho_full_mle","rho_full_wald","p_floored"),
             type = c("numeric","character","character","integer","integer","integer","numeric","numeric","logical"),
             permitted_values = c("0-1", "all3 | no861 | noNHD1 | noNHD2 | only861 | onlyNHD1 | onlyNHD2", "free text", "1-3", "10-11", "count", "-1 to 1", "-1 to 1", "TRUE | FALSE"),
             definition = c("Fraction of all shared genes agreeing in sign between the two cohorts (identical to baseline_concord; the threshold-free direction measure, and the quantity to read when n_both_sig = 0).",
                            "Which Zhou NHD donors form the NHD side of the pseudobulk DESeq2 contrast; the control side is the same in every arm.",
                            "Plain-language description of the arm; donor 861 = Zhou NHD3 (59-year-old female, PMI 36 h, TYROBP c.2T>C), the donor shared with this study. 1 v 11 arms have no NHD replicate (dispersion from controls) and are descriptive.",
                            "Number of Zhou NHD donors in the arm.",
                            "Number of Zhou control donors in the arm: 10 for Micro-PVM and OPC (one control has fewer than 10 nuclei of that type), otherwise 11.",
                            "Zhou genes with padj < 0.10 in this arm (the Zhou-side significance gate); inflated in 1 v 11 arms.",
                            "Spearman rho between our per-nucleus log2FC and Zhou's unshrunken DESeq2 log2FC over all shared genes — the metric comparable across arms (apeglm shrinkage attenuates lower-n arms more).",
                            "Spearman rho against the Zhou Wald statistic.",
                            "TRUE where any p value of the row (rho_full_p, rho_sig_p, tau_sig_p or binom_p) was below 2.2250739e-308 before flooring; such values are written at the floor and must be read as p < 2.2e-308, never as a measured value. rho_full_p treats genes as independent and is not evidence; binom_p is one-sided, nominal and uncorrected across rows; arms share controls and are not compared by test.")))
ST19B_COLS <- data.table(column = c("cell_type","donor_a","donor_b","n_genes","rho_apeglm","rho_mle"), type = c("character","character","character","integer","numeric","numeric"),
                         permitted_values = c("Zhou cell type", "861 | NHD1 | NHD2", "861 | NHD1 | NHD2", "count", "-1 to 1", "-1 to 1"),
                         definition = c("Cell type", "First Zhou NHD donor of the pair, profiled alone against the 11 controls (861 = Zhou NHD3, the donor shared with this study)", "Second Zhou NHD donor of the pair", "Genes tested in both 1 v 11 profiles",
                                        "Spearman rho between the two donors' apeglm log2FC profiles (within the Zhou cohort; no data from this study)", "Same on unshrunken log2FC"))

# ST20 = Visium microglia-proximity table (script 111; post-dates ST0). Constants quoted from 111:
# NBOOT = 2000 block-bootstrap resamples, NSHIFT = 200 toroidal shifts, BLOCK = 6 array rows x 12 array columns.
ST20_TITLE  <- unname(TITLE_PUBLIC[["ST20"]])
ST20_BACKS  <- "Fig. 6f"
ST20_SOURCE <- "111_visium_proximity_severity_FH.R (module scores: 65_visium_prep_depthmatched_FH.R; abundances: cell2location)"
ST20_COLS <- data.table(
  column = c("set","version","compartment","outcome","primary","n_proximal","n_distal","n_blocks","rank_median_proximal","rank_median_distal",
             "delta","ci_lo","ci_hi","null_lo","null_hi","null_percentile","delta_oligo_stratified","median_proximal","median_distal","condition"),
  type = c("character","character","character","character","logical","integer","integer","integer","numeric","numeric",
           "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","character"),
  permitted_values = c("NHD_Frontal1 | NHD_Frontal2 | NHD (both) | CON_Frontal1 | CON_Frontal2 | CON (both)",
                       "fraction | abundance | fraction_neighborhood | abundance_neighborhood",
                       "Grey matter | White matter",
                       "Astrocyte reactive | Complement / MHC-II | Antigen presentation | Cholesterol (sterol arm) | Structural myelin | Fatty-acid / sphingomyelin | Lipid uptake / salvage | Glycolytic shift | FTL/HAMP (iron) | SPP1 (DAM-2 arm) | Housekeeping (null) | Homeostatic microglia | Astrocyte abundance (c2l) | Oligodendrocyte abundance (c2l) | Excitatory-neuron abundance (c2l)",
                       "TRUE | FALSE", "count", "count", "count", ">= 0", ">= 0",
                       "-1 to 1", "-1 to 1", "-1 to 1", "-1 to 1 or NA", "-1 to 1 or NA", "0-1 or NA", "-1 to 1 or NA", "numeric", "numeric", "NHD | CON"),
  definition = c(
    "Visium section (NHD_Frontal1, NHD_Frontal2, CON_Frontal1, CON_Frontal2) or the two sections of one condition pooled ('NHD (both)', 'CON (both)'); quartiles are always formed within section x spatial domain before pooling.",
    "Ranking variable used to define microglia-rich versus microglia-poor spots: fraction = cell2location microglial q05 abundance divided by the total q05 abundance over the eight reference cell types (the primary ranking; absolute abundance in control tissue is a cellularity index); abundance = absolute microglial q05 abundance (sensitivity); fraction_neighborhood / abundance_neighborhood = the mean of the same quantity over the six hexagonal array neighbors of the spot (spot itself excluded; at least five of six neighbors present, otherwise NA and the spot is dropped), testing 'near microglia-rich tissue' at one-spot distance.",
    sprintf("Grey matter = the spatial domains flagged as grey matter in integrated_domain_levels.csv (%s); White matter = the %s domain. White matter rows exist only for the NHD sets that reach 60 depth-matched white-matter spots (NHD_Frontal1 and 'NHD (both)'): the control white-matter spots fall below the 3,000-UMI depth floor of the module scores, and the NHD_Frontal2 white matter falls below 60 spots once its sulcal layer I is assigned to L1/pia (set sizes: n_proximal + n_distal are the two quartiles of the set).",   # domain names from the levels file
            .list_and(DOM_GREY), paste(DOM_WHITE, collapse = " / ")),
    "Outcome measured per spot: an AddModuleScore module score on 3,000-UMI-thinned counts (script 65; depth removed) or a cell2location q05 cell-type abundance ('... abundance (c2l)' rows, the composition confounds). 'FTL/HAMP (iron)' and 'SPP1 (DAM-2 arm)' are the iron-handling and TREM2-dependent-lipid modules named for the single genes that carry them in Visium; 'Housekeeping (null)' is the negative-control module; 'Homeostatic microglia' (P2RY12, CX3CR1, CSF1R, TMEM119, SALL1) is the reference that separates more microglia from a different microglial state.",
    "TRUE for the four pre-specified primary contrasts (set = 'NHD (both)', version = fraction, compartment = Grey matter, outcome in {Astrocyte reactive, Complement / MHC-II, Cholesterol (sterol arm), Structural myelin}); every other row is exploratory.",
    "Spots in the microglia-rich group: the top quartile of the ranking variable (rank percentile >= 0.75), quartiles formed within section x spatial domain and then pooled across the set.",
    "Spots in the microglia-poor group: the bottom quartile of the ranking variable (rank percentile <= 0.25), formed as for n_proximal.",
    "Number of tissue blocks in the set x compartment used as the resampling unit of the block bootstrap: blocks of 6 array rows x 12 array columns within a section (BLOCK in script 111).",
    "Median of the ranking variable (see version) among the microglia-rich spots.",
    "Median of the ranking variable among the microglia-poor spots.",
    "Cliff's delta of the outcome, microglia-rich minus microglia-poor spots (probability that a rich spot exceeds a poor spot minus the reverse; midranks for ties). NA when either group has fewer than five spots.",
    "Lower bound of the 95 % block-bootstrap interval of delta: tissue blocks resampled with replacement 2,000 times (NBOOT), delta recomputed between the two quartile groups each time, 2.5th percentile. A within-donor sampling interval, not a donor-level confidence interval; NA when fewer than eight blocks.",
    "Upper bound of the 95 % block-bootstrap interval (97.5th percentile), as ci_lo.",
    "2.5th percentile of the toroidal-shift null: the ranking map is displaced on the array by 200 random parity-preserving shifts (NSHIFT; even row and column offsets, wrapped within the section's array extent), quartiles re-formed within section x domain and delta recomputed each time. Computed for the fraction and abundance versions only; NA for the neighborhood versions.",
    "97.5th percentile of the toroidal-shift null, as null_lo.",
    "Fraction of the 200 null deltas that are <= the observed delta (1 = the observed delta exceeds every shifted-map delta; 0 = it is below all of them). Reported as a percentile of the null, not as a P value: condition is perfectly aliased with donor, so no P value is given anywhere in this table.",
    "Cliff's delta recomputed within deciles of cell2location oligodendrocyte abundance (deciles over both groups pooled; strata with >= 5 spots per group, weighted by stratum size) — the composition control against grey/white-matter mixing. NA for the oligodendrocyte-abundance outcome itself and when fewer than two usable strata.",
    "Median of the outcome among the microglia-rich spots.",
    "Median of the outcome among the microglia-poor spots.",
    "Condition of the set: NHD or CON."))
stopifnot(nrow(ST20_COLS) == 20, !anyDuplicated(ST20_COLS$column))

# ST21 = Visium domain composition (script 114; post-dates ST0). Constants quoted from 114:
# NBOOT = 2000 block-bootstrap resamples, MIN_SPOTS = 30 spots per condition, BLOCK = 6 array rows x 12 array columns.
ST21_TITLE  <- unname(TITLE_PUBLIC[["ST21"]])
ST21_BACKS  <- "Fig. 6g"
ST21_SOURCE <- "114_visium_domain_composition_FH.R (abundances: cell2location; domains: integrated_cluster_identity.csv + integrated_domain_levels.csv)"
ST21_COLS <- data.table(
  column = c("cell_type","cell_type_label","reference_inflated","domain","domain_order","compartment","n_spots_CON","n_spots_NHD","n_blocks_CON","n_blocks_NHD",
             "mean_prop_CON","mean_prop_NHD","mean_prop_CON_Frontal1","mean_prop_CON_Frontal2","mean_prop_NHD_Frontal1","mean_prop_NHD_Frontal2",
             "log2_ratio","ci_lo","ci_hi","included","nboot","min_spots"),
  type = c("character","character","logical","character","integer","character","integer","integer","integer","integer",
           "numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","logical","integer","integer"),
  permitted_values = c("Astro | Endo | Micro-PVM | Neuron_Ex | Neuron_Inh | OPC | Oligo | Pericytes", "display name", "TRUE | FALSE",
                       paste(DOM_LEV, collapse = " | "), sprintf("1-%d", length(DOM_LEV)), "grey | white | other", "count", "count", "count", "count",   # from integrated_domain_levels.csv
                       "0-1", "0-1", "0-1 or NA", "0-1 or NA", "0-1 or NA", "0-1 or NA", "numeric or NA", "numeric or NA", "numeric or NA", "TRUE | FALSE", "2000", "30"),
  definition = c(
    "cell2location reference cell type (the eight types of the snRNA-seq atlas used as the deconvolution reference).",
    "Display name of the cell type (Micro-PVM, Astrocyte, Oligodendrocyte, OPC, Excitatory neuron, Inhibitory neuron, Endothelial, Pericyte).",
    "TRUE for Endo and Pericytes: the reference over-calls the vascular types (27.8 % + 32.9 % of control cortex), so their rows are shaded on the panel and are not interpreted.",
    "Spatial domain of the spot (integrated Banksy/Harmony cluster mapped through the curated identity table).",
    sprintf("Superficial-to-deep order of the domain (%s).", paste(sprintf("%d = %s", DOM_INFO$order, DOM_INFO$domain), collapse = ", ")),   # from the levels file
    sprintf("Compartment of the domain: grey (cortical layers: %s), white (%s) or other (%s).", paste(DOM_GREY, collapse = ", "), paste(DOM_WHITE, collapse = ", "), paste(DOM_OTHER, collapse = ", ")),
    "Control spots in the domain with a defined composition (total q05 abundance over the eight types > 0).",
    "NHD spots in the domain with a defined composition.",
    "Tissue blocks of 6 array rows x 12 array columns (within section) carrying the control spots of the domain; the resampling unit of the bootstrap.",
    "Tissue blocks carrying the NHD spots of the domain.",
    "Mean over control spots of the per-spot proportion of the type: q05 abundance of the type divided by the spot's total q05 abundance over the eight types.",
    "The same over NHD spots.",
    "The same over the spots of section CON_Frontal1 alone; NA when the section has no spot in the domain.",
    "The same for CON_Frontal2.",
    "The same for NHD_Frontal1.",
    "The same for NHD_Frontal2.",
    "log2(mean_prop_NHD / mean_prop_CON); NA when the domain is excluded (see included). Descriptive: condition is aliased with donor, so no P value is given.",
    "2.5th percentile of log2_ratio over 2,000 block-bootstrap draws (NBOOT): tissue blocks resampled with replacement within each condition, the two means and their ratio recomputed each draw. A within-donor sampling interval, not a donor-level confidence interval.",
    "97.5th percentile of the same bootstrap distribution.",
    "TRUE when the domain has at least 30 spots (MIN_SPOTS) in both conditions; only included domains are drawn in Fig. 6g.",
    "Number of bootstrap draws used (constant).",
    "Minimum spots per condition for a domain to carry a ratio (constant)."))
stopifnot(nrow(ST21_COLS) == 22, !anyDuplicated(ST21_COLS$column))

# ST22 = Visium spatial-domain annotation, one row per spot of the final object (script 125; post-dates ST0).
# Constants quoted from 120/119b/123: primary rule = BANKSY pial cluster within 300 um of the pial surface (1,039 spots moved),
# sensitivity rule = WM/L6 spots strictly inside the pathologists' layer-I band (151 spots); L4 IT signature over cortical spots.
ST22_TITLE  <- unname(TITLE_PUBLIC[["ST22"]])
ST22_BACKS  <- "Fig. 6a,b; Supplementary Fig. 4c"   # panel e (annotation overlay) retired
ST22_SOURCE <- "125_visium_domain_annotation_table_FH.R (domains: integrated_cluster_identity.csv + integrated_domain_levels.csv; layer I / pia refinement: 120_visium_L1_histology_refinement_FH.R)"
ST22_COLS <- data.table(
  column = c("spot_id","section","condition","seurat_cluster","domain","domain_order","compartment","tissue_call_source",
             "moved_to_L1pia","previous_domain","banksy_pial","interior_WM_L6","d_L1_line_um","d_WM_line_um","inside_L1_band"),
  type = c("character","character","character","integer","character","integer","character","character",
           "logical","character","logical","logical","numeric","numeric","logical"),
  permitted_values = c("<section>_<barcode>", "CON_Frontal1 | CON_Frontal2 | NHD_Frontal1 | NHD_Frontal2", "CON | NHD", "0-10",
                       paste(DOM_LEV, collapse = " | "), sprintf("1-%d", length(DOM_LEV)), "grey | white | other", "spaceranger | he_mask",   # domains from integrated_domain_levels.csv
                       "TRUE | FALSE", paste(DOM_LEV, collapse = " | "), "TRUE | FALSE", "TRUE | FALSE", ">= 0", ">= 0", "TRUE | FALSE"),
  definition = c(
    "Spot identifier: the Visium section followed by the 10x barcode; the column name of the integrated object.",
    "Visium section (capture area). CON_Frontal1 and CON_Frontal2 are from the control donor, NHD_Frontal1 and NHD_Frontal2 from the NHD donor; one donor per condition.",
    "Condition of the section: CON or NHD.",
    "Final integrated cluster of the spot (integrated Banksy/Harmony clustering, resolution 0.5; clusters 0-9), or 10 = the layer I / pia refinement cluster that received the spots moved by the primary rule.",
    "Spatial domain of the spot: the cluster mapped through the curated identity table (the domain used in every Visium panel and table).",
    sprintf("Superficial-to-deep order of the domain (%s).", paste(sprintf("%d = %s", DOM_INFO$order, DOM_INFO$domain), collapse = ", ")),   # from the levels file
    sprintf("Compartment of the domain: grey (cortical layers: %s), white (%s) or other (%s).", paste(DOM_GREY, collapse = ", "), paste(DOM_WHITE, collapse = ", "), paste(DOM_OTHER, collapse = ", ")),
    "Tier under which the spot entered the analysis: spaceranger = called in-tissue by Space Ranger; he_mask = outside the Space Ranger call but under tissue on the H&E mask (the pial rim the Space Ranger call had dropped). Spots under neither tier are not in the object.",
    "TRUE when the spot was re-assigned to L1/pia by the primary rule (BANKSY-primary): member of a pial BANKSY cluster (a cluster with at least 30 % of its spots within 100 um of the pia) and lying inside the pathologists' closed layer-I band or within 300 um of a layer-I line. 1,039 spots; the rule is applied identically to all four sections and is reversible through previous_domain.",
    "Domain the spot carried before the refinement (its expression-only cluster identity). Equal to domain for every spot that was not moved.",
    "TRUE when the spot satisfies the BANKSY-pial membership term of the primary rule. Under the shipped (BANKSY-primary) refinement this column equals moved_to_L1pia; it is kept so the two terms of the rule stay separable.",
    "TRUE under the histology-only sensitivity rule: a spot labeled WM or L6 before refinement that lies strictly inside the closed band between the pathologists' two layer-I lines. 151 spots; every domain-dependent value can be recomputed with these spots, rather than the moved_to_L1pia set, re-assigned to L1/pia.",
    "Distance in um from the spot center to the nearest pathologists' layer-I line of its section (H&E annotation drawn blind to the transcriptome).",
    "Distance in um from the spot center to the nearest pathologists' white-matter line of its section.",
    "TRUE when the spot lies inside the closed band between the pathologists' two layer-I lines (the band the sensitivity rule uses)."))
stopifnot(nrow(ST22_COLS) == 15, !anyDuplicated(ST22_COLS$column))
ST22B_COLS <- data.table(
  column = c("section","n_cortical","mean_L4IT_per10k","mean_RORB_per10k","frac_L4high","neighbor_L4high_given_L4high",
             "coherence_ratio","morans_I","median_depth_um_L4high","n_L4_domain"),
  type = c("character","integer","numeric","numeric","numeric","numeric","numeric","numeric","numeric","integer"),
  permitted_values = c("CON_Frontal1 | CON_Frontal2 | NHD_Frontal1 | NHD_Frontal2", "count", ">= 0", ">= 0", "0-1", "0-1", ">= 0", "-1 to 1", ">= 0 um", "count"),
  definition = c(
    "Visium section; one donor per condition, so every value on this sheet is descriptive and none is tested.",
    "Cortical spots of the section: all spatial domains except WM and Vasc/immune (the spots over which the signature, the graph and the statistics are computed).",
    "Mean over the cortical spots of the layer IV IT signature: RORB + TWIST2 + DPF3 + NGB + SCHLAP1 + LYPD6B UMI per 10,000 UMI of the spot (RORB plus the reference's L4 IT-over-L5 IT genes).",
    "Mean over the cortical spots of RORB alone, UMI per 10,000.",
    "Fraction of the section's cortical spots that are L4-high: signature at or above the 85th percentile of the pooled cortical distribution of all four sections (the top 15 %). Under random placement this is the expected fraction of an L4-high spot's neighbors that are also L4-high.",
    "Neighbor co-occurrence: the mean, over the L4-high spots, of the fraction of their six nearest cortical neighbors that are also L4-high, on a six-nearest-neighbor graph built over the cortical spots of the section and symmetrized.",
    "neighbor_L4high_given_L4high divided by frac_L4high: the neighbor co-occurrence relative to chance (1 = random placement; larger = the L4-high spots form a band).",
    "Moran's I of the signature on the same symmetrized six-nearest-neighbor graph over the cortical spots (spatial autocorrelation; 0 = no spatial structure).",
    "Median distance in um of the L4-high spots to the nearest WM spot of the section (cortical depth of the L4-high spots).",
    "Number of spots of the section in the L4 spatial domain (for comparison with the signature-based measures)."))
stopifnot(nrow(ST22B_COLS) == 10, !anyDuplicated(ST22B_COLS$column))

# ST23 = Visium program expression per spatial domain, three sheets (script 126; post-dates ST0). Copies of the panel
# tables: sheet Programs_per_domain = 67b (Fig. 6e; per spot the program UMIs per 10,000 UMIs as a ratio of sums over the spot and its
# present hexagonal array neighbors, six offsets within section, native depth; then mean/median/n per domain), sheet
# Genes_by_band = 73 (Fig. 6d; per-spot gene UMI per 10,000 at native depth, no neighbor pooling, mean per band x
# condition, log2((NHD + 1e-4) / (CON + 1e-4))), sheet Gene_sets_Visium = 65 (module membership as scored, `present`).
ST23_TITLE  <- unname(TITLE_PUBLIC[["ST23"]])
ST23_BACKS  <- "Fig. 6d,e"
ST23_SOURCE <- "126_visium_program_tables_FH.R (copies of the Fig. 6e, Fig. 6d and module-membership tables of 67b_visium_fig6e_alldepth_FH.R, 73_visium_supplychain_heatmap_FH.R and 65_visium_prep_depthmatched_FH.R; programs defined in 65_visium_prep_depthmatched_FH.R)"
# Gene lists quoted in the ST23 dictionary are read from the shipped Gene_sets sheet, never typed (the scored
# Complement / MHC-II module is not the Fig. 6d block; the two differ by member)
.f23c <- file.path(SRC, "ST23c_visium_program_gene_sets_FH.csv"); stopifnot("MISSING ST23c — run scripts/126_visium_program_tables_FH.R" = file.exists(.f23c))
.gs23 <- fread(.f23c); stopifnot(identical(names(.gs23), c("program", "gene", "detected_on_slide")), is.logical(.gs23$detected_on_slide))
.members <- function(p) paste(.gs23[program == p & detected_on_slide == TRUE]$gene, collapse = ", ")
ST23_COLS <- data.table(
  column = c("program","condition","domain","domain_order","mean_umi_per_10k","median_umi_per_10k","n_spots"),
  type = c("character","character","character","integer","numeric","numeric","integer"),
  permitted_values = c("Complement / MHC-II | Fatty-acid / sphingomyelin", "CON | NHD",
                       paste(DOM_LEV, collapse = " | "), sprintf("1-%d", length(DOM_LEV)), ">= 0", ">= 0", "count"),   # domains from integrated_domain_levels.csv
  definition = c(
    sprintf("One of the two programs mapped in Fig. 6e, scored with the members listed in sheet Gene_sets_Visium: Complement / MHC-II (%s) and Fatty-acid / sphingomyelin (%s).",
            .members("Complement / MHC-II"), .members("Fatty-acid / sphingomyelin")),
    "Condition of the spots: CON or NHD. One donor per condition, two sections each; the two sections of a condition are pooled, so every value is descriptive and none is tested.",
    paste0("Spatial domain of the spot (the final integrated Banksy/Harmony cluster mapped through the curated identity table; the domain of ", sd_ref("ST22"), "). All eight domains are listed, including the two that are not laminar bands (InN, Vasc/immune)."),
    sprintf("Superficial-to-deep order of the domain (%s).", paste(sprintf("%d = %s", DOM_INFO$order, DOM_INFO$domain), collapse = ", ")),   # from the levels file
    "Mean over the spots of the domain and condition of the per-spot program signal: UMIs of the program's genes per 10,000 UMIs of the spot at native sequencing depth (no UMI floor; every spot of the final object), pooled over the spot and its six hexagonal array neighbors as a ratio of sums (program UMIs over total UMIs of the spot plus its present neighbors, within section), so that a single UMI in a shallow white-matter spot cannot dominate. Not log-transformed. Comparable between conditions within a domain only in the sense of Fig. 6e (the domains are an outcome of the pathology; spot counts per domain differ between conditions, see n_spots).",
    "Median over the same spots of the same per-spot neighborhood signal; read beside the mean where shallow spots draw the mean up (control white matter).",
    paste0("Spots of the domain and condition entering the mean and median (every spot of the final object; the per-domain counts equal those of ", sd_ref("ST22"), ").")))
stopifnot(nrow(ST23_COLS) == 7, !anyDuplicated(ST23_COLS$column))
.f23b <- file.path(SRC, "ST23b_visium_program_genes_by_band_FH.csv"); stopifnot("MISSING ST23b — run scripts/126_visium_program_tables_FH.R" = file.exists(.f23b))
.gb23 <- fread(.f23b); stopifnot("program" %in% names(.gb23))
ST23B_COLS <- data.table(
  column = c("program","gene","band","band_order","mean_umi_per_10k_CON","mean_umi_per_10k_NHD","log2_ratio_NHD_vs_CON"),
  type = c("character","character","character","integer","numeric","numeric","numeric"),
  permitted_values = c(paste(unique(.gb23$program), collapse = " | "),   # the blocks drawn, read from the shipped sheet in panel order
                       "HGNC symbol", paste(DOM_BANDS, collapse = " | "), paste(DOM_INFO$order[match(DOM_BANDS, DOM_INFO$domain)], collapse = " | "), "> 0", ">= 0", "numeric"),
  definition = c(
    sprintf("Pathway block of the Fig. 6d heatmap the gene belongs to (the %d blocks drawn, in panel order). The Complement / MHC-II and Reactive astrocyte blocks of the heatmap are wider than, and differ in membership from, the corresponding scored programs of sheet Gene_sets_Visium (Complement / MHC-II, Astrocyte reactive); the membership here is the heatmap's.", uniqueN(.gb23$program)),
    "Gene (HGNC symbol); every gene appears once per band. Genes whose control mean was 0 in any band are not drawn and not listed.",
    paste0("Laminar band of the heatmap: the five cortical-layer domains and white matter (the six bands drawn in Fig. 6d; InN and Vasc/immune are not bands). Same domain vocabulary as ", sd_ref("ST22"), "."),
    paste0("Superficial-to-deep order of the band (the domain_order of ", sd_ref("ST22"), ")."),
    "Mean over the control spots of the band of the per-spot gene expression: UMIs of the gene per 10,000 UMIs of the spot at native sequencing depth (every spot of the final object; no UMI floor, no neighbor pooling).",
    "The same over the NHD spots of the band.",
    "log2((mean_umi_per_10k_NHD + 1e-4) / (mean_umi_per_10k_CON + 1e-4)): the absolute log2 NHD-to-control ratio of the two means (the pseudocount protects the log where a mean is 0; it is negligible at the means listed). This is the number the text quotes; the heatmap color itself is the centered, within-block-scaled value and is not shipped. Magnitudes are comparable within a pathway block, not between blocks. Descriptive: condition is aliased with donor, so no P value is given."))
stopifnot(nrow(ST23B_COLS) == 7, !anyDuplicated(ST23B_COLS$column))
ST23C_COLS <- data.table(
  column = c("program","gene","detected_on_slide"),
  type = c("character","character","logical"),
  permitted_values = c(paste(unique(.gs23$program), collapse = " | "),   # read from the shipped sheet, in scoring order
                       "HGNC symbol", "TRUE | FALSE"),
  definition = c(
    sprintf("Program scored on tissue (AddModuleScore module of the Visium analysis; %d programs — those of Fig. 6e, the proximity contrasts of %s and the null module).", uniqueN(.gs23$program), sd_ref("ST20")),
    "Member gene of the program as defined before scoring (HGNC symbol); each program x gene pair appears once.",
    paste0("TRUE when the gene was present on the slide and therefore entered the program score; FALSE = member dropped by the detection intersect (absent from the integrated object's features) and not scored.",
           if (all(.gs23$detected_on_slide)) " In the shipped build every member is TRUE." else sprintf(" %d of %d members are FALSE.", sum(!.gs23$detected_on_slide), nrow(.gs23)))))
stopifnot(nrow(ST23C_COLS) == 3, !anyDuplicated(ST23C_COLS$column))

# ST24 = program-level Cliff's delta per cell type and region, two sheets (script 127; post-dates ST0). Copies of the four
# panel statistics tables (Fig. 2d = 2B, Fig. 3e = 3D, Fig. 4d = 3Fp, Fig. 5d = 4C2): per program x region (x lineage /
# class) the per-nucleus Cliff's delta, the Wilcoxon p, BH q within the panel, the house gate (|delta| >= 0.15) and the
# grey q* mark (q < 0.05 & |delta| < 0.15), both recomputed by 127 and asserted equal to the panel's own columns. Nucleus
# counts are NA for Fig. 2d because that panel table does not record them. Source p / q written as 0 are floored to the
# ST0 convention (2.2250739e-308 = p < 2.2e-308) by 127, which logs the count.
ST24_TITLE  <- unname(TITLE_PUBLIC[["ST24"]])
ST24_BACKS  <- "Fig. 2d; Fig. 3e,f; Fig. 4d; Fig. 5d"
ST24_SOURCE <- paste0("127_program_delta_table_FH.R (copies of the panel tables of 2B_micro_program_violins_FH.R, 3D_astro_state_violins_FH.R, 3Fp_oligo_lineage_block_FH.R and 4C2_neuron_programme_dotmatrix_FH.R; the gene membership of every program is listed in ", sd_ref("ST6"), ", under the program name as it appears here)")
ST24_COLS <- data.table(
  column = c("figure_panel","cell_type","program","region","n_CON","n_NHD","cliffs_delta","p_wilcox","q_BH","effect_gate","q_star"),
  type = c("character","character","character","character","integer","integer","numeric","numeric","numeric","logical","logical"),
  permitted_values = c("Fig. 2d | Fig. 3e | Fig. 4d | Fig. 5d", "Micro-PVM | Astro | Oligo | OPC | Neuron_Ex | Neuron_Inh", "panel axis label",
                       "Frontal | Hippocampus", "count or NA", "count or NA", "-1 to 1 (3 decimals)", ">= 2.2250739e-308, <= 1", ">= 2.2250739e-308, <= 1", "TRUE | FALSE", "TRUE | FALSE"),
  definition = c(
    "Main-figure panel whose statistics table the row is copied from: Fig. 2d (microglial curated programs), Fig. 3e (astrocyte signature scores), Fig. 4d (oligodendrocyte-lineage programs, oligodendrocytes and OPC side by side) or Fig. 5d (neuronal programs, excitatory and inhibitory neurons side by side). Every value is the panel's own; nothing is recomputed for this file except the two logical gates, which are asserted equal to the panel's.",
    paste0("Atlas cell type of the nuclei scored (the vocabulary of ", sd_ref("ST1"), " and ", sd_ref("ST9"), "). Micro-PVM = the bona-fide microglial nuclei of Fig. 2 (CD163+/F13A1+ cluster held out); Oligo and OPC are scored separately on their own objects and read side by side in Fig. 4d; Neuron_Ex and Neuron_Inh are the two neuronal classes of Fig. 5d."),
    paste0("Program or signature as labeled on the panel axis. Gene membership: all sets are listed in ", sd_ref("ST6"), " (column signature), under the program name used here, with two exceptions of naming only: the five Fig. 5d programs are listed there as Fig5d_Presynaptic, Fig5d_Postsynaptic, Fig5d_GluR_ionotropic, Fig5d_OxPhos and Fig5d_IEG; and of the eleven astrocyte signatures, the five published sets are listed under their curated names (Complement/IFN-reactive (Liddelow) = A1_reactive, Ischemic/S100A10-reactive (Zamanian) = A2_reactive, DAA (Habib) = DAA_Habib, Zhou NHD astrocyte = Zhou_NHD_astrocyte_up, STAT3 targets = STAT3_targets), each of whose scoring_context names the program label used here; the other six astrocyte sets and the eight oligodendrocyte-lineage programs are listed under exactly these names. Per-nucleus scores are the panel's: Seurat AddModuleScore (seed 42) for Figs. 2d and 3e, the mean log-normalized expression of the detected members for Figs. 4d and 5d. The Fig. 5d Presynaptic, Postsynaptic and OxPhos sets are not the same-named sets of ", sd_ref("ST7"), " (see the conventions of ", sd_ref("ST0"), ")."),
    "Region of the nuclei: Frontal (frontal cortex) or Hippocampus. The hippocampal NHD arm rests on a single library, as in every panel.",
    paste0("Control nuclei entering the comparison for that program, cell type and region. NA for the Fig. 2d rows: the microglial panel table does not record the per-group counts (the counts per region and condition are those of the Micro-PVM rows of ", sd_ref("ST9"), " and of the selection chain of ", sd_ref("ST17b"), "), and no value is inferred here."),
    "NHD nuclei entering the comparison; NA for the Fig. 2d rows for the same reason.",
    "Cliff's delta of the per-nucleus program score, NHD versus control, computed as 2*U/(n_NHD*n_CON) - 1 with midranks (positive = higher in NHD), rounded to 3 decimals. The statistical unit is the nucleus and there is one donor per condition, so delta is an effect size, not donor-level inference; it is the number printed beside each violin / dot of the panel.",
    "Two-sided Mann-Whitney / Wilcoxon rank-sum p of the per-nucleus scores, NHD versus control, normal approximation with midranks (the panel's test). With one donor per condition it ranks programs only and carries no claim. Values the panel script wrote as 0 (double underflow) are floored to 2.2250739e-308 and must be read as p < 2.2e-308, never p = 0.",
    "Benjamini-Hochberg q of p_wilcox, adjusted within the panel over all its program x region (x lineage or class) tests (14 tests for Fig. 2d, 22 for Fig. 3e, 32 for Fig. 4d, 20 for Fig. 5d); the panel's own column, recomputed and asserted by script 127. Same underflow floor as p_wilcox.",
    "TRUE when |cliffs_delta| >= 0.15, the project's effect-size floor (the strict gate of every figure); the same rule as the panel's passes column, asserted equal. Together with q_star it reproduces the panel annotation: effect_gate = delta printed as a claim, q_star = grey mark, neither = below the floor and not significant.",
    "TRUE when q_BH < 0.05 and |cliffs_delta| < 0.15: the grey 'q*' mark of the panels (statistically significant across thousands of nuclei but below the effect-size floor; a kept annotation, never an exclusion). Asserted equal to the panel's mark column."))
stopifnot(nrow(ST24_COLS) == 11, !anyDuplicated(ST24_COLS$column))
ST24B_COLS <- data.table(
  column = c("figure_panel","cell_type","program","cliffs_delta_Frontal","cliffs_delta_Hippocampus","q_BH_Frontal","q_BH_Hippocampus","direction_consistent"),
  type = c("character","character","character","numeric","numeric","numeric","numeric","logical"),
  permitted_values = c("Fig. 2d | Fig. 3e | Fig. 4d | Fig. 5d", "Micro-PVM | Astro | Oligo | OPC | Neuron_Ex | Neuron_Inh", "panel axis label",
                       "-1 to 1 (3 decimals)", "-1 to 1 (3 decimals)", ">= 2.2250739e-308, <= 1", ">= 2.2250739e-308, <= 1", "TRUE | FALSE"),
  definition = c(
    "As sheet Program_deltas.", "As sheet Program_deltas.", "As sheet Program_deltas; one row per figure_panel x cell_type x program, in panel order.",
    "cliffs_delta of sheet Program_deltas for the frontal cortex.", "cliffs_delta of sheet Program_deltas for the hippocampus.",
    "q_BH of sheet Program_deltas for the frontal cortex (BH within the panel; same underflow floor).", "q_BH of sheet Program_deltas for the hippocampus.",
    "TRUE when the frontal and hippocampal deltas have the same sign (no delta is exactly 0). Direction consistency across the two regions is what the text reads as the program-level result; it is never called replication (one donor per condition, and the hippocampal NHD arm is one library)."))
stopifnot(nrow(ST24B_COLS) == 8, !anyDuplicated(ST24B_COLS$column))

# ST3 second sheet = the detection-gated Micro-PVM fgsea behind Supplementary Fig. 3b (script 110b): same
# recipe as sheet Data, ranking restricted to genes not on the tissue artifact list and detected in >= 10 %
# of nuclei in both conditions; both regions. Columns are the ST3 columns minus the power-tier pair, plus:
ST3B_COLS <- rbind(
  dict[table == "ST3" & column %in% c("subtype","pathway","pval","padj","ES","NES","size","leadingEdge","direction","db"), .(column, type, permitted_values, definition)],
  data.table(column = "metric", type = "character", permitted_values = "log2fc_winsorized_gated",
             definition = "The pre-ranking metric: winsorized signed log2 fold change, computed on the gated gene set only."),
  data.table(column = c("n_ranked", "gate"), type = c("integer", "character"),
             permitted_values = c(">= 1", "not is_artifact & min(pct) >= 0.10"),
             definition = c("Number of genes in the gated ranking for that comparison.",
                            "The detection gate applied before ranking: gene not on the curated ambient/artifact list for Micro-PVM and detected in >= 10 % of nuclei in both conditions (Supplementary Fig. 3b; Methods).")))

# ---- README helpers ---------------------------------------------------------
.sep <- function(nm) data.table(column = sprintf("— sheet %s —", nm), type = "", permitted_values = "", definition = "")
# the column dictionary of a former table, with the "— sheet X —" separators carrying the NEW sheet names
dict_rows <- function(id) {
  sh <- sd_sheets_of(id)
  switch(id,
    ST18 = { stopifnot(identical(sh, depth_sheet))
             rbindlist(lapply(seq_along(sh), function(k) rbind(
               .sep(sh[k]),
               data.table(column = "(sheet)", type = "", permitted_values = "", definition = ST18_DESC$description[k]),
               ST18_COLS[sheet == substr(sh[k], 1, 1), .(column, type = "", permitted_values = "", definition)]))) },
    ST19 = rbind(.sep(sh[1]), ST19_COLS, .sep(sh[2]), ST19B_COLS),
    ST20 = ST20_COLS, ST21 = ST21_COLS,
    ST22 = rbind(.sep(sh[1]), ST22_COLS, .sep(sh[2]), ST22B_COLS),
    ST23 = rbind(.sep(sh[1]), ST23_COLS, .sep(sh[2]), ST23B_COLS, .sep(sh[3]), ST23C_COLS),
    ST24 = rbind(.sep(sh[1]), ST24_COLS, .sep(sh[2]), ST24B_COLS),
    ST3  = rbind(.sep(sh[1]), dict[table == id, .(column, type, permitted_values, definition)], .sep(sh[2]), ST3B_COLS),
    dict[table == id, .(column, type, permitted_values, definition)])
}
backs_of  <- function(id) switch(id, ST18 = ST18_BACKS, ST19 = ST19_BACKS, ST20 = ST20_BACKS, ST21 = ST21_BACKS, ST22 = ST22_BACKS, ST23 = ST23_BACKS, ST24 = ST24_BACKS, titles[table == id]$backs_figures[1])
source_of <- function(id) switch(id, ST18 = ST18_SOURCE, ST19 = ST19_SOURCE, ST20 = ST20_SOURCE, ST21 = ST21_SOURCE, ST22 = ST22_SOURCE, ST23 = ST23_SOURCE, ST24 = ST24_SOURCE, titles[table == id]$source_script[1])
title_int <- function(id) switch(id, ST18 = ST18_TITLE, ST19 = ST19_TITLE, ST20 = ST20_TITLE, ST21 = ST21_TITLE, ST22 = ST22_TITLE, ST23 = ST23_TITLE, ST24 = ST24_TITLE, title_of(id))

CONVENTIONS_LINE <- sprintf("See %s (data dictionary) for analysis-wide conventions: discovery criteria, Cliff's delta, the q* mark, depth matching, the p-underflow floor (2.2250739e-308 means p < 2.2e-308, never p = 0), and NA = not tested.", sd_ref("ST0"))

# README sheet of workbook n: a 4-column grid. Rows flagged `bold` are section heads. No internal ids: the referee
# addresses everything by workbook + sheet (the former ST ids stay in _build/ and in the dictionary's `table` column).
readme_grid <- function(n, ids) {
  row4 <- function(a = "", b = "", c = "", d = "", bold = FALSE) data.table(c1 = a, c2 = b, c3 = c, c4 = d, bold = bold)
  g <- rbind(row4("Label", sprintf("Supplementary Data %d", n)),
             row4("Title", SD_WB$title[n]),
             row4("Sheets", paste(SD_SHEETS$sheet[SD_SHEETS$workbook == n], collapse = ", ")),
             row4("Conventions", CONVENTIONS_LINE),
             row4())
  for (id in ids) {
    sh <- sd_sheets_of(id)
    g <- rbind(g,
               row4(sprintf("— %s %s —", if (length(sh) > 1) "sheets" else "sheet", paste(sh, collapse = ", ")), bold = TRUE),
               row4("Holds", unname(TITLE_PUBLIC[[id]])),
               row4("Backs figures", backs_of(id)),
               row4("Source script", source_of(id)),
               row4("column", "type", "permitted_values", "definition", bold = TRUE),
               { d <- dict_rows(id); data.table(c1 = d$column, c2 = d$type, c3 = d$permitted_values, c4 = d$definition, bold = grepl("^— sheet .* —$", d$column)) },
               row4())
  }
  g
}

# ---- packaging-time re-keying + guard ----------------------------------------
# Every character cell of every sheet (README included) passes through sd_rekey_old(): an unambiguous
# pre-consolidation citation ("Supplementary Data 12-27", "Supplementary Table N") is rewritten to the new
# workbook/sheet and COUNTED; the count is expected to be 0 because every upstream source renders through
# sd_ref(). A non-zero count is reported per sheet and the run stops after packaging with the list — the csv
# and the xlsx would otherwise disagree. Numbers 1-11 exist in both numberings and cannot be told apart by
# pattern; they are only ever rendered through sd_ref() upstream (grep-guarded in the report).
REKEY_LOG <- list()
rekey_dt <- function(d, where) {
  d <- copy(d); n <- 0L
  for (j in names(d)) if (is.character(d[[j]])) { v <- sd_rekey_old(d[[j]]); n <- n + attr(v, "n_rekeyed"); set(d, j = j, value = as.vector(v)) }
  if (n > 0) REKEY_LOG[[where]] <<- n
  d
}
OLD_NUM_PAT <- "Supplementary Data (1[2-9]|2[0-7])\\b"
count_old <- function(d) sum(vapply(d, function(v) if (is.character(v)) sum(grepl(OLD_NUM_PAT, v, perl = TRUE)) else 0L, numeric(1)))

hs <- createStyle(textDecoration = "bold")
write_wb <- function(n, tables) {   # tables = named list (former id -> named list of data sheets, new names)
  wb <- createWorkbook()
  g  <- rekey_dt(readme_grid(n, names(tables)), sprintf("Supplementary Data %d / README", n))
  addWorksheet(wb, "README")
  writeData(wb, "README", g[, .(c1, c2, c3, c4)], colNames = FALSE)
  for (r in which(g$bold)) addStyle(wb, "README", hs, rows = r, cols = 1:4, gridExpand = TRUE)
  setColWidths(wb, "README", cols = 1:4, widths = c(30, 16, 40, 90))
  expected <- SD_SHEETS$sheet[SD_SHEETS$workbook == n]; got <- unlist(lapply(tables, names), use.names = FALSE)
  stopifnot("sheet set / order of the workbook differs from the fixed map" = identical(got, expected),
            all(nchar(got) <= 31L), !anyDuplicated(got))
  for (id in names(tables)) for (nm in names(tables[[id]])) {
    d <- rekey_dt(tables[[id]][[nm]], sprintf("Supplementary Data %d / %s", n, nm))
    addWorksheet(wb, nm)
    writeData(wb, nm, d, headerStyle = hs)
    freezePane(wb, nm, firstRow = TRUE)
    tables[[id]][[nm]] <- d
  }
  p <- file.path(OUT, sd_file_n(n))
  saveWorkbook(wb, p, overwrite = TRUE)
  dims <- c(list(README = c(nrow(g), 4L)), lapply(unlist(tables, recursive = FALSE, use.names = FALSE), function(d) c(nrow(d) + 1L, ncol(d))))
  fix_xlsx_rels(p, dims = dims)
  list(path = p, readme = g, tables = tables)
}

# openxlsx writes, for every worksheet, a relationship to a drawing part (and a VML drawing) that it
# never creates. Excel tolerates the dangling reference; openpyxl (and any strict reader) raises
# KeyError on the missing part. Post-process the zip: drop the drawing / vmlDrawing relationships
# from xl/worksheets/_rels/sheetN.xml.rels and the matching <drawing/> / <legacyDrawing/> tags
# from the sheet xml. Done in R with utils::unzip / utils::zip; the file is rebuilt in a temp dir
# and moved into place only when zip succeeds.
# openxlsx also writes <dimension ref="A1"/> on every sheet; Excel and pandas ignore it, but a plain
# openpyxl read_only reader trusts it and sees one row. `dims` (named list, sheet name -> c(nrow, ncol) including
# the header row, in the order the sheets were added) rewrites the tag with the true used range.
fix_xlsx_rels <- function(path, dims = NULL) {
  td <- tempfile("xlsxfix_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(path, exdir = td)
  rels <- list.files(file.path(td, "xl", "worksheets", "_rels"), pattern = "\\.rels$", full.names = TRUE)
  n_rel <- 0L; n_tag <- 0L
  for (rf in rels) {
    x <- readLines(rf, warn = FALSE, encoding = "UTF-8"); x <- paste(x, collapse = "\n")
    x2 <- gsub('<Relationship [^>]*Target="[^"]*(drawings|vmlDrawing)[^"]*"[^>]*/>', "", x)
    n_rel <- n_rel + (nchar(x) - nchar(x2) > 0)
    if (!identical(x, x2)) writeLines(x2, rf, useBytes = TRUE)
  }
  sheets <- list.files(file.path(td, "xl", "worksheets"), pattern = "^sheet[0-9]+\\.xml$", full.names = TRUE)
  sheets <- sheets[order(as.integer(sub("^sheet([0-9]+)\\.xml$", "\\1", basename(sheets))))]   # sheet1 = first added (README)
  if (!is.null(dims)) stopifnot("dims must name every sheet of the workbook, in order" = length(dims) == length(sheets))
  n_dim <- 0L
  for (k in seq_along(sheets)) { sf <- sheets[k]
    x <- paste(readLines(sf, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    x2 <- gsub("<drawing [^>]*/>|<legacyDrawing [^>]*/>", "", x)
    n_tag <- n_tag + (nchar(x) - nchar(x2) > 0)
    if (!is.null(dims) && grepl('<dimension ref="A1"/>', x2, fixed = TRUE)) {
      x2 <- sub('<dimension ref="A1"/>', sprintf('<dimension ref="A1:%s%d"/>', openxlsx::int2col(dims[[k]][2]), dims[[k]][1]), x2, fixed = TRUE); n_dim <- n_dim + 1L }
    if (!identical(x, x2)) writeLines(x2, sf, useBytes = TRUE)
  }
  # the content-type overrides for the never-written drawing parts are dangling too
  ct <- file.path(td, "[Content_Types].xml")
  x <- paste(readLines(ct, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  x2 <- gsub('<Override [^>]*PartName="/xl/drawings/[^"]*"[^>]*/>', "", x)
  if (!identical(x, x2)) writeLines(x2, ct, useBytes = TRUE)
  # empty rels files left behind are legal; rebuild the zip (mimetype ordering is not required for xlsx)
  tmpzip <- tempfile(fileext = ".xlsx")
  old <- setwd(td); on.exit(setwd(old), add = TRUE)
  st <- utils::zip(tmpzip, files = list.files(".", recursive = TRUE, all.files = TRUE, include.dirs = FALSE), flags = "-r9Xq")
  setwd(old)
  stopifnot("zip failed while repairing the workbook" = st == 0, file.exists(tmpzip))
  file.copy(tmpzip, path, overwrite = TRUE); unlink(tmpzip)
  say("   repaired %s: %d sheet rels / %d sheet xml stripped of drawing references; %d dimension tags written", basename(path), n_rel, n_tag, n_dim)
  invisible(path)
}

# ---- loaders: one former table -> its data sheet(s) under the NEW sheet names, every guard kept ----------
.blank <- function(v) data.table(column = v, type = "", permitted_values = "", definition = "")
load_st <- function(id) {
  sh <- sd_sheets_of(id)
  if (id == "ST18") {
    sheets <- setNames(lapply(depth_files, fread), depth_sheet)
    # sheet E column `programme` -> `program` at packaging time (the source table of 96 keeps its name;
    # the README dictionary above documents the shipped name)
    .e <- grep("^E_", names(sheets)); stopifnot(length(.e) == 1, "programme" %in% names(sheets[[.e]]))
    setnames(sheets[[.e]], "programme", "program")
    for (.sh in names(sheets)) stopifnot(setequal(names(sheets[[.sh]]), ST18_COLS[sheet == substr(.sh, 1, 1)]$column))
    stopifnot(identical(names(sheets), sh))
    return(sheets)
  }
  d <- fread(st_files[[id]], colClasses = if (id %in% c("ST0","ST12")) "character" else NULL)
  if (id == "ST0") {
    # make the dictionary readable under the NEW labels without touching its join keys: every row gets the
    # workbook number + sheet that holds the table it documents; the label is "Supplementary Data N, sheet X"
    d[, workbook := NA_integer_]; d[, sheet := NA_character_]
    .in <- d$table %in% ST_SHIP
    d[.in, workbook := sd_wb_of(table)]
    d[.in, sheet := sd_primary_sheet(table)]
    d[, supplementary_data_label := fifelse(.in, vapply(table, function(t) if (t %in% ST_SHIP) sd_ref(t) else "", character(1)),
                                            fifelse(table == "Table1", "Table 1",
                                            fifelse(table %in% SD_DROPPED$st, "Not a Supplementary Data file: ships with the code release (reference_tables/)", table)))]
    # sd_ref() drops the sheet where the table owns the whole workbook (ST3 -> "Supplementary Data 5"); the label
    # must still name the sheet the documented columns live in, so the sheet is appended where sd_ref omitted it
    d[.in & !grepl(", sheet", supplementary_data_label, fixed = TRUE) & !workbook %in% sd_single_sheet_wb(),
      supplementary_data_label := sprintf("%s, sheet %s", supplementary_data_label, sheet)]
    setcolorder(d, c("table", "supplementary_data_label", "workbook", "sheet"))
    stopifnot(all(d[.in]$supplementary_data_label != ""), !anyNA(d[.in]$workbook),
              all(d[table %in% SD_DROPPED$st]$workbook %in% NA))
    say("   ST0: %d rows document the two GTRD reference tables (ST4/ST5) that ship with the code release, labeled as such", sum(d$table %in% SD_DROPPED$st))
    ptr <- function(id2, column, definition) data.table(
      table = id2, supplementary_data_label = { r <- sd_ref(id2); s2 <- sd_sheets_of(id2)
        if (grepl(", sheet", r, fixed = TRUE) || sd_wb_of(id2) %in% sd_single_sheet_wb()) r                    # sheet already named / single-sheet workbook
        else if (length(s2) > 1) sprintf("%s (all %d sheets)", r, length(s2))                                  # ST18 owns all six sheets of workbook 10
        else sprintf("%s, sheet %s", r, s2) },
      workbook = sd_wb_of(id2), sheet = sd_primary_sheet(id2),
      table_title = title_int(id2), backs_figures = backs_of(id2), source_script = source_of(id2),
      column = column, type = "mixed", permitted_values = sprintf("see README sheet of Supplementary Data %d", sd_wb_of(id2)), definition = definition)
    d <- rbind(d,
      ptr("ST19", sprintf("(columns of %s + arm, arm_label, n_NHD_zhou, ...)", sd_ref("ST13")),
          sprintf("Cross-cohort statistic per Zhou-donor arm; columns as %s plus the arm identifiers.", sd_ref("ST13"))),
      ptr("ST20", "(20 columns)", "Visium microglia-rich versus microglia-poor spot contrasts (Cliff's delta, block-bootstrap CI, toroidal-shift null) per section/condition, ranking version, compartment and outcome."),
      ptr("ST21", "(22 columns)", "Cell2location composition per Visium spatial domain and condition (mean proportion per reference type) with the log2 NHD/CON ratio and its block-bootstrap interval."),
      ptr("ST22", sprintf("(15 columns + sheet %s)", sd_sheets_of("ST22")[2]), "Per-spot Visium spatial-domain annotation of the final object (domain, tissue-call tier, layer I / pia refinement flags and line distances) with the per-section layer IV signal-and-order statistics as a second sheet."),
      ptr("ST23", sprintf("(7 columns + sheets %s)", paste(sd_sheets_of("ST23")[-1], collapse = ", ")), "Visium program expression per spatial domain: per-domain mean and median program UMI per 10,000 per condition (Fig. 6e), the per-gene band means and log2 NHD/control ratios of the Fig. 6d heatmap, and the gene membership of every program scored on tissue; copies of the panel tables."),
      ptr("ST24", sprintf("(11 columns + sheet %s)", sd_sheets_of("ST24")[2]), "Program-level Cliff's delta, Wilcoxon p, BH q (within panel), effect gate and q* mark per program x cell type x region for Figs. 2d, 3e, 4d and 5d, long form and with both regions side by side; copies of the panel tables."),
      ptr("ST18", "(six sheets A-F)", "Sequencing-depth diagnostics; each sheet has its own schema, documented in the README sheet of its workbook."),
      fill = TRUE)
    .pub <- d$table %in% names(TITLE_PUBLIC)
    d[.pub, table_title := unname(TITLE_PUBLIC[table])]
    stopifnot(all(d[table == "ST0"]$table_title == TITLE_PUBLIC[["ST0"]]), all(d[table == "Table1"]$table_title != ""))
    # the dictionary must place every shipped table exactly where the fixed map does
    stopifnot(setequal(unique(d[!is.na(workbook)]$table), ST_SHIP),
              all(d[!is.na(workbook), .(ok = workbook[1] == sd_wb_of(table[1]) & sheet[1] == sd_primary_sheet(table[1])), by = table]$ok))
  }
  if (id %in% c("ST20","ST21","ST22","ST23")) {   # the Visium tables must post-date the final Visium object
    .vobj <- file.path(dirname(PROJ), "Visium", "integrated_harmony", "NHD_frontal_integrated_harmony.rds")
    stopifnot("Visium object missing" = file.exists(.vobj),
              "Visium supplementary table is OLDER than NHD_frontal_integrated_harmony.rds — re-run 111 / 114 / 125 / 126" = file.mtime(st_files[[id]]) > file.mtime(.vobj)) }
  if (id == "ST20") {   # the README dictionary is hand-written here; it must describe the shipped columns and outcome names exactly
    stopifnot("ST20 white-matter sets differ from the README sentence (NHD_Frontal1 + 'NHD (both)')" = setequal(unique(d[compartment == "White matter"]$set), c("NHD_Frontal1", "NHD (both)")))
    stopifnot("ST20 columns differ from the README dictionary (ST20_COLS)" = setequal(names(d), ST20_COLS$column),
              "ST20 outcome names outside the dictionary's permitted values" =
                all(unique(d$outcome) %in% trimws(strsplit(ST20_COLS[column == "outcome"]$permitted_values, "|", fixed = TRUE)[[1]]))) }
  if (id == "ST21") {
    d[cell_type_label == "Microglia / PVM", cell_type_label := "Micro-PVM"]   # display name aligned with the atlas vocabulary at packaging time (source table of 114 unchanged)
    stopifnot("ST21 columns differ from the README dictionary (ST21_COLS)" = setequal(names(d), ST21_COLS$column),
              "ST21 domain names outside the dictionary's permitted values" =
                all(unique(d$domain) %in% trimws(strsplit(ST21_COLS[column == "domain"]$permitted_values, "|", fixed = TRUE)[[1]])))
    # ST21 must carry every current domain, and its spot counts must agree with the per-spot
    # annotation shipped as ST22 (ST21 counts spots with a defined composition, so <= ST22 and within 5 %)
    stopifnot("ST21 domain set differs from integrated_domain_levels.csv — re-run 114" = setequal(unique(d$domain), DOM_LEV))
    .a <- fread(file.path(SRC, "ST22_visium_domain_annotation_FH.csv"))[, .(n22 = .N), by = .(domain, condition)]
    .m <- merge(unique(d[, .(domain, n_spots_CON, n_spots_NHD)]), dcast(.a, domain ~ condition, value.var = "n22"), by = "domain")
    stopifnot("ST21 spot counts disagree with ST22 (stale 114 run?)" = all(.m$n_spots_CON <= .m$CON & .m$n_spots_NHD <= .m$NHD & .m$n_spots_CON >= 0.95 * .m$CON & .m$n_spots_NHD >= 0.95 * .m$NHD)) }
  if (id == "ST22") {
    stopifnot("ST22 columns differ from the README dictionary (ST22_COLS)" = setequal(names(d), ST22_COLS$column),
              "ST22 domain names outside the dictionary's permitted values" =
                all(unique(c(d$domain, d$previous_domain)) %in% trimws(strsplit(ST22_COLS[column == "domain"]$permitted_values, "|", fixed = TRUE)[[1]]))) }
  if (id == "ST23") {
    stopifnot("ST23 columns differ from the README dictionary (ST23_COLS)" = setequal(names(d), ST23_COLS$column),
              "ST23 domain names outside the dictionary's permitted values" =
                all(unique(d$domain) %in% trimws(strsplit(ST23_COLS[column == "domain"]$permitted_values, "|", fixed = TRUE)[[1]])),
              "ST23 must carry every current domain" = setequal(unique(d$domain), DOM_LEV),
              "ST23 programs outside the dictionary's permitted values" =
                all(unique(d$program) %in% trimws(strsplit(ST23_COLS[column == "program"]$permitted_values, "|", fixed = TRUE)[[1]]))) }
  sheets <- list(d)
  if (id == "ST3") { g3 <- file.path(TAB, "micro_gsea_gated_FH.csv"); stopifnot("MISSING — run scripts/110b_micro_gsea_gated_FH.R" = file.exists(g3))
    g3 <- fread(g3); stopifnot(setequal(names(g3), ST3B_COLS$column)); g3[, leadingEdge := gsub(";", "|", leadingEdge, fixed = TRUE)]   # pipe-separated, as sheet Data
    sheets <- c(sheets, list(g3)) }
  if (id == "ST19") { b19 <- file.path(SRC, "ST19b_zhou_within_cohort_donor_concordance.csv"); stopifnot(file.exists(b19))
    b19 <- fread(b19); stopifnot("ST19b columns differ from the README dictionary (ST19B_COLS)" = setequal(names(b19), ST19B_COLS$column))   # parity with ST22b
    sheets <- c(sheets, list(b19)) }
  if (id == "ST22") { b22 <- file.path(SRC, "ST22b_visium_L4IT_signal_vs_order_FH.csv"); stopifnot("MISSING ST22b — run scripts/125_visium_domain_annotation_table_FH.R" = file.exists(b22))
    b22 <- fread(b22); if ("neighbour_L4high_given_L4high" %in% names(b22)) setnames(b22, "neighbour_L4high_given_L4high", "neighbor_L4high_given_L4high")   # American spelling at packaging time (source table of 125 unchanged)
    stopifnot("ST22b columns differ from the README dictionary (ST22B_COLS)" = setequal(names(b22), ST22B_COLS$column))
    sheets <- c(sheets, list(b22)) }
  if (id == "ST23") { b23 <- fread(.f23b); c23 <- fread(.f23c)   # existence asserted where the dictionary reads them
    stopifnot("ST23b columns differ from the README dictionary (ST23B_COLS)" = setequal(names(b23), ST23B_COLS$column),
              "ST23b band names outside the dictionary's permitted values" =
                all(unique(b23$band) %in% trimws(strsplit(ST23B_COLS[column == "band"]$permitted_values, "|", fixed = TRUE)[[1]])),
              "ST23b bands must be the six Fig. 6d bands" = setequal(unique(b23$band), DOM_BANDS), all(b23$band %in% DOM_LEV),
              "ST23c columns differ from the README dictionary (ST23C_COLS)" = setequal(names(c23), ST23C_COLS$column),
              "ST23c programs outside the dictionary's permitted values" =
                all(unique(c23$program) %in% trimws(strsplit(ST23C_COLS[column == "program"]$permitted_values, "|", fixed = TRUE)[[1]])),
              "ST23 Data programs must be scored programs of ST23c" = all(unique(d$program) %in% c23$program))
    sheets <- list(d, b23, c23) }
  if (id == "ST24") {   # snRNA tables: not under the Visium-object mtime guard above
    .pv <- function(col) trimws(strsplit(ST24_COLS[column == col]$permitted_values, "|", fixed = TRUE)[[1]])
    stopifnot("ST24 columns differ from the README dictionary (ST24_COLS)" = setequal(names(d), ST24_COLS$column),
              "ST24 figure_panel outside the dictionary's permitted values" = all(unique(d$figure_panel) %in% .pv("figure_panel")),
              "ST24 must carry all four panels" = setequal(unique(d$figure_panel), .pv("figure_panel")),
              "ST24 cell_type outside the dictionary's permitted values" = all(unique(d$cell_type) %in% .pv("cell_type")),
              "ST24 region outside the dictionary's permitted values" = setequal(unique(d$region), .pv("region")),
              "ST24 must be 88 rows (14 + 22 + 32 + 20)" = nrow(d) == 88L,
              "ST24 p / q below the ST0 underflow floor (README says >= 2.2250739e-308)" = all(d$p_wilcox >= 2.2250739e-308) && all(d$q_BH >= 2.2250739e-308),
              "ST24 n_CON / n_NHD NA outside Fig. 2d (README says NA only there)" = !anyNA(d[figure_panel != "Fig. 2d", .(n_CON, n_NHD)]) && all(is.na(d[figure_panel == "Fig. 2d"]$n_CON)),
              "ST24 gates must be the stated rules" = identical(d$effect_gate, abs(d$cliffs_delta) >= 0.15) && identical(d$q_star, d$q_BH < 0.05 & abs(d$cliffs_delta) < 0.15))
    # The README above promises that every program is listed in the Signatures sheet (ST6) — check it against the shipped ST6
    .st6 <- fread(st_files[["ST6"]])
    .p_need <- unique(d[figure_panel %in% c("Fig. 2d", "Fig. 3e", "Fig. 4d")]$program)
    .p_have <- c(unique(.st6$signature), trimws(sub("^.* as program ", "", grep(" as program ", unique(.st6$scoring_context), value = TRUE))))
    stopifnot("ST24 program(s) with no row in ST6 (sheet Signatures) — re-run 90 ST6" = all(.p_need %in% .p_have),
              "ST24 Fig. 5d programs must all have a Fig5d_ set in ST6" = length(grep("^Fig5d_", unique(.st6$signature))) == length(unique(d[figure_panel == "Fig. 5d"]$program)))
    say("   ST24: all %d Fig. 2d / 3e / 4d programs have a row in ST6 (%d by name, %d via a stated panel label)", length(.p_need),
        sum(.p_need %in% unique(.st6$signature)), sum(!.p_need %in% unique(.st6$signature)))
    b24 <- file.path(SRC, "ST24b_program_cliffs_delta_wide_FH.csv"); stopifnot("MISSING ST24b — run scripts/127_program_delta_table_FH.R" = file.exists(b24))
    b24 <- fread(b24)
    stopifnot("ST24b columns differ from the README dictionary (ST24B_COLS)" = setequal(names(b24), ST24B_COLS$column),
              "ST24b must be 44 rows (7 + 11 + 16 + 10)" = nrow(b24) == 44L,
              "ST24b programs must be the Data programs" = setequal(paste(b24$figure_panel, b24$cell_type, b24$program), paste(d$figure_panel, d$cell_type, d$program)),
              "ST24b direction_consistent must be the sign rule" = identical(b24$direction_consistent, sign(b24$cliffs_delta_Frontal) == sign(b24$cliffs_delta_Hippocampus)))
    sheets <- list(d, b24) }
  if (id == "ST6") {   # the scoring_context vocabulary must be the map's (a stale ST6 would carry old numbers)
    .ok <- unname(SD_ST6_CTX); .v <- unique(d$scoring_context)
    .base <- unique(unlist(lapply(.v, function(v) if (v %in% .ok) v else sub(" as (column|program) .*$", "", strsplit(v, "; ", fixed = TRUE)[[1]]))))
    stopifnot("ST6 scoring_context outside SD_ST6_CTX — re-run 90 ST6" = all(.base %in% .ok)) }
  stopifnot("number of data sheets differs from the fixed map" = length(sheets) == length(sh))
  setNames(sheets, sh)
}

# ---- build -----------------------------------------------------------------
index <- list(); built <- list()
for (n in SD_WB$workbook) {
  ids <- unique(SD_SHEETS$st[SD_SHEETS$workbook == n])
  say("== Supplementary Data %d — %s  [%s]", n, SD_WB$name[n], paste(ids, collapse = ", "))
  tables <- setNames(lapply(ids, load_st), ids)
  res <- write_wb(n, tables)
  built[[n]] <- res
  md5 <- unname(md5sum(res$path))
  for (id in ids) for (nm in names(res$tables[[id]])) {
    d <- res$tables[[id]][[nm]]
    index[[length(index) + 1]] <- data.table(
      supplementary_data = sprintf("Supplementary Data %d", n), number = n, workbook_title = SD_WB$title[n], workbook_name = SD_WB$name[n],
      file = basename(res$path), sheet = nm, internal_id = id, former_supplementary_data = old_num(id),
      sheet_holds = unname(TITLE_PUBLIC[[id]]), title_internal = title_int(id), backs_figures = backs_of(id), source_script = source_of(id),
      rows = nrow(d), cols = ncol(d), md5 = md5, bytes = file.size(res$path))
  }
  # guard: no pre-consolidation number in the README or in any data sheet of this workbook
  .old <- count_old(res$readme[, .(c1, c2, c3, c4)]) + sum(vapply(unlist(res$tables, recursive = FALSE), count_old, numeric(1)))
  stopifnot("a pre-consolidation 'Supplementary Data 12-27' citation survives in this workbook" = .old == 0)
}
index <- rbindlist(index)
stopifnot(nrow(index) == nrow(SD_SHEETS), identical(paste(index$number, index$sheet), paste(SD_SHEETS$workbook, SD_SHEETS$sheet)))
fwrite(index, file.path(BUILD, "Supplementary_Data_INDEX.csv"))

if (length(REKEY_LOG)) {
  for (w in names(REKEY_LOG)) say("!! %d pre-consolidation citation(s) re-keyed at packaging time in %s", REKEY_LOG[[w]], w)
  stop("The packaging-time re-key fired (see above): an upstream source still types an old 'Supplementary Data' number. ",
       "Render it through sd_ref() (scripts/_sd_map_FH.R) and re-run 90, so the csv and the xlsx agree.")
} else say("packaging-time re-key: 0 pre-consolidation citations found in any sheet (every source renders through sd_ref())")

# dictionary + every README: the old numbers must be gone (the guard the task asks for, run on the written workbooks)
.dict_out <- built[[1]]$tables[["ST0"]][["Data"]]
stopifnot("old numbers survive in the packaged dictionary" = count_old(.dict_out) == 0,
          "old numbers survive in a README" = sum(vapply(built, function(b) count_old(b$readme[, .(c1, c2, c3, c4)]), numeric(1))) == 0)
say("guard: 'Supplementary Data 12-27' occurs 0 times in the dictionary and in the %d READMEs", length(built))

# ---- rename map + find/replace list for the Word draft ----------------------
map <- as.data.table(sd_rename_map())
map[, former_supplementary_data := old_num(internal_id)]
stopifnot(sum(map$old_form == "Supplementary Data N") == 27L, all(1:27 %in% map[old_form == "Supplementary Data N"]$former_supplementary_data))
fwrite(map[, .(old_label, old_form, internal_id, former_supplementary_data, new_label, new_sheet)], file.path(BUILD, "SD_rename_map.csv"))

# A naive sequential find/replace
# CASCADES, so the recipe is two passes through a
# placeholder token that cannot occur in the text.
sd_old <- map[old_form == "Supplementary Data N"][order(-former_supplementary_data)]
sd_old[, token := fifelse(grepl("^DROPPED", new_label), "@@SD-DROPPED@@", sub("^Supplementary Data (\\d+)$", "@@SD\\1@@", new_label))]
sd_tab <- map[old_form == "Supplementary Table N"][order(-nchar(old_label))]
sd_tab[, token := fifelse(grepl("^DROPPED", new_label), "@@SD-DROPPED@@", sub("^Supplementary Data (\\d+)$", "@@SD\\1@@", new_label))]
folded <- map[old_form == "Supplementary Data N" & !grepl("^DROPPED", new_label), .(old = paste(sub("Supplementary Data ", "", old_label), collapse = ", "), n_old = .N), by = new_label][n_old > 1][order(as.integer(sub("Supplementary Data ", "", new_label)))]
md <- c("# Find / replace list for the Word draft — 27-file numbering -> 11-workbook numbering (2026-09-21 consolidation)",
        "",
        sprintf("Generated by scripts/100_package_supplementary_data_FH.R on %s.", format(Sys.time(), "%Y-%m-%d %H:%M")),
        "",
        "## Read first",
        "* **Two files were dropped from the package**: the former Supplementary Data 5 (GTRD transcription-factor survey) and 6 (GTRD network).",
        "  Every citation of them must be REMOVED from the text (or re-pointed to the code release, `reference_tables/`). They map to the",
        "  token `@@SD-DROPPED@@` below so that they stand out after pass 1; search for that token and rewrite each sentence by hand.",
        "* **Several old files now share one workbook**, so the text should cite the sheet where it matters:",
        sprintf("  * %s <- former %s", folded$new_label, folded$old),
        "  Column `new_sheet` of `SD_rename_map.csv` gives the sheet for every old number.",
        "* **Plurals and ranges** ('Supplementary Data 16–18', 'Supplementary Data 23, 26', 'Supplementary Data 2 and 10') are NOT caught",
        "  by the single-number rules: re-key each number by hand with the map (16–18 -> 7; 23, 26 -> 11; 2 and 10 -> 3 and 2).",
        "* **Do not run a one-pass sequential replace**: old 2 -> 3 followed by old 3 -> 4 would re-hit the number just written.",
        "  Use the two passes below ('Match case' ON, 'Whole words' OFF, 'Use wildcards' OFF).",
        "",
        "## Pass 1 — old label -> placeholder token (IN THIS ORDER, longest number first: 'Supplementary Data 2' must not run before 'Supplementary Data 27')",
        "| find | replace |", "|---|---|",
        sprintf("| %s | %s |", sd_old$old_label, sd_old$token),
        sprintf("| %s | %s |", sd_tab$old_label, sd_tab$token),
        "",
        "## Pass 2 — placeholder token -> new label (any order)",
        "| find | replace |", "|---|---|",
        sprintf("| @@SD%d@@ | Supplementary Data %d |", SD_WB$workbook, SD_WB$workbook),
        "| @@SD-DROPPED@@ | (rewrite by hand: the GTRD tables are no longer Supplementary Data; see the code release, reference_tables/) |",
        "",
        "Bare suffix citations to check by hand: '(17b)' -> Supplementary Data 2, sheet Myeloid_selection; '(15b)' -> Supplementary Data 7, sheet RIN_gap.",
        "",
        "## Final map (one row per sheet)",
        "", "| Supplementary Data | sheet | former file | former Supplementary Data | rows | cols |", "|---|---|---|---|---|---|",
        sprintf("| %d | %s | %s | %d | %s | %d |", index$number, index$sheet, index$internal_id, index$former_supplementary_data, format(index$rows, big.mark = ","), index$cols))
writeLines(md, file.path(BUILD, "FIND_REPLACE_for_docx.md"))

# ---- manifest ---------------------------------------------------------------
wb_sum <- index[, .(n_sheets = .N, rows = sum(rows), cols = paste(cols, collapse = "/"), bytes = bytes[1], md5 = md5[1], sheets = paste(sheet, collapse = ", "), title = workbook_title[1]), by = .(number, file)]
writeLines(c("# Supplementary Data package (11 themed workbooks, 2026-09-21 consolidation)",
             sprintf("Built %s from supplementary_tables/ (90_build outputs; ST19 from 106; ST20 from 111; ST21 from 114; ST22 from 125; ST23 from 126; ST24 from 127) and tables/depth_A-F.", format(Sys.time(), "%Y-%m-%d %H:%M")),
             "Each .xlsx: sheet README (label, title, sheet list, conventions, then per former table: sheet name(s), what it holds, figures backed, source, column dictionary) then the data sheets.",
             "Numbers are copied verbatim from the shipped tables; nothing is re-derived here.",
             "Not in the package: the GTRD transcription-factor survey and network (former Supplementary Data 5 and 6) — they ship with the code release (103, reference_tables/).",
             "",
             "`rows` = rows summed over all data sheets of the workbook (README sheet excluded); `cols` = columns per data sheet, in sheet order.",
             "",
             "| Supplementary Data | sheets | rows (all data sheets) | cols | bytes | md5 | title |", "|---|---|---|---|---|---|---|",
             sprintf("| %d | %s | %s | %s | %s | %s | %s |", wb_sum$number, wb_sum$sheets, format(wb_sum$rows, big.mark = ","), wb_sum$cols, format(wb_sum$bytes, big.mark = ","), wb_sum$md5, wb_sum$title),
             "",
             "## Sheet -> former file",
             "", "| Supplementary Data | sheet | former file | former Supplementary Data | rows | cols |", "|---|---|---|---|---|---|",
             sprintf("| %d | %s | %s | %d | %s | %d |", index$number, index$sheet, index$internal_id, index$former_supplementary_data, format(index$rows, big.mark = ","), index$cols)),
           file.path(BUILD, "MANIFEST.md"))

# ---- retire the superseded single-file workbooks (12-27) -------------------------
# regenerated artifacts
KEEP_OLD <- file.path(BUILD, "_pre_consolidation_20260921")
stale <- setdiff(list.files(OUT, pattern = "^Supplementary_Data_[0-9]+\\.xlsx$"), sd_file_n(SD_WB$workbook))
if (length(stale)) {
  stopifnot("the pre-consolidation copy is missing — refusing to delete the superseded workbooks" = dir.exists(KEEP_OLD) && all(file.exists(file.path(KEEP_OLD, stale))))
  unlink(file.path(OUT, stale)); say("removed %d superseded workbooks from Supplementary_Data/ (%s); copies in _build/%s", length(stale), paste(stale, collapse = ", "), basename(KEEP_OLD))
}

say("\nWrote %d Supplementary Data workbooks (%d data sheets) to %s", nrow(SD_WB), nrow(index), OUT)
print(index[, .(number, sheet, internal_id, former_supplementary_data, rows, cols)])
# the delivery folder must hold exactly the 11 workbooks (+ _build/)
.left <- setdiff(list.files(OUT, all.files = FALSE), c(sd_file_n(SD_WB$workbook), "_build"))
stopifnot("Supplementary_Data/ holds files other than the 11 xlsx + _build/" = length(.left) == 0)
say("Supplementary_Data/ holds exactly %d xlsx + _build/ (builder aids: %s)", nrow(SD_WB), paste(list.files(BUILD), collapse = ", "))
cat("\n== sessionInfo ==\n"); print(sessionInfo())
cat("\n=== DONE ===\n", file = stderr())
