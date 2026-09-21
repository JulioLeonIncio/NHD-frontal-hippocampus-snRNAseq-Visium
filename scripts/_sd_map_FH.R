# =============================================================================
# _sd_map_FH.R — The fixed map of the "Supplementary Data" package:
# 11 themed multi-sheet workbooks built from the former single-table ST files.
# -----------------------------------------------------------------------------
# Shared (sys.source'd) by
#   * 100_package_supplementary_data_FH.R  — builds the workbooks from this map,
#   * _ST0_data_dictionary_FH.R and 90_build_supplementary_tables_FH.R — every "Supplementary Data N"
#     citation typed into a dictionary cell or ST6 scoring_context is rendered by sd_ref(), never typed,
#     so the csv and the xlsx cannot disagree on a number,
#   * 101_build_supplementary_information_FH.R (SI contents) and 103_code_release_FH.R (the two
#     GTRD tables that left the package and ship as reference tables).
# Pure data + pure functions; no side effects; base R only (the dictionary is sourced into a bare env).
#
# History.
# CHANGING this map RENUMBERS every citation in the manuscript — the structure is asserted at the bottom.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

# ---- workbooks: number -> short name + one-sentence referee-facing title ------------------------
SD_WB <- data.frame(stringsAsFactors = FALSE,
  workbook = 1:11,
  name = c("Data dictionary",
           "Samples, nuclei and annotation",
           "Per-nucleus differential expression",
           "Pseudobulk differential expression",
           "Gene-set enrichment",
           "Gene sets and markers",
           "Program effect sizes and sensitivity",
           "Microglia: transcription-factor activity and myeloid pooling",
           "Cross-cohort contrast with Zhou et al. 2023",
           "Sequencing-depth diagnostics",
           "Visium spatial transcriptomics"),
  title = c(
    "Data dictionary: the columns of every sheet of every Supplementary Data workbook, with the workbook, sheet, figures backed and source script of each table, and the analysis-wide conventions.",
    "Samples, nuclei and annotation: per-library sequencing and quality-control metrics, post-QC nucleus counts per cell type, region and condition, frontal cortical neuron subclass counts, hippocampal neuron subtype assignment per nucleus, and the myeloid nucleus selection chain.",
    "Per-nucleus MAST differential expression, NHD versus control, per cell type and region and per frontal neuronal subclass.",
    "Pseudobulk DESeq2 differential expression per cell type and region (confirmatory tier; hippocampal P values withheld).",
    "Gene-set enrichment (fgsea) on genes ranked by trimmed log2 fold change, per cell type and region, with nucleus-count power tiers, and the detection-gated microglial enrichment behind Supplementary Fig. 3b.",
    "Gene sets and markers: the curated signature gene sets used for module scoring and program dot plots, with symbol-mapping provenance, and the top ten protein-coding positive markers per annotated cell type.",
    "Program effect sizes and sensitivity: program-level Cliff's delta per cell type and region for Figs. 2d, 3e, 4d and 5d, per-nucleus neuronal program scores, the hippocampal unresolved-cluster sensitivity, and the per-library RNA-integrity arm.",
    "Microglia: myeloid transcription-factor activity on the CollecTRI regulon behind Fig. 2h, with the admission gate and the down-sampling sensitivity arm, and the pooled Micro-PVM versus microglia-only comparison for eight activation genes.",
    "Cross-cohort contrast with Zhou et al. 2023: direction consistency per cell type and region against the empirical chance baseline, and the same statistic recomputed per Zhou NHD donor arm with and without the shared donor.",
    "Sequencing-depth diagnostics: per-compartment depth ratios, UMI-floor sweep, depth-only artifact by normalization, neuronal-program depth artifact and myeloid depth comparison (six sheets, A-F).",
    "Visium spatial transcriptomics: the spatial-domain annotation of every spot with the layer IV order statistics, program expression per domain and the per-gene band means of the Fig. 6d heatmap with the gene sets scored on tissue, the microglia-rich versus microglia-poor spot contrasts, and the cell2location composition of every domain."))

# ---- sheets: one row per data sheet, in shipping order --------------------------------------------
# key   = the sheet's id (the former table id; sub-sheets of a former multi-sheet table carry its id + suffix)
# st    = the former table id the sheet belongs to (its README block, dictionary and asserts)
SD_SHEETS <- data.frame(stringsAsFactors = FALSE, rbind(
  c("ST0",   "ST0",   1,  "Data",                      1),
  c("ST12",  "ST12",  2,  "QC_per_sample",             13),
  c("ST9",   "ST9",   2,  "Nucleus_counts",            10),
  c("ST10",  "ST10",  2,  "Frontal_subclasses",        11),
  c("ST11",  "ST11",  2,  "Hippocampal_subtypes",      12),
  c("ST17b", "ST17b", 2,  "Myeloid_selection",         20),
  c("ST1",   "ST1",   3,  "Data",                      2),
  c("ST2",   "ST2",   4,  "Data",                      3),
  c("ST3",   "ST3",   5,  "Data",                      4),
  c("ST3b",  "ST3",   5,  "Micro_PVM_detection_gated", 4),
  c("ST6",   "ST6",   6,  "Signatures",                7),
  c("ST8",   "ST8",   6,  "Cell_type_markers",         9),
  c("ST24",  "ST24",  7,  "Program_deltas",            27),
  c("ST24b", "ST24",  7,  "By_program",                27),
  c("ST7",   "ST7",   7,  "Neuron_scores_per_nucleus", 8),
  c("ST16",  "ST16",  7,  "Unresolved_sensitivity",    18),
  c("ST15",  "ST15",  7,  "RIN_per_library",           16),
  c("ST15b", "ST15b", 7,  "RIN_gap",                   17),
  c("ST14",  "ST14",  8,  "TF_activity_Fig2h",         15),
  c("ST17",  "ST17",  8,  "Pooled_vs_microglia_only",  19),
  c("ST13",  "ST13",  9,  "Direction_consistency",     14),
  c("ST19",  "ST19",  9,  "Donor_arms",                22),
  c("ST19b", "ST19",  9,  "Within_Zhou_donors",        22),
  c("ST18A", "ST18",  10, "A_relative_RNA_content",    21),
  c("ST18B", "ST18",  10, "B_threshold_sweep",         21),
  c("ST18C", "ST18",  10, "C_artifact_by_treatment",   21),
  c("ST18D", "ST18",  10, "D_condition_balance",       21),
  c("ST18E", "ST18",  10, "E_neuron_module_artifact",  21),
  c("ST18F", "ST18",  10, "F_myeloid_condition_depth", 21),
  c("ST22",  "ST22",  11, "Spot_annotation",           25),
  c("ST22b", "ST22",  11, "Layer_IV_order",            25),
  c("ST23",  "ST23",  11, "Programs_per_domain",       26),
  c("ST23b", "ST23",  11, "Genes_by_band",             26),
  c("ST23c", "ST23",  11, "Gene_sets_Visium",          26),
  c("ST20",  "ST20",  11, "Microglia_proximity",       23),
  c("ST21",  "ST21",  11, "Composition_per_domain",    24)))
names(SD_SHEETS) <- c("key", "st", "workbook", "sheet", "old_number")
SD_SHEETS$workbook   <- as.integer(SD_SHEETS$workbook)
SD_SHEETS$old_number <- as.integer(SD_SHEETS$old_number)

# ---- the former tables that left the package (ship in the code release, reference_tables/) --------
SD_DROPPED <- data.frame(stringsAsFactors = FALSE, st = c("ST4", "ST5"), old_number = c(5L, 6L),
                         where = "code release, reference_tables/ (not a Supplementary Data file)")

# ---- the pre-consolidation order (old_number = position), still needed for the rename map -----------
SD_OLD_ORDER <- c("ST0","ST1","ST2","ST3","ST4","ST5","ST6","ST7","ST8","ST9","ST10","ST11",
                  "ST12","ST13","ST14","ST15","ST15b","ST16","ST17","ST17b","ST18","ST19","ST20","ST21","ST22","ST23","ST24")

# ---- referee-facing titles of the former tables (the "what this sheet holds" line of every README block) ----
SD_TITLE_PUBLIC <- c(
  ST0   = "Data dictionary: the columns of every sheet of every Supplementary Data workbook and the analysis-wide conventions.",
  ST1   = "Per-nucleus MAST differential expression, NHD versus control, per cell type and region and per frontal neuronal subclass.",
  ST2   = "Pseudobulk DESeq2 differential expression per cell type and region (confirmatory tier; hippocampal P values withheld).",
  ST3   = "Gene-set enrichment (fgsea) on genes ranked by trimmed log2 fold change, per cell type and region, with nucleus-count power tiers.",
  ST4   = "Genome-wide transcription-factor activity survey (decoupleR ULM, GTRD regulon) across six cell types and two regions.",
  ST5   = "The GTRD transcription-factor-to-target network used for the genome-wide survey (reference_tables/ST4).",
  ST6   = "Curated signature gene sets used for module scoring and program dot plots, with symbol-mapping provenance.",
  ST7   = "Per-nucleus neuronal program scores for excitatory and inhibitory neurons.",
  ST8   = "Cell-type marker genes: top ten protein-coding positive markers per annotated cell type.",
  ST9   = "Post-QC nucleus counts per cell type, region and condition.",
  ST10  = "Frontal cortical neuron subclass counts per condition.",
  ST11  = "Hippocampal neuron subtype assignment per nucleus.",
  ST12  = "Per-library sequencing and quality-control metrics for the seven snRNA-seq libraries.",
  ST13  = "Direction consistency with the Zhou et al. 2023 cohort per cell type and region, against the empirical chance baseline.",
  ST14  = "Myeloid transcription-factor activity on the CollecTRI regulon behind Fig. 2h, with the admission gate and the down-sampling sensitivity arm.",
  ST15  = "Per-library mean excitatory-neuron program scores with RNA integrity number.",
  ST15b = "Per-region RNA-integrity gap beside the NHD-minus-control program loss, on library means.",
  ST16  = "Hippocampal excitatory-neuron program effect sizes with and without the unresolved cluster.",
  ST17  = "Pooled Micro-PVM versus microglia-only effect sizes for eight activation genes.",
  ST17b = "Myeloid nucleus selection chain per condition and region.",
  ST18  = "Sequencing-depth diagnostics: per-compartment depth ratios, UMI-floor sweep, depth-only artifact by normalization, neuronal-program depth artifact and myeloid depth comparison.",
  ST19  = "Cross-cohort direction consistency with the Zhou et al. 2023 cohort recomputed with and without the shared donor (donor 861 = Zhou NHD3): all three Zhou NHD donors, the two independent donors, donor 861 alone, and each independent donor alone; same pipeline and statistic as sheet Direction_consistency.",
  ST20  = "Visium spots ranked by cell2location microglial fraction within each section and spatial domain: Cliff's delta between the microglia-rich and microglia-poor quartiles for microglial, astrocytic and lipid module scores and cell-type abundances, with block-bootstrap confidence intervals and a toroidal-shift null; grey matter in both conditions, white matter in the two NHD sections; absolute-abundance and six-neighbor versions as sensitivity arms.",
  ST21  = "Cell2location composition of every Visium spatial domain: mean per-spot proportion of the eight reference cell types per domain and condition, and the log2 NHD-to-control ratio of the means with a block-bootstrap interval; no P values (one donor per condition).",
  ST22  = "Visium spatial-domain annotation of every spot: final domain, tissue-call tier, and the layer I / pia refinement flags (BANKSY-primary rule and the histology-only sensitivity rule), with the per-section layer IV signal-and-order statistics as a second sheet.",
  ST23  = "Visium program expression per spatial domain: mean and median UMI per 10,000 of the complement/MHC-II and fatty-acid/sphingomyelin programs per domain and condition (Fig. 6e), the per-gene means and log2 NHD/control ratios per band of the Fig. 6d heatmap, and the gene membership of every program scored on tissue.",
  ST24  = "Program-level Cliff's delta, NHD versus control, per cell type and region: the curated-program and signature scores of Figs. 2d, 3e, 4d and 5d (microglia, astrocytes, oligodendrocyte lineage, excitatory and inhibitory neurons), per nucleus Wilcoxon q and the effect gate, long form and by program with both regions side by side.")

# ---- lookups --------------------------------------------------------------------------------------
sd_workbook_of <- function(st) SD_SHEETS$workbook[match(st, SD_SHEETS$st)]            # former table id -> workbook number
sd_sheets_of   <- function(st) SD_SHEETS$sheet[SD_SHEETS$st == st]                    # former table id -> its sheet(s), in order
sd_primary_sheet <- function(st) SD_SHEETS$sheet[match(st, SD_SHEETS$st)]             # first sheet of the former table
sd_single_sheet_wb <- function() { n <- table(SD_SHEETS$workbook); as.integer(names(n)[n == 1L]) }   # workbooks with one data sheet: cite without a sheet
.sd_list_and <- function(x) if (length(x) <= 1) x else paste0(paste(head(x, -1), collapse = ", "), " and ", tail(x, 1))

# sd_ref(keys): render a citation to one or more sheets. keys = sheet keys ("ST13", "ST19b") or former table ids
# ("ST18" = all six depth sheets). Sheets of the same workbook are grouped; a workbook is cited without a sheet when
# it has a single data sheet or when the cited former table owns every sheet of that workbook (ST3 -> workbook 5,
# ST18 -> workbook 10); different workbooks are separated by "; ".
#   sd_ref("ST13")            -> "Supplementary Data 9, sheet Direction_consistency"
#   sd_ref("ST1")             -> "Supplementary Data 3"          (single data sheet)
#   sd_ref("ST3")             -> "Supplementary Data 5"          (ST3 owns both sheets of workbook 5)
#   sd_ref(c("ST20","ST23"))  -> "Supplementary Data 11, sheets Programs_per_domain and Microglia_proximity"
#   sd_ref("ST24", sheet = FALSE) -> "Supplementary Data 7"
sd_ref <- function(keys, sheet = TRUE) {
  res <- lapply(keys, function(k) {
    i <- match(k, SD_SHEETS$key)
    if (is.na(i)) { i <- match(k, SD_SHEETS$st); if (is.na(i)) stop("sd_ref: unknown sheet key / table id: ", k) }
    wb <- SD_SHEETS$workbook[i]; st <- SD_SHEETS$st[i]
    owns_wb <- (k == st) && all(SD_SHEETS$st[SD_SHEETS$workbook == wb] == st)   # the table id was cited and it owns the whole workbook
    list(wb = wb, sheet = SD_SHEETS$sheet[i], whole = owns_wb)
  })
  wb <- vapply(res, `[[`, integer(1), "wb"); sh <- vapply(res, `[[`, character(1), "sheet"); whole <- vapply(res, `[[`, logical(1), "whole")
  single <- sd_single_sheet_wb()
  out <- character(0)
  for (w in unique(wb)) {                       # unique() keeps first-appearance order = the caller's order
    s <- unique(sh[wb == w])
    s <- s[order(match(s, SD_SHEETS$sheet[SD_SHEETS$workbook == w]))]   # shipping order within the workbook
    out <- c(out, if (!sheet || w %in% single || any(whole[wb == w])) sprintf("Supplementary Data %d", w)
                  else sprintf("Supplementary Data %d, %s %s", w, if (length(s) > 1) "sheets" else "sheet", .sd_list_and(s)))
  }
  paste(out, collapse = "; ")
}

# ---- old label -> new (the rename map; also the packaging-time safety net) -------------------------
# Every label form the draft / dictionary ever used: "Supplementary Data N" (1-27), "Supplementary Table N"
# (0-19, 15b, 17b: the earlier labels) and the bare internal id "STn".
sd_rename_map <- function() {
  one <- function(st) {
    if (st %in% SD_DROPPED$st) return(c(new_label = sprintf("DROPPED — %s", SD_DROPPED$where[SD_DROPPED$st == st]), new_sheet = ""))
    c(new_label = sprintf("Supplementary Data %d", sd_workbook_of(st)), new_sheet = sd_primary_sheet(st))
  }
  rows <- list()
  for (st in SD_OLD_ORDER) {
    nn <- one(st); old_n <- match(st, SD_OLD_ORDER)
    rows[[length(rows) + 1]] <- data.frame(old_label = sprintf("Supplementary Data %d", old_n), old_form = "Supplementary Data N",
                                           internal_id = st, new_label = nn[["new_label"]], new_sheet = nn[["new_sheet"]], stringsAsFactors = FALSE)
  }
  old_tab <- setdiff(SD_OLD_ORDER, c("ST20","ST21","ST22","ST23","ST24"))   # tables that never had a "Supplementary Table" label
  for (st in old_tab) { nn <- one(st)
    rows[[length(rows) + 1]] <- data.frame(old_label = paste0("Supplementary Table ", sub("^ST", "", st)), old_form = "Supplementary Table N",
                                           internal_id = st, new_label = nn[["new_label"]], new_sheet = nn[["new_sheet"]], stringsAsFactors = FALSE) }
  for (st in SD_OLD_ORDER) { nn <- one(st)
    rows[[length(rows) + 1]] <- data.frame(old_label = st, old_form = "STn",
                                           internal_id = st, new_label = nn[["new_label"]], new_sheet = nn[["new_sheet"]], stringsAsFactors = FALSE) }
  do.call(rbind, rows)
}

# sd_rekey_old(text): rewrite the unambiguous pre-consolidation citations inside free text —
# "Supplementary Data N" with N in 12-27 (singles, lists "23, 26", ranges "16-18"/"16–18", "and") and every
# "Supplementary Table N" form — to the new workbook/sheet citation. Numbers 1-11 are not touched: they exist in
# both numberings, so a bare "Supplementary Data 3" cannot be told apart; the source
# text must render those through sd_ref(). Returns the text with attribute "n_rekeyed" = citations rewritten.
sd_rekey_old <- function(text) {
  n <- 0L
  old_to_st <- function(num) { num <- sub("b$", "", num); SD_OLD_ORDER[as.integer(num)] }
  expand <- function(spec) {          # "16-18, 27" -> c(16,17,18,27) as OLD numbers
    parts <- trimws(unlist(strsplit(gsub("\\s+and\\s+", ",", spec), ",")))
    unlist(lapply(parts, function(p) { r <- as.integer(trimws(strsplit(p, "[-–]")[[1]])); if (length(r) == 2) seq(r[1], r[2]) else r }))
  }
  # a citation of another paper's supplement ("Zhou et al. 2023 Supplementary Table 1") is not ours: skip when the
  # preceding text ends in an author-year cue
  external <- function(k, start) grepl("(et al\\.|Zhou|\\b(19|20)[0-9]{2})[,;:]?\\s*$", substr(text[k], max(1L, start - 40L), start - 1L), perl = TRUE)
  render <- function(sts) {
    keep <- setdiff(sts, SD_DROPPED$st)
    if (!length(keep)) return("the GTRD reference tables of the code release")
    sd_ref(keep)
  }
  # (1) "Supplementary Data <list/range>" where every number is >= 12
  pat <- "Supplementary Data (\\d+(?:\\s*[-–]\\s*\\d+)?(?:,\\s*\\d+(?:\\s*[-–]\\s*\\d+)?)*(?:,?\\s+and\\s+\\d+)?)\\b"
  m <- gregexpr(pat, text, perl = TRUE)
  for (k in seq_along(text)) {
    if (is.na(text[k]) || m[[k]][1] < 0) next   # NA cells (empty csv fields) pass through
    starts <- as.integer(m[[k]]); lens <- attr(m[[k]], "match.length")
    for (j in rev(seq_along(starts))) {          # right-to-left so earlier offsets stay valid
      tok  <- substr(text[k], starts[j], starts[j] + lens[j] - 1L)
      nums <- expand(sub("^Supplementary Data ", "", tok))
      if (!all(nums >= 12L & nums <= 27L) || external(k, starts[j])) next
      new  <- render(SD_OLD_ORDER[nums]); n <- n + 1L
      text[k] <- paste0(substr(text[k], 1L, starts[j] - 1L), new, substr(text[k], starts[j] + lens[j], nchar(text[k])))
    }
  }
  # (2) "Supplementary Table N[b]" — an earlier form, unambiguous in any numbering
  pat2 <- "Supplementary Tables? (\\d+b?)\\b"
  m2 <- gregexpr(pat2, text, perl = TRUE)
  for (k in seq_along(text)) {
    if (is.na(text[k]) || m2[[k]][1] < 0) next
    starts <- as.integer(m2[[k]]); lens <- attr(m2[[k]], "match.length")
    for (j in rev(seq_along(starts))) {
      tok <- substr(text[k], starts[j], starts[j] + lens[j] - 1L)
      st  <- paste0("ST", sub("^Supplementary Tables? ", "", tok))
      if (!st %in% SD_OLD_ORDER || external(k, starts[j])) next
      new <- render(st); n <- n + 1L
      text[k] <- paste0(substr(text[k], 1L, starts[j] - 1L), new, substr(text[k], starts[j] + lens[j], nchar(text[k])))
    }
  }
  attr(text, "n_rekeyed") <- n
  text
}

# ---- structure asserts (the map is load-bearing for every citation in the manuscript) ---------------
stopifnot(nrow(SD_WB) == 11L, identical(SD_WB$workbook, 1:11), !anyDuplicated(SD_WB$name),
          all(nchar(SD_SHEETS$sheet) <= 31L),                                       # Excel sheet-name limit
          !grepl("[][*?:/\\\\]", paste(SD_SHEETS$sheet, collapse = "")),           # characters Excel forbids in sheet names
          !anyDuplicated(paste(SD_SHEETS$workbook, SD_SHEETS$sheet)),              # unique within a workbook
          !anyDuplicated(SD_SHEETS$key),
          setequal(unique(SD_SHEETS$workbook), 1:11),
          setequal(c(unique(SD_SHEETS$st), SD_DROPPED$st), SD_OLD_ORDER),         # every former table is placed or dropped
          setequal(names(SD_TITLE_PUBLIC), SD_OLD_ORDER),
          sd_workbook_of("ST0") == 1L, sd_workbook_of("ST12") == 2L, sd_workbook_of("ST1") == 3L, sd_workbook_of("ST2") == 4L,
          sd_workbook_of("ST3") == 5L, sd_workbook_of("ST6") == 6L, sd_workbook_of("ST24") == 7L, sd_workbook_of("ST14") == 8L,
          sd_workbook_of("ST13") == 9L, sd_workbook_of("ST18") == 10L, sd_workbook_of("ST22") == 11L,
          identical(sd_single_sheet_wb(), c(1L, 3L, 4L)),
          identical(sd_ref("ST13"), "Supplementary Data 9, sheet Direction_consistency"),
          identical(sd_ref("ST1"), "Supplementary Data 3"), identical(sd_ref("ST3"), "Supplementary Data 5"), identical(sd_ref("ST18"), "Supplementary Data 10"),
          identical(sd_ref("ST3b"), "Supplementary Data 5, sheet Micro_PVM_detection_gated"),
          identical(sd_ref(c("ST15","ST15b","ST16")), "Supplementary Data 7, sheets Unresolved_sensitivity, RIN_per_library and RIN_gap"),
          identical(as.vector(sd_rekey_old("see Supplementary Data 27 and Supplementary Data 23, 26")),
                    "see Supplementary Data 7, sheet Program_deltas and Supplementary Data 11, sheets Programs_per_domain and Microglia_proximity"),
          identical(as.vector(sd_rekey_old("Supplementary Data 3 stays")), "Supplementary Data 3 stays"),
          identical(as.vector(sd_rekey_old(c(NA, "x"))), c(NA, "x")),
          identical(as.vector(sd_rekey_old("course from Zhou et al. 2023 Supplementary Table 1. 7 libraries")), "course from Zhou et al. 2023 Supplementary Table 1. 7 libraries"),
          identical(as.vector(sd_rekey_old("see Supplementary Table 17b")), "see Supplementary Data 2, sheet Myeloid_selection"),
          nrow(sd_rename_map()) == 27L + 22L + 27L)

# ---- ST6 scoring_context vocabulary (shared by 90, which writes the column, and the dictionary, which lists it) ----
SD_ST6_CTX <- c(nucleus = sprintf("scored per nucleus (%s)", sd_ref("ST7")),
                fig2d   = sprintf("Fig. 2d program scoring (%s)", sd_ref("ST24")),
                visium  = sprintf("Visium program scores (%s)", sd_ref(c("ST23", "ST20"))),
                astro   = sprintf("Fig. 3e,f signature scores (%s)", sd_ref("ST24")),
                fig4d   = sprintf("Fig. 4d program scoring (%s)", sd_ref("ST24")),
                fig5d   = sprintf("Fig. 5d program scoring (%s)", sd_ref(c("ST24", "ST15", "ST15b", "ST16"))),
                fig2g   = "Fig. 2g microglial-state enrichment set",
                refonly = "reference only; not scored in a shipped panel")
