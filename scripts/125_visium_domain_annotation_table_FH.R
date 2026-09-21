#!/usr/bin/env Rscript
# =============================================================================
# 125_visium_domain_annotation_table_FH.R — Supplementary Data 11, sheets Spot_annotation + Layer_IV_order (internal ST22): the per-spot Visium
# spatial-domain annotation of the final integrated object, plus the per-section layer IV signal-and-order
# statistics (script 123) as a second sheet.
# -----------------------------------------------------------------------------
# What. One row per spot of ../Visium/integrated_harmony/NHD_frontal_integrated_harmony.rds (the object after the
# layer I / pia refinement of script 120): the final cluster and its domain (integrated_cluster_identity.csv through
# _visium_domains_FH.R), the domain's order and compartment (integrated_domain_levels.csv), the tissue-call tier the
# spot entered the analysis under (117b: `spaceranger` = Space Ranger in_tissue, `he_mask` = the H&E mask tier), and
# the refinement flags so every domain-dependent number can be recomputed under the sensitivity rule:
#   moved_to_L1pia   the BANKSY-primary rule (120, rule "banksy"): spot re-assigned to L1/pia  (== obj$histology_refined_L1)
#   previous_domain  the domain the spot carried before refinement (equals `domain` when not moved)
#   banksy_pial      the BANKSY-pial membership itself (121/122; the primary rule's expression term)
#   interior_WM_L6   the histology-only sensitivity rule: WM/L6-labelled spot strictly inside the pathologists' closed
#                    layer-I band (119b)
#   d_L1_line_um / d_WM_line_um   distance to the nearest pathologists' layer-I / white-matter line (119)
#   inside_L1_band   inside the closed layer-I band (119b)
# Nothing is re-derived: every column is read from the object or from the tables the refinement scripts wrote.
#
# What it is not. Not a re-run of the refinement (120) and not a new domain call; the numbers asserted below
# (18,765 spots; 1,039 moved under the primary rule; 151 under the histology-only rule) are the ones the Methods quote.
#
# Inputs:  ../Visium/integrated_harmony/NHD_frontal_integrated_harmony.rds          (final object; meta.data only is used)
#          ../Visium/integrated_harmony/integrated_cluster_identity.csv + integrated_domain_levels.csv  (via _visium_domains_FH.R)
#          ../Visium/<section>_masked/outs/tissue_call.csv                          (117b; column tissue_call_source)
#          tables/visium_L1_refinement_20260918.csv                                  (120)
#          tables/visium_sulcus_check_perspot_FH.csv                                 (119; d_L1, d_WM)
#          tables/visium_L1band_membership_20260918.csv                              (119b; inside_L1_band)
#          tables/visium_L4IT_signal_vs_order_FH.csv                                 (123; copied unchanged)
# Outputs: supplementary_tables/ST22_visium_domain_annotation_FH.csv        (sheet Spot_annotation of Supplementary Data 11)
#          supplementary_tables/ST22b_visium_L4IT_signal_vs_order_FH.csv    (sheet Layer_IV_order; script 100 ships it
#                                                                             as the second sheet, never as its own file)
#          logs/125_visium_domain_annotation_table_FH.log                   (provenance: counts + sessionInfo)
#
# Column DICTIONARY (ST22 sheet Data; the referee-facing wording lives in scripts/100_* as ST22_COLS):
#   spot_id             <section>_<barcode>; the column name of the integrated object
#   section             Visium section (CON_Frontal1 | CON_Frontal2 | NHD_Frontal1 | NHD_Frontal2)
#   condition           CON | NHD (the section prefix; asserted equal to obj$condition)
#   seurat_cluster      final integrated cluster id (0-10; 10 = the refinement cluster)
#   domain              spatial domain (integrated_cluster_identity.csv)
#   domain_order        superficial -> deep order of the domain (integrated_domain_levels.csv)
#   compartment         grey | white | other
#   tissue_call_source  spaceranger | he_mask (117b tiers; `excluded` never enters the object — asserted)
#   moved_to_L1pia      logical; re-assigned to L1/pia by the primary (BANKSY) rule
#   previous_domain     domain before refinement
#   banksy_pial         logical; member of a pial BANKSY cluster within the positional constraint (the primary rule's term)
#   interior_WM_L6      logical; the histology-only sensitivity rule
#   d_L1_line_um        distance (um) to the nearest pathologists' layer-I line
#   d_WM_line_um        distance (um) to the nearest pathologists' white-matter line
#   inside_L1_band      logical; inside the closed layer-I band
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))   # DOM_MAP / DOM_LEV / DOM_INFO / DOM_COMP

ROOT <- dirname(PROJ); IH <- file.path(ROOT, "Visium", "integrated_harmony"); TAB <- file.path(PROJ, "tables")
OUTD <- file.path(PROJ, "supplementary_tables"); LOGD <- file.path(PROJ, "logs")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE); dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
OBJ  <- file.path(IH, "NHD_frontal_integrated_harmony.rds")
F_REF <- file.path(TAB, "visium_L1_refinement_20260918.csv")
F_SUL <- file.path(TAB, "visium_sulcus_check_perspot_FH.csv")
F_MEM <- file.path(TAB, "visium_L1band_membership_20260918.csv")
F_L4  <- file.path(TAB, "visium_L4IT_signal_vs_order_FH.csv")
SECS  <- c("CON_Frontal1", "CON_Frontal2", "NHD_Frontal1", "NHD_Frontal2")
F_TC  <- setNames(file.path(ROOT, "Visium", paste0(SECS, "_masked"), "outs", "tissue_call.csv"), SECS)
OUT_A <- file.path(OUTD, "ST22_visium_domain_annotation_FH.csv")
OUT_B <- file.path(OUTD, "ST22b_visium_L4IT_signal_vs_order_FH.csv")
LOG   <- file.path(LOGD, "125_visium_domain_annotation_table_FH.log")
# The Methods numbers this table must reproduce (rule banksy)
N_SPOTS <- 18765L; N_MOVED <- 1039L; N_INTERIOR <- 151L

for (f in c(OBJ, F_REF, F_SUL, F_MEM, F_L4, F_TC))
  if (!file.exists(f)) stop("MISSING ", f, " — upstream step (117b / 119 / 119b / 120 / 123) has not run")
# The run-consistency guard is the content check below —
# moved_primary must equal obj$histology_refined_L1 spot by spot — which is stronger than a timestamp.

cat("== read object meta ==\n")
obj <- readRDS(OBJ); md <- obj@meta.data; md$spot_id <- rownames(md); rm(obj); invisible(gc())
stopifnot("object must carry histology_refined_L1 (script 120)" = "histology_refined_L1" %in% names(md),
          is.logical(md$histology_refined_L1), all(md$sample_id %in% SECS))
cat("spots:", nrow(md), "\n")

cat("== domains ==\n")
d <- data.table(spot_id = md$spot_id, section = as.character(md$sample_id))
d[, condition := sub("_.*$", "", section)]
stopifnot(all(d$condition %in% c("CON", "NHD")), "obj$condition disagrees with the section prefix" = all(d$condition == as.character(md$condition)))
d[, seurat_cluster := as.integer(as.character(md$seurat_clusters))]
d[, domain := unname(DOM_MAP[as.character(seurat_cluster)])]
stopifnot("a cluster id has no domain in integrated_cluster_identity.csv" = !anyNA(d$domain), all(d$domain %in% DOM_LEV))
d[, domain_order := DOM_INFO$order[match(domain, DOM_INFO$domain)]]
d[, compartment  := unname(DOM_COMP[domain])]
stopifnot(!anyNA(d$domain_order), all(d$compartment %in% c("grey", "white", "other")))

cat("== tissue-call tier (117b) ==\n")
tc <- rbindlist(lapply(SECS, function(s) {
  x <- fread(F_TC[[s]], colClasses = list(character = "barcode"))
  stopifnot("tissue_call.csv lacks barcode / tissue_call_source" = all(c("barcode", "tissue_call_source") %in% names(x)))
  data.table(spot_id = paste0(s, "_", x$barcode), tissue_call_source = x$tissue_call_source) }))
stopifnot(!anyDuplicated(tc$spot_id), "spots missing from tissue_call.csv" = all(d$spot_id %in% tc$spot_id))
d[, tissue_call_source := tc$tissue_call_source[match(spot_id, tc$spot_id)]]
stopifnot("an analysed spot carries tissue call `excluded`" = !any(d$tissue_call_source == "excluded"),
          all(d$tissue_call_source %in% c("spaceranger", "he_mask")))

cat("== refinement flags (120) ==\n")
ref <- fread(F_REF)
stopifnot(all(c("spot_id", "previous_domain", "moved_primary", "banksy_pial", "interior_WM_L6") %in% names(ref)),
          !anyDuplicated(ref$spot_id), "spots missing from visium_L1_refinement_20260918.csv" = all(d$spot_id %in% ref$spot_id))
i <- match(d$spot_id, ref$spot_id)
d[, moved_to_L1pia  := md$histology_refined_L1]
d[, previous_domain := ref$previous_domain[i]]
d[, banksy_pial     := as.logical(ref$banksy_pial[i])]
d[, interior_WM_L6  := as.logical(ref$interior_WM_L6[i])]
stopifnot("obj$histology_refined_L1 differs from moved_primary in the refinement table" = all(d$moved_to_L1pia == as.logical(ref$moved_primary[i])),
          "a spot moved to L1/pia is not in the L1/pia domain" = all(d[moved_to_L1pia == TRUE]$domain == "L1/pia"),
          "previous_domain must equal domain for every spot that was not moved" = all(d[moved_to_L1pia == FALSE, previous_domain == domain]),
          all(d$previous_domain %in% DOM_LEV), !anyNA(d$banksy_pial), !anyNA(d$interior_WM_L6))

cat("== line distances (119) and band membership (119b) ==\n")
sul <- fread(F_SUL); stopifnot(all(c("spot_id", "d_L1", "d_WM") %in% names(sul)), !anyDuplicated(sul$spot_id),
                               "spots missing from visium_sulcus_check_perspot_FH.csv" = all(d$spot_id %in% sul$spot_id))
j <- match(d$spot_id, sul$spot_id)
d[, d_L1_line_um := round(sul$d_L1[j], 1)]; d[, d_WM_line_um := round(sul$d_WM[j], 1)]
mem <- fread(F_MEM); stopifnot(all(c("spot_id", "inside_L1_band") %in% names(mem)), !anyDuplicated(mem$spot_id),
                               "spots missing from visium_L1band_membership_20260918.csv" = all(d$spot_id %in% mem$spot_id))
d[, inside_L1_band := as.logical(mem$inside_L1_band[match(spot_id, mem$spot_id)])]
stopifnot(!anyNA(d$d_L1_line_um), !anyNA(d$d_WM_line_um), !anyNA(d$inside_L1_band), all(d$d_L1_line_um >= 0), all(d$d_WM_line_um >= 0))

# ---- the numbers the Methods quote ----------------------------------------
stopifnot("row count differs from the object" = nrow(d) == nrow(md),
          "expected 18,765 spots in the final object" = nrow(d) == N_SPOTS,
          "expected 1,039 spots moved to L1/pia (primary rule)" = sum(d$moved_to_L1pia) == N_MOVED,
          "expected 151 spots under the histology-only rule" = sum(d$interior_WM_L6) == N_INTERIOR,
          !anyDuplicated(d$spot_id))
setcolorder(d, c("spot_id", "section", "condition", "seurat_cluster", "domain", "domain_order", "compartment", "tissue_call_source",
                 "moved_to_L1pia", "previous_domain", "banksy_pial", "interior_WM_L6", "d_L1_line_um", "d_WM_line_um", "inside_L1_band"))
stopifnot(ncol(d) == 15)
setorder(d, section, spot_id)   # deterministic row order

cat("== sheet Layer_IV_order (123, copied unchanged) ==\n")
l4 <- fread(F_L4)
L4_COLS <- c("section", "n_cortical", "mean_L4IT_per10k", "mean_RORB_per10k", "frac_L4high", "neighbour_L4high_given_L4high",
             "coherence_ratio", "morans_I", "median_depth_um_L4high", "n_L4_domain")
stopifnot("visium_L4IT_signal_vs_order_FH.csv columns differ from the expected 10" = identical(names(l4), L4_COLS),
          nrow(l4) == 4, setequal(l4$section, SECS))

# ---- write (temp then rename; the shipped folder is never left half-written) ----
.wr <- function(x, p) { tmp <- paste0(p, ".tmp"); fwrite(x, tmp); invisible(file.rename(tmp, p)) }
.wr(d, OUT_A); .wr(l4, OUT_B)
# copy check: the shipped ST22b must be byte-equivalent in content to the 123 table
stopifnot(isTRUE(all.equal(as.data.frame(fread(OUT_B)), as.data.frame(l4), check.attributes = FALSE)))

# ---- summary + provenance -------------------------------------------------
summ <- c(
  sprintf("ST22 Data: %d spots x %d columns -> %s", nrow(d), ncol(d), OUT_A),
  sprintf("ST22b Layer_IV_order: %d sections x %d columns -> %s", nrow(l4), ncol(l4), OUT_B),
  "spots per section:", capture.output(print(table(d$section))),
  "domain x condition:", capture.output(print(table(d$domain, d$condition)[DOM_LEV, ])),
  "tissue_call_source x section:", capture.output(print(table(d$tissue_call_source, d$section))),
  sprintf("moved_to_L1pia = %d (primary, BANKSY); banksy_pial = %d; interior_WM_L6 = %d (histology-only); inside_L1_band = %d",
          sum(d$moved_to_L1pia), sum(d$banksy_pial), sum(d$interior_WM_L6), sum(d$inside_L1_band)),
  "previous_domain of moved spots:", capture.output(print(table(d[moved_to_L1pia == TRUE]$previous_domain, d[moved_to_L1pia == TRUE]$section))),
  sprintf("d_L1_line_um range %.1f-%.1f; d_WM_line_um range %.1f-%.1f", min(d$d_L1_line_um), max(d$d_L1_line_um), min(d$d_WM_line_um), max(d$d_WM_line_um)))
cat(summ, sep = "\n")
writeLines(c(sprintf("125_visium_domain_annotation_table_FH.R  run %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             sprintf("object: %s (mtime %s)", OBJ, format(file.mtime(OBJ))),
             sprintf("inputs: %s", paste(basename(c(F_REF, F_SUL, F_MEM, F_L4)), collapse = ", ")),
             "", summ, "", capture.output(sessionInfo())), LOG)
cat("log:", LOG, "\n=== DONE ===\n")
