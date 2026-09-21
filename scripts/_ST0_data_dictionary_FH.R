#!/usr/bin/env Rscript
# =============================================================================
# _ST0_data_dictionary_FH.R — _ST0_data_dictionary_FH.R Authored content for ST0, the supplementary-table data dictionary.
# Sourced (sys.source) by 90_build_supplementary_tables_FH.R, which then CHECKS
# this dictionary against the columns actually present in every shipped CSV and
# fails loudly on any column that is documented but absent, or present but
# undocumented.  That check is the whole point: a dictionary that can silently
# drift from the data is worse than none.
#
# Defines exactly two objects:
#   ST0_TABLES      data.frame(table, title, backs_figures, source_script)
#   ST0_COLUMNS     data.frame(table, column, type, permitted_values, definition)
#   ST0_CONVENTIONS data.frame(convention, definition)
#
# House rules honoured here:
#   * The Zhou 2023 comparison is direction consistency against an empirical
#     baseline; the word "replication" is never used for it.
#   * Our regions are Frontal and Hippo only.  "OCC / occipital" appears solely
#     inside Zhou's own panel labels (an external cohort), never as our tissue.
#   * American spelling.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

# 100 carry the same numbers and sheet names. The map is sourced into this environment (bare env under 90's sys.source).
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
local({
  .root <- if (exists("PROJ", inherits = TRUE)) get("PROJ", inherits = TRUE) else {
    .c <- Sys.getenv("NHD_PROJ")
    .c[dir.exists(.c)][1] }
  .m <- file.path(.root, "scripts", "_sd_map_FH.R")
  stopifnot("MISSING scripts/_sd_map_FH.R (the fixed Supplementary Data map)" = file.exists(.m))
  sys.source(.m, envir = parent.env(environment()))
})
stopifnot(exists("sd_ref"), exists("SD_ST6_CTX"))

ST0_TABLES <- data.frame(stringsAsFactors = FALSE, rbind(
c("Table1","Donor and sample characteristics of the Frontal + Hippocampus cohort.",
  "All figures (cohort description)","hand-curated"),
c("ST0","Data dictionary: the columns of every sheet of every Supplementary Data workbook (the tables built after the dictionary — sequencing-depth, Zhou donor arms, Visium and program effect sizes — are documented in full in the README sheet of their workbook) and the analysis-wide conventions.",
  "n/a","90_build_supplementary_tables_FH.R + _ST0_data_dictionary_FH.R"),
c("ST1","Per-nucleus MAST discovery differential expression, NHD vs control, for each cell type x region and for each frontal neuronal subclass (excitatory subclasses from the Jorstad 2023 dorsolateral prefrontal cortex label transfer, including L4 IT; inhibitory subclasses from Azimuth).",
  "Fig 1c-f; Fig 2b,e,f; Fig 3 astrocyte volcano; Fig 4 oligodendrocyte volcano; Fig 5 neuron volcano grid",
  "40_mast_percell_dual_FH.R (broad) + 45_mast_subclass_FH.R (subclass)"),
c("ST2","Pseudobulk DESeq2 differential expression per cell type x region, the confirmatory tier of the dual-method design.",
  "Fig 1 dual-method concordance; confirmatory tier for Figs 2-5",
  "30_pseudobulk_wald_FH.R"),
c("ST3","fGSEA pathway enrichment on the pre-ranked MAST gene list, per subtype x region, with nucleus-count power tiers.",
  "Fig 2g; Fig 3d; Supplementary Fig 3b",
  "44_gsea_MAST_cliffs_FH.R + 44b_gsea_subclass_FH.R"),
c("ST4","Genome-wide transcription-factor activity survey (decoupleR ULM on the GTRD regulon) across all six cell types x two regions. Not the source of Fig. 2h.",
  "None directly; broader survey supporting the TF section","02_TF_activity_FH.R"),
c("ST5","The GTRD transcription-factor -> target network used for ST4. Unsigned: mode_of_regulation is 1 for every edge. Not the network behind Fig. 2h.",
  "None directly; the regulon underlying ST4","02_TF_activity_FH.R"),
c("ST6","Curated signature gene sets used for module scoring and program dot plots.",
  "Fig 2d; Fig 3e,f; Fig 4d; Fig 5d","90_build_supplementary_tables_FH.R"),
c("ST7","Per-nucleus neuronal program module scores (8 programs) for excitatory and inhibitory neurons.",
  "Fig 5 neuron program violins and dot matrix","90_build_supplementary_tables_FH.R"),
c("ST8","Cell-type marker genes: top 10 protein-coding positive markers per annotated cell type.",
  "Fig 1 atlas / annotation support","90_build_supplementary_tables_FH.R"),
c("ST9","Post-QC nucleus counts per cell type x region x condition. Totals 55,354.",
  "Fig 1 composition; Supplementary composition panels","90_build_supplementary_tables_FH.R"),
c("ST10","Frontal cortical neuron subclass counts per condition: 15 subclasses (9 excitatory from the Jorstad 2023 dorsolateral prefrontal cortex label transfer, including L4 IT; 6 inhibitory from Azimuth).",
  "Fig 5 neuron subclass composition; Fig 5f",
  "90_build_supplementary_tables_FH.R"),
c("ST11","Hippocampal neuron subtype assignment per nucleus (de-novo re-annotation).",
  "Fig 5 hippocampal neuron panels","90_build_supplementary_tables_FH.R"),
c("ST12","Per-library sequencing and QC metrics for the seven snRNA-seq libraries.",
  "Methods; Supplementary QC figure","90_build_supplementary_tables_FH.R"),
c("ST13","Cross-cohort direction consistency against Zhou 2023, per cell type x region panel, measured against an empirical baseline. Not a replication test.",
  "Supplementary Fig. 3c (microglia); Fig 3c (astrocytes); Fig 5c (neurons); Supplementary Fig. 2e (all twelve comparisons)",
  "70_zhou_crosscohort_FH.R"),
c("ST14","Myeloid transcription-factor activity on the signed, curated CollecTRI regulon, scored on log-normalized counts across all 2,985 nuclei, with the pre-declared six-criterion admission gate and a 750-UMI down-sampling falsification test alongside. This is the table behind Fig. 2h.",
  "Fig 2h","68_micro_TF_collectri_FH.R + 70b_micro_TF_bipartite_FH.R"),
c("ST15","Per-library mean neuronal program scores in excitatory neurons, with RNA integrity number, for the seven libraries. Computed from the same per-nucleus scores that produced Fig. 5d.",
  "Fig 5d (per-library support); RNA-quality sensitivity",
  "90_build_supplementary_tables_FH.R (score cache of 4C_signature_violins_FH.R)"),
c("ST15b","Per-region summary of ST15: the control-minus-NHD RIN gap beside the NHD-minus-control program loss computed on library means.",
  "Fig 5d (RNA-quality sensitivity)","90_build_supplementary_tables_FH.R"),
c("ST16","Cliff's delta per neuronal program in hippocampal excitatory neurons, computed twice: over all nuclei, and excluding the nuclei whose ST11 subtype is Unresolved.",
  "Fig 5d (taxonomy sensitivity)","90_build_supplementary_tables_FH.R"),
c("ST17","Gene-level comparison of the pooled Micro-PVM analysis against a microglia-only re-analysis, for eight microglial activation genes x two regions. Shows that pooling CD163+/F13A1+ nuclei with microglia does not cost gene-level sensitivity.",
  "Fig 2 (pooling justification); Methods annotation paragraph",
  "92_micro_only_MAST_FH.R + 90_build_supplementary_tables_FH.R"),
c("ST17b","Myeloid nucleus selection chain per condition x region: Azimuth Micro-PVM, oligodendrocyte doublets removed, CD163+/F13A1+ removed, microglia analysed; with the separate Lymphocyte class counted alongside.",
  "Fig 2a/2d/2h; Supplementary Fig. 3","92_micro_only_MAST_FH.R")
))
names(ST0_TABLES) <- c("table", "title", "backs_figures", "source_script")

# -----------------------------------------------------------------------------
# ST14 family sizes quoted in the dictionary are read from the shipped table, never typed:
# the hand-typed "100 tests / 54 factors / 46 gate rows" went stale when the primary arm
# switched to the log-normalized all-nuclei run (118 rows / 64 factors). The BH family is
# asserted to be the whole table (p.adjust over every row reproduces q_BH).
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a   # NA-aware (base's is NULL-only)
.st14 <- local({
  .root <- if (exists("PROJ", inherits = TRUE)) get("PROJ", inherits = TRUE) else {
    .c <- Sys.getenv("NHD_PROJ")
    .c[dir.exists(.c)][1] }
  f <- file.path(.root, "supplementary_tables", "ST14_micro_TF_activity_CollecTRI_Fig2h.csv")
  stopifnot("MISSING ST14 — build ST14 before the dictionary" = file.exists(f))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  stopifnot(all(c("TF","p_value","q_BH","score_750","scored_both_regions") %in% names(d)),
            "ST14 q_BH is not BH over the whole table (family size differs from nrow)" =
              max(abs(stats::p.adjust(d$p_value, method = "BH") - d$q_BH)) < 1e-8)
  # the 750-UMI arm is BH-corrected within its own table (tables/micro_TF_collectri_thin750_FH.csv), whose
  # family can exceed the rows that join ST14; the family size is read from that table, not from the join
  f750 <- file.path(.root, "tables", "micro_TF_collectri_thin750_FH.csv")
  stopifnot("MISSING tables/micro_TF_collectri_thin750_FH.csv (the q_BH_750 family) — run 68 with DEPTH_MODE=thin750" = file.exists(f750))
  s750 <- utils::read.csv(f750, stringsAsFactors = FALSE)
  stopifnot(all(c("source","region","p_value") %in% names(s750)))
  s750$q <- stats::p.adjust(s750$p_value, method = "BH")
  j <- merge(d[!is.na(d$score_750), c("TF","region","q_BH_750")], s750[, c("source","region","q")], by.x = c("TF","region"), by.y = c("source","region"))
  stopifnot("ST14 q_BH_750 is not BH over the 750-UMI arm's own table" = nrow(j) == sum(!is.na(d$score_750)) && max(abs(j$q_BH_750 - j$q)) < 1e-8)
  # The 750-UMI tests that join no row of ST14 are named, with the
  # reason, and both are derived here so the sentence cannot go stale. Reason, verified at build
  # time: run_ulm scores a factor in a region only if >= MINSIZE (10) of its CollecTRI targets are
  # among the genes detected in >= 5% of that region's microglia (68_micro_TF_collectri_FH.R, gkeep).
  # At native depth each orphan factor has MINSIZE - 1 such targets in frontal microglia (it is
  # scored in the other region); the 750-UMI floor removes the shallow nuclei, per-gene detection
  # among the retained nuclei rises, and a further target clears the gene floor, so the arm can
  # score a test the primary estimator could not. This is the regulon-coverage floor, not the 10%
  # TF-detection floor (which is region-invariant by construction and admits all four factors).
  .orph <- s750[!paste(s750$source, s750$region) %in% paste(d$TF, d$region), c("source", "region")]
  .orph <- .orph[order(match(.orph$region, c("Frontal", "Hippo")), .orph$source), , drop = FALSE]
  orphan_sentence <- ""
  if (nrow(.orph)) {
    fnet <- file.path(.root, "data", "TF_regulon_collectri.rds")
    fgd  <- file.path(.root, "tables", "micro_gene_cliffsdelta_FH.csv")
    f68  <- file.path(.root, "scripts", "68_micro_TF_collectri_FH.R")
    stopifnot("MISSING data/TF_regulon_collectri.rds (needed to explain the unjoined 750-UMI tests)" = file.exists(fnet),
              "MISSING tables/micro_gene_cliffsdelta_FH.csv (primary per-gene table of 68)" = file.exists(fgd),
              "MISSING scripts/68_micro_TF_collectri_FH.R" = file.exists(f68))
    .l68 <- readLines(f68, warn = FALSE)
    .minsize <- suppressWarnings(as.integer(sub("^MINSIZE\\s*<-\\s*(\\d+)L?.*$", "\\1", grep("^MINSIZE\\s*<-", .l68, value = TRUE)[1])))
    .gfloor  <- suppressWarnings(as.numeric(sub("^.*gkeep <- detf >= ([0-9.]+).*$", "\\1", grep("gkeep <- detf >= ", .l68, value = TRUE)[1])))
    stopifnot("could not read MINSIZE / the per-gene detection floor from script 68" = is.finite(.minsize) && is.finite(.gfloor),
              "ST14 carries a factor below the MINSIZE read from script 68 — the reason text would be wrong" = min(d$n_target) >= .minsize)
    .net <- readRDS(fnet); .gd <- utils::read.csv(fgd, stringsAsFactors = FALSE)
    stopifnot(all(c("source", "target") %in% names(.net)), all(c("gene", "region") %in% names(.gd)))
    .orph$n_target_primary <- mapply(function(tf, rg) sum(.net$target[.net$source == tf] %in% .gd$gene[.gd$region == rg]),
                                     .orph$source, .orph$region)
    .orph$in_primary_other_region <- .orph$source %in% d$TF
    stopifnot("an unjoined 750-UMI test is NOT explained by the regulon-coverage floor (n_target >= MINSIZE at native depth) — re-derive the reason" =
                all(.orph$n_target_primary < .minsize),
              "an unjoined 750-UMI factor is absent from ST14 altogether — the 'scored in the other region' wording would be false" =
                all(.orph$in_primary_other_region))
    .reg_word <- c(Frontal = "frontal", Hippo = "hippocampal")
    .n_tgt <- sort(unique(.orph$n_target_primary))
    orphan_sentence <- sprintf(paste(
      "The %s 750-UMI tests without a primary row are %s; each of these factors is scored in the other region, but in %s microglia at native depth only %s of its CollecTRI targets are among the genes detected in >= %d%% of nuclei, short of the %d-target regulon-coverage floor of the primary estimator (the 10%% factor-detection floor is region-invariant and admits all of them);",
      "after the 750-UMI floor removes the shallow nuclei, per-gene detection among the retained nuclei rises and a further target clears the gene floor. Their sensitivity-arm scores are therefore not shipped in this table."),
      c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine")[nrow(.orph)] %||% as.character(nrow(.orph)),
      { .lab <- sprintf("%s (%s)", .orph$source, .orph$region); if (length(.lab) > 1) paste0(paste(head(.lab, -1), collapse = ", "), " and ", tail(.lab, 1)) else .lab },
      paste(unique(unname(.reg_word[.orph$region])), collapse = " or "),
      if (length(.n_tgt) == 1) as.character(.n_tgt) else sprintf("%d-%d", min(.n_tgt), max(.n_tgt)),
      as.integer(round(100 * .gfloor)), .minsize)
  }
  list(n_rows = nrow(d), n_tf = length(unique(d$TF)), n_gate = length(unique(d$TF[d$scored_both_regions %in% TRUE])),
       n_750 = nrow(s750), n_750_joined = sum(!is.na(d$score_750)), n_na750 = sum(is.na(d$score_750)),
       n_orphan = nrow(.orph), orphan = .orph, orphan_sentence = orphan_sentence)
})

.col <- function(tab, ...) {
  a <- list(...)
  do.call(rbind, lapply(a, function(x)
    data.frame(table = tab, column = x[1], type = x[2],
               permitted_values = x[3], definition = x[4],
               stringsAsFactors = FALSE)))
}

ST0_COLUMNS <- rbind(

# ---- Table 1 ----------------------------------------------------------------
.col("Table1",
 c("Characteristic","character","free text","Name of the donor or sample characteristic reported in that row."),
 c("NHD donor","character","free text","Value of the characteristic for the Nasu-Hakola disease donor."),
 c("Control donor","character","free text","Value of the characteristic for the control donor: female, 55, no neurological diagnosis, minimal Alzheimer-type change (Braak I / Thal 0 / CERAD 0), post-mortem interval 6 h.")),

# ---- ST1 --------------------------------------------------------------------
.col("ST1",
 c("comparison","character","<cell_type>_<region>, e.g. Astro_Frontal, L5_IT_Frontal","Identifier of the NHD vs control test; one MAST model was fitted per comparison."),
 c("cell_type","character","Astro, Micro-PVM, Neuron_Ex, Neuron_Inh, Oligo, OPC, or a frontal neuronal subclass (L2/3 IT, L4 IT, L5 IT, L6 IT, L6 IT Car3, L5/6 NP, L6 CT, L6b, Lamp5, Pvalb, Sncg, Sst, Vip)","Cell population in which the test was run. Excitatory subclasses are the Jorstad 2023 dorsolateral prefrontal cortex labels transferred per nucleus (Seurat TransferData, within-area subclass); inhibitory subclasses are the Azimuth human motor-cortex labels. L5 ET (2 control nuclei) and Sst Chodl (0 NHD nuclei) were not tested."),
 c("region","character","Frontal | Hippo","Brain region. Frontal pools the grey- and white-matter dissections; Hippo is hippocampus."),
 c("gene","character","HGNC symbol","Tested gene."),
 c("avg_log2FC","numeric","log2 units","Seurat average log2 fold change, NHD relative to control. Positive = higher in NHD."),
 c("p_val","numeric","0-1; floored at 2.2250739e-308","Raw two-sided p-value from the MAST hurdle model."),
 c("p_val_adj","numeric","0-1; floored at 2.2250739e-308","Bonferroni-adjusted p-value as returned by Seurat, corrected over all features in the assay (26,973), not over the features tested in this comparison. Conservative by roughly 7-14x; used to rank, not as a calibrated error rate."),
 c("pct_NHD","numeric","proportion 0-1","Fraction of NHD nuclei with at least one raw RNA UMI for the gene. Computed on raw RNA counts, never on SCT."),
 c("pct_CON","numeric","proportion 0-1","Fraction of control nuclei with at least one raw RNA UMI for the gene. Computed on raw RNA counts, never on SCT."),
 c("cliffs_delta","numeric","-1 to 1","Cliff's delta effect size on the raw RNA LogNormalize layer. Positive = higher in NHD."),
 c("n_NHD","integer","count of nuclei","Number of NHD nuclei in the comparison."),
 c("n_CON","integer","count of nuclei","Number of control nuclei in the comparison."),
 c("discovery","logical","TRUE | FALSE","Discovery feature flag: p_val_adj < 0.05 and |cliffs_delta| >= 0.15 and max(pct_NHD, pct_CON) >= 0.10. Unfiltered (before the artifact filter)."),
 c("artifact","logical","TRUE | FALSE","TRUE if the gene is removed by the curated ambient/artifact filter (unannotated lncRNA and antisense classes, clustered histones, housekeeping genes, sex-chromosome genes, unnamed clone contigs, and neuronal transcripts scored in a glial cell type)."),
 c("discovery_filtered","logical","TRUE | FALSE","discovery and not artifact. This is the set the figures plot: Fig 1c-f reproduce exactly from it (e.g. Micro-PVM 599, not the 654 given by `discovery` alone), as do the Fig 1e/1f correlations and their gene counts."),
 c("pb_confirmed","logical","TRUE | FALSE","Frontal only. The gene also clears the lane-level pseudobulk gate (DESeq2 padj < 0.01, |apeglm log2FC| > log2(2.5)) with the same sign as the per-nucleus estimate. Frontal is the only region with two libraries per condition. Note that those two libraries are the grey- and white-matter dissections of one donor, so this is a lane-consistency check, not biological replication."),
 c("xregion_supported","logical","TRUE | FALSE","Hippocampus only. The per-nucleus direction agrees with the frontal pseudobulk sign for the same gene and cell type. The hippocampus has one NHD library and no pseudobulk tier of its own, so this is cross-region direction agreement, a different and weaker quantity than pb_confirmed."),
 c("support_tier","character","lane_reproducible | xregion_supported | discovery_only | NA","Support tier for discovery_filtered features; NA otherwise. Fig 1c plots discovery_filtered as the pale bar and the non-discovery_only tiers as the solid bar. The solid bar therefore means different things in the two regions and the legend says so."),
 c("latent_vars","character","matter+log10nCount | hippo_area+log10nCount | log10nCount","The covariates in that row's MAST model, recorded per row. log10(nCount_RNA) is present in every model: NHD frontal nuclei carry roughly twice the post-QC median UMIs of control, so sequencing depth runs with diagnosis and is modelled rather than ignored. `matter` (grey/white) blocks the frontal dissection; `hippo_area` blocks the two control hippocampal captures.")),

# ---- ST2 --------------------------------------------------------------------
.col("ST2",
 c("gene","character","HGNC symbol","Tested gene."),
 c("baseMean","numeric","normalized counts","DESeq2 mean of normalized counts across libraries."),
 c("log2FoldChange","numeric","log2 units","DESeq2 Wald log2 fold change, NHD relative to control (unshrunken)."),
 c("lfcSE","numeric","log2 units","Standard error of log2FoldChange."),
 c("stat","numeric","z","DESeq2 Wald statistic."),
 c("pvalue","numeric","0-1, or NA; floored at 2.2250739e-308","Raw Wald p-value. NA = the gene was not tested (all-zero counts or a Cook's-distance outlier)."),
 c("padj","numeric","0-1, or NA; floored at 2.2250739e-308","BH adjusted p-value. NA means not tested (DESeq2 independent filtering) and is never recoded to 0 or 1."),
 c("log2FoldChange_apeglm","numeric","log2 units","apeglm shrunken log2 fold change. Reported as a bonus, not as a significance gate."),
 c("lfcSE_apeglm","numeric","log2 units","Posterior SD of the apeglm shrunken estimate."),
 c("comparison","character","<cell_type>_<region>","Cell type x region in which the pseudobulk model was fitted."),
 c("n_cells_total","integer","count of nuclei","Nuclei aggregated into the pseudobulk profiles for this comparison."),
 c("n_NHD_libraries_physical","integer","1 or 2","Physical NHD libraries contributing: 2 in frontal cortex (grey + white), 1 in hippocampus."),
 c("n_CON_libraries_physical","integer","2","Physical control libraries contributing: frontal grey + white, or the two hippocampal captures."),
 c("design_n_lanes","integer","count","Columns in the DESeq2 design matrix. For the hippocampus this is 4 because the single NHD library was randomly split into two pseudo-replicates; it is not a count of physical libraries."),
 c("design_n_NHD_lanes","integer","count","NHD columns in the design matrix. 2 in hippocampus refers to the two halves of one library."),
 c("design_n_CON_lanes","integer","count","Control columns in the design matrix."),
 c("dispersion_basis","character","within_donor_tissue_libraries | pseudo_split","How dispersion was estimated. within_donor_tissue_libraries = two physical libraries per condition, which in this design are the grey- and white-matter dissections of one donor (frontal only); dispersion is therefore estimated across tissue compartments of a single donor, not across donors. pseudo_split = the single NHD hippocampal library halved in silico; dispersion is not estimable at all and pvalue, padj and stat are NA by construction."),
 c("model_fit_reportable","logical","TRUE | FALSE","Whether the DESeq2 model fit for that row can be reported at all, i.e. whether dispersion had two physical libraries per condition to be estimated from. TRUE only for within_donor_tissue_libraries rows. It does not mean that a donor-level inference is licensed: with one donor per condition, the frontal pseudobulk tier is a descriptive compartment-consistency check — it asks whether an effect seen per nucleus also holds when the grey- and white-matter dissections of the same donor are aggregated — and it is used in ST1 only to set the support_tier label, never as evidence of replication. FALSE rows carry NA for pvalue, padj and stat and their fold changes are descriptive only."),
 c("note","character","free text or NA","Provenance note carried from the DESeq2 run; flags the pseudo-split rows."),
 c("region","character","Frontal | Hippo","Region, derived from `comparison`.")),

# ---- ST3 --------------------------------------------------------------------
.col("ST3",
 c("subtype","character","<cell_type>_<region> or <subclass>_<region>","Population whose ranked gene list was tested."),
 c("pathway","character","MSigDB set identifier","Gene-set name exactly as distributed by MSigDB. Set names containing REPLICATION (e.g. GOBP_DNA_REPLICATION) are database identifiers, not a claim of experimental replication."),
 c("pval","numeric","0-1, or NA; floored at 2.2250739e-308","fgsea permutation p-value. NA = fgsea could not estimate the enrichment (degenerate permutation null)."),
 c("padj","numeric","0-1, or NA; floored at 2.2250739e-308","BH adjusted p-value across pathways within the subtype x database. NA where pval is NA."),
 c("ES","numeric","-1 to 1","Enrichment score."),
 c("NES","numeric","unitless, or NA","Normalized enrichment score. Positive = enriched among genes higher in NHD. NA where the null was degenerate."),
 c("size","integer","count of genes","Number of set genes present in the ranked list after filtering."),
 c("leadingEdge","character","pipe-separated HGNC symbols","Leading-edge genes driving the enrichment."),
 c("direction","character","NHD | CON | empty","Sign of the enrichment: NHD = up in NHD, CON = up in control. Empty where NES is NA."),
 c("db","character","GO:BP | Reactome | Hallmark","MSigDB collection the pathway came from."),
 c("metric","character","log2fc_winsorized","The pre-ranking metric: winsorized signed log2 fold change. A signed-Cliff's-delta prerank is not used here."),
 c("n_nuclei_min","integer","count of nuclei","Smaller of the two condition nucleus counts for this subtype; drives power_tier."),
 c("power_tier","character","solid | low_conf | excluded","Nucleus-count power flag: solid >= 150, low_conf 50-149, excluded < 50.")),

# ---- ST4 --------------------------------------------------------------------
.col("ST4",
 c("statistic","character","ulm","Inference method: decoupleR univariate linear model."),
 c("source","character","HGNC symbol of a transcription factor","The transcription factor whose activity was inferred."),
 c("condition","character","avg_log2FC","The per-gene statistic that was scored; here the MAST average log2 fold change."),
 c("score","numeric","z-like, unitless","ULM activity score. Positive = the factor's GTRD targets move up in NHD. Because GTRD is unsigned, the sign reflects target direction only, not activation vs repression."),
 c("p_value","numeric","0-1","Two-sided p-value for the ULM score. No multiple-testing column is shipped for this survey; correct across the 5,184 rows if used inferentially."),
 c("comparison","character","<cell_type>_<region>, 12 levels","Cell type x region whose ranked statistic was scored.")),

# ---- ST5 --------------------------------------------------------------------
.col("ST5",
 c("source","character","HGNC symbol of a transcription factor","Regulator."),
 c("target","character","HGNC symbol","Target gene."),
 c("mode_of_regulation","numeric","1 for every edge","GTRD carries no sign, so this column is 1 throughout; the network is unsigned. Contrast CollecTRI (ST14), whose mode of regulation is +1 / -1."),
 c("weight","numeric","unitless","Edge weight passed to decoupleR.")),

# ---- ST6 --------------------------------------------------------------------
.col("ST6",
 c("signature","character","signature set name",paste0("Name of the curated gene set: the published glial/neuronal sets, the seven Fig. 2d microglial programs, the eleven Fig. 3e,f astrocyte sets (six defined for that panel under their panel labels; the five published ones under their curated names, with the panel label stated in scoring_context), the eight Fig. 4d oligodendrocyte-lineage programs, the Pan 2024 CSF1R states, and the Fig5d_* sets scored in Fig. 5d. Every program of ", sd_ref("ST24"), " is listed here.")),
 c("gene","character","HGNC symbol","Member gene. Rows are unique on (signature, gene)."),
 c("symbol_original","character","free text","The symbol as published in the source that defined the set, including mouse casing and legacy aliases. Equal to `gene` wherever the published symbol needed no mapping."),
 c("symbol_mapped","character","HGNC symbol, or empty","The current human HGNC symbol actually looked up in the expression matrix. Empty where no human symbol exists (mouse-only genes)."),
 c("in_scoring_matrix","logical","TRUE | FALSE","Whether symbol_mapped is a row of the expression matrix the signature was scored on (36,601 features)."),
 c("used_in_scoring","logical","TRUE | FALSE","Whether the gene actually contributed to a score. FALSE either because it is not in the matrix, or because the signature is reference-only and is never scored by any shipped panel - see scoring_context."),
 c("exclusion_reason","character","none | symbol_not_mapped | not_in_matrix | below_min_detection | alias_merged","Why the gene did not contribute. `none` = it did contribute. `alias_merged` = the published symbol is a legacy alias and was scored under symbol_mapped. `symbol_not_mapped` = no human ortholog symbol exists (mouse-only). `below_min_detection` is defined for completeness but is empty by construction: no scoring path in this project applies a per-gene detection floor."),
 c("scoring_context","character",paste(SD_ST6_CTX[c("nucleus","fig2d","visium","astro","fig4d","fig5d","fig2g","refonly")], collapse = " | "),
   paste0("Where the signature is scored: the eight neuronal programs and the four curated sets they are built from are scored per nucleus (", sd_ref("ST7"), "; the curated set names the column of that sheet it becomes); the seven microglial programs are the Fig. 2d programs (", sd_ref("ST24"), "), and the two of them that are also scored on tissue carry the Visium context as well (", sd_ref(c("ST23","ST20")), "); the eleven astrocyte sets are the Fig. 3e,f signature scores (", sd_ref("ST24"), ") - the six defined for that panel are listed under their panel labels, and each of the five published sets names the panel program it is reported as (for example A1_reactive as program Complement/IFN-reactive (Liddelow)); the eight oligodendrocyte-lineage programs are the Fig. 4d programs (", sd_ref("ST24"), "), and the five of them that are also scored on tissue carry the Visium context as well; the Fig5d_ sets are the Fig. 5d programs (", sd_ref(c("ST24","ST15","ST15b","ST16")), "); Zhou_NHD_microglia_up is a gene set of the Fig. 2g microglial-state enrichment. Every other set, including the Pan 2024 CSF1R states, is shipped for reference and is not scored in any panel."))),

# ---- ST7 --------------------------------------------------------------------
.col("ST7",
 c("cell","character","nucleus barcode","Unique nucleus identifier in the atlas."),
 c("Condition","character","NHD | CON","Donor condition."),
 c("Region","character","Frontal | Hippo","Brain region."),
 c("cell_type","character","Neuron_Ex | Neuron_Inh","Excitatory or inhibitory neuron."),
 c("orig.ident","character","TFHS000529-TFHS000535","Sequencing library of origin (7 libraries)."),
 c("Synaptic","numeric","mean log-normalized expression","Module score: union of the presynaptic and postsynaptic curated sets."),
 c("OxPhos","numeric","mean log-normalized expression","Module score: neuronal oxidative-phosphorylation set."),
 c("UPR","numeric","mean log-normalized expression","Module score: unfolded-protein response."),
 c("DNA_damage","numeric","mean log-normalized expression","Module score: DNA-damage response."),
 c("Apoptosis","numeric","mean log-normalized expression","Module score: apoptosis."),
 c("IEG","numeric","mean log-normalized expression","Module score: immediate-early genes."),
 c("Postsynaptic","numeric","mean log-normalized expression","Module score: postsynaptic set."),
 c("Presynaptic","numeric","mean log-normalized expression","Module score: presynaptic set.")),

# ---- ST8 --------------------------------------------------------------------
.col("ST8",
 c("gene","character","HGNC symbol","Marker gene."),
 c("cluster","character","cell-type label","Cell type for which the gene is a positive marker."),
 c("avg_log2FC","numeric","log2 units","Average log2 fold change of the cell type versus all other nuclei."),
 c("pct_in","numeric","proportion 0-1","Fraction of nuclei of that cell type expressing the gene.")),

# ---- ST9 --------------------------------------------------------------------
.col("ST9",
 c("cell_type","character","Astro, Endo, Lymphocyte, Micro-PVM, Micro-PVM_doublet, Neuron_Ex, Neuron_Inh, Oligo, OPC, Pericytes","Atlas cell-type annotation. Micro-PVM_doublet = myeloid nuclei co-expressing another lineage's markers at atlas annotation; retained in the post-QC total, excluded from every analysis."),
 c("Region","character","Frontal | Hippo","Brain region."),
 c("CON","integer","count of nuclei","Post-QC nuclei from the control donor."),
 c("NHD","integer","count of nuclei","Post-QC nuclei from the NHD donor.")),

# ---- ST10 -------------------------------------------------------------------
.col("ST10",
 c("neuron_subtype","character","L2/3 IT | L4 IT | L5 IT | L6 IT | L6 IT Car3 | L5 ET | L5/6 NP | L6 CT | L6b | Lamp5 | Pvalb | Sncg | Sst | Sst Chodl | Vip","Frontal neuronal subclass. Excitatory nuclei carry the label transferred from the Jorstad 2023 dorsolateral prefrontal cortex reference (which, unlike the motor-cortex reference, has an L4 IT class); inhibitory nuclei carry the Azimuth human motor-cortex subclass. The raw Azimuth call for every nucleus is retained in the deposited nucleus metadata as azimuth_subclass."),
 c("Region","character","Frontal","Frontal cortex only; hippocampal neurons are typed de novo in ST11."),
 c("Condition","character","NHD | CON","Donor condition."),
 c("n","integer","count of nuclei","Nuclei assigned to that subclass. Rows with n = 0 are omitted."),
 c("label_source","character","ctx_DFC_Jorstad2023_within | ctx_Azimuth_subclass","Reference that labeled the subclass: ctx_DFC_Jorstad2023_within = Jorstad 2023 dorsolateral prefrontal cortex within-area subclass transferred per nucleus (all nine excitatory subclasses); ctx_Azimuth_subclass = Azimuth human motor-cortex subclass (all six inhibitory subclasses).")),

# ---- ST11 -------------------------------------------------------------------
.col("ST11",
 c("barcode","character","nucleus barcode","Unique nucleus identifier in the atlas."),
 c("Region","character","Hippo","Hippocampus."),
 c("neuron_subtype","character","CA3 | Subiculum | Inh-MGE | Inh-CGE (VIP) | Inh-CGE (LAMP5) | Unresolved","De-novo hippocampal neuron subtype. Five subtypes are resolved; Unresolved is a cluster that failed the private-marker gate and is not a subtype (see taxonomy_flag). Dentate-gyrus granule neurons were not captured by this dissection, CA1/CA2 were not resolved, and no Cajal-Retzius population was recovered (TP73 is below the detection floor in every cluster)."),
 c("taxonomy_flag","character","private-marker supported | unresolved (failed private-marker gate)","Whether the label may be used as a subtype. Every call must be supported by at least one private marker of its own label at >= 10% detection; a group-mean argmax alone can be carried by one promiscuous gene. 'unresolved (failed private-marker gate)' marks the 190 nuclei of a low-complexity cluster that failed that gate: TP73, the private Cajal-Retzius marker, is detected in 1.1% of its nuclei, and the argmax call rested on RELN, which reaches 44-52% in other clusters; median nFeature 1,312 versus 2,543-5,458 for the resolved subtypes; 186 of the 190 nuclei come from one library. These nuclei are not a neuronal subtype and are excluded from every subtype-level comparison."),
 c("Condition","character","NHD | CON","Donor condition."),
 c("source_annot","character","hip_denovo_label","Provenance of the subtype call: the de-novo hippocampal clustering, not an external reference mapping.")),

# ---- ST12 -------------------------------------------------------------------
.col("ST12",
 c("Donor","character","NHD | Control | Total","Donor for the library; the final row is the cohort total."),
 c("Region","character","Frontal cortex | Hippocampus (capture 1) | Hippocampus (capture 2) | —","Long-form library label; — in the Total row. Frontal cortex corresponds to region token Frontal elsewhere; Hippocampus corresponds to Hippo."),
 c("Matter","character","Grey | White | —","Dissected matter for cortical libraries; — where not applicable (hippocampal libraries and the Total row)."),
 c("Reads","character","count, comma-grouped","Raw sequencing reads."),
 c("Valid barcodes","character","percent","Fraction of reads with a valid barcode."),
 c("Mean reads/nucleus","character","count, comma-grouped","Mean reads per called nucleus."),
 c("Sequencing saturation","character","percent","Cell Ranger sequencing saturation."),
 c("Reads mapped to transcriptome","character","percent","Fraction of reads confidently mapped to the transcriptome."),
 c("Median genes/nucleus","character","count, comma-grouped","Median genes detected per nucleus."),
 c("Median UMI/nucleus","character","count, comma-grouped","Median UMIs per nucleus."),
 c("Nuclei (Cell Ranger)","character","count, comma-grouped","Nuclei called by Cell Ranger before QC."),
 c("Nuclei (post-QC)","character","count, comma-grouped","Nuclei retained in the analysis atlas. The Total row sums to 55,354."),
 c("RIN","character","1-10, one decimal","RNA integrity number for that library (Agilent Bioanalyzer). Systematically lower in NHD (mean 5.6) than control (mean 7.2), and predictive of the transcriptome mapping rate (Pearson r = 0.96 across the seven libraries). Reported, not adjusted for: with one donor per condition it cannot be separated from condition.")),

# ---- ST13 -------------------------------------------------------------------
.col("ST13",
 c("panel","character","<our cell type> <our region> vs Zhou OCC","Panel identifier. OCC here is Zhou's occipital cohort, an external dataset; our own regions are only Frontal and Hippo."),
 c("n_shared","integer","count of genes","Genes tested in both cohorts for this panel; the threshold-free comparison uses all of them."),
 c("n_both_sig","integer","count of genes","Genes significant in both cohorts."),
 c("rho_full","numeric","-1 to 1","Spearman correlation of effect sizes across all n_shared genes (threshold-free)."),
 c("rho_full_p","numeric","0-1; floored at 2.2250739e-308","p-value for rho_full."),
 c("rho_sig","numeric","-1 to 1","Spearman correlation restricted to the n_both_sig genes."),
 c("rho_sig_p","numeric","0-1; floored at 2.2250739e-308","p-value for rho_sig."),
 c("tau_sig","numeric","-1 to 1","Kendall tau over the n_both_sig genes."),
 c("tau_sig_p","numeric","0-1; floored at 2.2250739e-308","p-value for tau_sig."),
 c("n_concord","integer","count of genes","Genes whose effect has the same sign in both cohorts, among n_both_sig."),
 c("concord_frac","numeric","proportion 0-1","Direction-consistency fraction, n_concord / n_both_sig. Must be read against baseline_concord, not against 0.5."),
 c("binom_p","numeric","0-1; floored at 2.2250739e-308","Binomial p-value for concord_frac."),
 c("rho_loo","numeric","-1 to 1","Leave-one-gene-out minimum of rho_full; a stability check that no single gene carries the correlation."),
 c("baseline_concord","numeric","proportion 0-1","Empirical direction-consistency baseline for this panel: the concordance expected by chance given the marginal sign imbalance of the two gene lists. This, not 0.5, is what concord_frac is compared with."),
 c("n_baseline_genes","integer","count of genes","Genes used to estimate baseline_concord."),
 c("n_dir_zero","integer","count of genes","Genes with a zero effect size in either cohort (no direction to compare)."),
 c("cell_type","character","Astro | Micro-PVM | Neuron_Ex | Neuron_Inh | Oligo | OPC","Our cell type in the panel."),
 c("region","character","Frontal | Hippo","Our region in the panel. Zhou's own tissue is occipital cortex.")),

# ---- ST14 -------------------------------------------------------------------
.col("ST14",
 c("TF","character","HGNC symbol of a transcription factor","Factor whose activity was inferred in microglia."),
 c("region","character","Frontal | Hippo","Region in which the factor was scored."),
 c("statistic","character","ulm","decoupleR univariate linear model."),
 c("ranking_metric","character","cliffs_delta","The per-gene statistic scored: microglial Cliff's delta, NHD vs control, on depth-matched nuclei."),
 c("score","numeric","z-like, unitless","ULM activity score on the signed CollecTRI regulon. Positive = inferred activity increased in NHD."),
 c("p_value","numeric","0-1; floored at 2.2250739e-308","Two-sided p-value for the ULM score."),
 c("q_BH","numeric","0-1; floored at 2.2250739e-308",sprintf("Benjamini-Hochberg q-value computed once over all %d TF x region tests in this table (%d unique factors). It is not computed over the %d factors of the gate audit.", .st14$n_rows, .st14$n_tf, .st14$n_gate)),
 c("score_750","numeric","ULM t statistic, or NA",sprintf("The same activity score recomputed on the 750-unique-molecular-identifier sensitivity arm (nuclei below the floor excluded, the remainder thinned to exactly that depth). NA where the factor does not clear the 10%% detection floor after thinning and so cannot be scored in that arm; %d of %d rows are NA for that reason, including PPARG, STAT1 and CEBPB. Down-sampling is a falsification test here, not the estimator - see Methods, Sequencing depth.", .st14$n_na750, .st14$n_rows)),
 c("q_BH_750","numeric","0-1, or NA",trimws(sprintf("Benjamini-Hochberg q for score_750, corrected within the sensitivity arm's own family of tests (%d TF x region tests scored after thinning, of which %d join a row of this table) rather than the primary family (%d), because the two arms score different numbers of factors. Compare the two columns for agreement, not for identical values. %s", .st14$n_750, .st14$n_750_joined, .st14$n_rows, .st14$orphan_sentence))),
 c("score_noHK","numeric","z-like, unitless","Housekeeping leave-out score: the same ULM score after removing housekeeping, heat-shock and nuclear-quality transcripts from the ranked list."),
 c("n_target","integer","count of genes","Target coverage: CollecTRI targets of this factor present in the scored gene list for this region."),
 c("scored_both_regions","logical","TRUE | FALSE","TRUE if the factor cleared the detection floor in both regions and therefore reached the admission gate. FALSE factors are scored in one region only and are ungated by construction."),
 c("gate_outcome","character","gated hit | excluded | declared control | declared null (DAM-2) | not gated","Outcome of the pre-declared six-criterion gate. 'declared control' is SPI1, admitted by declaration as a lineage control, not by passing; 'declared null (DAM-2)' marks the factor carried as the pre-declared DAM-2 null."),
 c("gate_first_failure","character","passes | C1 coverage | C2 reproducibility | C3 robustness | C4 coherence | C5 interpretability | C6 anchor support | not gated (scored in one region only)","First criterion the factor failed."),
 c("gate_mean_z","numeric","z-like, unitless","Mean of the frontal and hippocampal scores; the quantity plotted in Fig. 2h."),
 c("gate_coverage","integer","count of genes","Smaller of the two regional target counts."),
 c("gate_leaveout_ratio","numeric","ratio >= 0","Smaller of the two regional |score_noHK| / |score| ratios; how much of the score survives removing housekeeping transcripts."),
 c("gate_mor_purity","numeric","0.5-1","Fraction of the factor's CollecTRI edges sharing one sign; low purity makes the score uninterpretable."),
 c("gate_tf_own_delta","numeric","-1 to 1, or NA","The factor's OWN Cliff's delta as a gene; used to detect a score incoherent with the factor's own expression."),
 c("gate_anchor_detected","integer","count of genes, or NA","Declared literature anchor genes for this factor that were detected."),
 c("gate_anchor_total","integer","count of genes, or NA","Declared literature anchor genes for this factor."),
 c("gate_anchor_dir_frac","numeric","0-1, or NA","Fraction of detected anchors moving in the direction the inferred activity predicts."),
 c("gate_C1_coverage","logical","TRUE | FALSE | NA","C1: target coverage at or above the minimum."),
 c("gate_C2_reproducible","logical","TRUE | FALSE | NA","C2: both regional scores exceed the null band with the same sign."),
 c("gate_C3_robust","logical","TRUE | FALSE | NA","C3: the score survives the housekeeping leave-out in both regions."),
 c("gate_C4_coherent","logical","TRUE | FALSE | NA","C4: the score is not contradicted by the factor's own expression change."),
 c("gate_C5_interpretable","logical","TRUE | FALSE | NA","C5: the factor's regulon is sign-pure enough for the score to have a direction."),
 c("gate_C6_anchored","logical","TRUE | FALSE | NA","C6: enough declared literature anchors are detected and move in the predicted direction.")),

# ---- ST15 -------------------------------------------------------------------
.col("ST15",
 c("library_id","character","TFHS000529-TFHS000535","Sequencing library (10x lane) the row summarises."),
 c("Condition","character","CON | NHD","Donor condition of that library."),
 c("Region","character","Frontal | Hippo","Region the library was dissected from."),
 c("Library","character","Grey | White | capture 1 | capture 2","Dissection fraction (frontal grey or white matter) or which of the hippocampal captures the library is. The control hippocampus was captured twice; the NHD hippocampus once."),
 c("RIN","numeric","4.4-7.5","RNA integrity number for that library, joined from ST12; not re-measured here."),
 c("n_nuclei","integer","count","Excitatory neurons from that library contributing to the means."),
 c("Presynaptic","numeric","mean module score","Mean Fig. 5d Presynaptic score over that library's excitatory neurons. Gene set = ST6 signature Fig5d_Presynaptic, not ST6/ST7 Presynaptic."),
 c("OxPhos","numeric","mean module score","Mean Fig. 5d OxPhos score. Gene set = ST6 Fig5d_OxPhos (23 genes), not ST7 OxPhos (16)."),
 c("IEG","numeric","mean module score","Mean immediate-early-gene score. This set is identical to ST7's IEG and to IEG_core."),
 c("Postsynaptic","numeric","mean module score","Mean Fig. 5d Postsynaptic score. Gene set = ST6 Fig5d_Postsynaptic, not ST7 Postsynaptic.")),

# ---- ST15b ------------------------------------------------------------------
.col("ST15b",
 c("Region","character","Frontal | Hippo","Region summarised."),
 c("RIN_gap_CON_minus_NHD","numeric","RIN units","Mean control RIN minus mean NHD RIN across that region's libraries. Positive = control RNA is the better preserved."),
 c("Presynaptic_loss","numeric","module-score units","Mean NHD library score minus mean control library score. Negative = lower in NHD. Computed on library means, not nucleus means, so the deeper NHD libraries cannot dominate by nucleus count."),
 c("OxPhos_loss","numeric","module-score units","As Presynaptic_loss, for OxPhos."),
 c("IEG_loss","numeric","module-score units","As Presynaptic_loss, for IEG."),
 c("Postsynaptic_loss","numeric","module-score units","As Presynaptic_loss, for Postsynaptic.")),

# ---- ST16 -------------------------------------------------------------------
.col("ST16",
 c("set","character","all_nuclei | excluding_unresolved","Nucleus set. `excluding_unresolved` drops hippocampal excitatory nuclei whose ST11 neuron_subtype is Unresolved."),
 c("n_NHD","integer","count","NHD hippocampal excitatory nuclei in that set."),
 c("n_CON","integer","count","Control hippocampal excitatory nuclei in that set."),
 c("Presynaptic","numeric","-1 to 1","Cliff's delta, NHD vs control, for the Fig. 5d Presynaptic program. Same estimator as ST1's cliffs_delta. Negative = lower in NHD."),
 c("OxPhos","numeric","-1 to 1","Cliff's delta for the Fig. 5d OxPhos program."),
 c("IEG","numeric","-1 to 1","Cliff's delta for the IEG program."),
 c("Postsynaptic","numeric","-1 to 1","Cliff's delta for the Fig. 5d Postsynaptic program. This is the program that moves: the Unresolved exclusion is almost entirely an NHD exclusion, so the postsynaptic effect depends on those nuclei in a way the other three do not.")),

# ---- ST17 -------------------------------------------------------------------
.col("ST17",
 c("gene","character","HGNC symbol","One of the eight microglial activation genes compared between the pooled and microglia-only analyses."),
 c("region","character","Frontal | Hippo","Region of the contrast."),
 c("pooled_delta","numeric","-1 to 1","Cliff's delta from the pooled Micro-PVM analysis, copied verbatim from ST1."),
 c("pooled_pct_NHD","numeric","0-1","Fraction of NHD Micro-PVM nuclei with at least one raw UMI, from ST1."),
 c("pooled_pct_CON","numeric","0-1","Fraction of control Micro-PVM nuclei with at least one raw UMI, from ST1."),
 c("pooled_padj","numeric","0-1","Bonferroni-adjusted MAST p from ST1. Bonferroni, not BH: Seurat adjusts over all 26,973 assay features."),
 c("pooled_discovery","logical","TRUE | FALSE","ST1 discovery_filtered: padj < 0.05 AND |delta| >= 0.15 AND raw detection >= 10%, after the curated artifact filter."),
 c("micro_delta","numeric","-1 to 1","Cliff's delta from the microglia-only re-analysis on the 2,824 microglia, identical model."),
 c("micro_pct_NHD","numeric","0-1","Raw detection in NHD microglia."),
 c("micro_pct_CON","numeric","0-1","Raw detection in control microglia."),
 c("micro_padj","numeric","0-1 or NA","Bonferroni-adjusted MAST p from the microglia-only run. NA where MAST did not test the gene - see micro_tested."),
 c("micro_discovery","logical","TRUE | FALSE","Discovery call in the microglia-only run, same three criteria."),
 c("delta_difference","numeric","-2 to 2","pooled_delta minus micro_delta. The empirical size of the pooling effect for that gene; it needs no denominator."),
 c("pooled_tested","logical","TRUE | FALSE","Whether MAST returned the gene in the pooled run. FALSE means the gene fell below min.pct 0.1 or logfc.threshold 0.1 and was not tested - it does not mean no effect."),
 c("micro_tested","logical","TRUE | FALSE","As pooled_tested, for the microglia-only run. Where FALSE, delta and detection are still reported (computed directly on the raw RNA layer) so the row is never blank.")),

# ---- ST17b ------------------------------------------------------------------
.col("ST17b",
 c("Condition","character","CON | NHD","Donor condition."),
 c("Region","character","Frontal | Hippo","Region."),
 c("MicroPVM_azimuth","integer","count","Nuclei labelled Micro-PVM by Azimuth transfer, before any myeloid-specific filtering. Sums to 3,081."),
 c("minus_oligo_doublets","integer","count","Removed at myeloid re-clustering as oligodendrocyte doublets (module score > 0.5 over PLP1/MBP/MOG/MOBP/CNP/ST18/CLDN11/MAG). Sums to 95."),
 c("after_doublet_removal","integer","count","MicroPVM_azimuth minus minus_oligo_doublets. Sums to 2,986."),
 c("minus_CD163_F13A1","integer","count","Removed as the CD163+/F13A1+ cluster. Sums to 162."),
 c("microglia_analysed","integer","count","Microglia carried into Fig 2d/2h and into ST17's microglia-only arm. Sums to 2,824."),
 c("Lymphocyte_separate_class","integer","count","Nuclei in the Azimuth Lymphocyte class. These nuclei were never inside Micro-PVM: they are a separate class and are shown alongside only because Supplementary Fig. 3d re-embeds Micro-PVM and Lymphocyte together. Sums to 60, and they remain inside the 55,354 total.")),

# ---- ST0 (self) -------------------------------------------------------------
.col("ST0",
 c("table","character","Table1, ST0-ST24, or CONVENTIONS","Table the row documents (internal table id); CONVENTIONS rows document analysis-wide definitions."),
 c("supplementary_data_label","character","Supplementary Data N, sheet X (N = 1-11); Supplementary Data N for a single-sheet workbook; Table 1; CONVENTIONS; or the code-release location of a table that is not a Supplementary Data file","Where the table ships: the Supplementary Data workbook and the sheet inside it that holds the rows the row documents (the workbook and sheet columns carry the same in machine-readable form). Present in the Supplementary Data 1 workbook; the csv form of the dictionary carries the internal id only."),
 c("workbook","integer","1-11, or empty","Number of the Supplementary Data workbook that holds the table (empty for Table 1, for CONVENTIONS rows and for the two GTRD reference tables that ship with the code release). Present in the Supplementary Data 1 workbook only."),
 c("sheet","character","sheet name, or empty","Name of the sheet inside that workbook (its first sheet for a table that spans several sheets, whose further sheets are documented in the README sheet of the workbook). Present in the Supplementary Data 1 workbook only."),
 c("table_title","character","free text","Title of that table."),
 c("backs_figures","character","free text","Figure panels the table supplies numbers for."),
 c("source_script","character","script file name(s)","Script of the code release that produced the table."),
 c("column","character","column name, or convention name for CONVENTIONS rows","Column being defined."),
 c("type","character","character | numeric | integer | logical | mixed | convention","Storage type; mixed for the tables documented in full in the README sheet of their workbook."),
 c("permitted_values","character","free text","Units, range, or the closed set of permitted values."),
 c("definition","character","free text","What the column means."))
)

# -----------------------------------------------------------------------------
ST0_CONVENTIONS <- data.frame(stringsAsFactors = FALSE, rbind(
c("Cohort",
  "One NHD donor (donor 861; female, 59; Braak III / Thal 2 / CERAD 2, i.e. intermediate AD neuropathologic change; post-mortem interval 36 h) and one control donor with no neurological diagnosis (female, 55; Braak I / Thal 0 / CERAD 0; post-mortem interval 6 h). Two confounds run with diagnosis and are stated rather than adjusted for, as in Table 1: the AD co-pathology is in the NHD donor and not in the control, and the NHD donor has the longer post-mortem interval and lower RNA integrity (mean RIN 5.6 vs 7.2). Two further features of the NHD donor are stated as well: the NHD donor is one of the three NHD donors in Zhou et al. 2023 (their NHD3), and the NHD donor had epilepsy from age 47 with repeated status epilepticus, to which the reporting neuropathologist attributed her cerebellar neuronal loss. Seizure activity alters activity-dependent neuronal transcription and causes neuronal loss independently of NHD, so it cannot be separated from disease effect in one donor; it bears particularly on the immediate-early-gene program and on neuronal loss. Clinical course from Zhou et al. 2023 Supplementary Table 1. 7 snRNA-seq libraries (control 4, NHD 3), 55,354 post-QC nuclei: control 37,371 / NHD 17,983; frontal 37,034 / hippocampus 18,320."),
c("Region tokens",
  paste0("Our regions are only Frontal and Hippo. Frontal = frontal cortex, pooling the grey- and white-matter dissections. Hippo = hippocampus (dentate-gyrus granule neurons were not captured; CA1/CA2 were not resolved). ST12 uses the long-form library labels 'Frontal cortex' and 'Hippocampus (capture 1/2)' for the same tissue, and ST24 (", sd_ref("ST24"), ") writes the region as 'Hippocampus'. Any occurrence of OCC or occipital in these tables refers to Zhou 2023's external occipital cohort, never to our tissue.")),
c("Condition tokens",
  "NHD = the Nasu-Hakola donor; CON = the control donor. Effect signs are always NHD relative to control: positive = higher in NHD."),
c("Detection percentages",
  "pct_NHD / pct_CON / pct_in are the fraction of nuclei with at least one raw RNA UMI. They are computed on the raw RNA counts layer and never on SCT: the per-lane SCT model is aliased with condition in this design and fabricates detection for genes with zero raw UMIs."),
c("Cliff's delta",
  "Rank-based non-parametric effect size in [-1, 1], computed as 2*U1/(n1*n2) - 1 on the raw RNA LogNormalize layer, where U1 is the Mann-Whitney statistic for NHD. Positive = higher in NHD. |delta| >= 0.15 is the project's effect-size floor."),
c("Discovery feature",
  "A gene is a discovery feature in a comparison when p_val_adj < 0.05 (Bonferroni, see Multiple testing) and |cliffs_delta| >= 0.15 and max(pct_NHD, pct_CON) >= 0.10. All three conditions are required; the detection floor is applied on raw RNA."),
c("q* mark in figures",
  "A grey 'q*' printed beside an effect size in a figure marks a gene that is statistically significant (q < 0.05) but falls below the effect-size floor (|Cliff's delta| < 0.15). It is a kept annotation, not an exclusion flag."),
c("Statistical framing",
  "MAST-forward: with one donor per condition the per-nucleus tests are not donor-level inference, so effect sizes and direction consistency carry the claims. Pseudobulk DESeq2 (ST2) is the confirmatory tier. Results are described as direction-consistent across regions and cohorts, never as replicated."),
c("Cross-cohort comparison (Zhou 2023)",
  "Threshold-free direction consistency, not replication. concord_frac is compared with the panel's own empirical baseline (baseline_concord in ST13), never with 0.5. Zhou's tissue is occipital cortex; the comparison is region-matched where possible and otherwise stated as cross-regional."),
c("Depth matching",
  "Analyses sensitive to library depth carry a depth-matched arm. For the microglial transcription-factor analysis (ST14 / Fig 2h) the primary scores are computed at native depth on all nuclei (log-normalized counts) and the same scores recomputed after thinning every nucleus to 750 UMI ship beside them as the sensitivity arm (score_750, q_BH_750). The Visium module scores are computed on spots thinned to 3,000 UMI."),
c("GSEA power tiers",
  "Every GSEA result carries a nucleus-count power tier from the smaller condition arm: solid (>= 150 nuclei), low_conf (50-149), excluded (< 50). Low-confidence panels are flagged figure-wide for that cell type, not per gene."),
c("Multiple testing",
  sprintf("ST1: Bonferroni over all 26,973 assay features (Seurat default), not BH and not restricted to tested features. ST2: BH within each comparison, with DESeq2 independent filtering leaving padj = NA for untested genes. ST3: BH across pathways within subtype x database. ST14: BH once across all %d TF x region tests (%d unique factors), not across the %d factors of the gate audit.", .st14$n_rows, .st14$n_tf, .st14$n_gate)),
c("Missing values",
  "NA means not tested or not estimable and is preserved as NA throughout. A DESeq2 padj of NA (independent filtering) is never recoded to 0 or 1, and an fgsea NES of NA (degenerate permutation null) is never recoded to 0."),
c("p-value underflow",
  "A p-value printed as exactly 0 would be double-precision underflow, not a measured zero, and -log10(0) = Inf breaks rankings and color scales. Every p-value column in these tables is therefore floored at 2.2250739e-308: the smallest normalized double (.Machine$double.xmin = 2.2250738585072014e-308) rounded up to 8 significant digits, so that the floor is written to CSV exactly and a reader cannot find a value below the limit the tables claim to enforce. Any value equal to that floor must be read and reported as p < 2.2e-308, never as p = 0. The number of floored values per column is recorded in the build log."),
c("Transcription-factor regulons",
  "Two different regulons are used and must not be conflated. ST4/ST5 are a genome-wide survey across all six cell types on GTRD, which is unsigned (mode_of_regulation = 1 for every edge); they are not Supplementary Data files and ship as reference tables with the code release (reference_tables/). ST14 is the microglia-only analysis on the signed, curated CollecTRI network and is the only transcription-factor table that backs Fig. 2h."),
c("Neuronal program names: ST7 and Fig. 5d are not the same gene sets",
  "IEG and IEG_core are the same 14 genes under two labels. The other three same-named pairs are not identical, and this is the single most likely way to misread these tables. ST6/ST7's Presynaptic, Postsynaptic and OxPhos are the curated sets Synaptic_presynaptic (21 genes), Synaptic_postsynaptic (14) and Neuronal_OXPHOS (16), scored per nucleus on the atlas. Fig. 5d scores different sets of the same names, defined for that panel (19, 17 and 23 genes) on the full atlas; they ship in ST6 under the prefixed names Fig5d_Presynaptic, Fig5d_Postsynaptic, Fig5d_OxPhos and Fig5d_GluR_ionotropic. The two families differ by 7, 6 and 7 genes, and per-program effect sizes computed from ST7 differ from the Fig. 5d panel by up to 0.08. ST15, ST15b and ST16 are computed from the Fig. 5d definitions, because their purpose is to support that panel."),
c("Myeloid nucleus selection and the CD163+/F13A1+ call",
  "CD163+/F13A1+ nuclei were identified by cluster assignment, not by a marker threshold, and the two must be stated separately. Micro-PVM nuclei were re-embedded, oligodendrocyte doublets removed (module score > 0.5 over PLP1/MBP/MOG/MOBP/CNP/ST18/CLDN11/MAG), and the remainder clustered on the existing shared-nearest-neighbor graph at resolution 0.6. The cluster with the highest median lead of the perivascular-macrophage marker score (MRC1, CD163, LYVE1, F13A1, MS4A7, STAB1, DAB2, MSR1, MAF) over the microglial one (P2RY12, P2RY13, TMEM119, CX3CR1, GPR34, SELPLG, SLC2A5, OLFML3) was designated CD163+/F13A1+ (n = 162), under a pre-set requirement that its lead exceed the runner-up cluster's by more than 0.10. Marker detection is a description of that cluster, not the rule that made it: CD163 is detected in about 50-58% of its nuclei versus about 5% of the main microglial mass. The full chain is ST17b: 3,081 Micro-PVM minus 95 oligodendrocyte doublets = 2,986, minus 162 CD163+/F13A1+ = 2,824 microglia analysed. The 60 Lymphocyte nuclei are a separate Azimuth class that was never inside Micro-PVM, and they remain inside the 55,354 total."),
c("Signature provenance in ST6",
  "ST6's in_scoring_matrix is an annotation match against the atlas gene universe (36,601 symbols), not a detection test: TRUE means the symbol is a row of the expression matrix, never that it was detected in any nucleus. exclusion_reason separates a symbol problem from an expression problem for every gene. Its value below_min_detection is permitted but is empty by construction - no scoring path in this project applies a per-gene detection floor, because module scores are the mean over the member genes present (requiring at least three) and AddModuleScore has no per-gene floor. alias_merged marks published symbols scored under a current HGNC symbol (PSD95 as DLG4, DAP12 as TYROBP, MRE11A as MRE11). symbol_not_mapped marks mouse symbols with no human ortholog symbol in the matrix; all four (H2-T23, IIGP1, LIGP1, GGTA1) are in A1_reactive, so 4 of that signature's 10 genes cannot be scored in human tissue.")
))
names(ST0_CONVENTIONS) <- c("convention", "definition")
