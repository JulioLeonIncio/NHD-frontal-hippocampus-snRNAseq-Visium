#!/usr/bin/env Rscript
# =============================================================================
# 40b_mast_microglia_only_FH.R — Per-nucleus MAST for the split myeloid classes.
# -----------------------------------------------------------------------------
# Why: every Figure-2 downstream analysis was run on the pooled
# "Micro-PVM" atlas class. PVM are 2.5% of CON myeloid nuclei but 15% of NHD myeloid
# nuclei, so the mixing fraction changes across the very contrast being tested.
#   * The stats audit: the referee's question is not "is it robust" but "what cell is
#     this" — the claim's noun is *microglia*, so the analysis unit must match it.
#   * The biology audit adds the decisive point: DAM staging is defined as loss of the
#     parenchymal homeostatic checkpoint (P2RY12/TMEM119/CX3CR1/SALL1), which PVM never
#     express (Van Hove 2019, PMID 31061494). A cell with no homeostatic signature
#     cannot "lose" one, so DAM-1 is undefined for PVM and "DAM-1 arrest" on a pooled
#     object is a category error independent of effect size (Keren-Shaul 2017
#     PMID 28602351; Deczkowska 2018 PMID 29775591).
# The pooled result is not wrong (Cliff's-delta agreement r = 0.9986/0.9963, see
# 2Q2) — it is mis-named. This script produces the correctly-named primary.
#
# Model, covariates, gene universe, Cliff's delta and the |delta| >= 0.15 discovery
# call are identical to 40_mast_percell_dual_FH.R. The only change is the cell set:
# "Micro-PVM" is replaced by the split classes from the immune-compartment assignment
# (2N), which also removes the 95 oligo-doublets.
#
# Also runs the PVM-INTRINSIC NHD-vs-CON contrast. The composition-only argument
# assumes PVM per-cell expression is condition-invariant — an assumption that is
# questionable when TYROBP is null and PVM express TYROBP, so it is tested rather
# than assumed. It is UNDERPOWERED (CON 19 / NHD 47 Frontal; 53 / 43 Hippo) and must
# be reported as such.
#
# Output: tables/mast_dual/MAST_{Microglia,PVM}_{Frontal,Hippo}.csv
#         tables/mast_dual/MAST_microglia_only_discovery_counts.csv
# The pooled MAST_Micro-PVM_*.csv files are left in place as the supplementary
# sensitivity/disclosure tier — do not delete them.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(MAST); library(future)
})
set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)
plan("sequential")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
DIR <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
OUT <- file.path(DIR, "tables/mast_dual"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CLS <- file.path(DIR, "tables/micro_states/immune_umap_classes_FH.csv")
stopifnot("run 2N_immune_UMAP_FH.R first" = file.exists(CLS))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

msg("=== loading F+Hippo atlas %s ===", format(Sys.time(), "%H:%M:%S"))
obj <- readRDS(file.path(DIR, "atlas/NHD_FH_harmony.rds"))

# --- apply the split myeloid classes ---------------------------------------
cls <- read.csv(CLS)
obj$myeloid_class <- cls$class[match(colnames(obj), cls$barcode)]   # NA outside the compartment
msg("myeloid classes attached: %s",
    paste(sprintf("%s=%d", names(table(obj$myeloid_class)), table(obj$myeloid_class)), collapse = " "))

# --- covariates: identical derivation to script 40 -------------------------
samp <- as.character(obj$Sample)
obj$matter <- factor(ifelse(grepl("_g$", samp), "grey",
                     ifelse(grepl("_w$", samp), "white", NA_character_)),
                     levels = c("grey","white"))
obj$hippo_area <- factor(ifelse(grepl("Hippo1$", samp), "area1",
                         ifelse(grepl("Hippo2$", samp), "area2", NA_character_)),
                         levels = c("area1","area2"))
obj$Condition   <- factor(as.character(obj$Condition), levels = c("CON","NHD"))
obj$log10nCount <- log10(obj$nCount_RNA)

# RNA is needed for the raw detection / effect-size measurement below; join the
# per-lane layers on the full object first and normalise once (a subset with lanes
# contributing zero nuclei leaves empty layers that StitchMatrix cannot handle).
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"
msg("PrepSCTFindMarkers ... %s", format(Sys.time(), "%H:%M:%S"))
obj <- PrepSCTFindMarkers(obj, verbose = FALSE); invisible(gc())

genes_keep <- grep("^(MT-|RPL|RPS|MRPL|MRPS)", rownames(obj), value = TRUE, invert = TRUE)
msg("features kept (MT/ribo dropped): %d / %d", length(genes_keep), nrow(obj))

cliffs_delta <- function(seu, genes, grp) {
  if (!length(genes)) return(setNames(numeric(0), character(0)))
  expr <- as.matrix(GetAssayData(seu, assay = "SCT", layer = "data")[genes, , drop = FALSE])
  isN <- grp == "NHD"; n1 <- sum(isN); n2 <- sum(!isN)
  apply(expr, 1, function(x) { r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
}

MIN_PER_COND <- 10
run_one <- function(cls_name, rg) {
  label <- paste0(cls_name, "_", rg)
  cells <- colnames(obj)[which(!is.na(obj$myeloid_class) &
                               obj$myeloid_class == cls_name & obj$Region == rg)]
  if (length(cells) < 30) { msg("  SKIP %-22s (<30 nuclei: %d)", label, length(cells)); return(NULL) }
  sub <- subset(obj, cells = cells)
  tabc <- table(sub$Condition)
  if (any(tabc[c("CON","NHD")] < MIN_PER_COND)) {
    msg("  SKIP %-22s (per-cond %s)", label, paste(tabc, collapse = "/")); return(NULL) }
  msg("  RUN  %-22s  CON %d / NHD %d", label, tabc[["CON"]], tabc[["NHD"]])
  Idents(sub) <- "Condition"
  block <- if (rg == "Hippo") "hippo_area" else "matter"
  has_block <- length(unique(na.omit(sub@meta.data[[block]]))) >= 2
  lv <- if (has_block) c(block, "log10nCount") else "log10nCount"
  m <- tryCatch(
    FindMarkers(sub, ident.1 = "NHD", ident.2 = "CON", test.use = "MAST",
                latent.vars = lv, assay = "SCT", features = genes_keep,
                min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi = FALSE, verbose = FALSE),
    error = function(e) { msg("  MAST err %s: %s", label, conditionMessage(e)); NULL })
  if (is.null(m) || !nrow(m)) return(NULL)
  m$gene <- rownames(m)
  # ===========================================================================
  # Detection and effect size are measured on raw/RNA DATA, not on SCT.
  # KCNIP4, RBFOX1, CSMD1,
  # SYT1 and OPCML have raw max count = 0 across all 267 NHD and all 730 CON frontal
  # microglia, yet SCT reports 98.9% vs 11.0% detection — and they were the TOP
  # discovery features (delta 0.60-0.69). The per-lane SCT model is aliased with
  # Condition (NHD frontal lanes are ~48% neurons vs 13% in CON), so neuronal
  # transcripts get a large expected value and non-zero corrected counts.
  # 21.3% of the Frontal discovery set had raw max(pct) < 0.10 in both conditions.
  # MAST still supplies the p-value on SCT (MAST-forward policy unchanged); pct and
  # Cliff's delta now come from RNA, and the discovery gate additionally requires
  # real detection. Genuinely expressed genes are unaffected (P2RY12/F13A1/CD163
  # give identical raw and SCT pct).
  # ===========================================================================
  rna <- GetAssayData(sub, assay = "RNA", layer = "counts")[m$gene, , drop = FALSE]
  isN <- sub$Condition == "NHD"
  m$pct.1 <- round(Matrix::rowMeans(rna[, isN,  drop = FALSE] > 0), 3)
  m$pct.2 <- round(Matrix::rowMeans(rna[, !isN, drop = FALSE] > 0), 3)
  # Cliff's delta on RNA log-normalised data, for every tested gene (not only the
  # significant ones) so that null claims carry an effect size too.
  rnad <- GetAssayData(sub, assay = "RNA", layer = "data")[m$gene, , drop = FALSE]
  n1 <- sum(isN); n2 <- sum(!isN)
  m$cliffs_delta <- apply(as.matrix(rnad), 1, function(x) {
    r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
  m$max_pct_raw <- pmax(m$pct.1, m$pct.2)
  out <- m %>% transmute(gene, cell_type = cls_name, region = rg,
                         avg_log2FC, p_val, p_val_adj, pct.1, pct.2, cliffs_delta,
                         max_pct_raw,
                         n_CON = tabc[["CON"]], n_NHD = tabc[["NHD"]],
                         latent_vars = paste(lv, collapse = "+"),
                         discovery = p_val_adj < 0.05 & !is.na(cliffs_delta) &
                                     abs(cliffs_delta) >= 0.15 &
                                     max_pct_raw >= 0.10)
  # fail loud if a discovery feature is not actually detected in the raw data
  stopifnot("discovery feature below the raw detection floor" =
              all(out$max_pct_raw[out$discovery] >= 0.10))
  write.csv(out, file.path(OUT, sprintf("MAST_%s.csv", label)), row.names = FALSE)
  msg("       -> %d tested, %d discovery", nrow(out), sum(out$discovery))
  out
}

all_df <- bind_rows(lapply(c("Microglia","PVM"), function(ct)
  bind_rows(lapply(c("Frontal","Hippo"), function(rg) run_one(ct, rg)))))

summ <- all_df %>% group_by(cell_type, region) %>%
  summarise(n_tested = n(), n_discovery = sum(discovery),
            up_NHD = sum(discovery & cliffs_delta > 0),
            down_NHD = sum(discovery & cliffs_delta < 0),
            # NB use [1], not first(): MAST attaches S4Vectors, whose first() masks
            # dplyr's and errors on a plain integer vector.
            n_CON = n_CON[1], n_NHD = n_NHD[1], .groups = "drop")
write.csv(all_df, file.path(OUT, "MAST_microglia_only_all.csv"), row.names = FALSE)
write.csv(summ, file.path(OUT, "MAST_microglia_only_discovery_counts.csv"), row.names = FALSE)

msg("\n=== DISCOVERY COUNTS (split myeloid classes, |delta|>=0.15) ===")
print(as.data.frame(summ))
msg("\nNOTE: PVM rows are UNDERPOWERED and must be reported as such.")
msg("=== DONE: 40b_mast_microglia_only_FH ===")
