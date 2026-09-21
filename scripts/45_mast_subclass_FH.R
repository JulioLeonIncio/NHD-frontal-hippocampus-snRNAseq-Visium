#!/usr/bin/env Rscript
# =============================================================================
# 45_mast_subclass_FH.R — Per-nucleus MAST at neuron subclass resolution, frontal+hippo rebuild. Faithful port of manuscript_7fig/scripts/45_mast_subclass.R.
# -----------------------------------------------------------------------------
# MIRRORS the broad MAST (40_mast_percell_dual_FH.R) exactly — same model,
# covariates, Cliff's-delta helper, MT/ribo drop, min.pct/logfc.threshold,
# recorrect_umi=FALSE — but iterates over the cortical neuron subclasses in the
# frontal cortex, restricted to bona-fide neurons (new_annotation %in% {Neuron_Ex,
# Neuron_Inh}) so force-mapped glia are excluded.
#
# The per-nucleus label comes from scripts/_neuron_subclass_FH.R
#   (data/neuron_subtype_map_FH.rds), not from atlas `predicted.subclass`. Frontal
#   excitatory nuclei are labelled by the Jorstad 2023 DLPFC transfer (116b; adds
#   L4 IT), inhibitory nuclei keep the Azimuth subclass.
#
#   subclasses = L2/3 IT, L4 IT, L5 IT, L6 IT, L6 IT Car3, L6 CT, L6b, L5/6 NP,
#                Pvalb, Sst, Vip, Lamp5, Sncg          (full neocortical census, 13)
#   region     = Frontal only (the F+Hippo rebuild has no OCC; hippocampal
#                neurons use a de-novo taxonomy, not the cortical Azimuth bins,
#                and are handled as prerank-only comparisons in 47b/44b).
#
#   excluded: L5 ET (CON=2), Sst Chodl (NHD=0). These are
#   below the MIN_PER_COND floor and are also auto-skipped by the gate; they are
#   simply not in `subclasses` so the log shows them as census exclusions.
#
# Model:
#   Idents      = Condition (NHD vs CON)
#   latent.vars = matter(grey/white) + log10(nCount_RNA)   for Frontal
#                 (block term added only if >=2 of its levels are present)
# Gate: skip a comparison with <30 nuclei total or <MIN_PER_COND nuclei/condition.
#
# OUTPUT (schema == broad MAST_dual_discovery_all.csv so downstream is drop-in):
#   tables/mast_dual/MAST_<label>.csv                 (per comparison)
#   tables/mast_dual/MAST_subclass_discovery_all.csv  (combined)
#   tables/mast_dual/MAST_subclass_discovery_counts.csv
#   label = "<sanitized subclass>_Frontal", e.g. "L2_3_IT_Frontal".
#
# Heavy (reloads the 2.9 GB atlas + PrepSCTFindMarkers). Run detached; log tee'd.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(MAST); library(future)
})
set.seed(42)
options(future.globals.maxSize = 30 * 1024^3)   # PrepSCTFindMarkers exports a big global
plan("sequential")                              # memory-safe, no parallel workers

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot("Could not resolve PROJ root" = !is.na(PROJ), dir.exists(PROJ))
DIR  <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
OUT  <- file.path(DIR, "tables/mast_dual"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

# ---------------------------------------------------------------------------
# Back up the combined subclass outputs before this run overwrites them
# (idempotency / rollback safety; timestamped, never clobbered).
# ---------------------------------------------------------------------------
.ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
for (bn in c("MAST_subclass_discovery_all.csv", "MAST_subclass_discovery_counts.csv")) {
  fp <- file.path(OUT, bn)
  if (file.exists(fp)) {
    bak <- file.path(OUT, sub("\\.csv$", sprintf("_bak_%s.csv", .ts), bn))
    file.copy(fp, bak, overwrite = FALSE); msg("backed up %s -> %s", bn, basename(bak))
  }
}

# ---------------------------------------------------------------------------
# Load + prep  (identical prep to the broad MAST 40_FH)
# ---------------------------------------------------------------------------
ATLAS <- file.path(DIR, "atlas/NHD_FH_harmony.rds")
stopifnot("MISSING atlas RDS — prep step failed" = file.exists(ATLAS))
msg("=== loading F+Hippo atlas %s ===", format(Sys.time(), "%H:%M:%S"))
obj <- readRDS(ATLAS)

stopifnot("MISSING metadata column new_annotation"     = "new_annotation"     %in% colnames(obj@meta.data),
          "MISSING metadata column Region"             = "Region"             %in% colnames(obj@meta.data),
          "MISSING metadata column Condition"          = "Condition"          %in% colnames(obj@meta.data),
          "MISSING metadata column Sample"             = "Sample"             %in% colnames(obj@meta.data))
# Subclass label from the single-source map (DFC Ex + Azimuth Inh), by barcode.
source(file.path(DIR, "scripts", "_neuron_subclass_FH.R"))   # SUBCLASS_MAP, CTX_SUBCLASSES, subclass_of()
obj$neuron_subclass <- subclass_of(colnames(obj))
msg("subclass labels joined: %d nuclei (%d Frontal neurons); source: %s",
    sum(!is.na(obj$neuron_subclass)),
    sum(!is.na(obj$neuron_subclass) & obj$Region == "Frontal"), CTX_SUBCLASS_SOURCE)

samp <- as.character(obj$Sample)
obj$matter <- factor(ifelse(grepl("_g$", samp), "grey",
                     ifelse(grepl("_w$", samp), "white", NA_character_)),
                     levels = c("grey","white"))
obj$hippo_area <- factor(ifelse(grepl("Hippo1$", samp), "area1",
                         ifelse(grepl("Hippo2$", samp), "area2", NA_character_)),
                         levels = c("area1","area2"))
obj$Condition   <- factor(as.character(obj$Condition), levels = c("CON","NHD"))
obj$log10nCount <- log10(obj$nCount_RNA)

# restrict cortical-subclass tests to bona-fide neurons (the map only carries
# Neuron_Ex / Neuron_Inh, but keep the explicit flag so the mask is self-documenting).
obj$is_neuron <- obj$new_annotation %in% c("Neuron_Ex", "Neuron_Inh")
stopifnot("map label on a non-neuron nucleus" = all(is.na(obj$neuron_subclass) | obj$is_neuron))

DefaultAssay(obj) <- "SCT"
msg("PrepSCTFindMarkers ... %s", format(Sys.time(), "%H:%M:%S"))
obj <- PrepSCTFindMarkers(obj, verbose = FALSE); invisible(gc())
# RNA must be joined + normalized on the full object before any subset,
# otherwise GetAssayData(layer="data") hands back raw counts. Detection, Cliff's delta
# and the pct floor are all computed on RNA (SCT corrected counts fabricate expression
# for zero-UMI genes -- see 40_mast_percell_dual_FH.R).
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])
DefaultAssay(obj) <- "RNA"; obj <- NormalizeData(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"

genes_keep <- grep("^(MT-|RPL|RPS|MRPL|MRPS)", rownames(obj), value = TRUE, invert = TRUE)
msg("features kept (MT/ribo dropped): %d / %d", length(genes_keep), nrow(obj))

# ---------------------------------------------------------------------------
# Cliff's delta (identical midrank Mann-Whitney helper to the broad MAST)
# ---------------------------------------------------------------------------
cliffs_delta <- function(seu, genes, grp) {
  if (!length(genes)) return(setNames(numeric(0), character(0)))
  expr <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "data")[genes, , drop = FALSE])
  isN <- grp == "NHD"; n1 <- sum(isN); n2 <- sum(!isN)
  apply(expr, 1, function(x) { r <- rank(x); U1 <- sum(r[isN]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1 })
}

# ---------------------------------------------------------------------------
# Per (subclass x Frontal) MAST.
#   cell_type column = subclass label exactly as in the map (so downstream can
#                      subset by neuron_subtype == cell_type).
#   comparison       = "<sanitized subclass>_Frontal" (filename-safe key).
# EXCLUDED (unpowered, stated): L5 ET (CON=2), Sst Chodl (NHD=0).
# ---------------------------------------------------------------------------
# List = CTX_SUBCLASSES minus the two unpowered classes (adds L4 IT).
subclasses <- setdiff(CTX_SUBCLASSES, c("L5 ET","Sst Chodl"))   # 13: 8 Ex + 5 Inh
stopifnot("expected 13 subclasses incl. L4 IT" = length(subclasses) == 13 && "L4 IT" %in% subclasses)
regions      <- c("Frontal")         # cortical only; no OCC, hippo neurons = de-novo
MIN_PER_COND <- 10                   # per-condition nucleus floor (matches broad MAST)
MIN_CELLS    <- 30                   # total nucleus floor
FORCE        <- identical(Sys.getenv("FORCE_MAST"), "1")  # recompute even if cache exists
sanitize <- function(x) gsub("[ /]", "_", x)
msg("CENSUS excludes (unpowered, not run): L5 ET, Sst Chodl")

run_one <- function(sc, rg) {
  label <- paste0(sanitize(sc), "_", rg)
  of    <- file.path(OUT, sprintf("MAST_%s.csv", label))
  # a cache is reusable only if it postdates both the atlas and the subclass map
  if (!FORCE && file.exists(of) && file.mtime(of) > max(file.mtime(ATLAS), NEURON_MAP_MTIME)) {
    msg("  REUSE %-22s (existing MAST cache newer than atlas + map; FORCE_MAST=1 to recompute)", label)
    return(read.csv(of, stringsAsFactors = FALSE))
  }
  # Key on the map label (neuron_subclass), never on predicted.subclass
  cells <- colnames(obj)[obj$is_neuron & !is.na(obj$neuron_subclass) &
                         obj$neuron_subclass == sc & obj$Region == rg]
  if (length(cells) < MIN_CELLS) { msg("  SKIP %-22s (<%d nuclei)", label, MIN_CELLS); return(NULL) }
  sub <- subset(obj, cells = cells)
  tabc <- table(sub$Condition)
  if (any(tabc[c("CON","NHD")] < MIN_PER_COND) || length(unique(sub$Condition)) < 2) {
    msg("  SKIP %-22s (per-cond nuclei %s)", label, paste(tabc, collapse="/")); return(NULL) }
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
  sig <- m$gene[which(m$p_val_adj < 0.05)]
  # Detection + effect size on raw RNA (see 40_mast_percell_dual_FH.R).
  rna <- GetAssayData(sub, assay = "RNA", layer = "counts")[m$gene, , drop = FALSE]
  isN <- sub$Condition == "NHD"
  m$pct.1 <- round(Matrix::rowMeans(rna[,  isN, drop = FALSE] > 0), 3)
  m$pct.2 <- round(Matrix::rowMeans(rna[, !isN, drop = FALSE] > 0), 3)
  m$max_pct_raw <- pmax(m$pct.1, m$pct.2)
  d <- cliffs_delta(sub, m$gene, sub$Condition)
  m$cliffs_delta <- unname(d[m$gene])
  out <- m %>% transmute(
    comparison = label, cell_type = sc, region = rg, gene,
    avg_log2FC, p_val, p_val_adj, pct.1, pct.2, cliffs_delta,
    n_NHD = as.integer(tabc["NHD"]), n_CON = as.integer(tabc["CON"]),
    latent_vars = paste(lv, collapse = "+"),
    discovery = !is.na(p_val_adj) & p_val_adj < 0.05 &
                !is.na(cliffs_delta) & abs(cliffs_delta) >= 0.15 &
                max_pct_raw >= 0.10)
  tmp <- paste0(of, ".tmp"); write.csv(out, tmp, row.names = FALSE)
  if (!file.rename(tmp, of)) { file.copy(tmp, of, overwrite = TRUE); unlink(tmp) }
  msg("  %-22s nuclei=%d (NHD %d/CON %d) | sig(p<.05)=%d | discovery(|d|>=.15)=%d",
      label, ncol(sub), tabc["NHD"], tabc["CON"], length(sig), sum(out$discovery))
  invisible(gc()); out
}

all_res <- list()
for (rg in regions) for (sc in subclasses) {
  msg("=== %s_%s  START %s ===", sanitize(sc), rg, format(Sys.time(), "%H:%M:%S"))
  all_res[[paste0(sanitize(sc),"_",rg)]] <- run_one(sc, rg)
}

# Rebuild the combined table from all per-comparison files on disk (reproduces the
# discovery_all schema exactly; reuse-safe if a prior run left some on disk).
census_labels <- as.vector(outer(sanitize(subclasses), regions, paste, sep = "_"))
files <- file.path(OUT, sprintf("MAST_%s.csv", census_labels))
files <- files[file.exists(files)]
msg("\n=== rebuilding discovery_all from %d per-comparison file(s) on disk ===", length(files))
all_df <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
write.csv(all_df, file.path(OUT, "MAST_subclass_discovery_all.csv"), row.names = FALSE)

msg("\n=== SUBCLASS DISCOVERY COUNTS (per-nucleus MAST, |delta|>=0.15) ===")
if (nrow(all_df)) {
  summ <- all_df %>% filter(discovery) %>%
    group_by(comparison) %>%
    summarise(n_discovery = dplyr::n(),
              up_NHD = sum(avg_log2FC > 0), up_CON = sum(avg_log2FC < 0), .groups = "drop")
  print(as.data.frame(summ), row.names = FALSE)
  write.csv(summ, file.path(OUT, "MAST_subclass_discovery_counts.csv"), row.names = FALSE)
} else {
  msg("  (no comparison passed the gate)")
}
cat("\n=== DONE ===\n", file = stderr()); print(sessionInfo())
