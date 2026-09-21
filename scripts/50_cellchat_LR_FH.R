#!/usr/bin/env Rscript
# =============================================================================
# 50_cellchat_LR_FH.R — CellChat ligand-receptor signaling layer, frontal + hippocampus rebuild (OCC fully dropped).
# Faithful PORT (not a redesign) of two reference generators:
#   (A) glia-only per-region CellChat  == scripts/35_rebuild_glia_cellchat_
#       truncatedMean.R.  Produces the canonical aggregate tables the micro->glia
#       figure scripts consume:
#         tables/cellchat_all_LR_pairs.csv          (truncatedMean, trim=0.1)
#         tables/cellchat_all_LR_pairs_triMean.csv
#       Recipe: full CellChatDB.human; groups = Astro/Micro-PVM/OPC/Oligo (tiny
#       classes Lymphocyte/Endo/Pericytes/Micro-PVM_doublet excluded); RNA
#       LogNormalize data; drop MT/ribo; computeCommunProb(raw.use=TRUE,
#       population.size default FALSE); filterCommunication(min.cells=10).
#       Both methods are produced here directly.
#
#   (B) with-neurons per-region CellChat == scripts/_pending_from_tmp/
#       cellchat_withNeurons_3region_v2.R.  Produces the per-region merged objects
#       the micro->neuron / astro->neuron figure scripts consume:
#         diagnostics/09_cellchat_per_region/cellchat_{CON,NHD}_{region}_withNeurons_v2.rds
#         diagnostics/09_cellchat_per_region/merged_NHDvsCON_{region}_withNeurons_v2.rds
#       Recipe: subsetDB(Secreted Signaling + ECM-Receptor + Cell-Cell Contact);
#       truncatedMean trim=0.1 raw.use population.size=FALSE; min.cells=10;
#       liftCellChat -> mergeCellChat.
#
#       *** neuron grouping NOTE (deliberate, faithful choice) ***
#       The reference v2 script used fine neuron subtypes (Azimuth cortical
#       subclasses for Frontal + Tippani hippo subtypes) sourced from external
#       re-annotation maps (NHD_hippo_neurons_reannotated_v2.rds,
#       neuron_subtype_map.rds).
#       atlas barcodes/clusters and are not valid for the re-clustered FH atlas;
#       re-deriving neuron subtypes is a separate annotation/clustering task
#       (biology lane), not a data-layer port. So here neurons are grouped at
#       major-class resolution (Neuron_Ex / Neuron_Inh from new_annotation), which
#       needs no external map. This is sufficient for the pathway-level
#       micro->neuron panel (script 69), which sums over neuron receivers.
#       The finer per-subtype chord panels (10c/10e) will need a neuron
#       re-annotation before they can run on the FH atlas -> FLAGGED in report.
#
# REGIONS: c("Frontal","Hippo") only. Per (region x condition) fit, then merge
#          CON vs NHD within each region. No OCC anywhere.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(CellChat)
  library(tibble); library(tidyr)
})
set.seed(42)
options(future.globals.maxSize = 16 * 1024^3)
future::plan("sequential")   # memory-safe; deterministic on this machine

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot(!is.na(PROJ))
DIR      <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
OUT_RDS  <- file.path(DIR, "data", "cellchat_per_region"); dir.create(OUT_RDS, showWarnings = FALSE, recursive = TRUE)
TBL      <- file.path(DIR, "tables");                       dir.create(TBL,     showWarnings = FALSE, recursive = TRUE)
LOGDIR   <- file.path(DIR, "logs");                         dir.create(LOGDIR,  showWarnings = FALSE, recursive = TRUE)
ATLAS    <- file.path(DIR, "atlas", "NHD_FH_harmony.rds")
stopifnot(file.exists(ATLAS))

DROP_GENE  <- "^(MT-|RPL|RPS|MRPL|MRPS)"
REGIONS    <- c("Frontal", "Hippo")               # OCC removed
CONDITIONS <- c("CON", "NHD")
GLIA       <- c("Oligo", "OPC", "Micro-PVM", "Astro")  # tiny classes EXCLUDED (match ref 35)
NEURON     <- c("Neuron_Ex", "Neuron_Inh")             # major-class neurons for pipeline (B)

cat("== 50_cellchat_LR_FH ==  ", format(Sys.time()), "\n")
cat("PROJ:", PROJ, "\n")

# ---- STEP 1: load + prep atlas (RNA LogNormalize, exactly as both refs) ------
cat("\n== STEP 1: load atlas ==\n")
obj <- readRDS(ATLAS)
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj, assay = "RNA")
obj <- NormalizeData(obj, verbose = FALSE)
stopifnot(all(c("Region", "Condition", "new_annotation") %in% colnames(obj@meta.data)))

# Hard guard: no OCC, only the two expected regions.
reg_lv <- sort(unique(as.character(obj$Region)))
cat("Region levels present:", paste(reg_lv, collapse = ", "), "\n")
if (any(grepl("OCC|Occ|Occipital", reg_lv)))
  stop("OCC found in atlas Region — FH rebuild expects Frontal+Hippo only")
stopifnot(setequal(reg_lv, REGIONS))

cat("Full atlas new_annotation:\n"); print(table(obj$new_annotation))
genes_to_use <- setdiff(rownames(obj), grep(DROP_GENE, rownames(obj), value = TRUE))
cat("Genes retained (dropped MT/ribo):", length(genes_to_use), "of", nrow(obj), "\n")

# ---- shared CellChat build core ---------------------------------------------
# `db_mode`: "full"  -> full CellChatDB.human       (pipeline A, glia)
#            "subset"-> Secreted+ECM+Cell-Cell       (pipeline B, with-neurons)
# `cc_type`: "truncatedMean" or "triMean"
build_cc <- function(seu, cond, region, groups_col, db_mode, cc_type) {
  tag <- paste0(cond, "_", region)
  cells <- colnames(seu)[seu$Condition == cond & seu$Region == region]
  if (length(cells) < 50) { cat("  SKIP", tag, "(<50 cells)\n"); return(NULL) }
  # keep only cells with an assigned group in this layer
  g_all <- as.character(seu@meta.data[[groups_col]])
  names(g_all) <- colnames(seu)
  cells <- cells[!is.na(g_all[cells])]
  if (length(cells) < 50) { cat("  SKIP", tag, "(<50 grouped cells)\n"); return(NULL) }
  s <- subset(seu, cells = cells)
  grp <- as.character(s@meta.data[[groups_col]])
  ct  <- table(grp)
  keep_ct <- names(ct)[ct >= 10]
  dropped <- setdiff(names(ct), keep_ct)
  if (length(dropped) > 0) {
    cat("  ", tag, "drops group(s) <10 cells: ",
        paste(sprintf("%s(%d)", dropped, ct[dropped]), collapse = ", "), "\n", sep = "")
    keep_cells <- colnames(s)[grp %in% keep_ct]
    s   <- subset(s, cells = keep_cells)
    grp <- as.character(s@meta.data[[groups_col]])
  }
  s$cc_group <- factor(grp, levels = sort(unique(grp)))
  cat("    ", tag, "groups: ", paste(sprintf("%s=%d", levels(s$cc_group),
        as.integer(table(s$cc_group))), collapse = ", "), "\n", sep = "")

  data_use <- GetAssayData(s, assay = "RNA", layer = "data")[genes_to_use, ]
  meta <- data.frame(labels = as.character(s$cc_group),
                     samples = "single_donor",
                     row.names = colnames(data_use))
  cc <- createCellChat(object = data_use, meta = meta, group.by = "labels")
  cc@DB <- if (db_mode == "full") CellChatDB.human else
           subsetDB(CellChatDB.human,
                    search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = cc_type, trim = 0.1,
                          raw.use = TRUE, population.size = FALSE)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  cc
}

# extract per-object LR table with region_cond tag (matches CSV schema of ref 35)
extract_lr <- function(cc, tag) {
  if (is.null(cc)) return(NULL)
  df <- tryCatch(subsetCommunication(cc), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df %>% mutate(region_cond = tag)
}

# ============================================================================
# Pipeline A — GLIA-ONLY (full DB), both methods -> two canonical CSVs
# ============================================================================
cat("\n=========================================================\n")
cat("== PIPELINE A: glia-only CellChat (full DB) — 2 methods ==\n")
cat("=========================================================\n")

# annotate glia cc column once
obj$cc_glia <- ifelse(obj$new_annotation %in% GLIA, as.character(obj$new_annotation), NA_character_)
cat("Glia cells (Astro/Micro-PVM/OPC/Oligo) per region x condition:\n")
print(ftable(factor(obj$cc_glia, levels = GLIA), obj$Region, obj$Condition,
             exclude = NULL))

run_glia_method <- function(cc_type, csv_name) {
  cat("\n-- glia method:", cc_type, "->", csv_name, "--\n")
  cc_list <- list()
  for (cond in CONDITIONS) for (rg in REGIONS) {
    tag <- paste0(cond, "_", rg)
    t0 <- Sys.time()
    cc_list[[tag]] <- build_cc(obj, cond, rg, "cc_glia", "full", cc_type)
    cat("   ", tag, "done in",
        round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
  }
  # aggregate CSV (region_cond schema, intersect common cols like ref 35)
  lr_rows <- Filter(Negate(is.null),
                    lapply(names(cc_list), function(nm) extract_lr(cc_list[[nm]], nm)))
  common_cols <- Reduce(intersect, lapply(lr_rows, colnames))
  all_lr <- bind_rows(lapply(lr_rows, function(x) x[, common_cols]))
  write.csv(all_lr, file.path(TBL, csv_name), row.names = FALSE)
  cat("   wrote", nrow(all_lr), "LR rows ->", csv_name, "\n")
  # per-region merged glia objects (canonical no-suffix names, glia layer)
  merged <- list()
  for (rg in REGIONS) {
    con <- cc_list[[paste0("CON_", rg)]]; nhd <- cc_list[[paste0("NHD_", rg)]]
    if (is.null(con) || is.null(nhd)) { cat("   skip merge", rg, "\n"); next }
    merged[[rg]] <- mergeCellChat(list(CON = con, NHD = nhd),
                                  add.names = c("CON", "NHD"), cell.prefix = TRUE)
  }
  list(cc_list = cc_list, merged = merged, all_lr = all_lr)
}

glia_trunc <- run_glia_method("truncatedMean", "cellchat_all_LR_pairs.csv")
glia_tri   <- run_glia_method("triMean",       "cellchat_all_LR_pairs_triMean.csv")

# save glia rds (canonical = truncatedMean, matching ref 35 which makes trunc canonical)
saveRDS(glia_trunc$cc_list, file.path(OUT_RDS, "cellchat_per_region_list.rds"))
saveRDS(glia_tri$cc_list,   file.path(OUT_RDS, "cellchat_per_region_list_triMean.rds"))
for (rg in names(glia_trunc$merged))
  saveRDS(glia_trunc$merged[[rg]], file.path(OUT_RDS, sprintf("merged_NHDvsCON_%s.rds", rg)))
for (rg in names(glia_tri$merged))
  saveRDS(glia_tri$merged[[rg]], file.path(OUT_RDS, sprintf("merged_NHDvsCON_%s_triMean.rds", rg)))
cat("Glia rds objects written to", OUT_RDS, "\n")

# ============================================================================
# Pipeline B — WITH-NEURONS (subset DB), truncatedMean, major-class neurons
# ============================================================================
cat("\n=============================================================\n")
cat("== PIPELINE B: with-neurons CellChat (subset DB, trunc) ==\n")
cat("=============================================================\n")

obj$cc_neu <- ifelse(obj$new_annotation %in% c(GLIA, NEURON),
                     as.character(obj$new_annotation), NA_character_)
cat("With-neurons groups (glia + Neuron_Ex/Inh) per region x condition:\n")
print(ftable(factor(obj$cc_neu, levels = c(GLIA, NEURON)), obj$Region, obj$Condition,
             exclude = NULL))

neu_cc <- list()
for (rg in REGIONS) for (cond in CONDITIONS) {
  tag <- paste0(cond, "_", rg)
  t0 <- Sys.time()
  cc <- build_cc(obj, cond, rg, "cc_neu", "subset", "truncatedMean")
  neu_cc[[tag]] <- cc
  if (!is.null(cc))
    saveRDS(cc, file.path(OUT_RDS, sprintf("cellchat_%s_withNeurons_v2.rds", tag)))
  cat("   ", tag, "done in",
      round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
}

cat("\n== merge with-neurons per region (liftCellChat -> mergeCellChat) ==\n")
for (rg in REGIONS) {
  con <- neu_cc[[paste0("CON_", rg)]]; nhd <- neu_cc[[paste0("NHD_", rg)]]
  if (is.null(con) || is.null(nhd)) { cat("  skip merge", rg, "\n"); next }
  common_groups <- sort(union(levels(con@idents), levels(nhd@idents)))
  con_l <- liftCellChat(con, common_groups)
  nhd_l <- liftCellChat(nhd, common_groups)
  m <- mergeCellChat(list(CON = con_l, NHD = nhd_l), add.names = c("CON", "NHD"))
  saveRDS(m, file.path(OUT_RDS, sprintf("merged_NHDvsCON_%s_withNeurons_v2.rds", rg)))
  cat("  wrote merged_NHDvsCON_", rg, "_withNeurons_v2.rds (groups: ",
      length(common_groups), ")\n", sep = "")
}

# ============================================================================
# STEP: key-axis verification (report which named LR pairs survived)
# ============================================================================
cat("\n== KEY-AXIS CHECK: microglia-driven signaling in the new data ==\n")

# key ligand->receptor axes of the manuscript
KEY <- tribble(
  ~axis,          ~ligand,  ~receptor_grep,
  "SPP1",         "SPP1",   "CD44|ITGA|ITGB",
  "GAS6",         "GAS6",   "MERTK|TYRO3|AXL",
  "PSAP",         "PSAP",   "GPR37|GPR37L1|LRP1",
  "GRN",          "GRN",    "SORT1|TNFRSF1",
  "SEMA4D",       "SEMA4D", "PLXNB",
  "ENTPD1/purin", "ENTPD1", ".*")

# glia-layer net rows (truncatedMean canonical) for micro->glia
glia_net <- bind_rows(lapply(names(glia_trunc$merged), function(rg) {
  df <- subsetCommunication(glia_trunc$merged[[rg]], slot.name = "net")
  bind_rows(df$CON %>% mutate(Condition = "CON"),
            df$NHD %>% mutate(Condition = "NHD")) %>% mutate(Region = rg)
}))
# neuron-layer net rows for micro->neuron
neu_merged_paths <- file.path(OUT_RDS, sprintf("merged_NHDvsCON_%s_withNeurons_v2.rds", REGIONS))
neu_net <- bind_rows(lapply(REGIONS, function(rg) {
  fp <- file.path(OUT_RDS, sprintf("merged_NHDvsCON_%s_withNeurons_v2.rds", rg))
  if (!file.exists(fp)) return(NULL)
  m <- readRDS(fp); df <- subsetCommunication(m, slot.name = "net")
  bind_rows(df$CON %>% mutate(Condition = "CON"),
            df$NHD %>% mutate(Condition = "NHD")) %>% mutate(Region = rg)
}))

report_axis <- function(net, sender, receivers, label) {
  cat("\n--", label, "(", sender, "-> {", paste(receivers, collapse=","), "}) --\n")
  if (is.null(net) || nrow(net) == 0) { cat("  (no net rows)\n"); return(invisible()) }
  sub <- net %>% filter(source == sender, target %in% receivers)
  for (i in seq_len(nrow(KEY))) {
    lg <- KEY$ligand[i]; rgx <- KEY$receptor_grep[i]
    hit <- sub %>% filter(ligand == lg, grepl(rgx, receptor))
    if (nrow(hit) == 0) { cat(sprintf("  %-14s : ABSENT\n", KEY$axis[i])); next }
    smry <- hit %>% group_by(Region, Condition) %>%
      summarise(n = dplyr::n(), .groups = "drop") %>%
      mutate(rc = paste0(Region, "-", Condition, "(", n, ")"))
    cat(sprintf("  %-14s : PRESENT  %s\n", KEY$axis[i], paste(smry$rc, collapse=" ")))
  }
}
report_axis(glia_net, "Micro-PVM", c("Astro", "Oligo", "OPC"), "micro -> glia (glia layer)")
report_axis(neu_net,  "Micro-PVM", NEURON, "micro -> neuron (neuron layer)")

# significant-pair counts per region x condition (pval < 0.05)
cat("\n== Significant LR pairs (pval<0.05) per region x condition ==\n")
sig_counts <- function(net, layer) {
  if (is.null(net) || nrow(net) == 0) { cat("  (", layer, ": no rows)\n"); return(invisible()) }
  s <- net %>% filter(pval < 0.05) %>%
    group_by(Region, Condition) %>% summarise(n_sig = dplyr::n(), .groups = "drop")
  cat("  [", layer, "]\n"); print(as.data.frame(s))
}
sig_counts(glia_net, "glia layer, truncatedMean")
sig_counts(neu_net,  "neuron layer, truncatedMean")

# ---- provenance -------------------------------------------------------------
prov <- file.path(LOGDIR, "50_cellchat_LR_FH_sessionInfo.txt")
writeLines(c(paste("run:", format(Sys.time())),
             paste("PROJ:", PROJ),
             capture.output(sessionInfo())), prov)
cat("\nProvenance ->", prov, "\n")
cat("\n=== DONE ===\n")
