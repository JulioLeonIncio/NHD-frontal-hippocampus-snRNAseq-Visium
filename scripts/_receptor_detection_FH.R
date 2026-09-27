# =============================================================================
# _receptor_detection_FH.R — Single source for "is this receptor detected in the receiving population?", used by the ligand-receptor panels (Fig. 3g micro->astro, Fig. 3h astro->neuron).
# Why. Those panels drew a grey cross meaning "receptor below
# the 10 % detection floor", but computed it as `is.na(rec_delta)` — the receptor being absent
# from the receiver's MAST table. `40_mast_percell_dual_FH.R` runs FindMarkers(min.pct = 0.1,
# logfc.threshold = 0.1), so a gene is dropped from that table for failing either the detection
# floor or the fold-change floor. Stable, highly expressed receptors were therefore drawn as
# undetectable: LRRTM4 is detected in 84.1 % / 83.2 % of frontal inhibitory neurons (delta +0.03),
# CADM1 in 60-78 % of hippocampal neurons, CLSTN1 in 59.6 % / 86.5 % of frontal excitatory neurons.
# The mislabel inverted the panel's own message — a detected, unchanged receptor is evidence for
# preserved input, not missing data.
#
# What: detection fraction per (gene, cell class, region, condition) computed on raw RNA counts
# from the atlas, cached with an mtime guard against the atlas.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
receptor_detection <- function(proj, atlas_path, genes, classes, regions,
                               cache = file.path(proj, "tables", "_receptor_detection_FH.rds")) {
  genes <- sort(unique(genes))
  ok <- file.exists(cache) && file.mtime(cache) >= file.mtime(atlas_path)
  if (ok) { d <- readRDS(cache)
    ok <- all(genes %in% d$gene) && all(classes %in% d$class) && all(regions %in% d$region)
    if (ok) return(d[d$gene %in% genes & d$class %in% classes & d$region %in% regions, ]) }
  suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
  message("receptor_detection: building cache from ", basename(atlas_path))
  obj <- readRDS(atlas_path); DefaultAssay(obj) <- "RNA"
  if (length(SeuratObject::Layers(obj, assay = "RNA")) > 1) obj <- SeuratObject::JoinLayers(obj, assay = "RNA")
  md <- obj@meta.data; cnt <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
  gs <- intersect(genes, rownames(cnt))
  out <- do.call(rbind, lapply(classes, function(cl) do.call(rbind, lapply(regions, function(rg) {
    cells <- rownames(md)[md$new_annotation == cl & md$Region == rg]
    if (!length(cells) || !length(gs)) return(NULL)
    pct <- sapply(c("CON", "NHD"), function(cd) {
      cc <- cells[md[cells, "Condition"] == cd]
      if (!length(cc)) rep(NA_real_, length(gs)) else Matrix::rowMeans(cnt[gs, cc, drop = FALSE] > 0) })
    data.frame(gene = gs, class = cl, region = rg, pct_CON = pct[, "CON"], pct_NHD = pct[, "NHD"],
               n_CON = sum(md[cells, "Condition"] == "CON"), n_NHD = sum(md[cells, "Condition"] == "NHD"),
               row.names = NULL, stringsAsFactors = FALSE) }))))
  out$pct_max <- pmax(out$pct_CON, out$pct_NHD, na.rm = TRUE)
  saveRDS(out, cache); rm(obj, cnt); gc()
  out[out$gene %in% genes & out$class %in% classes & out$region %in% regions, ]
}
# TRUE where the receptor really is below the floor in that receiver (default 10 %, the paper's floor)
receptor_below_floor <- function(det, gene, class, region, floor = 0.10) {
  i <- match(paste(gene, class, region), paste(det$gene, det$class, det$region))
  ifelse(is.na(i), TRUE, det$pct_max[i] < floor)   # never measured -> treat as below the floor
}
