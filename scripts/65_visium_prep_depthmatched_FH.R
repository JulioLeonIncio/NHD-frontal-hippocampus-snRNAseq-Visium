#!/usr/bin/env Rscript
# =============================================================================
# 65_visium_prep_depthmatched_FH.R — Shared prep for the rebuilt Visium figure.
# Loads the Visium object once, thins every spot to a common UMI depth, and emits
# tidy tables the panel scripts plot from.
# -----------------------------------------------------------------------------
# Why depth matching is not optional HERE (REBUILD_PLAN entries 32-33):
# NHD spots carry half the depth of CON spots (order of magnitude: median ~4.5 k vs
# ~9 k UMI and ~2.3 k vs ~3.7 k genes per spot, consistent within condition; the live
# per-section values are written to tables/visium_dm_qc_FH.csv by this script).
# Detection prevalence is therefore never a valid cross-condition readout in this
# dataset: raw CUX2 prevalence "falls" to 0.699 of control and is entirely artifact
# (1.046 once matched), while PCP4 at 0.279 is entirely real (0.278 matched). The
# same panel contained both, and only matching separated them.
#
# Good news for Methods, measured here: at a matched 3,000 UMI, CON detects 1,956
# and NHD 1,835 median genes — only ~6% lower. The raw 3,689 vs 2,324 gene gap is
# almost all depth, not a biological complexity difference.
#
# Emits:
#   tables/visium_dm_marker_counts_FH.csv   per gene x section mean counts / 3,000 UMI
#   tables/visium_dm_module_spots_FH.csv    per spot module scores + c2l content
#   tables/visium_dm_qc_FH.csv              the depth-matching audit itself
#   tables/visium_dm_module_genes_FH.csv    module membership as scored (module, gene, present)
#
# House: spot-quality floor before any c2l use;
# gene-set key so a set change forces a rebuild (an mtime guard cannot see one).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(tidyr)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
TDIR <- file.path(PROJ, "tables"); LOGD <- file.path(PROJ, "logs")
for (d in c(TDIR, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

VIS  <- file.path(dirname(PROJ), "Visium", "integrated_harmony",
                  "NHD_frontal_integrated_harmony.rds")
PROP <- file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN",
                  "proportions_q05.csv")
stopifnot("MISSING Visium object" = file.exists(VIS),
          "MISSING c2l proportions" = file.exists(PROP))

DEPTH_T    <- 3000L   # every retained spot thinned to this
SPOT_FLOOR <- 1000L   # c2l quality floor

# ---------------------------------------------------------------------------
# Marker panel — grouped so the panel can show selectivity, not just loss
# ---------------------------------------------------------------------------
# The upper-layer group is the internal positive control: if CUX2 is unchanged in
# the same sections where PCP4 collapses, "selective deep-layer loss" is licensed;
# without it the honest reading degrades to "global degradation".
MARKERS <- list(
  "Deep-layer neuron" = c("PCP4","RORB","TLE4","FEZF2"),
  "Upper-layer neuron"= c("CUX2","CUX1","LAMP5"),
  "Pan-neuronal"      = c("SNAP25","RBFOX3","SYT1"),
  "Myelin"            = c("MBP","PLP1","MOBP","MAG","CNP"),
  # Pooling them with AIF1/C1QB/CD74 under one
  # "glial activation" heading let SERPINA3 (+3.54, the largest effect in the dataset) be
  # read as evidence about microglia. The two limbs are now named separately.
  "Microglial activation" = c("AIF1","C1QB","CD74","CSF1R"),   # TYROBP removed — the NHD donor is TYROBP-null, so it biases any NHD score downward by construction
  "Astrocyte reactive"    = c("GFAP","SERPINA3","VIM","CD44"))

# ---------------------------------------------------------------------------
# Microglial programmes
# ---------------------------------------------------------------------------
# The shipped Fig-6 modules are carried by pan-cellular genes: leave-one-out gives
# DAM-2 = FTH1 (-0.62) + SPP1 (-0.39) with CD9 -0.006 and LPL -0.074 contributing
# nothing, and inflammation = APOE (+0.27) + B2M (+0.21). At 55 um those genes
# report the compartment mixture, not a microglial state — the same standard that
# made us refuse to plot the oligodendrocyte stress module on tissue.
# FTH1 (98.6% of spots), B2M (87.9%), APOE (81.7%), CTSD (62%) are therefore OUT.
MODULES <- list(
  # TREM2-dependent lipid arm, restricted to what is both DAM-2 and detectable
  "TREM2-dep lipid (DAM-2)" = c("SPP1","CD9","LPL","GPNMB"),
  # TREM2-independent activation, microglia-restricted
  "Complement / MHC-II"     = c("C1QA","C1QB","C1QC","CD74","HLA-DRA","CD14"),
  # Neutral NULLS — the baseline that shows the axis does not move on its own
  "Housekeeping (null)"     = c("ACTB","GAPDH","TUBB","PGK1","RPL13A","RPS18",
                                "PPIA","B2M","TBP","UBC","YWHAZ","SDHA"))
POSCTRL <- c("CD74","C1QB")   # single genes that must rise with microglial content
# v3 — genes the Figure-6 cascade column needs per spot.
# Layer: the laminar-depth panel tests Figure 5's deep-layer vulnerability in tissue, so
# it needs deep and upper markers measured in the same spots at the same depth. The upper
# set is the internal positive control that licenses the word "selective" -- without it
# a deep-layer fall is indistinguishable from global degradation.
LAYER_GENES <- c(deep = c("PCP4","TLE4","FOXP2","SEMA3E","HS3ST4"),
                 upper = c("CUX2","RORB","LAMP5","CALB1"),
                 class = c("SLC17A7","GAD1","GAD2"),
                 myelin = c("MBP","PLP1","MOBP"),
                 glial = c("AIF1","SERPINA3","TREM2","GFAP"))          # TYROBP removed (TYROBP-null donor); it remains a row in the Fig 6d gene heatmap, where its absence is the point
LAYER_GENES <- unname(LAYER_GENES)
EXPORT_GENES <- unique(c(POSCTRL, LAYER_GENES))

# The spatial figure showed myelin alone, when the finding from
# Figure 4 is the supply chain — myelin protein held while the lipid arms that build the
# membrane fail. Every arm clears the 55 um detection floor on this slide (sterol 13/13
# genes, most >45% of spots; SCD 78%; UGT8 22%), so the dissociation is measurable in
# tissue and not only in nuclei. Sourced from _oligo_programmes_FH.R so the spatial and
# snRNA panels cannot drift apart.
source(file.path(PROJ, "scripts", "_oligo_programmes_FH.R"))
.oligo <- oligo_programme_sets(readRDS(file.path(PROJ, "data", "curated_signatures.rds")))
MODULES <- c(MODULES, .oligo[c("Structural myelin", "Cholesterol (sterol arm)",
                               "Fatty-acid / sphingomyelin", "Galactolipid",
                               "Lipid uptake / salvage")])

# v5 — the proximity and severity analyses need, per
# spot, the astrocyte-reactive pair and the Fig-2 microglial programmes beside the existing
# modules, and the TREM2-arm genes as detection. Same thinning, same spots, same AddModuleScore.
# Caveat carried into the legend: FTH1/FTL (iron) and SLC2A3 (glycolysis) are not
MODULES <- c(MODULES, list(
  "Astrocyte reactive"       = c("GFAP","SERPINA3","VIM","CD44"),                 # Fig 3 pan-reactive pair + vim/CD44
  "Antigen presentation"     = c("CD74","HLA-DRA","HLA-DRB1","HLA-DMB"),          # Fig 2b programme
  "Glycolytic shift"         = c("PFKFB3","SLC2A3"),                              # Fig 2b programme
  "Iron handling"            = c("FTH1","FTL","HAMP","TMEM163")))                 # Fig 2b programme (pan-cellular members)
EXPORT_GENES <- unique(c(EXPORT_GENES, "TREM2","ITGAX","CST7","LGALS3","CLEC7A","PPARG",  # TREM2-arm detection per spot
                         "VIM","CD44","HLA-DRB1","HLA-DMB","PFKFB3","HAMP","TMEM163"))

# v6 (Supplementary Fig. 6 review) — the homeostatic microglial module is the reference row
# of the proximity panel: it separates "more microglia" (homeostatic score rises with the
# microglial fraction exactly as the activation modules do) from "a different microglial
# state" (activation rises, homeostatic does not). Canonical homeostatic set (Butovsky 2014;
# Keren-Shaul 2017 stage-1 loss set); members absent from the slide are dropped by the
# intersect() below and reported in the log.
MODULES <- c(MODULES, list(
  "Microglia homeostatic"    = c("P2RY12","CX3CR1","CSF1R","TMEM119","SALL1")))

# v4 — export the module member genes too, for the gene-level heatmap
# that sits beside 6_supplychain_spatial. Same genes, same thinning, same spots as the
# scores, so the heatmap and the maps cannot disagree.
EXPORT_GENES <- unique(c(EXPORT_GENES, unlist(MODULES, use.names = FALSE)))
cat(sprintf("per-spot genes to export: %d\n", length(EXPORT_GENES)))

COORD <- file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN", "spot_coords.csv")
stopifnot("MISSING spot_coords.csv" = file.exists(COORD))
# Cache key = gene sets + parameters + the mtime and size of the three inputs (integrated
# object, c2l proportions, spot coordinates), so a re-integrated object or a fresh
# deconvolution invalidates the tables and every downstream Fig-6 panel rebuilds from the
# current spots.
.stamp <- function(f) sprintf("%s@%d:%d", basename(f), as.integer(file.mtime(f)), as.integer(file.size(f)))
WM_Q <- 0.85   # WM-tail quantile of the oligodendrocyte proportion (used below for is_wm / depth); part of the cache key
KEY <- paste0("v6|T", DEPTH_T, "|floor", SPOT_FLOOR, "|wmq", WM_Q, "|",
              paste(sort(unlist(MARKERS)), collapse = ","), "|",
              paste(sort(unlist(MODULES)), collapse = ","), "|",
              paste(sort(EXPORT_GENES), collapse = ","), "|",
              paste(vapply(c(VIS, PROP, COORD), .stamp, ""), collapse = "|"))
KEYF <- file.path(PROJ, "data", "_visium_dm.key")
MARK_CSV <- file.path(TDIR, "visium_dm_marker_counts_FH.csv")
MOD_CSV  <- file.path(TDIR, "visium_dm_module_spots_FH.csv")
QC_CSV   <- file.path(TDIR, "visium_dm_qc_FH.csv")
GENE_CSV <- file.path(TDIR, "visium_dm_module_genes_FH.csv")   # module membership: the one list downstream scripts read
stale <- !file.exists(KEYF) || !identical(readLines(KEYF, warn = FALSE)[1], KEY) ||
         !all(file.exists(MARK_CSV, MOD_CSV, QC_CSV, GENE_CSV))
if (!stale) { cat("depth-matched tables are current (gene sets, parameters AND input mtimes match) — nothing to do.\n"); quit(save = "no") }
cat("gene set, parameters or an INPUT (object / c2l / coords) changed — rebuilding\n")

cat("Loading Visium object...\n")
v <- readRDS(VIS)
DefaultAssay(v) <- "Spatial"
cts <- GetAssayData(v, assay = "Spatial", layer = "counts")
md  <- v@meta.data
sid <- if ("sample_id" %in% colnames(md)) as.character(md$sample_id) else as.character(md$orig.ident)
cond <- ifelse(grepl("^NHD", sid), "NHD", "CON")
tot  <- Matrix::colSums(cts)
nfeat<- Matrix::colSums(cts > 0)

# ---- the depth audit, written out so the figure can show its own control ----
qc_pre <- data.frame(sample_id = sid, Condition = cond, total_umi = tot, n_gene = nfeat) %>%
  group_by(Condition, sample_id) %>%
  summarise(n_spot = n(), med_umi = median(total_umi), med_gene = median(n_gene),
            .groups = "drop") %>% mutate(stage = "raw")
cat("\n== spot depth BEFORE matching ==\n"); print(as.data.frame(qc_pre))

keep <- tot >= DEPTH_T
for (cc in c("CON","NHD"))
  cat(sprintf("retained at T=%d: %s %d/%d (%.0f%%)\n", DEPTH_T, cc,
              sum(keep & cond == cc), sum(cond == cc),
              100 * sum(keep & cond == cc) / sum(cond == cc)))
stopifnot("depth matching drops >40% of a condition — lower DEPTH_T" =
            all(tapply(keep, cond, mean) > 0.6))

ds <- Seurat::SampleUMI(cts[, keep, drop = FALSE], max.umi = DEPTH_T,
                        upsample = FALSE, verbose = FALSE)
rownames(ds) <- rownames(cts); colnames(ds) <- colnames(cts)[keep]
sid_k <- sid[keep]; cond_k <- cond[keep]

qc_post <- data.frame(sample_id = sid_k, Condition = cond_k,
                      n_gene = Matrix::colSums(ds > 0)) %>%
  group_by(Condition, sample_id) %>%
  summarise(n_spot = n(), med_umi = DEPTH_T, med_gene = median(n_gene), .groups = "drop") %>%
  mutate(stage = "depth-matched")
cat("\n== AFTER matching (this is the Methods number) ==\n"); print(as.data.frame(qc_post))
write.csv(bind_rows(qc_pre, qc_post), QC_CSV, row.names = FALSE)

# ---- 1. per-gene mean counts per 3,000 UMI, per section --------------------
mk <- unique(unlist(MARKERS))
present <- intersect(mk, rownames(ds))
if (length(setdiff(mk, present)))
  cat("markers absent from the slide (dropped):",
      paste(setdiff(mk, present), collapse = ", "), "\n")
grp <- rep(names(MARKERS), lengths(MARKERS))
names(grp) <- unlist(MARKERS)

mark <- lapply(unique(sid_k), function(s) {
  i <- sid_k == s
  data.frame(sample_id = s, Condition = unique(cond_k[i]), gene = present,
             mean_count = as.numeric(Matrix::rowMeans(ds[present, i, drop = FALSE])),
             pct_spots  = as.numeric(Matrix::rowMeans(ds[present, i, drop = FALSE] > 0)))
}) %>% bind_rows() %>% mutate(group = unname(grp[gene]))
write.csv(mark, MARK_CSV, row.names = FALSE)
cat(sprintf("\nWrote %s (%d genes x %d sections)\n", basename(MARK_CSV),
            length(present), length(unique(sid_k))))

# ---- 2. per-spot module scores on the thinned data + c2l content -----------
# AddModuleScore on the thinned matrix: its expression-matched control bins absorb
# most of a global depth term (verified: housekeeping null -0.119, random +0.008
# after matching), but the bins are pooled across conditions, so thinning first and
# scoring second is belt-and-braces rather than either alone.
vs <- CreateSeuratObject(counts = ds)
vs <- NormalizeData(vs, verbose = FALSE)
for (nm in names(MODULES)) {
  g <- intersect(MODULES[[nm]], rownames(vs))
  cat(sprintf("module %-26s %d/%d genes\n", nm, length(g), length(MODULES[[nm]])))
  vs <- AddModuleScore(vs, features = list(g), name = paste0("M_", make.names(nm)),
                       seed = 42, ctrl = 100)
}
# module membership as scored (present = on the slide and therefore inside the score);
# script 111 reads this for the signature-overlap statement, so the two cannot drift apart
modgenes <- bind_rows(lapply(names(MODULES), function(nm)
  data.frame(module = nm, gene = MODULES[[nm]], present = MODULES[[nm]] %in% rownames(vs))))
write.csv(modgenes, GENE_CSV, row.names = FALSE)
cat(sprintf("Wrote %s (%d modules, %d genes, %d present)\n", basename(GENE_CSV),
            length(MODULES), nrow(modgenes), sum(modgenes$present)))
sc <- vs@meta.data[, grep("^M_", colnames(vs@meta.data)), drop = FALSE]
colnames(sc) <- names(MODULES)

pc <- intersect(EXPORT_GENES, rownames(vs))
cat(sprintf("per-spot genes exported: %d/%d (%s)\n", length(pc), length(EXPORT_GENES),
            paste(setdiff(EXPORT_GENES, pc), collapse = ", ")))
pcm <- as.matrix(GetAssayData(vs, layer = "data")[pc, , drop = FALSE])

prop <- read.csv(PROP, check.names = FALSE, stringsAsFactors = FALSE)
rownames(prop) <- prop$spot_id
common <- intersect(colnames(ds), rownames(prop))
stopifnot("c2l join lost >1% of spots — barcode mismatch" =
            length(common) >= 0.99 * ncol(ds))

xy <- read.csv(COORD, stringsAsFactors = FALSE); rownames(xy) <- xy$cell   # COORD defined with the key above
out <- data.frame(spot_id = common, sample_id = sid_k[match(common, colnames(ds))],
                  x = xy[common, "x"], y = xy[common, "y"],
                  Condition = cond_k[match(common, colnames(ds))],
                  micro_prop = prop[common, "Micro-PVM"],
                  oligo_prop = prop[common, "Oligo"],
                  total_umi_pre = tot[keep][match(common, colnames(ds))],
                  sc[common, , drop = FALSE], check.names = FALSE)
for (g in pc) out[[paste0("gene_", g)]] <- pcm[g, common]

# ---- all eight deconvolved cell types (composition panel) -------------------
C2L_TYPES <- intersect(c("Astro","Endo","Micro-PVM","Neuron_Ex","Neuron_Inh",
                         "OPC","Oligo","Pericytes"), colnames(prop))
for (ct in C2L_TYPES) out[[paste0("c2l_", ct)]] <- prop[common, ct]
cat("c2l cell types exported:", paste(C2L_TYPES, collapse = ", "), "\n")

# ---- cortical-depth proxy ---------------------------------------------------
# There is no layer annotation on these sections, so depth is derived, not assumed:
# white matter is the high-oligodendrocyte tail of each section's own spots, and every
# spot's depth is its distance to the nearest WM spot, scaled by that section's own
# 95th-percentile distance. Per section, so a difference in section size or orientation
# cannot masquerade as a difference in depth. Spots in the WM tail get depth 0 and are
# excluded from the laminar panel (they are the reference, not a cortical layer).
wm_q <- WM_Q
out$is_wm <- FALSE; out$depth <- NA_real_
for (sm in unique(out$sample_id)) {
  ii  <- which(out$sample_id == sm)
  thr <- as.numeric(quantile(out$oligo_prop[ii], wm_q, na.rm = TRUE))
  wm  <- ii[out$oligo_prop[ii] >= thr]
  gm  <- setdiff(ii, wm)
  out$is_wm[wm] <- TRUE
  if (!length(wm) || !length(gm)) next
  d <- apply(cbind(out$x[gm], out$y[gm]), 1, function(p)
        min(sqrt((out$x[wm] - p[1])^2 + (out$y[wm] - p[2])^2)))
  out$depth[gm] <- d / as.numeric(quantile(d, 0.95, na.rm = TRUE))
  out$depth[wm] <- 0
  cat(sprintf("  %-14s WM spots %4d (oligo >= %.3f), GM %4d, median depth %.2f\n",
              sm, length(wm), thr, length(gm), median(out$depth[gm], na.rm = TRUE)))
}
write.csv(out, MOD_CSV, row.names = FALSE)
cat(sprintf("Wrote %s (%d spots)\n", basename(MOD_CSV), nrow(out)))

writeLines(KEY, KEYF)
writeLines(capture.output(sessionInfo()),
           file.path(LOGD, "65_visium_prep_depthmatched_FH_sessionInfo.txt"))
cat("=== DONE ===\n")
