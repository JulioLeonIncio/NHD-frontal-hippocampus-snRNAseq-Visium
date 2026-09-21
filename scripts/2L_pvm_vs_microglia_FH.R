#!/usr/bin/env Rscript
# =============================================================================
# 2L_pvm_vs_microglia_FH.R — PVM / border-associated macrophages vs parenchymal
# microglia: identity, proportions, and markers.  Figure 2.
# -----------------------------------------------------------------------------
# "before we had a very small cluster of PVM, that cluster was
# very different. I want to include that one and at least show PVM vs Microglia
# proportions, including cell markers."
#
# Why this SURVIVES where the STATE-CLUSTERS did not ([[nhd-microglia-do-not-subcluster]]):
# PVM/BAM vs parenchymal microglia is a lineage distinction with a large canonical
# marker set, not a transcriptional state. The diagnostic (2L0) showed the split is
# reproducible and well-supported, unlike the DAM-style state clusters:
#   * the same population separates at res 0.6, 0.8 and 1.0 (n = 161 / 209 / 208)
#   * CD163 is detected in ~50-58% of its nuclei vs ~5% in the main microglial mass
#   * in the resolution sweep it carried 45 credible markers and was the least
#     donor-confounded cluster of any solution (max condition frac 0.57, max lane
#     frac 0.26) — i.e. it is not a one-donor/one-lane artefact.
# The PVM cluster is therefore taken from the res-0.6 clustering, chosen
# DATA-DRIVEN (highest median PVM-minus-microglia identity lead), never by a
# hardcoded cluster index, which would silently break on any re-embed.
#
# Caveat built into the panel: MRC1 is detected in ~23% of all myeloid nuclei here
# and is near-uniform across clusters, so it is not discriminating in this dataset
# (ambient / broad low-level nuclear signal). The dotplot therefore leads with the
# markers that actually separate (CD163, F13A1, MS4A7, LYVE1, STAB1, DAB2) and MRC1
# is shown but must not be described as the discriminating marker.
#
# Oligo-doublet nuclei (Oligo module score > 0.5) are excluded before proportions.
#
# Outputs: figures/Figure_2/panels/2L_pvm_vs_micro_proportions.{png,pdf}
#          figures/Figure_2/panels/2M_pvm_vs_micro_markers.{png,pdf}
#          tables/micro_states/pvm_vs_micro_{proportions,markers,assignment}_FH.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(ggh4x)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_artifact_genes.R"))

CACHE <- file.path(PROJ, "data", "_cache_micro_states_FH.rds")
PANEL <- file.path(PROJ, "figures", "Figure_2", "panels")
TBL   <- file.path(PROJ, "tables", "micro_states")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TBL, LOGS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
stopifnot("MISSING micro cache — run 2J0_micro_states_cache_FH.R" = file.exists(CACHE))

RES_PVM   <- 0.6     # resolution at which the PVM population resolves (see 2L0)
OLIGO_CUT <- 0.5     # doublet cut, as established in 2J2

PVM_MK   <- c("MRC1","CD163","LYVE1","F13A1","MS4A7","STAB1","DAB2","MSR1","MAF")
MICRO_MK <- c("P2RY12","P2RY13","TMEM119","CX3CR1","GPR34","SELPLG","SLC2A5","OLFML3")
OLIGO_MK <- c("PLP1","MBP","MOG","MOBP","CNP","ST18","CLDN11","MAG")

o <- readRDS(CACHE)
DefaultAssay(o) <- "RNA"
if (inherits(o[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(o[["RNA"]], search = "counts")) > 1)
  o[["RNA"]] <- SeuratObject::JoinLayers(o[["RNA"]])
o <- NormalizeData(o, verbose = FALSE)

keep_mk <- function(g) g[g %in% rownames(o)]
PVM_MK <- keep_mk(PVM_MK); MICRO_MK <- keep_mk(MICRO_MK); OLIGO_MK <- keep_mk(OLIGO_MK)

# ---- 1. drop oligo doublets (they are not myeloid cells at all) -------------
o <- AddModuleScore(o, features = list(OLIGO_MK), name = "OLSC", seed = 42, ctrl = 50)
n0 <- ncol(o)
o  <- o[, o$OLSC1 <= OLIGO_CUT]
cat(sprintf("dropped %d oligo-doublet nuclei (%.1f%%); %d myeloid nuclei remain\n",
            n0 - ncol(o), 100 * (n0 - ncol(o)) / n0, ncol(o)))

# ---- 2. identify the PVM cluster, data-driven ------------------------------
o <- AddModuleScore(o, features = list(PVM_MK, MICRO_MK), name = "IDSC", seed = 42, ctrl = 50)
o$pvm_lead <- o$IDSC1 - o$IDSC2
GRAPH <- grep("_snn$", names(o@graphs), value = TRUE)[1]
o <- FindClusters(o, resolution = RES_PVM, graph.name = GRAPH, verbose = FALSE)
o$cl <- factor(as.integer(as.character(o$seurat_clusters)))

lead_by_cl <- o@meta.data %>% group_by(cl) %>%
  summarise(n = n(), med_lead = median(pvm_lead), .groups = "drop") %>%
  arrange(desc(med_lead))
cat("\n=== PVM-identity lead by cluster (res ", RES_PVM, ") ===\n", sep = "")
print(as.data.frame(lead_by_cl), digits = 3)
PVM_CL <- as.character(lead_by_cl$cl[1])
cat(sprintf("-> PVM cluster = %s (n = %d)\n", PVM_CL, lead_by_cl$n[1]))
stopifnot("top cluster does not look PVM-like — inspect before shipping" =
            lead_by_cl$med_lead[1] > lead_by_cl$med_lead[2] + 0.10)

o$myeloid_class <- factor(ifelse(as.character(o$cl) == PVM_CL, "PVM", "Microglia"),
                          levels = c("Microglia", "PVM"))
o$Condition <- factor(o$Condition, levels = c("CON","NHD"))
o$Region    <- factor(o$Region,    levels = REGION_ORDER)
cat("\n=== class x condition x region ===\n")
print(ftable(table(o$myeloid_class, o$Condition, o$Region)))

write.csv(o@meta.data %>% select(Condition, Region, SampleID, cl, pvm_lead, myeloid_class),
          file.path(TBL, "pvm_vs_micro_assignment_FH.csv"))

# ---- 3. PROPORTIONS --------------------------------------------------------
# Descriptive only: n = 1 donor per condition, so the fraction shift is the effect
# size and no significance test is drawn [[composition-panels-effectsize-not-pseudorep]].
prop <- o@meta.data %>%
  count(Region, Condition, myeloid_class, name = "n") %>%
  complete(Region, Condition, myeloid_class, fill = list(n = 0L)) %>%
  group_by(Region, Condition) %>%
  mutate(n_myeloid = sum(n), pct = 100 * n / n_myeloid) %>% ungroup()
pooled <- o@meta.data %>%
  count(Condition, myeloid_class, name = "n") %>%
  complete(Condition, myeloid_class, fill = list(n = 0L)) %>%
  group_by(Condition) %>%
  mutate(n_myeloid = sum(n), pct = 100 * n / n_myeloid) %>% ungroup() %>%
  mutate(Region = "Pooled", .before = 1)
prop_out <- bind_rows(prop %>% mutate(Region = as.character(Region)), pooled)
write.csv(prop_out, file.path(TBL, "pvm_vs_micro_proportions_FH.csv"), row.names = FALSE)
cat("\n=== PVM as % of myeloid nuclei ===\n")
print(as.data.frame(prop_out %>% filter(myeloid_class == "PVM") %>%
                      select(Region, Condition, n, n_myeloid, pct)), digits = 3)

p_prop <- ggplot(prop %>% filter(myeloid_class == "PVM"),
                 aes(x = Condition, y = pct, fill = Condition)) +
  geom_col(width = 0.62, linewidth = 0) +
  facet_wrap2(~ Region, nrow = 1, labeller = region_full_lowconf_labeller(),
              strip = strip_region_x(REGION_ORDER, clip = "off")) +
  scale_fill_manual(values = PAL_COND, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "PVM (% of myeloid nuclei)") +
  theme_pub(base_size = 7) +
  theme(panel.spacing = unit(6, "pt"), plot.margin = margin(6, 6, 3, 3))
ggsave(file.path(PANEL, "2L_pvm_vs_micro_proportions.png"), p_prop,
       width = 2.5, height = 1.9, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "2L_pvm_vs_micro_proportions.pdf"), p_prop,
       width = 2.5, height = 1.9, bg = "white")
cat("wrote 2L_pvm_vs_micro_proportions.{png,pdf}\n")

# ---- 4. marker DGE (PVM vs microglia) --------------------------------------
Idents(o) <- o$myeloid_class
# Removed from the test universe, not merely from the display shortlist — ribosomal
# and mitochondrial transcripts track nuclear depth/quality rather than identity, and
# PVM nuclei here are systematically deeper, so leaving them in lets a depth
# difference be reported as a marker (RPLP1 had reached the shortlist that way).
# Covers RPL*/RPS* (incl. RPLP1, RPSA), mito-ribosomal MRPL*/MRPS*, and the MT- genes.
RIBO_MITO_RE <- "^RP[LS]|^MRP[LS]|^MT-|^MTRNR|^MTATP|^MTCO|^MTND|^MTCYB"

# p-value underflow FLOOR (reviewer-facing rigor). Wilcoxon p-values from thousands of
# nuclei routinely underflow to a literal 0 in double precision (ITK, themis here) --
# that is a representation limit, not evidence of an infinitely small p, and a bare
# "0" in a supplementary table is indefensible (it also becomes -log10(0) = Inf the
# moment anyone plots these columns). Zeros are floored to the smallest normalised
# double and the substitution is COUNTED in the LOG so it is auditable rather than a
# silent fudge. The floored rows must be reported as p < 2.2e-308, never as p = 0.
P_FLOOR <- .Machine$double.xmin      # 2.225074e-308
floor_p <- function(df, tag) {
  for (cl in intersect(c("p_val", "p_val_adj"), names(df))) {
    n_uf <- sum(df[[cl]] == 0, na.rm = TRUE)
    if (n_uf > 0) {
      cat(sprintf("[%s] %d %s values underflowed to 0; floored to %g (report as < %g)\n",
                  tag, n_uf, cl, P_FLOOR, P_FLOOR))
      df[[cl]][df[[cl]] == 0] <- P_FLOOR
    }
  }
  df
}
ribo_mito <- grep(RIBO_MITO_RE, rownames(o), value = TRUE)
feat_use  <- setdiff(rownames(o), ribo_mito)
cat(sprintf("\nexcluded %d ribosomal/mitochondrial genes from the DGE universe (%d tested)\n",
            length(ribo_mito), length(feat_use)))
mk <- FindMarkers(o, ident.1 = "PVM", ident.2 = "Microglia", features = feat_use,
                  logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox",
                  verbose = FALSE) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(is_artifact = is_artifact(gene, "Micro-PVM"),
         pct_diff = pct.1 - pct.2) %>%
  arrange(desc(avg_log2FC)) %>%
  floor_p("PVM-vs-microglia")
write.csv(mk, file.path(TBL, "pvm_vs_micro_markers_FH.csv"), row.names = FALSE)
cat(sprintf("\nPVM-vs-microglia markers: %d rows (%d artifact-flagged)\n",
            nrow(mk), sum(mk$is_artifact)))

# Confound check: if PVM nuclei are systematically deeper, every gene gets a higher
# detection rate in them and "markers" are really complexity. Report it explicitly.
cat("\n=== complexity by class (detection confound check) ===\n")
print(o@meta.data %>% group_by(myeloid_class) %>%
        summarise(n = n(), med_nFeature = median(nFeature_RNA),
                  med_nCount = median(nCount_RNA)) %>% as.data.frame())

# Display MARKERS are chosen from the DGE ("pick better gene
# markers after DGE analysis"), not from a hand-written canonical list.
# Ranked by detection DIFFERENCE (pct.1 - pct.2), not avg_log2FC: log2FC is inflated
# for genes detected in a handful of nuclei (the earlier top-15-by-log2FC surfaced
# EIF1/NCL/RPSA/CANX/LAMTOR1 — translation-machinery genes riding on depth), whereas
# a detection difference demands the gene actually be present in a large share of one
# subtype and absent from the other. That is exactly what a display marker must do.
N_MK      <- 8
PCT_FLOOR <- 0.25
# PVM nuclei are DEEPER than microglial nuclei (median ~721 vs ~512 genes, printed
# above). Detection rate therefore rises with depth for any broadly-expressed gene,
# so ranking on pct_diff alone lets depth masquerade as identity — that is how FTH1,
# FTL, APOE and FOXP1 (all broadly expressed, pct in microglia 0.20-0.27) reached the
# top. A display marker must be absent from the other subtype, not merely better
# detected in this one, so the contrast is capped as well as floored.
PCT_CEIL_OTHER <- 0.20
# (ribosomal/mitochondrial genes are already gone — excluded from the test universe
# above — so no further ribo filter is needed here.)
pick <- function(direction) {
  d <- mk %>% filter(!is_artifact, p_val_adj < 0.05)
  if (direction == "PVM")
    d %>% filter(pct.1 >= PCT_FLOOR, pct.2 <= PCT_CEIL_OTHER) %>%
      slice_max(pct_diff, n = N_MK, with_ties = FALSE)
  else
    d %>% filter(pct.2 >= PCT_FLOOR, pct.1 <= PCT_CEIL_OTHER + 0.15) %>%
      slice_min(pct_diff, n = N_MK, with_ties = FALSE)
}
PVM_TOP   <- pick("PVM");   MICRO_TOP <- pick("Microglia")
cat("\n=== PVM-enriched display markers (by detection difference) ===\n")
print(as.data.frame(PVM_TOP   %>% select(gene, avg_log2FC, pct.1, pct.2, pct_diff, p_val_adj)), digits = 3)
cat("\n=== Microglia-enriched display markers (by detection difference) ===\n")
print(as.data.frame(MICRO_TOP %>% select(gene, avg_log2FC, pct.1, pct.2, pct_diff, p_val_adj)), digits = 3)

# ---- 4b. T-LYMPHOCYTES ("we are also missing the T lymphocyte
# markers there") ------------------------------------------------------------
# additive by design: the Microglia/PVM split derived above is left exactly as it was,
# so the validated 2L proportions and the two myeloid marker sets do not move. T cells
# simply were not in this script's cache (_cache_micro_states_FH.rds is myeloid-only);
# they live in _cache_immune_compartment_FH.rds, which carries all three immune classes
# (Microglia 2824 / PVM 162 / T-lymphocyte 60).
# T cells are a separate lineage, not part of the pooled Micro-PVM compartment the DGE
# runs on -- they are shown here only so the reader can see what was set aside and why.
IMM <- file.path(PROJ, "data", "_cache_immune_compartment_FH.rds")
stopifnot("MISSING immune-compartment cache" = file.exists(IMM))
oi <- readRDS(IMM)
DefaultAssay(oi) <- "RNA"
if (inherits(oi[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(oi[["RNA"]], search = "counts")) > 1)
  oi[["RNA"]] <- SeuratObject::JoinLayers(oi[["RNA"]])
oi <- NormalizeData(oi, verbose = FALSE)
.myc <- setNames(as.character(o$myeloid_class), colnames(o))
oi$class3 <- ifelse(colnames(oi) %in% names(.myc), unname(.myc[colnames(oi)]),
                    ifelse(as.character(oi$immune_class) == "T-lymphocyte",
                           "T-lymphocyte", NA_character_))
oi <- oi[, !is.na(oi$class3)]
oi$class3 <- factor(oi$class3, levels = c("Microglia", "PVM", "T-lymphocyte"))
cat("\n=== three-class counts for the marker dotplot ===\n"); print(table(oi$class3))
Idents(oi) <- "class3"
feat_T <- setdiff(rownames(oi), grep(RIBO_MITO_RE, rownames(oi), value = TRUE))
mkT <- FindMarkers(oi, ident.1 = "T-lymphocyte", features = feat_T,
                   logfc.threshold = 0.25, min.pct = 0.10, test.use = "wilcox",
                   verbose = FALSE) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(is_artifact = is_artifact(gene, "Micro-PVM"), pct_diff = pct.1 - pct.2) %>%
  floor_p("T-vs-myeloid")
write.csv(mkT, file.path(TBL, "tcell_vs_myeloid_markers_FH.csv"), row.names = FALSE)
# Same display rule as the myeloid markers: rank on detection difference with a floor in
# the target class and a ceiling in the others, so depth cannot masquerade as identity.
T_TOP <- mkT %>% filter(!is_artifact, p_val_adj < 0.05,
                        pct.1 >= PCT_FLOOR, pct.2 <= PCT_CEIL_OTHER) %>%
  slice_max(pct_diff, n = N_MK, with_ties = FALSE)
cat("\n=== T-lymphocyte-enriched display markers (by detection difference) ===\n")
print(as.data.frame(T_TOP %>% select(gene, avg_log2FC, pct.1, pct.2, pct_diff, p_val_adj)), digits = 3)
stopifnot("no T-lymphocyte display markers survived the filters" = nrow(T_TOP) > 0)

# ---- 5. marker DOTPLOT -----------------------------------------------------
MICRO_SHOW <- MICRO_TOP$gene
PVM_SHOW   <- PVM_TOP$gene
T_SHOW     <- T_TOP$gene
GENES <- c(MICRO_SHOW, PVM_SHOW, T_SHOW)
# Dotplot now built on the three-class object so the T row is real, not padded
expr  <- FetchData(oi, vars = GENES, layer = "data")
expr$class <- oi$class3
dot <- expr %>%
  pivot_longer(-class, names_to = "gene", values_to = "e") %>%
  group_by(class, gene) %>%
  summarise(pct_exp = 100 * mean(e > 0), avg_exp = mean(e), .groups = "drop") %>%
  # Colour = mean expression, not a z-score. With only two groups a per-gene z-score
  # is degenerate — it can only ever be +-0.707, so every dot saturates and the colour
  # conveys nothing but "higher in this row", exaggerating genes that barely differ
  # (GPR34, SELPLG). Mean log-normalised expression keeps magnitude meaningful.
  ungroup() %>%
  mutate(panel = factor(dplyr::case_when(gene %in% MICRO_SHOW ~ "Microglia-enriched",
                                         gene %in% PVM_SHOW   ~ "CD163+/F13A1+-enriched",
                                         TRUE                 ~ "T-lymphocyte-enriched"),
                        levels = c("Microglia-enriched", "CD163+/F13A1+-enriched",
                                   "T-lymphocyte-enriched")),
         gene  = factor(gene, levels = GENES),
         class = factor(class, levels = c("T-lymphocyte", "PVM", "Microglia")))

write.csv(dplyr::arrange(dot, panel, gene, class),
          file.path(TBL, "immune_marker_dotplot_FH.csv"), row.names = FALSE)

# apparent-type parity. House rule:
# on-page pt = declared pt x placement scale (placed_width/saved_width in the
# assembled figure); main-figure band 4.8-5.3 pt, target 5.1.
# 2M is placed 6.50 -> 4.098 in (scale 0.630) in SuppFig3 -- the WIDEST canvas in the
# figure and therefore the most shrunken: 7 pt read only 4.41 pt on the page, the
# coldest panel there. TYPE_F is >1 (pre-compensating a wide-rendered panel, per
# [[apparent-font-parity-precompensate-wide-panels]]): 7 x 1.16 = 8.12 pt -> 5.11 pt.
# Text only -- the dot size range and legend.key.size unit() are unchanged.
TYPE_F_2M <- 0.90

# House dot-matrix grammar. Matched by eye to the reference block-
# diagonal marker matrix, SuppFig1 panel f (70_suppfig1_NHD_QC_FH.R): shape-21 dots
# with a grey35 outline, a thin black panel border in place of axis lines, no
# gridlines, two-line legend titles, and a compact colourbar (0.28 x 1.4 cm, three
# breaks) so the guide stops dominating a 2.05 in tall panel.
# The colour scale deliberately stays sequential. This panel plots mean
# log-normalised expression, a one-sided quantity (>= 0) -- not a z-score -- so panel
# f's diverging steel/rose ramp would imply a signed contrast that does not exist
# here, and its legend title "scaled mean expr" would be untrue. Hence "mean\nexpr"
# and a white -> brick ramp.
EXP_MAX <- max(dot$avg_exp)
EXP_BR  <- if (EXP_MAX >= 2) c(0, 1, 2) else signif(seq(0, EXP_MAX, length.out = 3), 2)
cat(sprintf("== 2M colour scale: mean expr [%.2f, %.2f]; breaks %s ==\n",
            min(dot$avg_exp), EXP_MAX, paste(EXP_BR, collapse = " / ")))

p_dot <- ggplot(dot, aes(x = gene, y = class)) +
  geom_point(aes(size = pct_exp, fill = avg_exp), shape = 21,
             colour = "grey35", stroke = 0.25) +
  # Row labels use the CD163+/F13A1+ marker name, matching SuppFig4 panel c
  # (2N). Bare "PVM" asserts a perivascular origin the data cannot establish; only the
  # joint-class label "Micro-PVM" is sanctioned. The facet titles stay
  # "PVM-enriched" -> "CD163+/F13A1+-enriched" for the same reason.
  scale_y_discrete(labels = IMM_LABEL) +
  # two-line strip titles: at 4.875 in the one-line "CD163+/F13A1+-enriched" overran its strip
  facet_grid2(~ panel, scales = "free_x", space = "free_x",
              labeller = as_labeller(c("Microglia-enriched" = "Microglia\nenriched",
                                       "CD163+/F13A1+-enriched" = "CD163+/F13A1+\nenriched",
                                       "T-lymphocyte-enriched" = "T-lymphocyte\nenriched"))) +
  scale_size_continuous(range = c(0.3, 3.2), name = "% cells\nexpressing",
                        breaks = c(10, 25, 50, 75)) +
  scale_fill_gradient(low = "grey95", high = "#B3132B", name = "mean\nexpr",
                      limits = c(0, EXP_MAX), breaks = EXP_BR) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(0.28, "cm"),
                                barheight = unit(1.4, "cm")),
         size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 7 * TYPE_F_2M) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                   face = "italic", colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.25),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.grid = element_blank(),
        panel.spacing.x = unit(1.5, "pt"),
        strip.background = element_rect(fill = "grey94", colour = NA),   # grey90 -> grey94: one strip tone across Supp Fig 3
        strip.text = element_text(colour = "black", lineheight = 0.82),
        strip.placement = "outside", strip.clip = "off",
        legend.position = "right", legend.box = "vertical",
        legend.key.size = unit(0.28, "cm"),
        legend.title = element_text(size = 6.2 * TYPE_F_2M, lineheight = 0.9),
        legend.text  = element_text(size = 6 * TYPE_F_2M),
        plot.margin = margin(6, 4, 3, 3))
W_2M <- 6.5 * 0.75
ggsave(file.path(PANEL, "2M_pvm_vs_micro_markers.png"), p_dot,
       width = W_2M, height = 2.05, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "2M_pvm_vs_micro_markers.pdf"), p_dot,
       width = W_2M, height = 2.05, bg = "white")
cat("wrote 2M_pvm_vs_micro_markers.{png,pdf}\n")

# ---- 6. UMAP with PVM shown as its own subtype -----------------------------
# Same embedding as the other Fig-2 microglia UMAPs (it comes from the shared cache),
# so this panel overlays the existing ones exactly. PVM are drawn last so that 162
# nuclei are not buried under 2,824 microglia. Point style matches 2A (Fig-2 sibling
# parity): size 0.55 with the thin dot-ring, which is the convention for the smaller
# per-cell-type UMAPs — the ring was only dropped on the dense whole-atlas UMAP.
PAL_MYELOID <- c(Microglia = unname(PAL_CELLTYPE[["Micro-PVM"]]), PVM = "#6D4C41")
ud <- as.data.frame(Embeddings(o, "umap"))[, 1:2]
colnames(ud) <- c("UMAP_1", "UMAP_2")
ud$class <- o$myeloid_class
ud <- rbind(ud[ud$class == "Microglia", ], ud[ud$class == "PVM", ])   # PVM on top

# Axis glyph outside the data (Fig-1 lesson: anchoring it inside the data range puts
# it on top of cells the moment the embedding changes).
xr <- range(ud$UMAP_1); yr <- range(ud$UMAP_2)
dx <- diff(xr); dy <- diff(yr)
gx <- xr[1] + 0.005 * dx; gy <- yr[1] - 0.135 * dy
p_umap <- ggplot(ud, aes(UMAP_1, UMAP_2, fill = class)) +
  geom_point(size = 0.55, alpha = 0.80, stroke = 0.04, shape = 21, colour = "grey20") +
  scale_fill_manual(values = PAL_MYELOID, name = NULL) +
  annotate("segment", x = gx, xend = gx + 0.11 * dx, y = gy, yend = gy,
           arrow = arrow(length = unit(1.3, "mm"), type = "closed"),
           linewidth = 0.3, colour = "grey25") +
  annotate("segment", x = gx, xend = gx, y = gy, yend = gy + 0.11 * dy,
           arrow = arrow(length = unit(1.3, "mm"), type = "closed"),
           linewidth = 0.3, colour = "grey25") +
  annotate("text", x = gx + 0.012 * dx, y = gy - 0.012 * dy, label = "UMAP 1",
           size = 1.9, colour = "grey25", hjust = 0, vjust = 1) +
  annotate("text", x = gx - 0.055 * dx, y = gy + 0.012 * dy, label = "UMAP 2",
           size = 1.9, colour = "grey25", hjust = 0, vjust = 0.5, angle = 90) +
  coord_cartesian(xlim = c(xr[1] - 0.11 * dx, xr[2] + 0.02 * dx),
                  ylim = c(yr[1] - 0.20 * dy, yr[2] + 0.02 * dy), clip = "off") +
  guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20",
                                                 stroke = 0.3, size = 2.2, alpha = 1))) +
  theme_pub(base_size = 8) +
  # Blank axis.title.x/.y explicitly: theme_pub() sets those two children directly,
  # so blanking only the parent `axis.title` leaves them in place and the default
  # "UMAP_1"/"UMAP_2" titles print alongside the corner glyph.
  theme(axis.title = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = c(0.84, 0.94),
        legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
        legend.key = element_blank(), legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 6.6),
        plot.margin = margin(3, 4, 2, 2))
ggsave(file.path(PANEL, "2N_myeloid_UMAP_PVM.png"), p_umap,
       width = 2.6, height = 2.5, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL, "2N_myeloid_UMAP_PVM.pdf"), p_umap,
       width = 2.6, height = 2.5, bg = "white")
cat("wrote 2N_myeloid_UMAP_PVM.{png,pdf}\n")

writeLines(c(sprintf("run: %s", format(Sys.time())), capture.output(sessionInfo())),
           file.path(LOGS, "2L_pvm_vs_microglia_FH_sessionInfo.txt"))
cat("\n=== DONE: 2L_pvm_vs_microglia_FH ===\n")
