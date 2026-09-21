#!/usr/bin/env Rscript
# =============================================================================
# 4_micro_neuron_subtype_LR_FH.R — F5g_micro_neuron_subtype_LR_FH.R — Figure 4 (Part B): subtype-resolved micro->neuron receptor readout, frontal + hippocampus.
# -----------------------------------------------------------------------------
# LIMITATION (stated, reviewer-facing): the FH withNeurons_v2 CellChat objects
# were fit at major-class resolution (Neuron_Ex / Neuron_Inh), so a subtype-level
# CellChat (per-subtype comm prob) is not available. This panel therefore resolves
# the micro->neuron axis at subtype resolution via receptor expression on neuron
# SUBTYPES (the receiver side of the LR pair), using the neuron subtype map
# (data/neuron_subtype_map_FH.rds; cortical Azimuth subclass in Frontal, de-novo
# label in Hippocampus) + the atlas. This is the honest subtype-resolvable readout
# of the micro->neuron LR axis; ligand-side comm prob stays at the class level
# (panels 4E). Noted in the log and intended for the figure legend.
#
# Receptors shown = the micro->neuron ligand receptors: SPP1 integrins
# (ITGAV/ITGA4/8/9-ITGB1) + GAS6/PROS1 TAM (TYRO3 down / MERTK / AXL undetected)
# + GRN-SORT1 + SEMA4D-PLXNB2 + ENTPD1-ADORA1. Direction = NHD-CON avg-expr.
# Valence split honoured: SPP1-integrin and GAS6-TAM are separate receptor blocks.
# No causal verbs.
#
# Cache: data/_cache_Fig4_subtype_micro_receptor_FH.rds (receptor-set + atlas-mtime
# invalidated). Output: figures/Figure_5/panels/F5g_micro_neuron_subtype_LR.{pdf,png}
# -----------------------------------------------------------------------------
# 20% SLIMMER + apparent-type parity. This panel now sits in the
# bottom-right slot of the assembled Fig5 (it replaced the astrocyte l-r panel there).
# Placed on-page width is 3.18 in; at 5.12 in saved that is a placement scale of 0.621,
# so declared type = apparent_target / 0.621 with the figure-wide targets:
#   tick / legend 5.2 pt | legend title 5.4 | axis title & strip 6.0 | in-panel 5.4.
# Deliberate deviation from the house full-name rule ([[nhd-canonical-region-order]]):
# "Hippocampus" is shortened to "Hippo" here, matching its sibling astro panel, because
# the hippocampal facet carries only 3 subtypes and the full word is wider than the data
# it labels. The legend must spell the region out once.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(patchwork)
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(tibble); library(ggh4x); library(scales); library(forcats)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
source(file.path(PROJ, "scripts", "_neuron_subclass_FH.R"))   # CTX_SUBCLASSES (fixed x-axis order)

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
SMAP  <- file.path(PROJ, "data", "neuron_subtype_map_FH.rds")
PANEL <- file.path(PROJ, "figures", "Figure_5", "panels")
DDIR  <- file.path(PROJ, "data")
TDIR  <- file.path(PROJ, "tables")
CACHE <- file.path(DDIR, "_cache_Fig4_subtype_micro_receptor_FH.rds")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, DDIR, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING atlas" = file.exists(ATLAS), "MISSING subtype map" = file.exists(SMAP))

REGIONS <- c("Frontal","Hippo")
LOW_CONF_REGIONS <- "Hippo"
MIN_CELLS <- 20          # subtype x cond floor for a stable mean (else dropped + logged)

# micro->neuron ligand receptors, grouped by ligand family (left strip). AXL kept
# to show its absence (uncoupled GAS6 arm).
RECEPTORS_BY_PATH <- list(
  "SPP1"  = c("ITGAV","ITGB1","ITGA4","ITGA8","ITGA9"),
  "GAS6/PROS1 (TAM)" = c("TYRO3","MERTK","AXL"),
  "GRN"   = c("SORT1"),
  "SEMA4D"= c("PLXNB2"),
  "ENTPD1"= c("ADORA1"))
ALWAYS_KEEP <- c("AXL")
REQ_REC <- sort(unique(unlist(RECEPTORS_BY_PATH)))
rec_sig <- paste0(paste(REQ_REC, collapse = ","), "|class=subtype_majority_v3")   # key bump (class fix; v3 adds minority_tbl attr)

.rebuild <- TRUE
if (file.exists(CACHE)) {
  .c <- readRDS(CACHE)
  .rebuild <- !identical(attr(.c, "rec_request"), rec_sig) ||
              file.mtime(CACHE) < file.mtime(ATLAS) ||
              file.mtime(CACHE) < file.mtime(SMAP)
  if (.rebuild) cat("Cache stale (receptor set / atlas / subtype-map changed) -> rebuild\n")
}
if (.rebuild) {
  cat("Building subtype receptor-expression cache...\n")
  smap <- readRDS(SMAP)   # barcode, Region, neuron_subtype, Condition, source_annot
  obj <- readRDS(ATLAS)
  DefaultAssay(obj) <- "RNA"
  if (length(SeuratObject::Layers(obj, assay = "RNA")) > 1)
    obj <- JoinLayers(obj, assay = "RNA")
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)

  keep <- obj$new_annotation %in% c("Neuron_Ex","Neuron_Inh") & obj$Region %in% REGIONS
  sub  <- subset(obj, cells = colnames(obj)[keep]); rm(obj); gc()

  # attach subtype label by barcode (map is authoritative for subtype identity)
  st <- smap$neuron_subtype[match(colnames(sub), smap$barcode)]
  sub$neuron_subtype <- st
  n_mapped <- sum(!is.na(st))
  cat(sprintf("  neurons in atlas subset: %d; subtype-mapped: %d (%.1f%%)\n",
              ncol(sub), n_mapped, 100 * n_mapped / ncol(sub)))
  sub <- subset(sub, cells = colnames(sub)[!is.na(sub$neuron_subtype)]); gc()

  rec_in  <- intersect(REQ_REC, rownames(sub))
  missing <- setdiff(REQ_REC, rec_in)
  if (length(missing)) cat("  Receptors not in assay:", paste(missing, collapse = ", "), "\n")

  mat <- GetAssayData(sub, assay = "RNA", layer = "data")[rec_in, , drop = FALSE]
  meta <- sub@meta.data %>%
    transmute(subtype = as.character(neuron_subtype),
              class = ifelse(new_annotation == "Neuron_Ex", "Ex", "Inh"),
              Condition = as.character(Condition), Region = as.character(Region)) %>%
    rownames_to_column("cell")
  # Hippocampal de-novo subtypes span both
  # new_annotation classes (e.g. CA3 carries 7 Neuron_Inh-annotated nuclei), and grouping by
  # the per-nucleus class produced (a) a phantom second dot column for CA3 computed from 3-4
  # nuclei and (b) an n-gate that read the tiny stratum (`first(n_cells)`), dropping
  # Inh-CGE (VIP) and Inh-MGE despite 600+/470+ real nuclei. Each subtype now takes its
  # majority class (used only for Ex-then-Inh ordering) and every nucleus of the subtype
  # counts once. Frontal labels are single-class by construction (DFC ran on Neuron_Ex;
  # Azimuth Inh labels sit on Neuron_Inh), so their numbers are unchanged.
  maj <- meta %>% count(subtype, Region, class) %>% group_by(subtype, Region) %>%
    mutate(n_sub = sum(n)) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(subtype, Region, class_major = class, minority_n = n_sub - n)
  mixed <- maj[maj$minority_n > 0, ]
  if (nrow(mixed)) cat(sprintf("  subtype x region spanning both classes -> majority class used: %s\n",
      paste(sprintf("%s/%s (%s; %d minority nuclei)", mixed$subtype, mixed$Region, mixed$class_major, mixed$minority_n), collapse = "; ")))
  meta <- meta %>% left_join(maj, by = c("subtype","Region")) %>% mutate(class = class_major) %>% select(-class_major, -minority_n)
  minority_tbl <- maj %>% transmute(subtype, Region, minority_class_n = minority_n)   # carried into the cache
  long <- as.matrix(mat) %>% as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expr") %>%
    left_join(meta, by = "cell")
  summ <- long %>%
    group_by(gene, subtype, class, Region, Condition) %>%
    summarise(pct_exp = mean(expr > 0) * 100, avg_exp = mean(expr),
              n_cells = dplyr::n(), .groups = "drop")
  attr(summ, "rec_request") <- rec_sig
  attr(summ, "minority_tbl") <- minority_tbl   # nuclei of the non-majority class kept per subtype x region
  tmp <- paste0(CACHE, ".tmp"); saveRDS(summ, tmp); file.rename(tmp, CACHE)
  rm(sub, mat, long); gc()
  cat("  Cache saved:", basename(CACHE), "\n")
}
summ <- readRDS(CACHE)

# Require MIN_CELLS in both conditions for a subtype to give a stable delta.
# (n_cells is identical across the receptor genes within a subtype x Region x
# Condition, so summarise to one value before the wide pivot.)
n_by <- summ %>%
  group_by(subtype, Region, Condition) %>%
  summarise(n_cells = dplyr::first(n_cells), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = n_cells,
              values_fill = 0, names_prefix = "n_")
if (!"n_CON" %in% names(n_by)) n_by$n_CON <- 0
if (!"n_NHD" %in% names(n_by)) n_by$n_NHD <- 0
ok_sub <- n_by %>% filter(n_CON >= MIN_CELLS, n_NHD >= MIN_CELLS) %>%
  distinct(subtype, Region)
dropped <- n_by %>% filter(!(n_CON >= MIN_CELLS & n_NHD >= MIN_CELLS))
if (nrow(dropped))
  cat(sprintf("Dropped %d subtype x region below %d cells/cond (unstable delta): %s\n",
              nrow(dropped), MIN_CELLS,
              paste(sprintf("%s/%s(%d,%d)", dropped$subtype, dropped$Region,
                            dropped$n_CON, dropped$n_NHD), collapse = "; ")))
# Persist the per-cell census + drop decision so the legend's "dropped for
# n < 20" list is read from a file, not from a scrolled log.
minority_tbl <- attr(summ, "minority_tbl")
stopifnot("cache lacks minority_tbl — delete the cache and re-run" = !is.null(minority_tbl))
census <- n_by %>%
  left_join(minority_tbl, by = c("subtype","Region")) %>%
  mutate(minority_class_n = ifelse(is.na(minority_class_n), 0L, as.integer(minority_class_n)),
         min_cells = MIN_CELLS, kept = n_CON >= MIN_CELLS & n_NHD >= MIN_CELLS) %>%
  arrange(Region, desc(kept), subtype)
write.csv(census, file.path(TDIR, "figF5g_micro_neuron_subtype_LR_census.csv"), row.names = FALSE)
mino <- census %>% filter(minority_class_n > 0)
cat(sprintf("Minority-class nuclei kept per subtype (non-majority new_annotation class): %s\n",
            if (nrow(mino)) paste(sprintf("%s/%s=%d", mino$subtype, mino$Region, mino$minority_class_n), collapse = "; ") else "none"))

path_lookup <- stack(RECEPTORS_BY_PATH) %>%
  transmute(gene = as.character(values), pathway = as.character(ind))
wide <- summ %>%
  semi_join(ok_sub, by = c("subtype","Region")) %>%
  left_join(path_lookup, by = "gene") %>% filter(!is.na(pathway)) %>%
  pivot_wider(id_cols = c(gene, pathway, subtype, class, Region),
              names_from = Condition, values_from = c(avg_exp, pct_exp),
              values_fill = list(avg_exp = 0, pct_exp = 0)) %>%
  mutate(delta_avg = avg_exp_NHD - avg_exp_CON,
         max_pct   = pmax(pct_exp_CON, pct_exp_NHD),
         max_avg   = pmax(avg_exp_CON, avg_exp_NHD))

# Order subtypes: fixed taxonomy order, not count order —
# Frontal = CTX_SUBCLASSES (helper; depth-ordered Ex then the shipped Inh order) and
# Hippo = HIP_LABEL_ORDER exactly as F5f1's legend (64_fig4_neuron_umap_FH.R), so this
# panel's x axis reads in the same order as the UMAP legends it sits beside.
HIP_LABEL_ORDER <- c("Inh-CGE (VIP)","Inh-CGE (LAMP5)","Inh-MGE","CA3","Subiculum","Unresolved")
.seen <- unique(summ$subtype)
.unk  <- setdiff(.seen, c(CTX_SUBCLASSES, HIP_LABEL_ORDER))
stopifnot("subtype(s) outside the fixed CTX/HIP orders" = length(.unk) == 0)
sub_order <- c(CTX_SUBCLASSES, HIP_LABEL_ORDER)
sub_order <- sub_order[sub_order %in% ok_sub$subtype]

wide <- wide %>%
  filter(gene %in% unlist(RECEPTORS_BY_PATH)) %>%
  # Keep genes expressed somewhere or ALWAYS_KEEP (AXL)
  group_by(gene) %>% mutate(any_expr = max(max_avg, na.rm = TRUE)) %>% ungroup() %>%
  filter(any_expr > 0.05 | gene %in% ALWAYS_KEEP) %>% select(-any_expr) %>%
  mutate(Region  = factor(Region, levels = REGIONS),
         subtype = factor(subtype, levels = sub_order),
         pathway = factor(pathway, levels = names(RECEPTORS_BY_PATH)),
         gene    = factor(gene, levels = unlist(RECEPTORS_BY_PATH)))

cap_abs <- as.numeric(quantile(abs(wide$delta_avg), 0.95, na.rm = TRUE))
cap_abs <- max(cap_abs, 0.05)
wide$delta_clip <- pmin(pmax(wide$delta_avg, -cap_abs), cap_abs)
cat(sprintf("Subtypes shown: %d; plot rows: %d; fill cap +/-%.3f (95th pct)\n",
            length(unique(wide$subtype)), nrow(wide), cap_abs))

# SPP1-integrin on frontal deep-layer subtypes (the requested callout)
cat("\n== SPP1-integrin delta on frontal deep-layer subtypes (L6 CT / L6 IT / L6b) ==\n")
deep <- wide %>% filter(Region == "Frontal", pathway == "SPP1",
                        subtype %in% c("L6 CT","L6 IT","L6b","L6 IT Car3")) %>%
  transmute(subtype, gene, delta = round(delta_avg, 3), max_pct = round(max_pct, 1)) %>%
  arrange(subtype, gene)
print(as.data.frame(deep), row.names = FALSE)

write.csv(wide %>% select(gene, pathway, subtype, class, Region,
                          avg_exp_CON, avg_exp_NHD, delta_avg, pct_exp_CON,
                          pct_exp_NHD, max_pct),
          file.path(TDIR, "figF5g_micro_neuron_subtype_LR_stats.csv"), row.names = FALSE)

DOWN_COL <- unname(PAL_DGE["Up in CON"]); UP_COL <- unname(PAL_DGE["Up in NHD"])
p <- ggplot(wide, aes(x = subtype, y = gene, size = max_pct, fill = delta_clip)) +
  geom_point(shape = 21, color = "grey25", stroke = 0.25) +
  facet_grid2(pathway ~ Region, scales = "free", space = "free", switch = "y",
              labeller = labeller(Region = ggplot2::as_labeller(
                c(Frontal = "Frontal", Hippo = "Hippo"))),
              strip = strip_region_x(REGIONS, fontsize = 9.7, low_conf = LOW_CONF_REGIONS,
                                     clip = "off")) +
  scale_fill_gradient2(low = DOWN_COL, mid = "grey96", high = UP_COL, midpoint = 0,
                       name = expression(Delta ~ "expr." ~ "(NHD" - "CON)"),
                       limits = c(-cap_abs, cap_abs), breaks = c(-cap_abs, 0, cap_abs), labels = function(x) sprintf("%.2f", x),   # unrounded breaks (rounding pushed ±cap outside the limits and dropped the end labels)
                       guide = guide_colorbar(order = 1, title.position = "left", title.vjust = 0.9,
                                              barheight = unit(0.26, "cm"),
                                              barwidth = unit(1.9, "cm"))) +
  scale_size_continuous(name = "% expr.", range = c(0.6, 3.2), limits = c(0, 100),
                        breaks = c(10, 30, 60),
                        guide = guide_legend(order = 2, title.position = "left", nrow = 1,
                                             override.aes = list(fill = "grey50"))) +
  # Sender declared on the y axis (propagated across every l-r panel).
  # The rows are ligand - receptor pairs, so the axis names both halves; the grey strips
  # beside them are CellChat pathway names, which are not always the ligand (TENM4 sits
  # under "ADGRL", named for the receptor; VSIR under "vista"; GAS6 under "GAS"). Naming
  # the axis for the pair — and the x axis for the receiver — lets the reader get sender
  # and receiver from the panel instead of the legend.
  # Never a unicode arrow here: base pdf() cannot encode U+2192 and it has already broken
  # a render in this project. ASCII hyphen only.
  labs(x = "Neuron subtype (receiver)", y = "Microglial ligand - receptor") +
  theme_pub(base_size = 12.1) +
  theme(axis.text.y = element_text(size = 8.4, face = "italic", color = "black"),
        axis.text.x = element_text(size = 8.4, color = "black", angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_text(size = 9.7, margin = margin(t = 3)),
        axis.title.y = element_text(size = 9.7),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x = unit(1.4, "lines"),   # room for clip='off' full region banners
        panel.spacing.y = unit(0.10, "lines"),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 9.7,
                                         face = "plain", color = "black",
                                         margin = margin(r = 2, l = 2)),
        strip.background.y = element_rect(fill = "grey97", color = "grey80", linewidth = 0.25),
        strip.placement = "outside",
        legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "center", legend.justification = "center",
        legend.margin = margin(t = 0, b = 0), legend.spacing.y = unit(0.02, "cm"),
        legend.box.spacing = unit(0.10, "cm"),
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 8.4), legend.title = element_text(size = 8.7),
        plot.margin = margin(6, 4, 4, 4))

W <- 5.12; H <- 4.02   # 20% taller (panel g)   # +0.15 in for the sender header
# Sender marker — width-aware so the arrow neither collides with the label nor
# drifts away from it (22_publication_theme_FH.R).
p <- with_sender(p, "Microglia", size = 3.41, panel_width = W, panel_height = H, arrow_x = 0.880)
BN <- "F5g_micro_neuron_subtype_LR"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nWrote %s.{pdf,png} + figF5g_micro_neuron_subtype_LR_stats.csv + _census.csv\n", BN))
cat("=== DONE ===\n", file = stderr())
