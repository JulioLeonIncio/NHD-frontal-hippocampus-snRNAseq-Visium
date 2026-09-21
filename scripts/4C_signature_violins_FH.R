#!/usr/bin/env Rscript
# =============================================================================
# 4C_signature_violins_FH.R — Figure 4 panel f: neuronal programme-module signature violins, frontal + hippocampus. Faithful port of manuscript_7fig/
# scripts/_run_5C_only.R to the FH rebuild (2-region).
# -----------------------------------------------------------------------------
# Scores 4 core neuronal-integrity programmes (all expected down in NHD deep-layer
# Ex) + IEG as an explicit snRNA capture-floor negative control, per nucleus, in
# excitatory neurons.  Layout: modules = facet COLUMNS (per-module independent y),
# regions = facet ROWS (pale PAL_REGION banners).  Cliff's delta is the sole
# annotation; no per-nucleus Wilcoxon stars (pseudoreplication at N=1 donor).
#
# Why PROGRAMMES (not genes): the gene-level NHD Neuron_Ex landscape is up-dominant
# (nominal padj<0.05 ~63% up both regions — see 4B); only the coordinated
# PROGRAMMES (synaptic/GluR/OXPHOS) fall.  This panel is the programme-down (not
# gene-down) evidence.
#
# Module re-validation vs the
# FH Neuron_Ex MAST tables — coverage logged by score_one():
#   * Presynaptic (SNARE)   19/19 present, down both regions   -> kept.
#   * Postsynaptic (PSD)    13/17 present (NLGN2/3, DLG2/3 absent), down -> kept.
#   * Ionotropic GluR        8/8  present, down                  -> kept.
#   * OXPHOS / etc           6/23 present (nuclear-encoded subset; MT-genome absent),
#                            down — module-average absorbs sparse coverage           -> kept.
#   * IEG (capture floor)    kept as an explicit negative control (under-captured in
#                            nuclei), never asserted as biology.
#   nothing dropped as biology; absent members simply never score (module = mean of
#   present genes, requires >=3).
#
# Cache: data/_cache_neuron_scores_FH.rds (mtime + module-content fingerprint guard).
# Output: figures/Figure_5/panels/4C_signature_violins.{pdf,png}
#         + tables/4C_signature_violin_stats_FH.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(Matrix); library(ggh4x); library(scales)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("Could not resolve PROJ root" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, PAL_REGION, REGION_FULL, REGION_ORDER

PANEL <- file.path(PROJ, "figures", "Figure_5", "panels")
DATA  <- file.path(PROJ, "data")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, DATA, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

CON_COL <- "#82B2D6"; NHD_COL <- "#D6837A"   # Fig-2 violin pale pair
REGIONS <- REGION_ORDER                      # c("Frontal","Hippo")

base_sigs <- readRDS(file.path(DATA, "curated_signatures.rds"))
stopifnot("curated_signatures.rds missing IEG_core" = "IEG_core" %in% names(base_sigs))

INCLUDE_IEG <- TRUE
neuron_sigs <- list(
  Presynaptic = c("SYN1","SYN2","SYN3","SYP","SYT1","STX1A","STX1B","SNAP25",
                  "VAMP2","RAB3A","SYNGR1","SLC17A7","STXBP1","CPLX1","SV2A","RIMS1","RIMS2","BSN","PCLO"),
  Postsynaptic = c("DLG4","HOMER1","SHANK1","SHANK2","SHANK3","CAMK2A","CAMK2B",
                   "DLGAP1","NLGN1","NLGN2","NLGN3","NRGN","DLG2","DLG3","SYNGAP1","GRIP1","BEGAIN"),
  GluR_ionotropic = c("GRIN1","GRIN2A","GRIN2B","GRIA1","GRIA2","GRIA3","GRIA4","GRIK2"),
  OxPhos = c("NDUFA1","NDUFA2","NDUFA4","NDUFA8","NDUFB1","NDUFB2","UQCRC1","UQCRC2",
             "COX4I1","COX5A","COX6C","ATP5F1A","ATP5F1B","ATP5F1C","ATP5MC2","ATP5PO",
             "NDUFS1","NDUFS2","NDUFV1","NDUFA9","SDHA","UQCRB","COX7A2"))
if (INCLUDE_IEG) neuron_sigs$IEG <- base_sigs$IEG_core

cache_v <- file.path(DATA, "_cache_neuron_scores_FH.rds")
ATLAS   <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
stopifnot("MISSING atlas NHD_FH_harmony.rds" = file.exists(ATLAS))

cache_ok <- FALSE
if (file.exists(cache_v) && file.mtime(cache_v) >= file.mtime(ATLAS)) {
  sc_df <- readRDS(cache_v)
  if (identical(attr(sc_df, "sigs"), neuron_sigs)) {
    cache_ok <- TRUE; cat("== cache == loaded FRESH score cache (module defs match).\n")
  } else cat("== cache == module defs DIFFER from current neuron_sigs -> rebuild.\n")
}
if (!cache_ok) {
  cat("== build == full Seurat load + scoring", length(neuron_sigs), "modules...\n")
  suppressPackageStartupMessages(library(Seurat))
  obj <- readRDS(ATLAS)
  DefaultAssay(obj) <- "RNA"
  if (length(SeuratObject::Layers(obj, assay = "RNA")) > 1) obj <- JoinLayers(obj, assay = "RNA")
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  neurons <- subset(obj, new_annotation %in% c("Neuron_Ex","Neuron_Inh"))
  rm(obj); gc()
  data_mat <- GetAssayData(neurons, assay = "RNA", layer = "data")   # sparse; never as.matrix()
  score_one <- function(genes) {
    g <- intersect(genes, rownames(data_mat))
    if (length(g) < 3) return(rep(NA_real_, ncol(data_mat)))   # too few genes: NA, not silent 0
    Matrix::colMeans(data_mat[g, , drop = FALSE])
  }
  cat("== build == per-module gene coverage (present / requested):\n")
  for (s in names(neuron_sigs)) {
    g <- intersect(neuron_sigs[[s]], rownames(data_mat))
    cat(sprintf("   %-16s %2d / %2d\n", s, length(g), length(neuron_sigs[[s]])))
  }
  score_mat <- sapply(neuron_sigs, score_one)
  sc_df <- data.frame(cell = colnames(neurons), Condition = neurons$Condition,
                      Region = neurons$Region, cell_type = neurons$new_annotation,
                      SampleID = neurons$SampleID, as.data.frame(score_mat), check.names = FALSE)
  attr(sc_df, "sigs") <- neuron_sigs
  tmp <- paste0(cache_v, ".tmp"); saveRDS(sc_df, tmp); file.rename(tmp, cache_v)
  cat("== build == wrote cache:", cache_v, "\n")
  rm(neurons, data_mat); gc()
}

# --- Filter to Ex only + valid cond/region ---
n0 <- nrow(sc_df)
sc_df <- sc_df %>%
  filter(cell_type == "Neuron_Ex", Condition %in% c("CON","NHD"), Region %in% REGIONS) %>%
  mutate(Condition = factor(Condition, levels = c("CON","NHD")),
         Region    = factor(Region, levels = REGIONS),
         donor     = sub("_.*$", "", cell))
cat(sprintf("== filter == Ex + valid cond/region: %d -> %d nuclei\n", n0, nrow(sc_df)))
stopifnot("No Ex nuclei after filtering" = nrow(sc_df) > 0)

sig_order  <- names(neuron_sigs)
sig_labels <- c(Presynaptic = "Presynaptic (SNARE)", Postsynaptic = "Postsynaptic (PSD)",
                GluR_ionotropic = "Ionotropic GluR", OxPhos = "OXPHOS / ETC",
                IEG = "IEG (capture floor)")[sig_order]

sc_long <- sc_df %>%
  pivot_longer(cols = all_of(sig_order), names_to = "signature", values_to = "score") %>%
  mutate(signature = factor(signature, levels = sig_order), Region = factor(Region, levels = REGIONS))

# --- Cliff's delta (rank-based, tie-correct midranks) ---
cliffs_d <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]; na <- length(a); nb <- length(b)
  if (!na || !nb) return(NA_real_)
  r <- rank(c(a, b)); U_b <- sum(r[(na + 1):(na + nb)]) - nb * (nb + 1) / 2
  2 * (U_b / (na * nb)) - 1
}
stat_df <- sc_long %>%
  group_by(signature, Region) %>%
  summarise(p = suppressWarnings(wilcox.test(score ~ Condition, exact = FALSE)$p.value),
            delta = cliffs_d(score[Condition == "CON"], score[Condition == "NHD"]),
            n_CON = sum(Condition == "CON" & !is.na(score)),
            n_NHD = sum(Condition == "NHD" & !is.na(score)),
            y_min = min(score, na.rm = TRUE), y_max = max(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(eff_tier = case_when(abs(delta) >= 0.474 ~ "large", abs(delta) >= 0.33 ~ "medium",
                              abs(delta) >= 0.15 ~ "small", TRUE ~ "negligible"),
         direction = case_when(is.na(delta) ~ NA_character_, abs(delta) < 0.15 ~ "negligible",
                               delta < 0 ~ "Down in NHD", TRUE ~ "Up in NHD"),
         q_BH = p.adjust(p, method = "BH"),
         # effect-gated grey "q*" mark only (q<0.05 & |delta|<0.15); no per-nucleus stars.
         star = ifelse(!is.na(q_BH) & q_BH < 0.05 & abs(delta) < 0.15, "q*", ""),
         delta_lab = sprintf("italic(delta) == '%+.2f'", delta),
         y_delta = y_max + 0.16 * (y_max - y_min))

# reviewer-facing rigor: p printed as exactly 0 is underflow, not true zero.
n_uf <- sum(stat_df$p == 0, na.rm = TRUE)
if (n_uf > 0) { cat(sprintf("== rigor == %d Wilcoxon p underflowed to 0; floored to %.3g\n", n_uf, .Machine$double.xmin))
  stat_df$p <- pmax(stat_df$p, .Machine$double.xmin) }

# --- persist stat table (legend keys to source) ---
stat_out <- stat_df %>%
  transmute(module = unname(sig_labels[as.character(signature)]), module_key = as.character(signature),
            region = unname(REGION_FULL[as.character(Region)]), cliffs_delta = round(delta, 4),
            direction, eff_tier, q_star = star, n_CON, n_NHD, n_total = n_CON + n_NHD,
            wilcox_p = signif(p, 3), q_BH = signif(q_BH, 3)) %>%
  arrange(match(module_key, sig_order), match(region, REGION_FULL))
write.csv(stat_out, file.path(TDIR, "4C_signature_violin_stats_FH.csv"), row.names = FALSE)
cat("== table == wrote tables/4C_signature_violin_stats_FH.csv\n")

# --- Plot (horizontal: modules = columns, regions = rows) ---
# no grey low-conf banner. Every region row strip gets the
# normal pale region tone (PAL_REGION_PALE) + black text, so Hippocampus reads
# pale green + black like its siblings (volcano grid / Zhou / LR). The single-
# NHD-lane caveat lives in the legend, never on the plot.
reg_fills <- unname(PAL_REGION_PALE[REGIONS])
reg_txt   <- rep("black", length(REGIONS))
strip_mr <- ggh4x::strip_themed(
  clip = "off",
  background_x = ggh4x::elem_list_rect(fill = NA, colour = NA),
  text_x = ggh4x::elem_list_text(colour = "black", size = 7, face = "plain"),
  background_y = ggh4x::elem_list_rect(fill = reg_fills, colour = NA),
  text_y = ggh4x::elem_list_text(colour = reg_txt, size = 7, face = "plain"))

p4C <- ggplot(sc_long, aes(x = Condition, y = score, fill = Condition)) +
  geom_violin(scale = "width", linewidth = 0.25, alpha = 0.85, trim = FALSE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "grey25", linewidth = 0.3, fatten = 1.8) +
  # Cliff's delta = sole annotation + effect-gated grey "q*" (no per-nucleus stars).
  geom_text(data = subset(stat_df, star == "q*"), aes(x = 1.5, y = y_delta, label = star),
            inherit.aes = FALSE, size = 2.0, vjust = 0.4, color = "grey50") +
  geom_text(data = stat_df, aes(x = 0.42, y = y_delta, label = delta_lab), inherit.aes = FALSE,
            parse = TRUE, hjust = 0, vjust = 0, size = 1.95, color = "black") +
  facet_grid2(Region ~ signature, scales = "free_y", independent = "y",
              labeller = labeller(signature = as_labeller(sig_labels),
                                  Region = region_full_lowconf_labeller(low_conf)),
              strip = strip_mr) +
  scale_fill_manual(values = c(CON = CON_COL, NHD = NHD_COL), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.28))) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_pub(base_size = 7.5) +
  theme(plot.title = element_blank(),
        strip.text.x = element_text(size = 7, face = "plain", margin = margin(b = 2, t = 1)),
        strip.text.y = element_text(size = 7, angle = -90, face = "plain", margin = margin(l = 3, r = 1)),
        strip.placement = "outside", strip.clip = "off",
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 5.5), axis.title.y = element_text(size = 7),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x = unit(0.5, "lines"), panel.spacing.y = unit(0.3, "lines"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(-2, 0, 0, 0), legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 7), plot.margin = margin(6, 4, 4, 4))

n_mod   <- length(sig_order)
panel_w <- (1.15 * n_mod + 1.6) * 0.9
panel_h <- 2.6                       # 2 region rows
BN <- "4C_signature_violins"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p4C, width = panel_w, height = panel_h, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p4C, width = panel_w, height = panel_h, dpi = 600, device = ragg::agg_png)

writeLines(c(sprintf("built: %s", Sys.time()),
             sprintf("modules: %s", paste(sig_order, collapse = ", ")),
             sprintf("INCLUDE_IEG: %s", INCLUDE_IEG), "", capture.output(sessionInfo())),
           file.path(LOGS, "4C_signature_violins_FH.sessionInfo.txt"))

cat(sprintf("\nWrote 4C signature violins (%d modules, Ex only, HORIZONTAL, 2 regions).\n", n_mod))
cat("\n=== delta summary (module x region) ===\n")
print(stat_out %>% select(module_key, region, cliffs_delta, direction, eff_tier, q_star, n_CON, n_NHD), row.names = FALSE)
cat("\n=== DONE ===\n", file = stderr())
