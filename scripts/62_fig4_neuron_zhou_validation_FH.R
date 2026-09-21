#!/usr/bin/env Rscript
# =============================================================================
# 62_fig4_neuron_zhou_validation_FH.R — 62_figF5c_neuron_zhou_validation_FH.R — Figure 4 panel i: cross-cohort generalization vs Zhou 2023, excitatory neurons, frontal + hippocampus.
# Ex sibling of 3_astro_zhou_validation_FH.R (same FH grammar).  atlas-free.
# -----------------------------------------------------------------------------
# Our OCC is gone in this rebuild, so this is cross-region GENERALIZATION (our
# Frontal / Hippocampus vs Zhou occipital) — an extra region difference on TOP of
# the cohort difference — never same-region "replication".  Two facets: Frontal |
# Hippocampus, each vs Zhou.
#   x = our avg_log2FC (NHD vs CON) in the facet's region
#   y = Zhou log2FoldChange (NHD vs CON, occipital)
#   plotted = MAST discovery set (is_disc) intersect Zhou-measured; grey = n.s. in
#             Zhou, black = Zhou-sig & concordant, orange = Zhou-sig & discordant.
#   headline (never
#             hand-typed): Spearman rho over all shared genes; direction-consistent %
#             over discovery vs the region-matched genome-wide BASELINE (binom p vs
#             baseline, not 0.5).
#   framing: direction-consistent / generalization, never "replicated".
#
# Honest NUMBERS (from ST13, reported in the log + panel): Neuron_Ex Frontal
# rho=+0.67, 100% direction-consistent but vs a high region-matched baseline (~76.5%)
# -> strong but the high baseline is disclosed on-panel; Hippocampus concordance
# (~80%) is near its baseline (~67%, n.s.) -> honestly weaker in hippocampus.
#
# No bold, no captions-in-panel, full region names, single-lane Hippocampus
#   grey "(low n)" figure-wide.  seed 42.
# Output: figures/Figure_5/panels/F5c_neuron_zhou_validation.{png,pdf}
# -----------------------------------------------------------------------------
# apparent-type parity PASS  ("adjust the font/panel size ratio
#   so everything looks more or less uniform, favouring the readability of the figure")
# Measured from the assembled Fig5.pdf: each panel is placed as a linked image with a
# known scale (parsed from the page's placement matrices). On-page type was 3.2-3.9 pt
# for most panels but 6.4-6.9 pt for the chord -- a 2x spread, and most of it below the
# print-legibility floor. Declared sizes here are therefore set to
#     declared = apparent_target / placement_scale
# with one set of apparent targets shared by every Figure-5 panel:
#     tick / legend text 5.2 pt | legend title 5.4 | axis title & facet strip 6.0
#     panel title 6.4 | in-panel gene labels & annotations 5.4
# canvas size and PIXEL DIMENSIONS are unchanged, so re-linking in Illustrator is a
# no-op: the frame keeps its box and only the type grows.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(ggh4x); library(ragg)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, REGION_FULL, strip_region_x, REGION_ORDER
source(file.path(PROJ, "scripts", "_artifact_genes.R"))           # is_artifact

CT      <- "Neuron_Ex"
REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
PANEL <- file.path(PROJ, "figures/Figure_5/panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ZHOU_CLS_COL <- PAL_ZHOU_CLS   # theme token (black / orange / grey); was hand-typed green / red
# single-lane Hippocampus grey "(low n)" banner figure-wide.
LOW_CONF <- "Hippo"

LOGF <- file.path(LOGS, "62_figF5c_neuron_zhou_validation_FH_provenance.txt")
logcon <- file(LOGF, open = "wt")
logln <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = logcon) }
logln("== 62_figF5c_neuron_zhou_validation_FH.R  run %s ==", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logln("PROJ = %s", PROJ)

# gene CSVs: fig4f_neuronEx_zhou_FH_{Frontal,Hippo}_genes.csv (Ex-specific set)
gene_files <- setNames(
  file.path(TDIR, sprintf("fig4f_neuronEx_zhou_FH_%s_genes.csv", REGIONS)), REGIONS)
st13_file  <- file.path(TDIR, "ST13_zhou_crosscohort_FH.csv")
for (f in c(gene_files, st13_file))
  if (!file.exists(f)) stop(sprintf("MISSING %s — run scripts/70_zhou_crosscohort_FH.R first", f))
for (f in c(gene_files, st13_file))
  logln("input %-44s  mtime %s", basename(f), format(file.mtime(f), "%Y-%m-%d %H:%M:%S"))

st13 <- read.csv(st13_file, stringsAsFactors = FALSE)
st13 <- st13[st13$cell_type == "Neuron_Ex", ]
stat_for <- function(rg) { r <- st13[st13$region == rg, ]; stopifnot(nrow(r) == 1); r }

load_region <- function(rg) {
  d <- read.csv(gene_files[[rg]], stringsAsFactors = FALSE)
  stopifnot(all(c("gene","our_lfc","zhou_lfc","zhou_padj","is_disc","zhou_sig",
                  "artifact","concordant") %in% names(d)))
  d %>% filter(is_disc) %>%
    mutate(region = factor(REGION_FULL[[rg]], levels = unname(REGION_FULL[REGIONS])),
           cls = factor(dplyr::case_when(
             zhou_sig &  concordant ~ "Zhou-sig, concordant",
             zhou_sig & !concordant ~ "Zhou-sig, discordant",
             TRUE                    ~ "n.s. in Zhou"), levels = names(ZHOU_CLS_COL)))
}
m <- bind_rows(lapply(REGIONS, load_region))
cat("== plotted (discovery) genes per region + concordance class ==\n")
print(as.data.frame(m %>% count(region, cls) %>% tidyr::pivot_wider(
  names_from = cls, values_from = n, values_fill = 0)))

annot_df <- bind_rows(lapply(REGIONS, function(rg) {
  s <- stat_for(rg)
  data.frame(region = factor(REGION_FULL[[rg]], levels = unname(REGION_FULL[REGIONS])),
             rho = s$rho_full, concord = s$concord_frac, baseline = s$baseline_concord,
             binom_p = s$binom_p, n_shared = s$n_shared, n_both = s$n_both_sig,
             stringsAsFactors = FALSE)
}))
for (rg in REGIONS) {
  s <- stat_for(rg)
  logln("%-8s: rho_full=%.2f | direction-consistent %.0f%% vs %.0f%% baseline (binom p %s) | n=%d | Zhou-both-sig %d",
        rg, s$rho_full, 100*s$concord_frac, 100*s$baseline_concord,
        ifelse(is.na(s$binom_p), "NA", sprintf("%.2g", s$binom_p)), s$n_shared, s$n_both_sig)
}

lim <- max(abs(c(m$our_lfc, m$zhou_lfc)), na.rm = TRUE) * 1.03

# labels: strongest Zhou-sig genes by joint magnitude (artifact excluded).
lab_df <- m %>% filter(zhou_sig, !artifact) %>% group_by(region) %>%
  slice_max(abs(our_lfc) + abs(zhou_lfc), n = 4, with_ties = FALSE) %>% ungroup()

fmt_p <- function(p) { if (is.na(p)) return("n.s.")
  if (p <= .Machine$double.xmin) return(sprintf("< %.0e", .Machine$double.xmin))
  if (p < 1e-3) sprintf("%.0e", p) else sprintf("%.2g", p) }
annot_df$l1 <- sprintf("rho == %.2f", annot_df$rho)
# Honest: lead with rho, then direction-consistency against its high baseline.
# SHORTENED for the one-column layout. The old line carried the baseline
# and the binomial p and n, which ran off the right edge of a 2 in panel — and a p-value
# printed inside a panel is caption text besides. The astro sibling's form is two short
# lines
annot_df$l2 <- sprintf("%.0f%% of %d both-sig\n(baseline %.0f%%)\nrho over n = %s",   # three short lines: the 2-in facet cannot hold the baseline clause on one line
                       100*annot_df$concord, annot_df$n_both,
                       100*annot_df$baseline,
                       trimws(format(annot_df$n_shared, big.mark = ",")))   # trimws: format() width-pads to a double space

p <- ggplot(m, aes(our_lfc, zhou_lfc)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_point(data = filter(m, cls == "n.s. in Zhou"), aes(colour = cls), size = 0.7, alpha = 0.45, stroke = 0) +
  geom_point(data = filter(m, cls != "n.s. in Zhou"), aes(colour = cls, shape = cls), size = 1.6, alpha = 0.9, stroke = 0.55) +
  geom_text_repel(data = lab_df, aes(label = gene), size = 1.90, fontface = "italic", colour = "black",
                  bg.color = "white", bg.r = 0.12, point.padding = 0.5, max.overlaps = Inf, force = 3, seed = 42, box.padding = 0.42, min.segment.length = 0,
                  ylim = c(-lim, 0.38 * lim),
                  segment.size = 0.2, segment.colour = "grey60") +
  geom_text(data = annot_df, aes(x = -lim, y = lim, label = l1), hjust = 0, vjust = 1, size = 2.05,
            colour = "black", parse = TRUE, inherit.aes = FALSE) +
  geom_text(data = annot_df, aes(x = -lim, y = lim - 0.17*lim, label = l2), hjust = 0, vjust = 1, size = 1.85,   # the baseline clause overran the frame
            colour = "black", inherit.aes = FALSE) +
  # Vertical / one column ("make the validation panel vertical, one
  # column as in the astrocyte figure"). REGIONS is Frontal -> Hippocampus and
  # facet_wrap fills by row, so ncol = 1 stacks them in the right order. This also fixes
  # the coord_equal problem the astro panel hit: with xlim = ylim the facets are square,
  # so a 2-across layout leaves dead space at the left and right of a wide box.
  ggh4x::facet_wrap2(~ region, ncol = 1,
                     strip = strip_region_x(REGIONS, low_conf = LOW_CONF, fontsize = 6.8, clip = "off")) +
  scale_colour_manual(values = ZHOU_CLS_COL, name = NULL, drop = FALSE,
                      guide = guide_legend(ncol = 1, byrow = TRUE,
                                           override.aes = list(size = 2.2, alpha = 1, stroke = 0.7))) +
  scale_shape_manual(values = SHAPE_ZHOU_CLS, name = NULL, drop = FALSE,
                     guide = guide_legend(ncol = 1, byrow = TRUE,
                                          override.aes = list(size = 2.2, alpha = 1, stroke = 0.7))) +
  scale_x_continuous(limits = c(-lim, lim)) + scale_y_continuous(limits = c(-lim, lim)) +
  coord_equal(clip = "off") +
  labs(x = expression("our excitatory-neuron log"[2]*"FC"),
       # the "(NHD vs CON, occipital)" parenthetical is caption text and clipped in a
       # one-column panel; it belongs in the legend (same fix as the astro sibling).
       y = expression("Zhou 2023 log"[2]*"FC")) +
  theme_pub(base_size = 6.3) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.line = element_blank(), strip.background = element_blank(),
        panel.spacing.x = unit(0.4, "lines"), legend.position = "bottom",
        legend.text = element_text(size = 4.8), legend.key.size = unit(0.22, "cm"),
        legend.spacing.x = unit(0.10, "cm"),
        plot.margin = margin(4, 6, 2, 4))
# Both panels are
# placed at ~61-62 % and their facet FRAMES are the same width on the page (the extra 0.3 in of this canvas is
# margin), so parity means the same on-page type as F3c, i.e. F3c's native sizes: base 8 x 0.79 = 6.3 pt,
# strip 7 -> 6.8, gene labels 1.9 mm, annotations 2.05 / 1.85 mm, legend 4.8 pt, key 0.22 cm. Verified by eye
# on a side-by-side crop of the two composites at equal magnification (first pass at 7.3 pt still read larger).
W <- 2.35; H <- 4.00   # portrait, one column
ragg::agg_png(file.path(PANEL, "F5c_neuron_zhou_validation.png"), width = W, height = H, units = "in", res = 600)
print(p); invisible(dev.off())
ggsave(file.path(PANEL, "F5c_neuron_zhou_validation.pdf"), p, width = W, height = H, useDingbats = FALSE)
logln("wrote F5c_neuron_zhou_validation.{png,pdf}")

cat("\n== sessionInfo ==\n", file = logcon)
capture.output(sessionInfo(), file = logcon, append = TRUE)
close(logcon)
cat("provenance:", LOGF, "\n")
cat("=== DONE ===\n", file = stderr())
