#!/usr/bin/env Rscript
# =============================================================================
# 2H_zhou_validation_FH.R — Figure 2 panel h (CLOSES): cross-cohort validation vs Zhou 2023, frontal + hippocampus. Port of 58_fig2_zhou_validation.R.
# -----------------------------------------------------------------------------
# MICROGLIA (Micro-PVM) cross-cohort direction concordance vs the independent Zhou
# et al. 2023 NHD cohort (occipital; PMID 36658241).  Our OCC is gone in this
# rebuild, so this is cross-region GENERALIZATION (our Frontal / Hippocampus vs
# Zhou occipital) — an extra region difference on TOP of the cohort difference —
# not same-region replication.  Two facets: Frontal | Hippocampus, each vs Zhou.
#
#   x = our avg_log2FC (NHD vs CON) in the facet's region
#   y = Zhou log2FoldChange (NHD vs CON, occipital)
#   genes = plotted = the MAST discovery set (is_disc) intersect Zhou-measured;
#           grey = n.s. in Zhou, black = Zhou-sig & concordant, orange = Zhou-sig &
#           discordant.  A discovery gene is never dropped for failing Zhou sig.
# Never
#           hand-typed): Spearman rho over all shared genes, direction-consistent
#           % over discovery vs the region-matched genome-wide baseline (binom p).
#   framing: direction-consistent / generalization, never "replicated".
#
# Honest narration of the magnitude-vs-direction split (represented, not hidden):
#   Frontal: weak-but-positive rho (+0.10), 59% concord vs 48.5% baseline.
#   Hippo:   rho slightly NEGATIVE (-0.07) but 86% sign-concordant vs 60.9%
#            baseline — the DIRECTIONS agree strongly while the magnitude ranking
#            does not.  Both numbers are printed on the panel so the split is
#            visible; the legend explains rho (magnitude) vs concord (direction).
#
# Atlas-free: reads the precomputed FH Zhou tables (tables/fig2d_micro_zhou_FH_*
# _genes.csv + ST13_zhou_crosscohort_FH.csv), built by scripts/70_zhou_crosscohort_FH.R.
# Output: 2_zhou_validation_<CELLTYPE>.{png,pdf} + provenance log.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(ggh4x); library(ragg)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, REGION_FULL, strip_region_x, REGION_ORDER
source(file.path(PROJ, "scripts", "_artifact_genes.R"))           # is_artifact (for label safety)

CT      <- "Micro-PVM"   # POOLED microglia+PVM; PVM = 5.4% of compartment
REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
PANEL <- file.path(PROJ, "figures/Figure_2/panels")

# Celltype selects which cross-cohort comparison is drawn; everything below — the
# classification, the annotation, the grammar, the sizing — is shared, so the added
# panels cannot drift in format from the microglial one they sit beside.
#   CELLTYPE=Micro-PVM (default) -> 2_zhou_validation_Micro-PVM  (SuppFig 3, panel b)
#   CELLTYPE=Oligo               -> 2_zhou_validation_Oligo
#   CELLTYPE=OPC                 -> 2_zhou_validation_OPC
CELLTYPE <- Sys.getenv("CELLTYPE", "Micro-PVM")
stopifnot("CELLTYPE must be Micro-PVM, Oligo or OPC" =
            CELLTYPE %in% c("Micro-PVM", "Oligo", "OPC"))
CT_LAB  <- c("Micro-PVM" = "Micro-PVM", "Oligo" = "oligodendrocytes", "OPC" = "OPC")[[CELLTYPE]]   # one class name figure-wide (panel a says Micro-PVM)
# The microglial panel used to be plain "2_zhou_validation", from when it was
# the only cross-cohort panel. Beside siblings that name their cell type, an unnamed one
# reads as the general case, and a filename that does not state its contents is what the
# publication-readiness check forbids. All three now carry the CELLTYPE token used in ST13,
# ST1 and this script's own parameter, so filename and table say the same word.
OUT_BN  <- paste0("2_zhou_validation_", CELLTYPE)
GENE_PAT <- if (CELLTYPE == "Micro-PVM") "fig2d_micro_zhou_FH_%s_genes.csv" else
              paste0("zhou_FH_", CELLTYPE, "_%s_genes.csv")
cat(sprintf("CELLTYPE = %s  -> %s.{png,pdf}\n", CELLTYPE, OUT_BN))
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# uniform Zhou concordance palette (black / orange / grey — PAL_ZHOU_CLS)
ZHOU_CLS_COL <- PAL_ZHOU_CLS   # theme token (black / orange / grey); was hand-typed green / red

# --- provenance log sink -----------------------------------------------------
LOGF <- file.path(LOGS, sprintf("2H_zhou_validation_FH_provenance%s.txt",
                                ifelse(CELLTYPE == "Micro-PVM", "", paste0("_", CELLTYPE))))
logcon <- file(LOGF, open = "wt")
logln <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = logcon) }
logln("== 2H_zhou_validation_FH.R  run %s ==", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logln("PROJ = %s", PROJ)

# --- load per-region Zhou gene tables + ST13 stats ---------------------------
gene_files <- setNames(file.path(TDIR, sprintf(GENE_PAT, REGIONS)), REGIONS)
st13_file  <- file.path(TDIR, "ST13_zhou_crosscohort_FH.csv")
for (f in c(gene_files, st13_file))
  if (!file.exists(f)) stop(sprintf("MISSING %s — run scripts/70_zhou_crosscohort_FH.R first", f))
for (f in c(gene_files, st13_file))
  logln("input %-44s  mtime %s", basename(f), format(file.mtime(f), "%Y-%m-%d %H:%M:%S"))

st13 <- read.csv(st13_file, stringsAsFactors = FALSE)
st13 <- st13[st13$cell_type == CELLTYPE, ]
stat_for <- function(rg) {
  r <- st13[st13$region == rg, ]
  stopifnot(nrow(r) == 1)
  r
}

# Plotted genes = discovery set (is_disc) per region; classify by Zhou-sig concordance.
load_region <- function(rg) {
  d <- read.csv(gene_files[[rg]], stringsAsFactors = FALSE)
  stopifnot(all(c("gene","our_lfc","zhou_lfc","zhou_padj","is_disc","zhou_sig",
                  "artifact","concordant") %in% names(d)))
  d %>% filter(is_disc) %>%
    mutate(region = factor(REGION_FULL[[rg]], levels = unname(REGION_FULL[REGIONS])),
           cls = factor(dplyr::case_when(
             zhou_sig &  concordant ~ "Zhou-sig, concordant",
             zhou_sig & !concordant ~ "Zhou-sig, discordant",
             TRUE                    ~ "n.s. in Zhou"),
             levels = names(ZHOU_CLS_COL)))
}
m <- bind_rows(lapply(REGIONS, load_region))
cat("== plotted (discovery) genes per region + concordance class ==\n")
print(as.data.frame(m %>% count(region, cls) %>% tidyr::pivot_wider(
  names_from = cls, values_from = n, values_fill = 0)))

# per-region ST13 headline (never hand-typed)
annot_df <- bind_rows(lapply(REGIONS, function(rg) {
  s <- stat_for(rg)
  data.frame(region = factor(REGION_FULL[[rg]], levels = unname(REGION_FULL[REGIONS])),
             rho = s$rho_full, concord = s$concord_frac, baseline = s$baseline_concord,
             binom_p = s$binom_p, n_shared = s$n_shared, n_both = s$n_both_sig,
             stringsAsFactors = FALSE)
}))
for (rg in REGIONS) {
  s <- stat_for(rg)
  logln("%-8s: rho_full=%.2f | direction-consistent %.0f%% vs %.0f%% baseline (binom p %.2g) | n=%d | Zhou-both-sig %d",
        rg, s$rho_full, 100*s$concord_frac, 100*s$baseline_concord, s$binom_p, s$n_shared, s$n_both_sig)
}

# --- panel (house style: coord_equal, shared symmetric limit; 2 region facets) -
# lim was max(|log2FC|), so one extreme gene (|log2FC| 14.3) stretched
# both axes to +-14.7 while the mass of the cloud sits inside +-5 -- about 88% of each
# facet rendered empty.
# Same rule as the astro sibling: the
# 97.5th percentile of the plotted values, symmetric, with a small pad.
lim <- as.numeric(stats::quantile(abs(c(m$our_lfc, m$zhou_lfc)), 0.975, na.rm = TRUE)) * 1.02
cat(sprintf("Zhou axis frame: +-%.2f  (max |log2FC| in the data was %.1f)\n",
            lim, max(abs(c(m$our_lfc, m$zhou_lfc)), na.rm = TRUE)))
# Choose the frame from the data, not by eye: report what each candidate would contain.
for (.q in c(0.90, 0.95, 0.975, 0.99)) {
  .l <- as.numeric(stats::quantile(abs(c(m$our_lfc, m$zhou_lfc)), .q, na.rm = TRUE)) * 1.02
  cat(sprintf("   frame +-%.2f (q%.3f) contains %.2f%% of points\n", .l, .q,
              100 * mean(abs(m$our_lfc) <= .l & abs(m$zhou_lfc) <= .l, na.rm = TRUE)))
}
# Labels: strongest Zhou-sig genes by joint magnitude, ambient-artifact excluded.
# Labels are chosen after the frame and restricted to genes inside it -- ranking by
# joint magnitude first would hand every label to the very outliers the frame excludes,
# and ggrepel would then draw a leader line to a point that is not on the panel.
lab_df <- m %>% filter(zhou_sig, !artifact,
                       abs(our_lfc) <= lim, abs(zhou_lfc) <= lim) %>%
  group_by(region) %>%
  slice_max(abs(our_lfc) + abs(zhou_lfc), n = 4, with_ties = FALSE) %>% ungroup()   # -> 4 per facet: the callouts live in the band below the cloud, and six ran into each other there
# NB lab_df$region carries the full region names (REGION_FULL), not REGIONS, so index
# the tally by its own levels -- keying it to REGIONS silently reported "Hippo NA".
cat("labelled genes inside the frame: ",
    paste(sprintf("%s %d", names(table(lab_df$region)), as.integer(table(lab_df$region))),
          collapse = "  |  "), "\n", sep = "")

# per-facet annotation text (2 lines). fmt_p never prints a literal 0.
fmt_p <- function(p) { if (is.na(p)) return("NA")
  if (p <= .Machine$double.xmin) return(sprintf("< %.0e", .Machine$double.xmin))
  if (p < 1e-3) sprintf("%.0e", p) else sprintf("%.2g", p) }
# This was a real mislabel, not a wording preference ====
# The annotation read "direction-consistent 75% vs 50% baseline ... n = 1,885", which
# puts a 72-GENE statistic beside an 1,885-GENE denominator and invites the reader to
# hear "75% of 1,885 genes agree". Checked against the source: ST13's `concord_frac`
# is n_concord / n_both_sig = 54/72 = 0.750 (Frontal) and 55/67 = 0.821 (Hippocampus),
# i.e. it is computed only over genes significant in both cohorts. `baseline_concord`
# is the separate 0.4997 / 0.5087 over all n_shared genes, and rho_full is over
# n_shared as well. Three different denominators were collapsed into one line.
# The astro and neuron Zhou panels already carry the corrected house wording
# ("X% of N genes sig. in both (baseline Y%)" + "rho over n = M"); this panel was the
# last one still on the old form. Direction still leads here;
# only the denominators are now stated.
# One headline template across the three placed Zhou
# panels — "N% of n both-sig (baseline B%)" + "rho over n". Direction still leads here
# ; the binomial p is caption text and lives in the legend / ST13.
# Same three-line
# headline as scripts 3 / 62 -- rho leads (plotmath, quoted to keep the trailing zero), then the direction
# statistic and the rho denominator on two plain lines.
annot_df$l1 <- sprintf("rho == '%.2f'", annot_df$rho)
# Dropped the explanatory parenthetical
# that is caption text
# and belongs in the figure legend, not baked into the panel (no-captions-in-panels).
# The on-plot annotation is now just the rho value.
# Rho has its own denominator (all shared genes), so it is stated with rho, exactly as
# the astro/neuron siblings do -- never left to inherit the line above it.
annot_df$l2 <- sprintf("%.0f%% of %d both-sig\n(baseline %.0f%%)\nrho over n = %s",   # three short lines: the one-column facet is 2 in wide
                       100*annot_df$concord, annot_df$n_both, 100*annot_df$baseline,
                       trimws(format(annot_df$n_shared, big.mark = ",")))

# apparent-type parity. House rule:
# on-page pt = declared pt x placement scale (placed_width/saved_width in the
# assembled figure); main-figure band 4.8-5.3 pt, target 5.1.
# Placed 5.40 -> 4.412 in (scale 0.817) in SuppFig3. axis.text was 8 pt = 6.54 pt on
# the page, the hottest panel in that figure; TYPE_F puts it at 5.07 pt.
# Exception, stated deliberately: the geom_text / geom_text_repel gene callouts and
# the headline annotation (1.9-2.0 mm = 5.41-5.69 pt) are not scaled. At this
# placement they already render 4.42-4.65 pt on the page -- the same on-page size as
# the identical 1.9 mm callouts in the main-figure sibling F1f_concordance_grid_MAST,
# which is placed at the same 0.817 scale. Riding them down by 0.78 would put them at
# 3.45 pt, below print legibility. Result: this panel's internal hierarchy now matches
# 2D's instead of being 1.48x wide.
TYPE_F <- 0.94

p <- ggplot(m, aes(our_lfc, zhou_lfc)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_point(data = filter(m, cls == "n.s. in Zhou"),
             aes(colour = cls), size = 0.7, alpha = 0.45, stroke = 0) +
  geom_point(data = filter(m, cls != "n.s. in Zhou"),
             aes(colour = cls, shape = cls), size = 1.6, alpha = 0.9, stroke = 0.55) +
  geom_text_repel(data = lab_df, aes(label = gene), size = 2.3, fontface = "italic",
                  colour = "black", bg.color = "white", bg.r = 0.12, point.padding = 0.5, max.overlaps = Inf, force = 4, seed = 42,
                  ylim = c(-Inf, -0.4), box.padding = 0.45,   # callouts below the cloud (which lives at y 0-2), never through it or the headline
                  min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey60") +
  # per-facet headline (top-left): l1 = direction statistic (leads, plain text);
  # l2 = rho with clarifying clause (plotmath -> base pdf safe).
  geom_text(data = annot_df, aes(x = -lim * 0.96, y = lim * 0.96, label = l1),
            hjust = 0, vjust = 1, size = 2.6, colour = "black", parse = TRUE, inherit.aes = FALSE) +
  geom_text(data = annot_df, aes(x = -lim * 0.96, y = lim * 0.78, label = l2),   # line 2 was touching the rho line at the larger type
            hjust = 0, vjust = 1, size = 2.3, colour = "black", inherit.aes = FALSE) +
  # Portrait, Frontal above Hippocampus, exactly as the Fig 3c / 5c siblings; the
  # GSEA dot plot (Supp Fig 3c) now sits beside it.
  ggh4x::facet_wrap2(~ region, ncol = 1,
                     strip = strip_region_x(REGIONS, clip = "off")) +
  scale_colour_manual(values = ZHOU_CLS_COL, name = NULL, drop = FALSE,
                      guide = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 2.2, alpha = 1, stroke = 0.7))) +   # one key per row: two keys no longer fit the 2.05-in canvas at the larger type
  scale_shape_manual(values = SHAPE_ZHOU_CLS, name = NULL, drop = FALSE,
                     guide = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 2.2, alpha = 1, stroke = 0.7))) +
  # The frame is set on the COORD, not on the scales. scale_*_continuous(limits=)
  # censors: it deletes out-of-range rows before any layer is drawn. coord_equal()
  # crops the view instead, so nothing is silently removed from the data. coord_equal
  # stays because the dashed y = x reference only reads as agreement at 45 degrees.
  # clip = "on": with clip = "off" the ~5% of points outside the frame rendered
  # beyond the panel border as loose dots floating in the margin. Clipping the panel
  # is the honest option -- the data is not censored (the scales are untouched), it is
  # simply not drawn outside its own frame. The region banner is unaffected: its
  # overflow is handled separately by strip_region_x(clip = "off").
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "on") +
  # "(NHD vs CON, occipital)" on the y axis is caption text -- Zhou being occipital-only
  # belongs in the legend, not inside the panel. Shortened to match the astro sibling.
  # Axis label follows CELLTYPE. It was hard-coded to "microglia" and silently
  # mislabelled the oligodendrocyte and OPC panels on their first render.
  labs(x = bquote("our" ~ .(CT_LAB) ~ "log"[2]*"FC (NHD vs CON)"),
       y = expression("Zhou 2023 log"[2]*"FC")) +
  theme_pub(base_size = 8 * TYPE_F) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.line = element_blank(),
        strip.background = element_blank(),
        panel.spacing.x = unit(0.4, "lines"),
        axis.title = element_text(size = 7.0),   # explicit: the x title must fit the 2.05-in canvas
        legend.position = "bottom", legend.text = element_text(size = 7.2 * TYPE_F),
        legend.key.size = unit(0.22, "cm"), legend.spacing.x = unit(0.10, "cm"),
        plot.margin = margin(4, 6, 2, 4))
# portrait, one column -- the same canvas as scripts 3 / 62 (2.05 x 3.99 in)
W <- 2.05; H <- 3.99
ragg::agg_png(file.path(PANEL, paste0(OUT_BN, ".png")), width = W, height = H,
              units = "in", res = 600); print(p); invisible(dev.off())
ggsave(file.path(PANEL, paste0(OUT_BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
logln("wrote %s.{png,pdf}  (cell type: %s)", OUT_BN, CELLTYPE)

cat("\n== sessionInfo ==\n", file = logcon)
capture.output(sessionInfo(), file = logcon, append = TRUE)
close(logcon)
cat("provenance:", LOGF, "\n")
cat("=== DONE ===\n", file = stderr())
