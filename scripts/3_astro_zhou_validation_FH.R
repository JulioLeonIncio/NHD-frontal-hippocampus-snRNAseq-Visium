#!/usr/bin/env Rscript
# =============================================================================
# 3_astro_zhou_validation_FH.R — F3c_astro_zhou_validation_FH.R — Figure 3: cross-cohort validation vs Zhou 2023, ASTROCYTES, frontal + hippocampus. ASTRO sibling of 2H_zhou_validation_FH.R.
# -----------------------------------------------------------------------------
# Our OCC is gone in this rebuild, so this is cross-region GENERALIZATION (our
# Frontal / Hippocampus vs Zhou occipital) — an extra region difference on TOP of
# the cohort difference — not same-region replication.  Two facets: Frontal |
# Hippocampus, each vs Zhou.  Generalization is strong for astrocytes.
#
#   x = our avg_log2FC (NHD vs CON) in the facet's region
#   y = Zhou log2FoldChange (NHD vs CON, occipital)
#   plotted = MAST discovery set (is_disc) intersect Zhou-measured; grey = n.s. in
#             Zhou, black = Zhou-sig & concordant, orange = Zhou-sig & discordant.  A
#             discovery gene is never dropped for failing Zhou sig.
#   headline (from tables/ST13_zhou_crosscohort_FH.csv, Astro rows
#             source, never hand-typed): Spearman rho over all shared genes,
#             direction-consistent % over discovery vs the region-matched
#             genome-wide baseline (binom p vs baseline, not 0.5).
#   framing: direction-consistent / generalization, never "replicated".
#
# Astro-specific label deviation (labels only, not stats): force-label the
# metallothionein programme + GFAP anchors, restricted to genes present +
# CONCORDANT (never mislabel a discordant/absent anchor).
#
# Atlas-free: reads tables/fig3f_astro_zhou_FH_{Frontal,Hippo}_genes.csv +
# ST13_zhou_crosscohort_FH.csv (built by scripts/70_zhou_crosscohort_FH.R).
# Output: F3c_astro_zhou_validation.{png,pdf} + provenance log.
# -----------------------------------------------------------------------------
# Harmonised to Figure 2's house type spec (base 8 / axis 6.3 /
# title 7 / legend 6.2) and rendered at its placed size so no assembly down-scaling is
# needed. Fig 2's panels total ~117 in2 on a 67.9 in2 page, i.e. they are shrunk ~0.77x
# at assembly, which is what pushed its type under the 5 pt floor. Fig 3 avoids that.
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
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))

# Font scale. Shrinking this panel while type stayed at absolute house sizes made it read
# text-heavy beside its Fig-3 siblings (the x-axis title alone spanned nearly the full
# width). Same rule as Fig 2's panel 2I: hold the house ratio but clamp so the smallest
# token lands on the 5 pt print floor. Strict ratio (2.92/3.90 = 0.75) would put axis text
# at 4.7 pt, under the floor; 5.0/6.3 = 0.79 is the largest factor that respects it.
FS <- 5.0 / 6.3
sz <- function(x) round(x * FS, 2)   # theme_pub, REGION_FULL, strip_region_x, REGION_ORDER
source(file.path(PROJ, "scripts", "_artifact_genes.R"))           # is_artifact

CT      <- "Astro"
REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
PANEL <- file.path(PROJ, "figures/Figure_3/panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ZHOU_CLS_COL <- PAL_ZHOU_CLS   # theme token (black / orange / grey); was hand-typed green / red

LOGF <- file.path(LOGS, "F3c_astro_zhou_validation_FH_provenance.txt")
logcon <- file(LOGF, open = "wt")
logln <- function(...) { s <- sprintf(...); cat(s, "\n"); cat(s, "\n", file = logcon) }
logln("== F3c_astro_zhou_validation_FH.R  run %s ==", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logln("PROJ = %s", PROJ)

gene_files <- setNames(
  file.path(TDIR, sprintf("fig3f_astro_zhou_FH_%s_genes.csv", REGIONS)), REGIONS)
st13_file  <- file.path(TDIR, "ST13_zhou_crosscohort_FH.csv")
for (f in c(gene_files, st13_file))
  if (!file.exists(f)) stop(sprintf("MISSING %s — run scripts/70_zhou_crosscohort_FH.R first", f))
for (f in c(gene_files, st13_file))
  logln("input %-44s  mtime %s", basename(f), format(file.mtime(f), "%Y-%m-%d %H:%M:%S"))

st13 <- read.csv(st13_file, stringsAsFactors = FALSE)
st13 <- st13[grepl("Astro", st13$cell_type), ]
stat_for <- function(rg) { r <- st13[st13$region == rg, ]; stopifnot(nrow(r) == 1); r }

# gene CSV schema: gene, our_lfc, our_delta, zhou_lfc, zhou_padj, is_disc, zhou_sig,
# both_sig, artifact, concordant  (note our_lfc name matches 2H's our_lfc).
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

# Lim was max(|log2FC|), so a single extreme gene stretched both axes to ~+-12
# while the mass of the data sits near +-5 -- roughly 60% of each facet was empty (same defect
# as the Fig-2 Zhou panel). Use the 99th percentile instead; points beyond it are squished to
# the frame by oob rather than dictating it. annot_df keys to `lim`, so the corner annotation
# follows automatically.
# Tightened further, 0.99 -> 0.975 quantile and a smaller pad. The concordant
# mass sits well inside +-5, so the previous +-7.1 frame left the top-left and bottom-right
# CORNERS empty (the cloud is diagonal). Points outside are squished to the frame by oob.
lim <- as.numeric(stats::quantile(abs(c(m$our_lfc, m$zhou_lfc)), 0.975, na.rm = TRUE)) * 1.02
cat(sprintf("Zhou axis limit: +-%.1f (99th pct); max |log2FC| was %.1f\n",
            lim, max(abs(c(m$our_lfc, m$zhou_lfc)), na.rm = TRUE)))
# Report how much of the cloud each candidate frame would actually contain, so the limit is
# chosen from the data rather than by eye. oob=squish pins anything outside onto the frame,
# so a limit that excludes much data shows up as a pile-up line on the border.
for (.q in c(0.90, 0.95, 0.975, 0.99)) {
  .l <- as.numeric(stats::quantile(abs(c(m$our_lfc, m$zhou_lfc)), .q, na.rm = TRUE)) * 1.02
  cat(sprintf("   frame +-%.2f (q%.3f) would contain %.2f%% of points\n", .l, .q,
              100 * mean(abs(m$our_lfc) <= .l & abs(m$zhou_lfc) <= .l, na.rm = TRUE)))
}

# Labels: strongest Zhou-sig genes by joint magnitude (artifact excluded) union
# curated astro anchors, anchors restricted to present + concordant.
ANCHORS_REQ <- c("MT2A","MT3","MT1G","GFAP","CD44","SERPINA3")
anchor_ok <- m %>% filter(gene %in% ANCHORS_REQ, concordant) %>% pull(gene) %>% unique()
logln("requested anchors: %s | labelled (present+concordant): %s | dropped: %s",
      paste(ANCHORS_REQ, collapse = ","),
      if (length(anchor_ok)) paste(anchor_ok, collapse = ",") else "none",
      paste(setdiff(ANCHORS_REQ, anchor_ok), collapse = ","))
lab_df <- m %>% filter(zhou_sig, !artifact) %>% group_by(region) %>%
  # Labels per facet. At 6.60 x 2.10 (two facets) ten italic gene labels
  # collided (UPP2 over the n= annotation, NAMPT/MYO5C, SERPINA3/GPSB1) and fought the corner
  # statistics for the same space.
  slice_max(abs(our_lfc) + abs(zhou_lfc), n = 6, with_ties = FALSE) %>% ungroup() %>%
  bind_rows(m %>% filter(gene %in% anchor_ok)) %>% distinct(region, gene, .keep_all = TRUE)

fmt_p <- function(p) { if (is.na(p)) return("NA")
  if (p <= .Machine$double.xmin) return(sprintf("< %.0e", .Machine$double.xmin))
  if (p < 1e-3) sprintf("%.0e", p) else sprintf("%.2g", p) }
# astro generalization is STRONG (rho + concord both high) -> lead with rho, then
# the direction statistic (both agree, no contradiction to disclaim).
# round-4 visual fix: the on-panel second line was a
# dense sentence ("direction-consistent NN% vs NN% baseline (p ...), n = ...") that
# drifted toward caption-in-panel and was wordier than the terse rho/delta annotations
# on the sibling volcano/violins. Trimmed to the effect-size ESSENTIALS on-panel
# ("NN% concordant, n = ..."); the genome-wide baseline % and the binomial p-value are
# moved to the figure LEGEND (still logged
# below for the legend author). not a data change — same numbers, fewer on-plot.
# round-5 visual fix: quote the number inside the
# plotmath expression so parse=TRUE keeps the trailing zero (0.30, not 0.3). Without
# the quotes, plotmath parses 0.30 as a numeric literal and strips the zero, so the two
# facets showed inconsistent precision (Frontal 0.37 vs Hippocampus 0.3). Number itself
# unchanged
annot_df$l1 <- sprintf("rho == '%.2f'", annot_df$rho)
annot_df$l2 <- sprintf("%.0f%% of %d both-sig (baseline %.0f%%)\nrho over n = %s",   # shortened: overran the frame
                       100*annot_df$concord, annot_df$n_both,
                       100*annot_df$baseline,
                       format(annot_df$n_shared, big.mark = ","))

p <- ggplot(m, aes(our_lfc, zhou_lfc)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_point(data = filter(m, cls == "n.s. in Zhou"), aes(colour = cls),
             size = 0.7, alpha = 0.45, stroke = 0) +
  geom_point(data = filter(m, cls != "n.s. in Zhou"), aes(colour = cls, shape = cls),
             size = 1.6, alpha = 0.9, stroke = 0.55) +
  geom_text_repel(data = lab_df, aes(label = gene), size = sz(2.4), fontface = "italic",
                  colour = "grey15", bg.color = "white", bg.r = 0.12, point.padding = 0.5, max.overlaps = Inf, force = 3, seed = 42,
                  box.padding = 0.28, min.segment.length = 0, segment.size = 0.2,
                  segment.colour = "grey60") +
  # Keep both annotation lines inside the frame -- nudging l1 to lim*1.02 pushed
  # the rho value off the panel entirely.
  geom_text(data = annot_df, aes(x = -lim, y = lim * 0.96, label = l1),
            hjust = 0, vjust = 1, size = 2.1, colour = "grey20", parse = TRUE, inherit.aes = FALSE) +
  geom_text(data = annot_df, aes(x = -lim, y = lim * 0.82, label = l2),
            hjust = 0, vjust = 1, size = 1.9, colour = "grey20", inherit.aes = FALSE) +
  # Vertical -- Frontal above Hippocampus. REGIONS is already in canonical
  # Frontal -> Hippocampus order, and facet_wrap fills by row, so ncol = 1 stacks them correctly.
  ggh4x::facet_wrap2(~ region, ncol = 1, strip = strip_region_x(REGIONS, clip = "off")) +
  scale_colour_manual(values = ZHOU_CLS_COL, name = NULL, drop = FALSE,
                      guide = guide_legend(nrow = 2, byrow = TRUE,
                                           override.aes = list(size = 2.2, alpha = 1, stroke = 0.7))) +
  scale_shape_manual(values = SHAPE_ZHOU_CLS, name = NULL, drop = FALSE,
                     guide = guide_legend(nrow = 2, byrow = TRUE,
                                          override.aes = list(size = 2.2, alpha = 1, stroke = 0.7))) +
  scale_x_continuous(limits = c(-lim, lim)) + scale_y_continuous(limits = c(-lim, lim)) +
  # Axes were running +-12 for data living in roughly -5..+5, leaving ~60% of
  # both facets empty (same defect as the Fig-2 Zhou panel). Crop to the 1st/99th percentile
  # of the plotted values, symmetric, with a small pad; coord_equal is kept so the y=x
  # reference line still reads at 45 degrees.
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
  labs(x = expression("our astrocyte log"[2]*"FC (NHD vs CON)"),
       # Was "...FC (NHD vs CON, occipital)" -- too long for a 2.10 in panel
       # (it clipped), and the parenthetical is caption text that belongs in the legend.
       y = expression("Zhou 2023 log"[2]*"FC")) +
  theme_pub(base_size = sz(8)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.line = element_blank(), strip.background = element_blank(),
        panel.spacing.x = unit(0.4, "lines"),
        # At 2.90 in wide the 3-item bottom legend ran off ("n.s. in Zh...").
        # Smaller text + tighter keys, and wrapped to 2 rows via guide_legend(nrow = 2) above.
        legend.position = "bottom", legend.text = element_text(size = sz(6.2)),
        legend.key.size = unit(0.22, "cm"), legend.spacing.x = unit(0.10, "cm"),
        plot.margin = margin(4, 6, 2, 4))
# The cause is
# coord_equal: with xlim = ylim = +-5 it forces square facets, so in a 6.60 x 2.10 box the
# two squares could only be ~2.1 in wide and ~2.4 in of white was left over, split down the
# left and right margins. coord_equal must stay -- the dashed y = x reference line is only
# interpretable at 45 degrees -- so the PANEL is resized to the aspect it actually wants:
# plot height ~1.8 in => two 1.8 in squares + ~1.0 in of axis/label furniture ~= 4.60 in.
# +20% wider (2.90 -> 3.48). height scaled by the same factor (4.60 -> 5.52)
# on purpose: coord_equal forces square facets, so their size is set by whichever dimension
# is binding. Widening alone would leave the squares at their height-limited size and simply
# re-open the left/right white margins we just removed. Both dimensions must move together.
# 30% smaller -> 2.44 x 3.86. both dimensions scale by the same 0.7 so the
# coord_equal aspect is preserved and no left/right white margin re-opens.
# 30% smaller, then 20% wider -- both dimensions scale together (x0.7 then
# x1.2 = x0.84 overall) because coord_equal forces square facets; moving width alone would
# re-open the left/right white margins.
# At 2.92 x 4.64 the rendered panel carried 0.36 in of white on
# the left and 0.38 in on the right -- 25% of its width. Cause: coord_equal makes the facets
# square, they were height-limited, and ggplot centres the plot area, so all the horizontal
# slack became equal side margins. Keeping the requested width and growing the height until
# the squares become width-limited removes it. Verified by measuring the ink bounding box of
# the rendered PNG, not by eye.
# a further 30% smaller -> 2.04 x 3.99. both dimensions scale (coord_equal), so
# the side margins shrink proportionally too and no new white opens up.
W <- 3.48 * 0.7 * 1.2 * 0.7; H <- 5.70 * 0.7
# Height chosen by measuring the rendered ink bbox, not by eye. Side margins vs height at
# this width:  4.64 -> 0.36/0.38 in |  5.30 -> 0.19/0.21 |  5.70 -> 0.09/0.12 |  6.00 -> 0.06/0.08.
# 5.70 is the knee: the remaining ~0.1 in is essentially the plot margin, while 6.00 would make
# the panel two thirds of the page height for a further 0.03 in. If the residual sliver still
# bothers, the alternative is narrowing to ~2.70 rather than growing taller.

ragg::agg_png(file.path(PANEL, "F3c_astro_zhou_validation.png"), width = W, height = H,
              units = "in", res = 600); print(p); invisible(dev.off())
ggsave(file.path(PANEL, "F3c_astro_zhou_validation.pdf"), p, width = W, height = H, useDingbats = FALSE)
logln("wrote F3c_astro_zhou_validation.{png,pdf}")

cat("\n== sessionInfo ==\n", file = logcon)
capture.output(sessionInfo(), file = logcon, append = TRUE)
close(logcon)
cat("provenance:", LOGF, "\n")
cat("=== DONE ===\n", file = stderr())
