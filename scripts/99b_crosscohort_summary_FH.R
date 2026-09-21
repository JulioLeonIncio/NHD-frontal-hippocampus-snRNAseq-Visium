#!/usr/bin/env Rscript
# =============================================================================
# 99b_crosscohort_summary_FH.R — Every cross-cohort comparison in one panel.
# Supplementary Fig. 2, panel e.
# -----------------------------------------------------------------------------
# Why a summary rather than more SCATTERPLOTS. The paper shows four cross-cohort
# scatter panels (astrocytes in Fig 3c, neurons in Fig 5c, microglia in Supplementary
# Fig 3b) and never shows the oligodendrocyte lineage at all. Adding two more scatters
# would have put the WEAKEST comparison in the set into a main figure: oligodendrocyte
# frontal agrees in 65% of 31 genes against a 48% baseline, an excess of 17 points at
# binomial p = 0.046, and the genes a reader looks for there -- ANLN, MBP, PLP1, CNP --
# are all non-significant in Zhou's data, so they would render as grey "n.s." dots.
# This panel instead shows all twelve comparisons at once, honestly ranked.
#
# The three THINGS A bare percentage HIDES, and how this panel shows them:
#   1. the baseline VARIES, from 48% to 77%. Concordance must be read against its own
#      empirical baseline, never against 50% [[nhd-crosscohort-validation-methodology]].
#      The bar is drawn from the baseline to the observed value, so the length of the bar
#      is the excess -- the quantity that carries the claim -- rather than the percentage.
#   2. power VARIES enormously, from 5 genes to 75. Inhibitory neurons reach 100% on five
#      genes; astrocytes reach 97% on sixty-four. Those are not the same statement, so n
#      is printed against every row and rows below 20 genes are drawn hollow.
#   3. direction and magnitude can disagree. Hippocampal microglia are 82% concordant
#      while their genome-wide rho is NEGATIVE (-0.12). rho was tried on the panel and cut:
#      at two facets wide it overran the edge and collided with points near 100%. It ships
#      in ST13 and in this panel's source CSV, and the legend states the microglial case.
#
# Reads ST13 only; no recomputation. Output:
#   figures/Supplementary/SuppFig2_Zhou_QC/panels/SuppFig2_e_crosscohort_summary.{pdf,png}
#   tables/suppfig2e_crosscohort_summary_FH.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggh4x); library(ragg)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Supplementary", "SuppFig2_Zhou_QC", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

LOW_N <- 20L   # below this many both-significant genes the point is drawn hollow

ST13 <- file.path(TDIR, "ST13_zhou_crosscohort_FH.csv")
stopifnot("MISSING ST13 — run scripts/70_zhou_crosscohort_FH.R" = file.exists(ST13))
d <- fread(ST13)
need <- c("cell_type","region","n_both_sig","concord_frac","baseline_concord","rho_full","binom_p","n_shared")
stopifnot("ST13 is missing a required column" = all(need %in% names(d)))
d <- d[, ..need]
stopifnot("expected 12 comparisons (6 cell types x 2 regions)" = nrow(d) == 12L)

# canonical cell-type order, matching the atlas ordering used elsewhere
CT_ORDER <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
CT_LAB   <- c("Micro-PVM" = "Micro-PVM", "Astro" = "Astrocytes", "Oligo" = "Oligodendrocytes",
              "OPC" = "OPC", "Neuron_Ex" = "Excitatory neurons", "Neuron_Inh" = "Inhibitory neurons")
stopifnot("ST13 carries an unexpected cell type" = all(d$cell_type %in% CT_ORDER))

d[, `:=`(pct       = 100 * concord_frac,
         base      = 100 * baseline_concord,
         excess    = 100 * (concord_frac - baseline_concord),
         low_n     = n_both_sig < LOW_N,
         Region    = factor(REGION_FULL[region], levels = unname(REGION_FULL[REGION_ORDER])),
         ct        = factor(CT_LAB[cell_type], levels = rev(unname(CT_LAB[CT_ORDER]))))]
# The label carries n, the quantity the percentage most badly hides. rho was tried here too
# and was dropped: at two facets wide it ran past the panel edge and collided with points
# near 100%. rho ships in ST13 and the legend calls out the one case where it matters --
# hippocampal microglia, 82% concordant with rho negative at -0.12.
# Carry the information the scatter panels show — rho over all shared genes and
# that gene count — into the row label, not only the both-significant n.
d[, lab := sprintf("n = %d; \u03c1 = %.2f\nover %s genes", n_both_sig, rho_full, trimws(format(n_shared, big.mark = ",")))]
setorder(d, Region, ct)
fwrite(d, file.path(TDIR, "suppfig2e_crosscohort_summary_FH.csv"))

say("\n== cross-cohort agreement, all twelve comparisons ==")
print(as.data.frame(d[order(-excess), .(cell_type, region, n_both_sig,
      concordance = round(pct), baseline = round(base), excess = round(excess),
      rho = round(rho_full, 2), low_power = low_n)]), row.names = FALSE)
say("\n  rows below %d both-significant genes (drawn hollow): %d of %d",
    LOW_N, sum(d$low_n), nrow(d))
say("  baselines span %.0f%% to %.0f%% — which is why the bar starts at the baseline",
    min(d$base), max(d$base))

p <- ggplot(d, aes(y = ct)) +
  # The bar is the excess over baseline: it starts at the empirical baseline, not at zero
  geom_segment(aes(x = base, xend = pct, yend = ct), colour = "grey72", linewidth = 1.05,
               lineend = "round") +
  geom_point(aes(x = base), shape = 124, size = 2.0, colour = "grey45") +
  geom_point(aes(x = pct, fill = Region, size = n_both_sig, alpha = !low_n),
             shape = 21, colour = "grey25", stroke = 0.3) +
  geom_text(aes(x = 105, label = lab), hjust = 0, size = 1.55, colour = "grey30", lineheight = 0.9) +
  facet_wrap2(~ Region, nrow = 1,
              strip = strip_region_x(REGION_ORDER, fontsize = 7.2, clip = "off")) +
  scale_x_continuous("Genes agreeing in direction (%)", limits = c(40, 158),
                     breaks = c(50, 75, 100), expand = expansion(mult = c(0.02, 0))) +
  scale_size_continuous(range = c(1.2, 4.0), breaks = c(10, 40, 75),
                        name = "genes significant\nin both cohorts") +
  scale_fill_manual(values = unname(PAL_REGION[REGION_ORDER]), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(y = NULL) +
  theme_pub(base_size = 8) +
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
        axis.text.y = element_text(size = 6.6, colour = "black"),
        axis.text.x = element_text(size = 6.2),
        axis.title.x = element_text(size = 6.8),
        strip.text = element_text(size = 7.2, face = "plain"),
        legend.position = "bottom", legend.direction = "horizontal", legend.key.size = unit(0.42, "lines"),
        legend.text = element_text(size = 5.8), legend.title = element_text(size = 6.0),
        panel.spacing.x = unit(3, "pt"), strip.clip = "off",
        plot.margin = margin(4, 4, 4, 4))

BN <- "SuppFig2_e_crosscohort_summary"; W <- 3.30; H <- 3.85   # e and f side by side, legends at the bottom
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
say("\nSaved %s.{pdf,png} (%.2f x %.2f in)", BN, W, H)
say("Wrote tables/suppfig2e_crosscohort_summary_FH.csv")
writeLines(capture.output(sessionInfo()),
           file.path(PROJ, "logs", "99b_crosscohort_summary_FH_sessionInfo.txt"))
say("=== DONE ===")
