#!/usr/bin/env Rscript
# =============================================================================
# 99c_crosscohort_donor861_panel_FH.R — Supplementary Fig. 2, panel f: the cross-cohort contrast with and without the shared donor (861 = Zhou NHD3).
# -----------------------------------------------------------------------------
# Why: the Zhou comparison could be the same brain
# measured twice. Panel e shows the twelve comparisons with all three Zhou donors; this
# panel repeats the three headline cell types in three arms — all three donors, the two
# independent donors only, donor 861 only — so the reader sees directly which agreement
# is external and which is the shared donor.
#
# One ROW (the space under panel e is 133 pt), same grammar as panel e. Left: genome-wide Spearman rho over all shared genes
# (threshold-free, so it is not penalised by the loss of DESeq2 power at n = 2). Bottom:
# genes agreeing in direction among those significant in both cohorts, drawn as a bar
# from the empirical baseline to the observed value (the bar is the excess), point size =
# n, hollow below 20 genes, and an explicit "n = 0" where the two-donor arm yields no
# Zhou-significant gene at all (excitatory neurons).
#
# Arms are SHAPES so no new hue is introduced: circle = all three, diamond = without 861,
# triangle = 861 only; fill = region (as in e). Per-donor arms (NHD1 alone, NHD2 alone)
# ship in the source table, not on the panel.
#
# Reads tables/ST13_zhou_crosscohort_arms_FH.csv (106) only; no recomputation.
# Output: figures/Supplementary/SuppFig2_Zhou_QC/panels/SuppFig2_f_crosscohort_donor861.{pdf,png}
#         tables/suppfig2f_crosscohort_donor861_FH.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggh4x); library(ragg); library(patchwork) })
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Supplementary", "SuppFig2_Zhou_QC", "panels")
say <- function(...) cat(sprintf(...), "\n")
LOW_N <- 20L

SRC <- file.path(TDIR, "ST13_zhou_crosscohort_arms_FH.csv")
stopifnot("MISSING — run scripts/106_zhou_crosscohort_arms_FH.R" = file.exists(SRC))
d <- fread(SRC)
ARMS <- c(all3 = "All three donors (3 v 11)", no861 = "Without 861 (2 v 11)", only861 = "861 only (1 v 11)", cmp = "Without NHD1 or NHD2 (2 v 11)")
CMP  <- c("noNHD1", "noNHD2")   # equal-power comparators: drop one of the other donors (2 v 11) — small open diamonds
CTS  <- c("Micro-PVM" = "Micro-PVM", "Astro" = "Astrocytes", "Oligo" = "Oligodendrocytes", "OPC" = "OPC",
          "Neuron_Ex" = "Excitatory neurons", "Neuron_Inh" = "Inhibitory neurons")   # all six, as panel e
d <- d[arm %in% c(names(ARMS), CMP) & cell_type %in% names(CTS)]
stopifnot("expected 5 arms x 6 cell types x 2 regions" = nrow(d) == 60L)
stopifnot("unshrunken rho column missing — re-run 106" = "rho_full_mle" %in% names(d))
d[, `:=`(pct    = 100 * concord_frac,
         base   = 100 * baseline_concord,
         low_n  = n_both_sig < LOW_N,
         Arm    = factor(ifelse(arm %in% CMP, ARMS[["cmp"]], ARMS[arm]), levels = unname(ARMS)),
         cmp    = arm %in% CMP,
         solo   = n_NHD_zhou == 1L,
         Region = factor(REGION_FULL[region], levels = unname(REGION_FULL[REGION_ORDER])),
         ct     = factor(CTS[cell_type], levels = rev(unname(CTS))))]
d[, lab := ifelse(n_both_sig == 0, "n = 0", sprintf("n = %d", n_both_sig))]
# explicit vertical offsets (dodge shifted only one end of each connector, drawing fans):
# all three donors on top, without 861 in the middle, 861 only at the bottom
OFF <- c(0.24, 0, -0.24, 0); names(OFF) <- unname(ARMS)   # comparators share the middle row
d[, yy := as.numeric(ct) + OFF[as.character(Arm)]]
YB <- seq_along(levels(d$ct)); YL <- levels(d$ct)
fwrite(d[order(Region, ct, Arm)], file.path(TDIR, "suppfig2f_crosscohort_donor861_FH.csv"))
say("== panel source =="); print(as.data.frame(d[order(Region, ct, Arm), .(cell_type, region, arm, n_NHD_zhou, n_CON_zhou,
  rho_mle = round(rho_full_mle, 2), rho_apeglm = round(rho_full, 2), n_both_sig, concordance = round(pct), baseline = round(base))]), row.names = FALSE)
main <- d[cmp == FALSE]

SHP <- c(21, 23, 24, 23); names(SHP) <- unname(ARMS)
# Colour in addition to shape — arms carry the hue (region is already the facet);
# Okabe-Ito, no red-green pair; equal-power comparators stay small and grey. All text black.
ARMCOL <- c("#000000", "#D55E00", "#0072B2", "#7F7F7F"); names(ARMCOL) <- unname(ARMS)
common <- list(
  scale_y_continuous(breaks = YB, labels = YL, limits = c(0.5, length(YB) + 0.5), expand = expansion(0)),
  facet_wrap2(~ Region, nrow = 1, strip = strip_region_x(REGION_ORDER, fontsize = 6.6, clip = "off")),
  scale_fill_manual(values = ARMCOL, name = NULL, drop = FALSE),
  scale_colour_manual(values = ARMCOL, name = NULL, drop = FALSE),
  scale_shape_manual(values = SHP, name = NULL, drop = FALSE),
  labs(y = NULL),
  theme_pub(base_size = 8),
  theme(panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25), panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 7.2, colour = "black"), axis.text.x = element_text(size = 6.8, colour = "black"),
        axis.title.x = element_text(size = 7.2, colour = "black"), strip.text = element_text(size = 7.4, face = "plain", colour = "black"),
        legend.key.size = unit(9, "pt"), legend.key.spacing.y = unit(2.5, "pt"),   # slot taller than the 2.1-size glyphs
        legend.text = element_text(size = 6.8, colour = "black", margin = margin(l = 3)),
        legend.title = element_text(size = 6.8, colour = "black", lineheight = 0.95, margin = margin(b = 4)),
        panel.spacing.x = unit(3, "pt"), strip.clip = "off",
        plot.margin = margin(4, 4, 4, 4)))

# top: genome-wide rho — the power-insensitive quantity
# left: rho on unshrunken DESeq2 log2FC.
# 1 v 11 arm hollow (descriptive); equal-power comparators (drop NHD1 / drop NHD2) as small open diamonds
p_rho <- ggplot(main, aes(y = yy)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = rho_full_mle, yend = yy), colour = "grey72", linewidth = 0.8, lineend = "butt") +
  geom_point(data = main[solo == FALSE], aes(x = rho_full_mle, fill = Arm, shape = Arm), size = 1.7, colour = "black", stroke = 0.3) +
  geom_point(data = main[solo == TRUE],  aes(x = rho_full_mle, shape = Arm, colour = Arm), fill = "white", size = 1.7, stroke = 0.55) +
  geom_point(data = d[cmp == TRUE], aes(x = rho_full_mle, y = yy, shape = Arm, colour = Arm), fill = "white", size = 1.0, stroke = 0.4) +
  scale_x_continuous(expression("Spearman " * rho * " (all shared genes)"), limits = c(-0.2, 0.85),   # unshrunken log2FC: stated in Methods, not on the axis
                     breaks = c(0, 0.4, 0.8)) +
  common + guides(shape = guide_legend(override.aes = list(size = c(2.1, 2.1, 2.1, 1.3), fill = c(ARMCOL[1], ARMCOL[2], "white", "white"),
                                                           colour = c("black", "black", ARMCOL[3], ARMCOL[4]), stroke = c(0.3, 0.3, 0.55, 0.4)), order = 1, ncol = 1),
                  fill = "none", colour = "none")

# bottom: direction agreement among both-significant genes, bar from the empirical baseline
p_conc <- ggplot(main, aes(y = yy)) +
  geom_segment(data = main[n_both_sig > 0], aes(x = base, xend = pct, yend = yy),
               colour = "grey72", linewidth = 0.8, lineend = "butt") +
  geom_point(aes(x = base), shape = 124, size = 1.8, colour = "black") +   # tick = genome-wide sign agreement (readable when n = 0)
  geom_point(data = main[n_both_sig > 0 & solo == FALSE], aes(x = pct, fill = Arm, shape = Arm, size = n_both_sig),
             colour = "black", stroke = 0.3) +      # no low-n fade: n is printed on every row, and alpha turned vermilion into a fourth hue
  geom_point(data = main[n_both_sig > 0 & solo == TRUE], aes(x = pct, shape = Arm, colour = Arm, size = n_both_sig),
             fill = "white", stroke = 0.55) +
  geom_text(aes(x = 108, label = lab), hjust = 0, size = 2.2, colour = "black") +
  scale_x_continuous("Direction agreement (%)", limits = c(42, 136), breaks = c(50, 75, 100),
                     expand = expansion(mult = c(0.02, 0))) +
  scale_size_continuous(range = c(0.9, 2.5), breaks = c(10, 50, 400), name = "Genes significant\nin both cohorts") +
  common + guides(shape = "none", fill = "none", colour = "none", size = guide_legend(order = 2, ncol = 1, override.aes = list(shape = 21, fill = "grey55", colour = "black")))

# This panel replaces the old e (99b summary) and takes the full row — the two
# blocks sit side by side (rho left, agreement right; the right block drops its duplicate y labels),
# legends collected at the bottom in one horizontal row.
p_conc <- p_conc + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p <- p_rho | p_conc
# a short plain title top-left so the reader knows what the panel tests
p <- p + plot_layout(guides = "collect", widths = c(1, 1.45)) +
  plot_annotation(title = "Replication in the Zhou et al. cohort by NHD donor set",
                  theme = theme(plot.title = element_text(size = 7.8, colour = "black", face = "plain", hjust = 0, margin = margin(0, 0, 3, 0)))) &
  theme(legend.position = "right", legend.box = "vertical", legend.spacing.y = unit(14, "pt"), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(8, 0, 0, -2), legend.justification = "top")
BN <- "SuppFig2_e_crosscohort_donor861"; W <- 7.05; H <- 2.85   # legends on the right; H 2.2 -> 2.7 so the sub-row pitch exceeds the n-label font (page has the room)
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600, device = ragg::agg_png)
say("Saved %s.{pdf,png} (%.2f x %.2f in)", BN, W, H)
writeLines(capture.output(sessionInfo()), file.path(PROJ, "logs", "99c_crosscohort_donor861_panel_FH_sessionInfo.txt"))
cat("\n=== DONE ===\n", file = stderr())
