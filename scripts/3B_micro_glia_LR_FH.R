#!/usr/bin/env Rscript
# =============================================================================
# 3B_micro_glia_LR_FH.R — Figure 3 panel b: microglia -> glial-receiver ligand-receptor rewiring, frontal + HIPPOCAMPUS (2 regions). Faithful port of
# 09a_v2_micro_to_glia_LR.R to the FH rebuild.
# -----------------------------------------------------------------------------
# Pulls the curated CellChat Micro-PVM -> {Astro, OPC, Oligo} l-r pairs per region
# (Frontal + Hippocampus), computes NHD - CON dprob, renders as a pathway-grouped
# dot plot.  Reads the aggregated subsetCommunication output directly
# (tables/cellchat_all_LR_pairs.csv) — atlas-light, no merged-object load.
#
# Critical gotcha: filter by ligand / interaction_name, not pathway_name (GAS6's
# pathway is "GAS", PROS1's is "PROS").
#
# FH re-validation of the curated pair whitelist
# against tables/cellchat_all_LR_pairs.csv.  A pair is shown only if it is present
# in the FH CellChat output for >=1 region/condition; each shown pair's per-region
# CON/NHD prob + p + BH-q + dprob is written to fig3b_LR_stats.csv so exactly what
# is significant per receiver is auditable.  Result (see log):
#   SPP1 -> CD44(astro) + ITGAV/ITGA8/ITGA9-integrins(OPC,oligo)  present CON+NHD.
#   GAS6 -> MERTK(astro) / TYRO3(OPC,oligo) / AXL(astro,low)      present CON+NHD;
#           GAS6->astro-MERTK up in NHD (F +0.026, H +0.005 dprob).
#   PROS1-> MERTK / TYRO3 / AXL   NHD-ONLY (2nd TAM ligand, added).
#   PSAP -> GPR37(oligo) / GPR37L1(astro,OPC)   nhd-only.
#   SEMA4D-> PLXNB1   present.
# Curated pairs absent from the FH output are dropped + reported (never shown as a
# blank row).  No pathway-level p-value cherry-picking; direction shown honestly.
#
# Full region names on the top banner; no in-panel title/caption (moved to
# legend); pale region strips.  Single-donor Hippocampus caveat -> manuscript text.
#
# Output: figures/Figure_3/panels/3b_micro_glia_LR.{pdf,png} + fig3b_LR_stats.csv
# -----------------------------------------------------------------------------
# Harmonised to figure 2. Fig 3 had a 1.5x internal spread in type
# (axis.text 5.6 in the LR/receptor dotplots vs 8.5 in the dotmatrix/gas6 panels) and none of
# it matched Fig 2's house spec. All Fig-3 panels now use the same tokens as Fig 2's
# 2D_secretome_MAST reference -- base_size 8 / axis.text 6.3 / axis.title 7 / strip 6.4-7 /
# legend 6.2 -- rendered at their placed width so 6.3 pt is 6.3 pt on the page.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(patchwork)
  library(dplyr); library(purrr); library(tidyr); library(ggplot2); library(scales)
  library(forcats); library(ggh4x)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, REGION_FULL, strip_region_x

LR    <- file.path(PROJ, "tables", "cellchat_all_LR_pairs.csv")
PANEL <- file.path(PROJ, "figures",
                   if (tolower(Sys.getenv("RECV_SET","astro")) == "oligo") "Figure_4" else "Figure_3",
                   "panels")
TDIR  <- file.path(PROJ, "tables")
# TDIR is defined before
# RECV_SET exists, so the suffix reads the env var directly, exactly as PANEL does above.
.stats_suffix <- if (tolower(Sys.getenv("RECV_SET", "astro")) == "oligo") "_oligo" else "_astro"
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(LR)) stop(sprintf("MISSING %s — CellChat prep step failed", LR))

# ---------------------------------------------------------------------------
# Receiver SET: "this figure only includes astrocytes, so why
# don't we include only astrocyte and use the micro OPC and Oligo for figure 4?"
# Correct -- Fig 3 is the astrocyte figure, yet this panel spent two thirds of its
# width on OPC and Oligo receivers. The panel is now parameterised rather than forked
# (a forked copy would drift from this one):
#   RECV_SET=astro (default) -> receivers = Astro          -> Figure_3/panels
#   RECV_SET=oligo           -> receivers = OPC + Oligo    -> Figure_4/panels
# This also gives Fig 4 the same "step out of Fig 2" opening that Fig 3 has: the
# microglial ligands of Fig 2's secretome arriving on the oligodendrocyte lineage.
RECV_SET <- tolower(Sys.getenv("RECV_SET", "astro"))
stopifnot("RECV_SET must be 'astro' or 'oligo'" = RECV_SET %in% c("astro","oligo"))
RECEIVERS <- if (RECV_SET == "astro") "Astro" else c("OPC","Oligo")
cat(sprintf("RECV_SET=%s -> receivers: %s\n", RECV_SET, paste(RECEIVERS, collapse=", ")))
REGIONS   <- c("Frontal","Hippo")
FDR_ALPHA <- 0.05

# Curated micro->glia programme pairs (FH-revalidated; PROS1->TAM added).
# interaction_name is the CellChat pair key.
CURATED_LR <- c(
  "SPP1_CD44", "SPP1_ITGAV_ITGB5", "SPP1_ITGAV_ITGB1",
  "SPP1_ITGA8_ITGB1", "SPP1_ITGA9_ITGB1",
  "GAS6_MERTK", "GAS6_TYRO3", "GAS6_AXL",
  # PROS1_AXL dropped — PROS1 (Protein S) does not activate
  # AXL; only GAS6 is an AXL agonist. PROS1 signals through MERTK/TYRO3 (biochem).
  "PROS1_MERTK", "PROS1_TYRO3",
  "PSAP_GPR37", "PSAP_GPR37L1",
  "SEMA4D_PLXNB1")
PATHWAY_ORDER <- c("SPP1","GAS","PROS","PSAP","SEMA4")   # ligand pathway groups (CellChat pathway_name values)
# The CellChat pathway_name is the ABBREVIATED ligand
# family ("GAS","PROS","SEMA4"). Relabel the row strips to the full ligand gene so a
# reader isn't left guessing (GAS -> GAS6, PROS -> PROS1, SEMA4 -> SEMA4D). Only the
# display changes; pathway_name factor levels stay the raw CellChat values.
PATHWAY_LABELS <- c(SPP1 = "SPP1", GAS = "GAS6", PROS = "PROS1",
                    PSAP = "PSAP", SEMA4 = "SEMA4D")

# ---- Load + per-(Region x Condition) BH-FDR --------------------------------
raw <- read.csv(LR, stringsAsFactors = FALSE) %>%
  mutate(Condition = sub("_.*$", "", region_cond),
         Region    = sub("^[^_]*_", "", region_cond)) %>%
  filter(source == "Micro-PVM", target %in% RECEIVERS,
         Region %in% REGIONS, interaction_name %in% CURATED_LR) %>%
  select(Region, Condition, source, target, ligand, receptor,
         interaction_name, interaction_name_2, pathway_name, prob, pval) %>%
  as_tibble()
stopifnot(nrow(raw) > 0)

# BH-FDR over each condition's per-region CellChat p-value pool (as reference).
# CellChat p-values here are per-region/condition over the full LR universe;
# we recompute BH on the curated subset's pool for a conservative "sig_any" flag
# consistent with the reference's disclosure column.
raw <- raw %>% group_by(Region, Condition) %>%
  mutate(q_BH = p.adjust(pval, method = "BH")) %>% ungroup()

cat(sprintf("Raw curated micro->glia rows: %d | unique pairs: %d\n",
            nrow(raw), length(unique(raw$interaction_name))))

# one row per (Region, target, LR, Condition)
agg <- raw %>%
  group_by(Region, target, interaction_name, interaction_name_2,
           ligand, receptor, pathway_name, Condition) %>%
  summarise(prob = max(prob, na.rm = TRUE), pval = min(pval, na.rm = TRUE),
            q_BH = min(q_BH, na.rm = TRUE), .groups = "drop")

wide <- agg %>%
  pivot_wider(id_cols = c(Region, target, interaction_name, interaction_name_2,
                          ligand, receptor, pathway_name),
              names_from = Condition, values_from = c(prob, pval, q_BH),
              values_fill = list(prob = 0, pval = 1, q_BH = 1)) %>%
  mutate(delta_prob = prob_NHD - prob_CON,
         max_prob   = pmax(prob_NHD, prob_CON),
         sig_any    = (q_BH_CON < FDR_ALPHA) | (q_BH_NHD < FDR_ALPHA))

n_tested <- nrow(wide); n_sig <- sum(wide$sig_any, na.rm = TRUE)
cat(sprintf("Condition-paired curated LR rows: %d | pass BH-FDR q<%g either cond: %d (%.1f%%)\n",
            n_tested, FDR_ALPHA, n_sig, 100 * n_sig / n_tested))

missing_lr <- setdiff(CURATED_LR, unique(wide$interaction_name))
if (length(missing_lr))
  cat("Curated LRs ABSENT from FH CellChat (dropped):",
      paste(missing_lr, collapse = ", "), "\n")

# per-pair mean dprob (for ordering) + pathway grouping
lr_rank <- wide %>%
  group_by(interaction_name, interaction_name_2, pathway_name) %>%
  summarise(mean_delta = mean(delta_prob, na.rm = TRUE), .groups = "drop") %>%
  mutate(pathway_name = factor(pathway_name, levels = PATHWAY_ORDER)) %>%
  arrange(pathway_name, desc(mean_delta))

cat(sprintf("\nFeatured %d curated LR pairs (%d NHD-up, %d CON-up mean dprob):\n",
            nrow(lr_rank), sum(lr_rank$mean_delta > 0), sum(lr_rank$mean_delta < 0)))
print(as.data.frame(lr_rank %>% mutate(mean_delta = round(mean_delta, 4))), row.names = FALSE)

# report SPP1 / GAS6 / PSAP status per receiver (the requested audit)
cat("\n== per-receiver programme status (max prob CON/NHD, is present) ==\n")
audit <- wide %>%
  transmute(Region, target, ligand, interaction_name_2,
            prob_CON = round(prob_CON, 4), prob_NHD = round(prob_NHD, 4),
            dprob = round(delta_prob, 4), sig_any) %>%
  arrange(match(ligand, c("SPP1","GAS6","PROS1","PSAP","SEMA4D")), target, Region)
print(as.data.frame(audit), row.names = FALSE)

# ---- Assemble plot frame ---------------------------------------------------
present_regions <- intersect(REGIONS, unique(as.character(wide$Region)))
plot_df <- wide %>%
  mutate(Region = factor(Region, levels = present_regions),
         target = factor(target, levels = RECEIVERS),
         pathway_name = factor(pathway_name, levels = PATHWAY_ORDER),
         interaction_name_2 = factor(interaction_name_2,
                                     levels = rev(lr_rank$interaction_name_2)))

# ---------------------------------------------------------------------------
# receptor-detection GATE. Two things were wrong in the shipped panel:
#  (1) It drew GAS6-AXL, GAS6-TYRO3 and PROS1-TYRO3 for the ASTRO receiver. Neither AXL nor
#      TYRO3 appears in MAST_Astro_{Frontal,Hippo}.csv at all -- both fail the >=10% detection
#      pre-filter -- so those dots implied receptor expression we cannot measure. The MANIFEST
#      itself forbids claiming anything about AXL/TYRO3 in astrocytes.
#  (2) CellChat communication probability rises when the ligand rises, even if the receptor is
#      being lost: GAS6-MERTK reads +0.026 (red, "up in NHD") in Frontal while astrocyte MERTK
#      is delta -0.113 there and -0.302 in Hippocampus. Showing only the probability made the
#      panel argue for INCREASED TAM signalling -- the opposite of the figure's claim.
# Fix: gate every pair on whether its receptor is measurable in the receiver, mark the
# undetected ones with the house keyed cross instead of a probability dot, and annotate each
# row with the receptor's own MAST direction so the probability is read in context.
.rec_stat <- dplyr::bind_rows(lapply(present_regions, function(rg) {
  f <- file.path(PROJ, "tables", "mast_dual", sprintf("MAST_%s_%s.csv",
        if (RECV_SET == "astro") "Astro" else "Oligo", rg))
  if (!file.exists(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE) %>%
    transmute(Region = rg, receptor = gene, rec_delta = cliffs_delta,
              rec_padj = p_val_adj, rec_disc = discovery %in% c(TRUE, "TRUE"))
}))
# CellChat receptor names can be heteromeric COMPLEXES ("ITGAV_ITGB5"), which never match a
# single MAST gene -- a naive join flags them all as undetected. Split on "_" and take the
# worst-detected member: a complex is only measurable if every subunit clears the floor, and its
# direction is that of the member that moves most.
.rec_lookup <- function(rg, rcp) {
  parts <- strsplit(rcp, "_", fixed = TRUE)[[1]]
  hit   <- .rec_stat[.rec_stat$Region == rg & .rec_stat$receptor %in% parts, , drop = FALSE]
  if (nrow(hit) < length(parts)) return(data.frame(rec_delta = NA_real_, rec_disc = NA))
  hit[which.min(hit$rec_delta), c("rec_delta", "rec_disc")]
}
plot_df <- plot_df %>%
  mutate(.rl = purrr::map2(as.character(Region), as.character(receptor), .rec_lookup),
         rec_delta = vapply(.rl, function(x) x$rec_delta[1], numeric(1)),
         rec_disc  = vapply(.rl, function(x) isTRUE(x$rec_disc[1]), logical(1)),
         rec_undetected = is.na(rec_delta)) %>%
  select(-.rl)
cat(sprintf("receptor gate: %d of %d pair x region cells have an UNDETECTED receptor -> %s\n",
            sum(plot_df$rec_undetected), nrow(plot_df),
            paste(unique(plot_df$receptor[plot_df$rec_undetected]), collapse = ", ")))
cat("receptor direction carried onto the panel (MAST delta in the receiver):\n")
print(as.data.frame(plot_df %>% filter(!rec_undetected) %>%
        distinct(Region, receptor, rec_delta, rec_disc) %>% arrange(Region, rec_delta)),
      row.names = FALSE, digits = 3)

# Drop pairs whose receptor is unmeasurable in every region ("if not
# measurable or absent in both regions lets just erase it"). A row that is a cross in both
# columns carries no information and only costs vertical space -- AXL and TYRO3 never clear the
# detection floor in astrocytes at all, so GAS6-AXL / GAS6-TYRO3 / PROS1-TYRO3 (and
# SPP1-(ITGAV+ITGB1), whose ITGB1 subunit is likewise absent) are removed outright.
# The cross is kept for pairs measurable in one region but not the other, where the contrast is
# real information; the key is drawn only if such a case survives.
.always_nd <- plot_df %>% group_by(interaction_name_2) %>%
  summarise(all_nd = all(rec_undetected), .groups = "drop")
.dropped <- .always_nd$interaction_name_2[.always_nd$all_nd]
if (length(.dropped))
  cat(sprintf("dropped %d L-R pair(s) whose receptor is undetected in EVERY region: %s\n",
              length(.dropped), paste(.dropped, collapse = ", ")))
plot_df <- plot_df %>% filter(!interaction_name_2 %in% .dropped) %>%
  mutate(interaction_name_2 = droplevels(interaction_name_2))
ANY_ND <- any(plot_df$rec_undetected)
cat(sprintf("partial-detection crosses remaining: %d (key %s)\n",
            sum(plot_df$rec_undetected), if (ANY_ND) "drawn" else "suppressed"))

# hard 98% symmetric cap on dprob (matches reference)
cap_abs <- as.numeric(quantile(abs(plot_df$delta_prob), 0.98, na.rm = TRUE))
plot_df$delta_clip <- pmin(pmax(plot_df$delta_prob, -cap_abs), cap_abs)

write.csv(plot_df %>% select(Region, target, interaction_name_2, pathway_name,
                             prob_CON, prob_NHD, delta_prob, pval_CON, pval_NHD,
                             q_BH_CON, q_BH_NHD, sig_any),
          file.path(TDIR, paste0("fig3b_LR_stats", .stats_suffix, ".csv")), row.names = FALSE)

pA <- ggplot(plot_df, aes(x = target, y = interaction_name_2,
                          size = max_prob, fill = delta_clip)) +
  geom_point(data = function(d) dplyr::filter(d, !rec_undetected),
             shape = 21, color = "grey25", stroke = 0.25) +
  # receptor not measurable in this receiver -> keyed grey cross, never a probability dot
  (if (ANY_ND) list(
     geom_point(data = function(d) dplyr::filter(d, rec_undetected),
                aes(shape = "receptor n.d."), inherit.aes = TRUE,
                colour = "grey55", size = 1.5, stroke = 0.5),
     scale_shape_manual(values = c("receptor n.d." = 4), name = NULL)) else NULL) +
  # Ring the dots whose receptor is significantly down in the receiver, so a red (ligand-driven)
  # probability cannot be misread as intact signalling
  geom_point(data = function(d) dplyr::filter(d, !rec_undetected, rec_disc, rec_delta < 0),
             shape = 21, fill = NA, colour = "black", stroke = 0.55,
             size = 3.1, inherit.aes = TRUE) +
  # These panels never
  # set scale_x_discrete, so ggplot's default discrete expansion (0.6 units either side) padded
  # every facet far wider than its dot column -- with a single receiver that is ~0.25 in of white
  # on each side of every dot. Tighten the padding to the dots.
  scale_x_discrete(expand = expansion(add = 0.30)) +
  facet_grid2(pathway_name ~ Region, scales = "free_y", space = "free_y",
              switch = "y",
              labeller = labeller(Region = as_labeller(REGION_FULL),
                                  pathway_name = as_labeller(PATHWAY_LABELS)),
              # clip="off": after the Astro / OPC-Oligo split these facets are
              # only ~0.4 in wide, so the CENTRED region banner truncated to "Hippocampu".
              # theme(strip.clip=) does not reach ggh4x strip_themed, which is why the shared
              # theme's strip_region_x(, clip = "off") carries this argument -- it lets the label overflow
              # the pale banner instead of being cut, keeping the full region name.
              strip = strip_region_x(present_regions, fontsize = 8.5, clip = "off")) +
  # Measured: the green fill was
  # 0.532 in but "Hippocampus" at 7 pt needs ~0.66 in, so the label overflowed 0.068 in past
  # its own banner (clip="off" turns truncation into overflow -- it is not a fix). Pin the
  # facet to 0.70 in so the banner contains its label, and so this panel and its oligo sibling
  # (RECV_SET=oligo -> 4b) share identical facet geometry.
  # Banner type 7 -> 8.5 pt and the
  # facet 0.70 -> 0.92 in. At 0.70 in "Hippocampus" fitted but filled the pale fill edge-to-edge
  # with no padding, and the panel is placed at ~0.82x in the composite, so 7 pt landed at ~5.7 pt.
  # 8.5 pt at 0.92 in gives a readable label with real margin inside its banner.
  ggh4x::force_panelsizes(cols = unit(0.92, "in")) +
  scale_fill_gradient2(low = "#3E7CB1", mid = "grey96", high = "#D1495B",
                       midpoint = 0, name = bquote(Delta ~ "prob"),
                       limits = c(-cap_abs, cap_abs),
                       # Breaks were round(c(-cap,0,cap), 2). cap_abs is
                       # ~0.025, so rounding pushed the outer breaks outward to +-0.03 -- outside
                       # limits = c(-cap_abs, cap_abs) -- and ggplot silently dropped them, leaving
                       # a colourbar whose only tick read "+0.00". The reader could not tell whether
                       # full crimson meant 0.005 or 0.05. Use the limits themselves as breaks and
                       # label to 3 dp (these probabilities are ~1e-3-1e-2).
                       breaks = c(-cap_abs, 0, cap_abs),
                       labels = function(x) sprintf("%+.3f", x),
                       guide = guide_colorbar(order = 1, barheight = unit(1.6, "cm"),
                                              barwidth = unit(0.28, "cm"))) +
  scale_size_continuous(name = "max prob.", range = c(0.8, 3.6),
                        # fixed breaks so the size key is identical in the astro and
                        # oligo runs (pretty_breaks resolved to one key in astro, two
                        # in oligo -- a one-key size legend conveys no scale).
                        breaks = c(0.01, 0.03, 0.05),
                        guide = guide_legend(order = 2, nrow = 3,
                                             override.aes = list(fill = "grey50"))) +
  # Sender declared on the y axis (propagated across every l-r panel).
  # The rows are ligand - receptor pairs, so the axis names both halves; the grey strips
  # beside them are CellChat pathway names, which are not always the ligand (TENM4 sits
  # under "ADGRL", named for the receptor; VSIR under "vista"; GAS6 under "GAS"). Naming
  # the axis for the pair — and the x axis for the receiver — lets the reader get sender
  # and receiver from the panel instead of the legend.
  # Never a unicode arrow here: base pdf() cannot encode U+2192 and it has already broken
  # a render in this project. ASCII hyphen only.
  labs(x = "Receiver cell", y = "Microglial ligand - receptor") +
  theme_pub(base_size = 8) +
  theme(
    axis.text.y       = element_text(size = 6.3, face = "italic", color = "black"),
    axis.text.x       = element_text(size = 6.3, color = "black",
                                     angle = 45, hjust = 1, vjust = 1),
    axis.title.x      = element_text(size = 7, margin = margin(t = 3)),
    axis.line         = element_blank(),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.spacing.x   = unit(0.25, "lines"),
    panel.spacing.y   = unit(0.10, "lines"),
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 6.4,
                                     face = "plain", color = "grey15",
                                     margin = margin(r = 3, l = 3, t = 1, b = 1)),
    # x region banner (fill + black text) comes from strip_region_x; leave it.
    # round-5 visual fix: unify the y-strip (ligand-family)
    # element_rect with the 3c sibling so the adjacently-placed dotplots read as an exact
    # pair. 3b was grey94/grey75/0.3 (darker/heavier); matched to 3c's lighter
    # grey97/grey80/0.25 so the black strip text stays clean on both.
    strip.background.y = element_rect(fill = "grey97", color = "grey80", linewidth = 0.25),
    strip.placement   = "outside",
    legend.position   = "right",
    legend.box        = "vertical",
    # The colourbar and size keys sat on the right panel border and
    # "-0.082" touched "max prob.". Positive box spacing + real spacing between the two guides.
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0.30, "cm"),
    legend.spacing.y  = unit(0.40, "cm"),
    legend.key.size   = unit(0.22, "cm"),
    legend.text       = element_text(size = 6.2),
    legend.title      = element_text(size = 6.2),
    plot.margin       = margin(6, 4, 4, 16))

# Width FLOOR: the left row-group strip column is a fixed width, so narrowing
# the panel starves the facets, not the strip. At 2.60 in the region banners clipped to
# "Fronta"/"ocam" -- the house rule is that facet strip titles must fit. 3.20 (1 receiver)
# and 3.45 (2 receivers) are the practical floors with full region names.
# W is measured from the gtable (fixed 0.92-in facets + strips + labels + legend). The hard-coded
# 4.20 / 4.35 in was narrower than the layout needs, and a gtable with fixed panel widths does not shrink —
# it overflowed, which put the colourbar and size keys on the right panel border (Fig 3g by eye).
.g <- ggplot2::ggplotGrob(pA)
W <- sum(grid::convertWidth(.g$widths, "in", valueOnly = TRUE)) + 0.08
H <- 2.55   # +0.15 in for the sender header
cat(sprintf("measured panel width: %.2f in\n", W))
# Sender marker — width-aware so the arrow neither collides with the label nor
# drifts away from it (22_publication_theme_FH.R).
pA <- with_sender(pA, "Microglia", panel_width = W,
                  arrow_x = if (RECV_SET == "astro") 0.97 else 0.956)
BN <- if (RECV_SET == "astro") "F3g2_micro_astro_LR" else "4b_micro_oligo_LR"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), pA, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), pA, width = W, height = H, dpi = 600)
cat(sprintf("\nWrote %s.{pdf,png} + fig3b_LR_stats.csv\n", BN))
cat("=== DONE ===\n", file = stderr())
