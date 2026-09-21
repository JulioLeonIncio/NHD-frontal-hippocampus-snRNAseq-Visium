#!/usr/bin/env Rscript
# =============================================================================
# 4G_astro_neuron_LR_FH.R — Figure 4 (Part B): astrocyte -> Ex/Inh neuron ligand-receptor pair bubble, frontal + hippocampus. Port of the 4E-grammar
# astro->neuron bubble in manuscript_7fig/scripts/09e_astro_neuron_3panel.R.
# -----------------------------------------------------------------------------
# Same visual grammar as F5e2_LR_pairs_ExInh (vertical: y = LR pair, x = Condition,
# facet class x Region, fill = signed delta prob, size = comm prob, thick border
# = BH-FDR q<0.1). The message: the curated astrocytic synapse-organiser + trophic
# inputs (NRXN1/NCAM1/CADM1 organisers, NRG3/NRG2->ERBB4 trophic) are broadly
# preserved in NHD — localizing the NHD perturbation upstream to microglia, not
# to astrocyte->neuron support. (Hippocampus is a single lane -> "(low n)".)
#
# FH data source: withNeurons_v2 CellChat objects (major-class Neuron_Ex/Inh
# targets — a direct Neuron_Ex->Ex / Neuron_Inh->Inh rename; no subtype roll-up).
#
# Critical gotcha: keyed by interaction_name_2 (ligand-RECEPTOR), not pathway_name.
# No causal verbs; the delta is the effect size.
#
# Output: figures/Figure_5/panels/4G_astro_neuron_LR.{pdf,png} + stats CSV.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(patchwork)
  library(CellChat); library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales); library(forcats); library(ggh4x)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))

CC    <- file.path(PROJ, "data", "cellchat_per_region")
# Astrocyte is the SENDER here, so the panel belongs with the astrocyte figure. Moved from
# Figure_5 to Figure_3.
# 14:51-14:56 -- after the FH atlas (14:08), inside the rebuild, Frontal+Hippo only. It is
# not the old 3-region analysis (manuscript_7fig/09e_v2_astro_to_ExInh_LR_bubble.R).
PANEL <- file.path(PROJ, "figures", "Figure_3", "panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

REGIONS <- c("Frontal","Hippo")
LOW_CONF_REGIONS <- "Hippo"
SENDER  <- "Astro"
NEURON  <- c("Neuron_Ex","Neuron_Inh")

# Curated astro->neuron "preserved support" whitelist (FH-audited present pairs):
# NRXN1/NCAM1/CADM1 synapse organisers + NRG3/NRG2->ERBB4 trophic (Inh-only).
LR_WHITELIST <- c(
  "NRG3 - ERBB4", "NRG2 - ERBB4",
  "NRXN1 - NLGN1", "NRXN1 - LRRTM4", "NRXN1 - CLSTN1",
  "NCAM1 - NCAM1",
  "CADM1 - CADM1", "CADM1 - NECTIN3")

load_region <- function(region) {
  path <- file.path(CC, sprintf("merged_NHDvsCON_%s_withNeurons_v2.rds", region))
  if (!file.exists(path))
    stop(sprintf("MISSING %s — CellChat withNeurons prep for %s failed", path, region))
  m  <- readRDS(path)
  df <- subsetCommunication(m, slot.name = "net")
  bind_rows(df$CON %>% mutate(Condition = "CON"),
            df$NHD %>% mutate(Condition = "NHD")) %>% mutate(Region = region)
}

raw <- bind_rows(lapply(REGIONS, load_region)) %>% as_tibble() %>%
  filter(source == SENDER, target %in% NEURON) %>%
  mutate(class = dplyr::recode(as.character(target),
                               Neuron_Ex = "Ex", Neuron_Inh = "Inh"))
stopifnot(nrow(raw) > 0)
cat(sprintf("Raw astro->neuron rows: %d\n", nrow(raw)))

agg <- raw %>%
  group_by(Region, class, interaction_name_2, pathway_name, Condition) %>%
  summarise(prob_sum = sum(prob, na.rm = TRUE),
            min_pval = suppressWarnings(min(pval, na.rm = TRUE)), .groups = "drop") %>%
  mutate(min_pval = ifelse(is.infinite(min_pval), 1, min_pval)) %>%
  pivot_wider(id_cols = c(Region, class, interaction_name_2, pathway_name),
              names_from = Condition, values_from = c(prob_sum, min_pval),
              values_fill = list(prob_sum = 0, min_pval = 1)) %>%
  mutate(min_p_either = pmin(min_pval_CON, min_pval_NHD),
         q_either     = p.adjust(min_p_either, method = "BH"),
         delta        = prob_sum_NHD - prob_sum_CON,
         peak         = pmax(prob_sum_CON, prob_sum_NHD),
         ok           = q_either < 0.1)

agg %>% group_by(interaction_name_2, pathway_name) %>%
  summarise(net_delta = sum(delta), max_abs = max(abs(delta)),
            min_q = min(q_either), n_ctx = dplyr::n(), .groups = "drop") %>%
  arrange(desc(max_abs)) %>%
  write.csv(file.path(TDIR, "fig4G_present_astro_neuron_LR.csv"), row.names = FALSE)

present_lr <- sort(unique(agg$interaction_name_2))
absent_lr  <- setdiff(LR_WHITELIST, present_lr)
if (length(absent_lr) > 0)
  cat("WARNING whitelist pair(s) absent, dropping:", paste(absent_lr, collapse = " | "), "\n")
LR_WHITELIST <- intersect(LR_WHITELIST, present_lr)
stopifnot(length(LR_WHITELIST) >= 1)
cat("Curated astro->neuron LR pairs used (", length(LR_WHITELIST), "):",
    paste(LR_WHITELIST, collapse = " | "), "\n")

# preservation audit (the requested "broadly preserved" check)
cat("\n== astro->neuron preserved-support status (prob CON/NHD, delta) ==\n")
audit <- agg %>% filter(interaction_name_2 %in% LR_WHITELIST) %>%
  transmute(Region, class, interaction_name_2,
            prob_CON = round(prob_sum_CON, 4), prob_NHD = round(prob_sum_NHD, 4),
            delta = round(delta, 4), q = round(q_either, 3)) %>%
  arrange(interaction_name_2, Region, class)
print(as.data.frame(audit), row.names = FALSE)
# summarise Frontal (2v2 lane) preservation vs Hippo (single lane)
frontal_preserved <- agg %>% filter(interaction_name_2 %in% LR_WHITELIST, Region == "Frontal") %>%
  summarise(n = dplyr::n(), n_up_or_stable = sum(delta >= -0.05, na.rm = TRUE))
cat(sprintf("\nFrontal astro->neuron pairs with delta >= -0.05 (preserved/up): %d of %d\n",
            frontal_preserved$n_up_or_stable, frontal_preserved$n))

bubble <- agg %>%
  filter(interaction_name_2 %in% LR_WHITELIST) %>%
  pivot_longer(cols = c(prob_sum_CON, prob_sum_NHD),
               names_to = "Condition", values_to = "prob") %>%
  mutate(Condition = sub("prob_sum_", "", Condition),
         class     = factor(class,     levels = c("Ex","Inh")),
         Region    = factor(Region,    levels = REGIONS),
         Condition = factor(Condition, levels = c("CON","NHD")))

PATHWAY_ORDER_AN <- c("NRXN", "NCAM", "CADM", "NRG")
bubble <- bubble %>%
  mutate(pathway_name = factor(pathway_name,
                               levels = intersect(PATHWAY_ORDER_AN, unique(pathway_name))))
# ---------------------------------------------------------------------------
# receptor-detection GATE — same treatment as the micro->astro l-r panel, for
# consistency. CellChat probability rises with the ligand even when the receptor is
# lost, so each pair is checked against the receiver's own MAST table:
#   * receptor not in the table at all (fails the >=10% floor)  -> keyed grey cross, no dot
#   * receptor down at discovery level (|delta|>=0.15, padj<.05) -> black ring on the dot
# The receiver here is class-specific (Neuron_Ex vs Neuron_Inh), so the lookup is per
# region and per class. Receptor is parsed from "ligand - receptor".
.rec_tab <- dplyr::bind_rows(lapply(c("Neuron_Ex","Neuron_Inh"), function(cl)
  dplyr::bind_rows(lapply(REGIONS, function(rg) {
    f <- file.path(PROJ, "tables", "mast_dual", sprintf("MAST_%s_%s.csv", cl, rg))
    if (!file.exists(f)) return(NULL)
    read.csv(f, stringsAsFactors = FALSE) %>%
      transmute(Region = rg, class = sub("^Neuron_", "", cl), receptor = gene,
                rec_delta = cliffs_delta, rec_disc = discovery %in% c(TRUE, "TRUE"))
  }))))
bubble <- bubble %>%
  mutate(receptor = trimws(sub("^.*\\s-\\s", "", as.character(interaction_name_2)))) %>%
  left_join(.rec_tab, by = c("Region", "class", "receptor")) %>%
  mutate(rec_undetected = is.na(rec_delta))
.all_nd <- bubble %>% group_by(interaction_name_2) %>%
  summarise(all_nd = all(rec_undetected), .groups = "drop")
.drop <- .all_nd$interaction_name_2[.all_nd$all_nd]
if (length(.drop)) {
  cat(sprintf("dropped %d pair(s) undetected in EVERY receiver: %s\n",
              length(.drop), paste(.drop, collapse = ", ")))
  bubble <- bubble %>% filter(!interaction_name_2 %in% .drop)
}
ANY_ND <- any(bubble$rec_undetected)
cat(sprintf("receptor gate: %d cross(es) [%s] | %d ring(s) for receptor DOWN at discovery\n",
            sum(bubble$rec_undetected),
            paste(unique(paste0(bubble$receptor[bubble$rec_undetected], "/",
                                bubble$class[bubble$rec_undetected])), collapse = ", "),
            sum(!bubble$rec_undetected & bubble$rec_disc & bubble$rec_delta < 0)))

y_order <- bubble %>% group_by(interaction_name_2) %>%
  summarise(net_delta = sum(delta, na.rm = TRUE), .groups = "drop") %>%
  arrange(net_delta) %>% pull(interaction_name_2)
bubble$interaction_name_2 <- factor(bubble$interaction_name_2, levels = unique(y_order))

# Signed-delta colour cap: data-driven 95th pct (stated). Raw delta kept in CSV.
delta_cells <- bubble %>% distinct(Region, class, interaction_name_2, delta)
delta_lim   <- as.numeric(stats::quantile(abs(delta_cells$delta), 0.95, na.rm = TRUE))
if (!is.finite(delta_lim) || delta_lim == 0) delta_lim <- max(abs(delta_cells$delta), na.rm = TRUE)
n_wins <- sum(abs(delta_cells$delta) > delta_lim, na.rm = TRUE)
cat(sprintf("\nbubble: %d dots / %d pairs; delta colour cap |%.3f| (95th pct), %d winsorized for colour\n",
            nrow(bubble), length(unique(bubble$interaction_name_2)), delta_lim, n_wins))

write.csv(bubble, file.path(TDIR, "fig4G_astro_neuron_LR_bubble.csv"), row.names = FALSE)

present_regions <- REGIONS[REGIONS %in% as.character(bubble$Region)]

p <- ggplot(bubble,
            aes(x = Condition, y = interaction_name_2, size = prob,
                fill = pmax(pmin(delta, delta_lim), -delta_lim))) +
  # probability dots only where the receptor is measurable in that receiver
  geom_point(data = function(d) dplyr::filter(d, !rec_undetected),
             shape = 21, stroke = 0.25, color = "grey20", alpha = 0.95) +
  # receptor absent from the receiver's MAST table -> keyed cross, never a probability dot
  (if (ANY_ND) list(
     geom_point(data = function(d) dplyr::filter(d, rec_undetected),
                aes(shape = "receptor n.d."), size = 1.5, stroke = 0.5,
                colour = "grey55", inherit.aes = TRUE),
     scale_shape_manual(values = c("receptor n.d." = 4), name = NULL)) else NULL) +
  # Receptor significantly down in the receiver -> ring, so a ligand-driven probability
  # cannot be misread as intact signalling (same encoding as the micro->astro panel)
  geom_point(data = function(d) dplyr::filter(d, !rec_undetected, rec_disc, rec_delta < 0),
             shape = 21, fill = NA, colour = "black", stroke = 0.55, size = 3.1,
             inherit.aes = TRUE, show.legend = FALSE) +
  # Source table check: ok is TRUE for 52/52 rows and q_either is exactly 0 for
  # every row, so the ring discriminated nothing, had no legend key, and its 0.55 stroke swamped
  # the pale fills. Restore only with a real gate and a key.
  # Horizontal. Was class (Ex/Inh) as ROWS stacked under each other,
  # which made a tall narrow panel; now Region x class runs across the top so the LR pairs
  # share one y column and the panel is landscape.
  facet_grid2(pathway_name ~ Region + class, scales = "free", space = "free",
              switch = "y",
              labeller = labeller(Region = region_full_lowconf_labeller(LOW_CONF_REGIONS)),
              # strip_region_x(, clip = "off") supplies one fill per region, but
              # `Region + class` builds a NESTED strip with four columns, so the two fills
              # cycled -- "Frontal Inh" rendered in the Hippocampus green and "Hippocampus Ex"
              # in the Frontal violet, i.e. the banner colour contradicted its own label.
              # Build the strip explicitly: each region tone repeated once per class for the
              # OUTER (region) layer, then a neutral pale grey for the INNER (class) layer.
              strip = ggh4x::strip_themed(
                clip = "off",
                background_x = ggh4x::elem_list_rect(
                  fill = c(rep(unname(PAL_REGION_PALE[present_regions]),
                               each = dplyr::n_distinct(bubble$class)),
                           rep("grey96", length(present_regions) *
                                         dplyr::n_distinct(bubble$class))),
                  colour = NA),
                text_x = ggh4x::elem_list_text(colour = "black", size = 7, face = "plain"))) +
  # House diverging tokens to match 3b/3c/4b/4c, which encode the same
  # quantity; PAL_DGE rendered brick-brown beside them. mid grey96 like the siblings.
  scale_fill_gradient2(low = "#3E7CB1", mid = "grey96",
                       high = "#D1495B", midpoint = 0,
                       limits = c(-delta_lim, delta_lim), oob = scales::squish,
                       name = expression(Delta ~ "prob (NHD" - "CON)"),
                       guide = guide_colorbar(order = 1, barheight = unit(1.6, "cm"),
                                              barwidth = unit(0.28, "cm"))) +
  scale_size_continuous(range = c(0.8, 3.6), name = "comm prob",
                        breaks = c(0.10, 0.25, 0.45),
                        guide = guide_legend(order = 2)) +
  # Sender declared on the y axis (propagated across every l-r panel).
  # The rows are ligand - receptor pairs, so the axis names both halves; the grey strips
  # beside them are CellChat pathway names, which are not always the ligand (TENM4 sits
  # under "ADGRL", named for the receptor; VSIR under "vista"; GAS6 under "GAS"). Naming
  # the axis for the pair — and the x axis for the receiver — lets the reader get sender
  # and receiver from the panel instead of the legend.
  # Never a unicode arrow here: base pdf() cannot encode U+2192 and it has already broken
  # a render in this project. ASCII hyphen only.
  labs(x = NULL, y = "Astrocytic ligand - receptor") +
  theme_pub(base_size = 8) +
  theme(  # y-strip styling now set below, with the pathway row-groups
        # Strip text was 8.0 pt
        # grey15 here vs 6.4 pt black in 3b -- a 25% type mismatch between two sibling l-r panels.
        axis.text.x        = element_text(size = 6.3),
        axis.text.y        = element_text(size = 6.3, face = "italic"),
        # this panel never blanked axis.line, so theme_pub's default axis line drew a heavy
        # black L down the left and across the bottom -- while its sibling l-r panels
        # (3b/3c/4b/4c) blank the axis line and rely on a thin black panel.border. It also
        # used a grey70 border where the siblings use black. Both now match the siblings.
        # row-group strips styled to match 3b/3c exactly (grey97 + grey80 hairline, 6.4 pt black)
        strip.background.y = element_rect(fill = "grey97", colour = "grey80", linewidth = 0.25),
        strip.text.y.left  = element_text(angle = 0, hjust = 0.5, size = 6.4,
                                          colour = "black", face = "plain"),
        strip.placement    = "outside",
        axis.line          = element_blank(),
        panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x    = unit(0.20, "lines"),
        panel.spacing.y    = unit(0.20, "lines"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.key.size    = unit(0.22, "cm"),
        legend.title       = element_text(size = 6.2),
        legend.text        = element_text(size = 6.2),
        plot.margin        = margin(4, 4, 4, 4),
        )

W <- 5.90; H <- 2.7   # + left row-group strip column and per-family row gaps   # +0.15 in for the sender header
# Sender marker — width-aware so the arrow neither collides with the label nor
# drifts away from it (22_publication_theme_FH.R).
p <- with_sender(p, "Astrocytes", panel_width = W, arrow_x = 0.542)
BN <- "F3h2_astro_neuron_LR"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nWrote %s.{pdf,png} + fig4G_astro_neuron_LR_bubble.csv\n", BN))
cat("=== DONE ===\n", file = stderr())
