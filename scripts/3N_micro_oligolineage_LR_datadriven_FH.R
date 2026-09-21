#!/usr/bin/env Rscript
# =============================================================================
# 3N_micro_oligolineage_LR_datadriven_FH.R — Figure 4, NEW panel:
# "F4j_micro_oligolineage_LR".  Microglia -> {OPC, Oligodendrocyte}
# ligand-receptor rewiring, frontal + HIPPOCAMPUS (no OCC).
# -----------------------------------------------------------------------------
# Why, given 4b already exists.  `3B_micro_glia_LR_FH.R` (RECV_SET=oligo) shows a
# curated whitelist — SPP1 / GAS6 / PROS1 / PSAP — carried over from the microglial
# secretome of Figure 2.  That answers "what happens to the Figure-2 ligands", not
# This
# panel is the unbiased companion: every Micro-PVM -> OPC/Oligo pair that CellChat
# calls significant in at least one condition, ranked by |NHD - CON| probability.
#
# It surfaces two axes the curated list could not:
#   * PDGF -> PDGFRA, the canonical OPC maintenance/proliferation signal: lost in
#     the frontal cortex (PDGFC-PDGFRA 0.026 -> 0.000) and GAINED in hippocampus
#     (PDGFB-PDGFRA 0.009 -> 0.045).  The receptor half moves with it — PDGFRA is
#     down in NHD OPC in both regions in the MAST tables (see 3g gene panel) — so
#     ligand and receptor fail together in the region that loses its myelin.
#   * TENM4 -> ADGRL3, lost in frontal (0.116 -> 0.003 on OPC).  TENM4 is required
#     for oligodendrocyte differentiation and small-axon myelination.
#
# Exclusion, stated not hidden.  CellChat's neurotransmitter entry
# "Glu-(SLC1A3+GLS)" dominates the raw ranking (11 of the top 20 frontal rows).
# Its "ligand" is inferred from SLC1A3 + GLS expression in the SENDER, and SLC1A3
# is an astrocyte gene: microglia scoring as glutamate senders is ambient/neighbour
# signal, not microglial biology.  Every Glu- row is dropped, counted, and written
# to the source CSV so the exclusion is auditable rather than silent.
#
# Sibling styling to 4b (pathway row strips, pale region banners, full region
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(forcats); library(ggh4x); library(stringr); library(patchwork)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))

LR    <- file.path(PROJ, "tables", "cellchat_all_LR_pairs.csv")
PANEL <- file.path(PROJ, "figures", "Figure_4", "panels")
TDIR  <- file.path(PROJ, "tables")
LOGD  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING CellChat L-R table — run 50_cellchat_LR_FH.R" = file.exists(LR))

SENDER    <- "Micro-PVM"
RECEIVERS <- c("OPC", "Oligo")
PVAL_CUT  <- 0.05
TOP_N     <- 9          # pairs shown, chosen on the union of both regions

raw <- read.csv(LR, stringsAsFactors = FALSE)
# The rebuild's CellChat run is Frontal + Hippocampus only; assert rather than assume.
stopifnot("CellChat table is NOT OCC-free — wrong input" =
            !any(grepl("OCC", raw$region_cond, ignore.case = TRUE)))
cat("region_cond present:", paste(sort(unique(raw$region_cond)), collapse = ", "), "\n")

lr <- raw %>%
  filter(source == SENDER, target %in% RECEIVERS) %>%
  separate(region_cond, into = c("Condition", "Region"), sep = "_", remove = FALSE)
cat(sprintf("%s -> %s rows: %d\n", SENDER, paste(RECEIVERS, collapse = "/"), nrow(lr)))

# ---- stated exclusion: CellChat's astrocyte-driven glutamate pseudo-ligand ---
is_glu <- grepl("^Glu-", lr$interaction_name_2)
n_glu  <- dplyr::n_distinct(lr$interaction_name_2[is_glu])
cat(sprintf("EXCLUDED %d glutamate pseudo-ligand pairs (Glu-(SLC1A3+GLS); SLC1A3 is astrocytic — ambient sender signal)\n", n_glu))
glu_dropped <- lr %>% filter(is_glu) %>%
  distinct(interaction_name_2, target, Region, Condition, prob, pval)
lr <- lr[!is_glu, ]

# ---- CON vs NHD per pair x receiver x region --------------------------------
wide <- lr %>%
  select(interaction_name_2, pathway_name, target, Region, Condition, prob, pval) %>%
  pivot_wider(names_from = Condition, values_from = c(prob, pval),
              values_fill = list(prob = 0, pval = 1)) %>%
  mutate(dprob   = prob_NHD - prob_CON,
         maxprob = pmax(prob_CON, prob_NHD),
         sig_any = pmin(pval_CON, pval_NHD) < PVAL_CUT,
         Region  = factor(Region, levels = REGION_ORDER),
         target  = factor(target, levels = RECEIVERS))

# Rank on the union of regions so both panels show the same rows (a pair that moves
# in one region only is still informative in the other, where it is flat)
rank_tab <- wide %>% filter(sig_any) %>%
  group_by(interaction_name_2, pathway_name) %>%
  summarise(score = max(abs(dprob)), .groups = "drop") %>%
  arrange(desc(score))
keep_pairs <- head(rank_tab$interaction_name_2, TOP_N)
cat("\nkept pairs (ranked on max |dprob| across regions/receivers):\n")
print(as.data.frame(rank_tab %>% filter(interaction_name_2 %in% keep_pairs) %>%
                    mutate(score = round(score, 4))))

dat <- wide %>% filter(interaction_name_2 %in% keep_pairs) %>%
  mutate(interaction_name_2 = factor(interaction_name_2, levels = rev(keep_pairs)),
         pathway_name = factor(pathway_name,
                               levels = unique(rank_tab$pathway_name[
                                 match(keep_pairs, rank_tab$interaction_name_2)])))
# a pair CellChat never scored in a given region/receiver is a true zero (both
# conditions absent) — keep it visible as an open cross rather than a blank cell
grid_full <- expand_grid(interaction_name_2 = factor(keep_pairs, levels = rev(keep_pairs)),
                         target = factor(RECEIVERS, levels = RECEIVERS),
                         Region = factor(REGION_ORDER, levels = REGION_ORDER))
dat <- grid_full %>%
  left_join(dat, by = c("interaction_name_2","target","Region")) %>%
  group_by(interaction_name_2) %>%
  mutate(pathway_name = first(na.omit(pathway_name))) %>% ungroup()
absent <- dat %>% filter(is.na(dprob) | maxprob == 0)
dat_pt <- dat %>% filter(!is.na(dprob), maxprob > 0)
cat(sprintf("\ncells: %d drawn, %d not scored in that region/receiver (open cross)\n",
            nrow(dat_pt), nrow(absent)))

# Colour scale clipped to a stated robust quantile (house idiom, same as 3G): the
# TENM4 collapse (-0.113) is ~4x the next largest change and would flatten every
# other pair to white if the scale ran to the maximum. Clipped values saturate.
LIMD <- max(quantile(abs(dat_pt$dprob), 0.85, na.rm = TRUE), 0.02)
n_clip <- sum(abs(dat_pt$dprob) > LIMD)
cat(sprintf("fill scale: symmetric +/- %.3f (85th pct |dprob|); %d cells clipped (max |dprob| %.3f)\n",
            LIMD, n_clip, max(abs(dat_pt$dprob))))
dat_pt$dprob_c <- pmax(pmin(dat_pt$dprob, LIMD), -LIMD)
p <- ggplot(dat_pt, aes(x = target, y = interaction_name_2)) +
  geom_point(data = absent, aes(x = target, y = interaction_name_2),
             inherit.aes = FALSE, shape = 4, size = 1.5, colour = "grey65", stroke = 0.5) +
  geom_point(aes(fill = dprob_c, size = maxprob), shape = 21, colour = "grey25", stroke = 0.3) +
  scale_fill_gradient2(low = "#3E7CB1", mid = "white", high = "#D1495B", midpoint = 0,
                       limits = c(-LIMD, LIMD),
                       breaks = c(-LIMD, 0, LIMD),
                       labels = function(b) sprintf("%+.3f", b),   # outer labels are the clip, not the max

                       name = expression(Delta*" prob (NHD - CON)")) +
  scale_size_continuous(range = c(0.9, 4.2), name = "max prob.") +
  facet_grid2(pathway_name ~ Region, scales = "free_y", space = "free_y", switch = "y",
              labeller = labeller(Region = region_full_lowconf_labeller()),
              strip = strip_region_x(REGION_ORDER, fontsize = 6.8, clip = "off")) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(0.22, "cm"),
                                barheight = unit(1.2, "cm")),
         size = guide_legend(order = 2)) +
  # Sender DECLARED ("we don't know the sender at first glance").
  # The panel stated its receivers ("Receiver cell") but never its sender, so the reader
  # had to get it from the legend. Naming the y axis for the sender closes that in the
  # panel's own grammar: microglial ligand (rows) -> receiver cell (columns). A drawn
  # microglia pictogram was considered and not used — no icon assets exist in the repo,
  # no other panel in Figs 2-4 uses pictograms, and a hand-rolled cell would read as
  # crude beside typographic siblings; a proper vector icon belongs in Illustrator.
  # Never a unicode arrow here — base pdf() cannot encode U+2192 (it has already
  # broken a render in this project); ASCII "to" only.
  # Sender declared on the y axis (propagated across every l-r panel).
  # The rows are ligand - receptor pairs, so the axis names both halves; the grey strips
  # beside them are CellChat pathway names, which are not always the ligand (TENM4 sits
  # under "ADGRL", named for the receptor; VSIR under "vista"; GAS6 under "GAS"). Naming
  # the axis for the pair — and the x axis for the receiver — lets the reader get sender
  # and receiver from the panel instead of the legend.
  # Never a unicode arrow here: base pdf() cannot encode U+2192 and it has already broken
  # a render in this project. ASCII hyphen only.
  labs(x = "Receiver cell", y = "Ligand - receptor") +
  theme_pub(base_size = 8) +
  # Only the legend and spacing were loose here — the
  # axis type already matched — so the panel can carry the astro geometry at 4.35 in.
  theme(axis.text.y  = element_text(size = 6.3, colour = "black", face = "italic"),
        axis.text.x  = element_text(size = 6.3, colour = "black"),
        axis.title.x = element_text(size = 7, margin = margin(t = 3)),
        axis.title.y = element_text(size = 7, margin = margin(r = 3)),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        strip.text.y.left = element_text(size = 6.4, angle = 0, hjust = 0.5, colour = "black"),
        strip.background.y = element_rect(fill = "grey97", colour = "grey80", linewidth = 0.25),
        strip.placement = "outside",
        panel.spacing.x = unit(0.25, "lines"), panel.spacing.y = unit(0.10, "lines"),
        legend.text  = element_text(size = 6.2),
        legend.title = element_text(size = 6.2),
        legend.key.size = unit(0.22, "cm"))

# Sender marker — now via the shared house helper (22_publication_theme_FH.R) so this
# panel and every other l-r panel are byte-identical in geometry.
p <- with_sender(p, "Microglia", panel_width = 4.35, arrow_x = 0.927)

BN <- "F4j_micro_oligolineage_LR"
W  <- 4.35; H <- 2.85   # astro l-r geometry (4.35 x 2.40); +0.45 in height for 9 rows vs 6
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("Saved %s.{pdf,png} (%.2f x %.2f in)\n", BN, W, H))

# ---- source CSV: shown rows and the stated exclusions ------------------------
out <- bind_rows(
  wide %>% mutate(status = ifelse(interaction_name_2 %in% keep_pairs,
                                  "shown", "significant but outside top N")) %>%
    select(interaction_name_2, pathway_name, target, Region,
           prob_CON, prob_NHD, pval_CON, pval_NHD, dprob, maxprob, sig_any, status),
  glu_dropped %>% transmute(interaction_name_2, pathway_name = "Glutamate",
                            target, Region, prob_CON = NA, prob_NHD = NA,
                            pval_CON = NA, pval_NHD = NA, dprob = NA, maxprob = prob,
                            sig_any = pval < PVAL_CUT,
                            status = "EXCLUDED: astrocytic pseudo-ligand (SLC1A3+GLS)"))
write.csv(out, file.path(TDIR, "figF4j_micro_oligolineage_LR_FH.csv"), row.names = FALSE)
cat("Wrote tables/figF4j_micro_oligolineage_LR_FH.csv (", nrow(out), " rows)\n", sep = "")

cat("\n== shown pairs, per region/receiver ==\n")
print(as.data.frame(dat_pt %>% arrange(Region, target, desc(abs(dprob))) %>%
                    mutate(across(c(prob_CON, prob_NHD, dprob), ~round(.x, 4))) %>%
                    select(Region, target, interaction_name_2, prob_CON, prob_NHD, dprob)))

writeLines(capture.output(sessionInfo()),
           file.path(LOGD, "3N_micro_oligolineage_LR_datadriven_FH_sessionInfo.txt"))
cat("=== DONE ===\n")
