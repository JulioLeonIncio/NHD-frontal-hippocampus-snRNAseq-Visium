#!/usr/bin/env Rscript
# =============================================================================
# 43_fig1_def_MAST_2region_FH.R — 43_fig1_def_MAST_2region_FH.R (frontal + hippocampus rebuild)
# -----------------------------------------------------------------------------
# Rebuilds NHD Figure 1 panels d / e / f using per-nucleus MAST avg_log2FC,
# redesigned for two regions (Frontal, Hippocampus). Adapted from
# manuscript_7fig/scripts/43_fig1_def_MAST_3region.R.
#
# With only 2 regions the 3-region design changes:
#   * panel d: a per-cell-type Frontal(x) vs Hippo(y) log2FC scatter (single
#     region pair) — NHD direction concordance across the two regions.
#   * panel e: a bar per cell type of the (single) Frontal-vs-Hippo Pearson r.
#   * panel f: the old region-specificity index is degenerate with one region
#     pair, so it is REPLACED by Frontal-unique / shared / Hippo-unique discovery
#     feature counts per cell type (a meaningful "where is the NHD response
#     shared vs region-restricted" summary).
#
# Metric throughout = per-nucleus MAST avg_log2FC.
# "Discovery DEG" in a region = discovery==TRUE there, after the ambient filter.
#
# Standalone + idempotent + portable. Panels only. All legend numbers are keyable to fig1_def_MAST_2region_FH.csv.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(patchwork); library(ragg); library(ggrepel)
})
set.seed(42)

# --- portable project root (julio.l origin/RIKEN or JulioLeon local) --------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
if (is.na(PROJ)) stop("MISSING project root — neither julio.l nor JulioLeon path exists")
DIR        <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
THEME_R    <- file.path(DIR, "scripts", "22_publication_theme_FH.R")
ARTIFACT_R <- file.path(DIR, "scripts", "_artifact_genes.R")
DATA_CSV   <- file.path(DIR, "tables", "mast_dual", "MAST_dual_discovery_all.csv")
PANEL_DIR  <- file.path(DIR, "figures", "Figure_1", "panels")
SUPP_DIR   <- file.path(PANEL_DIR, "_supp")   # demoted panels (e -> supplement)
TAB_DIR    <- file.path(DIR, "tables", "mast_dual")
OUT_CSV    <- file.path(TAB_DIR, "fig1_def_MAST_2region_FH.csv")
LOG_DIR    <- file.path(DIR, "logs")

# --- fail loud, fail early --------------------------------------------------
for (p in c(THEME_R, ARTIFACT_R, DATA_CSV)) {
  if (!file.exists(p)) stop(sprintf("MISSING %s — prep step failed", p))
}
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR,   recursive = TRUE, showWarnings = FALSE)

source(THEME_R)       # theme_pub(), PAL_CELLTYPE, PAL_REGION, REGION_FULL
source(ARTIFACT_R)    # is_artifact(gene, cell_type)

CELLTYPES <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
REGIONS   <- c("Frontal","Hippo")

# Panel d point categories (Frontal x Hippo discovery-DEG membership):
#   Shared NHD-up / Shared CON-up  (discovery DEG in both, sign-concordant)
#   Discordant                     (discovery DEG in both, opposite sign)
#   Frontal-only                   (discovery DEG in Frontal, not Hippo)
#   Hippo-only                     (discovery DEG in Hippo, not Frontal)
# Frontal-only / Hippo-only use the region UMAP tones so the reader maps them to
# the region colours across panels b and d. Shared/Discordant keep the CB-safe
# red/blue/grey scheme from the 3-region reference.
CATCOL <- c(`Shared NHD-up` = unname(PAL_COND[["NHD"]]), `Shared CON-up` = unname(PAL_COND[["CON"]]),
            `Discordant`         = "grey78",
            `Frontal-only`       = unname(PAL_REGION[["Frontal"]]),
            `Hippocampus-only`   = unname(PAL_REGION[["Hippo"]]))   # full name (house rule); banner fill still keyed by Hippo code above
CATLEV <- names(CATCOL)

# ============================================================================
# Load + ambient-filter (identical policy to script 43 / 41)
# ============================================================================
cat("== load ==\n")
dat <- read.csv(DATA_CSV, stringsAsFactors = FALSE)
stopifnot(all(c("cell_type","region","gene","avg_log2FC","discovery") %in% colnames(dat)))
n0 <- nrow(dat)
dat <- dat %>%
  filter(cell_type %in% CELLTYPES, region %in% REGIONS, !is.na(avg_log2FC)) %>%
  filter(!is_artifact(gene, cell_type))     # drop ambient/artifact rows entirely
cat(sprintf("rows: %d -> %d after celltype/region/ambient filter\n", n0, nrow(dat)))
stopifnot(nrow(dat) > 0)

# Accessors keyed by (cell_type, region). discovery is logical here but normalise
# defensively in case the CSV round-trips it as character.
isTRUE_vec <- function(x) (x %in% c(TRUE, "TRUE", "True", "true"))
.lfc_cache  <- list()
.disc_cache <- list()
.minpct_cache <- list()   # gene -> min(pct.1, pct.2) per (ct, region) — for callout floor
stopifnot(all(c("pct.1","pct.2") %in% colnames(dat)))
for (ct in CELLTYPES) for (rg in REGIONS) {
  sub <- dat[dat$cell_type == ct & dat$region == rg, , drop = FALSE]
  key <- paste(ct, rg, sep = "||")
  v <- sub$avg_log2FC; names(v) <- sub$gene
  .lfc_cache[[key]]  <- v
  .disc_cache[[key]] <- sub$gene[isTRUE_vec(sub$discovery)]
  mp <- pmin(sub$pct.1, sub$pct.2); names(mp) <- sub$gene
  .minpct_cache[[key]] <- mp
}
lfc_vec  <- function(ct, rg) .lfc_cache[[paste(ct, rg, sep = "||")]]
disc_set <- function(ct, rg) .disc_cache[[paste(ct, rg, sep = "||")]]
minpct_vec <- function(ct, rg) .minpct_cache[[paste(ct, rg, sep = "||")]]

# ============================================================================
# Frontal-vs-Hippo cross-region concordance  (single region pair)
# ----------------------------------------------------------------------------
# gene set = union of discovery DEGs in {Frontal U Hippo}, restricted to genes
#            with non-NA avg_log2FC in both regions (tested/present in both).
# r = Pearson(avg_log2FC_Frontal, avg_log2FC_Hippo) over those genes.
# Guard: < MIN_SHARED shared genes -> r = NA (no crash), logged as skipped.
# The gene-set choice matches script 43: union of the two discovery sets,
# intersected with genes tested in both regions.
# ============================================================================
MIN_SHARED <- 15L   # an r on <15 shared genes is unstable; shouldn't show to 2 dp
cat("== Frontal-vs-Hippo cross-region concordance ==\n")

series_rows <- list()   # per-gene points for panel d (x=Frontal, y=Hippo)
r_rows      <- list()   # per-cell-type r + n for panels d/e + CSV
skipped     <- character(0)

for (ct in CELLTYPES) {
  vf <- lfc_vec(ct, "Frontal"); vh <- lfc_vec(ct, "Hippo")
  union_genes <- union(disc_set(ct, "Frontal"), disc_set(ct, "Hippo"))
  shared <- intersect(union_genes, intersect(names(vf), names(vh)))
  n <- length(shared)
  if (n >= MIN_SHARED) {
    r <- suppressWarnings(stats::cor(vf[shared], vh[shared], method = "pearson"))
  } else {
    r <- NA_real_
    skipped <- c(skipped, sprintf("%s (n=%d < %d)", ct, n, MIN_SHARED))
  }
  r_rows[[length(r_rows) + 1L]] <- data.frame(
    cell_type = ct, r = r, n_genes = n, stringsAsFactors = FALSE)
  if (n >= MIN_SHARED) {
    mpf <- minpct_vec(ct, "Frontal"); mph <- minpct_vec(ct, "Hippo")
    series_rows[[length(series_rows) + 1L]] <- data.frame(
      cell_type = ct, gene = shared,
      x = as.numeric(vf[shared]), y = as.numeric(vh[shared]),
      in_fr = shared %in% disc_set(ct, "Frontal"),
      in_hp = shared %in% disc_set(ct, "Hippo"),
      minpct_fr = as.numeric(mpf[shared]), minpct_hp = as.numeric(mph[shared]),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("  %-11s n=%4d  r=%s\n", ct, n,
              ifelse(is.na(r), "NA", sprintf("%.3f", r))))
}
series_df <- bind_rows(series_rows)
r_df      <- bind_rows(r_rows)

if (length(skipped)) {
  cat("== SKIPPED (r set to NA, <", MIN_SHARED, "shared genes) ==\n", sep = "")
  for (s in skipped) cat("   -", s, "\n")
} else {
  cat("== no cell type skipped — all >= ", MIN_SHARED, " shared genes ==\n", sep = "")
}

# ============================================================================
# Panel f data: Frontal-unique / shared / Hippo-unique discovery features
# ----------------------------------------------------------------------------
# Over the discovery==TRUE (ambient-filtered) sets per region:
#   Shared        = gene is a discovery DEG in both Frontal and Hippo
#   Frontal-only  = discovery DEG in Frontal only
#   Hippo-only    = discovery DEG in Hippo only
# These are set-membership counts (not restricted to genes tested in both), so
# they answer "where is the NHD response shared vs region-restricted".
# ============================================================================
cat("== shared vs unique discovery features (panel f) ==\n")
share_rows <- list()
for (ct in CELLTYPES) {
  d_fr <- disc_set(ct, "Frontal"); d_hp <- disc_set(ct, "Hippo")
  shared_g <- intersect(d_fr, d_hp)
  fr_only  <- setdiff(d_fr, d_hp)
  hp_only  <- setdiff(d_hp, d_fr)
  share_rows[[length(share_rows) + 1L]] <- data.frame(
    cell_type = ct,
    n_frontal_only = length(fr_only),
    n_shared       = length(shared_g),
    n_hippo_only   = length(hp_only),
    n_frontal_disc = length(d_fr),
    n_hippo_disc   = length(d_hp),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-11s Frontal-only=%4d  shared=%4d  Hippo-only=%4d\n",
              ct, length(fr_only), length(shared_g), length(hp_only)))
}
share_df <- bind_rows(share_rows)

# ============================================================================
# Assemble report CSV (source of truth for all d/e/f legend numbers)
# ============================================================================
report <- r_df %>%
  left_join(share_df, by = "cell_type") %>%
  mutate(cell_type = factor(cell_type, levels = CELLTYPES)) %>%
  arrange(cell_type) %>%
  transmute(cell_type,
            r_Frontal_Hippo = r, n_shared_tested = n_genes,
            n_frontal_only, n_shared, n_hippo_only,
            n_frontal_disc, n_hippo_disc)
write.csv(report, OUT_CSV, row.names = FALSE)
cat("\n== REPORT ==\n")
print(report, row.names = FALSE, digits = 3)
cat("wrote:", OUT_CSV, "\n\n")

# ============================================================================
# Panel d — F1f_concordance_grid_MAST  (6 facets: Frontal x vs Hippo y)
# ============================================================================
cat("== render panel d ==\n")
series_df$cell_type <- factor(series_df$cell_type, levels = CELLTYPES)
# category per point from discovery membership + sign concordance:
series_df$cat5 <- factor(with(series_df,
  ifelse(in_fr & in_hp,
         ifelse(x > 0 & y > 0, "Shared NHD-up",
         ifelse(x < 0 & y < 0, "Shared CON-up", "Discordant")),
  ifelse(in_fr & !in_hp, "Frontal-only", "Hippocampus-only"))),
  levels = CATLEV)

# per-facet r annotation (bottom-right, empty in these scatters) — from r_df, with
# the shared-gene n so an r never appears without the count it rests on.
ann_d <- r_df %>%
  mutate(cell_type = factor(cell_type, levels = CELLTYPES),
         label = ifelse(is.na(r), "r = NA",
                        sprintf("r = %.2f (n=%d)", r, n_genes)))

# per-facet symmetric equal x/y ranges via invisible corner points so the dashed
# y=x diagonal lands at a true 45 deg under aspect.ratio=1 (ggplot2 4.0 forbids
# coord_fixed + free scales).
lims_d <- series_df %>%
  group_by(cell_type) %>%
  summarise(m = max(abs(c(x, y))) * 1.15, .groups = "drop")
blank_d <- bind_rows(transform(lims_d, x = m, y = m),
                     transform(lims_d, x = -m, y = -m))
blank_d$cell_type <- factor(blank_d$cell_type, levels = CELLTYPES)

# ----------------------------------------------------------------------------
# CALLOUT curation. The scatter
# points are untouched; only the labelled/ringed callouts are curated:
#   (4) expression floor: label-eligible only if min(pct.1,pct.2) >= PCT_FLOOR in
#       both regions -> excludes near-binary on/off dropout genes (huge log2FC
#       because pct~0 in one condition), which read as artifacts. e.g. the whole
#       Oligo APC/TANC2/CALN1/SACS/... set (pct.1~0-0.02) falls out here; we do
#       not force a per-cell oligo callout (oligo loss is compositional) — the
#       floor simply keeps whatever graded genes remain.
#   (5) extended non-coding / cross-lineage exclusion via is_callout_excluded().
#   (6) prefer-inject canonical headline genes when they are in the concordant
#       set and clear the floor (MT2A/MT3 astro; FTH1 micro; GPR17 OPC; etc.).
# ----------------------------------------------------------------------------
PCT_FLOOR   <- 0.05
N_LABEL     <- 8L
LAB_EXCLUDE <- c("CCDC26","F13A1")   # legacy ambient-adjacent (kept for parity)

# Canonical prefer-inject callouts per cell type (only used if present in the
# concordant + floor-passing eligible pool; never fabricated).
PREFER <- list(
  `Astro`      = c("MT2A","MT3"),
  `Micro-PVM`  = c("FTH1","HAMP","SLC2A3","SGK1","DUSP1","BCL6"),
  `OPC`        = c("GPR17"),
  `Neuron_Ex`  = c("SNAP25","SYT1"))

# Common gate for all callouts: shared sign-concordant, not an ambient-adjacent /
# extended-artifact / non-coding / cross-lineage symbol. (GLUL in Micro-PVM etc.
# Are removed here, so an excluded gene can never be rescued by prefer-inject.)
gate <- series_df %>%
  filter(cat5 %in% c("Shared NHD-up","Shared CON-up"),
         !gene %in% LAB_EXCLUDE) %>%
  filter(!is_callout_excluded(gene, cell_type))

# top-up pool: additionally clears the pct floor in both regions -> drops the
# near-binary on/off dropout genes (Oligo APC/CALN1/... pct~0). Auto-selected.
elig <- gate %>% filter(minpct_fr >= PCT_FLOOR, minpct_hp >= PCT_FLOOR)

# Prefer pool: the hand-vetted canonical callouts are exempt from the pct floor
# (but not from `gate`). Rationale: microglial disease genes
# HAMP/SLC2A3/SGK1/DUSP1/BCL6 are genuinely induced in NHD (pct.1 ~0.2-0.5) from a
# near-off CON baseline (pct.2 ~0.003-0.02); that asymmetry is the biology, not an
# artifact, so min(pct) would wrongly kill them. They stay only if `gate` passes.
pref_all <- unlist(PREFER)
prefer_pool <- gate %>%
  rowwise() %>% filter(gene %in% PREFER[[as.character(cell_type)]]) %>% ungroup()

# per facet: prefer-inject the canonical genes present in `prefer_pool`, then top
# up to N_LABEL by |x|+|y| from `elig` (floor-passing), avoiding duplicates.
pick_labels <- function(ct) {
  keep_pref <- prefer_pool %>% filter(as.character(cell_type) == ct)
  rest <- elig %>% filter(as.character(cell_type) == ct, !gene %in% keep_pref$gene) %>%
    slice_max(order_by = abs(x) + abs(y),
              n = max(0L, N_LABEL - nrow(keep_pref)), with_ties = FALSE)
  bind_rows(keep_pref, rest)
}
lab_d <- bind_rows(lapply(as.character(CELLTYPES), pick_labels))

# audit: report per-facet label counts + which prefer-injected genes landed / missed
.lc <- lab_d %>% count(cell_type, name = "n_labels")
cat("   panel d labels per facet (after pct-floor + exclusion + prefer-inject):\n")
for (i in seq_len(nrow(.lc))) {
  ct <- as.character(.lc$cell_type[i])
  labs_ct <- lab_d$gene[lab_d$cell_type == ct]
  # a prefer gene "missed" only if it failed the GATE (not merely the floor —
  # prefer genes are floor-exempt), so the log flags real biology drops only.
  pref_miss <- setdiff(PREFER[[ct]], gate$gene[gate$cell_type == ct])
  cat(sprintf("     %-11s %d  [%s]%s\n", ct, .lc$n_labels[i],
              paste(labs_ct, collapse = ", "),
              if (length(pref_miss)) sprintf("  prefer-miss(gate-fail): %s",
                                             paste(pref_miss, collapse = ",")) else ""))
}

# draw order: Discordant faint underneath, shared categories larger on top.
DRAW_ORDER <- c("Discordant","Hippocampus-only","Frontal-only","Shared CON-up","Shared NHD-up")
series_ord <- series_df %>%
  mutate(.ord = factor(as.character(cat5), levels = DRAW_ORDER)) %>%
  arrange(.ord)
disc_pts <- series_ord %>% filter(cat5 == "Discordant")
keep_pts <- series_ord %>% filter(cat5 != "Discordant")
# The journal pass: "Shared NHD-up" (house NHD red) and "Hippocampus-only"
# (house Hippocampus green) sit on the same scatter, which is the red-green contrast Nature
# Portfolio asks authors to avoid. The region tokens must stay (they key panels b/d), so a
# second, colour-free cue separates the two families: shared categories are filled discs,
# region-only categories are open rings. The legend carries both cues.
CATSHAPE <- c(`Shared NHD-up` = 16, `Shared CON-up` = 16, `Discordant` = 16,
              `Frontal-only` = 1, `Hippocampus-only` = 1)
shared_pts <- keep_pts %>% filter(cat5 %in% c("Shared NHD-up","Shared CON-up"))
region_pts <- keep_pts %>% filter(cat5 %in% c("Frontal-only","Hippocampus-only"))
# ring the labelled shared-concordant genes only (tie the highlight to the names).
lab_key  <- paste(lab_d$cell_type, lab_d$gene)
ring_pts <- keep_pts %>%
  filter(cat5 %in% c("Shared NHD-up","Shared CON-up"),
         paste(cell_type, gene) %in% lab_key)

# The single-NHD-Hippocampus-lane caveat lives in the figure legend, never on a
# panel. The old LOWN_CT / LOWCONF_MARK plumbing (panel-d facet strip + panel-f
# y-tick) has been deleted so the string can never resurface even if a constant
# is restored elsewhere. Facet strips now use the CT_DISPLAY recode (no "(low n)").
# ct_disp maps the coded facet key Neuron_Ex -> "Neuron Ex" etc. (house rule: no
# code tokens in labels); non-neuron types pass through unchanged.
ct_lowconf_labeller <- CT_LABELLER

# SE-corner obstacle for ggrepel (visual craft): the r/n annotation is
# pinned to the bottom-right of every facet (geom_text at Inf,-Inf). A near-origin
# callout occasionally gets flung into that same corner, overprinting "r=.. (n=..)"
# with a full-height leader line (Neuron_Ex). We add one invisible, blank-label
# obstacle point per facet at the SE corner (x = +max, y = -max, from lims_d) and
# feed it to the same geom_text_repel call: repel treats it as an occupied box and
# steers real gene labels away, so the corner stays clear for the r/n text. It is
# invisible (label=""), affects layout only, never prints, and is gene-agnostic
# (fixes the corner for all facets, not just whichever gene lands there today).
repel_obstacle <- dplyr::bind_rows(
  lims_d %>% transmute(cell_type, x =  m * 0.86, y = -m * 0.86, gene = ""),  # SE corner
  lims_d %>% transmute(cell_type, x =  m * 0.45, y = -m * 0.90, gene = ""),  # lower edge, right of centre
  lims_d %>% transmute(cell_type, x =  m * 0.90, y = -m * 0.45, gene = "")   # right edge, below centre
) %>% mutate(cell_type = factor(cell_type, levels = CELLTYPES))
lab_d_repel <- dplyr::bind_rows(lab_d, repel_obstacle)

p_d <- ggplot(series_df, aes(x = x, y = y, color = cat5)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, color = "grey45") +
  geom_blank(data = blank_d, aes(x = x, y = y), inherit.aes = FALSE) +
  geom_point(data = disc_pts, aes(shape = cat5), size = 0.50, alpha = 0.32, stroke = 0) +
  geom_point(data = region_pts, aes(shape = cat5), size = 1.15, alpha = 0.85, stroke = 0.38) +
  geom_point(data = shared_pts, aes(shape = cat5), size = 1.25, alpha = 0.80, stroke = 0) +
  geom_point(data = ring_pts, aes(x = x, y = y), inherit.aes = FALSE,
             shape = 21, fill = NA, colour = "black", size = 2.2, stroke = 0.5) +
  ggrepel::geom_text_repel(data = lab_d_repel, aes(x = x, y = y, label = gene),
            inherit.aes = FALSE, size = 1.9, fontface = "italic", color = "black",
            bg.color = "white", bg.r = 0.12, segment.color = "grey55",
            # box.padding/force bumped again (0.40->0.55, 3.5->6.0) + faint segments
            # (segment.alpha 0.6) + more point.padding so the two DENSEST facets fan
            # their leader lines out instead of crossing through the central cloud:
            # OPC (GPR17/FTH1/RASSF4/SIRT2/DNAJC6) and Neuron_Ex (YWHAG/YWHAB/SYT1/
            # VAMP2/CALM1 near the origin). Clean Micro-PVM/Oligo facets unaffected.
            segment.size = 0.16, segment.alpha = 0.6, min.segment.length = 0,
            box.padding = 0.55, point.padding = 0.28, force = 6.0, force_pull = 0.6,
            max.overlaps = Inf, max.time = 3, max.iter = 200000, seed = 42) +
  # r/n annotation, pinned to the bottom-right corner. The repel_obstacle above
  # keeps gene callouts out of this corner so the string is never overprinted.
  geom_text(data = ann_d, aes(x = Inf, y = -Inf, label = label),
            color = "black", hjust = 1.06, vjust = -0.7, size = 1.9,
            inherit.aes = FALSE, fontface = "plain") +
  facet_wrap(~ cell_type, scales = "free", nrow = 2, labeller = ct_lowconf_labeller) +
  scale_color_manual(values = CATCOL, name = NULL, breaks = names(CATCOL)) +
  scale_shape_manual(values = CATSHAPE, name = NULL, breaks = names(CATCOL)) +
  labs(x = "log2FC Frontal (MAST)", y = "log2FC Hippocampus (MAST)") +
  guides(color = guide_legend(override.aes = list(size = 2.0, alpha = 1, stroke = 0.6)),
         shape = guide_legend(override.aes = list(size = 2.0, alpha = 1, stroke = 0.6))) +
  theme_pub(base_size = 7) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        aspect.ratio = 1,
        panel.spacing = unit(0.6, "lines"),
        text        = element_text(colour = "black"),
        axis.text   = element_text(colour = "black"),
        axis.title  = element_text(colour = "black"),
        strip.text  = element_text(colour = "black"),
        legend.text = element_text(colour = "black"))

# 10% smaller (6.0 x 4.72 -> 5.40 x 4.25).
ggsave(file.path(PANEL_DIR, "F1f_concordance_grid_MAST.png"), p_d,
       width = 5.40, height = 4.25, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL_DIR, "F1f_concordance_grid_MAST.pdf"), p_d,
       width = 5.40, height = 4.25, bg = "white")
cat("  wrote F1f_concordance_grid_MAST.{png,pdf}\n")

# ============================================================================
# Panel e — F1e_concordance_bar_MAST  (one bar/cell type = Frontal-Hippo r)
# ----------------------------------------------------------------------------
# DEMOTED to supplement: it only re-plots the r values
# already printed in each panel-d facet. Rendered into panels/_supp/. Main Fig 1
# = a/b/c/d/f. No replacement burden panel (burden+direction live in c and f).
# ============================================================================
cat("== render panel e (SUPPLEMENT) ==\n")
# Remove any stale 2E_* left in the main panels dir from a prior run so it does
# not sit in the a-f set (idempotent demotion).
for (ext in c("png","pdf")) {
  stale_e <- file.path(PANEL_DIR, paste0("F1e_concordance_bar_MAST.", ext))
  if (file.exists(stale_e)) { file.remove(stale_e); cat("   removed stale main-dir", basename(stale_e), "\n") }
}
# rank ascending so coord_flip puts largest r at TOP; NA r sorts to bottom.
ord_e <- r_df %>% arrange(is.na(r), r) %>% pull(cell_type)
bar_e <- r_df %>% mutate(cell_type = factor(cell_type, levels = ord_e))

# shared bar style with panel f (matched sibling pair, per the reference).
BAR_W   <- 0.72
BAR_EXP <- expansion(mult = c(0.02, 0.20))
BAR_BRK <- c(0, 0.25, 0.5, 0.75, 1.0)
BAR_MAR <- margin(4, 12, 4, 4)
p_e <- ggplot(bar_e, aes(x = cell_type, y = r, fill = cell_type)) +
  geom_col(color = NA, width = BAR_W) +
  geom_text(aes(y = r + 0.02, label = ifelse(is.na(r), "NA", sprintf("%.2f", r))),
            hjust = 0, size = 2.0, color = "black") +
  scale_fill_manual(values = PAL_CELLTYPE, guide = "none") +
  # This axis printed the raw factor codes "Neuron_Ex"/"Neuron_Inh" while
  # every sibling panel (1C, 2D, 2F) uses the CT_DISPLAY recode.  House rule: no code
  # tokens with underscores on a printed axis.
  scale_x_discrete(labels = function(x) ct_disp(x)) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, 1), breaks = BAR_BRK, expand = BAR_EXP) +
  labs(x = NULL, y = "Frontal-Hippocampus concordance (Pearson r, MAST)") +
  theme_pub(base_size = 7) +
  theme(plot.margin = BAR_MAR)

ggsave(file.path(SUPP_DIR, "F1e_concordance_bar_MAST.png"), p_e,
       width = 3.2, height = 2.0, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(SUPP_DIR, "F1e_concordance_bar_MAST.pdf"), p_e,
       width = 3.2, height = 2.0, bg = "white")
cat("  wrote _supp/F1e_concordance_bar_MAST.{png,pdf}\n")

# ============================================================================
# Panel f — F1d_shared_vs_unique_MAST  (stacked bar: Frontal-only / shared / Hippo-only)
# ----------------------------------------------------------------------------
# Replaces the degenerate 2-region specificity index. Per cell type, a horizontal
# stacked bar of discovery-feature counts split into Frontal-unique, shared, and
# Hippo-unique. Ordered by total discovery burden (largest at top). Frontal-only
# and Hippo-only use the region UMAP tones; shared uses a neutral dark grey.
# ============================================================================
cat("== render panel f ==\n")
# Any direction) — deliberately not "Shared", which in panel d means SIGN-
# concordant-in-both. Same figure, so the words must not collide.
MID_LAB <- "Detected in both"
share_long <- share_df %>%
  mutate(total = n_frontal_only + n_shared + n_hippo_only) %>%
  arrange(total) %>%                                   # ascending -> largest at TOP after coord_flip
  mutate(cell_type_raw = as.character(cell_type),
         cell_type = factor(cell_type, levels = cell_type)) %>%
  select(cell_type, `Frontal-only` = n_frontal_only, !!MID_LAB := n_shared,
         `Hippocampus-only` = n_hippo_only) %>%
  pivot_longer(-cell_type, names_to = "cat", values_to = "n")
SEG_LEV <- c("Frontal-only", MID_LAB, "Hippocampus-only")   # stack order (bottom->top)
share_long$cat <- factor(share_long$cat, levels = SEG_LEV)
SEG_COL <- setNames(c(unname(PAL_REGION_UMAP[["Frontal"]]), "grey45",
                      unname(PAL_REGION_UMAP[["Hippo"]])), SEG_LEV)
# total label at the bar tip (keyable to n_frontal_only+n_shared+n_hippo_only).
tot_f <- share_long %>% group_by(cell_type) %>%
  summarise(total = sum(n), .groups = "drop")
xmax_f <- max(tot_f$total) * 1.12

# The single-NHD-
# Hippocampus-lane caveat is a legend note, never on the panel. y-ticks use the
f_ct_labs <- setNames(ct_disp(levels(share_long$cell_type)), levels(share_long$cell_type))

BAR_W   <- 0.72                                    # shared with panel e
BAR_MAR <- margin(4, 12, 4, 4)
p_f <- ggplot(share_long, aes(x = cell_type, y = n, fill = cat)) +
  geom_col(width = BAR_W, color = NA) +
  geom_text(data = tot_f, aes(x = cell_type, y = total + xmax_f * 0.01, label = total),
            hjust = 0, size = 2.0, color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = SEG_COL, name = NULL) +
  scale_x_discrete(labels = f_ct_labs) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(0, xmax_f),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Discovery features (per-nucleus MAST)") +
  theme_pub(base_size = 7) +
  theme(plot.margin = BAR_MAR,
        legend.position = "bottom",
        legend.key.size = unit(0.26, "cm"),
        # As part of matching sibling 2E's apparent font size; the key is
        # trimmed slightly so the 3-item legend still fits the narrower 3.2in canvas.
        legend.text = element_text(size = 6.8),
        legend.margin = margin(t = -2)) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE))

# Match sibling 2E's aesthetics — same canvas (3.2 x 2.0) and
# therefore the same apparent font size on the page. It had been widened to 3.9in
# to stop the 3-item legend clipping "Frontal-only"; the trimmed legend key above
# buys that room back at the sibling-matched width instead.
# [[apparent-font-parity-precompensate-wide-panels]]: a wider canvas at the same
# base_size renders smaller text once both panels are placed at one column width.
ggsave(file.path(PANEL_DIR, "F1d_shared_vs_unique_MAST.png"), p_f,
       width = 3.2, height = 2.0, dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(PANEL_DIR, "F1d_shared_vs_unique_MAST.pdf"), p_f,
       width = 3.2, height = 2.0, bg = "white")
cat("  wrote F1d_shared_vs_unique_MAST.{png,pdf}\n")

# ============================================================================
# Provenance
# ============================================================================
si_path <- file.path(LOG_DIR, "fig1_def_MAST_2region_FH_sessionInfo.txt")
writeLines(c(sprintf("run: %s", format(Sys.time())),
             sprintf("PROJ: %s", PROJ),
             capture.output(sessionInfo())), si_path)
cat("wrote provenance:", si_path, "\n")

cat("\n=== DONE: 43_fig1_def_MAST_2region_FH.R ===\n")
