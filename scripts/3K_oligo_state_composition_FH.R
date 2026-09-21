#!/usr/bin/env Rscript
# =============================================================================
# 3K_oligo_state_composition_FH.R — Figure 3k: oligo state-proportion shift, CON vs NHD, per region (Frontal, Hippocampus-grey-lowN). frontal + hippo.
# -----------------------------------------------------------------------------
# The load-bearing compositional-loss quantification: myelin loss in NHD is
# compositional — the oligo population redistributes away from the mature
# continuum toward the stress / high-myelin poles — not a per-cell gene collapse.
# One stacked bar per (Region x Condition) makes the redistribution undeniable:
#   Frontal   Mature   67.6% -> 21.0%   (the mature continuum collapses)
#             High-myelin 2.9% -> 32.0% (state expands)
#             Stress    0.2% -> 34.3%   (state expands)
#
# House rules honored:
#   * Effect size = the fraction shift itself. no significance stars (single-donor
#     composition; per-house-policy snRNA composition uses effect size, never
#     single-donor Fisher/pseudorep stars). The CON->NHD delta is annotated for the
#     two NHD states as the readable effect.
#   * Single-lane-Hippocampus honesty -> grey Hippo banner + "(low n)" (Hippo NHD =
#     one lane), same figure-wide treatment as every other astro/oligo Fig-3 panel.
#   * Region strips = pale banners + black text + full names; canonical Frontal ->
#     Hippocampus order. no bold, no in-panel captions.
#   * State palette = the shared PAL_STATE from 3F (same colours as 3h/3i).
#   * Legend numbers key to tables/oligo_state_composition_FH.csv (source of truth).
# Output: figures/Figure_3/panels/F4h_oligo_state_composition.{png,pdf}
# -----------------------------------------------------------------------------
# apparent-type lift 3K_oligo_state_composition_FH.R: every text size in this script scaled
# by 1.41 so the panel reads at ~5.0 pt on the assembled page, matching the
# Figure-5 dot panels. Point sizes, line widths and unit() dimensions are not
# touched. Canvas size unchanged, so re-linking is a no-op.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ggh4x)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # PAL_REGION*, REGION_FULL/ORDER, strip_region_x, region_full_lowconf_labeller

PANEL <- file.path(PROJ, "figures", "Figure_4", "panels")   # oligo -> Figure 4
TDIR  <- file.path(PROJ, "tables")
LOGD  <- file.path(PROJ, "logs")
for (d in c(PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

CSV <- file.path(TDIR, "oligo_state_composition_FH.csv")
if (!file.exists(CSV))
  stop("MISSING ", CSV, " — run 3F_oligo_embed_states_FH.R first (it writes this table).")

# ---- state palette (shared with 3F/3h/3i so colours match across the figure) --
PAL_STATE <- c(
  "OPC-proximal (early)"    = "#CCBB44",
  "Mature (continuum)"      = "#4477AA",
  "High-myelin (NHD)"       = "#228833",
  "Stress (NHD)"            = "#AA3377")
STATE_ORDER <- names(PAL_STATE)
REGIONS <- REGION_ORDER                              # c("Frontal","Hippo")

# ── single-lane-Hippocampus honesty (figure-wide) ────────────────────────────
HIPPO_NHD_LANES <- 1L                                # TFHS000535 is the only NHD Hippo lane
low_conf <- if (HIPPO_NHD_LANES < 2) "Hippo" else character(0)

# ---- load the composition table (single source of truth for the numbers) -----
comp <- read.csv(CSV, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("composition CSV missing expected columns" =
            all(c("Region","Condition","oligo_state","n","frac") %in% names(comp)))
stopifnot("unexpected states in CSV" = setequal(unique(comp$oligo_state), STATE_ORDER),
          "unexpected regions in CSV" = setequal(unique(comp$Region), REGIONS))

# validate the fractions sum to 1 within each (Region x Condition) group (fail loud
# if the CSV is malformed, rather than silently plotting a non-normalized bar).
grp_sum <- comp %>% group_by(Region, Condition) %>%
  summarise(s = sum(frac), .groups = "drop")
if (any(abs(grp_sum$s - 1) > 1e-6))
  stop("fractions do not sum to 1 within a Region x Condition group:\n",
       paste(capture.output(print(as.data.frame(grp_sum))), collapse = "\n"))
cat("== fraction sums per (Region x Condition) all ~1.0 (validated) ==\n")

comp <- comp %>% mutate(
  oligo_state = factor(oligo_state, levels = STATE_ORDER),
  Region      = factor(Region,      levels = REGIONS),
  Condition   = factor(Condition,   levels = c("CON","NHD")),
  pct         = 100 * frac)

# ---- CON->NHD fraction shift = the effect SIZE (written to CSV, not the panel) --
# No significance test: single-donor composition per house policy. The shift is the
# effect. The percentage-point shift is computed here so it is reported to the log
# and (already) carried in the source CSV, but is not printed inside the bars.
wide <- comp %>%
  select(Region, Condition, oligo_state, pct) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = pct) %>%
  mutate(d_pp = NHD - CON)                            # percentage-point shift
cat("\n== CON->NHD state shift (percentage points) ==\n")
print(as.data.frame(wide %>% mutate(across(where(is.numeric), ~round(.x, 1)))), row.names = FALSE)

# At 11.3 pt "Hippocampus" was wider than its facet strip and overflowed
# the pale banner on both sides (clip = "off" lets it spill rather than truncate).
# The panel is 4.0 in wide with two facets plus a y axis, so the strip is narrow;
# the fix is the strip type, not the region name, which stays spelled out per house
# style. 9.4 pt is the largest size at which "Hippocampus" sits inside its banner here.
STRIP_PT <- 9.4   # single source for both the helper and the theme pin below
strip_lo <- strip_region_x(REGIONS, low_conf = low_conf, fontsize = STRIP_PT, clip = "off")

p <- ggplot(comp, aes(x = Condition, y = pct, fill = oligo_state)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.2,
           position = position_stack(reverse = TRUE)) +
  facet_wrap2(~ Region, nrow = 1,
              labeller = region_full_lowconf_labeller(low_conf), strip = strip_lo) +
  scale_fill_manual(values = PAL_STATE, name = "Oligo state",
                    breaks = STATE_ORDER, limits = STATE_ORDER) +
  scale_y_continuous(name = "State proportion (%)", expand = expansion(mult = c(0, 0.02)),
                     breaks = seq(0, 100, 25)) +
  labs(x = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, reverse = FALSE)) +
  theme_pub(base_size = 11.3) +
  # This line used to set size = 11.3, which silently OVERRODE the
  # 9.4 pt passed to strip_region_x() six lines above -- so "Hippocampus" rendered at
  # 11.3 pt and overflowed its banner in every export, while the comment claimed 8 pt and
  # the line above claimed 9.4. The override dates from the era when elem_list_text()
  # styled only the first strip; strip_region_x now
  # sizes both banners itself, so the size is pinned to the same value here rather than a
  # different one. Keep these two numbers equal.
  theme(strip.text = element_text(size = STRIP_PT, face = "plain"),
        strip.clip = "off",
        axis.line = element_line(colour = "black", linewidth = 0.35),
        panel.spacing.x = unit(0.6, "lines"),
        # Legend moved right -> bottom. On the right it consumed roughly half
        # of the 4-inch panel, which squeezed each region banner down to 499 px of width
        # for 470 px of text -- a 2.6% margin that reads as overflow once the panel is
        # placed small on the page. At the bottom the facets take nearly the full width,
        # so "Hippocampus" clears its banner comfortably at the same 9.4 pt (no type is
        # made smaller) and the bars get wider too. Panel frame is unchanged, so the
        # placed artwork relinks 1:1.
        legend.position = "bottom", legend.key.size = unit(0.30, "cm"),
        legend.text = element_text(size = 8.6), legend.title = element_text(size = 9.0),
        legend.margin = margin(1, 0, 0, 0), legend.box.spacing = unit(0.10, "cm"),
        plot.margin = margin(5, 4, 3, 4))

BN <- "F4h_oligo_state_composition"
W <- 4.0; H <- 2.03   # 25% shorter; width and type unchanged
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600, device = ragg::agg_png)
cat(sprintf("Saved %s.{pdf,png}  (%.1f x %.1f in)\n", BN, W, H))

writeLines(capture.output(sessionInfo()), file.path(LOGD, "3K_oligo_state_composition_FH_sessionInfo.txt"))
cat("=== DONE ===\n", file = stderr())
