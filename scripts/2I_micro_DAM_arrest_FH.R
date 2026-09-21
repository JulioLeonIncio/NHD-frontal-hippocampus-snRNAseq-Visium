# =============================================================================
# 2I_micro_DAM_arrest_FH.R — F2e_micro_DAM_arrest_FH.R — Figure 2 panel I: the DAM-1 -> DAM-2 arrest contrast. frontal + HIPPOCAMPUS (2 regions). Microglia (Micro-PVM).
# -----------------------------------------------------------------------------
# Thesis panel.  States Figure 2 in one glance: NHD DAP12-null microglia enter
# the TREM2-INDEPENDENT DAM-1 program (homeostatic P2RY12 exit + APOE/FTH1/B2M/
# CD74/HLA-DRA/C1QB induction) but arrest before the TREM2-DEPENDENT DAM-2 lipid
# endpoint (PPARG/ITGAX/GPNMB/CLEC7A/TREM2 flat; LPL/CST7/LGALS3 never detected).
#
# Two ARMS (the contrast is the panel):
#   arm 1  "DAM-1 entry (TREM2-independent)"      -> reads mostly ON  (coloured bars + delta dots)
#   arm 2  "DAM-2 endpoint (TREM2-dependent)"     -> reads mostly OFF (flat / grey x)
#
# ENCODING (reused verbatim from 2D_secretome_MAST_FH.R — house style):
#   * avg log2FC bar, coloured by direction (NHD-up warm red / CON-up cool blue).
#   * effect-size = filled dot at the bar, SIZE = |Cliff's delta|, only for the
#     strict/discovery tier (padj<0.05 and |delta|>=0.15).  no bare per-nucleus
#     p-stars (single NHD donor -> pseudoreplication; effect size + q* only).
#   * grey "q*" text = the house effect-gated mark (padj<0.05 and |delta|<0.15,
#     i.e. significant-but-small).
#   * grey open "x" at 0 = below DETECTION (gene absent from the region's MAST
#     background, or max(pct.1,pct.2) < 0.10).  A gene with a real n.s. test is a
#     small/flat coloured bar; a gene absent from the table is a grey x.  These
#     are not plotted as biological zero.
#   * both region banners use the normal pale region tone (Frontal mauve /
#     Hippocampus green); the single-NHD-hippo-lane caveat is a legend note,
#     never an on-plot greyed banner.
#
# Detection semantics:
#   arm 2 splits into three visible states, all distinct from a real DEG:
#     tested-but-flat  : in table, max_pct>=0.10, padj n.s.  -> flat coloured bar
#                        (ITGAX both; CLEC7A both; TREM2 both; GPNMB Hippo;
#                         PPARG Hippo -1.55 n.s.)
#     below-detection  : not in the region's MAST table or max_pct<0.10 -> grey x
#                        (PPARG Frontal; GPNMB Frontal; LPL/CST7/LGALS3 both)
#   The brief's "PPARG ~-3.71 tested-DOWN Frontal" is not what the table shows:
#   PPARG is absent from the Frontal Micro-PVM background -> drawn as grey x
#   (below detection), Hippo PPARG is the tested-flat -1.55 (padj 0.51).  The
#   table governs (fail-loud): we do not invent a Frontal PPARG bar.
#
# Outputs:
#   figures/Figure_2/panels/F2e_micro_DAM_arrest.{png,pdf}
#   tables/mast_dual/fig2_DAM_arrest_FH.csv   (legend-keyable source values)
#   logs/F2e_micro_DAM_arrest_FH.log           (tee'd at run time)
#   diagnostics/10_microglia_core_signature/2I_DAM_arrest_sessionInfo_FH.txt
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(ggh4x); library(ragg)
})
set.seed(42)

# --- Portable project root (resolve whichever machine root exists) -----------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))

THEME     <- file.path(PROJ, "scripts", "22_publication_theme_FH.R")
ARTIFACT  <- file.path(PROJ, "scripts", "_artifact_genes.R")
PANEL     <- file.path(PROJ, "figures", "Figure_2", "panels")
OUT_DIAG  <- file.path(PROJ, "diagnostics", "10_microglia_core_signature")
TAB_DIR   <- file.path(PROJ, "tables", "mast_dual")
LOGS      <- file.path(PROJ, "logs")
for (d in c(PANEL, OUT_DIAG, TAB_DIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

CT      <- "Micro-PVM"   # POOLED microglia+PVM; PVM = 5.4% of compartment
REGIONS <- c("Frontal", "Hippo")               # canonical FH order
MAST <- setNames(
  file.path(TAB_DIR, sprintf("MAST_Micro-PVM_%s.csv", REGIONS)), REGIONS)

for (f in c(THEME, ARTIFACT, MAST))
  if (!file.exists(f)) stop("MISSING ", f, " — upstream prep step did not run.")
source(THEME)       # theme_pub, PAL_REGION_PALE, REGION_FULL, REGION_ORDER
source(ARTIFACT)    # is_artifact()
stopifnot(identical(REGIONS, REGION_ORDER))    # theme + this script agree on FH order

PCT_FLOOR <- 0.10   # per-region detection floor (matches 2D_secretome_MAST_FH.R)

# --- The two arms (the panel's structural thesis) ----------------------------
# Gene order within each arm is deliberate:
#   Arm 1: the homeostatic EXIT (P2RY12, down) first, then the DAM-1 induction set.
#   Arm 2: the tested-but-flat DAM-2 genes first, the never-detected lipid
#          endpoint (LPL/CST7/LGALS3) last — so the eye reads "fizzles out".
# ASCII hyphen (not U+2014 em-dash): the base pdf() device cannot encode the
# em-dash (mbcsToSbcs conversion), and these panels save both .png (ragg) and
# .pdf (base) — so the strip label must be plain ASCII to stay byte-consistent.
# Line-wrapped (\n) so the left arm strips stay a narrow column and do not steal
# horizontal space from the plotting region (a single long horizontal label made
# ggh4x widen the strip column and squeeze the bars).
arm_levels <- c("DAM-1 entry\n(TREM2-independent)",
                "DAM-2 endpoint\n(TREM2-dependent)\nnot engaged")
dam_map <- tibble::tribble(
  ~gene,       ~arm,
  "P2RY12",    arm_levels[1],   # homeostatic exit (down)
  "APOE",      arm_levels[1],
  "FTH1",      arm_levels[1],
  "B2M",       arm_levels[1],
  "CD74",      arm_levels[1],
  "HLA-DRA",   arm_levels[1],
  "C1QB",      arm_levels[1],
  "PPARG",     arm_levels[2],
  "ITGAX",     arm_levels[2],
  "GPNMB",     arm_levels[2],
  "CLEC7A",    arm_levels[2],
  "TREM2",     arm_levels[2],
  "LPL",       arm_levels[2],   # DAM-2 lipid endpoint — never detected
  "CST7",      arm_levels[2],
  "LGALS3",    arm_levels[2])
stopifnot(!any(duplicated(dam_map$gene)))
gene_levels <- dam_map$gene       # top-to-bottom display order (reversed at plot)
cat(sprintf("== DAM arrest panel: %d genes / 2 arms ==\n", nrow(dam_map)))

# --- Load MAST per region, restrict to the DAM gene set ----------------------
load_region <- function(reg) {
  d <- read.csv(MAST[[reg]], stringsAsFactors = FALSE)
  need <- c("cell_type","region","gene","avg_log2FC","p_val_adj",
            "pct.1","pct.2","cliffs_delta","discovery")
  stopifnot(all(need %in% names(d)))
  d <- d %>% filter(cell_type == CT, region == reg)
  stopifnot(nrow(d) > 0)
  d
}
mast_list <- lapply(REGIONS, load_region); names(mast_list) <- REGIONS
mp <- bind_rows(mast_list)

cat("== Microglia MAST background sizes per region ==\n")
for (r in REGIONS)
  cat(sprintf("  %-8s : %d genes tested\n", r,
              length(unique(mp$gene[mp$region == r]))))

# None of the DAM genes are ambient artifacts, but apply the filter for parity
# with the other Micro-PVM panels and to fail loud if that ever changes.
art <- dam_map$gene[is_artifact(dam_map$gene, CT)]
if (length(art)) stop("Unexpected: DAM gene(s) flagged as artifact: ",
                      paste(art, collapse = ", "))
mp <- mp %>% filter(!is_artifact(gene, CT))

# --- Build the full gene x region grid (absent genes = below detection) -------
grid <- expand_grid(gene = gene_levels, region = REGIONS)
df <- grid %>%
  left_join(mp %>% select(gene, region, avg_log2FC, p_val_adj,
                          pct.1, pct.2, cliffs_delta, discovery),
            by = c("gene", "region")) %>%
  mutate(
    in_table = !is.na(avg_log2FC),
    max_pct  = pmax(pct.1, pct.2),
    min_pct  = pmin(pct.1, pct.2),
    # Statistically SUPPORTED = the test actually returned an effect. Genes with
    # NA Cliff's delta / padj >= 0.05 are not supported and must not look like a result.
    supported = in_table & !is.na(cliffs_delta) & !is.na(p_val_adj) & p_val_adj < 0.05,
    # DETECTED = in the region's MAST background and above the pct floor. A gene
    # absent from the table (max_pct NA) is below detection, not a zero.
    # The house guardrail (min(pct) >= 0.10 & non-NA delta to
    # colour an effect glyph) is now applied to unsupported genes as well, so an
    # unsupported gene that is essentially absent from one condition renders as
    # below-detection rather than as a bar. Supported genes keep the max_pct rule, so
    # genuine on/off biology (F13A1: 33.7% NHD vs 0.3% CON, padj 2.5e-49) is unaffected.
    detected = in_table & !is.na(max_pct) & max_pct >= PCT_FLOOR &
               (supported | (!is.na(min_pct) & min_pct >= PCT_FLOOR)),
    dir = case_when(!detected      ~ NA_character_,
                    avg_log2FC > 0 ~ "NHD-up",
                    TRUE           ~ "CON-up"),
    abs_cd  = abs(cliffs_delta),                       # NA-safe (some tested genes have NA delta)
    is_disc = detected & !is.na(discovery) & as.logical(discovery),  # padj<0.05 & |delta|>=0.15
    # Same effect-size encoding as the secretome/volcano panels (no bare p-stars):
    #   delta  -> big filled dot, size = |delta|  (strict/discovery tier)
    #   qstar  -> grey "q*" text                  (significant-but-small)
    #   none   -> nothing (tested-but-flat, or below detection)
    eff_tier = case_when(
      !detected | is.na(p_val_adj)   ~ "none",
      is_disc                        ~ "delta",
      p_val_adj < 5e-2 & !is_disc    ~ "qstar",
      TRUE                           ~ "none"),
    star = ifelse(eff_tier == "qstar", "q*", ""),
    # display status for the source CSV / log (3 states as briefed)
    status = case_when(
      !detected  ~ "below_detection",
      is_disc    ~ "significant",
      eff_tier == "qstar" ~ "significant",   # significant-but-small still = significant test
      TRUE       ~ "ns"),
    arm  = factor(dam_map$arm[match(gene, dam_map$gene)], levels = arm_levels),
    gene = factor(gene, levels = rev(gene_levels)),      # rev -> top gene at top of panel
    region = factor(region, levels = REGIONS))

# --- Source CSV (legend-keyable) ---------------------------------------------
src <- df %>%
  transmute(gene = as.character(gene),
            arm = as.character(arm),
            region = as.character(region),
            avg_log2FC = round(avg_log2FC, 3),
            cliffs_delta = round(cliffs_delta, 3),
            padj = p_val_adj,
            status) %>%
  arrange(match(arm, arm_levels), match(gene, gene_levels), match(region, REGIONS))
csv_path <- file.path(TAB_DIR, "fig2_DAM_arrest_FH.csv")
write.csv(src, csv_path, row.names = FALSE)
cat(sprintf("Wrote source CSV: %s (%d rows)\n", basename(csv_path), nrow(src)))

cat("\n== Per-gene status (Frontal | Hippocampus) ==\n")
print(as.data.frame(src), row.names = FALSE)

cat("\n== Per-arm x region tally ==\n")
tally <- df %>% count(arm, region, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0)
print(as.data.frame(tally))

# --- Render ------------------------------------------------------------------
L2FC_CAP <- 5
df_det <- df %>% filter(detected) %>%
  mutate(l2fc_cap = pmax(pmin(avg_log2FC, L2FC_CAP), -L2FC_CAP),
         x_tip  = l2fc_cap + ifelse(l2fc_cap > 0, 0.12, -0.12),
         hj_tip = ifelse(l2fc_cap > 0, 0, 1))
# The single "below detection" glyph conflated two different facts and
# stated one of them falsely. Genes absent from the region's MAST background (LPL, CST7,
# LGALS3) really were never tested. But TREM2, GPNMB, PPARG and frontal CLEC7A are tested
# and are detected -- GPNMB reaches 14.5% of NHD nuclei in hippocampus -- and are held
# back only by the paired-detection guardrail (below 10% in the other condition and not
# significant), which exists so a large fold change cannot rest on near-absence in one
# arm. Calling those "below detection" contradicted Supplementary Table 1, and
# contradicted Fig 2f, which prints GPNMB hippocampus +3.8 in the same figure because it
# applies a max-pct rule rather than this panel's min-pct rule. The two classes are now
# drawn and named separately.
df_nt  <- df %>% filter(!detected, !in_table)   # never tested in this region
df_nd  <- df %>% filter(!detected,  in_table)   # tested, but not paired-estimable
df_eff <- df_det %>% filter(eff_tier == "delta")   # strict effect-size tier

# Font : PANEL ratio, CLAMPED to the print FLOOR.
# F2f_secretome_MAST is 3.90 in wide at base_size 8 / axis.text 6.3 -> 2.051 / 1.615 pt per
# inch of width. Holding that ratio literally at this panel's new 2.396 in gives
# FS = 0.614, which puts axis.text at 3.75 pt and every token under 5 pt -- below the
# print floor for Acta Neuropathologica (and any journal). 2Y survived the same treatment
# only because it started at 4.60 in; 2I starts at 3.33 in, so the identical relative
# shrink lands it far smaller.
# Resolution: keep the REQUESTED geometry exactly, and scale fonts by the largest factor
# that keeps the ratio's hierarchy while putting the smallest text on the 5 pt floor:
#     FS = 5.0 / 6.1 = 0.820   (vs 0.614 unclamped)
# Ratio spacing between tokens is preserved, so the panel still reads as the same design
# family as 2E; only the absolute floor is enforced. 0.820 also lands within 4% of 2Y's
# 0.849, so 2Y and 2I end up visually consistent with each other too.
FS_RATIO <- (5.2 * 0.8 * 0.8 * 0.8 * 0.9) / 3.90   # 0.614 -- documented, not used
MIN_PT   <- 5.0
FS       <- MIN_PT / 6.1                            # 0.820 -- smallest text hits the floor
sz <- function(x) round(x * FS, 2)
cat(sprintf("font scale: ratio-parity would be %.3f (axis.text %.2f pt, BELOW %.1f pt floor); using %.3f (axis.text %.2f pt)\n",
            FS_RATIO, 6.1*FS_RATIO, MIN_PT, FS, 6.1*FS))

p <- ggplot() +
  geom_blank(data = df, aes(x = 0, y = gene)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  # four-state encoding. Previously a
  # non-significant bar was drawn identically to a significant one and the only cue
  # was the absence of an effect-size dot — undefined in the key. APOE (padj 0.219
  # Frontal / 1.0 Hippo) was rendering as a solid NHD-red bar, which reads as the
  # canonical DAM-1 result it no longer is. Now: solid fill = statistically supported
  # (padj < 0.05); hollow outline = tested but not significant; grey x = below
  # detection. All three are defined in the key.
  geom_col(data = df_det %>% filter(supported),
           aes(x = l2fc_cap, y = gene, fill = dir),
           width = 0.72, colour = "grey30", linewidth = 0.15) +
  geom_col(data = df_det %>% filter(!supported),
           aes(x = l2fc_cap, y = gene, colour = dir),
           fill = NA, width = 0.72, linewidth = 0.35, linetype = "22",
           show.legend = FALSE) +
  # The hollow layer must not raise its own "dir" colour legend (it duplicates the
  # fill key with a meaningless title); the dashed = not-significant convention is
  # defined in the figure legend text, alongside the x = below-detection mark.
  scale_colour_manual(values = c("NHD-up" = unname(PAL_COND[["NHD"]]), "CON-up" = unname(PAL_COND[["CON"]])),
                      guide = "none") +
  # numeric avg log2FC at every detected bar tip (no p-stars).
  geom_text(data = df_det,
            aes(x = x_tip, y = gene, hjust = hj_tip,
                label = sprintf("%.1f", avg_log2FC)),
            size = sz(1.7), colour = "grey15") +
  # effect-size dot inside the bar near 0, size = |Cliff's delta|, strict tier only.
  geom_point(data = df_eff,
             aes(x = ifelse(l2fc_cap > 0, 0.55, -0.55), y = gene,
                 size = abs_cd, fill = dir),
             shape = 21, stroke = 0.2, colour = "white",
             show.legend = c(fill = FALSE, size = TRUE)) +
  # house effect-gated "q*" mark for significant-but-small, fixed outer column.
  geom_text(data = df_det %>% filter(eff_tier == "qstar"),
            aes(x = ifelse(l2fc_cap > 0, 4.6, -4.6), y = gene, hjust = 0.5,
                label = star),
            size = sz(1.7), colour = "grey50") +
  # not tested in this region -> grey "x" at 0; tested but not paired-estimable ->
  # grey open circle at 0. Neither is a flat coloured bar.
  geom_point(data = df_nt, aes(x = 0, y = gene, shape = "not tested"),
             colour = "grey60", size = sz(1.7), stroke = 0.5) +
  geom_point(data = df_nd, aes(x = 0, y = gene, shape = "detected, not estimable"),
             colour = "grey60", size = sz(1.7), stroke = 0.5) +
  facet_grid2(arm ~ region, scales = "free_y", space = "free_y",
              switch = "y", labeller = labeller(region = as_labeller(REGION_FULL)),
              strip = strip_themed(
                # clip="off": the narrow Hippocampus facet would otherwise truncate
                # the full region name ("ippocampu"); off lets the centred banner
                # label overflow the pale strip instead of being cut.
                clip = "off",
                # Arm strip = neutral pale grey with black text (structural group
                # label, not a data colour) — arms are the panel's thesis dividers.
                background_y = elem_list_rect(fill = NA,
                                              colour = "grey80",
                                              linewidth = 0.25),
                text_y = elem_list_text(angle = 0, hjust = 0.5, size = sz(6.2),
                                        face = "plain", colour = "black"),
                # Region banners: both regions use their normal pale region tone.
                # Matches 2D/2F/2_zhou siblings.
                background_x = elem_list_rect(
                  fill = unname(PAL_REGION_PALE[REGIONS]), colour = NA),
                # apparent-font parity: this
                # panel renders at W=5.2in, between the wide (2D volcano W=7.4) and
                # narrow (2F W=4.6) siblings, so its size-7 banner reads slightly
                # large on the assembled page.  Pre-scale by render_W / W_REF
                # (5.2/5.4=0.963, W_REF=2_zhou's not-flagged 5.4in render): 7 -> 6.7,
                # so all six region-faceted banners land at one on-page cap-height.
                text_x = elem_list_text(colour = "black", size = sz(6.7),
                                        face = "plain"))) +
  scale_fill_manual(values = c(`NHD-up` = unname(PAL_COND[["NHD"]]), `CON-up` = unname(PAL_COND[["CON"]])),
                    name = NULL) +
  scale_shape_manual(values = c("detected, not estimable" = 1,
                                "not tested" = 4), name = NULL) +
  # plotmath group("|", ., "|") renders TRUE vertical bars (the earlier literal
  # "|...|" string rendered as lowercase-l glyphs). ASCII apostrophe kept.
  scale_size_area(name = expression("|Cliff's " * delta * "|"),
                  max_size = sz(2.4), limits = c(0, 0.7),
                  breaks = c(0.2, 0.4, 0.6)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     expand = expansion(mult = c(0.25, 0.42))) +
  coord_cartesian(xlim = c(-L2FC_CAP, L2FC_CAP)) +
  guides(fill  = guide_legend(order = 1),
         size  = guide_legend(order = 2,
                   override.aes = list(fill = "grey55", colour = "white", shape = 21)),
         shape = guide_legend(order = 3, override.aes = list(size = sz(2.4)))) +
  labs(x = expression(avg~log[2]~"FC (NHD/CON)"), y = NULL) +
  theme_pub(base_size = sz(8)) +
  theme(
    plot.title         = element_blank(),
    # gene-row + x-axis labels pre-scaled by the same 5.2/5.4=0.963 factor as the
    # banner above so they land at the sibling on-page cap-height (6.3 -> 6.1).
    axis.text.y        = element_text(size = sz(6.1), face = "italic", colour = "black"),
    axis.text.x        = element_text(size = sz(6.1), colour = "black"),
    axis.title.x       = element_text(size = sz(7), margin = margin(t = 3)),
    axis.line          = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA,
                                       linewidth = 0.25),
    panel.grid.major.x = element_line(colour = "grey95", linewidth = 0.2),
    strip.placement    = "outside",
    panel.spacing.x    = unit(0.25, "lines"),
    panel.spacing.y    = unit(0.30, "lines"),   # a touch more between the two arms
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    # The single-row bottom legend (fill +
    # |Cliff's delta| size key + "below detection" shape key) was wider than the
    # canvas and clipped "below detection" -> "below de".  Wrap to 2 rows
    # (legend.box="vertical" stacks the three guides) and widen the canvas + right
    # margin so the last label sits fully inside the frame.
    legend.box         = "vertical",
    legend.box.just    = "left",
    legend.key.size    = unit(0.22, "cm"),
    legend.spacing.x   = unit(0.12, "cm"),
    # 0.02 cm + a negative top margin left the three stacked guide boxes touching:
    # the "|Cliff's delta|" title printed into the descenders of the CON-up/NHD-up row
    # above it.  Give the stack real vertical breathing room.
    legend.spacing.y   = unit(0.13, "cm"),
    legend.margin      = margin(1, 0, 0, 0),
    # Legend.text was pinned to 6.2 but legend.title was not, so it
    # inherited theme_pub's base size -- the size-key title rendered ~1.5x the rest of
    # the legend, overlapped the fill row above it and pushed the last size break off
    # the canvas.  Pin both, and drop the parenthetical from the title (it is caption
    # text; "|Cliff's delta|" already says "absolute").
    legend.text        = element_text(size = sz(6.2)),
    legend.title       = element_text(size = sz(6.2)),
    plot.margin        = margin(6, 10, 4, 4))

# 20% slimmer (width x0.8), then 20% smaller overall (both x0.8).
# 20% narrower, then 10% smaller, then adopt 2E's font:panel ratio.
#   3.328 x 3.860 -> 2.662 x 3.860 -> 2.396 x 3.474 in
W <- 5.2 * 0.8 * 0.8 * 0.8 * 0.9
# 15% shorter again -> 3.474 -> 2.953 in. Width is unchanged, so the
# width-derived font floor (FS = 5.0/6.1 = 0.820, axis.text 5.00 pt) is deliberately left
# alone -- shortening must not push type back under the 5 pt print floor.
H <- (4.6 * 0.8 + 0.18) * 0.9 * 0.85 * 0.9   # 10% shorter again -> 2.657 in
png_path <- file.path(PANEL, "F2e_micro_DAM_arrest.png")
pdf_path <- file.path(PANEL, "F2e_micro_DAM_arrest.pdf")
ragg::agg_png(png_path, width = W, height = H, units = "in", res = 600)
print(p); dev.off()
ggsave(pdf_path, p, width = W, height = H, useDingbats = FALSE)
cat(sprintf("\nWrote panel: %s\n         and: %s\n", png_path, pdf_path))

# --- Provenance --------------------------------------------------------------
si_path <- file.path(OUT_DIAG, "2I_DAM_arrest_sessionInfo_FH.txt")
writeLines(capture.output({
  cat("Rendered:", format(Sys.time()), "\n")
  cat("PROJ:", PROJ, "\n")
  cat("MAST sources:\n"); for (r in REGIONS) cat("  ", MAST[[r]], "\n")
  cat("ARM 1 (DAM-1 entry):", paste(dam_map$gene[dam_map$arm == arm_levels[1]], collapse=", "), "\n")
  cat("ARM 2 (DAM-2 endpoint):", paste(dam_map$gene[dam_map$arm == arm_levels[2]], collapse=", "), "\n")
  cat("PCT_FLOOR:", PCT_FLOOR, "  L2FC_CAP:", L2FC_CAP, "\n\n")
  print(sessionInfo())
}), si_path)
cat(sprintf("Wrote provenance: %s\n", basename(si_path)))

cat("\n=== DONE ===\n", file = stderr())
