#!/usr/bin/env Rscript
# =============================================================================
# 41_merge_dual_method_FH.R — Merge MAST discovery + pseudobulk, frontal+hippo.
# Faithful adaptation of manuscript_7fig/scripts/41_merge_dual_method.R:
#   * OCC removed everywhere (RGS, LOWN_OCC dropped — Hippo NHD n>>60, not low-n)
#   * pb_confirmed applies to frontal only (sole genuine 2v2 region); Hippo stays
#     discovery + xregion_supported (Hippo per-nucleus dir agrees w/ Frontal pb dir)
#   * dispersion_basis is written by 30_pseudobulk_wald_FH (biological_2v2 / pseudo_split)
#   * inputs/outputs under the rebuild folder
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
DIR <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
MD  <- file.path(DIR, "tables/mast_dual")
FIG <- file.path(DIR, "figures/_dual_method"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
source(file.path(DIR, "scripts/_artifact_genes.R"))
PADJ <- 0.01; LFC <- log2(2.5)
CTS <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
RGS <- c("Frontal","Hippo")     # OCC removed

## --- tier 1: MAST discovery (ambient filtered) ---
mast <- read.csv(file.path(MD, "MAST_dual_discovery_all.csv"))
mast$ambient <- is_artifact(mast$gene, mast$cell_type)
n_amb <- mast %>% filter(discovery %in% TRUE) %>% summarise(a = sum(ambient)) %>% pull(a)
cat(sprintf("ambient-filtered discovery features removed: %d\n", n_amb))
disc <- mast %>% filter(discovery %in% TRUE, !ambient) %>% mutate(mast_sign = sign(avg_log2FC))

## --- tier 2: pseudobulk (Hippo pseudo-split already padj-nulled by script 30_FH) ---
pb <- read.csv(file.path(DIR, "diagnostics/03_DGE_pseudobulk/pseudobulk_DESeq2_with_Wald_stat.csv"))
pb$genuine <- pb$dispersion_basis == "biological_2v2"
cat(sprintf("genuine 2v2 lanes: {%s}\n", paste(sort(unique(pb$comparison[pb$genuine])), collapse=", ")))
pb <- pb %>% mutate(pb_strict = !is.na(padj) & padj < PADJ & abs(log2FoldChange_apeglm) > LFC,
                    pb_sign   = sign(log2FoldChange_apeglm))
if (!"tier_detail" %in% names(pb)) pb$tier_detail <- NA_character_
if (!"low_n_flag"  %in% names(pb)) pb$low_n_flag  <- FALSE
front_conf <- pb %>% filter(grepl("_Frontal$", comparison), genuine, pb_strict) %>%
  transmute(cell_type = sub("_Frontal$", "", comparison), gene, front_pb_sign = pb_sign)

## --- merge ---
merged <- disc %>%
  left_join(pb %>% select(comparison, gene, genuine, pb_strict, pb_sign,
                          pb_lfc = log2FoldChange_apeglm, pb_padj = padj,
                          dispersion_basis, tier_detail, low_n_flag),
            by = c("comparison","gene")) %>%
  left_join(front_conf, by = c("cell_type","gene")) %>%
  mutate(
    pb_confirmed      = region == "Frontal" & genuine %in% TRUE & pb_strict %in% TRUE & pb_sign == mast_sign,
    xregion_supported = region == "Hippo" & !is.na(front_pb_sign) & front_pb_sign == mast_sign,
    support = case_when(pb_confirmed ~ "lane_reproducible",
                        xregion_supported ~ "xregion_supported",
                        TRUE ~ "discovery_only"))
write.csv(merged, file.path(MD, "dual_method_merged.csv"), row.names = FALSE)

## --- counts ---
counts <- merged %>% group_by(cell_type, region, comparison) %>%
  summarise(n_discovery = dplyr::n(),
            disc_up_NHD = sum(mast_sign > 0), disc_up_CON = sum(mast_sign < 0),
            n_confirmed = sum(pb_confirmed),
            conf_up_NHD = sum(pb_confirmed & mast_sign > 0), conf_up_CON = sum(pb_confirmed & mast_sign < 0),
            n_xregion   = sum(xregion_supported), .groups = "drop") %>%
  mutate(confirm_tier    = ifelse(region == "Hippo", "cross-region", "lane-reproducible"),
         n_confirm_shown = ifelse(region == "Hippo", n_xregion, n_confirmed)) %>%
  arrange(factor(cell_type, levels = CTS), factor(region, levels = RGS))
write.csv(counts, file.path(MD, "dual_method_counts.csv"), row.names = FALSE)

cat("\n=== DUAL-METHOD COUNTS: discovery (confirmed) per cell type x region ===\n")
cat("   confirmed = lane-reproducible (Frontal pseudobulk) | cross-region (Hippo)\n")
counts %>% transmute(comparison, discovery = n_discovery,
                     confirmed = as.character(n_confirm_shown), tier = confirm_tier) %>%
  as.data.frame() %>% print(row.names = FALSE)

## --- de-burden panel (2 regions) ---
src_theme <- file.path(DIR, "scripts/22_publication_theme_FH.R")
if (!file.exists(src_theme)) src_theme <- file.path(PROJ, "scripts/22_publication_theme.R")
if (file.exists(src_theme)) try(source(src_theme), silent = TRUE)
PAL <- c(`Up in NHD` = unname(PAL_COND[["NHD"]]), `Up in CON` = unname(PAL_COND[["CON"]]))
thm <- if (exists("theme_pub")) theme_pub(9) else theme_classic(base_size = 9)
# Region facet strips: pale region-tone banner + black plain text + full names
# (Frontal/Hippocampus) via the shared strip_region_x(, clip = "off") helper — matches the
# Fig 2/3 sibling panels (region-banner-pale-fill-black-text house rule). The
# `region` column is remapped to REGION_FULL for the strip label while the
# banner fill is keyed by the short codes RGS (strip_region_x fills in facet
# order, so the factor level order must equal RGS).
REG_FULL_LVL <- if (exists("REGION_FULL")) unname(REGION_FULL[RGS]) else RGS
base <- merged %>% mutate(dir = ifelse(mast_sign > 0, "Up in NHD", "Up in CON")) %>%
  group_by(cell_type, region, dir) %>%
  # region-aware confirmed overlay (round-4 fix): Frontal is the sole 2v2 region so
  # its confirmed tier = pb_confirmed (Frontal pseudobulk-reproducible). Hippo has a
  # single NHD lane, so its confirmed tier = xregion_supported (per-nucleus direction
  # agrees with the Frontal pseudobulk sign). Previously conf = sum(pb_confirmed) was
  # structurally 0 for every Hippo row (pb_confirmed hard-requires region=='Frontal'),
  # leaving the Hippo facet with no dark overlay. Both tiers now map to one plotted
  # layer so both facets read as the same two-tone encoding. Single-NHD-hippo-lane
  # caveat stays in the legend, not on the plot.
  summarise(Discovery = dplyr::n(),
            conf = sum(ifelse(region == "Hippo", xregion_supported, pb_confirmed)),
            .groups = "drop") %>%
  mutate(cell_type = factor(cell_type, levels = rev(CTS)),
         region    = factor(if (exists("REGION_FULL")) REGION_FULL[region] else region,
                            levels = REG_FULL_LVL),
         sgn       = ifelse(dir == "Up in CON", -1, 1))
# a shared +/-max(Discovery) axis left Hippocampus half-empty). geom_blank rows
# force each facet symmetric around 0 out to its own max |bar| + 6% headroom.
facet_lim <- base %>% group_by(region) %>%
  summarise(lim = max(Discovery, na.rm = TRUE) * 1.06, .groups = "drop")
blank_df <- facet_lim %>% tidyr::crossing(sgn2 = c(-1, 1)) %>%
  transmute(region, cell_type = factor(CTS[1], levels = rev(CTS)),
            dir = "Up in NHD", x = sgn2 * lim)
# round-4 fix: the solid overlay now covers both confirmation tiers (Frontal
# lane-reproducible + Hippo cross-region), so the swatch is relabelled
# "Cross-method confirmed" to describe both facets with one key.
SHADE <- c("Discovery (per-nucleus)" = 0.30, "Cross-method confirmed" = 1)
# The two
# alpha swatches now render as explicit distinct grey blocks (pale grey =
# discovery, solid grey = confirmed) via a fixed-fill override, so the
# solid-vs-pale meaning reads at swatch level instead of being nearly invisible.
p <- ggplot(base, aes(y = cell_type, fill = dir)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_col(aes(x = sgn * Discovery, alpha = "Discovery (per-nucleus)"), width = 0.7, linewidth = 0) +
  geom_col(aes(x = sgn * conf, alpha = "Cross-method confirmed"), width = 0.7, linewidth = 0) +
  geom_blank(data = blank_df, aes(x = x), inherit.aes = FALSE) +
  ggh4x::facet_wrap2(~ region, nrow = 1, scales = "free_x",
                     # Region banner labels +20% (7 -> 8.4pt),
                     # then a further +20% (8.4 -> 10.1pt).
                     strip = strip_region_x(RGS, fontsize = 10.1, clip = "off")) +
  scale_fill_manual(values = PAL, name = NULL) +
  scale_alpha_manual(values = SHADE, name = NULL, breaks = c("Cross-method confirmed","Discovery (per-nucleus)")) +
  # house rule: no code tokens on the y-axis -> Neuron_Ex -> "Neuron Ex", etc.
  # (ct_disp from the FH theme; factor LEVELS stay coded for ordering/joins).
  scale_y_discrete(labels = if (exists("ct_disp")) ct_disp else waiver()) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  labs(x = "per-nucleus discovery features", y = NULL) +
  thm +
  # Sibling banner parity: ggh4x facet_wrap2 free_x can size the two
  # region strip texts unequally (Hippocampus rendered ~1.7x Frontal cap-height
  # even though strip_region_x(RGS, clip = "off") sets a single fixed 7pt). Force one shared
  # strip text size in the panel theme after strip_themed so both banners render
  # at identical 7pt cap-height (house rule: region banners identical across
  # sibling facets). Fill/colour still come from strip_region_x (pale + black).
  # (size raised to 10.1 with strip_region_x above — this override is what actually
  # binds, so the two must be kept equal or the banners silently revert.)
  theme(strip.text = element_text(size = 10.1, colour = "black", face = "plain")) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(1,"pt")) +
  guides(fill = guide_legend(order=1, override.aes=list(alpha=1)),
         # explicit 2-swatch saturation key: solid grey vs pale grey blocks, both
         # drawn at fixed fill "grey25" so only the alpha (saturation) differs.
         alpha = guide_legend(order=2, override.aes=list(fill="grey25", colour = NA),
                              keywidth = unit(0.42, "cm"), keyheight = unit(0.30, "cm")))
# Panel 10% smaller overall (0.90 on both axes), then a further
# 10% slimmer (0.90 on width only).
PW <- 9 * 0.8 * 0.9 * 0.95 * 0.9 * 0.9; PH <- 4.4 * 0.8 * 0.9
ggsave(file.path(FIG, "deburden_dual_method.png"), p, width = PW, height = PH, dpi = 300, bg = "white")
ggsave(file.path(FIG, "deburden_dual_method.pdf"), p, width = PW, height = PH, bg = "white")
PANEL1 <- file.path(DIR, "figures", "Figure_1", "panels")
dir.create(PANEL1, recursive = TRUE, showWarnings = FALSE)
for (ext in c("png","pdf"))
  file.copy(file.path(FIG, paste0("deburden_dual_method.", ext)),
            file.path(PANEL1, paste0("F1c_deburden_dual_method.", ext)), overwrite = TRUE)
cat("copied panel c -> figures/Figure_1/panels/F1c_deburden_dual_method.{png,pdf}\n")
cat("\n=== DONE ===\n", file = stderr())
