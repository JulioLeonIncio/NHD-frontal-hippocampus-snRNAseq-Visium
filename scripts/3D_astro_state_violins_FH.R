#!/usr/bin/env Rscript
# =============================================================================
# 3D_astro_state_violins_FH.R — Figure 3 panel d: astrocyte reactive-state signature violins + stats, frontal + HIPPOCAMPUS (2 regions). Port of
# fig4_4b_v2_astro_violins.R.  Writes the delta/q/n/star stats table that the
# companion dot-matrix (3D_astro_reactive_dotmatrix_FH.R) reads (numbers match).
# -----------------------------------------------------------------------------
# Scores 11 astrocyte-state signatures per nucleus (Metallothionein-loss headline
# = MT2A/MT3, pan-reactive/DAA framing, never A1/A2 as biology).  Per (signature x
# Region): per-nucleus Cliff's delta + BH-FDR Wilcoxon q, and the house gated star
# (*/**/*** when q<alpha & |d|>=0.15; grey "q*" when q<alpha & |d|<0.15; n.s.
# otherwise).  no per-nucleus significance STARS as a donor claim — the effect
# size (Cliff's delta) is the finding; single-donor pseudoreplication caveat ->
# manuscript text.  The violin renders a 5-signature subset; the dot-matrix shows
# all 11.
#
# State gene-set re-validation — each set's present/dropped genes
# are logged by score_sig() (DROPPED = absent from the FH atlas).  Curation vs the
# 3-region reference (drops the FH-absent genes; keeps the biology):
#   * Metallothionein: MT2A/MT3 headline retained; MT1H absent-in-atlas dropped.
#   * Pan-reactive / DAA / A1 / A2 / STAT3 / Zhou-NHD: from curated_signatures.rds
#     (already vetted); absent genes dropped by score_sig, reported.
#   * Homeostatic astro: SLC1A3 kept in the module score (a module average absorbs
#     its Frontal direction quirk; the volcano handles the per-gene direction).
#   * Interferon / Heat-shock / Synaptogenic: curated inline; absent dropped.
#
# Cache: data/_cache_Fig3_astro_scores_FH.rds (mtime + signature-completeness guard).
# Output: figures/Figure_3/panels/F3e_astro_state_violins.{pdf,png}
#         + tables/fig3d_astro_signature_stats.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(scales); library(ggh4x)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, REGION_FULL, strip_region_x, REGION_ORDER

PANEL <- file.path(PROJ, "figures", "Figure_3", "panels")
TDIR  <- file.path(PROJ, "tables")
DDIR  <- file.path(PROJ, "data")
LOGS  <- file.path(PROJ, "logs")
CACHE <- file.path(DDIR, "_cache_Fig3_astro_scores_FH.rds")
ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
SIGP  <- file.path(DDIR, "curated_signatures.rds")
for (d in c(PANEL, TDIR, DDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING atlas" = file.exists(ATLAS), "MISSING curated_signatures.rds" = file.exists(SIGP))

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")

OR_CON_FILL <- "#82B2D6"; OR_NHD_FILL <- "#D6837A"

# The A1/A2 dichotomy is deprecated (Escartin 2021). Relabel
# the reactive rows to their provenance-anchored names so the panel never asserts the
# A1/A2 biology:  A1 reactive -> "Complement/IFN-reactive (Liddelow)";
#                 A2 reactive -> "Ischemic/S100A10-reactive (Zamanian)".
# Metallothionein carries the headline down direction -> "Metallothionein (lost)".
# These strings are the CSV `signature` keys, so the companion dot-matrix's SIG_LEVELS
# must match verbatim (it does).
# The six inline sets and the five published-set panel labels are declared once in
# scripts/_astro_inline_sets_FH.R, which 90_build_supplementary_tables_FH.R also sources so that
# ST6 (Supplementary Data 7) ships exactly what this panel scores. Edit them there, not here.
source(file.path(PROJ, "scripts", "_astro_inline_sets_FH.R"))   # astro_inline_sets(), astro_published_panel_names()
PUB_LAB <- astro_published_panel_names()
LAB_A1 <- unname(PUB_LAB[["A1_reactive"]])
LAB_A2 <- unname(PUB_LAB[["A2_reactive"]])
LAB_MT <- "Metallothionein (lost)"

# 6 new state sets (inline; absent genes dropped by score_sig + reported).
NEW_SIGS <- astro_inline_sets()
stopifnot("the shared helper no longer carries the Metallothionein set under its panel label" = LAB_MT %in% names(NEW_SIGS),
          "the shared helper must carry exactly six inline astrocyte sets" = length(NEW_SIGS) == 6L)

ALL_SIG_ORDER <- c(LAB_A1,"Pan-reactive",LAB_A2,"DAA (Habib)",
                   "Zhou NHD astrocyte","Interferon (Hasel)","STAT3 targets",
                   LAB_MT,"Homeostatic astro","Heat-shock (HSF1)","Synaptogenic")
VIOLIN_SIG_ORDER <- c("Pan-reactive","DAA (Habib)",LAB_MT,
                      "Homeostatic astro","Heat-shock (HSF1)")
# the eleven panel rows are exactly the five published labels + the six inline sets of the helper
stopifnot("ALL_SIG_ORDER drifted from the shared helper's names" =
            setequal(ALL_SIG_ORDER, c(unname(PUB_LAB), names(NEW_SIGS))) && length(ALL_SIG_ORDER) == 11L)

# ---- mtime + signature-completeness cache guard ----------------------------
cache_stale <- !file.exists(CACHE) || file.mtime(CACHE) < file.mtime(ATLAS)
if (!cache_stale) {
  have <- unique(readRDS(CACHE)$signature)
  if (!all(ALL_SIG_ORDER %in% have)) {
    cat(sprintf("== cache missing %d signature(s) -> rebuild ==\n", sum(!ALL_SIG_ORDER %in% have)))
    cache_stale <- TRUE
  }
}
if (cache_stale) {
  cat("== _cache_Fig3_astro_scores_FH.rds stale/missing -> rebuild from atlas ==\n")
  suppressPackageStartupMessages(library(Seurat))
  sigs <- readRDS(SIGP); obj <- readRDS(ATLAS)
  DefaultAssay(obj) <- "RNA"
  if (length(SeuratObject::Layers(obj, assay = "RNA")) > 1) obj <- JoinLayers(obj, assay = "RNA")
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  ast <- subset(obj, subset = new_annotation == "Astro" & Region %in% REGIONS)
  cat(sprintf("   astro nuclei: %d\n", ncol(ast))); print(table(Region = ast$Region, Condition = ast$Condition))
  score_sig <- function(so, genes, name) {
    gs <- intersect(genes, rownames(so)); dropped <- setdiff(genes, rownames(so))
    cat(sprintf("   [%-19s] %2d/%2d genes present%s\n", name, length(gs), length(genes),
                if (length(dropped)) paste0("; DROPPED (absent): ", paste(dropped, collapse = ", ")) else ""))
    if (length(gs) < 1L) stop("signature '", name, "' has 0 genes in the object")
    so <- AddModuleScore(so, features = list(gs), name = "sig", seed = 42, assay = "RNA")
    data.frame(cell = colnames(so), score = so$sig1, Condition = so$Condition,
               Region = so$Region, cell_type = "Astro", signature = name, stringsAsFactors = FALSE)
  }
  stopifnot("curated_signatures.rds lacks a published astrocyte set named by the shared helper" = all(names(PUB_LAB) %in% names(sigs)))
  score_calls <- c(
    # the five published sets, scored under their panel labels (A1_reactive -> LAB_A1, ...)
    lapply(names(PUB_LAB), function(k) score_sig(ast, sigs[[k]], unname(PUB_LAB[[k]]))),
    lapply(names(NEW_SIGS), function(nm) score_sig(ast, NEW_SIGS[[nm]], nm)))
  glial <- bind_rows(score_calls)
  tmp <- paste0(CACHE, ".tmp"); saveRDS(glial, tmp); file.rename(tmp, CACHE)
  cat(sprintf("   wrote %s (%d rows, %d signatures)\n", basename(CACHE), nrow(glial), length(unique(glial$signature))))
  rm(obj, ast, glial); gc()
}

astro <- readRDS(CACHE) %>% filter(cell_type == "Astro") %>%
  mutate(donor     = sub("_.*$", "", cell),
         Condition = factor(Condition, levels = c("CON","NHD")),
         Region    = factor(Region,    levels = REGIONS))
stopifnot("cache missing an expected signature" = all(ALL_SIG_ORDER %in% unique(astro$signature)))
astro <- astro %>% filter(signature %in% ALL_SIG_ORDER) %>%
  mutate(signature = factor(signature, levels = ALL_SIG_ORDER))
cat("Scoring stats for signatures (11):\n"); print(ALL_SIG_ORDER)

# ---- Cliff's delta (rank-based, vectorized) --------------------------------
cliff_delta <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]; nx <- length(x); ny <- length(y)
  if (nx < 3 || ny < 3) return(NA_real_)
  r <- rank(c(x, y)); Rx <- sum(r[seq_len(nx)]); U <- Rx - nx*(nx+1)/2
  (2*U/(nx*ny)) - 1
}

stat_df <- astro %>% group_by(signature, Region) %>%
  summarise(
    n_CON = sum(Condition == "CON"), n_NHD = sum(Condition == "NHD"),
    med_CON = median(score[Condition == "CON"], na.rm = TRUE),
    med_NHD = median(score[Condition == "NHD"], na.rm = TRUE),
    p_wilcox = tryCatch(wilcox.test(score[Condition == "CON"], score[Condition == "NHD"],
                                    exact = FALSE)$p.value, error = function(e) NA_real_),
    cliff_d = cliff_delta(score[Condition == "NHD"], score[Condition == "CON"]),
    .groups = "drop") %>%
  mutate(q_BH = p.adjust(p_wilcox, method = "BH"))   # BH over 11 sigs x 2 regions = 22 tests

stat_df <- stat_df %>% mutate(
  abs_d = abs(cliff_d), passes_effect = !is.na(abs_d) & abs_d >= 0.15,
  # Drop the per-nucleus
  # wilcox ***/**/* ladder — that q is a PSEUDOREPLICATED per-nucleus test with a
  # single NHD donor and cannot support a donor-level significance claim
  # (pi/stats policy: rare N=1 -> effect sizes only, no per-nucleus stars).
  # Canonical convention (matches 2B_micro_program_violins_FH.R ~L102): keep the
  # Cliff's delta label (the effect size) + only the grey effect-gated "q*" mark
  # (q<0.05 and |delta|<0.15, i.e. significant-but-small). q_BH stays in the CSV.
  star = case_when(q_BH < 0.05 & !passes_effect ~ "q*", TRUE ~ ""))

cat("\nPer-facet stats (BH-FDR + Cliff's delta):\n")
print(as.data.frame(stat_df %>% select(signature, Region, n_CON, n_NHD, med_CON, med_NHD,
        p_wilcox, q_BH, cliff_d, star) %>% mutate(across(where(is.numeric), ~signif(.x, 3)))),
      row.names = FALSE)
write.csv(stat_df, file.path(TDIR, "fig3d_astro_signature_stats.csv"), row.names = FALSE)
cat(sprintf("Wrote stats CSV: %d rows (%d signatures x %d regions)\n",
            nrow(stat_df), length(ALL_SIG_ORDER), length(REGIONS)))

# ---- violin scope: 5-signature subset --------------------------------------
astroV <- astro %>% filter(signature %in% VIOLIN_SIG_ORDER) %>%
  mutate(signature = factor(as.character(signature), levels = VIOLIN_SIG_ORDER))
statV  <- stat_df %>% filter(signature %in% VIOLIN_SIG_ORDER) %>%
  mutate(signature = factor(as.character(signature), levels = VIOLIN_SIG_ORDER))

donor_means <- astroV %>% group_by(signature, Region, Condition, donor) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop")
y_top <- bind_rows(
    astroV %>% mutate(src = "nucleus"),
    donor_means %>% mutate(src = "donor", cell = NA_character_, cell_type = "Astro")) %>%
  group_by(signature, Region) %>%
  summarise(y_star = quantile(score, 0.995, na.rm = TRUE),
            y_min = min(score, na.rm = TRUE), y_max = max(score, na.rm = TRUE), .groups = "drop")
# Directly over the NHD violin it annotates (x=2), at a small headroom above that
# group's own upper whisker/box top — not floating at the violin midpoint (x=1.5) in
# dead centre space. This mirrors how the volcano/dotmatrix/GSEA siblings key q* to
# its own dot. y_nhd = NHD-group 99th percentile (a robust "top of the box/whisker"
# that ignores the extreme violin tail) + a small headroom.
y_nhd_df <- astroV %>% filter(Condition == "NHD") %>%
  group_by(signature, Region) %>%
  summarise(y_nhd = quantile(score, 0.99, na.rm = TRUE), .groups = "drop")
statV <- statV %>% left_join(y_top, by = c("signature","Region")) %>%
  left_join(y_nhd_df, by = c("signature","Region")) %>%
  mutate(delta_lab = sprintf("italic(delta) == '%+.2f'", cliff_d),
         y_top   = y_max + 0.06 * (y_max - y_min),
         y_delta = y_max + 0.19 * (y_max - y_min),
         y_qstar = y_nhd + 0.05 * (y_max - y_min))

# Hippocampus astro rests on a single NHD lane -> flag the
# Hippo facet low-confidence FIGURE-WIDE (grey banner + "(low n)" suffix), the same
# treatment used on the astro volcano / dotmatrix / oligo panels. One consistent
# Hippo caveat across every astro and oligo panel of Figure 3.
low_conf <- "Hippo"

pB <- ggplot(astroV, aes(x = Condition, y = score, fill = Condition)) +
  geom_violin(scale = "width", linewidth = 0.25, alpha = 0.85, trim = FALSE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "grey25",
               linewidth = 0.3, fatten = 1.8) +
  geom_jitter(data = donor_means, aes(x = Condition, y = score), inherit.aes = FALSE,
              width = 0.08, size = 0.6, color = "grey25", fill = "white",
              shape = 21, stroke = 0.3, alpha = 0.95) +
  # Only the house effect-gated grey "q*" mark is drawn (no per-nucleus
  # */**/*** stars, no "n.s."). Cliff's delta below carries the effect.
  # placed over the NHD violin (x=2) at y_qstar (just above the NHD box/whisker), so it
  # keys to its own group like the volcano/dotmatrix/GSEA siblings — not floating at
  # the violin midpoint.
  geom_text(data = subset(statV, star == "q*"), aes(x = 2, y = y_qstar, label = star),
            inherit.aes = FALSE, size = 1.9, vjust = 0, color = "grey50") +
  geom_text(data = statV, aes(x = 0.42, y = y_delta, label = delta_lab), inherit.aes = FALSE,
            # delta labels were the largest text in the panel and read louder than the data;
            # keeps them legible but subordinate to the violins.
            parse = TRUE, hjust = 0, vjust = 0, size = 1.75, color = "black") +
  # Read horizontally -- transposed from signature-rows x region-columns
  # (tall, narrow) to region-ROWS x signature-COLUMNS (wide, short). Region therefore moves from
  # the x strip to the y strip, so the pale region banner has to be supplied as background_y;
  # strip_region_x(, clip = "off") only styles x strips and would have left the regions unstyled.
  facet_grid2(Region ~ signature, scales = "free_y", switch = "y",
              labeller = labeller(Region = region_full_lowconf_labeller(low_conf)),
              strip = ggh4x::strip_themed(
                clip = "off",
                background_y = ggh4x::elem_list_rect(
                  fill = unname(PAL_REGION_PALE[REGIONS]), colour = NA),
                text_y = ggh4x::elem_list_text(colour = "black", size = 7, face = "plain",
                                               angle = 90),
                background_x = ggh4x::elem_list_rect(fill = "grey96", colour = NA),
                text_x = ggh4x::elem_list_text(colour = "black", size = 6.2, face = "plain"))) +
  scale_fill_manual(values = c("CON" = OR_CON_FILL, "NHD" = OR_NHD_FILL), name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.30))) +
  labs(x = NULL, y = NULL) +
  # Base 7.5 -> 8 and axis text -> 6.3, matching the
  # house spec every other Fig-3 panel now uses (2D_secretome reference).
  theme_pub(base_size = 8) +
  theme(
    # strip styling now comes from strip_themed() above (region pale on y, signatures grey on x)
    strip.placement    = "outside", strip.clip = "off",
    axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 6.3), axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.spacing.x = unit(0.6, "lines"), panel.spacing.y = unit(0.25, "lines"),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.margin = margin(-2, 0, 0, 0), legend.key.size = unit(0.22, "cm"),
    legend.text = element_text(size = 6.2), plot.margin = margin(6, 4, 4, 4))

BN <- "F3e_astro_state_violins"
# Widened 3.0 -> 3.6 in so the "Hippocampus (low n)" banner
# (clip="off") is not clipped at the right panel edge.
# 20% slimmer (3.60 -> 2.88); height unchanged.
# 20% smaller (both dimensions) -> 4.46 x 1.84. Type is left at the house spec
# (base 8 / axis 6.3 / legend 6.2), which is absolute, so it stays legible; the delta labels
# were already reduced to 1.75 in the previous pass.
W <- 5.58 * 0.8; H <- 2.30 * 0.8
ggsave(file.path(PANEL, paste0(BN, ".pdf")), pB, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), pB, width = W, height = H, dpi = 600)
cat(sprintf("\nWrote %s.{pdf,png} + fig3d_astro_signature_stats.csv\n", BN))
cat("=== DONE ===\n", file = stderr())
