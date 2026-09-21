# =============================================================================
# 2B_micro_program_violins_FH.R — Figure 2 panel b (leads): microglial-PROGRAM signature violins across frontal + HIPPOCAMPUS (2 regions). Port of
# 57_fig2_micro_program_violins.R from the 3-region original.
# -----------------------------------------------------------------------------
# rows = regions, columns = 7 curated microglial PROGRAMS; CON vs NHD violins +
# boxplot; per-facet Wilcoxon q + Cliff's delta (effect-size-gated stars, house
# "q*" mark).  Signatures scored with AddModuleScore in bona-fide Micro-PVM.
#
# gene-set re-validation:
#   All 7 programme gene sets were re-checked against the FH Micro-PVM MAST tables
#   (Frontal + Hippo).  result: every gene in every programme is detected (max pct
#   >= 0.10) and moves in the programme's expected direction in >=1 region ->
#   nothing dropped, the 7 programmes are carried verbatim.  AddModuleScore uses
#   the genes actually present in the atlas (intersect), logged below.  This is the
#   DAM-1-arrest / TREM2-independent-inflammation signature at the programme level:
#   Antigen presentation / iron / phagolysosomal / glycolytic / inflammatory up;
#   homeostatic + immunoregulatory-brake lost.
#
# Atlas-based.  Output: F2d_micro_program_violins.{png,pdf} (+ stats CSV).
# -----------------------------------------------------------------------------
# apparent-type lift 2B_micro_program_violins_FH.R: every text size in this script scaled
# by 1.21 so the panel reads at ~5.0 pt on the assembled page, matching the
# Figure-5 dot panels. Point sizes, line widths and unit() dimensions are not
# touched. Canvas size unchanged, so re-linking is a no-op.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tidyr); library(ggplot2); library(ggh4x); library(ragg)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # PAL_REGION, REGION_FULL, strip_region_x, REGION_ORDER
ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
PANEL <- file.path(PROJ, "figures/Figure_2/panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
stopifnot("MISSING atlas" = file.exists(ATLAS))

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")

# --- the 7 microglial PROGRAMS (identical to the volcano MICRO_THEME; FH-revalidated) ---
PROGRAMS <- list(
  "Antigen presentation"           = c("CD74","HLA-DRA","HLA-DRB1","HLA-DMB"),
  "Metal handling (iron/zinc)"     = c("FTH1","FTL","HAMP","TMEM163"),
  "Phagocytic / lysosomal"         = c("FCGR3A","SCIN","ADAM28","CTSB"),
  "Glycolytic shift"               = c("PFKFB3","SLC2A3"),
  "Inflammatory / stress response" = c("CEBPD","ZFP36L1","RGS1","SRGN","HSPA1A","TNFRSF1B","NAIP"),
  "Homeostatic (lost)"             = c("P2RY12","MEF2C","PLXDC2","SORL1"),
  "Immunoregulatory brake (lost)"  = c("IRAK3","LDLRAD4","ZBTB16","HDAC9"))
PROG_ORDER <- names(PROGRAMS)

OR_CON_FILL <- unname(PAL_COND_PALE[["CON"]])   # theme single-source
OR_NHD_FILL <- unname(PAL_COND_PALE[["NHD"]])   # theme single-source

cat("loading atlas...\n"); atl <- readRDS(ATLAS)
# microglia-only — see 40b/41b. PVM are split out of the compartment
# because DAM/programme staging is undefined for them (Van Hove 2019 PMID 31061494);
# oligo-doublets are already excluded by the assignment. Fails loud if it is missing.
.cls_f <- file.path(PROJ, "tables", "micro_states", "immune_umap_classes_FH.csv")
if (!file.exists(.cls_f)) stop("MISSING ", .cls_f, " — run 2N_immune_UMAP_FH.R first")
.cls <- read.csv(.cls_f, stringsAsFactors = FALSE)
.keep_micro <- intersect(.cls$barcode[.cls$class == "Microglia"], colnames(atl))
.keep_pvm   <- intersect(.cls$barcode[.cls$class == "PVM"],       colnames(atl))
if (length(.keep_micro) < 100) stop("microglia assignment matched only ",
                                    length(.keep_micro), " cells")
o_pvm <- subset(atl, cells = c(.keep_micro, .keep_pvm))   # kept for the delta_PM audit
o <- subset(atl, cells = .keep_micro); rm(atl); gc()
cat(sprintf("microglia nuclei (PVM excluded): %d\n", ncol(o)))
DefaultAssay(o) <- "RNA"; o <- JoinLayers(o); o <- NormalizeData(o, verbose = FALSE)
o$Region    <- factor(o$Region, levels = REGIONS)
o$Condition <- factor(o$Condition, levels = c("CON","NHD"))
present <- lapply(PROGRAMS, function(g) intersect(g, rownames(o)))
cat("== programme genes present in atlas (used by AddModuleScore) ==\n")
for (nm in PROG_ORDER)
  cat(sprintf("  %-32s %d/%d: %s\n", nm, length(present[[nm]]), length(PROGRAMS[[nm]]),
              paste(present[[nm]], collapse = ", ")))
o <- AddModuleScore(o, features = present, name = "Prog_", seed = 42)
score_cols <- paste0("Prog_", seq_along(PROGRAMS))

md <- o@meta.data %>%
  select(Region, Condition, SampleID = tidyselect::any_of(c("SampleID","donor","orig.ident")),
         all_of(score_cols)) %>%
  rename_with(~ PROG_ORDER, all_of(score_cols)) %>%
  tidyr::pivot_longer(all_of(PROG_ORDER), names_to = "program", values_to = "score") %>%
  mutate(program = factor(program, levels = PROG_ORDER),
         Region  = factor(Region, levels = REGIONS))

# --- per-nucleus Wilcoxon + Cliff's delta per (program, Region) --------------
cliff_delta <- function(x, y) {          # delta = 2*AUC - 1, bounded [-1,1]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]; nx <- length(x); ny <- length(y)
  if (nx < 3 || ny < 3) return(NA_real_)
  r <- rank(c(x, y)); U <- sum(r[seq_len(nx)]) - nx*(nx+1)/2   # midranks (ties handled by rank())
  (2*U/(nx*ny)) - 1
}
stat_df <- md %>% group_by(program, Region) %>%
  summarise(p_wilcox = tryCatch(wilcox.test(score[Condition=="CON"],
                                            score[Condition=="NHD"], exact=FALSE)$p.value,
                                error = function(e) NA_real_),
            cliff_d = cliff_delta(score[Condition=="NHD"], score[Condition=="CON"]),
            y_star  = quantile(score, 0.995, na.rm = TRUE),
            y_min   = min(score, na.rm = TRUE), y_max = max(score, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(q_BH = p.adjust(p_wilcox, method = "BH"),
         passes = !is.na(cliff_d) & abs(cliff_d) >= 0.15,
         # Drop the per-nucleus
         # wilcox ***/**/* stars — that p is a PSEUDOREPLICATED per-nucleus test with a
         # single NHD donor, so it cannot support a donor-level significance claim
         # (pi/stats policy: rare N=1 -> effect sizes only, no per-nucleus stars).
         # We keep the Cliff's delta label (the effect size) and the house effect-gated
         # "q*" mark (grey: q<0.05 and |delta|<0.15, i.e. significant-but-small).  The
         # q_BH column is still written to the stats CSV for the record, but the only
         # on-panel annotations are delta + q*.
         star = case_when(q_BH < 0.05 & !passes ~ "q*", TRUE ~ ""),
         delta_lab = sprintf("italic(delta) == '%+.2f'", cliff_d),
         y_top   = y_max + 0.06 * (y_max - y_min),
         y_delta = y_max + 0.19 * (y_max - y_min))
write.csv(stat_df, file.path(TDIR, "fig2_micro_program_violin_stats_FH.csv"), row.names = FALSE)
cat("\n== per-facet stats (q_BH, Cliff's delta, star) ==\n")
print(as.data.frame(stat_df[, c("program","Region","q_BH","cliff_d","star")]), row.names = FALSE)

pC <- ggplot(md, aes(Condition, score, fill = Condition)) +
  geom_violin(scale = "width", linewidth = 0.25, alpha = 0.85, trim = FALSE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               colour = "grey25", linewidth = 0.3, fatten = 1.8) +
  # Only the house effect-gated "q*" mark is drawn on the panel; the
  # per-nucleus wilcox stars are gone.  Cliff's delta below carries the effect.
  geom_text(data = subset(stat_df, star == "q*"),
            aes(1.5, y_top, label = star), inherit.aes = FALSE,
            size = 2.2, vjust = 0.4, colour = "grey50") +
  geom_text(data = stat_df, aes(0.42, y_delta, label = delta_lab), inherit.aes = FALSE,
            parse = TRUE, hjust = 0, vjust = 0, size = 2.1, colour = "black") +
  # Regions = ROWS (pale region banners, black text), programs = columns.
  facet_grid2(Region ~ program, scales = "free_y",
              labeller = labeller(
                Region  = as_labeller(REGION_FULL),
                # Split before " (" and at " / " as before, then wrap any line that
                # is still wider than the facet. Without the wrap step
                # "Immunoregulatory brake" clipped to "mmunoregulatory brak"
                # [[facet-strip-titles-must-fit-panel-width]] — 20 chars is the widest
                # line that fits at this panel width ("Antigen presentation" is exactly
                # 20 and must stay on one line).
                program = as_labeller(function(x) {
                  s <- gsub(" / ", "/\n", gsub(" \\(", "\n(", x))
                  vapply(strsplit(s, "\n", fixed = TRUE), function(parts)
                    paste(vapply(parts, function(p)
                      paste(strwrap(p, width = 20), collapse = "\n"),
                      character(1)), collapse = "\n"),
                    character(1))
                })),
              strip = ggh4x::strip_themed(
                # Clip = "off": the right-hand banner was rendering
                # "lippocampu" -- both the H and the trailing s cut. theme(strip.clip=)
                # does not reach ggh4x strip constructors, so the argument has to be
                # given here for the centred label to overflow the pale banner instead
                # of being truncated.
                clip = "off",
                background_x = ggh4x::elem_list_rect(fill = "grey93", colour = NA),
                text_x = ggh4x::elem_list_text(size = 6.8, colour = "black", lineheight = 0.85),
                # The single-NHD-hippo-lane caveat lives in the legend, not
                # as a greyed-out on-plot banner.  both region banners now use their normal pale region
                # tone (PAL_REGION_PALE) so all six Fig-2 region-faceted panels read as
                # siblings (2D/2F/2_zhou already do this).
                background_y = ggh4x::elem_list_rect(
                  fill = unname(PAL_REGION_PALE[REGIONS]), colour = NA),
                text_y = ggh4x::elem_list_text(size = 8.5, colour = "black", angle = -90))) +
  scale_fill_manual(values = c("CON" = OR_CON_FILL, "NHD" = OR_NHD_FILL), name = NULL) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.30))) +
  labs(x = NULL, y = "module score") +
  theme_pub(base_size = 9.1) +
  theme(plot.title       = element_blank(),   # house rule: no in-panel title
        strip.clip       = "off",
        axis.text.x      = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y      = element_text(size = 6.5),
        axis.title.y     = element_text(size = 8.5),
        axis.line        = element_blank(),
        panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x  = unit(0.28, "lines"), panel.spacing.y = unit(0.28, "lines"),
        legend.position  = "bottom", legend.key.size = unit(0.24, "cm"),
        legend.text      = element_text(size = 8.5),
        plot.margin      = margin(6, 4, 4, 4))
# 2 region rows -> shorter.
# 10% slimmer (width only). 6.80 -> 6.12 in; height unchanged.
W <- 6.8 * 0.9; H <- 2.0
ragg::agg_png(file.path(PANEL, "F2d_micro_program_violins.png"), width = W, height = H,
              units = "in", res = 600); print(pC); invisible(dev.off())
ggsave(file.path(PANEL, "F2d_micro_program_violins.pdf"), pC, width = W, height = H, useDingbats = FALSE)
cat("wrote F2d_micro_program_violins.{png,pdf}\n")

# =============================================================================
# delta_PM AUDIT (risk 2 — the number that decides whether
# this panel is safe). The composition bias a programme score can absorb is
# c * delta_PM, where c is the mixing-fraction term (Frontal 0.124 / Hippo 0.097)
# and delta_PM is Cliff's delta of the score between PVM and microglia. A score
# built from co-enriched PVM genes can have a large delta_PM even when no single
# gene is confounded, so it must be measured per programme rather than assumed.
# Computed in CON CELLS only, so the condition contrast never enters it.
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
op <- o_pvm
op$.class <- ifelse(colnames(op) %in% .keep_pvm, "PVM", "Microglia")
op <- op[, op$Condition == "CON"]
DefaultAssay(op) <- "SCT"
prog_present <- lapply(PROGRAMS, function(g) intersect(g, rownames(op)))
op <- Seurat::AddModuleScore(op, features = prog_present, name = "PMAUD", seed = 42, ctrl = 50)
sc_cols <- paste0("PMAUD", seq_along(prog_present))

cliff <- function(x, isP) {
  n1 <- sum(isP); n2 <- sum(!isP)
  if (n1 < 5 || n2 < 5) return(NA_real_)
  r <- rank(x); U1 <- sum(r[isP]) - n1*(n1+1)/2; 2*U1/(n1*n2) - 1
}
isP <- op$.class == "PVM"
aud <- tibble::tibble(
  programme = names(PROGRAMS),
  n_genes   = vapply(prog_present, length, integer(1)),
  delta_PM  = vapply(sc_cols, function(cc) cliff(op@meta.data[[cc]], isP), numeric(1)),
  n_PVM     = sum(isP), n_micro = sum(!isP)) %>%
  mutate(max_bias_Frontal = round(0.1243 * delta_PM, 4),
         max_bias_Hippo   = round(0.0966 * delta_PM, 4),
         # a programme is flagged if composition alone could move it an appreciable
         # fraction of the 0.15 effect gate
         flag = ifelse(abs(max_bias_Frontal) >= 0.05, "CHECK", "ok"))
AUD_OUT <- file.path(PROJ, "tables", "micro_states", "programme_deltaPM_audit_FH.csv")
dir.create(dirname(AUD_OUT), recursive = TRUE, showWarnings = FALSE)
write.csv(aud, AUD_OUT, row.names = FALSE)
cat("\n== delta_PM audit (PVM vs microglia, CON cells only) ==\n")
print(as.data.frame(aud), digits = 3)
cat(sprintf("wrote %s\n", AUD_OUT))
cat("=== DONE ===\n", file = stderr())
