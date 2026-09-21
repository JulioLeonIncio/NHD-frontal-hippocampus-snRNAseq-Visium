# =============================================================================
# 2E_states_MAST_GSEA_FH.R — Figure 2 panel e: MAST-forward curated microglia state enrichment by effect-size preranked GSEA (signed winsorized avg_log2FC),
# frontal + HIPPOCAMPUS (2 regions).  Port of _run_3C_states_MAST_GSEA.R.
# -----------------------------------------------------------------------------
# For each region we prerank the full-expressed Micro-PVM genes by signed
# winsorized avg_log2FC (NHD vs CON; +delta = higher in NHD) and run fgsea against
# the 19 curated canonical microglia-state signatures.  NES>0 = state genes at the
# NHD-up end; NES<0 = CON-up end.  Signature/DIRECTION-level enrichment of our NHD
# effect sizes against published state labels.
#
# atlas-free: reads only the FH cached full-coverage Cliff's-delta CSVs (which
#   carry avg_log2FC), tables/mast_dual/cliffs_delta_full_Micro-PVM_<region>.csv.
#   Regenerates the 04_states_MAST_GSEA cache for FH.  If a cache is missing stop.
#
# gene-set re-validation:
#   The 19 curated canonical state sets are carried VERBATIM (published labels are
#   the point of the panel; genes are intersected with the FH ranked list by
#   fgsea).  FH testable-set sizes were re-checked:
#     * DAM/MGnD/ARM all have 9-22 testable genes/region and prerank NHD-up — but
#       the enrichment is carried by the DAM-1 / inflammatory members (FTH1, CTSB,
#       B2M, MHC), not the DAM-2 lipid endpoint: LPL/CST7/LGALS3/CD9 are below the
#       FH detection floor and absent from the ranked list; GPNMB Hippo-only; SPP1
#       flat Frontal.  So a positive DAM/MGnD NES here = DAM-1 arrest, not DAM-2.
#     * MHC-II strongly NHD-up; homeostatic CON-up — TREM2-independent inflammation.
#   Nothing is dropped from the curated sets; sets with <15 testable genes are
#   flagged LOW-CONFIDENCE (the honest universal state with a ~2-4k ranked list),
#   <3 testable -> "too few genes to test".
#
# fgsea params: gseaParam=1, maxSize=500, set.seed(42), nPermSimple=1000,
#   minSize=3.  padj = BH within each region.
#
# Hippo Micro-PVM is a single NHD lane -> within-region GSEA is descriptive only.
#
# Outputs:
#   diagnostics/10_microglia_core_signature/04_states_MAST_GSEA_FH.csv
#   figures/Figure_2/panels/F2g_states_MAST_GSEA.{png,pdf}
#   provenance: .../04_states_MAST_GSEA_FH.log
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
  library(scales); library(forcats); library(ggh4x); library(ragg)
  library(fgsea)
})
set.seed(42)
msg <- function(...) cat(sprintf(...), "\n")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))

OUT_DIAG  <- file.path(PROJ, "diagnostics", "10_microglia_core_signature")
PANEL     <- file.path(PROJ, "figures", "Figure_2", "panels")
LOGS      <- file.path(PROJ, "logs")
for (d in c(OUT_DIAG, PANEL, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

MAST_DIR  <- file.path(PROJ, "tables", "mast_dual")
ARTIFACT  <- file.path(PROJ, "scripts", "_artifact_genes.R")
THEME     <- file.path(PROJ, "scripts", "22_publication_theme_FH.R")
GSEA_MET  <- file.path(PROJ, "scripts", "_gsea_composite_metric.R")
CURSIG    <- file.path(PROJ, "data", "curated_signatures.rds")

stopifnot("MISSING _artifact_genes.R"           = file.exists(ARTIFACT),
          "MISSING 22_publication_theme_FH.R"    = file.exists(THEME),
          "MISSING _gsea_composite_metric.R"     = file.exists(GSEA_MET),
          "MISSING curated_signatures.rds"       = file.exists(CURSIG))
source(ARTIFACT)   # is_artifact()
source(THEME)      # theme_pub, PAL_DGE, PAL_REGION, REGION_FULL, REGION_ORDER
source(GSEA_MET)   # fc_rank, GSEA_METRIC_NAME

REGIONS  <- REGION_ORDER            # c("Frontal","Hippo")
CT      <- "Micro-PVM"   # POOLED microglia+PVM; PVM = 5.4% of compartment
MINSIZE  <- 3
MAXSIZE  <- 500
LOWCONF_N<- 15
TOOFEW_N <- 3
GSEA_PARAM <- 1
NPERM      <- 1000

# --- Cliff's-delta caches (the preranking metric) ----------------------------
cache_path <- function(r) file.path(MAST_DIR,
  sprintf("cliffs_delta_full_%s_%s.csv", CT, r))
missing <- REGIONS[!file.exists(vapply(REGIONS, cache_path, character(1)))]
if (length(missing))
  stop(sprintf("MISSING Cliff's-delta cache(s) for %s — run scripts/47_full_expressed_ranking_FH.R first.",
               paste(missing, collapse = ", ")))
msg("All %d Micro-PVM Cliff's-delta caches present (atlas not loaded).", length(REGIONS))

# -----------------------------------------------------------------------------
# 1. GENE_SETS — the 19 curated canonical sets (inlined verbatim from
#    diagnostics/10_microglia_core_signature/01_microglia_core_signature_ORA.R so
#    the FH rebuild is self-contained; ZHOU_NHD_WR read from curated_signatures.rds).
# -----------------------------------------------------------------------------
DAM_KS2017 <- toupper(c(
  "Trem2","Tyrobp","Apoe","Ctsb","Ctsd","Ctsz","Cst7","Lpl",
  "Cd9","Itgax","Spp1","Clec7a","Axl","Fth1","Lyz2","B2m",
  "H2-K1","Gpnmb","Lilrb4","Csf1","Igf1","Tnfrsf1b","Ccl6",
  "Ctsa","Npc2","Lgals3","Lamp1","Ccl3","Cst3","Mmp12",
  "Cxcr4","Cadm1","Itgb5","Ch25h","Gusb","Ank"))
DAM_KS2017 <- gsub("H2-K1","HLA-A",DAM_KS2017)
MGND_KRASEMANN <- c("APOE","CLEC7A","AXL","CD9","CST7","GPNMB","LPL","SPP1",
                    "ITGAX","TREM2","TYROBP","IGF1","LGALS3","CTSD","CTSB","CCL3")
ARM_SALA <- c("APOE","CST7","SPP1","LPL","CLEC7A","ITGAX","GPNMB","CD9",
              "CTSD","CTSB","CTSZ","TYROBP","TREM2","AXL","CCL6","GUSB","FTH1")
TREM2_LOCKED_HM <- c("P2RY12","TMEM119","CX3CR1","SELPLG","CSF1R","SALL1","MEF2C","SIGLECH",
                     "HEXB","FCRLS","OLFML3","MERTK","TGFBR1","TGFBR2","SPARC")
MIC1_MATHYS <- c("CD163","F13A1","LYVE1","PLD4","FCN1","FTL","RNASE6",
                 "SLC2A5","TREM2","C1QA","C1QB","C1QC","APOE","TYROBP")
ZHOU_NHD_WR <- readRDS(CURSIG)$Zhou_NHD_microglia_up
stopifnot(length(ZHOU_NHD_WR) > 0)
OLAH_C4 <- c("APOE","SPP1","LPL","FABP5","ITGAX","CD9","GPNMB","CTSB","CTSD","LAMP1",
             "PSAP","CHI3L1","CD68","CD163")
OLAH_C7 <- c("HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DQA2","HLA-DQB1",
             "CD74","HLA-A","HLA-B","HLA-C","B2M","CIITA","STAT1","IRF1","IFITM3")
HM <- c("P2RY12","TMEM119","CX3CR1","CSF1R","HEXB","SALL1","SELPLG","OLFML3","MERTK",
        "MEF2C","SIGLEC8","P2RY13","CX3CL1")
TREM2_AXIS <- c("TREM2","TYROBP","HCST","SYK","PLCG2","VAV1","VAV3",
                "PIK3CG","PIK3CD","PIK3AP1","CSF1R","CARD9",
                "NFATC1","NFATC2","PTPN6","PTPN11","LCP2")
MHC1 <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G",
          "B2M","TAP1","TAP2","TAPBP","PSMB8","PSMB9","PSMB10","ERAP1","ERAP2")
MHC2 <- c("HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5",
          "HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DPA1","HLA-DPB1",
          "CD74","CIITA","CTSS","IFI30","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB")
COMPL <- c("C1QA","C1QB","C1QC","C3","C4A","C4B","C5AR1","C5AR2","CFH","CFI","CFB","C3AR1")
LYSO <- c("CTSB","CTSD","CTSS","CTSZ","CTSL","CTSH","CTSA","CTSK",
          "LAMP1","LAMP2","LIPA","NPC1","NPC2","GRN","PSAP","HEXA","HEXB","GBA1",
          "GLA","GNS","IDS","GAA","MAN2B1","GALC","PPT1","TPP1","CLN3","ASAH1","ATP6V0A1")
LIPID <- c("APOE","LPL","PLIN2","ACSL1","NCEH1","CD36","FABP5","FABP3",
           "SOAT1","ACAT1","SREBF1","SREBF2","HMGCR","LDLR","ABCA1","ABCA7",
           "APOC1","SCD","DHCR24","DHCR7","MSMO1","SQLE","MVK","FDFT1")
PHAGO <- c("AXL","MERTK","TYRO3","GAS6","PROS1","MEGF10","STAB1","C1QA","C1QB","C1QC",
           "SCARB1","CD36","TREM2","ITGAV","ITGB5","ITGB3","CD14","TLR2","TLR4")
SENESC <- c("CDKN1A","CDKN2A","FERMT2","B2M","LMNB1","H2AFX","GLB1","IGFBP3","SERPINE1",
            "IL6","TNF","CXCL8")
IFN <- c("ISG15","IFIT1","IFIT2","IFIT3","IFITM3","MX1","OAS1","OAS2","OAS3","IRF1","IRF7",
         "STAT1","STAT2","B2M","HLA-A","HLA-B","HLA-C")
ACUTE_INFL <- c("IL1B","IL6","TNF","CXCL8","CXCL10","CCL2","CCL3","CCL4","NFKB1","NFKB2",
                "NLRP3","CASP1","CASP4","GBP1","GBP5","NOS2")
GENE_SETS <- list(
  "DAM (Keren-Shaul 2017)"             = DAM_KS2017,
  "MGnD (Krasemann 2017)"              = MGND_KRASEMANN,
  "ARM (Sala Frigerio 2019)"           = ARM_SALA,
  "Mic1 AD human (Mathys 2019)"        = MIC1_MATHYS,
  "Zhou NHD wound-repair"              = ZHOU_NHD_WR,
  "Olah c4 lipid-activated human"      = OLAH_C4,
  "Olah c7 HLA/IFN human"              = OLAH_C7,
  "TREM2-null homeostatic (Filipello)" = TREM2_LOCKED_HM,
  "Homeostatic microglia"              = HM,
  "TREM2-DAP12 signaling axis"         = TREM2_AXIS,
  "MHC class I"                        = MHC1,
  "MHC class II"                       = MHC2,
  "Complement"                         = COMPL,
  "Lysosomal"                          = LYSO,
  "Lipid / cholesterol"                = LIPID,
  "Phagocytosis / efferocytosis"       = PHAGO,
  "Senescence / dystrophic"            = SENESC,
  "Interferon / type-I IFN"            = IFN,
  "Acute inflammation (IL1/TNF/NFkB)"  = ACUTE_INFL)
stopifnot(length(GENE_SETS) == 19)
msg("Loaded %d curated canonical gene sets.", length(GENE_SETS))

ALL_SETS <- lapply(GENE_SETS, function(g) unique(as.character(g)))
stopifnot(!any(duplicated(names(ALL_SETS))))
GRP_CANON <- "Canonical states"
set_block <- setNames(rep(GRP_CANON, length(GENE_SETS)), names(ALL_SETS))

# Display labels (ASCII-only; base pdf-safe).
bio_rename <- c(
  "DAM (Keren-Shaul 2017)"             = "DAM (Keren-Shaul 2017)",
  "MGnD (Krasemann 2017)"              = "MGnD (Krasemann 2017)",
  "ARM (Sala Frigerio 2019)"           = "ARM (Sala Frigerio 2019)",
  "Mic1 AD human (Mathys 2019)"        = "Mic1 AD (Mathys 2019)",
  "Zhou NHD wound-repair"              = "Zhou NHD w-rep (Zhou 2023)",
  "Olah c4 lipid-activated human"      = "Olah c4 (Olah 2020, lipid)",
  "Olah c7 HLA/IFN human"              = "Olah c7 (Olah 2020, HLA/IFN)",
  "TREM2-null homeostatic (Filipello)" = "TREM2-null homeo (Filipello 2018)",
  "Homeostatic microglia"              = "Homeostatic (Hickman 2013)",
  "TREM2-DAP12 signaling axis"         = "TREM2-DAP12 (canonical axis)",
  "MHC class I"                        = "MHC-I (HLA-A/B/C)",
  "MHC class II"                       = "MHC-II (HLA-DR/DP/DQ)",
  "Complement"                         = "Complement (C1Q/C3/CR)",
  "Lysosomal"                          = "Lysosomal (CTSB/CTSD/LAMP)",
  # Biology 4: do not label this unqualified "lipid" — the
  # NHD enrichment here is APOE/SREBF cholesterol-handling (DAM-1), not the DAM-2
  # lipid-droplet endpoint (LPL/CST7/LGALS3 are below the FH detection floor).
  "Lipid / cholesterol"                = "Cholesterol handling (APOE/SREBF; DAM-1)",
  "Phagocytosis / efferocytosis"       = "Phagocytosis (GO:BP)",
  "Senescence / dystrophic"            = "Senescence (Saul 2022, SenMayo)",
  "Interferon / type-I IFN"            = "Interferon type-I (ISG)",
  "Acute inflammation (IL1/TNF/NFkB)"  = "Inflammation (IL1/TNF/NFkB)")
stopifnot(all(names(ALL_SETS) %in% names(bio_rename)))

# -----------------------------------------------------------------------------
# 2. fgsea per region on signed winsorized avg_log2FC (BH within region).
# -----------------------------------------------------------------------------
run_region <- function(r) {
  cf <- cache_path(r)
  cd <- read.csv(cf, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot("expanded cache missing avg_log2FC — re-run scripts/47_full_expressed_ranking_FH.R" =
            all(c("gene","avg_log2FC") %in% names(cd)))
  cd <- cd[!is.na(cd$avg_log2FC), ]
  cd <- cd[!is_artifact(cd$gene, CT), ]
  cd <- cd[!duplicated(cd$gene), ]
  comp  <- fc_rank(cd$avg_log2FC)
  ranks <- sort(setNames(comp, cd$gene), decreasing = TRUE)
  ranks <- ranks[!is.na(ranks)]

  n_in_rank <- vapply(ALL_SETS, function(g) length(intersect(g, names(ranks))), integer(1))
  sets_r <- ALL_SETS
  set.seed(42)
  fg <- suppressWarnings(fgsea(pathways = sets_r, stats = ranks,
                               minSize = MINSIZE, maxSize = MAXSIZE,
                               gseaParam = GSEA_PARAM, nPermSimple = NPERM))
  fg <- as_tibble(fg) %>% select(state = pathway, NES, padj, size)

  out <- tibble(state = names(ALL_SETS)) %>%
    mutate(block  = unname(set_block[state]),
           region = r,
           n_genes_in_rank = n_in_rank[state]) %>%
    left_join(fg, by = "state") %>%
    mutate(
      tested    = !is.na(NES),
      low_conf  = n_genes_in_rank < LOWCONF_N,
      too_few   = n_genes_in_rank < TOOFEW_N,
      direction = case_when(!tested ~ NA_character_,
                            NES > 0 ~ "NHD-up",
                            TRUE    ~ "CON-up"))
  msg("  %-8s: ranked genes=%4d | sig FDR<0.05=%2d (NHD-up=%2d, CON-up=%2d) | untestable(<%d)=%2d",
      r, length(ranks),
      sum(out$tested & out$padj < 0.05, na.rm = TRUE),
      sum(out$tested & out$padj < 0.05 & out$NES > 0, na.rm = TRUE),
      sum(out$tested & out$padj < 0.05 & out$NES < 0, na.rm = TRUE),
      MINSIZE, sum(!out$tested))
  out
}
msg("\n== fgsea per region (signed winsorized avg_log2FC, minSize=%d, BH within region) ==", MINSIZE)
gsea <- bind_rows(lapply(REGIONS, run_region))

# -----------------------------------------------------------------------------
# 3. CSV output.
# -----------------------------------------------------------------------------
csv_out <- gsea %>%
  transmute(state, block, region, NES, padj, n_genes_in_rank, low_conf, direction)
csv_path <- file.path(OUT_DIAG, "04_states_MAST_GSEA_FH.csv")
write.csv(csv_out, csv_path, row.names = FALSE)
msg("\nWrote %s (%d rows)", basename(csv_path), nrow(csv_out))

# -----------------------------------------------------------------------------
# 4. Console report.
# -----------------------------------------------------------------------------
lab <- function(k) unname(bio_rename[k])
report_dir <- function(r, dir, n = 4) {
  d <- gsea %>% filter(region == r, tested, direction == dir) %>%
    arrange(padj, desc(abs(NES))) %>% head(n)
  if (!nrow(d)) return(sprintf("    %-7s : (none)\n", dir))
  paste(sprintf("    %-7s : %s\n", dir,
                paste(sprintf("%s (NES=%+.2f, padj=%.1e)", lab(d$state), d$NES, d$padj),
                      collapse = "; ")), collapse = "")
}
cat("\n=== Top states by NES + padj, per region ===\n")
for (r in REGIONS) { cat(sprintf("  %s:\n", r)); cat(report_dir(r, "NHD-up")); cat(report_dir(r, "CON-up")) }
cat("\n=== DAM/MGnD/ARM (positive NES = DAM-1 arrest; lipid endpoint absent from rank) ===\n")
dm <- gsea %>% filter(state %in% c("DAM (Keren-Shaul 2017)","MGnD (Krasemann 2017)","ARM (Sala Frigerio 2019)"))
for (i in seq_len(nrow(dm)))
  with(dm[i, ], cat(sprintf("  %-24s %-8s n=%d  %s\n", state, region, n_genes_in_rank,
    if (tested) sprintf("NES=%+.2f padj=%.2g dir=%s", NES, padj, direction) else "NES=NA")))

# -----------------------------------------------------------------------------
# 5. Panel — dotplot. NES fill = direction, size = -log10 padj, black outline =
#    sig, low-conf faded, untestable = open grey glyph.
# -----------------------------------------------------------------------------
DIR_COL <- c("NHD-up" = unname(PAL_DGE["Up in NHD"]),
             "CON-up" = unname(PAL_DGE["Up in CON"]))
set_meta <- gsea %>% group_by(state, block) %>%
  summarise(maxn   = max(n_genes_in_rank),
            minpad = suppressWarnings(min(padj, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(minpad = ifelse(is.finite(minpad), minpad, Inf),
         disp   = lab(state))
canon_order <- set_meta %>% filter(block == GRP_CANON) %>%
  arrange(desc(minpad)) %>% pull(state)             # bottom -> top (least sig first)
lev_keys   <- canon_order
lev_labels <- set_meta$disp[match(lev_keys, set_meta$state)]

plot_df <- gsea %>%
  left_join(set_meta %>% select(state, disp), by = "state") %>%
  mutate(set_label = factor(disp, levels = lev_labels),
         region    = factor(region, levels = REGIONS),
         # floor padj at 1e-10 so a p=0/underflow never becomes -log10 = Inf and
         # blows up the dot-size scale (reviewer-facing edge case).
         neglogp   = -log10(pmax(padj, 1e-10)),
         sig       = ifelse(tested & padj < 0.05, "sig", "ns"),
         conf      = ifelse(low_conf, "low", "well"))
tested_df <- plot_df %>% filter(tested)
untest_df <- plot_df %>% filter(!tested)

PANEL_TITLE <- NULL   # house rule: no in-panel title

pC <- ggplot(tested_df, aes(x = set_label, y = region)) +
  geom_point(aes(size = neglogp, fill = direction, colour = sig, alpha = conf),
             shape = 21, stroke = 0.55) +
  geom_point(data = untest_df, aes(x = set_label, y = region, shape = "too few to test"),
             colour = "grey70", fill = NA, size = 1.7, stroke = 0.4,
             inherit.aes = FALSE) +
  scale_fill_manual(values = DIR_COL, name = "NES direction",
                    breaks = c("NHD-up","CON-up")) +
  scale_colour_manual(values = c(sig = "black", ns = "grey65"),
                      breaks = c("sig","ns"), labels = c("FDR<0.05","n.s."),
                      name = NULL) +
  scale_alpha_manual(values = c(well = 1, low = 0.55),
                     breaks = c("well","low"), drop = FALSE,
                     labels = c(">=15 genes", "<15 genes (dir-level)"),
                     name = "set coverage") +
  scale_size_continuous(name = expression(-log[10]~FDR),
                        range = c(0.6, 4.6),
                        breaks = -log10(c(0.05, 0.01, 0.001)),
                        labels = c("0.05","0.01","0.001"),
                        limits = c(0, NA)) +
  scale_shape_manual(values = c("too few to test" = 1), name = NULL) +
  # y-axis annotation — house rule: the single-NHD-hippo-lane power caveat lives in
  # the legend, not on any panel.  Row label reads just the full region name.
  # vertical-centering fix: the render
  # height (H=2.9) is set by the tall 5-block right-hand legend, not by the 2 data
  # rows, so with the default discrete expansion the Frontal/Hippocampus rows sat
  # top-anchored in the upper third and the lower half of the box was empty.  A
  # symmetric additive expansion pads equal whitespace above and below the two rows
  # so they read as vertically centered in the box.  Dot sizes/positions unchanged.
  scale_y_discrete(labels = c(Frontal = unname(REGION_FULL["Frontal"]),
                              Hippo   = unname(REGION_FULL["Hippo"])),
                   limits = rev(REGIONS),
                   expand = expansion(add = c(1.9, 1.9))) +   # Frontal top, rows centered
  guides(fill   = guide_legend(override.aes = list(size = 3, alpha = 1,
                                colour = "grey25"), order = 1),
         colour = guide_legend(override.aes = list(size = 3, fill = "grey80"),
                                order = 2),
         alpha  = guide_legend(override.aes = list(size = 3, fill = "grey30",
                                colour = "grey30"), order = 3),
         size   = guide_legend(order = 4),
         shape  = guide_legend(override.aes = list(size = 2.4), order = 5)) +
  labs(x = NULL, y = NULL, title = PANEL_TITLE) +
  theme_pub(base_size = 8) +
  theme(axis.text.y       = element_text(size = 7.5, colour = "black"),
        axis.text.x       = element_text(size = 5.8, colour = "black", angle = 45,
                                         hjust = 1, margin = margin(t = 1)),
        axis.line         = element_blank(),
        panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.35),
        axis.ticks.x      = element_line(colour = "black", linewidth = 0.3),
        axis.ticks.y      = element_line(colour = "black", linewidth = 0.3),
        panel.grid.major  = element_line(colour = "grey92", linewidth = 0.25),
        panel.spacing.y   = unit(0.55, "lines"),
        plot.title        = element_blank(),
        legend.position   = "right",
        legend.key.size   = unit(0.30, "cm"),
        legend.title      = element_text(size = 6.4),
        legend.text       = element_text(size = 6.2),
        legend.spacing.y  = unit(0.05, "cm"),
        # The right-hand legend column stacks
        # five blocks (fill/colour/alpha/size/shape) with size=3 override keys, so it
        # is taller than the 1.9in panel and the topmost swatch (above "n.s.") was
        # clipped off the top canvas edge.  Widen the top plot.margin and the render
        # height so the full legend column fits inside the frame.
        plot.margin       = margin(10, 4, 2, 4))

# 2 region rows -> shorter.  Height bumped 1.9 -> 2.2 -> 2.9 so the tall
# right-hand legend column (five stacked blocks: fill/colour/alpha/size/shape) fits
# fully inside the frame; at 2.2 the top "NES direction" title + NHD-up swatch were
# still clipped off the top canvas edge.
W <- 6.58; H <- 2.9
png_path <- file.path(PANEL, "F2g_states_MAST_GSEA.png")
pdf_path <- file.path(PANEL, "F2g_states_MAST_GSEA.pdf")
ragg::agg_png(png_path, width = W, height = H, units = "in", res = 600)
print(pC); dev.off()
ggsave(pdf_path, pC, width = W, height = H, useDingbats = FALSE)
msg("\nWrote panel: %s\n         and: %s", png_path, pdf_path)

# -----------------------------------------------------------------------------
# 6. Provenance log.
# -----------------------------------------------------------------------------
log_path <- file.path(OUT_DIAG, "04_states_MAST_GSEA_FH.log")
writeLines(c(
  sprintf("run: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("PROJ: %s", PROJ),
  sprintf("metric: %s (signed avg_log2FC, winsorized at 1/99th pct; FC-only)", GSEA_METRIC_NAME),
  "pathway pre-filter: NONE (curated state signatures — all tested)",
  sprintf("params: minSize=%d maxSize=%d gseaParam=%g nPermSimple=%d seed=42",
          MINSIZE, MAXSIZE, GSEA_PARAM, NPERM),
  sprintf("sets: %d canonical", length(GENE_SETS)),
  "", "--- sessionInfo ---", capture.output(sessionInfo())), log_path)
msg("Wrote provenance: %s", log_path)

cat("\n=== DONE ===\n")
