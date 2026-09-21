# =============================================================================
# 2G_micro_region_venn_FH.R — Figure 2 panel g: microglia regional-overlap Venn, frontal vs HIPPOCAMPUS (2 circles). Port of 54_fig2_micro_region_venn.R
# from the 3-circle (Frontal/OCC/Hippo) original.
# -----------------------------------------------------------------------------
# 2-way Venn of the Micro-PVM per-nucleus MAST discovery hits (padj<0.05 &
# |Cliff's d|>=0.15) between Frontal and Hippocampus, to show which microglial DGE
# genes are region-unique vs shared.  Atlas-free (reads the FH discovery CSV only).
#
# Manifest note: panel g is SUPP-RECOMMENDED (redundant with the volcano,
#   panel c) — the region-unique/shared split is already visible there.  Curated here
#   nonetheless so it is submission-ready if kept in the main figure.
#
# Gene LISTS (re-derived for FH here):
#   The displayed region-unique gene lists are data-driven from the FH discovery
#   table — the top microglia-credible unique genes by |Cliff's delta| per lobe,
#   with ambient/soup genes stripped via is_artifact() and lncRNA/contigs removed.
#   (The old hardcoded credible lists were OCC/old-atlas-keyed and are not ported.)
#
# Outputs (Figure_2/panels):
#   2_micro_region_venn_genes.{png,pdf}   fixed 2-circle Venn with gene lists
# Tables (tables/mast_dual):
#   micro_region_venn_membership_FH.csv   every gene x region-set x direction
#   micro_region_unique_genes_FH.csv      the region-EXCLUSIVE gene lists
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggforce); library(ragg)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # PAL_REGION, REGION_FULL, REGION_ORDER
source(file.path(PROJ, "scripts", "_artifact_genes.R"))           # is_artifact()

MAST  <- file.path(PROJ, "tables/mast_dual/MAST_dual_discovery_all.csv")
PANEL <- file.path(PROJ, "figures/Figure_2/panels")
TBL   <- file.path(PROJ, "tables/mast_dual")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TBL, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
CT      <- "Micro-PVM"   # POOLED

# --- load Microglia discovery hits ------------------------------------------
d <- read.csv(MAST, stringsAsFactors = FALSE) %>%
  filter(cell_type == CT, region %in% REGIONS,
         discovery %in% c(TRUE, "TRUE")) %>%
  mutate(dir = ifelse(avg_log2FC > 0, "NHD-up", "CON-up"),
         max_pct = pmax(pct.1, pct.2))
cat("Microglia discovery hits per region:\n"); print(table(d$region))

# The displayed region-unique callouts were neuronal-ambient soup (FRMD4A, PDE4D,
# KCNIP4/RBFOX1/CSMD1/NRG3/SYT1/LRRTM4/DLGAP1/OPCML) + lncRNAs (SNHG14, MEG3,
# NUTM2B-AS1, LINC01578) — not microglial.  Root causes fixed here:
#   (a) removed the MAX_LNC lncRNA re-injection slot (lncRNA are never re-added).
#   (b) a bare top-|delta| ranking surfaces globally-abundant neuronal ambient that
#       is not caught by is_artifact() (FRMD4A/PDE4D/OXR1/SORL1/PDE3B are not in the
#       neuronal list) -> add a positive-identity requirement: a displayed callout
#       must be a microglial/myeloid-credible gene (MICRO_IDENTITY below), and clear
#       a min-pct detection floor, and not be an explicitly suppressed survivor.
# The Venn COUNTS (289/173/215) are unchanged — this only curates the displayed
# labels.  Genuine microglial region-unique genes present in the FH discovery table
# (see tables/mast_dual/micro_region_unique_genes_FH.csv) now surface.
MICRO_IDENTITY <- c(
  # complement / Fc / phagocytic-receptor myeloid core
  "C1QA","C1QB","C1QC","C3","CD163","CD14","CD68","CD86","FCGR1A","FCGR2A","FCGR2B",
  "FCGR3A","FCGRT","MRC1","MSR1","MERTK","AXL","GAS6","TREM2","TYROBP","ITGAM","ITGAX",
  # antigen presentation / MHC
  "CD74","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DMA","HLA-DMB","HLA-DPA1","HLA-DPB1",
  "HLA-DQA1","HLA-DQB1","CIITA","B2M",
  # iron / lysosomal / lipid microglial effectors
  "FTL","FTH1","HAMP","SLC11A1","HMOX1","CTSB","CTSD","CTSS","LAMP1","GPNMB","APOE",
  "PSAP","GRN","LPL","CST7","SPP1","CD9",
  # inflammatory / signalling / TF microglial effectors
  "CEBPD","RGS1","SRGN","TNFRSF1B","NAIP","IL1B","TLR2","TLR5","TLR8","NFKB1","IRF8",
  "SPI1","MITF","MAP3K8","GPR183","PLAUR","SERPINE1","CHI3L1","F13A1",
  # microglial GTPase/adaptor/kinase identity genes seen region-unique in FH
  "SRGAP3","RAPGEF1","SIPA1L1","PEAK1","MGAT1","P2RY12","MEF2C","PLXDC2","SORL1",
  "IRAK3","HDAC9","PPARG","PPARD","VSIR","MAMDC2","FAM49A","ADGRG6","CD163L1")
# explicit suppression (belt-and-braces vs the allowlist): named neuronal/ambient
# survivors + lncRNAs the reviewer flagged, must never surface as callouts.
CALLOUT_SUPPRESS <- c("FRMD4A","PDE4D","OXR1","SORL1","PDE3B","MEF2C","PLXDC2",
                      "SNHG14","MEG3","NUTM2B-AS1","LINC01578")
PCT_FLOOR_LABEL <- 0.15   # min-pct detection floor for a displayed callout
# Per (gene, region) max_pct so the detection floor is applied in the same region
# whose exclusive lobe the gene sits in.
pct_gr <- d %>% distinct(gene, region, max_pct)

# --- membership + counts ----------------------------------------------------
LNC_RE    <- "^LINC[0-9]|-AS[0-9]*$|-DT$|^MIR[0-9].*HG$|^SNHG[0-9]|^MEG3$|^MALAT1$|^NEAT1$|^PVT1$|^XIST$|^JPX$|^KCNQ1OT1$|^CCDC26$"
CONTIG_RE <- "^(AC|AL|AP|FP|Z|BX|CR)[0-9]+\\.[0-9]+$"
memb <- d %>%
  mutate(is_lnc    = grepl(LNC_RE, gene, ignore.case = TRUE),
         is_contig = grepl(CONTIG_RE, gene),
         hard_amb  = is_artifact(gene, CT) & !grepl(LNC_RE, gene, ignore.case = TRUE)) %>%
  group_by(gene) %>%
  summarise(regions = paste(sort(unique(region)), collapse = "+"),
            nreg = dplyr::n_distinct(region),
            maxcd = max(abs(cliffs_delta), na.rm = TRUE),
            is_lnc = any(is_lnc), is_contig = any(is_contig),
            hard_amb = any(hard_amb), .groups = "drop")
combo_count <- table(memb$regions)
cc <- function(k) { v <- combo_count[k]; if (length(v)==0 || is.na(v)) 0L else as.integer(v) }
n_Frontal <- cc("Frontal"); n_Hippo <- cc("Hippo"); n_both <- cc("Frontal+Hippo")
cat(sprintf("Venn lobes: Frontal-only=%d, Hippo-only=%d, shared=%d\n",
            n_Frontal, n_Hippo, n_both))

# --- positive-identity-gated region-unique display lists (FH) -----------------
# A displayed callout must be: (1) exclusive to the lobe, (2) not ambient/contig/
# lncRNA, (3) not on CALLOUT_SUPPRESS, (4) on the MICRO_IDENTITY allowlist (positive
# microglial/myeloid identity), AND (5) clear the min-pct detection floor in this
# region.  no lncRNA reserve slot (MAX_LNC removed).  Ranked by |Cliff's delta|.
N_SHOW  <- 8L      # genes listed per exclusive lobe
lobe_display <- function(reg) {
  memb %>%
    left_join(pct_gr %>% filter(region == reg) %>% select(gene, max_pct),
              by = "gene") %>%
    filter(regions == reg, !hard_amb, !is_contig, !is_lnc,
           !gene %in% CALLOUT_SUPPRESS,
           gene %in% MICRO_IDENTITY,
           !is.na(max_pct), max_pct >= PCT_FLOOR_LABEL) %>%
    arrange(desc(maxcd)) %>% head(N_SHOW) %>% pull(gene)
}
FRONTAL_LIST <- lobe_display("Frontal")
HIPPO_LIST   <- lobe_display("Hippo")
cat("Frontal-unique microglial callouts (identity-gated, top):",
    paste(FRONTAL_LIST, collapse=", "), "\n")
cat("Hippo-unique microglial callouts (identity-gated, top):",
    paste(HIPPO_LIST,   collapse=", "), "\n")
stopifnot("no Frontal callouts survived the identity+pct gate" = length(FRONTAL_LIST) > 0,
          "no Hippo callouts survived the identity+pct gate"   = length(HIPPO_LIST)   > 0)

# --- fixed 2-circle geometry (overlapping vertically, Frontal on TOP) --------
# "make it vertical, frontal on top, and make sure we read well
# all in it". Previously the circles overlapped left-right and the flanking gene
# lists were anchored at x = -+2.35 — close enough that the longest names (MITF,
# TLR5, APOE / SIPA1L1, MGAT1, FAM49A) printed on TOP of the circle outlines.
# Now: Frontal circle above, Hippocampus below, and each list is anchored outside
# the circles' x-range (they span x = -+R = -+1.55, lists start at -+1.85 and grow
# away from the figure), so no label can touch a circle regardless of name length.
R   <- 1.55
cen <- data.frame(region = factor(REGIONS, levels = REGIONS),
                  x0 = c(0, 0), y0 = c(0.85, -0.85))       # Frontal TOP, Hippo below
# counts inside: Frontal-only (top) | shared (centre lens) | Hippo-only (bottom)
in_count <- data.frame(
  x = c(0, 0, 0), y = c(1.35, 0.00, -1.35),
  lab = c(as.character(n_Frontal), as.character(n_both), as.character(n_Hippo)))
# Gene lists outside the circles: Frontal upper-LEFT, Hippocampus lower-RIGHT
# (diagonal placement keeps each list beside its own lobe without stacking the
# panel taller than the circles).
out <- data.frame(region = factor(REGIONS, levels = REGIONS),
                  x = c(-1.85, 1.85), y = c(2.15, 0.55), hjust = c(1, 0))
out$genes <- c(paste(head(FRONTAL_LIST, N_SHOW), collapse = "\n"),
               paste(head(HIPPO_LIST,   N_SHOW), collapse = "\n"))
# invisible reserve box tight to content (circles + outside lists).
# The symmetric box therefore reserved ~1.3 units (~0.5 in) of dead white on the
# left. Box is now sized to the actual content on each side; some headroom is kept
# because shrinking the canvas makes fixed-pt text occupy more data units.
reserve <- data.frame(x = c(-3.25, 3.85), y = c(2.55, -2.55))

# apparent-type fix. This panel was saved 3.22 in wide but placed about
# 1.41 in wide in SuppFig3 -- a scale of 0.44, which drove the gene lists from 7.1 pt
# declared down to roughly 3.1 pt on the page, below the print floor. That is why they
# read as illegible; the declared size was never the problem, the placement scale was.
# The panel is therefore re-rendered at its intended placed width (2.55 in), so the
# scale is ~1.0 and declared pt == on-page pt, and the type is set to the house band
# (~5.1 pt body). Nothing about the panel's content changes: same circles, same counts,
# same gene lists. If the panel is placed at a width other than 2.55 in, these numbers
# must be re-derived.
VEN_GENE_MM  <- 2.50   # gene lists
VEN_LABEL_MM <- 2.60   # region names + counts (2.25 -> 2.60)
VEN_TITLE_PT <- 8.0    # panel title (7 -> 8)
p_genes <- ggplot() +
  geom_blank(data = reserve, aes(x = x, y = y)) +
  ggforce::geom_circle(data = cen,
    aes(x0 = x0, y0 = y0, r = R, fill = region, colour = region),
    alpha = 0.15, linewidth = 0.7) +
  scale_fill_manual(values = PAL_REGION, guide = "none") +
  scale_colour_manual(values = PAL_REGION, guide = "none") +
  geom_label(data = in_count, aes(x = x, y = y, label = lab),
             fontface = "plain", size = VEN_LABEL_MM, colour = "grey15", fill = "white", alpha = 0.55,
             label.size = 0, label.padding = unit(0.4, "mm")) +
  geom_text(data = out, aes(x = x, y = y, label = REGION_FULL[as.character(region)],
                            colour = region, hjust = hjust),
            fontface = "plain", size = VEN_LABEL_MM, show.legend = FALSE) +
  geom_text(data = out, aes(x = x, y = y - 0.30, label = genes, hjust = hjust),
            vjust = 1, size = VEN_GENE_MM, lineheight = 0.98, fontface = "italic", colour = "grey15") +
  coord_equal(clip = "off", expand = FALSE) +
  # The panel never said which cell type it was about. Name it, using
  # the paper-wide joint-class label "Micro-PVM" and the same title treatment as the lead
  # UMAP (F2a_micro_UMAP_composite), so the two panels agree. This is a TITLE (what the
  # panel contains), not caption text, so it does not breach the no-captions-in-panels rule.
  ggtitle("Micro-PVM") +
  theme_void(base_size = 8) +
  theme(plot.title = element_text(size = VEN_TITLE_PT, hjust = 0.5, face = "plain", colour = "black",
                                  margin = margin(t = 1, b = 2)),
        plot.margin = margin(2, 3, 2, 3))
# canvas matches the reserve-box aspect (7.10 x 5.10 units = 1.392) so coord_equal
# leaves no dead margin, and is 10% smaller than the first vertical version
# (height 2.57 -> 2.31; width follows from the aspect).
# +0.17 in of height for the new title so coord_equal does not shrink the circles to
# make room (the width stays on the reserve-box aspect).
# canvas keeps the reserve-box aspect (7.10/5.10) so coord_equal leaves no dead margin;
# k is chosen so the panel is saved at its intended placed width of 2.55 in.
VEN_K <- 2.55 / (7.10 / 5.10)
VW <- VEN_K * (7.10 / 5.10); VH <- VEN_K + 0.17
ggsave(file.path(PANEL, "2_micro_region_venn_genes.pdf"), p_genes, width = VW, height = VH)
agg_png(file.path(PANEL, "2_micro_region_venn_genes.png"), width = VW, height = VH,
        units = "in", res = 600); print(p_genes); invisible(dev.off())
cat("wrote 2_micro_region_venn_genes.{png,pdf}\n")

# --- membership table + region-EXCLUSIVE gene lists (deliverable) ------------
gene_region <- d %>%
  group_by(gene) %>%
  summarise(regions   = paste(sort(unique(region)), collapse = "+"),
            n_regions = n_distinct(region),
            dirs      = paste(sprintf("%s:%s", region, dir)[order(region)], collapse = "; "),
            dir_conflict = n_distinct(dir) > 1,
            max_abs_cliffs = max(abs(cliffs_delta), na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(n_regions), regions, desc(max_abs_cliffs))
write.csv(gene_region, file.path(TBL, "micro_region_venn_membership_FH.csv"), row.names = FALSE)

uniq <- gene_region %>% filter(n_regions == 1) %>%
  mutate(unique_to = regions) %>%
  select(unique_to, gene, dirs, max_abs_cliffs) %>%
  arrange(unique_to, desc(max_abs_cliffs))
write.csv(uniq, file.path(TBL, "micro_region_unique_genes_FH.csv"), row.names = FALSE)

cat("\n== region-set sizes (all discovery hits) ==\n")
print(gene_region %>% count(regions, name = "n_genes") %>% arrange(desc(n_genes)))
cat("\n== region-UNIQUE gene counts ==\n")
print(uniq %>% count(unique_to, name = "n_unique"))
cat("\n== shared-by-both core (n =", sum(gene_region$n_regions == 2), ") ==\n")
core <- gene_region %>% filter(n_regions == 2) %>% arrange(desc(max_abs_cliffs))
cat(paste(head(core$gene, 40), collapse = ", "), "\n")
if (any(gene_region$dir_conflict))
  cat("\nNOTE:", sum(gene_region$dir_conflict),
      "gene(s) change direction across regions (flagged dir_conflict in CSV):\n  ",
      paste(gene_region$gene[gene_region$dir_conflict], collapse = ", "), "\n")

cat("\n=== DONE ===\n", file = stderr())
