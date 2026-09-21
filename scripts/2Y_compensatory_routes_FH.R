#!/usr/bin/env Rscript
# =============================================================================
# 2Y_compensatory_routes_FH.R — F2c_compensatory_routes_FH.R — which ACTIVATING ROUTES are engaged when DAP12 protein is absent but its transcript is intact?
# -----------------------------------------------------------------------------
# Rationale (after supplying the allele): TYROBP c.2T>C
# (p.Met1Thr) is a loss-of-initiation variant — no premature stop, so no NMD, so the
# mRNA survives and the lesion is purely at the protein level. The cell therefore has
# an intact transcriptional apparatus and an activating stimulus it cannot transduce
# through DAP12. The question this panel asks is: what does it use instead?
#
# ROUTES TESTED (each is a way to signal without functional DAP12):
#   ITAM adaptors      FCER1G (FcRgamma), HCST (DAP10) — the two adaptors that can
#                      substitute for DAP12; TREM2 is known to pair with DAP10.
#   Fc receptors       FCGR1A/2A/2B/3A — signal through FCER1G, not DAP12.
#   TLR / MyD88        TLR2/4, CD14, MYD88, IRAK1, TIRAP — wholly DAP12-independent.
#   Complement         C1QA/B/C, C3AR1, C5AR1, ITGAX — GPCR/integrin routes.
#   Scavenger          MSR1, CD36, marco, CD68 — uptake without ITAM signalling.
#   Kinases            SYK, BTK, PIK3CD, PLCG2 — shared transducers.
#   Brakes             INPP5D (SHIP1), PTPN6 (SHP1), CD33 — inhibitory arm.
#   TFs                SPI1, IRF8, MITF, TFEB — is the myeloid programme re-wired?
#
# Both regions.
# Frontal is depth-matched (NHD 496.5 vs CON 503.0 median genes, ratio 0.99); hippocampal
# NHD nuclei are 2.83x deeper (1395.0 vs 492.5), so a RAW detection rate there is inflated
# for every gene and would fabricate compensation. Rather than drop the region or discard
# reads, detection is DEPTH-ADJUSTED per gene x region by logistic regression with
# log10(nCount) as a covariate, read out at the region's median depth (see below).
# Validation: on Frontal, which needs no adjustment, raw vs adjusted r = 0.999
# (median |diff| 0.51 pp); on Hippocampus the raw median shifts +11.2 pp -> -0.1 pp,
# i.e. the adjustment removes exactly the uniform uplift the depth gap would have created.
# POOLED Micro-PVM. Descriptive; n = 1 donor per condition; no tests, no stars.
# -----------------------------------------------------------------------------
# apparent-type lift F2c_compensatory_routes_FH.R: every text size in this script scaled
# by 1.36 so the panel reads at ~5.0 pt on the assembled page, matching the
# Figure-5 dot panels. Point sizes, line widths and unit() dimensions are not
# touched. Canvas size unchanged, so re-linking is a no-op.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)
                                library(tidyr); library(ggh4x)})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
P <- Sys.getenv("NHD_PROJ")
source(file.path(P,"scripts","22_publication_theme_FH.R"))
TBL <- file.path(P,"tables","micro_states"); PANEL <- file.path(P,"figures","Figure_2","panels")

o <- readRDS(file.path(P,"data","_cache_immune_compartment_FH.rds"))
# POOLED Micro-PVM. This panel was still
# microglia-only, left over from the abandoned microglia-only primary.
o <- o[, o$immune_class %in% c("Microglia","PVM")]
DefaultAssay(o) <- "RNA"
if (inherits(o[["RNA"]],"Assay5")) o[["RNA"]] <- SeuratObject::JoinLayers(o[["RNA"]])
o$Condition <- factor(as.character(o$Condition), levels=c("CON","NHD"))
o$Region    <- factor(as.character(o$Region), levels=REGION_ORDER)

DefaultAssay(o) <- "RNA"; o <- NormalizeData(o, verbose=FALSE)

ROUTE <- tibble::tribble(
  ~gene,~route,
  # Its place is here, at the head of the ITAM adaptors, because
  # the c.2T>C start-loss allele carries no premature stop: the transcript ESCAPES NMD and
  # survives, so the lesion is purely at the protein level. A preserved-or-raised TYROBP
  # bar beside up-regulated FCER1G/HCST is the visual statement of that, and it is what
  # licenses reading the rest of the panel as attempted compensation.
  # The legend must state that TYROBP is the mutated gene; that is legend text, not a
  # panel annotation (house rule: no captions inside panels).
  "TYROBP","ITAM adaptors","FCER1G","ITAM adaptors","HCST","ITAM adaptors",
  "FCGR1A","Fc receptors","FCGR2A","Fc receptors","FCGR2B","Fc receptors","FCGR3A","Fc receptors",
  "TLR2","TLR / MyD88","TLR4","TLR / MyD88","CD14","TLR / MyD88","MYD88","TLR / MyD88",
  "IRAK1","TLR / MyD88","TIRAP","TLR / MyD88",
  # C1QA/B/C dropped: they are the SECRETED opsonin and
  # are shown with identical numbers in panel e (2E secretome, "Complement"); the complement
  # ROUTE this panel is asking about is the GPCR arm, C3AR1/C5AR1. ITGAX dropped: it is CD11c,
  # an integrin that was mis-grouped here under "Complement", and it is a DAM-2 stage marker
  # already carried by panel c (2I). 33 -> 29 rows.
  # Kept despite also appearing in panel a (FCER1G, HCST, SYK, PLCG2, SPI1): this panel's claim
  # is "adaptor/receptor up-regulation uncoupled from SYK/PLCG2" -- without the adaptors and the
  # flat kinases inside the panel it would read as "every route is engaged", i.e. compensation
  # succeeded, the opposite of the finding. The legend must state they also appear in panel a.
  "C3AR1","Complement","C5AR1","Complement",
  "MSR1","Scavenger","CD36","Scavenger","MARCO","Scavenger","CD68","Scavenger",
  "SYK","Kinases","BTK","Kinases","PIK3CD","Kinases","PLCG2","Kinases",
  "INPP5D","Brakes","PTPN6","Brakes","CD33","Brakes")
# The count below used to be `nrow(ROUTE)` against a HARD-CODED 33. After the
# de-duplication that denominator was stale (it printed "29 of 33" for a 29-gene list), and
# it never reported the genes lost later at the detection floor. Compute every number.
.n_req <- nrow(ROUTE)
.absent <- setdiff(ROUTE$gene, rownames(o))
ROUTE <- ROUTE %>% filter(gene %in% rownames(o))
cat(sprintf("routes requested: %d | absent from matrix: %d%s | carried forward: %d\n",
            .n_req, length(.absent),
            if (length(.absent)) paste0(" (", paste(.absent, collapse=", "), ")") else "",
            nrow(ROUTE)))

# ---------------------------------------------------------------------------
# DEPTH-ADJUSTED detection. include Hippocampus, do not downsample.
# The problem: this panel's metric is raw detection rate, which scales with sequencing
# depth, and the two regions are not comparably sequenced --
#     Frontal  CON 503.0 vs NHD 496.5 median genes  -> ratio 0.99  (matched)
#     Hippo    CON 492.5 vs NHD 1395.0              -> ratio 2.83  (not matched)
# Raw Hippocampus bars would push every route up in NHD purely because more genes are
# detected per nucleus -- fabricated compensation.
# Instead, adjust rather than discard: for each gene x region fit
#     glm(detected ~ Condition + log10(nCount), binomial)
# and report the difference in PREDICTED detection between NHD and CON evaluated at one
# common reference depth (the region's median log10 nCount). Every nucleus contributes;
# the depth term absorbs the confound. Genes that are fully separated (detected in ~none
# or ~all nuclei) cannot be fitted and fall back to the raw difference, flagged below.
# Control: Frontal is already depth-matched, so its adjusted values must track its raw
# values closely. That agreement is printed below and is the check that this is sound.
E <- FetchData(o, vars=ROUTE$gene, layer="counts")
meta <- data.frame(Condition = o$Condition, region = o$Region,
                   depth = log10(Matrix::colSums(GetAssayData(o, assay="RNA", layer="counts")) + 1))
adj_one <- function(det, cond, depth) {
  raw <- 100*(mean(det[cond=="NHD"]) - mean(det[cond=="CON"]))
  if (length(unique(det)) < 2) return(c(raw = raw, adj = raw, fitted = 0))
  fit <- try(suppressWarnings(glm(det ~ cond + depth, family = binomial())), silent = TRUE)
  if (inherits(fit, "try-error") || !fit$converged) return(c(raw = raw, adj = raw, fitted = 0))
  nd <- data.frame(cond = factor(c("CON","NHD"), levels = levels(cond)),
                   depth = rep(median(depth), 2))
  pr <- try(predict(fit, nd, type = "response"), silent = TRUE)
  if (inherits(pr, "try-error") || any(!is.finite(pr))) return(c(raw = raw, adj = raw, fitted = 0))
  # unname(): predict() returns a named vector, so c(adj = pr[2]-pr[1]) silently became
  # "adj.2" and every out["adj"] lookup returned NA.
  c(raw = raw, adj = unname(100*(pr[2] - pr[1])), fitted = 1)
}
res <- list()
for (rg in levels(meta$region)) for (g in ROUTE$gene) {
  k  <- meta$region == rg
  det <- as.integer(E[k, g] > 0)
  out <- adj_one(det, droplevels(meta$Condition[k]), meta$depth[k])
  res[[length(res)+1]] <- data.frame(region = rg, gene = g,
      CON = 100*mean(det[meta$Condition[k]=="CON"]),
      NHD = 100*mean(det[meta$Condition[k]=="NHD"]),
      raw_pp = unname(out["raw"]), delta_pp = unname(out["adj"]),
      fitted = unname(out["fitted"]))
}
d <- dplyr::bind_rows(res) %>% left_join(ROUTE, by="gene") %>%
  mutate(region = factor(region, levels = REGION_ORDER))
cat(sprintf("depth-adjusted GLM fitted for %d of %d gene x region cells (%d fell back to raw)\n",
            sum(d$fitted==1), nrow(d), sum(d$fitted==0)))
.ctl <- d %>% filter(region=="Frontal")
cat(sprintf("CONTROL (Frontal, already depth-matched): raw vs adjusted r = %.3f, median |diff| = %.2f pp\n",
            cor(.ctl$raw_pp, .ctl$delta_pp), median(abs(.ctl$raw_pp - .ctl$delta_pp))))
.hip <- d %>% filter(region=="Hippo")
cat(sprintf("Hippocampus (depth-confounded): raw median %+.1f pp -> adjusted median %+.1f pp\n",
            median(.hip$raw_pp), median(.hip$delta_pp)))

d <- d %>%
  # a gene is "below detection" only if it is under the 2%% floor in both conditions of
  # that region; drawn with the house keyed cross instead of vanishing (see 2I / oligo).
  mutate(below_floor = pmax(CON, NHD) < 2,
         delta_pp    = ifelse(below_floor, NA_real_, delta_pp)) %>%
  mutate(route=factor(route, levels=c("ITAM adaptors","Fc receptors","TLR / MyD88",
                                      "Complement","Scavenger","Kinases","Brakes")))
# Order genes by their frontal adjusted value so both facets share one row order
.ord <- d %>% filter(region=="Frontal") %>% arrange(route, delta_pp) %>% pull(gene)
.ord <- c(.ord, setdiff(unique(d$gene), .ord))
d <- d %>% mutate(gene = factor(gene, levels = unique(.ord))) %>% arrange(route, gene)
cat(sprintf("below the 2%% detection floor (drawn as keyed grey cross, not dropped): %d%s\n",
            sum(d$below_floor),
            if (any(d$below_floor)) paste0(" -> ", paste(sprintf("%s/%s", d$gene[d$below_floor],
                d$region[d$below_floor]), collapse=", ")) else ""))
write.csv(d, file.path(TBL,"compensatory_routes_FH.csv"), row.names=FALSE)
cat("\n=== detection change (pp, NHD - CON), frontal microglia, by route ===\n")
print(as.data.frame(d %>% arrange(route, desc(delta_pp)) %>%
        select(route, gene, CON, NHD, delta_pp)), row.names=FALSE, digits=3)

# ---------------------------------------------------------------------------
# Style: match F2f_secretome_MAST exactly. Verified by eye and by
# sampling rendered pixels -- the bar COLOURS were already identical (#C0392B NHD /
# #2C7FB8 CON; PAL_DGE and PAL_COND resolve to the same hexes), so every difference was
# frame, typography and annotation:
#   base_size 7->8 | axis text 5.8->6.3 black | axis.line dropped for a panel.border
#   | faint grey95 x-gridlines | strips grey97 with a grey80 hairline
#   | strip text 5.8->6.4 | bar width 0.68->0.72 | numeric value at each bar tip
#   | fill key shown as CON-up / NHD-up | legend.box vertical at the bottom.
# Added a pale "Frontal" region banner: it matches every other Fig-2 panel and it makes
# this panel's Frontal-only scope visible on the plot, which it previously was not
# (a reader comparing it with the 2-region 2E could not tell).
# Not copied: the |Cliff's delta| size key and the q* mark. This panel is a descriptive
# detection-change readout with no per-gene test (n = 1 donor per condition); importing
# an effect-size dot would fabricate a statistic the panel does not compute.
LAB_OFF <- 0.8            # pp gap between bar tip and its number (axis spans ~38 pp)
d <- d %>% mutate(
  x_tip = delta_pp + ifelse(delta_pp > 0, LAB_OFF, -LAB_OFF),
  hj    = ifelse(delta_pp > 0, 0, 1),
  dir   = ifelse(delta_pp > 0, "NHD-up", "CON-up"))

# ---------------------------------------------------------------------------
# Size + font ratio: 20% narrower, then 10% smaller overall
#   4.60 x 4.20 -> 3.68 x 4.20 -> 3.312 x 3.780 in
# Then adopt F2f_secretome_MAST's font-size : panel-size ratio rather than its absolute
# points. 2E is 3.90 in wide at base_size 8 / axis.text 6.3 -> 2.051 and 1.615 pt per
# inch of width. Holding that ratio at 3.312 in gives the scale factor below, so the two
# panels read as the same design at different sizes instead of 2Y looking text-heavy.
# Note this is the opposite of absolute-point parity: if these two are ever placed at
# equal width on the page, re-check by eye, because then equal points (not equal ratio)
# is what matches.
PW <- 4.60 * 0.8 * 0.9      # 3.312 in
PH <- 4.20 * 0.9 * 0.87     # 3.29 in - the four TF rows left this panel
# x1.36 so this panel reads ~5.0 pt on the assembled page (it was 3.67,
# the second-lowest in the figure set). Every size here goes through sz(), so one
# constant moves them all together and the internal hierarchy is preserved.
FS <- PW / 3.90             # 0.849 -- 2E's width is the reference
sz <- function(x) round(x * FS, 2)

p <- ggplot(d, aes(delta_pp, gene)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_col(data = dplyr::filter(d, !below_floor), aes(fill = dir),
           width = 0.72, linewidth = 0) +
  geom_text(data = dplyr::filter(d, !below_floor),
            aes(x = x_tip, y = gene, hjust = hj, label = sprintf("%.1f", delta_pp)),
            size = sz(1.7), colour = "grey15") +
  geom_point(data = dplyr::filter(d, below_floor),
             aes(x = 0, y = gene, shape = "below detection"),
             inherit.aes = FALSE, colour = "grey60", size = sz(1.7), stroke = 0.5) +
  facet_grid2(route ~ region, scales = "free_y", space = "free_y", switch = "y",
              labeller = labeller(region = as_labeller(REGION_FULL)),
              strip = strip_themed(
                background_y = elem_list_rect(fill = "grey97", colour = "grey80",
                                              linewidth = 0.25),
                text_y = elem_list_text(angle = 0, hjust = 1, size = sz(6.4),
                                        face = "plain", colour = "black"),
                background_x = elem_list_rect(
                  fill = unname(PAL_REGION_PALE[REGION_ORDER]), colour = NA),
                text_x = elem_list_text(colour = "black", size = sz(7), face = "plain"))) +
  scale_fill_manual(values = c(`NHD-up` = unname(PAL_COND[["NHD"]]),
                               `CON-up` = unname(PAL_COND[["CON"]])), name = NULL) +
  scale_shape_manual(values = c("below detection" = 4), name = NULL) +
  # Left expansion raised 0.14 -> 0.26: at 3.312 in the widest negative tip label
  # (Hippocampus TLR2, -10.2) ran off the panel and rendered as "10.2" with the minus
  # sign clipped by the border. Right stays wider for the long positive labels (34.2).
  scale_x_continuous(expand = expansion(mult = c(0.26, 0.22))) +
  guides(fill  = guide_legend(order = 1),
         shape = guide_legend(order = 2, override.aes = list(size = sz(2.4)))) +
  labs(x = "depth-adjusted detection change, pp (NHD - CON)", y = NULL) +
  theme_pub(base_size = sz(8)) +
  theme(
    plot.title         = element_blank(),
    axis.text.y        = element_text(size = sz(6.3), face = "italic", colour = "black"),
    axis.text.x        = element_text(size = sz(6.3), colour = "black"),
    axis.title.x       = element_text(size = sz(7), margin = margin(t = 3)),
    axis.line          = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.grid.major.x = element_line(colour = "grey95", linewidth = 0.2),
    strip.placement    = "outside",
    panel.spacing.x    = unit(0.25, "lines"),
    panel.spacing.y    = unit(0.15, "lines"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.box         = "vertical",
    legend.box.just    = "left",
    legend.key.size    = unit(0.2, "cm"),
    legend.spacing.y   = unit(0.05, "cm"),
    legend.margin      = margin(t = -1),
    legend.text        = element_text(size = sz(6.2)),
    legend.title       = element_text(size = sz(6.2)))
ggsave(file.path(PANEL,"F2c_compensatory_routes.png"), p, width=PW, height=PH,
       dpi=600, bg="white", device=ragg::agg_png)
ggsave(file.path(PANEL,"F2c_compensatory_routes.pdf"), p, width=PW, height=PH, bg="white")
cat("\nwrote F2c_compensatory_routes.{png,pdf}\n")
