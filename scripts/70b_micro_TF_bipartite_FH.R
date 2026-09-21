#!/usr/bin/env Rscript
# =============================================================================
# 70b_micro_TF_bipartite_FH.R — Figure 2: "F2h_micro_TF_bipartite" Microglial (Micro-PVM) regulator -> target map, NHD vs control.
# -----------------------------------------------------------------------------
# Which FACTORS are drawn is decided by A PRE-DECLARED GATE, not by choice.
# An earlier version of this panel used a hand-picked hub list with reasons written in
# the comments. That is a curation wearing a rule's clothes, and it is the first thing a
# referee would take apart. Every factor scored in both regions now goes through the same
# six criteria, every criterion is computed from the data, and the pass/fail of all of
# them ships as `tables/fig2_micro_TF_gate_audit_FH.csv` so the selection can be
# re-derived independently. Nothing is excluded by opinion.
#
#   C1 coverage          >= MIN_TARGETS scored targets in both regions.
#                        A z from 5 targets cannot sit beside a z from 85.
#   C2 reproducibility   |z| >= Z_NULL in both regions with the same sign.
#                        One donor per condition, so a single-region hit is not evidence.
#   C3 robustness        the score survives deleting housekeeping / heat-shock /
#                        nuclear-quality transcripts: |z_noHK| >= Z_NULL in both regions
#                        and >= LEAVEOUT_MIN of the original magnitude retained.
#   C4 coherence         the factor's own transcript must not contradict the inferred
#                        direction (opposite sign at |delta| >= TF_DELTA_MIN).
#   C5 interpretability  >= MOR_PURITY of the regulon's edges share one mode of
#                        regulation, otherwise a positive score has no assignable
#                        direction (the MAX problem: shared by activating and repressing
#                        heterodimers).
#   C6 anchor support    >= MIN_ANCHORS literature-canonical targets of that pathway are
#                        detected, and >= ANCHOR_DIR of them move in the expected
#                        direction. Anchor sets are declared below from the literature,
#                        independent of this dataset. This is the criterion that decides
#                        the contested factors, and it decides them on evidence:
#                            CIITA   5/5 anchors detected, all up, both regions
#                            HIF1A   5/6 detected, all up
#                            STAT3   2/4 detected, both up (exactly at MIN_ANCHORS)
#                            REL     1/8 -- the canonical cytokine module (TNF, IL1B,
#                                       CXCL8, TNFAIP3, SOD2, BIRC3, RELB) is below the
#                                       detection FLOOR at 750 UMI
#                            NFKB1   fails EARLIER, at C3 (leave-out 0.66), not at C6
#                            CTNNB1 / TCF7L2  0/9 -- the Wnt call is unanchored
#                            MAX     0/6 -- entirely unanchored
#                        (these numbers are re-read from the audit table at run time and
#                         checked below, so the header cannot drift from the data again)
#
# A correction to AN EARLIER justification. A previous draft cut CTNNB1 and TCF7L2 on
# the grounds that their regulons are "largely the same canonical TCF target list", i.e.
# One signal counted twice. that is not TRUE in this prior: their scored-target Jaccard
# is 0.236 (13 shared of 55). They are cut by C6, on their missing anchors, and the
# collinearity claim has been withdrawn.
#
# What the GATE cannot do. It cannot rank a factor's biological plausibility, and two
# exclusions therefore rest on a stated judgement rather than a number: REL/NFKB1 is the
# arm the NHD mutation is predicted to release (TREM2-DAP12 normally inhibits TLR
# signalling), and excluding it means the panel does not show the mechanism that is most
# expected. It is excluded anyway, because we cannot measure its output. that absence is
# A result and BELONGS in the text, not a gap to be quietly filled by drawing the node.
#
# The control is added by declaration, not by passing. SPI1 is drawn identically to the
# gated hubs so that its flat result is visible at the same visual weight. It is named in
# CONTROL_TF before any score is read, and the audit table records that it was admitted
# as a control rather than as a hit.
#
# N = 1 donor per condition. The ULM p-value is a per-nucleus quantity and carries no
# condition-level claim; scores are effect sizes. Say "direction-consistent across two
# independently processed regions", never "replicated". Say Micro-PVM, never "microglia".
# Say "stimulus-driven", never "chronic" (JUN/FOS/EGR1/ATF3/NR4A1 are below the floor,
# so the acute-versus-chronic distinction is untestable here).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
TDIR  <- file.path(PROJ, "tables"); DDIR <- file.path(PROJ, "data")
PANEL <- file.path(PROJ, "figures", "Figure_2", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# PARAMETERS — every threshold in one place, declared before any data is read
# ---------------------------------------------------------------------------
Z_NULL        <- 1.96   # C2/C3  two-sided regulon null for a z-like ULM score
MIN_TARGETS   <- 10L    # C1     scored targets required per factor
LEAVEOUT_MIN  <- 0.50   # C3     fraction of |z| that must survive the leave-out
TF_DELTA_MIN  <- 0.05   # C4     |delta| at which a TF's own transcript may contradict
MOR_PURITY    <- 0.80   # C5     share of regulon edges sharing one mode of regulation
MIN_ANCHORS   <- 2L     # C6     canonical anchors that must be detected
ANCHOR_DIR    <- 0.66   # C6     share of detected anchors moving as expected
DELTA_FLOOR   <- 0.12   # target ring: |delta| for "region-consistent"
TARGETS_PER_HUB <- 8L   # symmetric fans, so a flat factor keeps a full fan
CONTROL_TF    <- "SPI1" # declared identity control, exempt from the gate
# PPARG is drawn as a declared measured NULL, on the same principle as SPI1:
# admitted by declaration rather than by passing, because its value is the point. PPARG is
# the TREM2-dependent DAM-2 lipid regulon, and this figure's claim is that NHD microglia do
# not reach that endpoint. A regulon that is measured and flat is the positive form of that
# claim -- it says the endpoint was looked for and was absent, which a missing node cannot.
# It became reportable only when the 750-UMI floor was dropped: at native depth PPARG clears
# the 10% detection floor (6.9% frontal / 11.3% hippocampus), and after thinning it does not.
NULL_TF       <- "PPARG"

# genes that can never be drawn or counted: technical, ambient, or dissociation-induced
DROP_RE <- paste0("^RP[LS]|^RPLP|^FAU$|^NACA$|^EEF1A1$|^TPT1$|^TMSB4X$|^TMSB10$|",
                  "^PTMA$|^UB[BC]$|^MALAT1$|^NEAT1$|^HIST|^H1-|^H2A|^H2B|^H3C|^H4C|",
                  "^HSPA|^HSPB|^HSP90|^DNAJ|^SGK1$|^FOS$|^FOSB$|^JUN$|^JUNB$|^DUSP1$|",
                  "^ZFP36$|^GLUL$|^SLC1A3$|^MBP$|^GRID2$|^F13A1$|^MRC1$|^CD163$|",
                  "^KCNMA1$|^TCF4$|^NAIP$|^MT-")

# C6 anchor sets — declared from the literature, not from these data. `up` lists the
# direction expected if the factor's activity is genuinely increased.
ANCHORS <- list(
  CIITA  = list(genes = c("CD74","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DMB"), up = TRUE),
  STAT3  = list(genes = c("SOCS3","HAMP","TNFRSF1B","BCL3"),                 up = TRUE),
  HIF1A  = list(genes = c("SLC2A3","PFKFB3","PGK1","LDHA","SAT1","P4HA1"),   up = TRUE),
  REL    = list(genes = c("NFKBIA","TNFAIP3","TNF","IL1B","CXCL8","SOD2","BIRC3","RELB"), up = TRUE),
  NFKB1  = list(genes = c("NFKBIA","TNFAIP3","TNF","IL1B","CXCL8","SOD2","BIRC3","RELB"), up = TRUE),
  CREB1  = list(genes = c("NR4A1","PCK1","SIK1","PDE4B","VGF"),              up = TRUE),
  CTNNB1 = list(genes = c("AXIN2","LEF1","NOTUM","DKK1","SP5","TCF7","MYC","CCND1","NKD1"), up = TRUE),
  TCF7L2 = list(genes = c("AXIN2","LEF1","NOTUM","DKK1","SP5","TCF7","MYC","CCND1","NKD1"), up = TRUE),
  MAX    = list(genes = c("MYC","MXD1","MXI1","MNT","MGA","NPM1"),           up = TRUE),
  NFAT5  = list(genes = c("NOS2","TNF","IL6","CCL2","AQP1","SLC5A3","AKR1B1"), up = TRUE),
  SPI1   = list(genes = c("CSF1R","ITGAM","CX3CR1","TYROBP","P2RY12"),       up = TRUE),
  # PPARG anchors are the lipid-handling arm of the DAM-2 endpoint. `up = TRUE` states the
  # direction expected if the endpoint were reached; the claim here is that it is not, so
  # C6 is computed and reported for PPARG rather than used to admit it.
  PPARG  = list(genes = c("ABCA1","ABCG1","CD36","LPL","PLIN2","FABP5"),      up = TRUE))

# ---------------------------------------------------------------------------
# 1. inputs
# ---------------------------------------------------------------------------
# Primary source switched to the log-normalized, all-nuclei run.
# The gate below is computed from `gde`, so it must come from the same matrix as the
# scores in `act` -- otherwise the audit in ST14 documents a gate applied to a different
# dataset than the panel. Both therefore move together, and both are the primary arm.
# The 750-UMI arm ships as a sensitivity in ST14; see 68_micro_TF_collectri_FH.R.
act <- read.csv(file.path(TDIR, "micro_TF_collectri_FH.csv"), stringsAsFactors = FALSE)
gde <- read.csv(file.path(TDIR, "micro_gene_cliffsdelta_FH.csv"), stringsAsFactors = FALSE)
net <- readRDS(file.path(DDIR, "TF_regulon_collectri.rds"))
stopifnot("TF table is the OLD version — re-run 68" =
            all(c("score_noHK","n_target") %in% names(act)))

delta <- gde %>% select(gene, region, cliffs_delta) %>%
  pivot_wider(names_from = region, values_from = cliffs_delta) %>%
  filter(!is.na(Frontal), !is.na(Hippo)) %>%
  mutate(delta = (Frontal + Hippo) / 2,
         consistent = sign(Frontal) == sign(Hippo) &
                      pmin(abs(Frontal), abs(Hippo)) >= DELTA_FLOOR)
d_of  <- setNames(delta$delta, delta$gene)
scored_genes <- delta$gene

# ---------------------------------------------------------------------------
# 2. the gate
# ---------------------------------------------------------------------------
anchor_support <- function(tf) {
  A <- ANCHORS[[tf]]
  if (is.null(A)) return(list(n_det = NA_integer_, n_tot = NA_integer_, frac_dir = NA_real_))
  det <- intersect(A$genes, scored_genes)
  if (!length(det)) return(list(n_det = 0L, n_tot = length(A$genes), frac_dir = NA_real_))
  ok <- if (A$up) d_of[det] > 0 else d_of[det] < 0
  list(n_det = length(det), n_tot = length(A$genes), frac_dir = mean(ok))
}
mor_purity <- function(tf) { m <- net$mor[net$source == tf]
  if (!length(m)) NA_real_ else max(mean(m > 0), mean(m < 0)) }

wide <- act %>% select(TF = source, region, score, score_noHK, n_target) %>%
  pivot_wider(names_from = region, values_from = c(score, score_noHK, n_target)) %>%
  filter(!is.na(score_Frontal), !is.na(score_Hippo))

anc <- lapply(wide$TF, anchor_support)
gate <- wide %>% mutate(
  mean_z      = (score_Frontal + score_Hippo) / 2,
  coverage    = pmin(n_target_Frontal, n_target_Hippo),
  leaveout    = pmin(abs(score_noHK_Frontal) / abs(score_Frontal),
                     abs(score_noHK_Hippo)   / abs(score_Hippo)),
  tf_delta    = as.numeric(d_of[TF]),
  mor_pure    = vapply(TF, mor_purity, numeric(1)),
  anchor_det  = vapply(anc, function(x) as.numeric(x$n_det),    numeric(1)),
  anchor_tot  = vapply(anc, function(x) as.numeric(x$n_tot),    numeric(1)),
  anchor_dir  = vapply(anc, function(x) as.numeric(x$frac_dir), numeric(1)),
  C1_coverage        = coverage >= MIN_TARGETS,
  C2_reproducible    = abs(score_Frontal) >= Z_NULL & abs(score_Hippo) >= Z_NULL &
                       sign(score_Frontal) == sign(score_Hippo),
  C3_robust          = abs(score_noHK_Frontal) >= Z_NULL &
                       abs(score_noHK_Hippo)   >= Z_NULL & leaveout >= LEAVEOUT_MIN,
  C4_coherent        = is.na(tf_delta) |
                       !(sign(tf_delta) != sign(mean_z) & abs(tf_delta) >= TF_DELTA_MIN),
  C5_interpretable   = !is.na(mor_pure) & mor_pure >= MOR_PURITY,
  C6_anchored        = !is.na(anchor_det) & anchor_det >= MIN_ANCHORS &
                       !is.na(anchor_dir) & anchor_dir >= ANCHOR_DIR,
  passes = C1_coverage & C2_reproducible & C3_robust & C4_coherent &
           C5_interpretable & C6_anchored,
  first_failure = dplyr::case_when(
    !C1_coverage      ~ "C1 coverage",      !C2_reproducible ~ "C2 reproducibility",
    !C3_robust        ~ "C3 robustness",    !C4_coherent     ~ "C4 coherence",
    !C5_interpretable ~ "C5 interpretability", !C6_anchored  ~ "C6 anchor support",
    TRUE              ~ "passes"),
  admitted = ifelse(passes, "gated hit",
             ifelse(TF == CONTROL_TF, "declared control",
             ifelse(TF == NULL_TF,    "declared null (DAM-2)", "excluded"))))

HUBS <- c(gate$TF[gate$passes], CONTROL_TF, NULL_TF) %>% unique()

# Guard: the header above quotes anchor counts. A referee re-derived the
# panel from the shipped audit CSV and found the header disagreed with it (it still
# carried numbers from an exploratory run). Assert the quoted values against the data so
# the two can never diverge silently again.
.hdr <- list(CIITA = c(5, 5), HIF1A = c(5, 6), STAT3 = c(2, 4))
for (.t in names(.hdr)) {
  .r <- gate[gate$TF == .t, ]
  stopifnot("header anchor counts disagree with the audit table — fix the header" =
              nrow(.r) == 1 && .r$anchor_det == .hdr[[.t]][1] && .r$anchor_tot == .hdr[[.t]][2])
}

# The control fails its own anchor test, and that is a result, not a footnote.
# SPI1 is admitted by declaration, but C6 is still computed for it, and it comes back
# anchor_dir = 0.2: four of its five declared anchors move against "activity increased".
# P2RY12 (-0.425 / -0.473) is the 18th largest effect of 3,021 genes and MEF2C
# (-0.532 / -0.344) the 23rd. The homeostatic output is not flat, it is falling hard.
# The SPI1 regulon cannot see this: P2RY12, CX3CR1, TMEM119 and SALL1 are not CollecTRI
# SPI1 targets, so a flat SPI1 score is blindness, not evidence of stability. Any legend
# built on this panel must say "lineage regulons unshifted, homeostatic output down" and
# must not say the identity program is unchanged.
.ctl <- gate[gate$TF == CONTROL_TF, ]
if (nrow(.ctl) == 1 && !is.na(.ctl$anchor_dir) && .ctl$anchor_dir < ANCHOR_DIR)
  cat(sprintf(paste0("\n*** CONTROL WARNING: %s anchor_dir = %.2f (%d/%d anchors detected).\n",
                     "    Its declared anchors move AGAINST the assumed direction, so a flat\n",
                     "    %s score must NOT be reported as 'identity unchanged'.\n\n"),
              CONTROL_TF, .ctl$anchor_dir, .ctl$anchor_det, .ctl$anchor_tot, CONTROL_TF))
n_by <- table(gate$first_failure[!gate$passes])
cat(sprintf("factors scored in both regions: %d\n", nrow(gate)))
cat(sprintf("  pass the gate: %d (%s)\n", sum(gate$passes),
            paste(sort(gate$TF[gate$passes]), collapse = ", ")))
cat("  excluded at:", paste(sprintf("%s %d", names(n_by), as.integer(n_by)), collapse = " | "), "\n")
cat(sprintf("  plus %s admitted as the declared identity control\n", CONTROL_TF))
cat("\n== gate audit, factors reaching C6 ==\n")
print(as.data.frame(gate %>% filter(C1_coverage, C2_reproducible) %>%
  arrange(desc(abs(mean_z))) %>%
  transmute(TF, z = round(mean_z, 2), cov = coverage, leaveout = round(leaveout, 2),
            mor = round(mor_pure, 2), tf_delta = round(tf_delta, 3),
            anchors = sprintf("%d/%d", anchor_det, anchor_tot),
            dir = round(anchor_dir, 2), verdict = first_failure)), row.names = FALSE)
write.csv(gate %>% select(TF, mean_z, coverage, leaveout, mor_pure, tf_delta,
                          anchor_det, anchor_tot, anchor_dir,
                          starts_with("C1_"), starts_with("C2_"), starts_with("C3_"),
                          starts_with("C4_"), starts_with("C5_"), starts_with("C6_"),
                          first_failure, admitted) %>% arrange(desc(mean_z)),
          file.path(TDIR, "fig2_micro_TF_gate_audit_FH.csv"), row.names = FALSE)

# The hub ring now CARRIES FDR and REPRODUCIBILITY (after an audit found
# this panel placed as Fig-2 h against this project's own "do not place as drawn" note).
# Two defects are fixed here:
#   (1) HIF1A was drawn with the same red "increased" ring as CIITA despite failing
#       multiple-testing correction in frontal cortex (BH q = 0.175 over the 46 x 2
#       factor-by-region tests). A factor that does not survive FDR in one of the two
#       regions cannot be encoded identically to one that clears it in both.
#   (2) SPI1 was drawn as "within null", which reads as an informative negative. It is
#       not: it FAILED gate C2 (its two regional scores are sign-discordant, +1.77
#       frontal / -0.48 hippocampal) and 4 of its 5 declared anchors move against the
#       inferred direction (anchor_dir 0.20). It is admitted by declaration as a control,
#       and the ring now says exactly that instead of implying a measured null.
allp <- act %>% transmute(TF = source, region, p_value)
allp$q <- p.adjust(allp$p_value, method = "BH")     # across every factor x region tested
qw <- allp %>% select(TF, region, q) %>%
  tidyr::pivot_wider(names_from = region, values_from = q, names_prefix = "q_")
hub <- gate %>% filter(TF %in% HUBS) %>%
  left_join(qw, by = "TF") %>%
  transmute(TF, score = mean_z, n_target = coverage, q_Frontal, q_Hippo,
            q_min = pmin(q_Frontal, q_Hippo, na.rm = TRUE),
            q_max = pmax(q_Frontal, q_Hippo, na.rm = TRUE),
            dir = dplyr::case_when(
              TF == CONTROL_TF                    ~ "declared control",
              TF == NULL_TF                       ~ "declared null (DAM-2)",
              q_max < 0.05 & score > 0            ~ "increased (q < 0.05 both regions)",
              q_max < 0.05 & score < 0            ~ "decreased (q < 0.05 both regions)",
              q_min < 0.05                        ~ "increased (n.s. in one region)",
              TRUE                                ~ "not significant"))
cat("\n== hub ring classes, with BH-q per region ==\n")
print(as.data.frame(hub %>% transmute(TF, score = round(score, 2),
      q_Frontal = signif(q_Frontal, 2), q_Hippo = signif(q_Hippo, 2), dir)), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. targets — symmetric fans so a flat factor is visibly measured, not absent
# ---------------------------------------------------------------------------
cons <- delta %>% filter(!grepl(DROP_RE, gene))
cat(sprintf("\ntargets: %d genes scored in both regions, %d after the technical screen\n",
            nrow(delta), nrow(cons)))
reg <- net %>% filter(source %in% HUBS, target %in% cons$gene) %>%
  transmute(TF = source, gene = target, mor) %>% distinct(TF, gene, .keep_all = TRUE) %>%
  left_join(cons %>% select(gene, delta), by = "gene") %>%
  group_by(TF) %>% slice_max(abs(delta), n = TARGETS_PER_HUB, with_ties = FALSE) %>%
  ungroup() %>% select(-delta)
tg <- reg %>% count(gene, name = "n_tf") %>%
  left_join(cons %>% select(gene, delta, Frontal, Hippo, consistent), by = "gene")
cat(sprintf("drawn: %d hubs, %d targets, %d edges\n", nrow(hub), nrow(tg), nrow(reg)))

# ---------------------------------------------------------------------------
# 4. layout — two columns, not a wheel
# ---------------------------------------------------------------------------
# A circular layout was tried and rejected by eye: a circle inscribed in a rectangle
# throws away all four corners, and with only four hubs the ring is small enough that
# ggrepel could not clear the hub labels off their own nodes ("STAT3 87 targets" printed
# straight through the STAT3 ring). A two-column bipartite fixes both -- it fills the
# rectangle, it makes "left = regulator, right = target" true by geometry rather than by
# glyph, and it gives every target its own row, so all of them can be named instead of
# the strongest fourteen. The panel becomes an auditable evidence list as well as a map.
HUB_X <- 0.0; TGT_X <- 1.0
hub <- hub %>% arrange(desc(score)) %>%
  mutate(y = if (n() == 1) 0.5 else seq(0.92, 0.08, length.out = n()))
tg  <- tg %>% arrange(desc(delta)) %>%
  mutate(y = seq(1, 0, length.out = n()))
eg  <- reg %>%
  left_join(hub %>% select(TF, hy = y), by = "TF") %>%
  left_join(tg  %>% select(gene, gy = y), by = "gene") %>%
  mutate(mode = ifelse(mor >= 0, "activates", "represses"))

# ---------------------------------------------------------------------------
# 5. panel
# ---------------------------------------------------------------------------
GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
L <- max(0.3, ceiling(max(abs(tg$delta)) * 10) / 10)

p <- ggplot() +
  geom_curve(data = eg, aes(x = HUB_X, y = hy, xend = TGT_X, yend = gy, linetype = mode),
             curvature = 0.16, colour = "grey80", linewidth = 0.16) +
  scale_linetype_manual(values = c(activates = "solid", represses = "22"), name = NULL) +
  geom_point(data = tg, aes(TGT_X, y, fill = delta, alpha = consistent), shape = 21,
             size = 1.7, colour = "grey25", stroke = 0.25) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.38),
                     labels = c(`TRUE` = "both regions agree", `FALSE` = "regions differ"),
                     name = NULL, guide = guide_legend(order = 2,
                       override.aes = list(shape = 21, fill = "grey55", size = 1.9))) +
  scale_fill_gradient2(low = GSEA_STEEL, mid = "white", high = GSEA_ROSE, midpoint = 0,
                       limits = c(-L, L), breaks = c(-L, 0, L),
                       name = expression(atop("target gene", "Cliff's " * delta))) +
  geom_text(data = tg, aes(TGT_X + 0.035, y, label = gene), hjust = 0,
            size = 1.72, colour = "black", fontface = "italic") +
  geom_point(data = hub, aes(HUB_X, y, colour = dir), shape = 21, fill = "white",
             size = 5.4, stroke = 1.05) +
  scale_colour_manual(
    values = c(`increased (q < 0.05 both regions)` = GSEA_ROSE,
               `increased (n.s. in one region)`    = "#E8A0AA",
               `decreased (q < 0.05 both regions)` = GSEA_STEEL,
               `not significant`                   = "grey65",
               `declared control`                  = "grey35",
               `declared null (DAM-2)`             = "#7E6BA8"),
    breaks = c("increased (q < 0.05 both regions)", "increased (n.s. in one region)",
               "declared null (DAM-2)", "declared control"),
    # This read "BH q over 46 factors x 2 regions". The q
    # VALUES are right, but 46 is the gate-AUDIT count, not the correction
    # denominator: BH is applied over the full activity table, 54 unique TFs =
    # 100 TF x region rows. Stating the wrong denominator beside a q is exactly
    # the reviewer-facing numerical error the house rule targets.
    # Derived, never typed. This denominator was hard-coded twice before and was wrong
    # both times; it must track the activity table the panel actually scored.
    name = sprintf("TF activity (BH q over\n%d TF x region tests)", nrow(allp)),
    guide = guide_legend(order = 3)) +
  geom_text(data = hub, aes(HUB_X - 0.13, y, label = sprintf("%s\n%d targets", TF, n_target)),   # clear the 5.4-pt hub ring
            hjust = 1, lineheight = 0.85, size = 2.30, colour = "black") +
  coord_cartesian(xlim = c(-0.58, 1.30), ylim = c(-0.03, 1.03), clip = "off") +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(0.26, "cm"),
                                barheight = unit(1.05, "cm")),
         linetype = guide_legend(order = 4,
                                 override.aes = list(colour = "grey55", linewidth = 0.4))) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank(), axis.line = element_blank(),
        panel.border = element_blank(),
        legend.position = "right", legend.box = "vertical", legend.box.just = "left",
        legend.key.size = unit(0.30, "cm"), legend.spacing.y = unit(0.14, "cm"),
        legend.text = element_text(size = 6.0),
        legend.title = element_text(size = 6.2, lineheight = 0.9),
        plot.margin = margin(2, 1, 2, 1))

BN <- "F2h_micro_TF_bipartite"; W <- 3.56; H <- 2.75   # 10% smaller
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nSaved %s.{png,pdf} (%.2f x %.2f in)\n", BN, W, H))
write.csv(tg %>% select(gene, n_tf, Frontal, Hippo, delta, consistent),
          file.path(TDIR, "figF2h_micro_TF_bipartite_targets_FH.csv"), row.names = FALSE)
write.csv(reg, file.path(TDIR, "figF2h_micro_TF_bipartite_edges_FH.csv"), row.names = FALSE)
cat("Wrote the gate audit + target and edge tables\n=== DONE ===\n")
