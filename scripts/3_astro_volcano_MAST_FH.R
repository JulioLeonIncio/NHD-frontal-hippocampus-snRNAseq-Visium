#!/usr/bin/env Rscript
# =============================================================================
# 3_astro_volcano_MAST_FH.R — F3b_astro_volcano_MAST_FH.R — Figure 3: per-nucleus MAST astrocyte volcano, frontal + HIPPOCAMPUS (2 regions). Port of _run_3D_astro_volcano_MAST.R.
# -----------------------------------------------------------------------------
# Sibling of 2C_micro_volcano_MAST_FH.R: identical construction (facet, two-tier
# size, threshold key, cap guard) — only cell type (Astro) + programme colours
# differ.  MAST-forward, atlas-free (reads MAST_dual_discovery_all.csv only).
#
#   x = avg_log2FC (NHD vs CON)   y = -log10(padj) [descriptive ranking axis only]
#   colour = astro pathway program; two-tier size (nominal padj<0.05 small /
#            strict + |Cliff's d|>=0.15 big).
# Honest FRAMING (project policy): y-axis padj is per-nucleus + pseudoreplication-
# inflated -> descriptive only.  No donor-level significance claim; the p=0.05
# h-guide is a threshold marker.  No in-panel title/caption.
#
# gene-set re-validation
# against the astro MAST full-expressed tables.  kept genes are detected + moving
# in their programme direction in >=1 region.  CHANGES vs the 3-region reference:
#   Metallothionein (lost): kept MT2A,MT3,MT1G,MT1E,MT1F,MT1M,MT1X,SLC30A1.
#       dropped SLC39A14 (up in Frontal, wrong direction), MT1H (only Hippo at the
#       0.10 pct floor).  Headline is MT2A (d-0.62/-0.43) + MT3 (brain-canonical).
#   Reactive/DAA (up): kept GFAP,CD44,TNC,OSMR,SERPINA3,EMP1,STAT3.
#       dropped C3,LCN2,TIMP1,S100A6,A2M (absent in FH astro), VIM (flat d~0),
#       CHI3L1 (sign-FLIPS Frontal+/Hippo-).
#   Heat-shock (up, Hippo-led): kept HSPA1A,HSPA1B,HSPB1,CRYAB,HSP90AA1,DNAJB1,
#       HSPH1,BAG3.  dropped DNAJA1,HSPA6 (absent).  (Frontal heat-shock is flat;
#       the programme is Hippocampus-specific — direction-match handles this.)
#   Homeostatic (lost): kept SLC1A2,GLUL,AQP4,ATP1A2,GJA1,NDRG2,GPC5.
#       dropped SLC1A3 (up in Frontal d+0.57, direction-inconsistent), GLIS3 (flat).
#
# Output: figures/Figure_3/panels/F3b_astro_volcano_MAST.{png,pdf} + provenance.
# -----------------------------------------------------------------------------
# Harmonised to Figure 2's house type spec (base 8 / axis 6.3 /
# title 7 / legend 6.2) and rendered at its placed size so no assembly down-scaling is
# needed. Fig 2's panels total ~117 in2 on a 67.9 in2 page, i.e. they are shrunk ~0.77x
# at assembly, which is what pushed its type under the 5 pt floor. Fig 3 avoids that.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggrepel); library(ggh4x); library(ragg)
  library(patchwork)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))

MAST_PATH <- file.path(PROJ, "tables", "mast_dual", "MAST_dual_discovery_all.csv")
ARTIFACT  <- file.path(PROJ, "scripts", "_artifact_genes.R")
THEME     <- file.path(PROJ, "scripts", "22_publication_theme_FH.R")
PANEL     <- file.path(PROJ, "figures", "Figure_3", "panels")
LOGS      <- file.path(PROJ, "logs")
for (f in c(MAST_PATH, ARTIFACT, THEME))
  if (!file.exists(f)) stop(sprintf("MISSING %s — prep step failed", f))
for (d in c(PANEL, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
source(ARTIFACT); source(THEME)

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
CT      <- "Astro"
YCAP    <- 60
FC_MIN  <- 1.20
LFC_MIN <- log2(FC_MIN)
CD_CUT   <- 0.15
PADJ_CUT <- 0.05

# --- themed astrocyte programs (FH-revalidated; see header) ------------------
ASTRO_THEME <- c(
  # Up: reactive / DAA
  GFAP="Reactive / DAA", CD44="Reactive / DAA", TNC="Reactive / DAA",
  OSMR="Reactive / DAA", SERPINA3="Reactive / DAA", EMP1="Reactive / DAA",
  STAT3="Reactive / DAA",
  # Up: heat-shock / proteostasis (Hippo-led)
  HSPA1A="Heat-shock / proteostasis", HSPA1B="Heat-shock / proteostasis",
  HSPB1="Heat-shock / proteostasis", CRYAB="Heat-shock / proteostasis",
  HSP90AA1="Heat-shock / proteostasis", DNAJB1="Heat-shock / proteostasis",
  HSPH1="Heat-shock / proteostasis", BAG3="Heat-shock / proteostasis",
  # Lost: metallothionein / metal-detox
  MT2A="Metallothionein (lost)", MT3="Metallothionein (lost)",
  MT1G="Metallothionein (lost)", MT1E="Metallothionein (lost)",
  MT1F="Metallothionein (lost)", MT1M="Metallothionein (lost)",
  MT1X="Metallothionein (lost)", SLC30A1="Metallothionein (lost)",
  # Lost: homeostatic glutamate / water transport
  SLC1A2="Homeostatic (lost)", GLUL="Homeostatic (lost)", AQP4="Homeostatic (lost)",
  ATP1A2="Homeostatic (lost)", GJA1="Homeostatic (lost)", NDRG2="Homeostatic (lost)",
  GPC5="Homeostatic (lost)")

THEME_UP   <- c("Reactive / DAA", "Heat-shock / proteostasis")
THEME_LOST <- c("Metallothionein (lost)", "Homeostatic (lost)")
THEME_LEVELS <- c(THEME_UP, THEME_LOST)
THEME_PAL <- c(
  "Reactive / DAA"            = "#C1272D",
  "Heat-shock / proteostasis" = "#D81B60",
  "Metallothionein (lost)"    = "#2E86AB",
  "Homeostatic (lost)"        = "#1B9E77")

# -----------------------------------------------------------------------------
cat("== load FH MAST discovery table ==\n")
raw <- read.csv(MAST_PATH, stringsAsFactors = FALSE)
need <- c("cell_type","region","gene","avg_log2FC","p_val_adj",
          "cliffs_delta","discovery","pct.1","pct.2")
stopifnot(all(need %in% colnames(raw)))

d <- raw %>%
  filter(cell_type == CT, region %in% REGIONS) %>%
  mutate(
    region    = factor(region, levels = REGIONS),
    discovery = as.logical(discovery),
    negl10    = -log10(pmax(p_val_adj, .Machine$double.xmin)),
    abs_cd    = abs(cliffs_delta),
    low_pct   = ifelse(avg_log2FC > 0, pct.2, pct.1),
    dir_col   = dplyr::case_when(
                  discovery & avg_log2FC >  LFC_MIN ~ "Up in NHD",
                  discovery & avg_log2FC < -LFC_MIN ~ "Up in CON",
                  TRUE                              ~ "ns"),
    artifact  = is_artifact(gene, CT),
    theme     = unname(ASTRO_THEME[gene]),
    theme_dir = dplyr::case_when(theme %in% THEME_UP  ~ "up",
                                 theme %in% THEME_LOST ~ "lost",
                                 TRUE                  ~ NA_character_),
    dir_match = !is.na(theme) &
                ((theme_dir == "up"   & avg_log2FC >  LFC_MIN) |
                 (theme_dir == "lost" & avg_log2FC < -LFC_MIN)),
    themed_strict = dir_match & discovery,
    themed        = dir_match & (p_val_adj < PADJ_CUT),
    tier      = factor(ifelse(themed_strict, "strict",
                       ifelse(themed, "nominal", "ns")),
                       levels = c("ns", "nominal", "strict")),
    theme_f   = factor(ifelse(themed, theme, NA), levels = THEME_LEVELS),
    gene_disp = gene)
stopifnot(nrow(d) > 0)

n_uf <- sum(d$p_val_adj <= 0 | d$p_val_adj < .Machine$double.xmin, na.rm = TRUE)
if (n_uf > 0)
  cat(sprintf("NOTE: %d padj values <= xmin floored to %.3e (-log10 capped, no Inf)\n",
              n_uf, .Machine$double.xmin))

n0 <- nrow(d); cat(sprintf("Astro rows across %d regions: %d\n", length(REGIONS), n0))
stopifnot(sum(d$discovery & is.na(d$abs_cd)) == 0)
d$y_disp <- pmin(d$negl10, YCAP); d$capped <- d$negl10 > YCAP

counts <- d %>% group_by(region) %>%
  summarise(NHD_up = sum(dir_col == "Up in NHD"), CON_up = sum(dir_col == "Up in CON"),
            theme_nominal = sum(tier == "nominal"), theme_strict = sum(tier == "strict"),
            ns = sum(dir_col == "ns"), capped = sum(capped), .groups = "drop")
cat("\n== discovery points + two-tier theme counts per region ==\n"); print(as.data.frame(counts))

lab_df <- d %>% filter(themed)
cat("\n== themed (coloured + labelled) genes per region ==\n")
for (r in REGIONS) {
  gg <- lab_df %>% filter(region == r) %>% arrange(theme_f, desc(abs(avg_log2FC)))
  cat(sprintf("  %-7s [%d]: %s\n", r, nrow(gg),
              paste(sprintf("%s(%s)", gg$gene_disp, abbreviate(gg$theme, 4)), collapse = ", ")))
}
cat("\n== program membership present in the volcano (any region) ==\n")
print(lab_df %>% distinct(gene_disp, theme) %>% count(theme, name = "n_genes") %>%
      arrange(match(theme, THEME_LEVELS)) %>% as.data.frame())

chk <- lab_df %>% filter(gene %in% c("MT2A","MT3","MT1G","GFAP","CD44","HSP90AA1","SLC1A2","GLUL")) %>%
  transmute(region, gene, avg_log2FC, cliffs_delta, dir = ifelse(avg_log2FC > 0, "NHD-up", "NHD-down"))
cat("\n== headline-gene direction check ==\n"); print(as.data.frame(chk))

# -----------------------------------------------------------------------------
# Build volcano (2-region facet).
# -----------------------------------------------------------------------------
d <- d[order(match(d$tier, c("ns","nominal","strict"))), ]
d$region <- factor(as.character(d$region), levels = REGIONS)
df_ns  <- dplyr::filter(d, !themed, dir_col == "ns")
df_d0  <- dplyr::filter(d, !themed, dir_col != "ns")
df_nom <- dplyr::filter(d, tier == "nominal")
df_str <- dplyr::filter(d, tier == "strict")
df_th  <- dplyr::filter(d, themed)
xl <- max(abs(d$avg_log2FC)) * 1.04

YCEIL <- YCAP * 1.22
reg_upper <- vapply(REGIONS, function(r) {
  yv <- d$y_disp[d$region == r & d$themed]
  tm <- if (length(yv)) as.numeric(quantile(yv, 0.90, names = FALSE))
        else max(d$y_disp[d$region == r], na.rm = TRUE)
  pinned <- any(d$y_disp[d$region == r & d$themed] >= YCAP - 1e-6)
  cap <- min(tm + max(2.5, 0.12 * tm), YCEIL)
  if (pinned) min(YCEIL, YCAP + 0.20 * YCAP) else cap
}, numeric(1))
cat(sprintf("\n== per-facet y upper crop (topmost label + headroom; ceiling %.1f) ==\n", YCEIL))
print(round(reg_upper, 1))
y_scales <- lapply(REGIONS, function(r)
  scale_y_continuous(limits = c(0, reg_upper[[r]]), oob = scales::squish,
                     breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0.02, 0.03))))

# Matching the fix applied to Fig 2's microglial volcano. Above ~11
# labels per facet ggrepel's force_pull drags them back together whatever the start
# positions, and max.overlaps=Inf means a non-converged solve DRAWS the overlap rather than
# dropping a label. Selection order is unchanged (strict tier, then headline genes) and
# theme colour still encodes every point that loses its text.
# Fig 2's volcano converges at 11 because its facets are
# ~3.3 in wide; this panel is 4.00 in for two facets (~2.0 in each), so the same label count
# still overlapped (MT1F/MT2A, MT1M/MT1E, ATP1A2/GFAP). Label capacity scales with facet
# width, not with the gene list.
# The remaining collisions (MT1M/MT1E, SLC1A2/GFAP) are all
# in the y-CAP band at 60, where many genes pin at the same height and their labels compete
# for one horizontal strip -- extra labels there cost more than they add. Six per facet
# converges; programme colour still encodes every themed point.
# Panel widened 4.00 -> 4.80 (+20%). I tried raising the cap 6 -> 8 on the
# assumption that label capacity scales with facet width -- it does not here, and 8 brought
# back MT1M/MT1E and SLC1A2/GFAP. Those pairs pin in the y-CAP band at 60, where every capped
# gene shares one horizontal strip; extra width gives them nowhere new to go, only extra
# height would. Cap stays at 6, verified clean; the extra width goes to the point cloud.
LAB_CAP  <- 6L
HEADLINE <- c("MT2A","MT3","MT1G","MT1E","GFAP","CD44","HSP90AA1","HSPA1A","SLC1A2","GLUL","AQP4")
df_lab <- df_th %>% group_by(region) %>%
  arrange(desc(tier == "strict"), desc(gene %in% HEADLINE), desc(abs(avg_log2FC)), .by_group = TRUE) %>%
  slice_head(n = LAB_CAP) %>% ungroup() %>% as.data.frame()
cat(sprintf("\n== labels drawn per facet (cap %d of themed) ==\n", LAB_CAP))
print(df_lab %>% count(region, name = "n_labels") %>% as.data.frame())
# halve the outward nudge: at 0.10*xl it pinned right-hand labels against the panel edge
# where they had nowhere left to resolve.
df_lab$nudge_x <- sign(df_lab$avg_log2FC) * 0.05 * xl

# Deterministic y PRE-STAGGER (same algorithm as Fig 2's volcano). ggrepel starts every
# label at its point, so labels at nearly the same y on the same side start stacked and the
# solver can settle in that overlapping local minimum. Walk each region x side group in y
# order, separate any pair closer than SEP, then fit the result into the panel: shrink SEP to
# what the axis can hold and re-centre on the group's original centre. Pushing blindly in one
# direction jams the ceiling; a uniform slide jams the floor.
df_lab$nudge_y <- 0
for (r in unique(df_lab$region)) for (sd in c(-1, 1)) {
  ii <- which(df_lab$region == r & sign(df_lab$avg_log2FC) == sd)
  if (length(ii) < 2) next
  top <- reg_upper[[as.character(r)]] * 0.97; bot <- 0.02 * top
  ii  <- ii[order(df_lab$y_disp[ii])]
  yp  <- df_lab$y_disp[ii]; n_l <- length(yp)
  SEP <- min(0.055 * top, (top - bot) / max(n_l - 1, 1))
  for (k in seq_len(n_l)[-1]) if (yp[k] - yp[k-1] < SEP) yp[k] <- yp[k-1] + SEP
  ctr <- mean(range(df_lab$y_disp[ii]))
  yp  <- yp - (mean(range(yp)) - ctr)
  if (max(yp) > top) yp <- yp - (max(yp) - top)
  if (min(yp) < bot) yp <- yp + (bot - min(yp))
  df_lab$nudge_y[ii] <- yp - df_lab$y_disp[ii]
}
cat(sprintf("label pre-stagger: %d of %d shifted in y\n",
            sum(abs(df_lab$nudge_y) > 1e-9), nrow(df_lab)))

p <- ggplot() +
  geom_vline(xintercept = c(-LFC_MIN, LFC_MIN), linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey60", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25, linetype = "22") +
  geom_point(data = df_ns, aes(avg_log2FC, y_disp), colour = "#949494", size = 0.5, alpha = 0.5, shape = 16, stroke = 0) +
  geom_point(data = df_d0, aes(avg_log2FC, y_disp), colour = "#ADADAD", size = 0.75, alpha = 0.45, shape = 16, stroke = 0) +
  geom_point(data = df_nom, aes(avg_log2FC, y_disp, colour = theme_f), size = 0.9, alpha = 0.5, shape = 16, stroke = 0) +
  geom_point(data = df_str, aes(avg_log2FC, y_disp, colour = theme_f), size = 2.3, alpha = 0.90, shape = 16, stroke = 0) +
  geom_text_repel(data = df_lab, aes(avg_log2FC, y_disp, label = gene_disp),
                  colour = "black", fontface = "italic", size = 2.0, seed = 42,
                  nudge_x = df_lab$nudge_x, nudge_y = df_lab$nudge_y,
                  max.overlaps = Inf, force = 14, force_pull = 0.12,
                  max.iter = 50000, max.time = 3,
                  box.padding = 0.7, point.padding = 0.3, min.segment.length = 0,
                  segment.size = 0.2, segment.colour = "grey55", segment.alpha = 0.7,
                  bg.color = "white", bg.r = 0.12, show.legend = FALSE) +
  facet_wrap2(~ region, nrow = 1, scales = "free_y",
              # The Hippo banner was already grey (low_conf);
              # add the matching "(low n)" suffix so the astro volcano carries the same
              # figure-wide single-lane-Hippo caveat as every other astro/oligo panel.
              labeller = region_full_lowconf_labeller(low_conf = "Hippo"),
              strip = strip_region_x(REGIONS, low_conf = "Hippo", clip = "off")) +
  ggh4x::facetted_pos_scales(y = y_scales) +
  scale_colour_manual(values = THEME_PAL, name = NULL, drop = FALSE, na.translate = FALSE,
                      guide = guide_legend(override.aes = list(size = 2.6, alpha = 1), ncol = 2, byrow = TRUE)) +
  scale_x_continuous(limits = c(-xl, xl), oob = scales::squish, expand = expansion(mult = 0.04),
                     breaks = scales::pretty_breaks(n = 5)) +
  labs(x = expression(log[2]~"fold change (NHD vs CON)"),
       y = expression(-log[10]~"(adjusted "*italic(p)*")")) +
  theme_classic(base_size = 8) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.35),
        axis.text = element_text(colour = "black"), axis.title = element_text(colour = "black"),
        strip.background = element_blank(), panel.spacing.x = unit(0.32, "lines"),
        legend.position = "bottom", legend.text = element_text(size = 7.4),
        legend.key.size = unit(0.34, "cm"), legend.key.spacing.x = unit(0.25, "cm"),
        legend.box.spacing = unit(0.1, "cm"), legend.margin = margin(t = 1, b = 0),
        plot.margin = margin(2, 3, 1, 3))

# right-side threshold KEY (verbatim from the micro volcano)
kt <- function(x, y, lab, size = 2.35, face = "plain", col = "black", parse = FALSE)
  annotate("text", x = x, y = y, label = lab, hjust = 0, size = size, colour = col, fontface = face, parse = parse)
UP_COL <- unname(PAL_DGE["Up in NHD"]); DN_COL <- unname(PAL_DGE["Up in CON"]); NS_COL <- unname(PAL_DGE["ns"])
key <- ggplot() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  kt(0.00, 0.985, "Coloured = theme gene:", size = 2.5) +
  kt(0.06, 0.915, "adjusted~italic(p)~\"< 0.05\"", parse = TRUE) +
  kt(0.06, 0.855, "\"fold change \" >= \"1.20x\"", parse = TRUE) +
  annotate("point", x = 0.03, y = 0.775, colour = "grey45", size = 0.9, shape = 16) +
  kt(0.10, 0.775, "nominal (small)") +
  annotate("point", x = 0.03, y = 0.700, colour = "grey45", size = 2.3, shape = 16) +
  kt(0.10, 0.700, "\"strict: |\"*delta*\"| \" >= 0.15", parse = TRUE) +
  annotate("segment", x = 0.02, xend = 0.02, y = 0.585, yend = 0.635, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  kt(0.10, 0.610, "fold-change floor") +
  annotate("segment", x = 0.00, xend = 0.05, y = 0.520, yend = 0.520, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  kt(0.10, 0.520, "italic(p)~\"= 0.05\"", parse = TRUE) +
  annotate("point", x = 0.03, y = 0.390, colour = UP_COL, size = 2.4, shape = 16) +
  kt(0.10, 0.390, "up in NHD (warm)") +
  annotate("point", x = 0.03, y = 0.320, colour = DN_COL, size = 2.4, shape = 16) +
  kt(0.10, 0.320, "lost in NHD (cool)") +
  annotate("point", x = 0.03, y = 0.250, colour = NS_COL, size = 2.4, shape = 16) +
  kt(0.10, 0.250, "not significant") +
  theme_void() + theme(plot.margin = margin(2, 2, 2, 1))

combo <- p + key + patchwork::plot_layout(widths = c(1, 0.34))

W <- 4.80; H <- 2.8
png_path <- file.path(PANEL, "F3b_astro_volcano_MAST.png")
pdf_path <- file.path(PANEL, "F3b_astro_volcano_MAST.pdf")
# Res was 300 -- the only 300-dpi panel in Fig 3, and it would render soft beside its
# 600-dpi siblings once placed at the same on-page size.
ragg::agg_png(png_path, width = W, height = H, units = "in", res = 600)
print(combo); invisible(dev.off())
ggsave(pdf_path, combo, width = W, height = H, useDingbats = FALSE)
cat(sprintf("\nWrote panel: %s\n         and: %s\n", png_path, pdf_path))

prov <- file.path(PANEL, "F3b_astro_volcano_MAST.provenance.txt")
sink(prov)
cat("F3b_astro_volcano_MAST (FH) — provenance\n")
cat("generated:", format(Sys.time()), "\n")
cat("script   : scripts/F3b_astro_volcano_MAST_FH.R\n")
cat("input    :", MAST_PATH, "\n")
cat("filter   : cell_type==Astro, region in {Frontal,Hippo}\n")
cat("YCAP:", YCAP, " | padj underflow floored to xmin:", n_uf, "rows\n\n")
cat("== discovery counts per region ==\n"); print(as.data.frame(counts))
cat("\n== headline-gene direction check ==\n"); print(as.data.frame(chk))
cat("\n== labelled (themed) genes ==\n")
print(as.data.frame(lab_df[, c("region","dir_col","theme","gene","avg_log2FC","cliffs_delta")]))
cat("\n== sessionInfo ==\n"); print(sessionInfo())
sink()
cat(sprintf("Wrote provenance: %s\n", prov))
cat("\n=== DONE ===\n", file = stderr())
