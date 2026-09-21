#!/usr/bin/env Rscript
# =============================================================================
# 4_oligolineage_volcano_MAST_FH.R — Figure 4: thematic volcano for the oligodendrocyte lineage. One script, two cell types:
#     VOLC_CT=Oligo (default) -> 4_oligo_volcano_MAST
#     VOLC_CT=OPC             -> 4_opc_volcano_MAST
# -----------------------------------------------------------------------------
# Faithful port of 3_astro_volcano_MAST_FH.R, which is the house gold standard for
# this panel type (2-region facet, two-tier significance encoding, plotmath key,
# capped ggrepel labels). Only the theme map, the palette and the output names
# differ — the plotting grammar is copied so Figs 2, 3 and 4 read identically.
#
# THEMES come from scripts/_oligo_programmes_FH.R so the volcano, the programme
# dot-matrix and the gene heatmap cannot drift apart. Each theme carries an EXPECTED
# direction and a gene is only coloured when it MOVES that WAY (dir_match), which is
# what keeps a "lost" theme from being coloured by an up-gene and vice versa.
#
# The Oligo theme set is the supply-chain argument: myelin protein up, every lipid
# arm that has to supply the membrane down. TF/ANLN are their own theme because they
# fall while the membrane proteins rise (and ANLN replicates Zhou 2023).
# The OPC set is the lineage-arrest argument: identity/commitment markers down.
# Galactolipid is deliberately absent from both: it is flat (UGT8 +0.07/-0.08,
# FA2H 0.00/+0.05, GALC -0.01/+0.04), so no gene would satisfy dir_match and an
# always-empty legend entry would mislead. Its flatness is carried by the heatmap.
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
PANEL     <- file.path(PROJ, "figures", "Figure_4", "panels")
LOGS      <- file.path(PROJ, "logs")
for (f in c(MAST_PATH, ARTIFACT, THEME))
  if (!file.exists(f)) stop(sprintf("MISSING %s — prep step failed", f))
for (d in c(PANEL, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
source(ARTIFACT); source(THEME)

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
CT <- Sys.getenv("VOLC_CT", "Oligo")
stopifnot("VOLC_CT must be Oligo or OPC" = CT %in% c("Oligo", "OPC"))
BN <- if (CT == "Oligo") "4_oligo_volcano_MAST" else "4_opc_volcano_MAST"
cat(sprintf("== thematic volcano for %s -> %s ==\n", CT, BN))
source(file.path(PROJ, "scripts", "_oligo_programmes_FH.R"))
.sets <- oligo_programme_sets(readRDS(file.path(PROJ, "data", "curated_signatures.rds")))
mk <- function(genes, label) setNames(rep(label, length(genes)), genes)
YCAP    <- 60
FC_MIN  <- 1.20
LFC_MIN <- log2(FC_MIN)
CD_CUT   <- 0.15
PADJ_CUT <- 0.05

THEME_PAL_MASTER <- c(
  "Myelin protein (up)"               = "#C1272D",  # red
  "Stress / reactive (up)"            = "#D81B60",  # magenta
  "OPC state shift (up)"              = "#E8A33D",  # amber
  "Sterol synthesis (lost)"           = "#2E86AB",  # blue
  "Fatty-acid / sphingomyelin (lost)" = "#1B9E77",  # green
  "Lipid import (lost)"               = "#7E57A6",  # purple
  "Iron / scaffold (lost)"            = "#8C564B",  # brown
  "OPC identity / commitment (lost)"  = "#16456B")  # navy

# --- themed oligodendrocyte-lineage programs (from _oligo_programmes_FH.R) ----
if (CT == "Oligo") {
  ASTRO_THEME <- c(
    mk(.sets[["Structural myelin"]],           "Myelin protein (up)"),
    mk(.sets[["Stress / reactive"]],           "Stress / reactive (up)"),
    mk(.sets[["Cholesterol (sterol arm)"]],    "Sterol synthesis (lost)"),
    mk(.sets[["Fatty-acid / sphingomyelin"]],  "Fatty-acid / sphingomyelin (lost)"),
    mk(.sets[["Lipid uptake / salvage"]],      "Lipid import (lost)"),
    mk(oligo_myelin_support_genes,             "Iron / scaffold (lost)"))
  THEME_UP   <- c("Myelin protein (up)", "Stress / reactive (up)")
  THEME_LOST <- c("Sterol synthesis (lost)", "Fatty-acid / sphingomyelin (lost)",
                  "Lipid import (lost)", "Iron / scaffold (lost)")
} else {
  # OPC: the lineage-arrest argument. PDGFRA and GPR17 are the load-bearing genes
  # (down in both regions); GPR17 marks differentiation-committed pre-oligodendrocytes,
  # so its loss says the pipeline that would replace the lost myelinating cells is not
  # being fed. LHFPL3 and XYLT1 go up and are given their own theme rather than being
  # averaged against the losses.
  ASTRO_THEME <- c(
    mk(c("PDGFRA","GPR17","SOX6","BCAN","MEGF11","DSCAM","OLIG2","OLIG1"),
       "OPC identity / commitment (lost)"),
    mk(c("LHFPL3","XYLT1"),                    "OPC state shift (up)"),
    mk(.sets[["Stress / reactive"]],           "Stress / reactive (up)"),
    mk(.sets[["Cholesterol (sterol arm)"]],    "Sterol synthesis (lost)"),
    mk(.sets[["Fatty-acid / sphingomyelin"]],  "Fatty-acid / sphingomyelin (lost)"))
  THEME_UP   <- c("OPC state shift (up)", "Stress / reactive (up)")
  THEME_LOST <- c("OPC identity / commitment (lost)", "Sterol synthesis (lost)",
                  "Fatty-acid / sphingomyelin (lost)")
}
THEME_LEVELS <- c(THEME_UP, THEME_LOST)
stopifnot("a theme has no master colour" = all(THEME_LEVELS %in% names(THEME_PAL_MASTER)))
THEME_PAL <- THEME_PAL_MASTER[THEME_LEVELS]
if (Sys.getenv("VOLC_UNION_LEVELS", "0") == "1") {
  THEME_LEVELS <- names(THEME_PAL_MASTER)
  THEME_PAL    <- THEME_PAL_MASTER
  cat("union theme levels ON (composite mode): both panels carry all",
      length(THEME_LEVELS), "themes so their legends merge\n")
}

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

n0 <- nrow(d); cat(sprintf("%s rows across %d regions: %d\n", CT, length(REGIONS), n0))
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
LAB_CAP  <- 8L
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

if (Sys.getenv("VOLC_SAVE_GG", "0") == "1") {
  gg_path <- file.path(PROJ, "data", paste0(BN, "_gg.rds"))
  saveRDS(list(p = p, key = key, CT = CT), gg_path)
  cat("stashed plot objects for the composite ->", basename(gg_path), "\n")
}
combo <- p + key + patchwork::plot_layout(widths = c(1, 0.34))

W <- 4.80; H <- 2.8
png_path <- file.path(PANEL, paste0(BN, ".png"))
pdf_path <- file.path(PANEL, paste0(BN, ".pdf"))
# Res was 300 -- the only 300-dpi panel in Fig 3, and it would render soft beside its
# 600-dpi siblings once placed at the same on-page size.
ragg::agg_png(png_path, width = W, height = H, units = "in", res = 600)
print(combo); invisible(dev.off())
ggsave(pdf_path, combo, width = W, height = H, useDingbats = FALSE)
cat(sprintf("\nWrote panel: %s\n         and: %s\n", png_path, pdf_path))

prov <- file.path(PANEL, paste0(BN, ".provenance.txt"))
sink(prov)
cat(sprintf("%s (FH) — provenance\n", BN))
cat("generated:", format(Sys.time()), "\n")
cat("script   : scripts/4_oligolineage_volcano_MAST_FH.R\n")
cat("input    :", MAST_PATH, "\n")
cat(sprintf("filter   : cell_type==%s, region in {Frontal,Hippo}\n", CT))
cat("YCAP:", YCAP, " | padj underflow floored to xmin:", n_uf, "rows\n\n")
cat("== discovery counts per region ==\n"); print(as.data.frame(counts))
cat("\n== headline-gene direction check ==\n"); print(as.data.frame(chk))
cat("\n== labelled (themed) genes ==\n")
print(as.data.frame(lab_df[, c("region","dir_col","theme","gene","avg_log2FC","cliffs_delta")]))
cat("\n== sessionInfo ==\n"); print(sessionInfo())
sink()
cat(sprintf("Wrote provenance: %s\n", prov))
cat("\n=== DONE ===\n", file = stderr())
