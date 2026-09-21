#!/usr/bin/env Rscript
# =============================================================================
# 3J_oligo_dissociation_curves_FH.R — Figure 4 (main): curated 3x4 pseudotime dissociation panel. CON-vs-NHD loess curves along the oligodendrocyte maturation
# continuum for 12 fixed genes covering the four arms of the white-matter lesion:
#   structural myelin : PLP1, CNP (flat), OPALIN (up)   — NHD tracks at-or-above CON
#   cholesterol/SREBP2: LSS, DHCR24, FDFT1 (down)       — NHD sits below CON
#   myelin LIPID      : SCD, SGMS1 (down), UGT8 (flat) — the sterol and fatty-acid/
#                       sphingomyelin arms fail; the galactolipid arm is preserved
#   stress            : CRYAB, FTL, APOD (up)           — NHD rises above CON
#
# The one-glance reading, and the reason the panel was extended from 6 genes to 12
# ("add more relevant genes since we are exploring different aspects
# of the white matter dysregulation"): the NHD oligodendrocyte is not losing its myelin
# programme per cell — PLP1/CNP overlap and OPALIN rises — while both lipid supply arms
# that have to build the membrane fall away beneath it, sterol synthesis and myelin
# lipid metabolism alike. The cell keeps ordering myelin protein it can no longer
# lipidate. The myelin deficit is therefore COMPOSITIONAL (fewer oligodendrocytes;
# carried by the UMAP/density/composition panels and by 4_myelin_scales_bridge), while
# the per-cell failure is the lipid supply chain.
#
# The lipid arm is NEW here and its genes were chosen by detection in these
# nuclei — SCD 0.65, UGT8 0.69, SGMS1 0.53 detected — after the audit that showed the old
# canonical lipid list (ELOVL4/ENPP1/GBA/CERS5) sits at the floor and made the programme
# read flat. OPALIN is added because it is the one myelin gene
# genuinely induced rather than merely share-inflated (detection 0.22 CON -> 0.55 NHD).
# -----------------------------------------------------------------------------
# Curated companion to the dense 24-gene supp grid (3J_pt_pathway_curves_FH.R,
# now in _supp/). same cache/CDS + same loess-along-pseudotime machinery, but a
# fixed 6-gene selection (no data-driven ranking) so the panel is deterministic
# and reads the thesis directly. Reads the shared oligo cache + CDS built by
# 3F_oligo_embed_states_FH.R (stale-cache guarded).
#
# House rules: no bold; no in-panel caption (the myelin/cholesterol/stress role is
# the row label, a plain axis title, not a caption); CON=blue / NHD=red (PAL_COND);
# canonical no-unicode-arrow text. Pooled across regions (this is the per-cell
# programme readout along the continuum, not a region contrast).
#
# Output: figures/Figure_4/panels/F4i_oligo_dissociation_curves.{png,pdf} + prov.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2); library(patchwork); library(Matrix)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot("PROJ root not found" = !is.na(PROJ) && dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, PAL_COND
DDIR  <- file.path(PROJ, "data")
TDIR  <- file.path(PROJ, "tables")
PANEL <- file.path(PROJ, "figures", "Figure_4", "panels")         # oligo -> Figure 4
LOGD  <- file.path(PROJ, "logs")
for (d in c(TDIR, PANEL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
CACHE <- file.path(DDIR, "_cache_oligo_umap.rds")
CDS_C <- file.path(DDIR, "_cds_oligo_FH.rds")
stopifnot("MISSING atlas" = file.exists(ATLAS),
          "MISSING oligo cache — run 3F first" = file.exists(CACHE),
          "MISSING oligo CDS — run 3F first"   = file.exists(CDS_C))
if (file.mtime(CACHE) < file.mtime(ATLAS) || file.mtime(CDS_C) < file.mtime(CACHE))
  stop("STALE oligo cache/CDS: re-run 3F_oligo_embed_states_FH.R.")

CON_COL <- unname(PAL_COND["CON"]); NHD_COL <- unname(PAL_COND["NHD"])

# --- the curated 6-gene dissociation set (fixed; not data-driven) ------------
# role = the dissociation axis each gene stands for; order fixes the 2x3 grid
# (row 1 = Myelin, row 2 = Cholesterol, row 3 = Stress -> but laid out as 2 cols x
# 3 rows via facet_wrap ncol=2, so the pairs sit side-by-side within each role).
gene_role <- tibble::tribble(
  ~gene,     ~role,                        ~expect,
  "PLP1",    "Structural myelin",          "flat",
  "CNP",     "Structural myelin",          "flat",
  "OPALIN",  "Structural myelin",          "up",
  "LSS",     "Cholesterol / SREBP2",       "down",
  "DHCR24",  "Cholesterol / SREBP2",       "down",
  "FDFT1",   "Cholesterol / SREBP2",       "down",
  "SCD",     "Myelin lipid",               "down",
  "SGMS1",   "Myelin lipid",               "down",
  # UGT8 is expected flat, not down, and that is a result rather than a weak gene: the
  # galactolipid arm of myelin lipid synthesis is PRESERVED (mean NHD-CON -0.08) while the
  # sterol arm and the fatty-acid/sphingomyelin arm both fall. Keeping it in the panel
  # states which parts of the lipid supply chain fail and which do not, and pre-empts the
  # reviewer question of whether only the moving genes were shown.
  "UGT8",    "Myelin lipid",               "flat",
  "CRYAB",   "Stress / reactive",          "up",
  "FTL",     "Stress / reactive",          "up",
  "APOD",    "Stress / reactive",          "up")
ROLE_ORDER <- c("Structural myelin","Cholesterol / SREBP2","Myelin lipid","Stress / reactive")
GENE_ORDER <- gene_role$gene
GRID_N <- 60L; LOESS_SUBS <- 4000L; DOTS_PER <- 1200L; SPAN <- 0.7

# ── load cache + CDS, attach pseudotime (identical prep to 3J supp) ───────────
cat("== load oligo cache + CDS, attach pseudotime ==\n")
obj <- readRDS(CACHE)
cds <- readRDS(CDS_C)
DefaultAssay(obj) <- "RNA"
if (length(SeuratObject::Layers(obj, assay = "RNA")) > 1) obj <- SeuratObject::JoinLayers(obj)
obj <- NormalizeData(obj, verbose = FALSE)
pt <- monocle3::pseudotime(cds); pt[is.infinite(pt)] <- NA_real_
common <- intersect(colnames(obj), names(pt))
obj$pt <- NA_real_; obj$pt[common] <- pt[common]

# depth-matched expression. The panel previously plotted
# LogNormalize values, i.e. each gene's share of its nucleus — the exact quantity
# that made the pooled myelin claim need correcting (REBUILD_PLAN entries 27/31).
# NHD frontal oligodendrocytes carry half the transcriptome of a control
# oligodendrocyte (median 1520 -> 799 UMI), so on a share axis "NHD myelin tracks
# at or above CON" is partly the denominator, and a reviewer can say so.
# Every nucleus is now binomially thinned to an identical depth before scoring, so
# the curves compare like with like and the dissociation cannot be dismissed as
# normalisation. The thinning removes nuclei below the floor, which costs NHD more
# than CON; retention is printed and belongs in the legend.
DEPTH_T <- 750L
cts_all <- GetAssayData(obj, assay = "RNA", layer = "counts")
tot_all <- Matrix::colSums(cts_all)
keep_d  <- tot_all >= DEPTH_T
for (cc in c("CON","NHD"))
  cat(sprintf("depth-match retention %s: %d/%d nuclei (%.0f%%)\n", cc,
              sum(keep_d & obj$Condition == cc), sum(obj$Condition == cc),
              100 * sum(keep_d & obj$Condition == cc) / sum(obj$Condition == cc)))
set.seed(42)
.ds <- Seurat::SampleUMI(cts_all[, keep_d, drop = FALSE], max.umi = DEPTH_T,
                         upsample = FALSE, verbose = FALSE)
dimnames(.ds) <- dimnames(cts_all[, keep_d, drop = FALSE])
data_mat <- log1p(.ds / DEPTH_T * 1e4)   # CP10K on an identical denominator
obj <- obj[, keep_d]
cells_keep <- which(!is.na(obj$pt) & obj$Condition %in% c("CON","NHD"))
meta <- data.frame(cell = colnames(obj)[cells_keep], pt = obj$pt[cells_keep],
                   Condition = factor(obj$Condition[cells_keep], levels = c("CON","NHD")))
cat(sprintf("cells with pseudotime + condition: %d (CON %d, NHD %d)\n",
            nrow(meta), sum(meta$Condition=="CON"), sum(meta$Condition=="NHD")))

# --- fail loud if any curated gene is absent (no silent drop from a main panel) --
missing <- setdiff(GENE_ORDER, rownames(data_mat))
if (length(missing))
  stop(sprintf("MISSING curated gene(s) from the oligo matrix: %s — cannot build the fixed panel.",
               paste(missing, collapse = ", ")))

# small gene x cell dense slice (12 genes only — no OOM risk)
expr_dense <- as.matrix(data_mat[GENE_ORDER, meta$cell, drop = FALSE])

# ── long frame for the curves + a per-condition dot subsample ────────────────
expr_long <- as.data.frame(t(expr_dense))
expr_long$cell <- rownames(expr_long)
expr_long <- expr_long %>%
  pivot_longer(cols = all_of(GENE_ORDER), names_to = "gene", values_to = "expr") %>%
  left_join(meta, by = "cell") %>%
  left_join(gene_role %>% select(gene, role), by = "gene") %>%
  mutate(gene = factor(gene, levels = GENE_ORDER),
         role = factor(role, levels = ROLE_ORDER),
         Condition = factor(Condition, levels = c("CON","NHD")))

set.seed(42)
dots_sub <- expr_long %>% group_by(gene, Condition) %>%
  slice_sample(n = DOTS_PER) %>% ungroup()

# ── direction sanity: compare NHD vs CON mean loess on the shared pt support ──
# reviewer-facing: confirm the curated thesis holds (myelin ~flat, chol down, stress
# up) before committing to the panel; report the mean NHD-CON curve gap per gene.
fit_dir <- function(g) {
  d <- expr_long %>% filter(gene == g, !is.na(pt))
  dc <- d %>% filter(Condition == "CON"); dn <- d %>% filter(Condition == "NHD")
  rng <- c(max(min(dc$pt), min(dn$pt)), min(max(dc$pt), max(dn$pt)))
  if (diff(rng) <= 0) return(NA_real_)
  grid <- seq(rng[1], rng[2], length.out = GRID_N)
  subs <- function(df) if (nrow(df) <= LOESS_SUBS) df else df[sample.int(nrow(df), LOESS_SUBS), ]
  fc <- predict(loess(expr ~ pt, data = subs(dc), span = SPAN),
                newdata = data.frame(pt = grid))
  fn <- predict(loess(expr ~ pt, data = subs(dn), span = SPAN),
                newdata = data.frame(pt = grid))
  mean(fn - fc, na.rm = TRUE)   # +ve = NHD above CON
}
dir_tab <- gene_role %>% rowwise() %>%
  mutate(mean_NHD_minus_CON = round(fit_dir(gene), 3)) %>% ungroup() %>%
  mutate(reads = case_when(
    expect == "flat" ~ ifelse(abs(mean_NHD_minus_CON) <= 0.20, "OK flat/>=CON",
                              ifelse(mean_NHD_minus_CON > 0, "NHD>CON (still not a per-cell loss) OK","CHECK: NHD<CON")),
    expect == "down" ~ ifelse(mean_NHD_minus_CON < 0, "OK NHD<CON", "CHECK: not down"),
    expect == "up"   ~ ifelse(mean_NHD_minus_CON > 0, "OK NHD>CON", "CHECK: not up")))
cat("\n== curated dissociation direction check (mean NHD-CON along shared pseudotime) ==\n")
print(as.data.frame(dir_tab), row.names = FALSE)

# write the source table (legend numbers key here)
src_path <- file.path(TDIR, "fig4_oligo_dissociation_curves_FH.csv")
write.csv(dir_tab, src_path, row.names = FALSE)
cat(sprintf("wrote %s\n", src_path))

# per-facet gap annotation, positioned in each facet's own free-y space
gap_lab <- expr_long %>% group_by(gene, role) %>%
  summarise(x = min(pt, na.rm = TRUE), y = max(expr, na.rm = TRUE) * 1.02, .groups = "drop") %>%
  left_join(dir_tab %>% select(gene, mean_NHD_minus_CON), by = "gene") %>%
  mutate(lab = sprintf("Delta == '%+.2f'", mean_NHD_minus_CON))

# ── build the 3x4 curated panel ──────────────────────────────────────────────
# facet_wrap by gene, ncol=3 so each role triple sits on one row (gene order fixes
# the layout). Gene names italic; the role is not captioned in-panel (it lives in
# the legend, per house rule). free_y so each gene's own dynamic range is visible.
p <- ggplot(expr_long, aes(pt, expr, colour = Condition, fill = Condition)) +
  geom_point(data = dots_sub, aes(colour = Condition), size = 0.10,
             alpha = 0.10, stroke = 0) +
  geom_smooth(method = "loess", span = SPAN, se = TRUE, linewidth = 0.9,
              alpha = 0.18, na.rm = TRUE) +
  # NESTED strips: the four ROLES are the argument, and with a plain gene facet the
  # reader had to hold the grouping in their head from the legend. ggh4x nests the
  # role above its three genes; because each role has exactly 3 genes and ncol = 3,
  # each role gets its own labelled row. This is a label, not an in-panel caption.
  ggh4x::facet_nested_wrap(dplyr::vars(role, gene), ncol = 3, scales = "free_y",
                           nest_line = element_line(colour = "grey70", linewidth = 0.3),
                           # by_layer_x = TRUE is required: without it elem_list_* recycles
                           # its values across individual STRIPS, so the faces alternated
                           # (CNP plain, OPALIN italic, PLP1 plain...) and the role strips
                           # alternated with them. by_layer_x applies one setting per NESTING
                           # layer, which is what "roles plain, genes italic" means. Same
                           # recycling trap as strip_region_x's fontsize.
                           strip = ggh4x::strip_nested(
                             background_x = ggh4x::elem_list_rect(
                               fill = c("grey92", "#EDEEF2"), colour = NA),
                             text_x = ggh4x::elem_list_text(
                               size = c(7.6, 8), face = c("plain", "italic")),
                             by_layer_x = TRUE)) +
  # mean NHD-CON curve gap per gene — the same number the direction check writes to
  # the source CSV, so the panel and the legend cannot disagree. House idiom: an
  # effect size annotated on the panel, never a p-value.
  # Size 2.05 -> 2.9 and grey35 -> grey20 ("very hard to see"). The
  # panel had been reduced 30% since these were added, which shrank them on the page as
  # well; they carry the panel's effect size, so they have to be readable at placed size.
  geom_text(data = gap_lab, aes(x = x, y = y, label = lab), inherit.aes = FALSE,
            hjust = 0, vjust = 1, size = 2.9, colour = "black", parse = TRUE) +
  scale_colour_manual(values = c(CON = CON_COL, NHD = NHD_COL), name = NULL) +
  scale_fill_manual(values = c(CON = CON_COL, NHD = NHD_COL), name = NULL, guide = "none") +
  scale_y_continuous(breaks = scales::pretty_breaks(2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(3)) +
  # ASCII "to" (no unicode arrow) — base pdf() cannot encode U+2192.
  labs(x = "Pseudotime  (OPC-proximal to Mature to Stress / High-myelin)",
       y = "Expression at matched depth (750 UMI)") +
  theme_pub(base_size = 8) +
  # Pale house gene-strip tone: a soft neutral pale so the 3j facet
  # strips share the Fig4 pale-banner idiom instead of the default ggplot grey. Gene
  # strips are not region/state strips, so a neutral pale (not a region/state colour)
  # is used; text stays plain black italic.
  theme(strip.text = element_text(colour = "black", margin = margin(t = 1.5, b = 1.5)),
        axis.title.x = element_text(size = 7.6, margin = margin(t = 3)),
        axis.title.y = element_text(size = 8.4, margin = margin(r = 3)),
        axis.text = element_text(size = 6.6),
        panel.spacing = unit(0.35, "lines"),
        # Legend right: on top it consumed a full row of height in a
        # 4-row panel that had just been reduced 30%, and the CON/NHD key is only two
        # entries — it costs far less width on the right than it did height above.
        legend.position = "right", legend.key.size = unit(0.30, "cm"),
        legend.margin = margin(0, 0, 0, 2),
        legend.text = element_text(size = 8),
        plot.margin = margin(3, 4, 3, 3))

W <- 4.59; H <- 4.46   # 30% smaller, 20% longer, then 10% slimmer
BN <- "F4i_oligo_dissociation_curves"
png_path <- file.path(PANEL, paste0(BN, ".png"))
pdf_path <- file.path(PANEL, paste0(BN, ".pdf"))
ggsave(pdf_path, p, width = W, height = H, useDingbats = FALSE)
ggsave(png_path, p, width = W, height = H, dpi = 600, device = ragg::agg_png)
cat(sprintf("Saved %s.{pdf,png} (%.2f x %.2f in)\n", BN, W, H))

prov <- file.path(PANEL, paste0(BN, ".provenance.txt"))
sink(prov)
cat("F4i_oligo_dissociation_curves (FH) — provenance\n")
cat("generated:", format(Sys.time()), "\n")
cat("script   : scripts/3J_oligo_dissociation_curves_FH.R\n")
cat("inputs   : cache", CACHE, "\n           CDS", CDS_C, "\n")
cat("source CSV:", src_path, "\n")
cat(sprintf("curated genes (fixed): %s\n", paste(GENE_ORDER, collapse = ", ")))
cat(sprintf("LOESS span=%.2f, grid=%d, dots/cond=%d, loess subsample=%d\n",
            SPAN, GRID_N, DOTS_PER, LOESS_SUBS))
cat("\n== direction check ==\n"); print(as.data.frame(dir_tab))
cat("\n== sessionInfo ==\n"); print(sessionInfo())
sink()
cat(sprintf("Wrote provenance: %s\n", prov))
cat("=== DONE ===\n", file = stderr())
