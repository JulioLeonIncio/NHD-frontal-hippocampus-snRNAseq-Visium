#!/usr/bin/env Rscript
# =============================================================================
# 109_suppfig5_split_umaps_FH.R — Supplementary Fig. 5: the five main-figure embeddings shown separately by condition x region (comment 2; Task A).
# -----------------------------------------------------------------------------
# Why. Overlaid UMAPs hide group differences. Each embedding (Fig 1b atlas, 2a Micro-PVM,
# 3a astrocytes, 4a oligodendrocytes, 5a neurons) is re-drawn on its own joint coordinates in
# a 1 x 4 strip — region-major: CON frontal | NHD frontal | CON hippocampus | NHD hippocampus
#  — with identical limits, point size and the same colour scale as
# the main-figure panel (cell type / myeloid class / oligodendrocyte state / neuron class;
# astrocytes are one continuous population and are drawn in the astrocyte cell-type colour).
#
# View LIMITS = the main panel's. Where the main-figure script crops its view (radial
# quantile about the median: astro 0.97 / oligo 0.99 / neuron 0.97, 2 % pad — see
# 3_astro_UMAP_FH.R, 3F_oligo_embed_states_FH.R, 64_fig4_neuron_umap_FH.R) the same rule is
# applied here on the full embedding, so e.g. the 47-nucleus ependymal island that Fig 3a
# crops out is cropped here too; nuclei outside the frame are counted in the log and CSV
# (n_outside_frame) and are retained in the embedding. Fig 1b / 2a use the full range.
#
# Points match the main-figure siblings: shape 21, thin grey20 outline (stroke 0.04), fixed
# size (abundant populations must look abundant — [[nhd-umap-abundance-pointsize]]).
#
# Density is not biology. Each panel is DOWN-SAMPLED to the smallest of the four groups for
# that embedding (seed 42, re-seeded per embedding), so a denser cloud means a different
# distribution, not more nuclei. The full-n version is rendered beside it (figures/_diagnostics)
# and compared by eye so that down-sampling did not erase structure; n (drawn / total) is
# printed in each panel title. Axis arrows are drawn once per embedding (first panel).
#
# Coordinates and metadata come from the caches the main figures read (never recomputed):
#   1b  data/_cache_Fig1/umap_df.rds            (06d)
#   2a  data/_cache_2A_myeloid_umap_FH.rds      (2A, mtime-guarded; = Fig 2a's embedding)
#   3a  data/_cache_astro_umap_FH.rds           (3_astro_UMAP)
#   4a  data/_cache_oligo_umap.rds              (3F)
#   5a  data/_cache_neuron_umap_FH.rds          (64)
# Output: figures/Supplementary/SuppFig5_split_UMAPs/SuppFig5_<row>_<embedding>.{png,pdf}
#         figures/Supplementary/SuppFig5_split_UMAPs/SuppFig5_split_UMAPs.{pdf,png}  (5 rows, portrait)
#         tables/suppfig5_split_umap_counts_FH.csv   (n per panel, keyed for the legend)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork); library(ragg) })
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
OUT  <- file.path(PROJ, "figures", "Supplementary", "SuppFig5_split_UMAPs"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DIAG <- file.path(PROJ, "figures", "_diagnostics", "suppfig5_full_n"); dir.create(DIAG, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJ, "logs"), showWarnings = FALSE); dir.create(file.path(PROJ, "tables"), showWarnings = FALSE)
say <- function(...) cat(sprintf(...), "\n")
# region-major column order: the reader compares CON vs NHD within a region first.
GROUPS <- c("CON Frontal", "NHD Frontal", "CON Hippocampus", "NHD Hippocampus")
PT <- 0.55; PT_ALPHA <- 0.80          # fixed across rows (house: abundant must look abundant)

# Radial view crop, identical to the main-figure scripts (computed on the full embedding)
crop_lims <- function(df, q, pad = 0.02) {
  cx <- median(df$UMAP_1); cy <- median(df$UMAP_2)
  r  <- sqrt((df$UMAP_1 - cx)^2 + (df$UMAP_2 - cy)^2)
  keep <- r <= quantile(r, q, na.rm = TRUE)
  qx <- range(df$UMAP_1[keep]); qy <- range(df$UMAP_2[keep])
  list(x = c(qx[1] - pad * diff(qx), qx[2] + pad * diff(qx)), y = c(qy[1] - pad * diff(qy), qy[2] + pad * diff(qy)))
}
must <- function(f, hint) { if (!file.exists(f)) stop("MISSING ", f, " — ", hint); f }

# ---- assemble the five embeddings as data frames: UMAP_1, UMAP_2, Condition, Region, col ----
emb <- list()
# 1b atlas (cache of 06d; 8 canonical classes, as the main panel; full range like Fig 1b)
u <- readRDS(must(file.path(PROJ, "data", "_cache_Fig1", "umap_df.rds"), "run 06d_figure1_umap_compact_FH.R"))
CT_LEVELS <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh","Endo","Pericytes")
u <- u[u$cell_type %in% CT_LEVELS, ]
emb[["1b"]] <- list(title = "All nuclei (Fig. 1b)", crop = NA,
                    df = data.frame(UMAP_1 = u$UMAP_1, UMAP_2 = u$UMAP_2, Condition = u$Condition, Region = u$Region,
                                    col = factor(u$cell_type, levels = CT_LEVELS)),
                    pal = PAL_CELLTYPE[CT_LEVELS], lab = CT_DISPLAY[CT_LEVELS])
# 2a Micro-PVM: the exact embedding Fig 2a draws (2A's mtime-guarded cache; pooled myeloid,
# T-lymphocytes excluded as in 2A). Coloured by myeloid class with the immune palette that
# 2N/2O use (Microglia = Micro-PVM cell-type colour; CD163+/F13A1+ = brown), labels from IMM_LABEL.
ATLAS <- file.path(PROJ, "atlas", "NHD_FH_harmony.rds")
CLS_F <- file.path(PROJ, "tables", "micro_states", "immune_umap_classes_FH.csv")
C2A   <- must(file.path(PROJ, "data", "_cache_2A_myeloid_umap_FH.rds"), "run 2A_micro_UMAP_FH.R first (writes the Fig 2a embedding cache)")
if (file.exists(ATLAS) && file.mtime(C2A) < file.mtime(ATLAS)) stop("STALE ", C2A, " (older than the atlas) — re-run 2A_micro_UMAP_FH.R")
if (file.exists(CLS_F) && file.mtime(C2A) < file.mtime(CLS_F)) stop("STALE ", C2A, " (older than immune_umap_classes_FH.csv) — re-run 2A_micro_UMAP_FH.R")
m <- readRDS(C2A)$df
stopifnot(all(c("UMAP_1","UMAP_2","Condition","Region","myeloid_class") %in% names(m)), nrow(m) > 1000)
PAL_IMM <- c("Microglia" = unname(PAL_CELLTYPE[["Micro-PVM"]]), "PVM" = "#6D4C41")   # as 2N/2O (T-lymphocyte slot unused here)
emb[["2a"]] <- list(title = "Micro-PVM, pooled myeloid nuclei (Fig. 2a)", crop = NA,
                    df = data.frame(UMAP_1 = m$UMAP_1, UMAP_2 = m$UMAP_2, Condition = m$Condition, Region = m$Region,
                                    col = factor(m$myeloid_class, levels = names(PAL_IMM))),
                    pal = PAL_IMM, lab = IMM_LABEL[names(PAL_IMM)])
say("2a: %d myeloid nuclei from the Fig 2a cache (%s)", nrow(m), format(file.mtime(C2A), "%Y-%m-%d %H:%M"))
# 3a astrocytes (Seurat cache; one continuous population -> single colour; radial-0.97 crop as Fig 3a)
a <- readRDS(must(file.path(PROJ, "data", "_cache_astro_umap_FH.rds"), "run 3_astro_UMAP_FH.R")); ae <- Seurat::Embeddings(a, "umap")
emb[["3a"]] <- list(title = "Astrocytes (Fig. 3a)", crop = 0.97,
                    df = data.frame(UMAP_1 = ae[, 1], UMAP_2 = ae[, 2], Condition = a$Condition, Region = a$Region, col = factor("Astrocytes")),
                    pal = c(Astrocytes = unname(PAL_CELLTYPE[["Astro"]])), lab = c(Astrocytes = "Astrocytes"))
# 4a oligodendrocytes (states as in Fig. 4e's state palette, 3F; radial-0.99 crop as Fig 4a)
o <- readRDS(must(file.path(PROJ, "data", "_cache_oligo_umap.rds"), "run 3F_oligo_embed_states_FH.R")); oe <- Seurat::Embeddings(o, "umap")
PAL_STATE <- c("OPC-proximal (early)" = "#CCBB44", "Mature (continuum)" = "#4477AA", "High-myelin (NHD)" = "#228833", "Stress (NHD)" = "#AA3377")
emb[["4a"]] <- list(title = "Oligodendrocytes (Fig. 4a)", crop = 0.99,
                    df = data.frame(UMAP_1 = oe[, 1], UMAP_2 = oe[, 2], Condition = o$Condition, Region = o$Region,
                                    col = factor(o$oligo_state, levels = names(PAL_STATE))),
                    pal = PAL_STATE, lab = setNames(names(PAL_STATE), names(PAL_STATE)))
# 5a neurons (class, teal/magenta; radial-0.97 crop as Fig 5a composite A)
n <- readRDS(must(file.path(PROJ, "data", "_cache_neuron_umap_FH.rds"), "run 64_fig4_neuron_umap_FH.R"))
cls <- ifelse(n$meta$new_annotation == "Neuron_Ex", "Excitatory", "Inhibitory")
emb[["5a"]] <- list(title = "Neurons (Fig. 5a)", crop = 0.97,
                    df = data.frame(UMAP_1 = n$embedding[, 1], UMAP_2 = n$embedding[, 2], Condition = n$meta$Condition, Region = n$meta$Region,
                                    col = factor(cls, levels = c("Excitatory", "Inhibitory"))),
                    pal = c(Excitatory = "#006D77", Inhibitory = "#C42E7B"), lab = c(Excitatory = "Excitatory", Inhibitory = "Inhibitory"))
rm(a, o, u, m, n); invisible(gc())

# ---- one 1 x 4 strip per embedding ----------------------------------------------------
counts <- list(); rows <- list(); rows_full <- list()
grid_of <- function(df, pal, lab, title, crop = NA, downsample = TRUE, seed = 42) {
  df$Region <- ifelse(df$Region %in% c("Hippo", "Hippocampus"), "Hippocampus", "Frontal")
  df$group  <- factor(paste(df$Condition, df$Region), levels = GROUPS)
  df <- df[!is.na(df$col) & !is.na(df$group), ]
  # Limits from the full embedding (main-panel rule), identical across the four panels
  lim <- if (is.na(crop)) list(x = range(df$UMAP_1), y = range(df$UMAP_2)) else crop_lims(df, crop)
  n_out <- sum(df$UMAP_1 < lim$x[1] | df$UMAP_1 > lim$x[2] | df$UMAP_2 < lim$y[1] | df$UMAP_2 > lim$y[2])
  n_tot <- table(df$group); n_min <- min(n_tot)
  set.seed(seed)
  if (downsample) df <- df %>% group_by(group) %>% slice_sample(n = n_min) %>% ungroup()
  n_drawn <- table(df$group)
  lev <- sprintf("%s\nn = %s of %s", GROUPS, trimws(format(n_drawn[GROUPS], big.mark = ",")), trimws(format(n_tot[GROUPS], big.mark = ",")))
  df$panel <- factor(lev[match(as.character(df$group), GROUPS)], levels = lev)
  df <- df[sample(nrow(df)), ]                                  # shuffled draw order: no group paints over another
  # Axis arrows once per embedding, in the first panel (bottom-left, inside the 5 % expansion)
  dx <- diff(lim$x); dy <- diff(lim$y); fr <- 0.16
  x0 <- lim$x[1] - 0.02 * dx; y0 <- lim$y[1] - 0.02 * dy
  arr <- data.frame(panel = factor(lev[1], levels = lev), x = x0, y = y0, xend = c(x0 + fr * dx, x0), yend = c(y0, y0 + fr * dy))
  arl <- data.frame(panel = factor(lev[1], levels = lev), x = c(x0 + 0.01 * dx, x0 - 0.02 * dx), y = c(y0 - 0.015 * dy, y0 + 0.01 * dy),
                    label = c("UMAP 1", "UMAP 2"), angle = c(0, 90), hjust = 0, vjust = c(1, 0))
  p <- ggplot(df, aes(UMAP_1, UMAP_2, fill = col)) +
    geom_point(size = PT, alpha = PT_ALPHA, stroke = 0.04, shape = 21, colour = "grey20") +
    geom_segment(data = arr, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
                 arrow = grid::arrow(length = grid::unit(1.1, "mm"), type = "closed"), linewidth = 0.3, colour = "black") +
    geom_text(data = arl, aes(x = x, y = y, label = label, angle = angle, hjust = hjust, vjust = vjust), inherit.aes = FALSE,
              size = 1.7, colour = "black") +
    scale_fill_manual(values = pal, labels = lab, name = NULL, drop = FALSE) +
    facet_wrap(~ panel, nrow = 1) +
    # square panels: UMAP axes are arbitrary; identical limits on all four (main-panel crop rule above)
    coord_cartesian(xlim = lim$x, ylim = lim$y, expand = TRUE, clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_pub(base_size = 7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          panel.border = element_blank(), aspect.ratio = 1,
          strip.text = element_text(size = 6.4, colour = "black", face = "plain"), strip.background = element_blank(),
          plot.title = element_text(size = 7.6, face = "plain", colour = "black", hjust = 0),
          legend.position = "right", legend.key.size = unit(0.35, "lines"), legend.text = element_text(size = 6, colour = "black"),
          panel.spacing = unit(6, "pt"), plot.margin = margin(2, 2, 4, 2)) +
    guides(fill = guide_legend(override.aes = list(shape = 21, colour = "grey20", stroke = 0.3, size = 2.4, alpha = 1), ncol = 1))
  list(p = p, n_tot = n_tot, n_drawn = n_drawn, n_out = n_out, lim = lim)
}
for (k in names(emb)) {
  e <- emb[[k]]
  g  <- grid_of(e$df, e$pal, e$lab, e$title, crop = e$crop, downsample = TRUE)
  gf <- grid_of(e$df, e$pal, e$lab, paste(e$title, "— all nuclei"), crop = e$crop, downsample = FALSE)
  rows[[k]] <- g$p; rows_full[[k]] <- gf$p
  counts[[k]] <- data.frame(embedding = k, group = GROUPS, n_total = as.integer(g$n_tot[GROUPS]), n_drawn = as.integer(g$n_drawn[GROUPS]),
                            crop_quantile = e$crop, n_outside_frame = g$n_out,
                            xlim_lo = g$lim$x[1], xlim_hi = g$lim$x[2], ylim_lo = g$lim$y[1], ylim_hi = g$lim$y[2])
  bn <- sprintf("SuppFig5_%s_%s", k, gsub("[^A-Za-z0-9]+", "_", tolower(sub(" \\(.*$", "", e$title))))
  ggsave(file.path(OUT, paste0(bn, ".png")), g$p, width = 7.0, height = 1.85, dpi = 500, device = ragg::agg_png)
  ggsave(file.path(OUT, paste0(bn, ".pdf")), g$p, width = 7.0, height = 1.85, useDingbats = FALSE)
  ggsave(file.path(DIAG, paste0(bn, "_FULL_N.png")), gf$p, width = 7.0, height = 1.85, dpi = 300, device = ragg::agg_png)
  say("%s: n per group %s (drawn %d each); crop %s -> %d/%d nuclei outside frame (retained in embedding)",
      k, paste(GROUPS, g$n_tot[GROUPS], collapse = "; "), min(g$n_tot), ifelse(is.na(e$crop), "none", e$crop), g$n_out, sum(g$n_tot))
}
counts <- bind_rows(counts); write.csv(counts, file.path(PROJ, "tables", "suppfig5_split_umap_counts_FH.csv"), row.names = FALSE)
# stale-twin guard: earlier renders named row b after the SuppFig 3c embedding; remove so no stale panel lingers
for (f in c(list.files(OUT, pattern = "^SuppFig5_2a_myeloid_and_lymphoid_nuclei\\.(png|pdf)$", full.names = TRUE),
            list.files(DIAG, pattern = "^SuppFig5_2a_myeloid_and_lymphoid_nuclei_FULL_N\\.png$", full.names = TRUE))) { file.remove(f); say("removed superseded %s", basename(f)) }

# ---- composite: 5 rows, portrait 509 x 691 pt (7.07 x 9.6 in) like every other figure ------
comp <- wrap_plots(rows, ncol = 1) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold", colour = "black"))   # panel letters are the one sanctioned bold (house)
ggsave(file.path(OUT, "SuppFig5_split_UMAPs.pdf"), comp, width = 7.07, height = 9.6, useDingbats = FALSE)
ggsave(file.path(OUT, "SuppFig5_split_UMAPs.png"), comp, width = 7.07, height = 9.6, dpi = 400, device = ragg::agg_png)
say("wrote composite SuppFig5_split_UMAPs.{pdf,png} and per-row panels in %s", OUT)
writeLines(c(sprintf("run: %s", format(Sys.time())), sprintf("2a cache: %s (%s)", C2A, format(file.mtime(C2A))), "", capture.output(sessionInfo())),
           file.path(PROJ, "logs", "109_suppfig5_split_umaps_FH_sessionInfo.txt"))
cat("\n=== DONE ===\n", file = stderr())
