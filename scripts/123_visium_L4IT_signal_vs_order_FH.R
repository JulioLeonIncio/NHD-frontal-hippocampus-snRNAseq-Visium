#!/usr/bin/env Rscript
# =============================================================================
# 123_visium_L4IT_signal_vs_order_FH.R — Diagnostic: is layer IV missing from the NHD Visium map because L4 IT neurons are gone, or because they no longer form a band? snRNA-seq says the L4 IT share is unchanged
#   (21.5 % vs 21.0 %); the Visium L4 domain is 1,660 vs 77 spots. Test on the Visium object:
#   (i) per-spot L4 IT signature score (RORB + the Jorstad-DFC L4 IT-over-L5 IT genes: TWIST2, DPF3, NGB, SCHLAP1, LYPD6B,
#       plus RORB) as UMI per 10,000 and as a z within cortical spots of each section;
#   (ii) amount: mean score over cortical spots (L1/pia–L6, InN) per section, and the fraction of "L4-high" spots
#        (top 15 % of the pooled cortical distribution);
#   (iii) order: spatial coherence of the L4-high spots per section — Moran's I of the score on the hex graph (6 nearest
#        neighbours) and the mean fraction of a L4-high spot's neighbours that are also L4-high (vs the expectation under
#        random placement = the section's L4-high fraction). Control: banded (high coherence); NHD: same amount but
#        dispersed would support "order, not class".
#   Also reports, per section, the median cortical depth (distance to the nearest WM spot) of L4-high spots, and the width
#   of the RORB band in control.
# OUT: tables/visium_L4IT_signal_vs_order_FH.csv; figures/_diagnostics/L4IT_signal_vs_order.png (maps + summary).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr); library(ggplot2); library(patchwork) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R")); source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))
IH <- file.path(dirname(PROJ), "Visium", "integrated_harmony"); C2L <- file.path(dirname(PROJ), "Visium", "cell2location", "c2l_MAIN")
obj <- readRDS(file.path(IH, "NHD_frontal_integrated_harmony.rds"))
xy <- read.csv(file.path(C2L, "spot_coords.csv")); names(xy)[1] <- "spot_id"; rownames(xy) <- xy$spot_id
md <- obj@meta.data; md$spot_id <- rownames(md); md$domain <- unname(DOM_MAP[as.character(md$seurat_clusters)]); md$x <- xy[md$spot_id, "x"]; md$y <- xy[md$spot_id, "y"]
cnt <- GetAssayData(obj, assay = "Spatial", layer = "counts"); tot <- Matrix::colSums(cnt)
L4G <- intersect(c("RORB","TWIST2","DPF3","NGB","SCHLAP1","LYPD6B"), rownames(cnt)); cat("L4 IT signature genes:", paste(L4G, collapse = " "), "\n")
md$l4 <- 1e4 * Matrix::colSums(cnt[L4G, , drop = FALSE]) / tot
md$rorb <- 1e4 * cnt["RORB", ] / tot
ctx <- md$domain %in% c(DOM_GREY, "InN")
thr <- quantile(md$l4[ctx], 0.85); md$l4high <- ctx & md$l4 >= thr
UM <- 1000 / 81.5
res <- list(); maps <- list()
for (sec in levels(factor(md$sample_id))) {
  s <- md[md$sample_id == sec & !is.na(md$x), ]; sc <- s[ctx[match(s$spot_id, md$spot_id)], ]
  ok <- s$spot_id %in% sc$spot_id
  # 6-nearest-neighbour graph built on cortical spots only and symmetrised (a graph over all spots,
  # subset afterwards, under-connects cortical spots at the white-matter / pial border, differently by condition)
  Dc <- as.matrix(dist(s[ok, c("x","y")])) * UM; diag(Dc) <- Inf
  nnc <- t(apply(Dc, 1, function(r) order(r)[1:6])); Wc <- matrix(0, nrow(Dc), nrow(Dc)); for (i in seq_len(nrow(Dc))) Wc[i, nnc[i, ]] <- 1
  Wc <- pmax(Wc, t(Wc))
  zc <- s$l4[ok] - mean(s$l4[ok])
  moran <- (length(zc) / sum(Wc)) * (t(zc) %*% Wc %*% zc) / sum(zc^2)
  hi <- s$l4high; hic <- hi[ok]; frac_hi <- mean(hic); nb_same <- mean(sapply(which(hic), function(i) mean(hic[nnc[i, ]])))
  D <- as.matrix(dist(s[, c("x","y")])) * UM; diag(D) <- Inf   # all-spot distances only for the depth-to-WM statistic
  wm <- s$domain == "WM"; depth <- if (any(wm)) apply(D[, wm, drop = FALSE], 1, min) else rep(NA, nrow(s))
  res[[sec]] <- data.frame(section = sec, n_cortical = sum(ok), mean_L4IT_per10k = round(mean(s$l4[ok]), 2), mean_RORB_per10k = round(mean(s$rorb[ok]), 2),
    frac_L4high = round(frac_hi, 3), neighbour_L4high_given_L4high = round(nb_same, 3), coherence_ratio = round(nb_same / frac_hi, 2), morans_I = round(as.numeric(moran), 3),
    median_depth_um_L4high = round(median(depth[hi & ok], na.rm = TRUE)), n_L4_domain = sum(s$domain == "L4"))
  maps[[sec]] <- ggplot(s, aes(x, y, colour = pmin(l4, quantile(md$l4[ctx], .98)))) + geom_point(size = 0.9, stroke = 0) + scale_colour_viridis_c(option = "B", name = "L4 IT /10k") +
    coord_equal() + theme_void(base_size = 8) + labs(title = sprintf("%s  I = %.2f  coh %.1f", sec, moran, nb_same / frac_hi))
}
out <- bind_rows(res); print(out, row.names = FALSE); write.csv(out, file.path(PROJ, "tables", "visium_L4IT_signal_vs_order_FH.csv"), row.names = FALSE)
p <- wrap_plots(maps, nrow = 1) + plot_annotation(title = "L4 IT signature per spot (RORB TWIST2 DPF3 NGB SCHLAP1 LYPD6B); I = Moran's I over cortical spots; coh = neighbour co-occurrence of L4-high spots / chance")
ggsave(file.path(PROJ, "figures", "_diagnostics", "L4IT_signal_vs_order.png"), p, width = 16, height = 4.8, dpi = 150, bg = "white")
cat("=== DONE ===\n")
