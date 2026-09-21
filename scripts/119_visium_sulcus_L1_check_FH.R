#!/usr/bin/env Rscript
# =============================================================================
# 119_visium_sulcus_L1_check_FH.R — Triple check of the pathologists' anatomical annotation against the Visium domain labels, spot by spot, in all four sections.
# The dashed L1 / WM lines drawn on figures/Figure_6/Fig6_260918_anatomy.pdf were digitised from the
# domain-map row of that page -> tables/annotation_dashes_20260918_mapunits.csv.
# Dashed pixels are grouped into connected components and each component is assigned to the label written
# beside it on the page (L1 or WM) by POSITION (table ANNOT below; the assignment is printed so it is auditable).
#
# For every spot: distance to the nearest L1-line pixel and to the nearest WM-line pixel (um). Zones:
#   L1 zone  = within 100 um of an L1 line (the band between two lines two spots apart, plus one spot of margin)
#   WM zone  = within 100 um of a WM line, or enclosed by the WM annotation (NHD1 U-band / corners are drawn as
#              bands, so "within 100 um" catches them; the CON WM wedges are large -> also use the domain itself)
# Outputs
#   tables/visium_sulcus_check_zone_composition_FH.csv   domain composition of the L1 zone and WM zone per section
#   tables/visium_sulcus_check_markers_FH.csv            UMI per 10k of myelin / glia limitans / leptomeningeal /
#                                                         pan-neuronal / HBB for: WM-labelled spots in the L1 zone,
#                                                         L1/pia spots in the L1 zone, WM spots in the WM zone, WM spots
#                                                         far from both (> 300 um), per section
#   figures/_diagnostics/sulcus_check_<sec>.png          H&E + spots (domain) + the digitised lines; zoom on the L1 zone
#   figures/_diagnostics/sulcus_check_<sec>_L1zone.png   only WM / L1-pia / Vasc-immune spots in the L1 zone, big points
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr); library(ggplot2); library(patchwork); library(png); library(jsonlite) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R")); source(file.path(PROJ, "scripts", "_visium_domains_FH.R")); source(file.path(PROJ, "scripts", "_visium_L1_annotation_rule_FH.R"))
ROOT <- dirname(PROJ); IH <- file.path(ROOT, "Visium", "integrated_harmony"); C2L <- file.path(ROOT, "Visium", "cell2location", "c2l_MAIN")
HEP <- file.path(ROOT, "Visium", "_HE_fullres", "panels"); OUTD <- file.path(PROJ, "figures", "_diagnostics"); TAB <- file.path(PROJ, "tables")
SEC <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")
UM_PER_UNIT <- 1000 / 81.5   # 1 mm = 81.5 map units (HE_frames.json)

obj <- readRDS(file.path(IH, "NHD_frontal_integrated_harmony.rds"))
stopifnot(file.mtime(file.path(C2L, "spot_coords.csv")) > file.mtime(file.path(IH, "NHD_frontal_integrated_harmony.rds")) - 3600)
xy <- read.csv(file.path(C2L, "spot_coords.csv")); names(xy)[1] <- "spot_id"
stopifnot(all(colnames(obj) %in% xy$spot_id))
md <- obj@meta.data; md$spot_id <- rownames(md)
md$domain <- unname(DOM_MAP[as.character(md$seurat_clusters)]); md$domain[is.na(md$domain)] <- paste0("g", md$seurat_clusters[is.na(md$domain)])
md <- md %>% left_join(xy, by = "spot_id")
dash <- read.csv(file.path(TAB, "annotation_dashes_20260918_mapunits.csv"))

# ---- gene sets per 10k ------------------------------------------------------------------------
SETS <- list(myelin = c("MBP","PLP1","MOBP","MOG"), glia_limitans = c("GFAP","AQP4","ID3"), leptomeningeal = c("PTGDS","DCN","COL1A2","SLC6A13","CEMIP"),
             pan_neuronal = c("SNAP25","RBFOX3","SLC17A7"), HBB = "HBB", vascular = c("CLDN5","PECAM1","PDGFRB"), immune = c("CD74","C1QB","TYROBP"))
cnt <- GetAssayData(obj, assay = "Spatial", layer = "counts"); tot <- Matrix::colSums(cnt)
per10k <- sapply(SETS, function(g) { g <- intersect(g, rownames(cnt)); 1e4 * Matrix::colSums(cnt[g, , drop = FALSE]) / tot })
md <- cbind(md, per10k[md$spot_id, ])

# ---- group dashed pixels into components and label them by position -----------------------------------
comp_label <- function(d, sec) {
  d <- d[d$section == sec, ]; if (!nrow(d)) return(d)
  # single-linkage components with a 3-spot-pitch gap (25 map units) on a thinned pixel set
  d <- d[!duplicated(round(d[, c("x","y")] / 2)), ]
  hc <- hclust(dist(d[, c("x","y")]), method = "single"); d$comp <- cutree(hc, h = H_LINK)
  d
}
# hand assignment of components to the page labels, by the centroid position written on the page (audited in the log)
assign_comp <- function(d, sec) {
  ce <- d %>% group_by(comp) %>% summarise(n = n(), cx = mean(x), cy = mean(y), .groups = "drop")
  xr <- range(md$x[md$sample_id == sec]); yr <- range(md$y[md$sample_id == sec])
  fx <- function(v) (v - xr[1]) / diff(xr); fy <- function(v) (v - yr[1]) / diff(yr)    # 0..1 within the section frame (y up)
  ce$fx <- fx(ce$cx); ce$fy <- fy(ce$cy)
  ce$label <- ifelse(is_L1(sec, ce$fx, ce$fy), "L1", "WM")   # rule from _visium_L1_annotation_rule_FH.R
  cat("\n", sec, "dashed-line components (centroid as fraction of the frame, x right / y up):\n"); print(as.data.frame(ce[order(-ce$n), ]), digits = 3)
  d$label <- ce$label[match(d$comp, ce$comp)]; d
}
nn_dist <- function(px, py, qx, qy) { if (!length(qx)) return(rep(Inf, length(px))); D <- sqrt(outer(px, qx, "-")^2 + outer(py, qy, "-")^2); apply(D, 1, min) * UM_PER_UNIT }

zone_rows <- list(); marker_rows <- list()
for (sec in SEC) {
  d <- assign_comp(comp_label(dash, sec), sec)
  s <- md[md$sample_id == sec, ]
  s$d_L1 <- nn_dist(s$x, s$y, d$x[d$label == "L1"], d$y[d$label == "L1"])
  s$d_WM <- nn_dist(s$x, s$y, d$x[d$label == "WM"], d$y[d$label == "WM"])
  s$zone <- ifelse(s$d_L1 <= 100, "L1 zone (<=100 um of an L1 line)", ifelse(s$d_WM <= 100, "WM zone (<=100 um of a WM line)", ifelse(s$d_L1 > 300 & s$d_WM > 300, "far (>300 um from both)", "between")))
  md$d_L1[md$sample_id == sec] <- s$d_L1; md$d_WM[md$sample_id == sec] <- s$d_WM; md$zone[md$sample_id == sec] <- s$zone
  zone_rows[[sec]] <- s %>% count(zone, domain) %>% group_by(zone) %>% mutate(frac = round(n / sum(n), 3), section = sec) %>% ungroup()
  grp <- function(df, name) if (nrow(df)) data.frame(section = sec, group = name, n = nrow(df), t(round(colMeans(df[, names(SETS)]), 2)), median_UMI = median(df$nCount_Spatial), median_d_L1_um = round(median(df$d_L1)), check.names = FALSE) else NULL
  marker_rows[[sec]] <- bind_rows(
    grp(s[s$zone == "L1 zone (<=100 um of an L1 line)" & s$domain == "WM", ], "WM-labelled spots INSIDE the L1 zone"),
    grp(s[s$zone == "L1 zone (<=100 um of an L1 line)" & s$domain == "L1/pia", ], "L1/pia spots inside the L1 zone"),
    grp(s[s$zone == "L1 zone (<=100 um of an L1 line)" & s$domain == "Vasc/immune", ], "Vasc/immune spots inside the L1 zone"),
    grp(s[s$zone == "L1 zone (<=100 um of an L1 line)" & s$domain %in% DOM_GREY & s$domain != "L1/pia", ], "cortical-layer spots inside the L1 zone"),
    grp(s[s$zone == "WM zone (<=100 um of a WM line)" & s$domain == "WM", ], "WM spots inside the WM zone"),
    grp(s[s$d_L1 > 300 & s$domain == "WM", ], "WM spots > 300 um from any L1 line (true WM)"),
    grp(s[s$d_L1 > 300 & s$domain == "L1/pia", ], "L1/pia spots > 300 um from any L1 line (check: should be few)"),
    grp(s[s$domain == "L2/3" & s$d_L1 > 300, ], "L2/3 spots far from L1 lines (reference)"))
  cat("\n== ", sec, " — domain composition inside the annotated L1 zone ==\n"); print(as.data.frame(zone_rows[[sec]] %>% filter(grepl("L1 zone", zone)) %>% arrange(-n)))
  cat("\n== ", sec, " — marker profile (UMI per 10k) ==\n"); print(marker_rows[[sec]], row.names = FALSE)

  # ---- figures: H&E + spots + lines; zoom on the L1 zone -----------------------------------------
  fr <- fromJSON(file.path(HEP, "HE_frames.json"))[[sec]]
  img <- readPNG(file.path(HEP, paste0(sec, "_HE_fig6b_frame.png")))
  ras <- annotation_custom(grid::rasterGrob(img, interpolate = TRUE, width = unit(1, "npc"), height = unit(1, "npc")), xmin = fr$x_min, xmax = fr$x_max, ymin = fr$y_min, ymax = fr$y_max)
  s$domain <- factor(s$domain, levels = c(DOM_LEV, setdiff(unique(s$domain), DOM_LEV)))
  pal <- c(DOM_PAL, setNames(rep("black", length(setdiff(levels(s$domain), DOM_LEV))), setdiff(levels(s$domain), DOM_LEV)))
  base <- function(p, sz) p + ras + geom_point(aes(x, y, colour = domain), size = sz, stroke = 0, alpha = 0.85) +
    geom_point(data = d, aes(x, y), colour = ifelse(d$label == "L1", "black", "grey20"), size = 0.25, inherit.aes = FALSE) +
    scale_colour_manual(values = pal, name = NULL, drop = TRUE) + coord_equal() + theme_void(base_size = 9) + theme(legend.position = "bottom")
  p_full <- base(ggplot(s), 1.3) + labs(title = paste(sec, "— domains over H&E; dashed lines = pathologists' annotation (black L1, grey WM)"))
  zone_sp <- s[s$d_L1 <= 250, ]; if (nrow(zone_sp) > 10) {
    xl <- range(zone_sp$x) + c(-15, 15); yl <- range(zone_sp$y) + c(-15, 15)
    p_zoom <- base(ggplot(s), 3.2) + coord_equal(xlim = xl, ylim = yl, expand = FALSE) + labs(title = "zoom: annotated L1 band (all domains)")
    sub <- s %>% filter(d_L1 <= 250, domain %in% c("WM","L1/pia","Vasc/immune"))
    p_zoom2 <- ggplot(sub) + ras + geom_point(aes(x, y, colour = domain), size = 3.6, stroke = 0) +
      geom_point(data = d[d$label == "L1", ], aes(x, y), colour = "black", size = 0.3, inherit.aes = FALSE) +
      geom_point(data = d[d$label == "WM", ], aes(x, y), colour = "grey30", size = 0.3, inherit.aes = FALSE) +
      scale_colour_manual(values = DOM_PAL, name = NULL) + coord_equal(xlim = xl, ylim = yl, expand = FALSE) + theme_void(base_size = 9) +
      theme(legend.position = "bottom") + labs(title = "zoom: only WM / L1-pia / Vasc-immune spots")
    ggsave(file.path(OUTD, paste0("sulcus_check_", sec, ".png")), p_full / (p_zoom | p_zoom2), width = 12, height = 15, dpi = 150, bg = "white")
    ggsave(file.path(OUTD, paste0("sulcus_check_", sec, "_L1zone.png")), p_zoom2, width = 8, height = 8, dpi = 200, bg = "white")
  } else ggsave(file.path(OUTD, paste0("sulcus_check_", sec, ".png")), p_full, width = 9, height = 9, dpi = 150, bg = "white")
}
write.csv(bind_rows(zone_rows), file.path(TAB, "visium_sulcus_check_zone_composition_FH.csv"), row.names = FALSE)
write.csv(bind_rows(marker_rows), file.path(TAB, "visium_sulcus_check_markers_FH.csv"), row.names = FALSE)
write.csv(md[, c("spot_id","sample_id","seurat_clusters","domain","x","y","d_L1","d_WM","zone", names(SETS))], file.path(TAB, "visium_sulcus_check_perspot_FH.csv"), row.names = FALSE)
cat("\n=== DONE ===\n")
