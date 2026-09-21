# assemble_figure_Visium.R — Apparent-type lift: all text sizes x1.25. These panels are placed at scale 0.43-0.59 in the assembled Figure 6, so their 7-9 pt type was reading at only
# 3.5-4.1 pt on the page -- the lowest in the whole figure set. Point sizes, line
# widths and unit() dimensions are untouched; canvas sizes unchanged.
# !/usr/bin/env Rscript
# NHD manuscript — VISIUM "last figure" composite (NHD house style: no bold,
# italic gene symbols, CON=#B0BEC5 / NHD=#FFB27A, clean bordered panels).
# Panels: a QC violins | b integrated UMAP | c annotated spatial clusters |
#         d cluster marker dotplot | e per-condition canonical-marker depth
#         profile (CON ordered vs NHD flat = headline) | f marker+ prevalence.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork); library(dplyr) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
BASE <- file.path(dirname(Sys.getenv("NHD_PROJ")), "Visium")
stopifnot(!is.na(BASE), dir.exists(BASE))
IH <- file.path(BASE,"integrated_harmony"); PC <- file.path(BASE,"per_condition_layers")
OUT <- file.path(BASE,"figure_Visium"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
# Condition palette = the microglia module-violin pair used across all condition
# UMAPs (CON pale blue / NHD pale red)
# the pair is read from the theme (PAL_COND_PALE = dense-UMAP/violin tint), never typed here
source(file.path(dirname(BASE), "NHD_frontal_hippo_rebuild", "scripts", "22_publication_theme_FH.R"))
PAL_COND <- PAL_COND_PALE
# QC violins: colour per donor/section but within the house condition families so the
# two CON sections read as one family and the two NHD sections as another (no off-palette
# drift). Shades derived from PAL_COND: lighter/darker of #96BFDB (CON) and #DF9C95 (NHD).
PAL_DONOR <- c(CON_Frontal1 = colorspace::lighten(PAL_COND[["CON"]], 0.25), CON_Frontal2 = colorspace::darken(PAL_COND[["CON"]], 0.20),
               NHD_Frontal1 = colorspace::lighten(PAL_COND[["NHD"]], 0.25), NHD_Frontal2 = colorspace::darken(PAL_COND[["NHD"]], 0.20))   # derived, not typed
th <- theme_bw(base_size=11.2) + theme(
  text=element_text(family="Helvetica", face="plain"),
  plot.title=element_text(size=11.2, hjust=0), panel.grid=element_blank(),
  axis.text=element_text(size=8.8), legend.key.size=unit(.32,"cm"), legend.text=element_text(size=8.8))
ital <- element_text(size=8.1, face="italic", family="Helvetica")

obj <- readRDS(file.path(IH,"NHD_frontal_integrated_harmony.rds"))
# Supplement keeps cluster-level labels (N:identity, one row per cluster of the current integration)
# read from integrated_cluster_identity.csv; the cluster count is taken from the file, never typed.
idc <- read.csv(file.path(IH,"integrated_cluster_identity.csv"))
idc$lab <- paste0(sub("^g","",idc$cluster),":",idc$identity)
lab_map <- setNames(idc$lab, sub("^g","",idc$cluster))
obj$cl_lab <- factor(unname(lab_map[as.character(obj$seurat_clusters)]),
                     levels=unname(lab_map[order(as.numeric(names(lab_map)))]))
# Main figure: collapse same-identity clusters into named domains; no cluster numbers.
# The cluster -> domain map is derived from the same identity CSV, so it follows the object's
# cluster numbering by construction and is never hard-coded here.
DOM_MAP <- setNames(as.character(idc$identity), sub("^g","",idc$cluster))
# Domain order and palette are read from integrated_domain_levels.csv (the same file the
# NHD_frontal_hippo_rebuild scripts read through _visium_domains_FH.R), never typed here
dlv <- read.csv(file.path(IH,"integrated_domain_levels.csv"), stringsAsFactors=FALSE)
dlv <- dlv[order(dlv$order),]
DOM_LEVELS <- dlv$domain
PAL_DOM <- setNames(dlv$colour, DOM_LEVELS)   # = Fig 6 DOM_PAL
stopifnot("identity CSV carries a domain outside DOM_LEVELS" = all(DOM_MAP %in% DOM_LEVELS),
          "object has a cluster with no identity row — refresh integrated_cluster_identity.csv" =
            all(as.character(unique(obj$seurat_clusters)) %in% names(DOM_MAP)))
cat("clusters:", nlevels(obj$cl_lab), "| spots:", ncol(obj), "\n"); print(table(obj$cl_lab, obj$condition))
obj$domain <- factor(unname(DOM_MAP[as.character(obj$seurat_clusters)]), levels=DOM_LEVELS)
# cluster palette = the Fig-6 domain hue of each cluster, shaded when a domain spans several
# clusters (e.g. three L6 clusters = three reds), so the supplement reads on the main figure's key
.shade <- function(col, k, n) { if (n == 1) return(col); colorspace_ok <- requireNamespace("colorspace", quietly = TRUE)
  f <- seq(-0.28, 0.28, length.out = n)[k]; if (colorspace_ok) colorspace::lighten(col, f) else col }
PAL_CL <- vapply(seq_len(nrow(idc)), function(i) { dom <- idc$identity[i]; sib <- which(idc$identity == dom)
  .shade(PAL_DOM[[dom]], match(i, sib), length(sib)) }, character(1))
names(PAL_CL) <- idc$lab; PAL_CL <- PAL_CL[levels(obj$cl_lab)]
# Clean display NAMES for the section/sample tick + strip labels (don't leak the
# internal underscore token "CON_Frontal1" into visible axis text / facet titles).
SEC_LEVELS <- c("CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2")
SEC_LAB    <- c(CON_Frontal1="CON Frontal 1", CON_Frontal2="CON Frontal 2",
                NHD_Frontal1="NHD Frontal 1", NHD_Frontal2="NHD Frontal 2")
.sd <- factor(SEC_LAB[as.character(obj$sample_id)], levels=unname(SEC_LAB[SEC_LEVELS]))
names(.sd) <- colnames(obj)                        # name by cell barcode for AddMetaData alignment
obj$sample_disp <- .sd
md <- obj@meta.data

## a) QC violins (UMIs / Features / %Mito), CON vs NHD
qcv <- function(y,lab,logy=FALSE){ p<-ggplot(md, aes(sample_disp,.data[[y]],fill=sample_id))+
  geom_violin(scale="width",linewidth=.25)+geom_boxplot(width=.12,outlier.size=.15,fill="white",linewidth=.25)+
  scale_fill_manual(values=PAL_DONOR)+th+labs(x=NULL,y=lab)+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=7.5),legend.position="none"); if(logy)p<-p+scale_y_log10(); p }
pa <- qcv("nCount_Spatial","UMIs",TRUE)/qcv("nFeature_Spatial","Genes",TRUE)/qcv("percent.mt","% Mito")

## b) integrated UMAP (domain + condition) — both sub-UMAPs share one builder.
##    Legends at BOTTOM (not right) so each blob uses the full panel width, and
##    near-zero axis expansion so the blob nearly touches the frame; coord_equal
##    keeps the embedding's natural ~1.47:1 (x14.2 : y9.6) shape (no distortion).
##    This kills the large top/bottom whitespace bands the right-legend+coord_equal
##    squeeze produced. Manual ggplot (not DimPlot) so we can shuffle the draw order
##    for fair overplotting (DimPlot draws in factor order, hiding lower levels).
.red2 <- if ("umap" %in% Reductions(obj)) "umap" else tail(Reductions(obj), 1)
.emb2 <- as.data.frame(Embeddings(obj, .red2))[, 1:2]; colnames(.emb2) <- c("umap_1","umap_2")
.umap <- data.frame(.emb2, domain=obj$domain, cl_lab=obj$cl_lab,
                    condition=factor(obj$condition, levels=names(PAL_COND)))
umap_panel <- function(data, grp, pal, legend_nrow=1, pt=.35){
  # Thin dark outline on every dot: shape 21 fillable circle,
  # data colour on `fill`, grey20 border via colour + thin stroke 0.04.
  ggplot(data, aes(umap_1, umap_2, fill=.data[[grp]])) +
    geom_point(size=pt, stroke=0.04, shape=21, colour="grey20") +
    scale_fill_manual(values=pal, name=NULL) +
    scale_x_continuous(name="UMAP 1", expand=expansion(mult=0.02)) +   # presentation label, not raw umap_1 token
    scale_y_continuous(name="UMAP 2", expand=expansion(mult=0.02)) +
    coord_equal(clip="off") + th + ggtitle(NULL) +
    theme(legend.position="bottom", legend.margin=margin(t=-3,b=0),
          legend.box.spacing=unit(2,"pt"), legend.spacing.x=unit(2,"pt"),
          plot.margin=margin(1,2,1,2)) +
    guides(fill=guide_legend(nrow=legend_nrow, byrow=TRUE,
                             override.aes=list(shape=21, colour="grey20", stroke=0.3, size=2)))
}
# Set.seed(43) before each shuffle -> the CON/NHD draw order of pb2 (and the domain draw order
# of pb1) is a reproducible permutation of the object's rows.
set.seed(43); pb1 <- umap_panel(.umap[sample(nrow(.umap)), ], "domain", PAL_DOM, legend_nrow=2)
set.seed(43); pb2 <- umap_panel(.umap[sample(nrow(.umap)), ], "condition", PAL_COND, legend_nrow=1)

## c) annotated spatial clusters on tissue (4 samples) -- main composite only (Fig 6b itself is script 74)
sp <- SpatialDimPlot(obj, group.by="domain", combine=FALSE, pt.size.factor=1.5, image.alpha=0.55, stroke=0)
sp <- lapply(seq_along(sp), function(i) sp[[i]] + scale_fill_manual(values=PAL_DOM) +
  ggtitle(unname(SEC_LAB[levels(factor(obj$sample_id))[i]])) +
  theme(legend.position="none", plot.title=element_text(size=10.0,hjust=.5)))  # paired w/ UMAP -> no legend
pc <- wrap_plots(sp, nrow=1)   # main figure: spatial in one row

## d) spatial-domain marker-annotation dotplot — SCZ Fig1-e style: marker genes
##    faceted by the cell-type/layer they mark (group headers on top); size = %
##    expressing, fill = z-scored mean expression. Shows the canonical-marker
##    basis of the domain annotation (transparency/QC), not independent validation
##    (domains were annotated using these markers). NEFH = L5 not L4.
MARK_GROUPS <- list(
  "L1/pia"=c("DCN","COL1A2","SLC6A13"),   # leptomeningeal / VLMC genes (Yang 2022; Kearns 2023); AQP4 sits under Astro
  "L2/3"=c("CUX2","LAMP5"), "L4"=c("RORB"), "L5"=c("PCP4","FEZF2","BCL11B","NEFH"),
  "L6"=c("TLE4","FOXP2"), "ExN"=c("RBFOX3","SLC17A7","SNAP25"), "InN"=c("GAD1","GAD2"),
  "WM"=c("MBP","PLP1","MOBP","MOG"), "OPC"=c("PDGFRA","OLIG1"),
  "Astro"=c("GFAP","AQP4","SLC1A2"), "Micro"=c("TYROBP","C1QB","P2RY12"),
  "Vasc"=c("CLDN5","PECAM1","PDGFRB"))
MARK_GROUPS <- lapply(MARK_GROUPS, function(g) intersect(g, rownames(obj)))
MARK_GROUPS <- MARK_GROUPS[lengths(MARK_GROUPS) > 0]
g2grp <- stack(MARK_GROUPS); g2grp <- setNames(as.character(g2grp$ind), g2grp$values)
mk_genes <- unlist(MARK_GROUPS, use.names=FALSE)
dd <- DotPlot(obj, features=mk_genes, group.by="domain")$data
dd$group <- factor(g2grp[as.character(dd$features.plot)], levels=names(MARK_GROUPS))
dd$gene  <- factor(dd$features.plot, levels=mk_genes)
dd$z     <- pmax(pmin(dd$avg.exp.scaled, 2), -1)
pd <- ggplot(dd, aes(gene, id)) +
  geom_point(aes(size=pct.exp, fill=z), shape=21, colour="grey30", stroke=0.2) +
  facet_grid(. ~ group, scales="free_x", space="free_x") +   # group headers at TOP
  scale_fill_gradient2(low="#4575B4", mid="grey95", high="#D73027", midpoint=0,
                       limits=c(-1,2), name="scaled\nexpr.") +
  scale_size_continuous(range=c(0.3,4), breaks=c(25,50,75), name="% expr.") +
  labs(x=NULL, y=NULL) + th +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=7.5,face="italic"),
        axis.text.y=element_text(size=10.0),
        strip.text.x=element_text(size=8.1,face="plain"), strip.placement="outside",   # no-bold house rule
        strip.background=element_rect(fill="grey92", colour=NA),   # pale, borderless banner as in Fig 6d
        panel.spacing.x=unit(2,"pt"), panel.border=element_rect(color="grey88",fill=NA,linewidth=.3))

## e) laminar ORGANIZATION (strong, depth-axis-free): Moran's I of layer markers
##    (banded in CON -> dispersed in NHD; MBP/GM-WM boundary preserved = positive
##    control) + neighbourhood domain coherence. Computed in laminar_organization.R.
mo <- read.csv(file.path(PC,"laminar_moran.csv")); en <- read.csv(file.path(PC,"laminar_entropy.csv"))
ml <- mo[mo$class=="laminar",]; ml$gene <- factor(ml$gene, levels=c("CUX2","RORB","PCP4","FEZF2","TLE4"))
mc <- mo[mo$gene=="MBP",]                                   # GM/WM boundary positive control
pe1 <- ggplot(ml, aes(gene, I_resid, color=condition)) +
  geom_boxplot(aes(group=interaction(gene,condition)), width=.55, outlier.shape=NA, linewidth=.3, position=position_dodge(.62)) +
  geom_point(position=position_dodge(.62), size=1.4) +
  geom_hline(yintercept=0, linetype=2, color="grey75", linewidth=.3) +
  scale_color_manual(values=PAL_COND) + th + theme(legend.position="top") +
  labs(x=NULL, y="layer-marker Moran's I\n(UMI-residualized)")
pe2 <- ggplot(en, aes(condition, coherence_gap, color=condition)) +
  geom_boxplot(width=.45, outlier.shape=NA, linewidth=.3) + geom_jitter(width=.07, size=1.8) +
  scale_color_manual(values=PAL_COND) + th + theme(legend.position="none") +
  labs(x=NULL, y="domain coherence gap\n(null - obs entropy)")
pe <- (pe1 | pe2) + plot_layout(widths=c(1.9,1))

## e (support) — updated marker-vs-depth gradient: CON depth-ordered -> NHD flat (intuitive companion to the Moran's I)
dep <- rbind(transform(read.csv(file.path(PC,"CON_marker_depth.csv")), cond="CON"),
             transform(read.csv(file.path(PC,"NHD_marker_depth.csv")), cond="NHD"))
gms <- c("CUX2","RORB","PCP4","TLE4","MBP")   # L2/3,L4,L5,L6,WM (drop noisy FEZF2)
depz <- do.call(rbind, lapply(gms, function(g) data.frame(cond=dep$cond, depth=dep$depth, gene=g,
          z=ave(dep[[g]], dep$cond, FUN=function(x) scale(x)[,1]))))   # z within condition
depz$gene <- factor(depz$gene, levels=gms)
# label the depth axis by cortical LAYER (data-driven marker-peak positions), not numbers
LAY_BREAKS <- c(-1.8,-0.9,0.0,0.9,2.2); LAY_LABS <- c("L2/3","L4","L5","L6","WM")
pe_depth <- ggplot(depz, aes(depth, z, color=gene)) +
  geom_hline(yintercept=0, linetype=2, color="grey75", linewidth=.3) +
  geom_smooth(se=FALSE, method="loess", span=0.7, formula=y~x, linewidth=.9) +
  facet_wrap(~cond, nrow=1) + th + scale_color_brewer(palette="Dark2") +
  scale_x_continuous(breaks=LAY_BREAKS, labels=LAY_LABS) +
  labs(x="cortical layer (superficial -> deep)", y="z(expr)")   # ASCII arrow (no unicode in base-pdf text)

## f) marker+ prevalence CON vs NHD
fr <- rbind(read.csv(file.path(PC,"CON_marker_prevalence.csv")),
            read.csv(file.path(PC,"NHD_marker_prevalence.csv")))
fr$cond <- ifelse(grepl("^CON",fr$sample_id),"CON","NHD")
prev <- fr %>% tidyr::pivot_longer(c(RORB_pos,CUX2_pos,PCP4_pos,TLE4_pos,MBP_pos), names_to="marker", values_to="frac") %>%
  group_by(cond,marker) %>% summarise(frac=mean(frac), .groups="drop")
prev$marker <- sub("_pos","",prev$marker)
pf <- ggplot(prev, aes(marker,frac,fill=cond)) + geom_col(position="dodge", width=.7) +
  scale_fill_manual(values=PAL_COND, name="Condition") + th + labs(x=NULL,y="fraction spots +") +   # not raw "cond" token
  theme(axis.text.x=ital)

## ---- apparent-type parity for supplementary figure 4 -----------
## House rule: on-page pt = declared pt x placement scale (placed_width/saved_width);
## main-figure band ~4.8-5.3 pt, target 5.1. Placements measured from the assembled
## figures/Supplementary/SuppFig4_Visium_QC/Fig4_suppl.pdf:
##     a  4.00 -> 2.08 in (0.521)   b  4.60 -> 2.08 in (0.453)
##     c  6.00 -> 2.08 in (0.347)   d 10.00 -> 6.94 in (0.694)
## `th` above is tuned for the main composite and is not touched here -- each SUPP
## export gets its own text-only multiplier, so Figure_Visium_composite and every
## main panel_*.{png,pdf} render byte-identically to before. text only: geom_point
## size, linewidth and unit() dimensions are never multiplied.
SUPP_F <- c(a = 1.11, b = 1.28, d = 0.86)   # panel c handled separately (canvas rebuilt)
supp_type <- function(f) theme(
  text         = element_text(size = 11.2 * f),
  plot.title   = element_text(size = 11.2 * f, hjust = 0),
  axis.title   = element_text(size = 11.2 * f),
  axis.text    = element_text(size =  8.8 * f),
  legend.title = element_text(size = 11.2 * f),
  legend.text  = element_text(size =  8.8 * f))

## ---- composite ----
top <- (pa | (pb1/pb2)) + plot_layout(widths=c(1,1.6))
botrow <- (pe1 | pe2 | pf) + plot_layout(widths=c(1.7,0.9,1.1))   # laminar Moran's I | coherence | prevalence
fig <- top / pc / pd / botrow + plot_layout(heights=c(1.28,1.0,0.63,1.2))  # a/b 20% smaller; c one row; d 30%+10% less tall
# VISIUM_SUPP_ONLY=1 skips the two slow 14x17 in composite renders (and, below, the
# main per-panel exports) so the supplement panels can be iterated in ~1 min instead
# of ~8. Unset/0 => the full original behaviour, byte-for-byte.
SUPP_ONLY <- identical(Sys.getenv("VISIUM_SUPP_ONLY"), "1")
if (!SUPP_ONLY) {
ggsave(file.path(OUT,"Figure_Visium_composite.png"), fig, width=14, height=17, dpi=600, limitsize=FALSE, device=ragg::agg_png)
ggsave(file.path(OUT,"Figure_Visium_composite.pdf"), fig, width=14, height=17, limitsize=FALSE)
}
## ---- per-panel editable files: every panel (a-f) as PNG + PDF ----
pb <- pb1 / pb2                       # panel b = clusters UMAP over condition UMAP
# panel_a_QC is consumed only by Supplementary Figure 4 (verified: no image with its
# 2400x3360 px footprint is placed in figures/Figure_6/panels/Fig6.pdf), so it carries
# the SUPP type multiplier while the composite's own `pa` stays untouched.
pa_s <- pa & supp_type(SUPP_F[["a"]]) &
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7.5*SUPP_F[["a"]]))
panels <- list(panel_a_QC=list(p=pa_s, w=4,   h=5.6),   # 20% smaller
               panel_b_UMAP=list(p=pb, w=3.9, h=4.7),  # bottom legends + tight expand -> blobs fill frame, no white bands
               panel_c_spatial=list(p=pc, w=11.2, h=3.2),  # main: spatial in one row, 20% smaller -> bigger fonts
               panel_d_dotplot=list(p=pd, w=10, h=2.52),  # 30% + a further 10% less tall
               panel_e_laminar=list(p=pe, w=8,  h=3.0),   # strong, depth-free (Moran's I + coherence)
               panel_e_depth=list(p=pe_depth, w=8, h=2.30),  # layer-depth gradient, 20%+10% less tall
               panel_f_prevalence=list(p=pf, w=4, h=3))
for (nm in names(panels)) {
  # panel_a_QC is a supplement panel -> always written; the rest are main-only.
  if (SUPP_ONLY && nm != "panel_a_QC") next
  el <- panels[[nm]]
  ggsave(file.path(OUT, paste0(nm,".png")), el$p, width=el$w, height=el$h, dpi=600, limitsize=FALSE, device=ragg::agg_png)
  ggsave(file.path(OUT, paste0(nm,".pdf")), el$p, width=el$w, height=el$h, limitsize=FALSE)
}

## ---- supplement panels: CLUSTER-LEVEL (N:identity, every cluster of the current integration) ----
set.seed(43); pb1c <- umap_panel(.umap[sample(nrow(.umap)), ], "cl_lab", PAL_CL, legend_nrow=3)
pbc  <- pb1c / pb2
# Panel c (supp) -- drawn in the Fig-6b FRAME (script 74): spot coordinates from
# cell2location/c2l_MAIN/spot_coords.csv (x = pxl_row-based, y = pxl_col-based, the frame
# every Figure-6 map uses), each section centred and scaled by its own larger extent so
# the four tiles share one unit square, coord_equal. Colour = the cluster palette above
# (domain hue shaded per sibling cluster). No H&E underlay: the tissue image is Fig 6b's
# top row, and the supplement's job is the cluster-level map on the same frame.
C2L_DIR <- file.path(BASE, "cell2location", "c2l_MAIN")
xy <- read.csv(file.path(C2L_DIR, "spot_coords.csv"), stringsAsFactors=FALSE)
stopifnot("spot_coords.csv does not cover the object's spots" = mean(colnames(obj) %in% xy$cell) > 0.99)
spc_d <- data.frame(spot_id=colnames(obj), cl_lab=obj$cl_lab, sample_disp=obj$sample_disp, sample_id=as.character(obj$sample_id))
spc_d <- merge(spc_d, data.frame(spot_id=xy$cell, x=xy$x, y=xy$y), by="spot_id")
spc_d <- do.call(rbind, lapply(split(spc_d, spc_d$sample_id), function(z) {
  sp <- max(diff(range(z$x)), diff(range(z$y)))
  z$xr <- (z$x - mean(range(z$x))) / sp; z$yr <- (z$y - mean(range(z$y))) / sp; z }))
cat("supp panel c spots per section (Fig-6b frame):\n"); print(table(spc_d$sample_disp))
pcc <- ggplot(spc_d, aes(xr, yr, colour=cl_lab)) +
  geom_point(size=0.85, stroke=0) +   # closes the hex grid at the 1.7-in facet width (dots were too small to see)
  scale_colour_manual(values=PAL_CL, name=NULL, drop=FALSE) +
  facet_wrap(~sample_disp, nrow=2) + coord_equal() +
  guides(colour=guide_legend(nrow=3, byrow=TRUE, override.aes=list(size=2.2))) +
  th + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank(),
             axis.line=element_blank(), panel.border=element_blank(),
             strip.background=element_blank(), strip.text=element_text(size=9.0, colour="black"),
             panel.spacing=unit(2,"pt"), plot.margin=margin(2,1,0,1),
             legend.position="bottom", legend.text=element_text(size=8.8),
             legend.key.size=unit(0.30,"cm"), legend.margin=margin(t=-2,b=0),
             legend.box.spacing=unit(2,"pt"), legend.spacing.x=unit(2,"pt"))
ddc <- DotPlot(obj, features=mk_genes, group.by="cl_lab")$data
ddc$group <- factor(g2grp[as.character(ddc$features.plot)], levels=names(MARK_GROUPS))
ddc$gene  <- factor(ddc$features.plot, levels=mk_genes); ddc$z <- pmax(pmin(ddc$avg.exp.scaled,2),-1)
pdc <- ggplot(ddc, aes(gene, id)) + geom_point(aes(size=pct.exp, fill=z), shape=21, colour="grey30", stroke=0.2) +
  facet_grid(. ~ group, scales="free_x", space="free_x") +
  scale_fill_gradient2(low="#4575B4", mid="grey95", high="#D73027", midpoint=0, limits=c(-1,2), name="scaled\nexpr.") +
  scale_size_continuous(range=c(0.3,4), breaks=c(25,50,75), name="% expr.") + labs(x=NULL,y=NULL) + th +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=7.5,face="italic"), axis.text.y=element_text(size=10.0),
        strip.text.x=element_text(size=8.1,face="plain"), strip.placement="outside",   # no-bold house rule
        strip.background=element_rect(fill="grey92", colour=NA),   # pale, borderless banner as in Fig 6d
        panel.spacing.x=unit(2,"pt"), panel.border=element_rect(color="grey88",fill=NA,linewidth=.3))
# apparent-type lift, text only (see SUPP_F above). Legend keys are TRIMMED (not
# scaled) so the bigger type does not push the legend block off the fixed b canvas.
pbc <- pbc & supp_type(SUPP_F[["b"]]) &
  # axis.title pinned to the axis.text size (flat hierarchy) rather than the theme
  # base: on a blank-ish UMAP the 11.2*1.28 = 14.3 pt title read 6.5 pt on the page
  # (out of the 4.8-5.3 band) and stole the height that coord_equal converts into
  # panel width. At 8.8*1.28 it lands at 5.10 pt and the blobs get wider.
  theme(axis.title=element_text(size=8.8*SUPP_F[["b"]]),
        legend.key.size=unit(0.24,"cm"), legend.margin=margin(t=-3,b=0),
        legend.box.spacing=unit(1,"pt"), legend.spacing.x=unit(1,"pt"),
        plot.margin=margin(3,2,1,2))   # top pad: the "5.0" y tick was flush to row 0
pdc <- pdc + supp_type(SUPP_F[["d"]]) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=7.5*SUPP_F[["d"]],face="italic"),
        axis.text.y=element_text(size=10.0*SUPP_F[["d"]]),
        strip.text.x=element_text(size=8.1*SUPP_F[["d"]],face="plain"),
        legend.key.size=unit(0.30,"cm"))
PANEL_C_W <- 3.60; PANEL_C_H <- 4.55
for (e in list(list(p=pbc,n="panel_b_UMAP_clusters",w=4.6,h=5.2),
               list(p=pcc,n="panel_c_spatial_clusters",w=PANEL_C_W,h=PANEL_C_H),
               list(p=pdc,n="panel_d_dotplot_clusters",w=10,h=2.8))) {
  ggsave(file.path(OUT, paste0(e$n,".png")), e$p, width=e$w, height=e$h, dpi=600, limitsize=FALSE, device=ragg::agg_png)
  ggsave(file.path(OUT, paste0(e$n,".pdf")), e$p, width=e$w, height=e$h, limitsize=FALSE)
}
cat("\n=== DONE ===\n", file=stderr())
cat("== MAIN (collapsed-named): Figure_Visium_composite + panel_b_UMAP/panel_c_spatial/panel_d_dotplot\n")
cat("== SUPP (cluster-level): panel_b_UMAP_clusters/panel_c_spatial_clusters/panel_d_dotplot_clusters  in:", OUT, "\n")
