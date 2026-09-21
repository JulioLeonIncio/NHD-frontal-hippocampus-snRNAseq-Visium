#!/usr/bin/env Rscript
# =============================================================================
# 116b_L4_label_transfer_DFC_FH.R — Re-annotate frontal excitatory nuclei against the Jorstad 2023 dorsolateral prefrontal cortex reference (has an L4 IT class).
# Why: Fig 5 / Table subclasses are Azimuth vs the human motor cortex (no L4 IT). Visium
#   shows an L4 domain in CON (1,546 spots) and almost none in NHD (118). Question: does the
#   single-nucleus data contain L4 IT neurons, and is their share lower in NHD?
# Method: reference = data/_reference_Jorstad2023_DFC (116a; excitatory, WithinArea_subclass
#   capped at 6,000/class). Both sides LogNormalize on shared genes; reference PCA (30);
#   Seurat FindTransferAnchors (project.query = FALSE, reference.reduction = "pca") ->
#   TransferData for WithinArea_subclass and CrossArea_subclass. Query = data/_cache_116_frontalEx_FH.rds
#   (frontal Neuron_Ex, 116). Outputs the per-nucleus prediction + score, the CON/NHD
#   share of every subclass per lane, the overlap with the 116 clusters and the Azimuth
#   labels, and a marker sanity check (RORB / CUX2 / L4-vs-L5 IT DE genes from the reference).
# Outputs:
#   tables/L4_DFC_transfer_percell_FH.csv, tables/L4_DFC_transfer_summary_FH.csv,
#   tables/L4_DFC_reference_L4vsL5_markers_FH.csv
#   figures/_diagnostics/L4_neuron_search/L4_DFC_transfer.png
#   logs/116b_L4_label_transfer_DFC_FH.log
# Run: OMP_NUM_THREADS=4 KMP_DUPLICATE_LIB_OK=TRUE Rscript scripts/116b_L4_label_transfer_DFC_FH.R
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr); library(ggplot2); library(patchwork) })
set.seed(42); Sys.setenv(OMP_NUM_THREADS = "4", KMP_DUPLICATE_LIB_OK = "TRUE")
options(future.globals.maxSize = 16 * 1024^3)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
REFD  <- file.path(PROJ, "data", "_reference_Jorstad2023_DFC")
QRY   <- file.path(PROJ, "data", "_cache_116_frontalEx_FH.rds")
PERC  <- file.path(PROJ, "tables", "L4_neuron_search_percell_FH.csv")
OUTD  <- file.path(PROJ, "figures", "_diagnostics", "L4_neuron_search"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
LOG   <- file.path(PROJ, "logs", "116b_L4_label_transfer_DFC_FH.log"); logcon <- file(LOG, "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
say("== 116b == ", format(Sys.time()))
stopifnot(file.exists(file.path(REFD, "matrix.mtx.gz")), file.exists(QRY))

# ---- reference ------------------------------------------------------------------
rm_ <- ReadMtx(file.path(REFD, "matrix.mtx.gz"), file.path(REFD, "barcodes.tsv.gz"), file.path(REFD, "features.tsv.gz"),
               feature.column = 2, skip.feature = 0, skip.cell = 0)
rmeta <- read.csv(file.path(REFD, "meta.csv"), row.names = 1)
rm_ <- rm_[!duplicated(rownames(rm_)), ]
ref <- CreateSeuratObject(rm_, meta.data = rmeta[colnames(rm_), ], min.cells = 10)
say("reference: ", ncol(ref), " nuclei x ", nrow(ref), " genes; WithinArea_subclass: ",
    paste(names(table(ref$WithinArea_subclass)), table(ref$WithinArea_subclass), sep = "=", collapse = " "))

# ---- query ----------------------------------------------------------------------
q <- readRDS(QRY); q <- JoinLayers(q)
shared <- intersect(rownames(ref), rownames(q)); say("shared genes: ", length(shared))
ref <- subset(ref, features = shared); q <- subset(q, features = shared)
ref <- NormalizeData(ref, verbose = FALSE) |> FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
q <- NormalizeData(q, verbose = FALSE)
anch <- FindTransferAnchors(reference = ref, query = q, normalization.method = "LogNormalize",
                            reference.reduction = "pca", dims = 1:30, features = VariableFeatures(ref), verbose = FALSE)
say("anchors: ", nrow(anch@anchors))
pw <- TransferData(anch, refdata = ref$WithinArea_subclass, dims = 1:30, verbose = FALSE)
pc <- TransferData(anch, refdata = ref$CrossArea_subclass, dims = 1:30, verbose = FALSE)
q$dfc_within <- pw$predicted.id; q$dfc_within_score <- pw$prediction.score.max
q$dfc_cross  <- pc$predicted.id; q$dfc_cross_score  <- pc$prediction.score.max

# ---- join the 116 clusters / UMAP ----------------------------------------------------
per <- read.csv(PERC); rownames(per) <- per$barcode
q$cl116 <- per[colnames(q), "cl"]; q$umap_1 <- per[colnames(q), "umap_1"]; q$umap_2 <- per[colnames(q), "umap_2"]
md <- q@meta.data

# ---- summaries ------------------------------------------------------------------
share_tab <- function(lab) {
  a <- md %>% count(Condition, .data[[lab]]) %>% group_by(Condition) %>% mutate(share = n / sum(n)) %>% ungroup() %>%
    rename(subclass = all_of(lab)) %>% tidyr::pivot_wider(names_from = Condition, values_from = c(n, share), values_fill = 0) %>%
    mutate(log2_share_NHD_vs_CON = round(log2((share_NHD + 1e-4) / (share_CON + 1e-4)), 2), taxonomy = lab)
  ln <- md %>% count(SampleID, Condition, .data[[lab]]) %>% group_by(SampleID) %>% mutate(share = round(n / sum(n), 4)) %>% ungroup() %>%
    rename(subclass = all_of(lab)) %>% select(subclass, SampleID, share) %>%
    tidyr::pivot_wider(names_from = SampleID, values_from = share, values_fill = 0, names_prefix = "lane_")
  left_join(a, ln, by = "subclass")
}
S <- bind_rows(share_tab("dfc_within"), share_tab("dfc_cross"))
say("\nDFC within-area taxonomy — share of each condition's frontal excitatory nuclei:")
print(as.data.frame(S %>% filter(taxonomy == "dfc_within") %>% select(subclass, n_CON, n_NHD, share_CON, share_NHD, log2_share_NHD_vs_CON, starts_with("lane_"))), digits = 3)
say("\nDFC cross-area taxonomy:")
print(as.data.frame(S %>% filter(taxonomy == "dfc_cross") %>% select(subclass, n_CON, n_NHD, share_CON, share_NHD, log2_share_NHD_vs_CON)), digits = 3)
say("\nmedian prediction score (within): ", round(median(md$dfc_within_score), 3), "; L4 IT: ",
    round(median(md$dfc_within_score[md$dfc_within == "L4 IT"]), 3), "; n L4 IT with score >= 0.5: ", sum(md$dfc_within == "L4 IT" & md$dfc_within_score >= 0.5))
say("\nL4 IT (within) x Azimuth motor-cortex label:"); print(table(md$predicted.subclass[md$dfc_within == "L4 IT"]))
say("\nL4 IT (within) x 116 cluster:"); print(table(md$cl116[md$dfc_within == "L4 IT"]))
say("\nDFC within x Azimuth (all):"); print(table(DFC = md$dfc_within, Azimuth = md$predicted.subclass))
write.csv(S, file.path(PROJ, "tables", "L4_DFC_transfer_summary_FH.csv"), row.names = FALSE)
write.csv(cbind(barcode = rownames(md), md[, c("Condition","SampleID","predicted.subclass","cl116","dfc_within","dfc_within_score","dfc_cross","dfc_cross_score")]),
          file.path(PROJ, "tables", "L4_DFC_transfer_percell_FH.csv"), row.names = FALSE)

# ---- marker sanity: L4 IT vs L5 IT genes in the reference, checked in our predicted classes ----
Idents(ref) <- "WithinArea_subclass"
mk <- FindMarkers(ref, ident.1 = "L4 IT", ident.2 = "L5 IT", only.pos = FALSE, logfc.threshold = 0.5, min.pct = 0.2, verbose = FALSE)
mk$gene <- rownames(mk); mk <- mk[order(-mk$avg_log2FC), ]
write.csv(mk, file.path(PROJ, "tables", "L4_DFC_reference_L4vsL5_markers_FH.csv"), row.names = FALSE)
top_up <- head(mk$gene[mk$avg_log2FC > 0 & mk$p_val_adj < 0.01], 8); top_dn <- head(rev(mk$gene[mk$avg_log2FC < 0 & mk$p_val_adj < 0.01]), 8)
say("reference L4 IT > L5 IT: ", paste(top_up, collapse = " "), " | L5 IT > L4 IT: ", paste(top_dn, collapse = " "))
Idents(q) <- "dfc_within"
genes <- unique(c("CUX2","LAMP5","RORB", top_up, top_dn, "THEMIS","FEZF2","SYT6","TSHZ2"))
genes <- genes[genes %in% rownames(q)]
pd <- DotPlot(q, features = genes, group.by = "dfc_within", assay = "RNA") + RotatedAxis() + labs(title = "our frontal excitatory nuclei by DFC-transferred subclass")
ggsave(file.path(OUTD, "L4_DFC_transfer_markers.png"), pd, width = 12, height = 5, dpi = 200, bg = "white")

# ---- panels ------------------------------------------------------------------------
emb <- md; emb$U1 <- emb$umap_1; emb$U2 <- emb$umap_2
p1 <- ggplot(emb, aes(U1, U2, colour = dfc_within)) + geom_point(size = 0.3, stroke = 0) + theme_pub() +
  labs(title = "Jorstad 2023 DFC within-area subclass") + guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1))
p2 <- ggplot(emb, aes(U1, U2, colour = dfc_within == "L4 IT")) + geom_point(size = 0.3, stroke = 0) +
  scale_colour_manual(values = c(`TRUE` = "#DD8452", `FALSE` = "grey85"), labels = c("other", "L4 IT")) + theme_pub() + labs(title = "L4 IT", colour = NULL)
p3 <- ggplot(emb, aes(U1, U2, colour = dfc_within_score)) + geom_point(size = 0.3, stroke = 0) + scale_colour_viridis_c() + theme_pub() + labs(title = "prediction score", colour = NULL)
sw <- S %>% filter(taxonomy == "dfc_within") %>% select(subclass, share_CON, share_NHD) %>% tidyr::pivot_longer(-subclass, names_to = "cond", values_to = "share") %>%
  mutate(cond = sub("share_", "", cond), subclass = factor(subclass, levels = c("L2/3 IT","L4 IT","L5 IT","L6 IT","L6 IT Car3","L5 ET","L5/6 NP","L6 CT","L6b")))
p4 <- ggplot(sw, aes(share, subclass, fill = cond)) + geom_col(position = position_dodge(width = 0.8), width = 0.7) + scale_fill_manual(values = PAL_COND) +
  scale_y_discrete(limits = rev) + theme_pub() + labs(x = "share of the condition's frontal excitatory nuclei", y = NULL, fill = NULL)
ln <- md %>% count(SampleID, Condition, dfc_within) %>% group_by(SampleID) %>% mutate(share = n / sum(n)) %>% ungroup() %>% filter(dfc_within == "L4 IT")
p5 <- ggplot(ln, aes(Condition, share, colour = Condition)) + geom_point(size = 2.5, position = position_jitter(width = 0.08, seed = 1)) +
  geom_text(aes(label = SampleID), size = 2.2, nudge_x = 0.22, colour = "black") + scale_colour_manual(values = PAL_COND) + theme_pub() + theme(legend.position = "none") +
  labs(y = "L4 IT share per lane", x = NULL, title = "L4 IT per lane")
ggsave(file.path(OUTD, "L4_DFC_transfer.png"), (p1 | p2 | p3) / (p4 | p5), width = 15, height = 9, dpi = 200, bg = "white")
saveRDS(list(anchors_n = nrow(anch@anchors), summary = S), file.path(PROJ, "data", "_cache_116b_L4_transfer_FH.rds"))
say("done ", format(Sys.time())); close(logcon)
