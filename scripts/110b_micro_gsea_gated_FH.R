#!/usr/bin/env Rscript
# =============================================================================
# 110b_micro_gsea_gated_FH.R — Microglial GO:BP + Reactome fgsea with the detection GATE, for Supplementary Fig. 3g.
# -----------------------------------------------------------------------------
# Why. ST3's frontal "lost in NHD microglia" terms (synaptic / calcium signalling) are carried by
# genes with pct_NHD ~ 0 and pct_CON 0.2-0.5 in the frontal lane (DSCAM, CACNA1A, CACNB4, LRRK2,
# PRKCB, ITPR1, PPP3CA): log2FC -3 to -8 on zero detection in 330 nuclei — the detection-floor
# pattern the project forbids colouring — and three of them sit on ARTIFACT_NEURONAL. The same
# genes are detected at 30-45 % in hippocampal NHD microglia and the terms run the other way
# there. So the ranking is rebuilt with the same recipe as script 44 (signed avg_log2FC
# winsorized 1/99 %, GO:BP + Reactome minus the tissue blacklist, minSize 10, maxSize 500,
# BH within comparison, seed 42) but restricted to genes that (i) are not is_artifact() for
# Micro-PVM and (ii) are detected in >= 10 % of nuclei in both conditions. Whatever survives is
# shown; if no "lost" term survives, the null is stated.
# Output: tables/micro_gsea_gated_FH.csv  (ST3 schema subset + gate provenance), read by 110.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(fgsea); library(msigdbr); library(data.table) })
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "_artifact_genes.R"))
source(file.path(PROJ, "scripts", "_gsea_composite_metric.R"))    # winsorize_lfc / fc_rank, as 44
say <- function(...) cat(sprintf(...), "\n")
MIN_PCT_BOTH <- 0.10

blacklist_keywords <- c("KERATINOCYTE","EPIDERM","HAIR_FOLLICLE","HAIR_CYCLE","TOOTH","ODONTOGENESIS",
  "AMELOBLAST","LENS","CORNIF","CARDIAC_MUSCLE","CARDIAC_CHAMBER","HEART_MORPHO","B_CELL_ACTIV",
  "THYMO","SPERMAT","OOCYTE","OVULATION","RETINA_LAYER","PHOTORECEPTOR_CELL_DIFF","ALCOHOL","ETHANOL")
go_bp    <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>% { split(.$gene_symbol, .$gs_name) }
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>% { split(.$gene_symbol, .$gs_name) }
sets <- c(go_bp, reactome); sets <- sets[!grepl(paste(blacklist_keywords, collapse = "|"), names(sets), ignore.case = TRUE)]
dbof <- ifelse(grepl("^GOBP_", names(sets)), "GO:BP", "Reactome"); names(dbof) <- names(sets)

out <- list()
for (rg in c("Frontal", "Hippo")) {
  f <- file.path(PROJ, "tables", "mast_dual", sprintf("cliffs_delta_full_Micro-PVM_%s.csv", rg))
  tbl <- fread(f)[!is.na(avg_log2FC)]
  n0 <- nrow(tbl)
  tbl <- tbl[pmax(pct.1, pct.2) >= 0.10]                                        # script-44 floor
  tbl[, artifact := is_artifact(gene, "Micro-PVM")]
  tbl[, gated_out := artifact | pmin(pct.1, pct.2) < MIN_PCT_BOTH]
  say("%s: %d genes -> %d after the 44 floor -> %d after artifact purge + min(pct) >= %.2f (removed %d artifact, %d low-detection)",
      rg, n0, nrow(tbl), sum(!tbl$gated_out), MIN_PCT_BOTH, sum(tbl$artifact), sum(!tbl$artifact & tbl$gated_out))
  g <- tbl[gated_out == FALSE]
  ranks <- sort(setNames(fc_rank(g$avg_log2FC), g$gene), decreasing = TRUE)
  set.seed(42)
  fg <- fgsea(pathways = sets, stats = ranks, minSize = 10, maxSize = 500, gseaParam = 1, nPermSimple = 1000)
  fg <- as.data.table(fg)[!is.na(padj)]
  fg[, `:=`(subtype = paste0("Micro-PVM_", rg), db = dbof[pathway], direction = ifelse(NES > 0, "NHD", "CON"),
            leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"), metric = "log2fc_winsorized_gated",
            n_ranked = length(ranks), gate = sprintf("not is_artifact & min(pct) >= %.2f", MIN_PCT_BOTH))]
  say("  %s: %d sets tested; FDR < 0.05 & size >= 15: NHD-up %d, lost %d", rg, nrow(fg),
      fg[padj < 0.05 & size >= 15 & NES > 0, .N], fg[padj < 0.05 & size >= 15 & NES < 0, .N])
  out[[rg]] <- fg
}
out <- rbindlist(out)
fwrite(out[, .(subtype, pathway, pval, padj, ES, NES, size, leadingEdge, direction, db, metric, n_ranked, gate)],
       file.path(PROJ, "tables", "micro_gsea_gated_FH.csv"))
say("\n== frontal, FDR < 0.05, size >= 15 ==")
print(out[subtype == "Micro-PVM_Frontal" & padj < 0.05 & size >= 15][order(-NES), .(pathway = substr(pathway, 1, 60), NES = round(NES, 2), padj = signif(padj, 2), size)], nrows = 60)
writeLines(capture.output(sessionInfo()), file.path(PROJ, "logs", "110b_micro_gsea_gated_FH_sessionInfo.txt"))
cat("\n=== DONE ===\n", file = stderr())
