#!/usr/bin/env Rscript
# =============================================================================
# 44_gsea_MAST_cliffs_FH.R — GSEA preranked on signed winsorized avg_log2FC, frontal+hippo rebuild, broad cell types. Faithful copy of
# manuscript_7fig/scripts/44_gsea_MAST_cliffs.R with rebuild paths + no OCC.
# Reads the full-expressed caches from 47_full_expressed_ranking_FH.R (broad only).
# Subclass + hippo-subtype GSEA are deferred to the neuron-figure stage.
# Metric = signed avg_log2FC winsorized 1/99 pct (p-value-free; pseudorep-safe).
# BH within each comparison family. Power tiers: excluded<50 / low_conf 50-149 / solid>=150.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr); library(fgsea); library(msigdbr)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
DIR      <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
MAST_DIR <- file.path(DIR, "tables/mast_dual")
OUT_DIR  <- file.path(DIR, "diagnostics/05_GSEA/fGSEA_NHD_MAST_cliffs"); dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(DIR, "scripts/_gsea_composite_metric.R"))
msg <- function(...) { cat(sprintf(...), "\n"); flush.console() }

BROAD_TYPES <- c("Micro-PVM","Astro","Oligo","OPC","Neuron_Ex","Neuron_Inh")
PCT_FLOOR <- 0.10; HIDE_N <- 50; FADE_N <- 150
POWER_TIERS <- c("excluded","low_conf","solid"); FDR_THRESH <- 0.05

## STEP 1: MAST discovery table (broad)
# microglia-only primary: the myeloid compartment is the split class, so
# the GSEA prerank is built on Microglia (and PVM) rather than the pooled "Micro-PVM".
# Comparisons are derived from this table, so pointing it here is sufficient; the
# matching cliffs_delta_full_{Microglia,PVM}_* caches come from 47 (same change).
broad_path <- file.path(MAST_DIR, "MAST_dual_discovery_all.csv")
stopifnot("MISSING MAST_dual_discovery_all.csv — run 40b + 41b" = file.exists(broad_path))
mast_all <- read_csv(broad_path, show_col_types = FALSE) %>% mutate(family = "broad")
cmp_meta <- mast_all %>% group_by(comparison) %>%
  summarise(cell_type = first(cell_type), region = first(region), family = first(family),
            n_NHD = first(n_NHD), n_CON = first(n_CON), .groups = "drop") %>%
  mutate(n_nuclei_min = pmin(n_NHD, n_CON),
         power_tier = cut(n_nuclei_min, breaks = c(-Inf, HIDE_N, FADE_N, Inf),
                          labels = POWER_TIERS, right = FALSE, ordered_result = TRUE),
         low_n = n_nuclei_min < FADE_N)
comparisons <- cmp_meta$comparison
DISC_GENES <- mast_all %>% filter(discovery %in% c(TRUE, "TRUE")) %>%
  group_by(comparison) %>% summarise(g = list(unique(gene)), .groups = "drop")
disc_genes_for <- function(cmp) { i <- match(cmp, DISC_GENES$comparison); if (is.na(i)) character(0) else DISC_GENES$g[[i]] }

## STEP 2: full-expressed caches (47_FH)
cache_path <- function(cmp) file.path(MAST_DIR, sprintf("cliffs_delta_full_%s.csv", cmp))
missing <- comparisons[!file.exists(vapply(comparisons, cache_path, character(1)))]
if (length(missing)) stop(sprintf("MISSING %d cache(s) — run 47_full_expressed_ranking_FH.R: %s",
                                   length(missing), paste(missing, collapse=", ")))
full_delta <- bind_rows(lapply(comparisons, function(cmp)
  read_csv(cache_path(cmp), show_col_types = FALSE) %>% mutate(comparison = cmp)))
msg("loaded %d cache rows across %d comparisons", nrow(full_delta), length(comparisons))

## STEP 3: gene sets
go_bp    <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP") %>% {split(.$gene_symbol, .$gs_name)}
reactome <- msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME") %>% {split(.$gene_symbol, .$gs_name)}
broad_sets <- c(go_bp, reactome)
blacklist_keywords <- c("KERATINOCYTE","EPIDERM","HAIR_FOLLICLE","HAIR_CYCLE","TOOTH","ODONTOGENESIS",
  "AMELOBLAST","LENS","CORNIF","CARDIAC_MUSCLE","CARDIAC_CHAMBER","HEART_MORPHO","B_CELL_ACTIV",
  "THYMO","SPERMAT","OOCYTE","OVULATION","RETINA_LAYER","PHOTORECEPTOR_CELL_DIFF","ALCOHOL","ETHANOL")
bl <- paste(blacklist_keywords, collapse="|")
broad_sets <- broad_sets[!grepl(bl, names(broad_sets), ignore.case = TRUE)]
msg("GO:BP=%d Reactome=%d -> broad %d", length(go_bp), length(reactome), length(broad_sets))
clean_pathway <- function(x, wrap=55) { out <- sub("^(GOBP_|REACTOME_|HALLMARK_)","",x); out <- gsub("_"," ",out); str_wrap(str_to_sentence(out), width=wrap) }

## STEP 4: fgsea per comparison (BH within comparison)
run_one_gsea <- function(cmp) {
  info <- cmp_meta[cmp_meta$comparison == cmp, ]
  tbl  <- full_delta %>% filter(comparison == cmp) %>%
    filter(pmax(pct.1, pct.2) >= PCT_FLOOR, !is.na(avg_log2FC)) %>% distinct(gene, .keep_all = TRUE)
  if (nrow(tbl) < 50) { msg("  %-22s SKIP (%d ranked genes)", cmp, nrow(tbl)); return(NULL) }
  comp  <- fc_rank(tbl$avg_log2FC)
  ranks <- sort(setNames(comp, tbl$gene), decreasing = TRUE); ranks <- ranks[!is.na(ranks)]
  sets <- broad_sets
  dge_g <- disc_genes_for(cmp)
  set.seed(42)
  fg <- tryCatch(fgsea(pathways=sets, stats=ranks, minSize=10, maxSize=500, gseaParam=1, nPermSimple=1000),
                 error=function(e){ msg("  fgsea err %s: %s", cmp, conditionMessage(e)); NULL })
  if (is.null(fg) || !nrow(fg)) return(NULL)
  has_dge <- if (length(dge_g)) vapply(sets[fg$pathway], function(g) any(g %in% dge_g), logical(1)) else rep(FALSE, nrow(fg))
  fg <- as_tibble(fg) %>% mutate(
      subtype=cmp, direction=if_else(NES>0,"NHD","CON"),
      db=if_else(grepl("^REACTOME_",pathway),"Reactome","GO:BP"),
      pathway_clean=clean_pathway(pathway), metric=GSEA_METRIC_NAME,
      n_nuclei_min=info$n_nuclei_min, power_tier=as.character(info$power_tier), low_n=info$low_n,
      dge_member=has_dge, leadingEdge=vapply(leadingEdge, paste, collapse="|", FUN.VALUE=character(1))) %>%
    arrange(padj, desc(abs(NES)))
  msg("  %-22s %5d tested | sig(FDR<%.2g)=%d | DGE-member=%d | ranked=%d",
      cmp, nrow(fg), FDR_THRESH, sum(fg$padj<FDR_THRESH,na.rm=TRUE), sum(fg$dge_member,na.rm=TRUE), length(ranks))
  fg %>% select(subtype, pathway, pval, padj, ES, NES, size, leadingEdge, direction, db,
                pathway_clean, metric, n_nuclei_min, power_tier, low_n, dge_member)
}
gsea_all <- bind_rows(lapply(comparisons, run_one_gsea))
write_csv(gsea_all, file.path(OUT_DIR, "fGSEA_NHD_MAST_cliffs_all.csv"))
msg("Saved fGSEA_NHD_MAST_cliffs_all.csv (%d rows)", nrow(gsea_all))

## STEP 6: summary
s2 <- gsea_all %>% group_by(subtype) %>%
  summarise(pathways=dplyr::n(), sig_FDR005=sum(padj<FDR_THRESH,na.rm=TRUE),
            sig_up_NHD=sum(padj<FDR_THRESH & NES>0,na.rm=TRUE),
            sig_up_CON=sum(padj<FDR_THRESH & NES<0,na.rm=TRUE),
            n_nuclei_min=first(n_nuclei_min), power_tier=first(power_tier), .groups="drop")
print(as.data.frame(s2), row.names = FALSE)
write_csv(s2, file.path(OUT_DIR, "fGSEA_NHD_MAST_cliffs_summary.csv"))
cat("\n=== DONE (broad; concordance-vs-pseudobulk-GSEA deferred) ===\n", file = stderr())
