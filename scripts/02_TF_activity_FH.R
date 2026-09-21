# =============================================================================
# 02_TF_activity_FH.R — Frontal+Hippocampus-only rebuild (Occipital dropped)
# Faithful port of manuscript_7fig/scripts/02_TF_activity.R to the new FH DGE
# inputs. not a redesign — same tool (decoupleR::run_ulm), same GTRD regulon,
# same output schema (statistic, source, condition, score, p_value, comparison),
# same region-stripping regex `_(Frontal|OCC|Hippo)$` (OCC simply never appears).
#
# input-matrix CHOICE (the one substantive port decision) -- documented here so a
# reviewer can see exactly what changed and why:
#   The original fed the pseudobulk DESeq2 Wald `stat` (per gene x comparison) to
#   run_ulm (confirmed: old ST4 has condition == "stat"). In the FH rebuild the
#   pseudobulk table has the Wald `stat` NULLED for every Hippo comparison
#   (dispersion_basis == "pseudo_split", stat/pvalue/padj retired). Feeding that
#   table would yield zero output for all 6 Hippo comparisons — half the TF layer.
#   Per the manuscript's MAST-forward policy we therefore use the MAST
#   avg_log2FC as the ranking statistic, taken from the full-expressed caches
#   (cliffs_delta_full_<ct>_<region>.csv), which give avg_log2FC over the full
#   detectable transcriptome per comparison — the direct analog of a per-gene
#   Wald stat vector, and the one choice that keeps all 12 comparisons computable.
#   `condition` in the output is set to "avg_log2FC" so the
#   provenance of the ranking metric is self-documenting in the CSV.
#
# Outputs (drop-in for the FH consumers):
#   tables/TF_main_celltype_activity_per_comparison.csv
#   supplementary_tables/ST4_TF_activity_decoupleR_ULM_GTRD.csv (== main, canonical copy)
#   supplementary_tables/ST5_TF_regulon_GTRD_unsigned.csv     (regulon, renamed cols)
#   data/TF_regulon_GTRD.rds                                 (copied from manuscript_7fig)
#   logs/02_TF_activity_FH.log                               (this run's provenance)
#
# The second original output (cortical-neuron-subclass TF activity) is
# deferred — FH subclass-level DGE has not been regenerated yet. It is not
# fabricated here; only the main-celltype TF layer is produced.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(decoupleR)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

# --- project root (NHD_PROJ) ------------------------------
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
stopifnot("MISSING project root" = !is.na(PROJ) && dir.exists(PROJ))

RB   <- file.path(PROJ, "NHD_frontal_hippo_rebuild")
DATA <- file.path(RB, "data")
TBL  <- file.path(RB, "tables")
STBL <- file.path(RB, "supplementary_tables")
MAST <- file.path(TBL, "mast_dual")
LOGD <- file.path(RB, "logs")
for (d in c(DATA, TBL, STBL, LOGD)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- tee a provenance log next to the outputs -------------------------------
LOGF <- file.path(LOGD, "02_TF_activity_FH.log")
con  <- file(LOGF, open = "wt"); sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)
cat("== 02_TF_activity_FH.R ==\n", format(Sys.time()), "\n")
cat("PROJ:", PROJ, "\n\n")

# ============================================================================
# 1. TF regulon network (GTRD) — data-independent; copy the frozen artifact
#    from manuscript_7fig rather than regenerating (avoids msigdbr version
#    drift; the original was built once and saved as TF_regulon_GTRD.rds).
# ============================================================================
cat("== step: TF regulon (GTRD) ==\n")
REG_SRC <- file.path(PROJ, "manuscript_7fig", "data", "TF_regulon_GTRD.rds")
REG_DST <- file.path(DATA, "TF_regulon_GTRD.rds")
stopifnot("MISSING source regulon manuscript_7fig/data/TF_regulon_GTRD.rds" =
            file.exists(REG_SRC))
file.copy(REG_SRC, REG_DST, overwrite = TRUE)
net <- readRDS(REG_DST)
stopifnot(all(c("source", "target", "mor") %in% colnames(net)), nrow(net) > 0)
cat(sprintf("Copied GTRD regulon: %d interactions, %d TFs (source: manuscript_7fig)\n",
            nrow(net), length(unique(net$source))))

# ============================================================================
# 2. Assemble the per-comparison ranking matrix from the MAST full-expressed
#    caches (avg_log2FC over the full detectable transcriptome per comparison).
# ============================================================================
cat("\n== step: load MAST full-expressed ranking caches ==\n")
cell_types <- c("Micro-PVM", "Astro", "Oligo", "OPC", "Neuron_Ex", "Neuron_Inh")
regions    <- c("Frontal", "Hippo")           # no OCC in the FH rebuild
comparisons <- as.vector(t(outer(cell_types, regions, paste, sep = "_")))

STAT_COL <- "avg_log2FC"                       # ranking metric (see header)
cache_path <- function(cmp) file.path(MAST, paste0("cliffs_delta_full_", cmp, ".csv"))

# fail loud & early if any expected FH cache is missing
missing <- comparisons[!file.exists(vapply(comparisons, cache_path, character(1)))]
if (length(missing))
  stop("MISSING full-expressed cache(s): ", paste(missing, collapse = ", "),
       " — the MAST-dual prep step must run first.")

dge_list <- lapply(comparisons, function(cmp) {
  d <- read.csv(cache_path(cmp), stringsAsFactors = FALSE)
  stopifnot(all(c("gene", STAT_COL) %in% colnames(d)))
  # strip region exactly as the original does (regex kept verbatim; OCC unused)
  d$comparison <- cmp
  d$cell_type  <- sub("_(Frontal|OCC|Hippo)$", "", cmp)
  d$region     <- sub("^.*_", "", cmp)
  d[, c("gene", STAT_COL, "comparison", "cell_type", "region")]
})
dge <- bind_rows(dge_list)
cat(sprintf("Loaded %d comparisons, %d total gene rows.\n",
            length(comparisons), nrow(dge)))

# ============================================================================
# 3. Run decoupleR ULM per comparison (identical call to the original;
#    minsize = 5, mor from the network).
# ============================================================================
cat("\n== step: decoupleR run_ulm per comparison ==\n")
results_list <- list()
for (cmp in comparisons) {
  sub <- dge[dge$comparison == cmp, ]
  # gene must be non-missing & unique (original did the same before matrix build)
  sub <- sub[!is.na(sub$gene) & !duplicated(sub$gene) & !is.na(sub[[STAT_COL]]), ]
  if (nrow(sub) < 5) { cat(sprintf("  %-20s SKIP (<5 genes)\n", cmp)); next }
  mat <- matrix(sub[[STAT_COL]], nrow = nrow(sub), ncol = 1,
                dimnames = list(sub$gene, STAT_COL))
  res <- tryCatch(
    decoupleR::run_ulm(mat = mat, network = net,
                       .source = "source", .target = "target",
                       .mor = "mor", minsize = 5),
    error = function(e) { cat("  ERROR for", cmp, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  res$comparison <- cmp
  results_list[[cmp]] <- res
  n_sig <- sum(res$p_value < 0.05, na.rm = TRUE)
  cat(sprintf("  %-20s genes=%-6d TFs tested=%-4d  p<0.05=%d\n",
              cmp, nrow(sub), nrow(res), n_sig))
}
out_df <- bind_rows(results_list)
stopifnot("run_ulm produced no output" = nrow(out_df) > 0)

# Guard: the OCC region must appear nowhere in the output
stopifnot("OCC leaked into TF output" =
            !any(grepl("_OCC$", out_df$comparison)) &&
            !any(grepl("OCC", out_df$comparison)))

# ============================================================================
# 4. Write outputs (drop-in schema: statistic, source, condition, score,
#    p_value, comparison — run_ulm's `condition` field carries the matrix
#    column name = "avg_log2FC" here, was "stat" in the original).
# ============================================================================
cat("\n== step: write outputs ==\n")
main_path <- file.path(TBL, "TF_main_celltype_activity_per_comparison.csv")
write.csv(out_df, main_path, row.names = FALSE)
cat("Wrote:", main_path, "\n")

# ST4 == canonical copy of the main table (same schema, supplementary home)
# ST4/ST5 are the GTRD genome-wide
# survey across all cell types; they do not back Figure 2h, which is microglia-only on
# the SIGNED CollecTRI network (scripts 68 + 70b, shipped as ST14).
st4_path <- file.path(STBL, "ST4_TF_activity_decoupleR_ULM_GTRD.csv")
write.csv(out_df, st4_path, row.names = FALSE)
cat("Wrote:", st4_path, "\n")

# ST5 == regulon with the supplementary column names (source,target,MoR,weight)
st5 <- net %>%
  transmute(source, target, mode_of_regulation = mor,
            weight = if ("likelihood" %in% colnames(net)) likelihood else 1)
st5_path <- file.path(STBL, "ST5_TF_regulon_GTRD_unsigned.csv")
write.csv(st5, st5_path, row.names = FALSE)
cat("Wrote:", st5_path, "\n")

# ============================================================================
# 5. Verification summary in the log
# ============================================================================
cat("\n== VERIFY ==\n")
cat("Unique TFs scored:", length(unique(out_df$source)), "\n")
cat("Comparisons present (expect 12, 0 OCC):\n"); print(sort(unique(out_df$comparison)))
# per-comparison significance (same column the consumers ultimately threshold:
# p_value here; the heatmap/network scripts apply BH within comparison downstream)
sig_tab <- out_df %>%
  group_by(comparison) %>%
  summarise(n_TF = n(),
            n_p05 = sum(p_value < 0.05, na.rm = TRUE),
            n_q10_BH = sum(p.adjust(p_value, "BH") < 0.10, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(comparison)
cat("\nPer-comparison TF significance:\n"); print(as.data.frame(sig_tab), row.names = FALSE)

# --- neuron-TF new-data check (report direction; do not assume old list holds) --
cat("\n== NEURON-TF NEW-DATA CHECK (score sign: - = NHD-down) ==\n")
neuron_tfs_watch <- c("CUX1","EMX1","MEF2C","MEF2D","PPARGC1A",
                      "TFEB","NFE2L1","ATF5")
chk <- out_df %>%
  filter(source %in% neuron_tfs_watch,
         grepl("^Neuron_(Ex|Inh)_", comparison)) %>%
  select(source, comparison, score, p_value) %>%
  mutate(score = round(score, 3), p_value = signif(p_value, 3)) %>%
  arrange(source, comparison)
if (nrow(chk)) print(as.data.frame(chk), row.names = FALSE) else
  cat("(none of the watched neuron TFs carry a regulon in this collection)\n")
absent <- setdiff(neuron_tfs_watch, out_df$source)
if (length(absent))
  cat("\nWatched TFs with NO regulon in GTRD (never scorable):",
      paste(absent, collapse = ", "), "\n")

cat("\n=== DONE ===\n")
