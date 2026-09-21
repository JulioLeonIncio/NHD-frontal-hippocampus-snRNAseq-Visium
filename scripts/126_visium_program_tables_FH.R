#!/usr/bin/env Rscript
# =============================================================================
# 126_visium_program_tables_FH.R — Supplementary Data 11, sheets Programs_per_domain / Genes_by_band / Gene_sets_Visium (internal ST23): the Visium program expression
# per spatial domain, three sheets, copied from the tables behind Fig. 6d,e with referee-facing column names.
# -----------------------------------------------------------------------------
# What. Nothing is re-derived. The three sheets are the panel tables already on disk, renamed and asserted:
#   sheet Data           tables/fig6e_program_umi_per10k_by_domain_FH.csv   (67b; Fig. 6e per-domain means/medians)
#                        per spot: UMIs of the program's genes per 10,000 UMIs at native depth (no 3,000-UMI floor),
#                        as a ratio of sums over the spot and its present hexagonal array neighbours (six offsets,
#                        within section); then mean / median / n over the spots of each program x condition x domain,
#                        the two sections of a condition pooled. Domain = the final spatial domain (as ST22).
#   sheet Genes_by_band  tables/figF6d_supplychain_genes_FH.csv             (73; Fig. 6d per-gene means + log2 ratio)
#                        per spot: gene UMIs per 10,000 UMIs at native depth, no neighbour pooling; mean per band x
#                        condition over the six bands drawn in Fig. 6d (DOM_BANDS = grey layers + WM; InN and
#                        Vasc/immune are not bands); lfc = log2((NHD + 1e-4) / (CON + 1e-4)) of those means, genes with
#                        CON mean = 0 dropped upstream (73). Gene membership = the Fig. 6d blocks (73's SETS), which for
#                        Complement / MHC-II and Reactive astrocyte are wider than the scored modules of sheet Gene_sets.
#   sheet Gene_sets      tables/visium_dm_module_genes_FH.csv                (65; module membership as scored)
#                        present = the gene was on the slide and therefore inside the AddModuleScore score.
#
# What it is not. Not a re-run of 65 / 67b / 73; if a source table is stale the guards below fail (6e spot counts
# per domain x condition must equal the shipped ST22 annotation spot by spot count).
#
# Inputs:  tables/fig6e_program_umi_per10k_by_domain_FH.csv        (67b)
#          tables/figF6d_supplychain_genes_FH.csv                  (73)
#          tables/visium_dm_module_genes_FH.csv                    (65)
#          supplementary_tables/ST22_visium_domain_annotation_FH.csv (125; the spot-count reference)
#          ../Visium/integrated_harmony/integrated_domain_levels.csv (via _visium_domains_FH.R: DOM_LEV / DOM_BANDS)
# Outputs: supplementary_tables/ST23_visium_program_by_domain_FH.csv        (sheet Programs_per_domain of Supplementary Data 11)
#          supplementary_tables/ST23b_visium_program_genes_by_band_FH.csv   (sheet Genes_by_band; script 100 ships
#          supplementary_tables/ST23c_visium_program_gene_sets_FH.csv        sheet Gene_sets;  b/c never as own files)
#          logs/126_visium_program_tables_FH.log                             (provenance: counts + sessionInfo)
#
# Column DICTIONARY (the referee-facing wording lives in scripts/100_* as ST23_COLS / ST23B_COLS / ST23C_COLS):
#   Data:           program, condition, domain, domain_order, mean_umi_per_10k, median_umi_per_10k, n_spots
#   Genes_by_band:  program, gene, band, band_order, mean_umi_per_10k_CON, mean_umi_per_10k_NHD, log2_ratio_NHD_vs_CON
#   Gene_sets:      program, gene, detected_on_slide
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "_visium_domains_FH.R"))   # DOM_LEV / DOM_INFO / DOM_BANDS

TAB  <- file.path(PROJ, "tables")
OUTD <- file.path(PROJ, "supplementary_tables"); LOGD <- file.path(PROJ, "logs")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE); dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
F_E   <- file.path(TAB, "fig6e_program_umi_per10k_by_domain_FH.csv")
F_D   <- file.path(TAB, "figF6d_supplychain_genes_FH.csv")
F_G   <- file.path(TAB, "visium_dm_module_genes_FH.csv")
F_22  <- file.path(OUTD, "ST22_visium_domain_annotation_FH.csv")
OUT_A <- file.path(OUTD, "ST23_visium_program_by_domain_FH.csv")
OUT_B <- file.path(OUTD, "ST23b_visium_program_genes_by_band_FH.csv")
OUT_C <- file.path(OUTD, "ST23c_visium_program_gene_sets_FH.csv")
LOG   <- file.path(LOGD, "126_visium_program_tables_FH.log")
# the numbers the Fig. 6e table must reproduce: every spot of the final object (18,765 = 9,359 CON + 9,406 NHD), two programs
N_CON <- 9359L; N_NHD <- 9406L
PROGS_E <- c("Complement / MHC-II", "Fatty-acid / sphingomyelin")
BANDS_6D <- c("L1/pia", "L2/3", "L4", "L5", "L6", "WM")

for (f in c(F_E, F_D, F_G, F_22))
  if (!file.exists(f)) stop("MISSING ", f, " — upstream step (65 / 67b / 73 / 125) has not run")
stopifnot("Fig. 6d bands (DOM_BANDS from integrated_domain_levels.csv) are not the expected six" = identical(DOM_BANDS, BANDS_6D),
          all(BANDS_6D %in% DOM_LEV))

# ---- sheet Data: Fig. 6e per-domain program signal (67b) -------------------
cat("== sheet Data (67b) ==\n")
e <- fread(F_E)
stopifnot("fig6e table columns differ from 67b" = identical(names(e), c("arm", "condition", "spatial_domain", "mean_per10k", "median_per10k", "n")))
setnames(e, c("arm", "condition", "spatial_domain", "mean_per10k", "median_per10k", "n"),
            c("program", "condition", "domain", "mean_umi_per_10k", "median_umi_per_10k", "n_spots"))
e[, domain_order := match(domain, DOM_LEV)]
stopifnot("a 6e domain is not in integrated_domain_levels.csv" = !anyNA(e$domain_order),
          "6e must cover every domain of the vocabulary" = setequal(unique(e$domain), DOM_LEV),
          setequal(unique(e$program), PROGS_E), setequal(unique(e$condition), c("CON", "NHD")),
          "expected 2 programs x 2 conditions x 8 domains = 32 rows" = nrow(e) == 2L * 2L * length(DOM_LEV),
          !anyDuplicated(e[, .(program, condition, domain)]),
          !anyNA(e$mean_umi_per_10k), !anyNA(e$median_umi_per_10k), all(e$mean_umi_per_10k >= 0), all(e$median_umi_per_10k >= 0),
          is.integer(e$n_spots) || all(e$n_spots == round(e$n_spots)), all(e$n_spots > 0))
e[, n_spots := as.integer(n_spots)]
.tot <- e[, .(n = sum(n_spots)), by = .(program, condition)]
stopifnot("6e spot totals differ from the final object (9,359 CON / 9,406 NHD per program) — re-run 67b" =
            all(.tot[condition == "CON"]$n == N_CON) && all(.tot[condition == "NHD"]$n == N_NHD))
# consistency with the shipped per-spot annotation: n per domain x condition must equal the ST22 counts (same object, all spots)
n22 <- fread(F_22)[, .(n22 = .N), by = .(domain, condition)]
.m <- merge(e[program == PROGS_E[1], .(domain, condition, n_spots)], n22, by = c("domain", "condition"), all = TRUE)
stopifnot("6e spot counts per domain x condition disagree with ST22 (stale 67b or 125 run?)" = !anyNA(.m$n_spots), !anyNA(.m$n22), all(.m$n_spots == .m$n22))
setcolorder(e, c("program", "condition", "domain", "domain_order", "mean_umi_per_10k", "median_umi_per_10k", "n_spots"))
setorder(e, program, condition, domain_order)
stopifnot(ncol(e) == 7)

# ---- sheet Genes_by_band: Fig. 6d per-gene means + log2 ratio (73) ---------
cat("== sheet Genes_by_band (73) ==\n")
g <- fread(F_D)
stopifnot("figF6d table columns differ from 73" = identical(names(g), c("programme", "gene", "band", "CON", "NHD", "lfc")))
setnames(g, c("programme", "gene", "band", "CON", "NHD", "lfc"),
            c("program", "gene", "band", "mean_umi_per_10k_CON", "mean_umi_per_10k_NHD", "log2_ratio_NHD_vs_CON"))
g[, band_order := match(band, DOM_LEV)]
stopifnot("6d bands must be exactly the six drawn (L1/pia, L2/3, L4, L5, L6, WM)" = setequal(unique(g$band), BANDS_6D), !anyNA(g$band_order),
          "NA in the 6d table" = !anyNA(g),
          !anyDuplicated(g[, .(program, gene, band)]),
          all(g$mean_umi_per_10k_CON > 0), all(g$mean_umi_per_10k_NHD >= 0), all(is.finite(g$log2_ratio_NHD_vs_CON)))
.per <- g[, .N, by = .(program, gene)]
stopifnot("every gene must appear in all six bands" = all(.per$N == length(BANDS_6D)),
          "a gene appears under two programs" = !anyDuplicated(.per$gene))
# the shipped ratio must be the stated formula of 73: log2((NHD + 1e-4) / (CON + 1e-4)) of the per-band means
stopifnot("log2 ratio is not log2((NHD+1e-4)/(CON+1e-4)) of the shipped means — 73 changed its formula" =
            isTRUE(all.equal(g$log2_ratio_NHD_vs_CON, log2((g$mean_umi_per_10k_NHD + 1e-4) / (g$mean_umi_per_10k_CON + 1e-4)), tolerance = 1e-8)))
setcolorder(g, c("program", "gene", "band", "band_order", "mean_umi_per_10k_CON", "mean_umi_per_10k_NHD", "log2_ratio_NHD_vs_CON"))
# row order = the panel's block order (73 writes programme block -> gene -> band); kept, not re-sorted
stopifnot(ncol(g) == 7)

# ---- sheet Gene_sets: module membership as scored (65) ---------------------
cat("== sheet Gene_sets (65) ==\n")
s <- fread(F_G)
stopifnot("visium_dm_module_genes columns differ from 65" = identical(names(s), c("module", "gene", "present")))
setnames(s, c("module", "gene", "present"), c("program", "gene", "detected_on_slide"))
stopifnot(is.logical(s$detected_on_slide), !anyNA(s$detected_on_slide), !anyNA(s$program), !anyNA(s$gene), all(nzchar(s$gene)),
          "expected >= 10 scored programs" = uniqueN(s$program) >= 10,
          "duplicate program x gene" = !anyDuplicated(s[, .(program, gene)]),
          "the two Fig. 6e programs must be scored programs" = all(PROGS_E %in% s$program))
stopifnot(ncol(s) == 3)

# ---- write (temp then rename; the shipped folder is never left half-written) ----
.wr <- function(x, p) { tmp <- paste0(p, ".tmp"); fwrite(x, tmp); invisible(file.rename(tmp, p)) }
.wr(e, OUT_A); .wr(g, OUT_B); .wr(s, OUT_C)
# copy check: the shipped values must equal the source values (only names / order / the *_order column were added)
.e0 <- fread(F_E); .g0 <- fread(F_D); .s0 <- fread(F_G)
.eb <- fread(OUT_A); .gb <- fread(OUT_B); .sb <- fread(OUT_C)
stopifnot(isTRUE(all.equal(.eb[order(program, condition, domain), .(program, condition, domain, mean_umi_per_10k, median_umi_per_10k, n_spots)],
                           .e0[order(arm, condition, spatial_domain), .(program = arm, condition, domain = spatial_domain, mean_umi_per_10k = mean_per10k, median_umi_per_10k = median_per10k, n_spots = as.integer(n))],
                           check.attributes = FALSE)),
          isTRUE(all.equal(.gb[, .(program, gene, band, mean_umi_per_10k_CON, mean_umi_per_10k_NHD, log2_ratio_NHD_vs_CON)],
                           .g0[, .(program = programme, gene, band, mean_umi_per_10k_CON = CON, mean_umi_per_10k_NHD = NHD, log2_ratio_NHD_vs_CON = lfc)],
                           check.attributes = FALSE)),
          isTRUE(all.equal(.sb, .s0[, .(program = module, gene, detected_on_slide = present)], check.attributes = FALSE)))

# ---- summary + provenance -------------------------------------------------
summ <- c(
  sprintf("ST23  Data:          %d rows x %d columns -> %s", nrow(e), ncol(e), OUT_A),
  sprintf("ST23b Genes_by_band: %d rows x %d columns (%d genes, %d programs, %d bands) -> %s", nrow(g), ncol(g), uniqueN(g$gene), uniqueN(g$program), uniqueN(g$band), OUT_B),
  sprintf("ST23c Gene_sets:     %d rows x %d columns (%d programs, %d genes, %d detected on slide) -> %s", nrow(s), ncol(s), uniqueN(s$program), nrow(s), sum(s$detected_on_slide), OUT_C),
  "Data: n_spots per program x condition:", capture.output(print(.tot)),
  "Data: mean UMI per 10,000 by domain (rows) x program / condition:",
  capture.output(print(dcast(e, domain_order + domain ~ program + condition, value.var = "mean_umi_per_10k"))),
  "Genes_by_band: genes per program:", capture.output(print(g[, .(genes = uniqueN(gene)), by = program])),
  sprintf("Genes_by_band: log2_ratio_NHD_vs_CON range %.2f to %.2f", min(g$log2_ratio_NHD_vs_CON), max(g$log2_ratio_NHD_vs_CON)),
  "Gene_sets: members per program:", capture.output(print(s[, .(members = .N, detected = sum(detected_on_slide)), by = program])))
cat(summ, sep = "\n")
writeLines(c(sprintf("126_visium_program_tables_FH.R  run %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             sprintf("inputs: %s", paste(sprintf("%s (mtime %s)", basename(c(F_E, F_D, F_G, F_22)), format(file.mtime(c(F_E, F_D, F_G, F_22)))), collapse = "; ")),
             "", summ, "", capture.output(sessionInfo())), LOG)
cat("log:", LOG, "\n=== DONE ===\n")
