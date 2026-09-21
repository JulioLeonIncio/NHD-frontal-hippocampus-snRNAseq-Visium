#!/usr/bin/env Rscript
# =============================================================================
# 127_program_delta_table_FH.R — Supplementary Data 7, sheets Program_deltas + By_program (internal ST24): program-level Cliff's delta, NHD versus
# control, per cell type and region, two sheets, copied from the four panel statistics tables with referee-facing
# column names.
# -----------------------------------------------------------------------------
# What. Nothing is re-derived. The four tables behind the program / signature-score panels are stacked, renamed and
# asserted (the statistical unit is the nucleus; one donor per condition, so delta is the effect size and q ranks only):
#   Fig. 2d  microglia (Micro-PVM)        tables/fig2_micro_program_violin_stats_FH.csv   (2B;  7 programs x 2 regions = 14)
#   Fig. 3e  astrocytes                   tables/fig3d_astro_signature_stats.csv          (3D; 11 signatures x 2 regions = 22;
#                                         the file keeps its historical "3d" name — the violins are panel e of the assembled
#                                         Fig. 3, the dot-matrix is panel f; IMPROVEMENTS_20260907.md)
#   Fig. 4d  oligodendrocyte lineage      tables/oligo_lineage_deltas_FH.csv              (3Fp; 8 programs x 2 lineages x 2 regions = 32)
#   Fig. 5d  neurons (Ex / Inh)           tables/fig5_neuron_programme_dotmatrix_FH.csv   (4C2; 5 programs x 2 classes x 2 regions = 20)
#   not tables/4C_signature_violin_stats_FH.csv: that is the superseded excitatory-only neuron violin (4C, Jul 13) whose
#   10 rows duplicate the Fig. 5d Ex rows; it backs no assembled panel.
#   sheet Data        one row per figure_panel x cell_type x program x region (88 rows)
#   sheet By_program  one row per figure_panel x cell_type x program, both regions side by side (44 rows)
#
# What is RE-DERIVED (checks only, never shipped in place of the source):
#   * q_BH is recomputed as p.adjust(p, "BH") within each panel and must equal the source column (asserts the
#     "BH within panel" statement of the README).
#   * effect_gate = |delta| >= 0.15 (house strict gate) and q_star = q < 0.05 & |delta| < 0.15 (house grey mark) are
#     computed here and asserted against the source's own passes / passes_effect / label / star columns.
#   * p-value underflow: a source p or q written as exactly 0 is floored to the ST0 convention floor 2.2250739e-308
#     (read as p < 2.2e-308, never p = 0); the count per column is logged. BH is checked on the raw values first.
#
# Inputs:  the four tables above
# Outputs: supplementary_tables/ST24_program_cliffs_delta_FH.csv         (sheet Program_deltas of Supplementary Data 7)
#          supplementary_tables/ST24b_program_cliffs_delta_wide_FH.csv   (sheet By_program; script 100 ships it as sheet 2, never as its own file)
#          logs/127_program_delta_table_FH.log                            (provenance: counts + sessionInfo)
#
# Column DICTIONARY (the referee-facing wording lives in scripts/100_* as ST24_COLS / ST24B_COLS):
#   Data:       figure_panel, cell_type, program, region, n_CON, n_NHD, cliffs_delta, p_wilcox, q_BH, effect_gate, q_star
#   By_program: figure_panel, cell_type, program, cliffs_delta_Frontal, cliffs_delta_Hippocampus, q_BH_Frontal,
#               q_BH_Hippocampus, direction_consistent
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))

TAB  <- file.path(PROJ, "tables")
OUTD <- file.path(PROJ, "supplementary_tables"); LOGD <- file.path(PROJ, "logs")
dir.create(OUTD, showWarnings = FALSE, recursive = TRUE); dir.create(LOGD, showWarnings = FALSE, recursive = TRUE)
F_MI  <- file.path(TAB, "fig2_micro_program_violin_stats_FH.csv")
F_AS  <- file.path(TAB, "fig3d_astro_signature_stats.csv")
F_OL  <- file.path(TAB, "oligo_lineage_deltas_FH.csv")
F_NE  <- file.path(TAB, "fig5_neuron_programme_dotmatrix_FH.csv")
OUT_A <- file.path(OUTD, "ST24_program_cliffs_delta_FH.csv")
OUT_B <- file.path(OUTD, "ST24b_program_cliffs_delta_wide_FH.csv")
LOG   <- file.path(LOGD, "127_program_delta_table_FH.log")
for (f in c(F_MI, F_AS, F_OL, F_NE))
  if (!file.exists(f)) stop("MISSING ", f, " — upstream panel script (2B / 3D / 3Fp / 4C2) has not run")

# ---- constants (house rules; ST0 CONVENTIONS) --------------------------------
GATE    <- 0.15                     # |delta| >= 0.15 = the project's effect-size floor (strict gate)
ALPHA   <- 0.05                     # q threshold of the grey q* mark
P_FLOOR <- 2.2250739e-308           # ST0 "p-value underflow": .Machine$double.xmin rounded up to 8 s.f.
stopifnot(P_FLOOR > .Machine$double.xmin)
REGION_LEV <- c("Frontal", "Hippocampus")
REGION_MAP <- c(Frontal = "Frontal", Hippo = "Hippocampus", Hippocampus = "Hippocampus")   # normalise the source tokens
PANEL_LEV  <- c("Fig. 2d", "Fig. 3e", "Fig. 4d", "Fig. 5d")
CT_LEV     <- c("Micro-PVM", "Astro", "Oligo", "OPC", "Neuron_Ex", "Neuron_Inh")           # atlas cell-type vocabulary (ST1 / ST9)
# panel display labels of the Fig. 5d y axis (4C2 MOD_LAB / ORD); the CSV keys are the ST7 / ST6 "Fig5d_<key>" names
NEURON_LAB <- c(Presynaptic = "Presynaptic (SNARE)", Postsynaptic = "Postsynaptic (PSD)",
                GluR_ionotropic = "Ionotropic GluR", OxPhos = "OxPhos", IEG = "IEG")
DATA_COLS <- c("figure_panel", "cell_type", "program", "region", "n_CON", "n_NHD", "cliffs_delta", "p_wilcox", "q_BH", "effect_gate", "q_star")

# the panel mark columns hold "q*" or ""; an all-empty column (Fig. 2d) is typed logical NA by fread, so both are read as
# character and NA / "" = no mark (any other token is an error, never silently "no mark")
.is_qstar <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; stopifnot("unexpected token in a panel mark column" = all(x %in% c("", "q*"))); x == "q*" }
.norm_region <- function(x) { r <- unname(REGION_MAP[as.character(x)]); stopifnot("unknown region token in a source table" = !anyNA(r)); r }
# each panel block -> the common long schema, keeping the source's own gate / mark columns for the assertions below
.block <- function(panel, cell_type, program, region, n_CON, n_NHD, delta, p, q, src_gate, src_qstar, prog_order) {
  d <- data.table(figure_panel = panel, cell_type = cell_type, program = program, region = .norm_region(region),
                  n_CON = as.integer(n_CON), n_NHD = as.integer(n_NHD), cliffs_delta_raw = as.numeric(delta),
                  p_raw = as.numeric(p), q_raw = as.numeric(q), src_gate = src_gate, src_qstar = src_qstar)
  stopifnot(all(d$program %in% prog_order))
  d[, prog_ord := match(program, prog_order)]
  d
}

# ---- Fig. 2d microglia (2B) --------------------------------------------------
cat("== Fig. 2d microglia (2B) ==\n")
mi <- fread(F_MI, colClasses = list(character = "star"))
stopifnot("fig2 micro table columns differ from 2B" = identical(names(mi), c("program", "Region", "p_wilcox", "cliff_d", "y_star", "y_min", "y_max", "q_BH", "passes", "star", "delta_lab", "y_top", "y_delta")),
          nrow(mi) == 14L, is.logical(mi$passes))
b_mi <- .block("Fig. 2d", "Micro-PVM", mi$program, mi$Region, NA_integer_, NA_integer_, mi$cliff_d, mi$p_wilcox, mi$q_BH,
               mi$passes, .is_qstar(mi$star), unique(mi$program))     # 2B writes its rows in PROG_ORDER; n per group is not in the table

# ---- Fig. 3e astrocytes (3D) -------------------------------------------------
cat("== Fig. 3e astrocytes (3D) ==\n")
as <- fread(F_AS, colClasses = list(character = "star"))
stopifnot("fig3d astro table columns differ from 3D" = identical(names(as), c("signature", "Region", "n_CON", "n_NHD", "med_CON", "med_NHD", "p_wilcox", "cliff_d", "q_BH", "abs_d", "passes_effect", "star")),
          nrow(as) == 22L, is.logical(as$passes_effect), isTRUE(all.equal(as$abs_d, abs(as$cliff_d))))
b_as <- .block("Fig. 3e", "Astro", as$signature, as$Region, as$n_CON, as$n_NHD, as$cliff_d, as$p_wilcox, as$q_BH,
               as$passes_effect, .is_qstar(as$star), unique(as$signature))   # 3D writes ALL_SIG_ORDER

# ---- Fig. 4d oligodendrocyte lineage (3Fp) -----------------------------------
cat("== Fig. 4d oligodendrocyte lineage (3Fp) ==\n")
ol <- fread(F_OL, colClasses = list(character = "label"))
stopifnot("oligo lineage table columns differ from 3Fp" = identical(names(ol), c("signature", "Region", "lineage", "p", "q_BH", "cliff_d", "cliff_d_depth_matched", "ratio_depth_matched", "n_depth_matched", "passes", "label", "n_CON", "n_NHD")),
          nrow(ol) == 32L, is.logical(ol$passes), setequal(unique(ol$lineage), c("Oligodendrocyte", "OPC")))
# The panel plots cliff_d (all nuclei); cliff_d_depth_matched is the 3Fp sensitivity arm and is not shipped here
b_ol <- .block("Fig. 4d", c(Oligodendrocyte = "Oligo", OPC = "OPC")[ol$lineage], ol$signature, ol$Region, ol$n_CON, ol$n_NHD, ol$cliff_d, ol$p, ol$q_BH,
               ol$passes, .is_qstar(ol$label), unique(ol$signature))          # 3Fp writes oligo_programme_order()

# ---- Fig. 5d neurons (4C2) ---------------------------------------------------
cat("== Fig. 5d neurons (4C2) ==\n")
ne <- fread(F_NE, colClasses = list(character = "label"))
stopifnot("fig5 neuron table columns differ from 4C2" = identical(names(ne), c("signature", "Region", "class", "p", "cliff_d", "n_CON", "n_NHD", "q_BH", "passes", "label", "dotsize")),
          nrow(ne) == 20L, is.logical(ne$passes), setequal(unique(ne$class), c("Ex", "Inh")), setequal(unique(ne$signature), names(NEURON_LAB)))
# 4C2's group_by() wrote the CSV in alphabetical key order; the panel's y axis is ORD (Presynaptic, Postsynaptic, GluR, OxPhos, IEG) — used here
b_ne <- .block("Fig. 5d", c(Ex = "Neuron_Ex", Inh = "Neuron_Inh")[ne$class], unname(NEURON_LAB[ne$signature]), ne$Region, ne$n_CON, ne$n_NHD, ne$cliff_d, ne$p, ne$q_BH,
               ne$passes, .is_qstar(ne$label), unname(NEURON_LAB))

# ---- stack, check, derive ----------------------------------------------------
cat("== stack + checks ==\n")
d <- rbindlist(list(b_mi, b_as, b_ol, b_ne))
d[, figure_panel := factor(figure_panel, PANEL_LEV)]; d[, cell_type := factor(cell_type, CT_LEV)]; d[, region := factor(region, REGION_LEV)]
stopifnot(!anyNA(d$figure_panel), !anyNA(d$cell_type), !anyNA(d$region))
setorder(d, figure_panel, cell_type, prog_ord, region)
d[, `:=`(figure_panel = as.character(figure_panel), cell_type = as.character(cell_type), region = as.character(region))]

# (1) BH within panel: recomputing p.adjust over the raw p of each panel must give the source q (asserts the README statement)
.bh <- d[, .(ok = isTRUE(all.equal(p.adjust(p_raw, method = "BH"), q_raw, tolerance = 1e-6))), by = figure_panel]
stopifnot("a source q_BH is not BH over its own panel's p values — README statement would be false" = all(.bh$ok))

# (2) p / q underflow floor (ST0 convention); logged, never silent
n_p0 <- sum(d$p_raw < P_FLOOR); n_q0 <- sum(d$q_raw < P_FLOOR)
d[, p_wilcox := pmax(p_raw, P_FLOOR)]; d[, q_BH := pmax(q_raw, P_FLOOR)]
cat(sprintf("p-value underflow: %d p_wilcox and %d q_BH values written as 0 by the panel scripts, floored to %.7e (read as p < 2.2e-308)\n", n_p0, n_q0, P_FLOOR))
stopifnot(all(d$p_wilcox >= P_FLOOR), all(d$q_BH >= P_FLOOR), all(d$p_wilcox <= 1), all(d$q_BH <= 1), all(d$q_BH >= d$p_wilcox - 1e-12))

# (3) the house gate and mark, computed here on the unrounded delta and asserted against every source's own column
d[, effect_gate := abs(cliffs_delta_raw) >= GATE]
d[, q_star := q_BH < ALPHA & abs(cliffs_delta_raw) < GATE]
stopifnot("effect_gate disagrees with the source table's passes / passes_effect column" = identical(d$effect_gate, d$src_gate),
          "q_star disagrees with the source table's q* label / star column" = identical(d$q_star, d$src_qstar),
          "q_star and effect_gate must be mutually exclusive" = !any(d$q_star & d$effect_gate))
d[, cliffs_delta := round(cliffs_delta_raw, 3)]
# the shipped 3-dp delta must tell the same gate story as the unrounded one (no row sits on the rounding edge of 0.15)
stopifnot("a delta rounds across the 0.15 gate — ship more decimals" = identical(abs(d$cliffs_delta) >= GATE, d$effect_gate))

# (4) shape
stopifnot("expected 14 + 22 + 32 + 20 = 88 rows" = nrow(d) == 88L,
          !anyDuplicated(d[, .(figure_panel, cell_type, program, region)]),
          all(d$region %in% REGION_LEV), all(d$cliffs_delta >= -1 & d$cliffs_delta <= 1), !anyNA(d$cliffs_delta),
          "every program must have both regions" = all(d[, .N, by = .(figure_panel, cell_type, program)]$N == 2L),
          "n columns: NA only where the source did not record them (Fig. 2d)" =
            all(is.na(d[figure_panel == "Fig. 2d"]$n_CON)) && !anyNA(d[figure_panel != "Fig. 2d", .(n_CON, n_NHD)]),
          all(d[!is.na(n_CON)]$n_CON > 0), all(d[!is.na(n_NHD)]$n_NHD > 0))
# per-panel n per region x cell type must be constant across programs (same nuclei scored for every program)
.nchk <- d[!is.na(n_CON), .(k = uniqueN(paste(n_CON, n_NHD))), by = .(figure_panel, cell_type, region)]
stopifnot("n_CON / n_NHD vary across programs within a panel x cell type x region" = all(.nchk$k == 1L))
e <- d[, ..DATA_COLS]
stopifnot(identical(names(e), DATA_COLS), ncol(e) == 11)

# ---- sheet By_program (wide): both regions side by side ----------------------
cat("== sheet By_program ==\n")
w <- dcast(d, figure_panel + cell_type + program + prog_ord ~ region, value.var = c("cliffs_delta", "q_BH"))
stopifnot("dcast did not emit the <value>_<region> names expected" = all(c("cliffs_delta_Frontal", "cliffs_delta_Hippocampus", "q_BH_Frontal", "q_BH_Hippocampus") %in% names(w)))
# same sign in both regions; a delta of exactly 0 (none here) would be neither sign, so it is asserted absent
stopifnot("a shipped delta is exactly 0 — direction undefined" = all(w$cliffs_delta_Frontal != 0 & w$cliffs_delta_Hippocampus != 0))
w[, direction_consistent := sign(cliffs_delta_Frontal) == sign(cliffs_delta_Hippocampus)]
w[, figure_panel := factor(figure_panel, PANEL_LEV)]; w[, cell_type := factor(cell_type, CT_LEV)]
setorder(w, figure_panel, cell_type, prog_ord)
w[, `:=`(figure_panel = as.character(figure_panel), cell_type = as.character(cell_type), prog_ord = NULL)]
setcolorder(w, c("figure_panel", "cell_type", "program", "cliffs_delta_Frontal", "cliffs_delta_Hippocampus", "q_BH_Frontal", "q_BH_Hippocampus", "direction_consistent"))
stopifnot("expected 7 + 11 + 16 + 10 = 44 rows" = nrow(w) == 44L, ncol(w) == 8, !anyNA(w),
          !anyDuplicated(w[, .(figure_panel, cell_type, program)]), is.logical(w$direction_consistent),
          "By_program rows must follow the Data sheet order" =
            isTRUE(all.equal(w[, .(figure_panel, cell_type, program)], unique(e[, .(figure_panel, cell_type, program)]), check.attributes = FALSE)))

# ---- write (temp then rename; the shipped folder is never left half-written) ----
.wr <- function(x, p) { tmp <- paste0(p, ".tmp"); fwrite(x, tmp); invisible(file.rename(tmp, p)) }
.wr(e, OUT_A); .wr(w, OUT_B)
# copy check: every shipped delta / p / q must equal the source value (only names, region tokens, rounding and the floor changed)
.eb <- fread(OUT_A); .wb <- fread(OUT_B)
.chk <- function(src_prog, src_reg, src_d, src_p, src_q, panel, ct) {
  s <- data.table(program = src_prog, region = .norm_region(src_reg), cell_type = ct, d0 = as.numeric(src_d), p0 = as.numeric(src_p), q0 = as.numeric(src_q))
  m <- merge(.eb[figure_panel == panel], s, by = c("cell_type", "program", "region"))
  stopifnot(nrow(m) == nrow(s), isTRUE(all.equal(m$cliffs_delta, round(m$d0, 3))),
            isTRUE(all.equal(m$p_wilcox, pmax(m$p0, P_FLOOR))), isTRUE(all.equal(m$q_BH, pmax(m$q0, P_FLOOR))))
}
.chk(mi$program, mi$Region, mi$cliff_d, mi$p_wilcox, mi$q_BH, "Fig. 2d", "Micro-PVM")
.chk(as$signature, as$Region, as$cliff_d, as$p_wilcox, as$q_BH, "Fig. 3e", "Astro")
.chk(ol$signature, ol$Region, ol$cliff_d, ol$p, ol$q_BH, "Fig. 4d", c(Oligodendrocyte = "Oligo", OPC = "OPC")[ol$lineage])
.chk(unname(NEURON_LAB[ne$signature]), ne$Region, ne$cliff_d, ne$p, ne$q_BH, "Fig. 5d", c(Ex = "Neuron_Ex", Inh = "Neuron_Inh")[ne$class])
stopifnot(nrow(.eb) == 88L, nrow(.wb) == 44L, identical(names(.eb), DATA_COLS),
          is.logical(.eb$effect_gate), is.logical(.eb$q_star), is.logical(.wb$direction_consistent))   # fread round-trips TRUE/FALSE as logical

# ---- summary + provenance -----------------------------------------------------
summ <- c(
  sprintf("ST24  Data:       %d rows x %d columns -> %s", nrow(e), ncol(e), OUT_A),
  sprintf("ST24b By_program: %d rows x %d columns -> %s", nrow(w), ncol(w), OUT_B),
  sprintf("p-value underflow floor %.7e applied to %d p_wilcox and %d q_BH values (source wrote 0)", P_FLOOR, n_p0, n_q0),
  "Data: rows / programs / effect_gate / q_star per panel x cell type:",
  capture.output(print(e[, .(rows = .N, programs = uniqueN(program), gated = sum(effect_gate), q_star = sum(q_star),
                             n_CON = if (anyNA(n_CON)) NA_character_ else paste(unique(n_CON), collapse = "/"),
                             n_NHD = if (anyNA(n_NHD)) NA_character_ else paste(unique(n_NHD), collapse = "/")), by = .(figure_panel, cell_type)])),
  sprintf("Data: cliffs_delta range %+.3f to %+.3f", min(e$cliffs_delta), max(e$cliffs_delta)),
  "By_program: direction_consistent per panel x cell type:",
  capture.output(print(w[, .(programs = .N, consistent = sum(direction_consistent)), by = .(figure_panel, cell_type)])),
  "By_program: the programs whose sign flips between regions:",
  capture.output(print(w[direction_consistent == FALSE, .(figure_panel, cell_type, program, cliffs_delta_Frontal, cliffs_delta_Hippocampus)])))
cat(summ, sep = "\n")
writeLines(c(sprintf("127_program_delta_table_FH.R  run %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             sprintf("inputs: %s", paste(sprintf("%s (mtime %s)", basename(c(F_MI, F_AS, F_OL, F_NE)), format(file.mtime(c(F_MI, F_AS, F_OL, F_NE)))), collapse = "; ")),
             "", summ, "", capture.output(sessionInfo())), LOG)
cat("log:", LOG, "\n=== DONE ===\n")
