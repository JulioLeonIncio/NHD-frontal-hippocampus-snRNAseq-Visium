#!/usr/bin/env Rscript
# =============================================================================
# 132_supp_data_12_communication_FH.R — Source data for the cell-cell communication and composition panels, which no other Supplementary Data workbook carries:
#   Fig. 2c, Fig. 3g, Fig. 3h, Fig. 4h, Fig. 4j, Fig. 5e and Fig. 5g.
#   Reads the panel tables written by the panel scripts themselves (never recomputed
#   here, so a number in a sheet is the number the panel drew) and writes them as
#   ST25-ST31 for the packager.
# OUT: supplementary_tables/ST25..ST31_*.csv
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
TAB <- file.path(PROJ, "tables"); OUT <- file.path(PROJ, "supplementary_tables")
rd <- function(...) { f <- file.path(...); stopifnot("missing panel table" = file.exists(f)); fread(f) }
say <- function(...) cat(sprintf(...), "\n")

# ---- Fig. 2c compensatory routes --------------------------------------------------
cr <- rd(TAB, "micro_states", "compensatory_routes_FH.csv")
stopifnot(identical(names(cr), c("region","gene","CON","NHD","raw_pp","delta_pp","fitted","route","below_floor")))
setnames(cr, "region", "Region")
fwrite(cr, file.path(OUT, "ST25_compensatory_routes.csv")); say("ST25 compensatory routes: %d rows", nrow(cr))

# ---- Fig. 3g / 3h / 5e chords: pathway totals per region x condition ----------------
# The three chord tables differ only in column spelling (src/tgt) and in whether the panel
# kept the share it drew. share is therefore recomputed for all three on one rule and
# ASSERTED against the one table that already carries it, so the sheet is internally consistent.
# 3A writes <panel>_stats.csv, so the Fig. 3g table is F3g1_micro_astro_chord_stats.csv after the next
# 3A run; the legacy 3a_ name is read only while the new one does not exist. Either way the table must be
# no older than the panel it backs, or a re-render would ship stale numbers without warning.
.f3g <- file.path(TAB, c("F3g1_micro_astro_chord_stats.csv", "3a_micro_astro_chord_stats.csv"))
.f3g <- .f3g[file.exists(.f3g)][1]
.p3g <- file.path(PROJ, "figures", "Figure_3", "panels", "F3g1_micro_astro_chord.pdf")
stopifnot("Fig. 3g chord table missing" = !is.na(.f3g),
          "Fig. 3g chord table is older than the F3g1 panel (re-run 3A or check the table name)" =
            !file.exists(.p3g) || difftime(file.mtime(.p3g), file.mtime(.f3g), units = "mins") < 10)
chord <- rbindlist(list(
  cbind(figure_panel = "Fig. 3g", fread(.f3g)),
  cbind(figure_panel = "Fig. 3h", rd(TAB, "F3h1_astro_neuron_chord_stats.csv")[, .(Region, Condition, source, target, pathway_name, prob_sum, weight, share_panel = share)]),
  cbind(figure_panel = "Fig. 5e", { x <- rd(TAB, "fig4E_micro_neuron_chord_stats.csv"); setnames(x, c("src","tgt"), c("source","target")); x })
), fill = TRUE)
chord[, total := sum(prob_sum), by = .(figure_panel, Region, Condition)]
chord[, share := prob_sum / total]
.chk <- chord[!is.na(share_panel)]
stopifnot("recomputed share does not reproduce the share the Fig. 3h chord drew" =
            nrow(.chk) > 0 && max(abs(.chk$share - .chk$share_panel)) < 1e-9)
say("chord share check: %d Fig. 3h rows reproduce to %.1e", nrow(.chk), max(abs(.chk$share - .chk$share_panel)))
chord[, c("weight","share_panel") := NULL]          # weight is a log1p ribbon-width transform, not data
setcolorder(chord, c("figure_panel","Region","Condition","source","target","pathway_name","prob_sum","total","share"))
fwrite(chord, file.path(OUT, "ST26_chord_pathway_totals.csv")); say("ST26 chord pathway totals: %d rows", nrow(chord))

# ---- Fig. 3g / 4j microglial ligand -> glial receptor -------------------------------
LRCOLS <- c("Region","target","interaction_name_2","pathway_name","prob_CON","prob_NHD",
            "delta_prob","pval_CON","pval_NHD","q_BH_CON","q_BH_NHD","sig_any",
            "drawn_as","rec_pct_max","rec_below_floor")
# Fig. 4j is the data-driven panel (3N_micro_oligolineage_LR_datadriven_FH.R, F4j in the rename
# map), not the curated SPP1/GAS6/PROS1/PSAP table 3B writes for RECV_SET=oligo — that table backs
# no assembled panel. 3N selects on the raw permutation P (< 0.05 in either condition) and computes
# no BH adjustment, so q_BH is NA there.
# Fig. 4j ships every cell of its grid (pairs x OPC/Oligo x region), read from the grid
# table 3N writes: the dots are the "shown" rows of the reported-pair table, and the open crosses
# (keyed "receptor not detected") are cells CellChat did not report at P < 0.05 in either condition
# whose receptor is below the 10 % raw-count detection floor in the receiver. Cross rows carry the
# existing fill rule (prob 0 and P 1 = not reported), delta_prob NA and sig_any FALSE.
f4j <- rd(TAB, "figF4j_micro_oligolineage_LR_FH.csv")[status == "shown"]
g4j <- rd(TAB, "figF4j_micro_oligolineage_LR_grid_FH.csv")
stopifnot("Fig. 4j grid table and reported-pair table were not written by the same 3N run (re-run 3N)" =
            abs(as.numeric(difftime(file.mtime(file.path(TAB, "figF4j_micro_oligolineage_LR_grid_FH.csv")),
                                    file.mtime(file.path(TAB, "figF4j_micro_oligolineage_LR_FH.csv")), units = "mins"))) < 10)
stopifnot("a shipped Fig. 4j row is not drawn in the panel" = all(f4j$maxprob > 0))
stopifnot("Fig. 4j shown rows are not unique per pair x receiver x region" =
            !anyDuplicated(f4j[, .(interaction_name_2, target, Region)]))
stopifnot("Fig. 4j grid cells are not unique per pair x receiver x region" =
            !anyDuplicated(g4j[, .(interaction_name_2, target, Region)]),
          "Fig. 4j grid drawn_as outside dot | cross" = all(g4j$drawn_as %in% c("dot", "cross")),
          "Fig. 4j grid is not every pair x receiver x region" =
            nrow(g4j) == uniqueN(g4j$interaction_name_2) * uniqueN(g4j$target) * uniqueN(g4j$Region),
          "a Fig. 4j cross has its receptor above the 10 % floor" = all(g4j[drawn_as == "cross"]$rec_below_floor),
          "a Fig. 4j cross carries a probability" = all(g4j[drawn_as == "cross", prob_CON == 0 & prob_NHD == 0 & is.na(dprob)]))
# the dots of the grid are exactly the shown rows, number for number
.dm <- merge(g4j[drawn_as == "dot"], f4j, by = c("interaction_name_2", "target", "Region"), suffixes = c("", ".shown"))
stopifnot("Fig. 4j grid dots are not exactly the shown rows" = nrow(.dm) == nrow(f4j) && nrow(.dm) == sum(g4j$drawn_as == "dot"),
          "Fig. 4j grid dot probabilities differ from the shown rows" =
            identical(.dm$prob_CON, .dm$prob_CON.shown) && identical(.dm$prob_NHD, .dm$prob_NHD.shown) &&
            identical(.dm$dprob, .dm$dprob.shown) && identical(.dm$pathway_name, .dm$pathway_name.shown))
.gk <- function(x) paste(x$interaction_name_2, x$target, x$Region)
f4j <- f4j[, .(Region, target, interaction_name_2, pathway_name, prob_CON, prob_NHD,
               delta_prob = dprob, pval_CON = as.numeric(pval_CON), pval_NHD = as.numeric(pval_NHD),
               q_BH_CON = NA_real_, q_BH_NHD = NA_real_, sig_any, drawn_as = "dot",
               rec_pct_max = g4j$rec_pct[match(.gk(f4j), .gk(g4j))],
               rec_below_floor = g4j$rec_below_floor[match(.gk(f4j), .gk(g4j))])]
x4j <- g4j[drawn_as == "cross", .(Region, target, interaction_name_2, pathway_name,
               prob_CON = 0, prob_NHD = 0, delta_prob = NA_real_, pval_CON = 1, pval_NHD = 1,
               q_BH_CON = NA_real_, q_BH_NHD = NA_real_, sig_any = FALSE, drawn_as = "cross",
               rec_pct_max = rec_pct, rec_below_floor)]
# Fig. 3g (3B) draws only pairs CellChat reported (no crosses survive the all-regions drop; an
# unreported cell is left blank), so every shipped 3g row is a dot. 3B computes receptor detection
# for its gate but does not write it to fig3b_LR_stats_astro.csv, so rec_* is NA for 3g rows.
f3g <- rd(TAB, "fig3b_LR_stats_astro.csv")
mg <- rbindlist(list(
  cbind(figure_panel = "Fig. 3g", f3g[, `:=`(drawn_as = "dot", rec_pct_max = NA_real_, rec_below_floor = NA)]),
  cbind(figure_panel = "Fig. 4j", f4j),
  cbind(figure_panel = "Fig. 4j", x4j)))
stopifnot(identical(names(mg), c("figure_panel", LRCOLS)),
          "a shipped microglia-to-glia DOT row is not significant in either condition (the README says every dot row is)" =
            all(mg[drawn_as == "dot"]$sig_any),
          "a Fig. 4j cross row is flagged significant" = !any(mg[drawn_as == "cross"]$sig_any),
          "Fig. 4j rows are not the full grid" = sum(mg$figure_panel == "Fig. 4j") == nrow(g4j),
          "Fig. 4j rows are not unique per pair x receiver x region" =
            !anyDuplicated(mg[figure_panel == "Fig. 4j", .(interaction_name_2, target, Region)]),
          "Fig. 3g rows are not the 3B table unchanged" =
            isTRUE(all.equal(mg[figure_panel == "Fig. 3g", names(fread(file.path(TAB, "fig3b_LR_stats_astro.csv"))), with = FALSE],
                             fread(file.path(TAB, "fig3b_LR_stats_astro.csv")), check.attributes = FALSE)))
# the two values the Results cite to this sheet must be present and match the panel table
.cite <- function(p, tg, rg) mg[figure_panel == "Fig. 4j" & interaction_name_2 == p & target == tg & Region == rg, delta_prob]
stopifnot(isTRUE(all.equal(round(.cite("TENM4 - ADGRL3", "OPC", "Frontal"), 3), -0.113)),
          isTRUE(all.equal(round(.cite("PDGFB - PDGFRA", "OPC", "Hippo"), 3),  0.036)))
say("Fig. 4j rows: %d = %d dots + %d crosses (TENM4-ADGRL3 OPC Frontal %+.4f; PDGFB-PDGFRA OPC Hippo %+.4f; max cross rec_pct_max %.4f)",
    sum(mg$figure_panel == "Fig. 4j"), nrow(f4j), nrow(x4j),
    .cite("TENM4 - ADGRL3", "OPC", "Frontal"), .cite("PDGFB - PDGFRA", "OPC", "Hippo"), max(x4j$rec_pct_max))
fwrite(mg, file.path(OUT, "ST27_LR_microglia_to_glia.csv")); say("ST27 microglia->glia L-R: %d rows", nrow(mg))

# ---- Fig. 3h astrocytic ligand -> neuronal receptor ---------------------------------
an <- rd(TAB, "fig4G_astro_neuron_LR_bubble.csv")
fwrite(an, file.path(OUT, "ST28_LR_astrocyte_to_neuron.csv")); say("ST28 astrocyte->neuron L-R: %d rows, %d pairs", nrow(an), uniqueN(an$interaction_name_2))

# ---- Fig. 5e microglial ligand -> neuronal receptor ---------------------------------
mn <- rd(TAB, "fig4E_micro_neuron_LR_bubble.csv")
fwrite(mn, file.path(OUT, "ST29_LR_microglia_to_neuron.csv")); say("ST29 microglia->neuron L-R: %d rows", nrow(mn))

# ---- Fig. 5g the same pairs resolved by neuronal subtype ----------------------------
sb <- rd(TAB, "figF5g_micro_neuron_subtype_LR_stats.csv")
fwrite(sb, file.path(OUT, "ST30_LR_microglia_to_neuron_subtype.csv")); say("ST30 subtype-resolved L-R: %d rows", nrow(sb))

# ---- Fig. 4h oligodendrocyte state composition --------------------------------------
oc <- rd(TAB, "oligo_state_composition_FH.csv")
stopifnot(identical(names(oc), c("Region","Condition","oligo_state","n","frac")))
stopifnot("fractions do not sum to 1 within a region x condition" =
            all(abs(oc[, sum(frac), by = .(Region, Condition)]$V1 - 1) < 1e-8))
fwrite(oc, file.path(OUT, "ST31_oligo_state_composition.csv")); say("ST31 oligodendrocyte state composition: %d rows", nrow(oc))
cat("=== DONE ===\n")
