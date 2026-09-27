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
chord <- rbindlist(list(
  cbind(figure_panel = "Fig. 3g", rd(TAB, "3a_micro_astro_chord_stats.csv")),
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
            "delta_prob","pval_CON","pval_NHD","q_BH_CON","q_BH_NHD","sig_any")
mg <- rbindlist(list(
  cbind(figure_panel = "Fig. 3g", rd(TAB, "fig3b_LR_stats_astro.csv")),
  cbind(figure_panel = "Fig. 4j", rd(TAB, "fig3b_LR_stats_oligo.csv"))))
stopifnot(identical(names(mg), c("figure_panel", LRCOLS)))
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
