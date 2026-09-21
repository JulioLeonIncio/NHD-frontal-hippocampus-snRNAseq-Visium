# =============================================================================
# _astro_inline_sets_FH.R — Single source for the astrocyte-state gene sets that Figure 3e,f scores but that do not live in data/curated_signatures.rds.
# -----------------------------------------------------------------------------
# Sourced by:
#   3D_astro_state_violins_FH.R          (scores them per nucleus -> Fig. 3e,f, ST24 Fig. 3e rows)
#   90_build_supplementary_tables_FH.R   (ships them in ST6 / Supplementary Data 7)
# The panel and the supplementary table must carry the same membership, so the sets are
# declared once here rather than typed into each script (
# "Supplementary Data 7 must list every gene set that Data 27 scores").
#
# Names are the panel axis labels exactly as ST24 (Supplementary Data 7, sheet Program_deltas; figure_panel
# "Fig. 3e") carries them, so a referee can join ST6.signature to ST24.program directly.
#
# Curation notes: absent genes are dropped at scoring
# time by score_sig() and reported in that script's log. Metallothionein: MT2A/MT3
# headline retained, MT1H (absent in the FH atlas) dropped. Homeostatic astro: SLC1A3
# kept in the module average despite its Frontal per-gene direction quirk.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
astro_inline_sets <- function() {
  list(
    "Metallothionein (lost)" = c("MT2A","MT3","MT1G","MT1E","MT1F","MT1M","MT1X","SLC30A1","SLC39A14"),
    "Homeostatic astro"      = c("SLC1A2","SLC1A3","GLUL","AQP4","GPC5","NDRG2","GJA1","ATP1A2"),
    "Heat-shock (HSF1)"      = c("HSPA1A","HSPA1B","HSPB1","CRYAB","HSP90AA1","DNAJB1","HSPH1","BAG3"),
    "Pan-reactive"           = c("LCN2","STEAP4","S1PR3","TIMP1","GFAP","VIM","CP","SERPINA3","CD44","OSMR","CD109","HSPB1"),
    "Interferon (Hasel)"     = c("ISG15","IFIT1","IFIT3","IFITM3","STAT1","B2M","IRF7","OASL","USP18","BST2"),
    "Synaptogenic"           = c("SPARCL1","THBS1","THBS2","GPC4","GPC6","CHRDL1","NRXN1","CADM1"))
}

# The five PUBLISHED astrocyte sets of curated_signatures.rds that Fig. 3e,f also scores,
# keyed by their curated_signatures / ST6 name, valued by the panel label under which ST24
# reports them. The A1/A2 dichotomy is deprecated (Escartin 2021), so the reactive rows are
# labelled by provenance, never as A1/A2 biology.
astro_published_panel_names <- function()
  c(A1_reactive           = "Complement/IFN-reactive (Liddelow)",
    A2_reactive           = "Ischemic/S100A10-reactive (Zamanian)",
    DAA_Habib             = "DAA (Habib)",
    STAT3_targets         = "STAT3 targets",
    Zhou_NHD_astrocyte_up = "Zhou NHD astrocyte")
