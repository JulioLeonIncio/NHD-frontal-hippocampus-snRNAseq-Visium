#!/usr/bin/env Rscript
# =============================================================================
# _pan2024_csf1r_signatures.R — Pan et al. 2024 (Acta Neuropathol Commun, PMC11365264) CSF1R-RD microglia / oligo / OPC state SIGNATURES, curated for the NHD cross-microgliopathy
# comparison. Marker symbols transcribed from the paper's Results +
# Fig 6/8 dot plots + trajectory heatmaps.
#
# Usage: source() this to get PAN2024_SIGS (named list of character vectors).
# Intended to be scored as added rows in the existing curated-state ora over our
# Micro-PVM pseudobulk (signature/direction level only — never cross-disease DEG
# correlation; both studies are single disease donors).
#
# KEY design choice: the CSF1R-RD M3 "autophagy/lipid-laden"
# state is split into two sub-signatures, because in our NHD data these arms
# behave OPPOSITELY (lysosomal partially up in Hippo; lipid/PPARG arm down in both
# regions) — the split is the headline of the comparison.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

PAN2024_SIGS <- list(

  # --- Microglia states (the focus) ---------------------------------------
  # M1 homeostatic (lost in disease). Pan's reported set:
  M1_homeostatic = c("CX3CR1","SOX5","FOXP2","GRID2","KHDRBS3","RASGEF1C","SYNDIG1"),

  # M2 pro-inflammatory. note: HSPH1/HSP90AA1 are heat-shock (generic stress) —
  # flagged; consider excluding when testing inflammation specificity.
  M2_proinflammatory = c("FCGBP","FCGR3A","IL2RA","IL18","VSIR","ADGRG6",
                         "MAMDC2","ARHGAP15","STARD13"),   # HSPH1, HSP90AA1 dropped (HSP)

  # M3 split — lysosomal / autophagy arm (partially induced in NHD):
  M3_lysosomal_autophagy = c("GPNMB","SQSTM1","ATG7","SCARB2","HEXA"),

  # M3 split — lipid-laden / PPARG arm (fails / reverses in NHD = the headline):
  M3_lipid_laden = c("PPARG","CPM","MGLL","ABHD3","PLPP3","LGALS1","LGALS3",
                     "HS3ST2","PTPRG"),

  # M4 peripheral monocyte-derived. caution: F13A1/MARCO overlap perivascular /
  # border-associated macrophage identity (Van Hove 2019, PMID 30894658) — in a
  # Micro-PVM cluster this is a PVM-compartment signal, not proof of monocyte
  # infiltration. Interpret accordingly.
  M4_monocyte_PVM = c("CD44","CCDC88C","MARCO","F13A1","ITGA4"),

  # M0 AD-associated/phagocytic (Pan's CSF1R-RD donor had AD co-pathology) —
  # included to control for AD-driven (not NHD-specific) signal.
  M0_AD_phagocytic = c("ARHGAP18","CDC42","CD163L1","CLEC7A","SLC11A1",
                       "CEBPD","IFI44L","FAM110B"),

  # --- Bonus: glial states (for the Fig 3 astro/oligo comparison, optional) --
  # O3 CSF1R-RD oligodendrocyte stress/apoptosis (HSP-heavy — mostly chaperone):
  O3_oligo_stress = c("DNAJB1","DNAJB6","CRYAB","PCBP1","BAG3","FOS",
                      "TNFRSF12A","SPP1","BEX1"),   # HSP90AA1/AB1/HSPA1A dropped (HSP)

  # OP1 CSF1R-RD OPC differentiation-ARREST — parallels our NHD oligo maturation
  # block (Fig 3). Up = arrest/immature; Down = pro-myelination/maturation.
  OP1_arrest_up   = c("LRRN1","LRRTM3","SCN1A","SCN3A","CDH10","LRP1",
                      "PDGFRA","SOX5","SOX6","NFIA"),
  OP1_arrest_down = c("NKX2-2","NKX6-2","SOX4","SOX8","TCF7L2","YY1","ZNF488",
                      "MBP","MOG","S100B","ITGB1")
)

# Canonical homeostatic anchors used in our manuscript (for the divergence
# scatter's "shared collapse" quadrant; not Pan-specific but the same biology):
PAN2024_HOMEOSTATIC_CANONICAL <- c("P2RY12","CX3CR1","TMEM119","CSF1R","SELPLG")

# Genes Pan reports as down in a state (sign matters for the divergence scatter):
PAN2024_DOWN <- list(
  M2_proinflammatory = c("ST6GALNAC3"),
  OP1_arrest         = c("NKX2-2","NKX6-2","SOX4","SOX8","TCF7L2","YY1",
                         "ZNF488","MBP","MOG","S100B","ITGB1")
)

if (sys.nframe() == 0) {
  cat("Pan 2024 CSF1R-RD signatures loaded:\n")
  for (nm in names(PAN2024_SIGS))
    cat(sprintf("  %-24s %2d genes: %s\n", nm, length(PAN2024_SIGS[[nm]]),
                paste(PAN2024_SIGS[[nm]], collapse=", ")))
}
