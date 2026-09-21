# =============================================================================
# _oligo_programmes_FH.R — Single source for the oligodendrocyte-lineage programme gene sets used by Figure 4.
# -----------------------------------------------------------------------------
# Sourced by:
#   3M_prep_opc_and_scales_FH.R        (per-nucleus / per-spot myelin quantities)
#   3Fp_oligo_lineage_block_FH.R       (per-nucleus dot-matrix + supplement violins)
#   3P_visium_programme_reconcile_FH.R (the same programmes scored on Visium spots)
# The snRNA and spatial panels must measure the same gene sets or the reconciliation
# is meaningless, so they are declared once here rather than copied into each script.
#
# Every member's Cliff's delta
# was printed before shipping, and the old `Myelination` set failed badly: it carried
# two strongly down genes and OMITTED the one gene that is genuinely induced.
#
#   old set member   delta Frontal / Hippo   detection   verdict
#   PLP1                +0.638 / +0.224        0.98      on-direction
#   CNP                 +0.369 / +0.174        0.69      on-direction
#   MOG                 +0.296 / +0.211        0.63      on-direction
#   MBP                 +0.217 / +0.047        0.99      on-direction
#   MAG                 +0.160 / +0.235        0.43      on-direction
#   CLDN11              +0.158 / +0.045        0.58      on-direction
#   MOBP                +0.083 / +0.245        0.72      weak/on
#   mal                 -0.026 / -0.005        0.47      flat (kept: canonical)
# removed
# removed
#   GJB1                +0.018 / -0.016        0.11      near floor — removed
#   GJC2                    absent             0.00      not detected — removed
# added
#
# Pooling TF and ANLN against PLP1 and omitting OPALIN drove the module to +0.13/+0.03,
# so the figure's headline panel rendered the myelin arm as a near-white dot carrying
# the house "q*" triviality mark — asserting the opposite of what the data show.
# TF (iron/transferrin) and ANLN (paranodal cytoskeleton) are not noise: their fall is
# its own result and they are shown as individual genes in the gene-level panel, never
# averaged into structural myelin.
#
# The LIPID supply chain is split into its arms, not pooled. The arms behave
# differently and that specificity is the result:
#   sterol           fails in both regions          (depth-matched delta -0.38 / -0.38)
#   sphingolipid/FA  SCD -0.28/-0.61, SGMS1 -0.19/-0.32, ACER3 -0.26/-0.30
#   galactolipid     UGT8 +0.07/-0.08, FA2H 0.00/+0.05, GALC -0.01/+0.04  = preserved
#   uptake/salvage   SLC44A1 -0.17/-0.43 (choline for phosphatidylcholine, 99% detected),
#                    NPC1 -0.12/-0.13; LDLR -0.12 hippo (below detection in frontal)
# A pooled "myelin lipid" row averages a robust both-region deficit with a
# region-specific one and a null, and hides which parts of the chain break. The
# uptake arm closes the loop: the cell cannot synthesise sterol and is not
# compensating by importing it.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
oligo_programme_sets <- function(sigs) {
  list(
    # membrane proteins only — what the cell is ORDERING
    "Structural myelin"       = c("PLP1","MBP","MOBP","MOG","MAG","CNP","MAL",
                                  "CLDN11","OPALIN"),
    # --- the myelin-lipid supply chain, one row per arm ----------------------
    "Cholesterol (sterol arm)"= c("HMGCR","HMGCS1","SQLE","DHCR7","DHCR24","MSMO1",
                                  "FDFT1","LSS","CYP51A1","SC5D","IDI1","MVD","FDPS"),
    # It carries
    # no signal and its presence implies a ceramide-supply claim the data do not make.
    # ELOVL1, the VLCFA elongase relevant to myelin sphingolipids, is flat here
    # (-0.013/+0.029), so this arm is carried by MUFA desaturation (SCD) and
    # sphingomyelin synthesis (SGMS1); ACER3 is a ceramidASE (catabolic) and is kept
    # only because it moves with them — the arm is named for what it measures.
    "Fatty-acid / sphingomyelin" = c("SCD","SGMS1","ACER3","ELOVL5"),
    # GAL3ST1 (sulfatide synthase) is undetected in these nuclei and cannot be added.
    "Galactolipid"            = c("UGT8","FA2H","GALC"),
    "Lipid uptake / salvage"  = c("SLC44A1","NPC1","LDLR"),
    # The discriminator row. A coordinated fall of SREBF2 + INSIG1 +
    # LDLR + the sterol enzymes is also the canonical signature of sterol repletion sensed
    # at the er — which in a DAP12-null brain, where phagocytes cannot clear myelin
    # cholesterol, is a live alternative to synthesis failure. If the cell were
    # sterol-LOADED, LXR targets would rise. Measured: ABCG1 -0.01/-0.10, NPC1
    # -0.12/-0.13, SREBF1 -0.02/-0.07, MYLIP +0.02 — flat to down, i.e. the LXR arm is
    # not activated, which argues against repletion. Thin (ABCA1/APOE/NR1H3 undetected in
    # oligodendrocytes), so it is a control row, not a claim.
    "Sterol efflux (LXR)"     = c("ABCG1","NPC1","SREBF1","MYLIP","NPC2"),
    # --- context -------------------------------------------------------------
        # BCAS1 added: newly-formed / early-myelinating oligodendrocyte marker (Fard 2017),
    # up here (+0.175/+0.067, detection 0.48/0.54). With OPALIN and MYRF it makes a
    # three-gene "attempted terminal differentiation" phenotype rather than one gene.
    "Oligo differentiation"   = c("MYRF","BCAS1","SOX10","OLIG1","OLIG2","NKX2-2","ZEB2"),
    "Stress / reactive"       = c("CRYAB","HSPA1A","HSPB1","DNAJB1","FTL","FTH1",
                                  "APOD","CLU","B2M")
  )
}

# Genes deliberately not in "Structural myelin" but reported at gene level, with the
# direction that disqualified them from the module (see the audit above).
oligo_myelin_support_genes <- c("TF", "ANLN")

# The row order the figure reads in: the order is placed, then every arm that has to
# supply the membrane, then context. Consumers should use this rather than names().
oligo_programme_order <- function()
  c("Structural myelin",
    "Cholesterol (sterol arm)", "Fatty-acid / sphingomyelin", "Galactolipid",
    "Lipid uptake / salvage", "Sterol efflux (LXR)",
    "Oligo differentiation", "Stress / reactive")

# The three arms that form the bracketed "myelin lipid supply" group in the panels.
oligo_lipid_arms <- function()
  c("Cholesterol (sterol arm)", "Fatty-acid / sphingomyelin", "Galactolipid",
    "Lipid uptake / salvage")
