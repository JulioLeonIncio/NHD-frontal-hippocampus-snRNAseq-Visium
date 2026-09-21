#!/usr/bin/env Rscript
# =============================================================================
# _artifact_genes.R — Ambient / artifact gene filter
# -----------------------------------------------------------------------------
# A per-nucleus reviewer
# (Cobos) flags this immediately. Rationale: SoupX (Young & Behjati 2020,
# PMID 33367645) / CellBender (Fleming 2023). This is an interim curated filter
# (a full SoupX re-estimation is future work); applied to glial comparisons only
# (neuronal genes are real biology in neurons).
#
# is_artifact(gene, cell_type) -> logical:
#   - lncRNA normalization/ambient artifacts + unannotated clone contigs: all cells
#   - curated neuronal/synaptic ambient markers: glial cells only
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

# Micro-PVM / Oligo / OPC discovery lists.
ARTIFACT_NEURONAL <- c(
  "RBFOX1","RBFOX3","KCNIP4","KCNIP1","CSMD1","SYT1","SNAP25","SYP","NRGN",
  "SLC17A7","GAD1","GAD2","LRRTM4","LRRN1","LRFN5","OPCML","NRG3","NRXN1",
  "NRXN3","NLGN1","DLGAP1","DLGAP2","FGF14","FGF12","DPP10","FAM155A",
  "STXBP5L","CTNNA2","ROBO2","RALYL","RIMS2","AGBL4","TENM2","SNTG1","ATRNL1",
  "CCSER1","NEGR1","ANKS1B","DSCAM","MDGA2","NEDD4L","ZFPM2","CELF2","CACNA1A",
  "CACNA1C","CNTNAP2","KALRN","NTNG1","PTPRN2","SYNDIG1","KHDRBS3","PLCB1",
  # Added R2: clean soup pattern (CON-glia pct ~0), no glial program.
  # GRIA4 deliberately not added — genuine OPC AMPA-receptor biology.
  "GRID2","GRIK2","ADGRL3","OTUD7A","GALNT13","PPP2R2B","NAV3","CACNB4")

# lncRNA whose cross-condition "DE" is a library-size / ambient normalization
# artifact (MALAT1/NEAT1 near-saturation; MEG3 neuronal imprinted lncRNA).
ARTIFACT_LNCRNA <- c("MALAT1","NEAT1","MEG3","SNHG14","DLEU1","MIR4300HG","KCNQ1OT1",
                     "ZFPM2-AS1",
                     "MIR181A1HG","PVT1","LINC00578","XACT")
# lncRNA host-gene / antisense PATTERNS (catch the long tail of mir*HG / LINC###### /
# SNHG# / *-AS# that individual symbols miss).
ARTIFACT_LNCRNA_REGEX <- "^MIR[0-9].*HG$|^LINC[0-9]{4,6}$|^SNHG[0-9]+$|-AS[0-9]+$"

# Replication-dependent cluster histones (legacy HIST1/2/3* and new H_C_ names).
# Non-polyadenylated, S-phase-restricted -> ambient/internal-priming artifact in
# poly-dT snRNA-seq absent a proliferation signature (Marzluff 2008 PMID 18927579;
# nomenclature Seal 2022 PMID 36180920). Deliberately does not match the real
# replication-INDEPENDENT variants H1-0 / H2AZ1 / H2AZ2 / H2AX / MACROH2A* / H3-3A/B.
ARTIFACT_HISTONE_REGEX <- "^(HIST[0-9]|H1-[1-9][0-9]?$|H2A[CJ][0-9]|H2B[CU][0-9]|H3C[0-9]|H4C[0-9])"

# Ultra-abundant translation/ribosomal-adjacent housekeeping (RNA-content artifact).
# NOTE (Lesson 1): excludes FTL/FTH1 (real NHD iron biology), ACTB/TUBA1A (possible
# morphology biology), B2M (real DAM-1 marker) -> those are never auto-killed.
ARTIFACT_HOUSEKEEPING <- c("EEF1A1","EEF1A2","EEF1B2","EEF1G","EEF2",
                           "FAU","NACA","BTF3","TPT1","TMSB4X","TMSB10",
                           "H3F3A","H3F3B")
# Ribosomal protein + mitochondrial-encoded (abundance artifacts). ^RP[LS][0-9]
# spares RPSA (laminin receptor); ^MT- is mito-encoded (percent.mt already regressed).
ARTIFACT_RIBO_MITO_REGEX <- "^(RPL|RPS|MRPL|MRPS)[0-9]|^MT-"

# Sex-chromosome genes: in an N=1 vs N=1 sex-matched design (both donors same sex,
# ~55 yo) these are donor/library noise, never NHD biology. XIST scored as a
# "discovery" gene in 7 comparisons -> remove from all cells.
ARTIFACT_SEX <- c("XIST","RPS4Y1","DDX3Y","UTY","USP9Y","KDM5D","EIF1AY","NLGN4Y",
                  "TXLNGY","ZFY")

# Ensembl unannotated clone/contig symbols (e.g. AC012150.1, AL359091.1, FP236383.3)
ARTIFACT_CONTIG_REGEX <- "^(AC|AL|AP|AF|FP|BX|CR|CT|Z)[0-9]{5,}\\.[0-9]+$"

# Silently disabled the neuronal-ambient filter for those scripts: is_artifact()
# only applies ARTIFACT_NEURONAL when cell_type is in .GLIAL, and the new names were
# not in the list. Verified: is_artifact("KCNIP4","Micro-PVM") was TRUE while
# is_artifact("KCNIP4","Microglia") was FALSE, and the Frontal microglial discovery
# set went from 76 flagged to 32. Any future rename must be added here.
.GLIAL <- c("Micro-PVM","Microglia","PVM","Astro","Oligo","OPC")

is_artifact <- function(gene, cell_type) {
  gene <- as.character(gene); cell_type <- as.character(cell_type)
  amb_all   <- gene %in% ARTIFACT_LNCRNA | gene %in% ARTIFACT_SEX |
               gene %in% ARTIFACT_HOUSEKEEPING |
               grepl(ARTIFACT_CONTIG_REGEX, gene)   |
               grepl(ARTIFACT_HISTONE_REGEX, gene)  |
               grepl(ARTIFACT_RIBO_MITO_REGEX, gene)|
               grepl(ARTIFACT_LNCRNA_REGEX, gene)
  amb_glial <- cell_type %in% .GLIAL & gene %in% ARTIFACT_NEURONAL
  amb_all | amb_glial
}

# ============================================================================
# is_callout_excluded(gene, cell_type) — label-only exclusion for figure callouts
# ----------------------------------------------------------------------------
# Stricter than is_artifact(): applied only to which genes are labelled/ringed on
# a scatter (Fig 1 panel d), never to the plotted points or to any DE set. Drops
# non-coding / antisense / divergent-transcript symbols and specific lncRNAs that
# a reviewer reads as artifacts even though they pass the ambient point filter,
# plus cross-lineage ambient (a canonical marker of another lineage leaking into
# this cluster). Keeps everything is_artifact() already removes.
# ============================================================================
# non-coding / antisense / divergent-transcript symbol PATTERNS (label-only)
CALLOUT_EXCLUDE_REGEX <- "-DT$|-AS[0-9]*$|^MEG[0-9]|^PWRN|^LINC[0-9]|MIR[0-9].*HG"
CALLOUT_EXCLUDE_GENES <- c("GMDS-DT","PWRN1","MEG8")
# cross-lineage ambient: canonical marker of lineage X leaking into cluster Y ->
# never a headline callout for Y. keyed cell_type -> genes.
CALLOUT_EXCLUDE_CROSSLIN <- list(`Micro-PVM` = c("GLUL"))   # GLUL = astrocyte-canonical

is_callout_excluded <- function(gene, cell_type) {
  gene <- as.character(gene); cell_type <- as.character(cell_type)
  ex <- is_artifact(gene, cell_type) |
        grepl(CALLOUT_EXCLUDE_REGEX, gene) |
        gene %in% CALLOUT_EXCLUDE_GENES
  # cross-lineage (vectorised over the keyed list)
  for (ct in names(CALLOUT_EXCLUDE_CROSSLIN))
    ex <- ex | (cell_type == ct & gene %in% CALLOUT_EXCLUDE_CROSSLIN[[ct]])
  ex
}

if (sys.nframe() == 0) {
  cat(sprintf("ARTIFACT_NEURONAL: %d | ARTIFACT_LNCRNA: %d | contig regex: %s\n",
              length(ARTIFACT_NEURONAL), length(ARTIFACT_LNCRNA), ARTIFACT_CONTIG_REGEX))
  cat(sprintf("CALLOUT_EXCLUDE_REGEX: %s | genes: %s\n",
              CALLOUT_EXCLUDE_REGEX, paste(CALLOUT_EXCLUDE_GENES, collapse = ",")))
}
