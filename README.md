# Nasu-Hakola disease frontal cortex and hippocampus — analysis code

Code for Leon et al. (manuscript under review): single-nucleus RNA-seq of prefrontal cortex and hippocampus from one
Nasu-Hakola disease (NHD, *TYROBP* c.2T>C homozygous) donor and one control donor, with Visium spatial
transcriptomics of the prefrontal cortex. The scripts take the processed count matrices to every figure,
supplementary figure and Supplementary Data table of the paper.

## Layout

* `atlas_build/` — per-library quality control, doublet removal, Azimuth annotation and the first Harmony
  integration of the seven single-nucleus libraries (`NHD_QC_Azimuth+Harmony.R`).
* `scripts/` — the frontal + hippocampus analysis: the scripts behind every figure panel, supplementary figure
  and Supplementary Data table, with the pipeline steps they read from (nothing exploratory). `01_build_FH_atlas.R`
  builds the working atlas; `22_publication_theme_FH.R` holds palettes and the figure theme; `30`–`47` differential
  expression (per-nucleus MAST, pseudobulk DESeq2, the dual-method merge, signed rankings); `44*`, `110b` gene-set
  enrichment; `02_`, `68_`, `70b_` transcription-factor activity; `50_` ligand–receptor analysis; `60_`, `116*` neuron
  re-annotation and the L4 IT label transfer; `70_`, `105_`, `106_`, `99b/c` the cross-cohort comparison with Zhou et al.
  2023; `65_`–`74_`, `107`–`127` the Visium analyses; figure scripts are named by figure and panel (`06d`/`43` Fig. 1,
  `2*` Fig. 2, `3*` Figs 3–4, `4*`/`62`/`64` Fig. 5, `70_suppfig1`, `109` supplementary figures); `90_`, `100_`,
  `102_`, `125_`–`127_` the Supplementary Data tables and the data-deposition bundle.
  Files starting with `_` are helpers that other scripts `source()`.
* `visium/` — Harmony integration of the four sections, cluster annotation, the cell2location deconvolution
  (Python) and the Supplementary Fig. 4 assembler.
* `zhou_reprocessing/` — the Zhou et al. 2023 occipital-cortex cohort (GEO GSE190015): quality control, doublet
  removal, Azimuth annotation and Harmony integration (`01_zhou_QC_Azimuth_Harmony.R`, exported from the R
  Markdown notebook it was run as), the per-cell-type pseudobulk validation object (`04_zhou_cross_cohort_validation.R`)
  and the Supplementary Fig. 2a–d panels (`70_suppfig2_zhou_QC.R`).
* `reference_tables/` — two tables that are not Supplementary Data: the genome-wide GTRD transcription-factor
  survey (decoupleR ULM) and the unsigned GTRD regulon it used, with a README. The transcription-factor
  result of the paper (Fig. 2h) is the CollecTRI analysis shipped as Supplementary Data 8.
* `SCRIPT_INDEX.csv` — every script with its one-line purpose, in run order.
* `package_versions.csv`, `R_version.txt` — the R environment (R 4.5.1); `python_env_cell2loc_env.txt` and
  `python_env_ctm_env.txt` — the two Python environments as pip freezes.

## Requirements

* R 4.5.1 with the packages in `package_versions.csv` (the rows flagged `core` are the ones the analysis
  depends on directly: Seurat 5, harmony, Azimuth, scDblFinder, MAST, DESeq2 + apeglm, fgsea, msigdbr,
  decoupleR, CellChat, monocle3, Banksy, SpatialExperiment, ComplexHeatmap, ggplot2, patchwork, data.table,
  openxlsx, officer).
* Python: `cell2loc_env` (cell2location 0.1.5 / scvi-tools 1.3.3, scanpy, anndata; used for deconvolution and
  the Jorstad reference export) and `ctm_env` (numpy, scipy, pandas, Pillow; used for the H&E
  tissue mask, the registration check and the annotation digitisation). Recreate them from the freeze files,
  e.g. `conda create -n cell2loc_env python=3.10 && pip install -r python_env_cell2loc_env.txt`.
* Space Ranger 1.3.1 outputs for the four Visium sections and the CellBender-filtered h5 files for the seven
  single-nucleus libraries (see Data).

## Running

Every script resolves its paths from one environment variable:

```
export NHD_PROJ=/path/to/NHD_frontal_hippo_rebuild      # the project folder (holds atlas/, data/, tables/, figures/, scripts/)
export NHD_REFERENCES=/path/to/references                # holds Jorstad2023/Jorstad2023_DFC.h5ad (script 116a)
```

The parent of `NHD_PROJ` must hold `Visium/` (one folder per section with its Space Ranger `outs/`, plus
`VISIUM_OUTS_MANIFEST.json` written by `117b`) and `NHD_QC_harmony.rds`, the first-pass atlas written by
`atlas_build/NHD_QC_Azimuth+Harmony.R`. That first script reads the CellBender h5 files from
`$NHD_CELLBENDER` and writes to `$NHD_ATLAS_OUT` (both default to the cluster paths it was run with).
Scripts stop with a message if `NHD_PROJ` is unset. Copy `scripts/` into `$NHD_PROJ/scripts/` and
`visium/` next to the Visium data, then run in this order:

1. **Atlas** — `atlas_build/NHD_QC_Azimuth+Harmony.R`, then `scripts/01_build_FH_atlas.R`, `60_neuron_reannotation_FH.R`,
   `116a_export_Jorstad_DFC_reference_FH.py` and `116b_L4_label_transfer_DFC_FH.R`.
2. **Differential expression** — `40_mast_percell_dual_FH.R`, `45_mast_subclass_FH.R`, `30_pseudobulk_wald_FH.R`,
   `41_merge_dual_method_FH.R`, `47_full_expressed_ranking_FH.R`, `47b_full_expressed_subclass_FH.R`;
   enrichment with `44_gsea_MAST_cliffs_FH.R` and `44b_gsea_subclass_FH.R`.
3. **Programs and networks** — the per-cell-type state and program caches (`2J0_`, `3F_`, `3M_`, `4C2_`), transcription
   factors (`02_`, `68_`), ligand–receptor (`50_`), the Zhou et al. comparison (`105_zhou_pseudobulk_cache_FH.R`,
   `70_zhou_crosscohort_FH.R`, `106_zhou_crosscohort_arms_FH.R`).
4. **Figures** — the panel scripts; each reads only the tables and caches written above and can be re-run alone.
   Per-cell caches carry modification-time guards and rebuild themselves when the atlas changes.
5. **Tables** — `90_build_supplementary_tables_FH.R`, `126_`, `127_`, then `100_package_supplementary_data_FH.R`.

Visium, in order: tissue call and registration (`107_`, `108a`, `108c`–`108e`, `117a`, `117b`) → integration
(`visium/integrated_harmony.R`, `visium/annotate_clusters.R`, `visium/refresh_identity.R`, `visium/propose_identity_20260918.R`)
→ BANKSY layer I / pia (`119_`, `119b_`, `121_`, `120_`, `118_`) → cell2location (`visium/01e_prep_NHD.R`, `01f_prep_visium_query.R`,
`02e_run_NHD.py`, `02f_run_NHD_spatial.py`, `03e_export_NHD.py`) → figures and tables (`65_`–`74_`, `108b`, `111_`, `114_`, `115_`,
`123_`, `125_`, `126_`). `112_visium_downstream_rebuild_FH.sh` runs the downstream Visium steps in dependency order.

Zhou et al. cohort (`zhou_reprocessing/`, run with `NHD_PROJ` set and `NHD_ZHOU_DIR` pointing at the folder that holds the
GEO downloads): `01_zhou_QC_Azimuth_Harmony.R`, then `04_zhou_cross_cohort_validation.R`, then `scripts/105_`, `70_`, `106_`;
`70_suppfig2_zhou_QC.R` draws Supplementary Fig. 2a–d.

Stochastic steps (clustering, UMAP, bootstraps, permutations, subsampling) set their seeds in the scripts.

## Data

Processed data — the single-nucleus count matrices with per-nucleus metadata (55,354 nuclei), the four Visium
Space Ranger outputs with spot metadata, and the raw sequencing reads — are deposited in Synapse, project
[syn77518970](https://www.synapse.org/#!Synapse:syn77518970) (doi:10.7303/syn77518970). The raw reads carry a
controlled-access requirement; the processed data are released with the paper. The Supplementary Data
workbooks that the table scripts produce are distributed with the article.

External references (downloaded, not shipped): the Jorstad et al. 2023 dorsolateral prefrontal cortex
snRNA-seq reference from CELLxGENE, collection `d17249d2-0e6e-4500-abb8-e6c93fa1ac6f` (raw counts;
`scripts/116a_export_Jorstad_DFC_reference_FH.py` documents the subsampling), and the Zhou et al. 2023 NHD
occipital-cortex atlas from GEO (GSE190015).

## Licence

MIT — see `LICENSE`.

## How to cite

Please cite the article (Leon et al.) and this repository; `CITATION.cff` carries the
reference and the archive DOI.
