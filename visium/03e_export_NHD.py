#!/usr/bin/env python
# =============================================================================
# 03e_export_NHD.py — NHD cell2location step 3 — join per-spot abundance to spatial metadata.
# Port of SCZ 03e_export_MAIN_BANKSY.py.
# Output: Visium/cell2location/c2l_MAIN/abundance_q05.csv
#   columns: spot_id, <one col per cell type>, sample_id, condition,
#            seurat_clusters (domain), coords. Also writes proportions.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
import os, sys
import pandas as pd
import scanpy as sc

if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
NHD = os.path.dirname(os.environ["NHD_PROJ"])
C2L = os.path.join(NHD, "Visium/cell2location/c2l_MAIN")

q05 = pd.read_csv(os.path.join(C2L, "q05_cell_abundance.csv"), index_col=0)
adata_vis = sc.read_h5ad(os.path.join(C2L, "visium_c2l.h5ad"))

# Spot annotations come from the current query export (visium_meta.csv, written by the prep step from
# the integrated object) so that cluster labels always follow the object; the h5ad obs is the copy
# the model was trained with and carries no information the abundances depend on
vm = pd.read_csv(os.path.join(C2L, "visium_meta.csv"), index_col=0)
vm.index = vm.index.astype(str)
meta_cols = [c for c in ["sample_id", "condition", "seurat_clusters",
                         "Spatial_snn_res.0.5", "percent.mt",
                         "nCount_Spatial", "nFeature_Spatial",
                         "imagerow", "imagecol", "row", "col"]
             if c in vm.columns]
meta = vm.loc[adata_vis.obs_names, meta_cols].copy()

ct = list(q05.columns)
out = q05.join(meta, how="left")
out.index.name = "spot_id"
out.to_csv(os.path.join(C2L, "abundance_q05.csv"))

# per-spot proportions (sum-to-1 across cell types)
prop = q05[ct].div(q05[ct].sum(1), axis=0)
prop = prop.join(meta, how="left")
prop.index.name = "spot_id"
prop.to_csv(os.path.join(C2L, "proportions_q05.csv"))

print(f"[c2l] abundance_q05.csv written: {out.shape}", flush=True)
print(f"[c2l] cell types: {ct}", flush=True)
if {"condition"}.issubset(out.columns):
    print("\nMean proportion by condition:", flush=True)
    pc = prop[ct].copy(); pc["condition"] = prop["condition"].values
    print(pc.groupby("condition").mean().round(4).T, flush=True)
if {"seurat_clusters", "condition"}.issubset(out.columns):
    print("\nSpots per domain x condition:", flush=True)
    print(out.groupby(["seurat_clusters", "condition"]).size())
print("\n=== DONE ===", file=sys.stderr)
