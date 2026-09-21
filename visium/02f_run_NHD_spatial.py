#!/usr/bin/env python
# =============================================================================
# 02f_run_NHD_spatial.py — NHD cell2location.
# Reference NB regression (8 broad types, batch=SampleID) + spatial mapping of
# the 4 integrated frontal Visium sections -> per-spot q05 cell abundances.
# CPU run on Apple Silicon (no GPU). ~12-16 h (3000 spatial epochs).
# Inputs : Visium/cell2location/c2l_MAIN/{sc_ref_*, visium_*} (mtx bundles)
# 02f: spatial mapping step on the current query (01f export; NHD_Frontal2 from
# NHD_Frontal2_registered/outs) with the reference signatures from the reference-regression
# step (02e). Same recipe as 02e, spatial model only.
# Outputs: c2l_MAIN/{sp_model/, q05_cell_abundance.csv, visium_c2l.h5ad}
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
import os, sys
# Headless fix: cell2location.filter_genes draws a QC scatter; an interactive
# backend blocks forever with no display. Force Agg before pyplot import.
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import coo_matrix
import cell2location as c2l
from cell2location.models import RegressionModel, Cell2location

def load_mtx_bundle(prefix):
    """MatrixMarket coordinate (genes x cells, 1-indexed) + _genes/_barcodes/_meta.
    PERF FIX: scipy.io.mmread is pure-Python and hangs on tens of M nonzeros;
    parse triplets with the pandas C tokenizer and build the sparse matrix directly."""
    mtx = prefix + "_matrix.mtx"
    with open(mtx) as fh:
        _hdr = fh.readline()                      # %%MatrixMarket ...
        dims = fh.readline().split()              # nrow ncol nnz
    nrow, ncol = int(dims[0]), int(dims[1])
    trip = pd.read_csv(mtx, delim_whitespace=True, skiprows=2, header=None,
                       names=["i", "j", "v"],
                       dtype={"i": np.int32, "j": np.int32, "v": np.float32},
                       engine="c")
    Xgc = coo_matrix((trip["v"].values,
                      (trip["i"].values - 1, trip["j"].values - 1)),
                     shape=(nrow, ncol))
    X = Xgc.T.tocsr()                             # cells x genes
    del trip, Xgc
    genes = [l.strip() for l in open(prefix + "_genes.tsv")]
    bcs   = [l.strip() for l in open(prefix + "_barcodes.tsv")]
    meta  = pd.read_csv(prefix + "_meta.csv", index_col=0)
    meta.index = meta.index.astype(str)
    a = ad.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [str(b) for b in bcs]
    meta = meta.reindex(a.obs_names)
    for c in meta.columns:
        a.obs[c] = meta[c].values
    a.var_names_make_unique()
    return a

if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
NHD = os.path.dirname(os.environ["NHD_PROJ"])
C2L = os.environ.get("C2L_DIR", os.path.join(NHD, "Visium/cell2location/c2l_MAIN"))   # RIKEN hedge: C2L_DIR=<bundle dir>
os.makedirs(os.path.join(C2L, "sp_model"), exist_ok=True)

LABEL_KEY = "c2l_label"
BATCH_KEY = "SampleID"
N_CELLS_PER_SPOT = 10          # cortical Visium ~5-15 cells/spot
DETECTION_ALPHA  = 20          # cell2location default for controlled tech variation
REF_EPOCHS = 250
SP_EPOCHS  = 3000              # ~12-16 h CPU; ELBO plateau well before 3k for ~18k spots

def log(m): print(f"[c2l] {m}", flush=True)

# ---- 1. reference signatures: from the reference-regression step ------------
# The snRNA reference is the same atlas; ref_signatures.csv from the regression step is
# loaded instead of re-training the regression model (identical inf_aver, ~40 min saved).
log("loading reference signatures (reused from the 2026-06-13 run)")
inf_aver = pd.read_csv(os.path.join(C2L, "ref_signatures.csv"), index_col=0)
log(f"reference signatures: {inf_aver.shape} | types: {list(inf_aver.columns)}")

# ---- 2. spatial mapping -----------------------------------------------------
log("loading visium mtx bundle")
adata_vis = load_mtx_bundle(os.path.join(C2L, "visium"))
adata_vis.var_names_make_unique()
shared = [g for g in adata_vis.var_names if g in inf_aver.index]
log(f"shared genes ref∩visium: {len(shared)}")
adata_vis = adata_vis[:, shared].copy()
inf_aver  = inf_aver.loc[shared, :]

Cell2location.setup_anndata(adata=adata_vis)
mod_sp = Cell2location(adata_vis, cell_state_df=inf_aver,
                       N_cells_per_location=N_CELLS_PER_SPOT,
                       detection_alpha=DETECTION_ALPHA)
log(f"training spatial model ({SP_EPOCHS} epochs, CPU) — long pole")
mod_sp.train(max_epochs=SP_EPOCHS, batch_size=None, train_size=1, accelerator="cpu")
adata_vis = mod_sp.export_posterior(
    adata_vis, sample_kwargs={"num_samples": 1000, "batch_size": mod_sp.adata.n_obs,
                              "accelerator": "cpu"})

mod_sp.save(os.path.join(C2L, "sp_model"), overwrite=True)
q05 = adata_vis.obsm["q05_cell_abundance_w_sf"].copy()
q05.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in q05.columns]
q05.to_csv(os.path.join(C2L, "q05_cell_abundance.csv"))
adata_vis.write(os.path.join(C2L, "visium_c2l.h5ad"))
log(f"q05 abundance written: {q05.shape}")
print("\n=== DONE ===", file=sys.stderr)
