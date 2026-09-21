#!/usr/bin/env python3
# 116a_export_Jorstad_DFC_reference_FH.py — Export the Jorstad et al. 2023 dorsolateral prefrontal cortex (DFC, 10x 3' v3, 5 donors) excitatory nuclei as a 10x-style MTX bundle for Seurat label transfer.
# Why: our cortical neuron subclasses come from Azimuth's human motor-cortex reference, which has no
# L4 IT class (M1 is agranular). Jorstad 2023 DFC carries a within-area and a cross-area taxonomy that
# both contain "L4 IT" (6,734 / 5,059 of 63,366 excitatory nuclei), region-matched to our prefrontal
# samples. Source: cellxgene collection d17249d2-0e6e-4500-abb8-e6c93fa1ac6f, dataset
# "Dissection: Dorsolateral prefrontal cortex (DFC)" (raw counts in raw/X, gene symbols in feature_name).
# The 2.7 GB h5ad is read from $NHD_REFERENCES/Jorstad2023/Jorstad2023_DFC.h5ad.
#
# What: excitatory nuclei only, stratified subsample (cap CAP per WithinArea_subclass, seed 42) so the
# reference stays light (~27k nuclei) while every class incl. L4 IT is complete or near-complete.
# OUT: data/_reference_Jorstad2023_DFC/{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz, meta.csv, README.md}
# Python environment: cell2loc_env (python_env_cell2loc_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import os, gzip, numpy as np, h5py, scipy.sparse as sp, scipy.io as sio
if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
if "NHD_REFERENCES" not in os.environ:
    raise SystemExit("Set the NHD_REFERENCES environment variable to the folder that holds Jorstad2023/ (see README.md)")
H5   = os.path.join(os.environ["NHD_REFERENCES"], "Jorstad2023", "Jorstad2023_DFC.h5ad")
PROJ = os.environ["NHD_PROJ"]
OUT  = os.path.join(PROJ, "data", "_reference_Jorstad2023_DFC"); os.makedirs(OUT, exist_ok=True)
CAP  = 6000; rng = np.random.default_rng(42)

f = h5py.File(H5); obs = f["obs"]
def cat(k):
    g = obs[k]; return g["categories"][:].astype(str)[g["codes"][:]]
Class, wa, ca, donor, layer = cat("Class"), cat("WithinArea_subclass"), cat("CrossArea_subclass"), cat("donor_id"), cat("Layer")
idx_all = obs["_index"][:].astype(str)
exc = np.where(Class == "excitatory")[0]
keep = []
for s in np.unique(wa[exc]):
    i = exc[wa[exc] == s]
    keep.append(i if len(i) <= CAP else rng.choice(i, CAP, replace=False))
keep = np.sort(np.concatenate(keep))
print(f"excitatory {len(exc)} -> kept {len(keep)}")
for s in np.unique(wa[keep]): print(f"  {s:12s} {int((wa[keep]==s).sum()):6d}  (of {int((wa[exc]==s).sum())})")

X = f["raw/X"]; indptr = X["indptr"][:]
rows = []
for i in keep:
    a, b = indptr[i], indptr[i + 1]
    rows.append((X["indices"][a:b], X["data"][a:b]))
ngene = f["raw/X"].attrs["shape"][1]
ip = np.concatenate([[0], np.cumsum([len(r[0]) for r in rows])])
M = sp.csr_matrix((np.concatenate([r[1] for r in rows]).astype(np.int32), np.concatenate([r[0] for r in rows]), ip), shape=(len(keep), ngene))
assert M.data.min() >= 1 and np.all(M.data == np.round(M.data)), "raw/X is not integer counts"
genes_ens = f["raw/var/_index"][:].astype(str)
fn = f["raw/var/feature_name"]; genes_sym = fn["categories"][:].astype(str)[fn["codes"][:]] if isinstance(fn, h5py.Group) else fn[:].astype(str)

with gzip.open(os.path.join(OUT, "matrix.mtx.gz"), "wb") as g: sio.mmwrite(g, M.T.tocsc(), field="integer")
with gzip.open(os.path.join(OUT, "features.tsv.gz"), "wt") as g:
    for e, s in zip(genes_ens, genes_sym): g.write(f"{e}\t{s}\tGene Expression\n")
with gzip.open(os.path.join(OUT, "barcodes.tsv.gz"), "wt") as g:
    for b in idx_all[keep]: g.write(b + "\n")
with open(os.path.join(OUT, "meta.csv"), "w") as g:
    g.write("barcode,WithinArea_subclass,CrossArea_subclass,donor_id,Layer\n")
    for i in keep: g.write(f"{idx_all[i]},{wa[i]},{ca[i]},{donor[i]},{layer[i]}\n")
with open(os.path.join(OUT, "README.md"), "w") as g:
    g.write("# Jorstad 2023 DFC excitatory reference (subsampled)\n\nSource: cellxgene d17249d2-0e6e-4500-abb8-e6c93fa1ac6f, "
            "'Dissection: Dorsolateral prefrontal cortex (DFC)', raw counts. Excitatory nuclei, WithinArea_subclass capped at "
            f"{CAP} per class (seed 42): {len(keep)} of {len(exc)}. Written by scripts/116a_export_Jorstad_DFC_reference_FH.py.\n")
print("wrote", OUT, M.shape)
