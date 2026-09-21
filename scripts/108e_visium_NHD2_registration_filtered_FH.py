#!/usr/bin/env python3
# 108e_visium_NHD2_registration_filtered_FH.py — Visium barcode-map registration: filtered matrix.
# filtered_feature_bc_matrix (h5 + MTX) for NHD_Frontal2_registered: the raw Space Ranger h5 subset to the
# in_tissue barcodes of the registered map (108d), written in Space Ranger's own h5 layout
# (matrix/{barcodes,data,indices,indptr,shape,features/*} + root attrs) so that Seurat::Load10X_Spatial
# reads it exactly like any other section's outs. 108f (R) round-trips it.
# Python environment: ctm_env (python_env_ctm_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import os, gzip, io, h5py, numpy as np, scipy.sparse as sps, scipy.io
if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
ROOT = os.path.dirname(os.environ["NHD_PROJ"])
OUTS = os.path.join(ROOT, "Visium", "NHD_Frontal2_registered", "outs")
tp = [l.rstrip("\n").split(",") for l in open(os.path.join(OUTS, "spatial", "tissue_positions_list.csv"))]
keep = set(r[0] for r in tp if r[1] == "1")
with h5py.File(os.path.join(OUTS, "raw_feature_bc_matrix.h5"), "r") as f:
    m = f["matrix"]; bcs = m["barcodes"][:].astype(str); shape = m["shape"][:]
    M = sps.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=(shape[0], shape[1]))
    feats = {k: m["features"][k][:] for k in m["features"].keys()}
    attrs = dict(f.attrs)
col = np.array([i for i, b in enumerate(bcs) if b in keep]); F = M[:, col].tocsc(); F.sort_indices()
print(f"raw {M.shape} -> filtered {F.shape} (nnz {F.nnz})")
out = os.path.join(OUTS, "filtered_feature_bc_matrix.h5")
if os.path.exists(out): os.remove(out)
with h5py.File(out, "w") as g:
    for k, v in attrs.items():
        if k == "library_ids": v = np.array([b"NHD_Frontal2_registered"], dtype="|S18")
        g.attrs[k] = v
    mm = g.create_group("matrix")
    mm.create_dataset("barcodes", data=bcs[col].astype("S18"))
    mm.create_dataset("data", data=F.data.astype(np.int32), compression="gzip")
    mm.create_dataset("indices", data=F.indices.astype(np.int64), compression="gzip")
    mm.create_dataset("indptr", data=F.indptr.astype(np.int64))
    mm.create_dataset("shape", data=np.array(F.shape, dtype=np.int32))
    ff = mm.create_group("features")
    for k, v in feats.items(): ff.create_dataset(k, data=v)
# MTX folder
d = os.path.join(OUTS, "filtered_feature_bc_matrix"); os.makedirs(d, exist_ok=True)
buf = io.BytesIO(); scipy.io.mmwrite(buf, F.astype(np.int32), field="integer")
with gzip.open(os.path.join(d, "matrix.mtx.gz"), "wb") as z: z.write(buf.getvalue())
with gzip.open(os.path.join(d, "barcodes.tsv.gz"), "wt") as z: z.write("\n".join(bcs[col]) + "\n")
with gzip.open(os.path.join(d, "features.tsv.gz"), "wt") as z:
    for i in range(F.shape[0]): z.write(f"{feats['id'][i].decode()}\t{feats['name'][i].decode()}\t{feats['feature_type'][i].decode()}\n")
print("wrote", out, "and", d, "\n=== DONE ===")
