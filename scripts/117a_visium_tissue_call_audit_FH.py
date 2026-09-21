#!/usr/bin/env python3
# 117a_visium_tissue_call_audit_FH.py — Audit of the Space Ranger tissue call (in_tissue) per section.
# Why: the pathologists' H&E annotation places layer 1 on the pale pial rim of the control
# sections, and that rim carries no spots in Fig 6c: Space Ranger's automatic tissue detection called it
# out-of-tissue (CON_Frontal1 379, CON_Frontal2 1,002 barcodes with median 4,334 / 3,297 UMI). The object
# was built from filtered_feature_bc_matrix.h5, so layer 1 never entered the analysis in CON, while in NHD
# the sulcal layer 1 is inside the mask and is labelled WM / Vasc-immune. This script shows, per section,
# every barcode of the capture area coloured by log10 UMI (from raw_feature_bc_matrix.h5) over the hires
# H&E, marks the barcodes Space Ranger excluded, and tabulates the UMI / gene distribution of excluded vs
# included spots, so the recovery mask (117b) can be set on evidence and checked by eye.
#
# OUT: figures/_diagnostics/visium_tissue_call/<section>_tissue_call_audit.png
#      tables/visium_tissue_call_audit_FH.csv  (per section: n in/out, UMI quantiles, n out with >= 500 / 1000 UMI)
# Python environment: ctm_env (python_env_ctm_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import os, json, csv, numpy as np, h5py, scipy.sparse as sp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
ROOT = os.path.dirname(os.environ["NHD_PROJ"])
PROJ = os.path.join(ROOT, "NHD_frontal_hippo_rebuild")
OUTD = os.path.join(PROJ, "figures", "_diagnostics", "visium_tissue_call"); os.makedirs(OUTD, exist_ok=True)
SECS = {"CON_Frontal1": "CON_Frontal1", "CON_Frontal2": "CON_Frontal2", "NHD_Frontal1": "NHD_Frontal1", "NHD_Frontal2": "NHD_Frontal2_registered"}

def read_raw(h5):
    f = h5py.File(h5); g = f["matrix"]
    bc = g["barcodes"][:].astype(str); shape = g["shape"][:]
    M = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=(shape[0], shape[1]))
    genes = g["features/name"][:].astype(str)
    umi = np.asarray(M.sum(0)).ravel(); ngene = np.asarray((M > 0).sum(0)).ravel()
    mt = np.asarray(M[np.char.startswith(genes, "MT-"), :].sum(0)).ravel()
    return bc, umi, ngene, np.where(umi > 0, 100 * mt / np.maximum(umi, 1), 0)

rows_out = []
for sec, folder in SECS.items():
    outs = os.path.join(ROOT, "Visium", folder, "outs"); spd = os.path.join(outs, "spatial")
    sf = json.load(open(os.path.join(spd, "scalefactors_json.json"))); hs = sf["tissue_hires_scalef"]
    pos = {r[0]: r for r in csv.reader(open(os.path.join(spd, "tissue_positions_list.csv")))}
    pos.pop("barcode", None)
    bc, umi, ng, pmt = read_raw(os.path.join(outs, "raw_feature_bc_matrix.h5"))
    keep = np.array([b in pos for b in bc]); bc, umi, ng, pmt = bc[keep], umi[keep], ng[keep], pmt[keep]
    it = np.array([int(pos[b][1]) for b in bc]); pr = np.array([float(pos[b][4]) for b in bc]) * hs; pc = np.array([float(pos[b][5]) for b in bc]) * hs
    img = Image.open(os.path.join(spd, "tissue_hires_image.png")).convert("RGB")
    inn, out = umi[it == 1], umi[it == 0]
    q = lambda v, p: int(np.quantile(v, p)) if len(v) else 0
    rows_out.append(dict(section=sec, outs=folder, n_in=int((it == 1).sum()), n_out=int((it == 0).sum()),
                         umi_in_q10=q(inn, .1), umi_in_med=q(inn, .5), umi_out_q10=q(out, .1), umi_out_med=q(out, .5), umi_out_q90=q(out, .9),
                         n_out_ge200=int((out >= 200).sum()), n_out_ge500=int((out >= 500).sum()), n_out_ge1000=int((out >= 1000).sum()),
                         genes_out_med=q(ng[it == 0], .5), mito_out_med=round(float(np.median(pmt[it == 0])), 1) if (it == 0).any() else 0,
                         mito_in_med=round(float(np.median(pmt[it == 1])), 1)))
    print(rows_out[-1])
    fig, ax = plt.subplots(1, 3, figsize=(19, 6.5))
    ax[0].imshow(img); ax[0].scatter(pc[it == 1], pr[it == 1], s=3, c="tab:red", alpha=.35, lw=0); ax[0].scatter(pc[it == 0], pr[it == 0], s=6, facecolors="none", edgecolors="black", lw=.4)
    ax[0].set_title(f"{sec}: Space Ranger in_tissue (red) / excluded (black rings, n={int((it==0).sum())})"); ax[0].axis("off")
    ax[1].imshow(img); s = ax[1].scatter(pc, pr, s=5, c=np.log10(umi + 1), cmap="magma", vmin=1.5, vmax=4.3, lw=0)
    ax[1].scatter(pc[it == 0], pr[it == 0], s=9, facecolors="none", edgecolors="cyan", lw=.35)
    plt.colorbar(s, ax=ax[1], fraction=.03, label="log10 UMI (all barcodes)"); ax[1].set_title("UMI of every barcode; cyan ring = excluded"); ax[1].axis("off")
    b = np.linspace(0, 4.6, 47)
    ax[2].hist(np.log10(inn + 1), bins=b, color="tab:red", alpha=.5, label=f"in tissue (n={len(inn)})")
    ax[2].hist(np.log10(out + 1), bins=b, color="black", alpha=.5, label=f"excluded (n={len(out)})")
    for v in (200, 500, 1000): ax[2].axvline(np.log10(v), ls=":", c="grey")
    ax[2].set_xlabel("log10 UMI"); ax[2].set_ylabel("barcodes"); ax[2].legend(frameon=False); ax[2].set_title("UMI distribution")
    plt.tight_layout(); plt.savefig(os.path.join(OUTD, f"{sec}_tissue_call_audit.png"), dpi=130); plt.close()
with open(os.path.join(PROJ, "tables", "visium_tissue_call_audit_FH.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)
print("wrote", OUTD)
