#!/usr/bin/env python3
# 108d_visium_NHD2_registration_outs_FH.py — Visium barcode-map registration: reflection + tissue mask.
# Builds Visium/NHD_Frontal2_registered/outs, the NHD_Frontal2 Space Ranger output with the barcode->pixel
# map registered to the scan and the tissue call derived from the H&E.
#
# Why this STEP. For NHD_Frontal2 the Space Ranger --reorient-images run (NHD_Frontal2_reorient) fits
# the fiducial frame but places the barcode layout in mirror image relative to the scan (a known
# failure mode of fiducial alignment on mirrored images): barcode (array_row r, array_col c)
# physically sits where the layout puts (77 - r, c), a reflection about the array's centre line
# (rows are the horizontal image axis in this scan; a quarter-row scan of the reflection axis peaks
# at 38.5-38.75 for Spearman(UMI, tissue brightness)). On the unreflected map the tissue call is
# mirrored as well: 356 barcodes on the tear / background carry in_tissue = 1 (median 1,318 UMI) and
# 373 barcodes over intact tissue carry in_tissue = 0 (median 2,984 UMI) — the "notch" that is not
# present in the H&E. 108c quantifies the mis-registration for every section.
#
# What it does. positions: (pxl_row, pxl_col) = lattice fit evaluated at (77 - r, c) (fit residual
# 0.5 px); in_tissue: >= 50 % of the 55-um footprint on H&E tissue (grey < background - 12,
# closed/opened), a mask that reproduces Space Ranger's own call at the delivered positions for
# 98.9 % of spots. Images: native scan orientation (unchanged). Matrices: raw copied verbatim;
# filtered rebuilt from raw with the registered in_tissue (MTX + h5 written by R, see 108e).
# Everything downstream of Space Ranger (integration, domains, cell2location, Figs 6/7,
# Supplementary Data) reads NHD_Frontal2 from this folder.
# Python environment: ctm_env (python_env_ctm_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import os, json, shutil, gzip, numpy as np
from PIL import Image
from scipy.ndimage import uniform_filter, binary_closing, binary_opening
if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
ROOT = os.path.dirname(os.environ["NHD_PROJ"])
SRC  = os.path.join(ROOT, "Visium", "NHD_Frontal2_reorient", "outs")
DST  = os.path.join(ROOT, "Visium", "NHD_Frontal2_registered", "outs")
os.makedirs(os.path.join(DST, "spatial"), exist_ok=True)
sp = os.path.join(SRC, "spatial")
sf = json.load(open(os.path.join(sp, "scalefactors_json.json"))); hs = sf["tissue_hires_scalef"]
rows = [l.rstrip("\n").split(",") for l in open(os.path.join(sp, "tissue_positions_list.csv"))]
bc = [r[0] for r in rows]; it = np.array([int(r[1]) for r in rows]); ar = np.array([int(r[2]) for r in rows]); ac = np.array([int(r[3]) for r in rows])
pr = np.array([float(r[4]) for r in rows]); pc = np.array([float(r[5]) for r in rows])
A = np.vstack([ar, ac, np.ones_like(ar)]).T
cx = np.linalg.lstsq(A, pc, rcond=None)[0]; cy = np.linalg.lstsq(A, pr, rcond=None)[0]
assert np.abs(A @ cx - pc).max() < 1 and np.abs(A @ cy - pr).max() < 1, "lattice fit is not exact"
X = cx[0] * (77 - ar) + cx[1] * ac + cx[2]; Y = cy[0] * (77 - ar) + cy[1] * ac + cy[2]
img = np.asarray(Image.open(os.path.join(sp, "tissue_hires_image.png")).convert("RGB")).astype(float)
grey = img.mean(2); bg = np.quantile(grey, 0.995)
mask = binary_opening(binary_closing(grey < bg - 12, iterations=2), iterations=1)
r = int(sf["spot_diameter_fullres"] * hs / 2); mf = uniform_filter(mask.astype(float), size=2 * r + 1)
def frac(Xp, Yp):
    xs = np.clip((Xp * hs).astype(int), 0, mf.shape[1] - 1); ys = np.clip((Yp * hs).astype(int), 0, mf.shape[0] - 1); return mf[ys, xs]
agree = ((frac(pc, pr) > 0.5) == (it == 1)).mean()
it_new = (frac(X, Y) > 0.5).astype(int)
print(f"mask vs Space Ranger call at delivered positions: {agree:.4f} agreement")
print(f"in_tissue: delivered {it.sum()} -> corrected {it_new.sum()} (newly in {((it_new==1)&(it==0)).sum()}, newly out {((it_new==0)&(it==1)).sum()})")
# write spatial/
for f in ["tissue_hires_image.png", "tissue_lowres_image.png", "detected_tissue_image.jpg", "aligned_fiducials.jpg", "scalefactors_json.json"]:
    shutil.copy2(os.path.join(sp, f), os.path.join(DST, "spatial", f))
with open(os.path.join(DST, "spatial", "tissue_positions_list.csv"), "w") as f:
    for b, i, a, c, y, x in zip(bc, it_new, ar, ac, Y, X):
        f.write(f"{b},{i},{a},{c},{int(round(y))},{int(round(x))}\n")
# raw matrices verbatim; filtered rebuilt by 108e (R) from raw + this in_tissue
for f in ["raw_feature_bc_matrix.h5", "metrics_summary.csv", "web_summary.html"]:
    if os.path.exists(os.path.join(SRC, f)): shutil.copy2(os.path.join(SRC, f), os.path.join(DST, f))
if os.path.isdir(os.path.join(SRC, "raw_feature_bc_matrix")):
    shutil.copytree(os.path.join(SRC, "raw_feature_bc_matrix"), os.path.join(DST, "raw_feature_bc_matrix"), dirs_exist_ok=True)
with open(os.path.join(DST, "in_tissue_registered.csv"), "w") as f:
    f.write("barcode,in_tissue_delivered,in_tissue_registered,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n")
    for b, i0, i1, a, c, y, x in zip(bc, it, it_new, ar, ac, Y, X): f.write(f"{b},{i0},{i1},{a},{c},{int(round(y))},{int(round(x))}\n")
json.dump(dict(reflection="array_row -> 77 - array_row (array_col unchanged), positions from the lattice fit",
               lattice_fit=dict(pxl_col=cx.tolist(), pxl_row=cy.tolist()), mask_threshold_below_background=12,
               mask_agreement_with_spaceranger_at_delivered_positions=float(agree), n_in_tissue_delivered=int(it.sum()),
               n_in_tissue_registered=int(it_new.sum()), newly_in=int(((it_new==1)&(it==0)).sum()), newly_out=int(((it_new==0)&(it==1)).sum())),
          open(os.path.join(DST, "REGISTRATION_provenance.json"), "w"), indent=1)
print("wrote", DST, "\n=== DONE ===")
