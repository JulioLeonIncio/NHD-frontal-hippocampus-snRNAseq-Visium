#!/usr/bin/env python3
# 108a_visium_HE_crop_FH.py — Crop each full-resolution Visium H&E to the spot extent used by Fig. 6b and put it into Fig. 6b's plotting frame.
# Fig. 6b (script 74) plots x = image row and y = image column of each spot at low-res scale
# (spot_coords.csv from cell2location = tissue_positions pxl_row/pxl_col x tissue_lowres_scalef),
# so the maps are the image transposed. The mapping pixel -> map coordinate is fitted per section
# from the spots themselves (a linear fit of x on pxl_row and y on pxl_col; a negative slope would
# be a mirror), so the script is agnostic to how any section's outs were produced. With the
# barcode-map registration step applied to NHD_Frontal2 (108d/108e), the fit reports mirror=False
# for all four sections — recorded per section in HE_frames.json.
# The crop is then transposed/flipped to sit in the same frame as the domain map.
#
# Panel SOURCE = the Space Ranger hires image of the same outs that produced the spot coordinates.
# Which outs folder each section is read from comes from Visium/VISIUM_OUTS_MANIFEST.json (
# <sec>_masked/outs with the three-tier tissue call; NHD_Frontal2's inherits the reflected barcode map of
# NHD_Frontal2_registered, 108d/108e). 1,968 px across the capture area is 3.5x what a 0.95-in panel needs
# at 600 dpi. All four hires images are in the native scan orientation, so the full-resolution tiff
# is in the same frame and the annotation copy is cropped with the same box for every section.
# Registration is verified data-side by 108c (spots on the tear are UMI-poor; WM domain on the pale
# WM) for every section, NHD_Frontal2 included.
#
# Outputs (Visium/_HE_fullres/panels/):
#   <section>_HE_fig6b_frame.png   the hires crop in Fig-6b orientation
#   <section>_HE_for_annotation.png the same crop from the full-resolution tiff, native orientation,
#                                   4,000 px long side (for Takao-sensei / Satoshi; lossless PNG)
#   HE_frames.json                  per section: map-coordinate extent of the crop, µm per map unit,
#                                   scale-bar length in map units, fit residuals
# Python environment: ctm_env (python_env_ctm_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import json, os, sys
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
ROOT = os.path.dirname(os.environ["NHD_PROJ"])
HE   = os.path.join(ROOT, "Visium", "_HE_fullres")
OUT  = os.path.join(HE, "panels"); os.makedirs(OUT, exist_ok=True)
COORDS = os.path.join(ROOT, "Visium", "cell2location", "c2l_MAIN", "spot_coords.csv")
SECTIONS = ["CON_Frontal1", "CON_Frontal2", "NHD_Frontal1", "NHD_Frontal2"]
# Outs folder per section from the manifest (same source as R's _visium_outs_FH.R); never hard-code it here
MANIFEST = os.path.join(ROOT, "Visium", "VISIUM_OUTS_MANIFEST.json")
assert os.path.exists(MANIFEST), f"MISSING {MANIFEST} — run scripts/117b_visium_tissue_call_outs_FH.py"
OUTS_MAP = json.load(open(MANIFEST))["outs"]
assert all(s in OUTS_MAP for s in SECTIONS), "VISIUM_OUTS_MANIFEST.json lacks a section"
def vis_outs(sec):
    d = os.path.join(ROOT, "Visium", OUTS_MAP[sec], "outs")
    assert os.path.isdir(d), f"outs folder missing for {sec}: {d}"
    return d
SPOT_UM = 55.0          # Visium spot diameter

coords = {}
with open(COORDS) as f:
    next(f)
    for line in f:
        cell, x, y = line.rstrip("\n").replace('"', "").split(",")
        sec, bc = cell.rsplit("_", 1)
        coords.setdefault(sec, {})[bc] = (float(x), float(y))

frames = {}
for sec in SECTIONS:
    outs = os.path.join(vis_outs(sec), "spatial")   # manifest outs
    pos = {}
    with open(os.path.join(outs, "tissue_positions_list.csv")) as f:
        for line in f:
            bc, it, ar, ac, pr, pc = line.rstrip("\n").split(",")
            pos[bc] = (int(pr), int(pc))
    sf = json.load(open(os.path.join(outs, "scalefactors_json.json")))
    px_per_um = sf["spot_diameter_fullres"] / SPOT_UM
    hs = sf["tissue_hires_scalef"]
    # fit map x ~ pxl_row, map y ~ pxl_col on the spots Fig 6b plots
    bcs = [b for b in coords[sec] if b in pos]
    X = np.array([coords[sec][b] for b in bcs]); P = np.array([pos[b] for b in bcs], float)
    ax, bx = np.polyfit(P[:, 0], X[:, 0], 1); ay, by = np.polyfit(P[:, 1], X[:, 1], 1)
    resid = max(np.abs(ax * P[:, 0] + bx - X[:, 0]).max(), np.abs(ay * P[:, 1] + by - X[:, 1]).max())
    mirror = ay < 0
    # crop box in full-res pixels = spot extent + half a spot spacing margin
    pad = 1.5 * sf["spot_diameter_fullres"]
    r0, r1 = int(max(0, P[:, 0].min() - pad)), int(P[:, 0].max() + pad)
    c0, c1 = int(max(0, P[:, 1].min() - pad)), int(P[:, 1].max() + pad)
    # annotation copy from the full-resolution TIFF (native orientation = the hires frame for all four)
    tif = os.path.join(HE, f"{sec}_HE_fullres.tif")
    if os.path.exists(tif):
        im = Image.open(tif).convert("RGB")
        a = im.crop((c0, r0, c1 + 1, r1 + 1))
        a.thumbnail((4000, 4000), Image.LANCZOS); a.save(os.path.join(OUT, f"{sec}_HE_for_annotation.png"), optimize=True); im.close()
    # panel crop from the hires image of the same outs (coordinates scaled by tissue_hires_scalef)
    hi = Image.open(os.path.join(outs, "tissue_hires_image.png")).convert("RGB")
    crop = hi.crop((int(c0 * hs), int(r0 * hs), int(c1 * hs) + 1, int(r1 * hs) + 1))
    # Fig-6b frame: horizontal axis = row (increasing to the right), vertical axis = col with
    # map y increasing up. Transpose (rows -> columns), then the sign of ay decides the flip.
    arr = np.asarray(crop)                                # [row, col, 3]
    fr = np.transpose(arr, (1, 0, 2))                    # [col, row, 3]: image row -> horizontal
    # after transpose, vertical index = col; map y increases with col if ay > 0 -> image top must be
    # the high col end (because image row 0 is at the top): flip vertically when ay > 0
    if ay > 0: fr = fr[::-1, :, :]
    if ax < 0: fr = fr[:, ::-1, :]
    fim = Image.fromarray(np.ascontiguousarray(fr))
    fim.save(os.path.join(OUT, f"{sec}_HE_fig6b_frame.png"), optimize=True)
    # extent of the crop in map units
    xs = sorted([ax * r0 + bx, ax * r1 + bx]); ys = sorted([ay * c0 + by, ay * c1 + by])
    um_per_map = 1.0 / (abs(ax) * px_per_um)            # map unit = |ax| px... -> µm per map unit
    frames[sec] = dict(x_min=xs[0], x_max=xs[1], y_min=ys[0], y_max=ys[1], mirror=bool(mirror),
                       um_per_map_unit=um_per_map, scalebar_1mm_map_units=1000.0 / um_per_map,
                       fit_max_residual_map_units=float(resid), crop_px=[r0, r1, c0, c1], px_per_um=px_per_um)
    print(f"{sec}: crop rows {r0}-{r1} cols {c0}-{c1}; mirror={mirror}; fit residual {resid:.3f} map units; 1 mm = {1000/um_per_map:.1f} map units")
json.dump(frames, open(os.path.join(OUT, "HE_frames.json"), "w"), indent=1)
print("=== DONE ===")
