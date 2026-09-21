#!/usr/bin/env python3
# 108c_visium_HE_registration_check_FH.py — Visium barcode-map registration: check. Is each section's barcode->pixel map registered to its H&E? Data-driven, independent of Space Ranger's
# tissue detection (Task C).
#
# For every section and each of the 8 dihedral transforms of the spot positions (as delivered,
# mirrors, rotations), two statistics:
#   (1) spots whose 55-um footprint is mostly BACKGROUND (white) — their median UMI vs the median on
#       tissue, and Spearman(UMI, white fraction). Under the true registration, only barcodes over a
#       tear / hole can sit on background and they must be UMI-poor. (Space Ranger's in_tissue call
#       is image-based, so 'as delivered' has ~0 background spots by construction; the test is whether a
#       candidate transform finds a UMI-poor set on background.)
#   (2) saturation of WM-domain spots vs cortical-domain spots (WM is paler/greyer in these H&Es;
#       not decisive alone — CON_Frontal1's spongiotic superficial cortex is also grey).
# Also writes an overlay per section (as delivered | mirrored) with the lowest-UMI decile in black.
#
# Validation. CON_Frontal1, CON_Frontal2 and NHD_Frontal1 read as registered from
# their Space Ranger outs (5 / 0 / 9 background spots). For NHD_Frontal2 the Space Ranger --reorient-images
# run places the barcode layout in mirror image relative to the scan (a known failure mode of fiducial
# alignment on mirrored images): on that run's map the horizontal-mirror hypothesis puts 238 spots on the
# tear with median 1,084 UMI vs 4,779 on tissue (rho = -0.43), and WM lands on the pale wedge and
# Vasc/immune on the vessel. The pipeline therefore read NHD_Frontal2 from NHD_Frontal2_registered/outs,
# where 108d/108e reflect the barcode map about the array centre and derive the tissue mask from the H&E;
# on those outs NHD_Frontal2 reads as registered like the other three sections.
#
# Every section is now read from the outs named in
# Visium/VISIUM_OUTS_MANIFEST.json (<sec>_masked/outs; NHD_Frontal2's inherits the reflected barcode map).
# The analysis tissue call has three tiers (tissue_call.csv, column tissue_call_source): spaceranger
# (Space Ranger's image-based call), he_mask (relaxed H&E mask) and expression (an expression-supported
# contiguous rim of barcodes at the tissue edge that the H&E reads as pale/white). The expression tier is
# therefore on background by construction under the white-fraction test, so:
#   * the registered threshold (<= max(10, 1 %) background spots as delivered) is applied to the
#     spaceranger + he_mask spots only (n_on_background_core);
#   * the 8-transform minimum test (as delivered must put no more spots on background than any other
#     dihedral transform) is kept over all spots — mirror detection is the point of the check;
#   * the expression tier is reported separately per transform: n_expression_tier_on_white (expected to be
#     about the tier's size as delivered) with its median UMI (median_UMI_expression_tier_on_white).
# If an outs folder has no tissue_call.csv, every spot is treated as tier 'spaceranger' (a printed warning).
# This script prints a registered / not registered verdict per section and exits non-zero if any
# section fails, so the downstream rebuild driver (112) uses it as its final gate. Run again after any
# change to the section outs (spatial/ or tissue_call.csv) to confirm every section reads 'registered'.
# Python environment: ctm_env (python_env_ctm_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import json, os, sys, numpy as np
from PIL import Image, ImageDraw, ImageOps
if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
ROOT = os.path.dirname(os.environ["NHD_PROJ"])
PROJ = os.path.join(ROOT, "NHD_frontal_hippo_rebuild")
C2L  = os.path.join(ROOT, "Visium", "cell2location", "c2l_MAIN")
DIAG = os.path.join(PROJ, "figures", "_diagnostics"); os.makedirs(DIAG, exist_ok=True)
# cluster -> domain map is read from integrated_cluster_identity.csv (the curated identity file of the
# current integration); never hard-code it here.
DOM = {}
with open(os.path.join(ROOT, "Visium", "integrated_harmony", "integrated_cluster_identity.csv")) as f:
    hdr = next(f).rstrip("\n").replace('"', "").split(","); ic_, ii_ = hdr.index("cluster"), hdr.index("identity")
    for line in f:
        p = line.rstrip("\n").replace('"', "").split(","); DOM[p[ic_].lstrip("g")] = p[ii_]
assert DOM, "integrated_cluster_identity.csv is empty"
# Outs folder per section from the manifest (same source as R's _visium_outs_FH.R); never hard-code it here
MANIFEST = os.path.join(ROOT, "Visium", "VISIUM_OUTS_MANIFEST.json")
assert os.path.exists(MANIFEST), f"MISSING {MANIFEST} — run scripts/117b_visium_tissue_call_outs_FH.py"
OUTS = json.load(open(MANIFEST))["outs"]          # section -> folder under Visium/ (its outs/ is read)
SECTIONS = ["CON_Frontal1", "CON_Frontal2", "NHD_Frontal1", "NHD_Frontal2"]
assert all(s in OUTS for s in SECTIONS), "VISIUM_OUTS_MANIFEST.json lacks a section"
OUTS = {s: OUTS[s] for s in SECTIONS}             # fixed section order
CORE_TIERS = ("spaceranger", "he_mask")           # tiers the registered threshold is judged on
def read_tissue_call(outs_dir, sec):
    """barcode -> tissue_call_source from tissue_call.csv (117b); all-'spaceranger' if the file is absent."""
    f = os.path.join(outs_dir, "tissue_call.csv")
    if not os.path.exists(f):
        print(f"  WARNING {sec}: no tissue_call.csv in {outs_dir} — every spot treated as tier 'spaceranger'"); return {}
    src = {}
    with open(f) as fh:
        hdr = next(fh).rstrip("\n").replace('"', "").split(","); ib, isrc = hdr.index("barcode"), hdr.index("tissue_call_source")
        for line in fh:
            q = line.rstrip("\n").replace('"', "").split(","); src[q[ib]] = q[isrc]
    assert src, f"{f} is empty"
    return src
meta = {}
with open(os.path.join(C2L, "visium_meta.csv")) as f:
    hdr = next(f).rstrip("\n").replace('"', "").split(","); ic = hdr.index("seurat_clusters"); iu = hdr.index("nCount_Spatial")
    for line in f:
        p = line.rstrip("\n").replace('"', "").split(",")
        if p[ic] not in DOM: sys.exit(f"cluster {p[ic]} in visium_meta.csv has no identity — integrated_cluster_identity.csv is stale vs the object")
        meta[p[0]] = (DOM[p[ic]], float(p[iu]))
def sat_feats(img, X, Y, r):
    H, W, _ = img.shape; out = np.full(len(X), np.nan)
    for i, (x, y) in enumerate(zip(X, Y)):
        x0, x1, y0, y1 = int(max(0, x - r)), int(min(W, x + r + 1)), int(max(0, y - r)), int(min(H, y + r + 1))
        if x1 <= x0 or y1 <= y0: continue
        w = img[y0:y1, x0:x1].reshape(-1, 3).astype(float); mx, mn = w.max(1), w.min(1); out[i] = ((mx - mn) / np.maximum(mx, 1)).mean()
    return out
rows = []
for sec, var in OUTS.items():
    outs_dir = os.path.join(ROOT, "Visium", var, "outs"); assert os.path.isdir(outs_dir), f"outs folder missing for {sec}: {outs_dir}"
    outs = os.path.join(outs_dir, "spatial")
    tier = read_tissue_call(outs_dir, sec)   # per-barcode tissue-call tier
    sf = json.load(open(os.path.join(outs, "scalefactors_json.json"))); hs = sf["tissue_hires_scalef"]
    im = Image.open(os.path.join(outs, "tissue_hires_image.png")).convert("RGB"); img = np.asarray(im); H, W, _ = img.shape
    grey = np.asarray(im.convert("L")).astype(float); white = grey > np.quantile(grey, 0.995) - 12
    P, D, U, SRC = [], [], [], []
    for line in open(os.path.join(outs, "tissue_positions_list.csv")):
        bc, it, ar, ac, pr, pc = line.rstrip("\n").split(","); k = f"{sec}_{bc}"
        if k in meta: P.append((float(pc) * hs, float(pr) * hs)); D.append(meta[k][0]); U.append(meta[k][1]); SRC.append(tier.get(bc, "spaceranger"))
    P = np.array(P); D = np.array(D); U = np.array(U); SRC = np.array(SRC); r = sf["spot_diameter_fullres"] * hs / 2
    core = np.isin(SRC, CORE_TIERS); expr = SRC == "expression"   # tiers judged / reported separately
    if tier and (SRC == "excluded").any(): sys.exit(f"{sec}: {(SRC == 'excluded').sum()} object spots are tier 'excluded' in tissue_call.csv — object and outs disagree")
    x = P[:, 0] - (W - 1) / 2; y = P[:, 1] - (H - 1) / 2
    T = {"as_delivered": (x, y), "mirrorH": (-x, y), "flipV": (x, -y), "rot180": (-x, -y), "transpose": (y, x), "rot90": (-y, x), "rot270": (y, -x), "antitranspose": (-y, -x)}
    wm = D == "WM"; ctx = np.isin(D, ["L2/3", "L4", "L5", "L6"])
    print(f"== {sec} ({var})  n = {len(U)}  median UMI = {np.median(U):.0f}  tiers: " + " ".join(f"{t}={int((SRC == t).sum())}" for t in ("spaceranger", "he_mask", "expression")))
    for name, (tx, ty) in T.items():
        X = tx + (W - 1) / 2; Y = ty + (H - 1) / 2; ri = int(r)
        wf = np.array([white[int(max(0, yy - ri)):int(yy + ri + 1), int(max(0, xx - ri)):int(xx + ri + 1)].mean() if 0 <= xx < W and 0 <= yy < H else np.nan for xx, yy in zip(X, Y)])
        ok = ~np.isnan(wf); onbg = ok & (wf > 0.5)
        rho = np.corrcoef(np.argsort(np.argsort(U[ok])), np.argsort(np.argsort(wf[ok])))[0, 1]
        S = sat_feats(img, X, Y, r); a = S[wm & ~np.isnan(S)]; b = S[ctx & ~np.isnan(S)]; auc = (a[:, None] < b[None, :]).mean()
        med_bg = np.median(U[onbg]) if onbg.sum() else np.nan; med_t = np.median(U[ok & ~onbg])
        onbg_core = onbg & core; expr_white = onbg & expr   # core tiers carry the threshold; expression tier reported apart
        med_ew = np.median(U[expr_white]) if expr_white.sum() else np.nan
        rows.append(dict(section=sec, outs=var, transform=name, n_on_background=int(onbg.sum()), n_on_background_core=int(onbg_core.sum()), n_core=int(core.sum()),
                         n_expression_tier=int(expr.sum()), n_expression_tier_on_white=int(expr_white.sum()), median_UMI_expression_tier_on_white=med_ew,
                         median_UMI_on_background=med_bg, median_UMI_on_tissue=med_t, spearman_UMI_vs_whitefrac=rho, AUC_WM_greyer_than_cortex=auc))
        print(f"  {name:14s} on-background n={onbg.sum():4d} (core {onbg_core.sum():4d}; expression tier {expr_white.sum():4d}/{expr.sum()} on white, median UMI {med_ew:6.0f})  "
              f"UMI bg/tissue = {med_bg:7.0f} / {med_t:6.0f}  rho(UMI,white) = {rho:+.3f}  AUC(WM greyer) = {auc:.3f}")
    # overlay: as delivered | mirrorH, lowest-UMI decile black, WM domain purple
    r2 = r; panels = []
    for mode in ["as_delivered", "mirrorH"]:
        ov = ImageOps.autocontrast(im.copy(), cutoff=1); d = ImageDraw.Draw(ov); q = np.quantile(U, 0.1)
        for (px, py), dom, u in zip(P, D, U):
            xx = px if mode == "as_delivered" else (W - 1 - px)
            if dom == "WM": d.ellipse([xx - r2, py - r2, xx + r2, py + r2], fill=(150, 130, 220))
            if u < q: d.ellipse([xx - r2 * 0.6, py - r2 * 0.6, xx + r2 * 0.6, py + r2 * 0.6], fill=(0, 0, 0))
        panels.append(ov)
    out = Image.new("RGB", (W * 2 + 20, H), (255, 255, 255)); out.paste(panels[0], (0, 0)); out.paste(panels[1], (W + 20, 0))
    out.resize((out.width // 2, out.height // 2)).save(os.path.join(DIAG, f"HE_registration_{sec}_asdelivered_vs_mirror.png"))
import csv
with open(os.path.join(PROJ, "tables", "visium_HE_registration_check_FH.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
# ---- verdict per section (the rebuild driver 112 gates on it) ------------------------------------------
# A section reads registered when (i) the as-delivered barcode map puts no more in-tissue barcodes on
# background than any of the 7 other dihedral transforms (all spots — mirror detection) AND (ii) the
# as-delivered background count over the spaceranger + he_mask tiers is <= max(10, 1 %) of those spots.
# (ii) was over all spots before the three-tier tissue call; the expression tier sits on
# H&E-white pixels by construction and is reported apart (n_expression_tier_on_white), not gated on.
# Validation: this criterion separates the three sections registered by Space Ranger
# (5 / 0 / 9 background spots) from the NHD_Frontal2 --reorient-images map (243 background spots, median
# 1,085 vs 4,783 UMI), and passes the reflected NHD_Frontal2_registered map. The UMI-poor background set
# and the WM-saturation AUC stay in the table as supporting evidence but are not gated on: a rotation can
# land genuine low-UMI edge spots on background in a registered section (CON_Frontal1 rot180).
verdict_ok = True
print("\n== registration verdict ==")
for sec in OUTS:
    rs = [r for r in rows if r["section"] == sec]; n = sum(1 for k in meta if k.startswith(sec + "_"))
    asd = next(r for r in rs if r["transform"] == "as_delivered"); best = min(rs, key=lambda r: r["n_on_background"])
    n_core = asd["n_core"]
    ok_thr = asd["n_on_background_core"] <= max(10, 0.01 * n_core)      # (ii) core tiers only
    ok_min = asd["n_on_background"] <= best["n_on_background"]          # (i) all spots, 8 transforms
    ok = ok_thr and ok_min
    verdict_ok &= ok
    print(f"  {sec:13s} ({asd['outs']:22s}) {'REGISTERED' if ok else 'NOT REGISTERED'}: as-delivered background spots = "
          f"{asd['n_on_background']} of {n} (core tiers {asd['n_on_background_core']} of {n_core}; expression tier "
          f"{asd['n_expression_tier_on_white']} of {asd['n_expression_tier']} on white, median UMI {asd['median_UMI_expression_tier_on_white']:.0f}); "
          f"fewest under '{best['transform']}' = {best['n_on_background']}")
print("wrote tables/visium_HE_registration_check_FH.csv + figures/_diagnostics/HE_registration_*.png")
if not verdict_ok:
    sys.exit("=== FAILED: at least one section is NOT registered — do not ship Fig 6b ===")
print("=== DONE ===")
