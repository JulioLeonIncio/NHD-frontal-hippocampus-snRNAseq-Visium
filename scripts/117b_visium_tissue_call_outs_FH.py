#!/usr/bin/env python3
# 117b_visium_tissue_call_outs_FH.py — Evidence-based tissue call for all four Visium sections.
# Why. Space Ranger's automatic tissue detection called the pale pial rim of the control sections
# out-of-tissue (CON_Frontal1 379, CON_Frontal2 1,002 barcodes) although those barcodes carry median
# 4,334 / 3,297 UMI — 10-20x the UMI of the control white matter that was kept (median ~200). The rim is
# layer I / glia limitans. NHD_Frontal1 (197) and
# NHD_Frontal2 (394, registered call) have smaller excluded rims of the same kind. An H&E grey threshold
# cannot recover them (the NHD_Frontal1 rim is as pale as glass at full resolution yet carries 2,989 UMI),
# so the call is made in three TIERS, recorded per barcode:
#   T1 spaceranger : the delivered call (NHD_Frontal2: the 108d registered call, tear already removed);
#   T2 he_mask     : >= 50 % of the 55-um footprint on a relaxed H&E mask (grey < background - 8 on the
#                    hires image, closing 2 / opening 1 as in 108d); only barcodes outside TEARS are eligible
#                    (tear = interior hole of the mask >= HOLE_MIN_PX on the hires image; the NHD_Frontal2 tear);
#   T3 expression  : nCount_raw >= EXPR_UMI and nFeature_raw >= EXPR_GENES (the cell2location spot floor of
#                    script 65), outside interior holes, and contiguous with T1 u T2 by hex-neighbour flood
#                    fill on the array grid — RNA can only come from tissue that was on the slide.
# Positions are never changed (NHD_Frontal2 inherits the 108d reflection).
#
# What IT WRITES (staged in Visium/_staging/, then renamed into place; an existing target is renamed
# _superseded_<ts>, never deleted): Visium/<section>_masked/outs/
#   spatial/ (images + scalefactors copied; tissue_positions_list.csv with the analysis call),
#   raw_feature_bc_matrix.h5 (verbatim copy), filtered_feature_bc_matrix.h5 + filtered_feature_bc_matrix/
#   (rebuilt from raw with the analysis call, Space Ranger h5 layout as in 108e), metrics_summary.csv,
#   tissue_call.csv (barcode, in_tissue_spaceranger, in_tissue_analysis, tissue_call_source, he_fraction,
#   nCount_raw, nFeature_raw, pct_mt_raw, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres),
#   TISSUE_CALL_provenance.json (thresholds, tier counts, mask agreement; "approved_by" left empty for the pi).
# Plus Visium/VISIUM_OUTS_MANIFEST.json (section -> outs folder) read by every downstream script, the
# review figure figures/_diagnostics/visium_tissue_call/<section>_tissue_call_tiers.png and
# tables/visium_tissue_call_summary_FH.csv.
# Run: scripts/117b_visium_tissue_call_outs_FH.py   [env TIERS=spaceranger,he_mask,expression]
# Python environment: ctm_env (python_env_ctm_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import os, sys, json, csv, gzip, io, shutil, tempfile, time, numpy as np, h5py, scipy.sparse as sps, scipy.io
from PIL import Image
from scipy.ndimage import uniform_filter, binary_closing, binary_opening, binary_fill_holes, label
DARK_GREY = 70; DARK_MIN_PX = 1000; DARK_FOOT = 0.05   # exogenous dark particulate contaminant (near-black object >= 1000 hires px,
                                                         # nothing in a section stains that dark): barcodes with > 5 % of their footprint on it are excluded
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
ROOT = os.path.dirname(os.environ["NHD_PROJ"])
VIS  = os.path.join(ROOT, "Visium"); PROJ = os.path.join(ROOT, "NHD_frontal_hippo_rebuild")
FIGD = os.path.join(PROJ, "figures", "_diagnostics", "visium_tissue_call"); os.makedirs(FIGD, exist_ok=True)
SRC  = {"CON_Frontal1": "CON_Frontal1", "CON_Frontal2": "CON_Frontal2", "NHD_Frontal1": "NHD_Frontal1", "NHD_Frontal2": "NHD_Frontal2_registered"}
TIERS = os.environ.get("TIERS", "spaceranger,he_mask,expression").split(",")
HE_DELTA = 8; EXPR_UMI = 1000; EXPR_GENES = 500; FOOT = 0.5
NEUROPIL = ["SNAP25", "RBFOX3", "SLC17A7"]; NEUROPIL_MIN = 12.0   # recovered rim spots must contain cortical neuropil
# (>= 12 per 10,000; layer I reads 16-40, the outermost ring over detached arachnoid / subarachnoid blood reads < 10). Applied to the
# H&E-mask tier only: control white-matter libraries (~200 UMI) are too small for a transcript floor and stay under the Space Ranger call.
HOLE_MIN_PX = 10000   # an interior hole of the mask counts as a TEAR (never eligible) only above this size on the hires image
                      # (NHD_Frontal2's tear = two components of 14.6k / 13.9k px; every other section's largest hole is <= 7.4k px:
                      # small holes are pale patches inside tissue and stay eligible)
TS = time.strftime("%Y%m%d_%H%M%S")
STAGE = os.path.join(VIS, "_staging"); os.makedirs(STAGE, exist_ok=True)

def read_positions(p):
    rows = [l.rstrip("\n").split(",") for l in open(p)]
    rows = [r for r in rows if r[0] != "barcode"]
    return {r[0]: (int(r[1]), int(r[2]), int(r[3]), float(r[4]), float(r[5])) for r in rows}

def read_raw(h5):
    with h5py.File(h5, "r") as f:
        m = f["matrix"]; bcs = m["barcodes"][:].astype(str); shape = m["shape"][:]
        M = sps.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=(shape[0], shape[1]))
        feats = {k: m["features"][k][:] for k in m["features"].keys()}; attrs = dict(f.attrs)
    names = feats["name"].astype(str)
    umi = np.asarray(M.sum(0)).ravel(); ng = np.asarray((M > 0).sum(0)).ravel()
    mt = np.asarray(M[np.char.startswith(names, "MT-"), :].sum(0)).ravel()
    neu = np.asarray(M[np.isin(names, NEUROPIL), :].sum(0)).ravel()            # cortical neuropil transcripts
    return bcs, M, feats, attrs, umi, ng, np.where(umi > 0, 100 * mt / np.maximum(umi, 1), 0.0), 1e4 * neu / np.maximum(umi, 1)

def hex_neighbours(r, c):
    return [(r, c - 2), (r, c + 2), (r - 1, c - 1), (r - 1, c + 1), (r + 1, c - 1), (r + 1, c + 1)]

def write_filtered(dst, bcs, M, feats, attrs, keep_idx, libname):
    F = M[:, keep_idx].tocsc(); F.sort_indices()
    out = os.path.join(dst, "filtered_feature_bc_matrix.h5")
    with h5py.File(out, "w") as g:
        for k, v in attrs.items():
            if k == "library_ids": v = np.array([libname.encode()], dtype=f"|S{len(libname)}")
            g.attrs[k] = v
        mm = g.create_group("matrix")
        mm.create_dataset("barcodes", data=bcs[keep_idx].astype("S18"))
        mm.create_dataset("data", data=F.data.astype(np.int32), compression="gzip")
        mm.create_dataset("indices", data=F.indices.astype(np.int64), compression="gzip")
        mm.create_dataset("indptr", data=F.indptr.astype(np.int64))
        mm.create_dataset("shape", data=np.array(F.shape, dtype=np.int32))
        ff = mm.create_group("features")
        for k, v in feats.items(): ff.create_dataset(k, data=v)
    d = os.path.join(dst, "filtered_feature_bc_matrix"); os.makedirs(d, exist_ok=True)
    buf = io.BytesIO(); scipy.io.mmwrite(buf, F.astype(np.int32), field="integer")
    with gzip.open(os.path.join(d, "matrix.mtx.gz"), "wb") as z: z.write(buf.getvalue())
    with gzip.open(os.path.join(d, "barcodes.tsv.gz"), "wt") as z: z.write("\n".join(bcs[keep_idx]) + "\n")
    with gzip.open(os.path.join(d, "features.tsv.gz"), "wt") as z:
        for i in range(F.shape[0]): z.write(f"{feats['id'][i].decode()}\t{feats['name'][i].decode()}\t{feats['feature_type'][i].decode()}\n")
    return F.shape

summary = []; manifest = {"written_by": "117b_visium_tissue_call_outs_FH.py", "date": TS, "tiers": TIERS,
                          "thresholds": {"he_delta": HE_DELTA, "footprint_fraction": FOOT, "expr_umi": EXPR_UMI, "expr_genes": EXPR_GENES}, "outs": {}}
for sec, srcname in SRC.items():
    src = os.path.join(VIS, srcname, "outs"); sp = os.path.join(src, "spatial")
    sf = json.load(open(os.path.join(sp, "scalefactors_json.json"))); hs = sf["tissue_hires_scalef"]
    pos = read_positions(os.path.join(sp, "tissue_positions_list.csv"))
    bcs, M, feats, attrs, umi, ng, pmt, neu10k = read_raw(os.path.join(src, "raw_feature_bc_matrix.h5"))
    inpos = np.array([b in pos for b in bcs]); assert inpos.all(), f"{sec}: barcodes missing from tissue_positions"
    it_sr = np.array([pos[b][0] for b in bcs]); ar = np.array([pos[b][1] for b in bcs]); ac = np.array([pos[b][2] for b in bcs])
    pr = np.array([pos[b][3] for b in bcs]); pc = np.array([pos[b][4] for b in bcs])
    # --- T2: relaxed H&E mask on the hires image (same frame as pxl * hires scalef) ---
    img = np.asarray(Image.open(os.path.join(sp, "tissue_hires_image.png")).convert("RGB")).astype(float)
    grey = img.mean(2); bg = np.quantile(grey, 0.995)
    mask = binary_opening(binary_closing(grey < bg - HE_DELTA, iterations=2), iterations=1)
    holes_all = binary_fill_holes(mask) & ~mask
    lab_h, n_h = label(holes_all); sz = np.bincount(lab_h.ravel()); sz[0] = 0
    holes = sz[lab_h] >= HOLE_MIN_PX
    dk_lab, _ = label(grey < DARK_GREY); dk_sz = np.bincount(dk_lab.ravel()); dk_sz[0] = 0; dark = dk_sz[dk_lab] >= DARK_MIN_PX
    r = max(1, int(round(sf["spot_diameter_fullres"] * hs / 2)))
    mf = uniform_filter(mask.astype(float), size=2 * r + 1); dkf = uniform_filter(dark.astype(float), size=2 * r + 1)
    xs = np.clip((pc * hs).astype(int), 0, mf.shape[1] - 1); ys = np.clip((pr * hs).astype(int), 0, mf.shape[0] - 1)
    hef = mf[ys, xs]; in_hole = holes[ys, xs]; on_dark = dkf[ys, xs] > DARK_FOOT
    mask108 = binary_opening(binary_closing(grey < bg - 12, iterations=2), iterations=1)
    agree108 = ((uniform_filter(mask108.astype(float), size=2 * r + 1)[ys, xs] > FOOT) == (it_sr == 1)).mean()
    t1 = (it_sr == 1) & ~on_dark
    t2 = (hef > FOOT) & ~in_hole & ~t1 & ~on_dark & (neu10k >= NEUROPIL_MIN)
    n_dark = int(((it_sr == 1) & on_dark).sum())   # Space Ranger in-tissue barcodes excluded as contaminant
    n_menin = int(((hef > FOOT) & ~in_hole & ~t1 & (neu10k < NEUROPIL_MIN)).sum())   # H&E-mask spots rejected as meninges / blood
    # --- T3: expression tier, contiguous with T1 u T2 by hex flood fill ---
    elig = (umi >= EXPR_UMI) & (ng >= EXPR_GENES) & ~in_hole & ~t1 & ~t2
    grid = {(int(a), int(c)): i for i, (a, c) in enumerate(zip(ar, ac))}
    seed = t1 | (t2 if "he_mask" in TIERS else False)
    reach = seed.copy(); frontier = list(np.where(seed)[0])
    while frontier:
        i = frontier.pop()
        for nb in hex_neighbours(int(ar[i]), int(ac[i])):
            j = grid.get(nb)
            if j is not None and not reach[j] and elig[j]:
                reach[j] = True; frontier.append(j)
    t3 = reach & elig
    source = np.full(len(bcs), "excluded", dtype=object)
    source[t1] = "spaceranger"
    if "he_mask" in TIERS: source[t2] = "he_mask"
    if "expression" in TIERS: source[t3] = "expression"
    source[on_dark & (it_sr == 1)] = "contaminant"
    it_new = ((source != "excluded") & (source != "contaminant")).astype(int)   # contaminant spots are OUT of the analysis
    n = {k: int((source == k).sum()) for k in ["spaceranger", "he_mask", "expression", "excluded", "contaminant"]}
    med = lambda m: int(np.median(umi[m])) if m.any() else 0
    print(f"{sec}: contaminant-excluded {n_dark} | spaceranger {n['spaceranger']} | he_mask +{n['he_mask']} (median UMI {med(source=='he_mask')}; {n_menin} H&E-mask spots rejected as meninges, neuropil < {NEUROPIL_MIN}/10k) | expression +{n['expression']} "
          f"(median UMI {med(source=='expression')}) | excluded {n['excluded']} (median UMI {med(source=='excluded')}) | "
          f"108d-rule agreement with SR {agree108:.3f} | holes {int(in_hole.sum())} barcodes")
    # --- stage the outs ---
    dst_final = os.path.join(VIS, f"{sec}_masked")
    stage = tempfile.mkdtemp(prefix=f"{sec}_masked_", dir=STAGE); outs = os.path.join(stage, "outs"); os.makedirs(os.path.join(outs, "spatial"))
    for f in ["tissue_hires_image.png", "tissue_lowres_image.png", "detected_tissue_image.jpg", "aligned_fiducials.jpg", "scalefactors_json.json"]:
        if os.path.exists(os.path.join(sp, f)): shutil.copy2(os.path.join(sp, f), os.path.join(outs, "spatial", f))
    with open(os.path.join(outs, "spatial", "tissue_positions_list.csv"), "w") as f:
        for b, i, a, c, y, x in zip(bcs, it_new, ar, ac, pr, pc): f.write(f"{b},{i},{a},{c},{int(round(y))},{int(round(x))}\n")
    for f in ["raw_feature_bc_matrix.h5", "metrics_summary.csv", "in_tissue_registered.csv", "REGISTRATION_provenance.json"]:
        if os.path.exists(os.path.join(src, f)): shutil.copy2(os.path.join(src, f), os.path.join(outs, f))
    shp = write_filtered(outs, bcs, M, feats, attrs, np.where(it_new == 1)[0], f"{sec}_masked")
    with open(os.path.join(outs, "tissue_call.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["barcode", "in_tissue_spaceranger", "in_tissue_analysis", "tissue_call_source", "he_fraction", "nCount_raw", "nFeature_raw", "pct_mt_raw", "neuropil_per10k", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
        for i, b in enumerate(bcs): w.writerow([b, int(it_sr[i]), int(it_new[i]), source[i], round(float(hef[i]), 3), int(umi[i]), int(ng[i]), round(float(pmt[i]), 2), round(float(neu10k[i]), 2), int(ar[i]), int(ac[i]), int(round(pr[i])), int(round(pc[i]))])
    prov = {"section": sec, "source_outs": src, "written": TS, "tiers": TIERS, "he_delta": HE_DELTA, "footprint_fraction": FOOT, "expr_umi": EXPR_UMI, "expr_genes": EXPR_GENES,
            "he_mask_neuropil_min_per10k": NEUROPIL_MIN, "he_mask_rejected_as_meninges": n_menin,
            "contaminant_rule": {"dark_grey_max": DARK_GREY, "min_component_px": DARK_MIN_PX, "footprint_fraction": DARK_FOOT}, "spaceranger_excluded_as_contaminant": n_dark,
            "n_by_source": n, "n_in_tissue_analysis": int(it_new.sum()), "n_in_tissue_spaceranger": int(t1.sum()), "median_umi_by_source": {k: med(source == k) for k in n},
            "mask108_agreement_with_spaceranger": round(float(agree108), 4), "n_barcodes_in_interior_holes": int(in_hole.sum()), "hole_min_px": HOLE_MIN_PX, "filtered_shape": [int(shp[0]), int(shp[1])],
            "approved_by": "", "approved_on": ""}
    json.dump(prov, open(os.path.join(outs, "TISSUE_CALL_provenance.json"), "w"), indent=2)
    if os.path.exists(dst_final): os.rename(dst_final, dst_final + f"_superseded_{TS}")
    os.replace(stage, dst_final)
    manifest["outs"][sec] = f"{sec}_masked"
    # --- review figure ---
    col = {"spaceranger": "#BBBBBB", "he_mask": "#2CA02C", "expression": "#D62728", "excluded": "black", "contaminant": "#FF7F0E"}
    fig, ax = plt.subplots(1, 3, figsize=(20, 6.8))
    ax[0].imshow(img.astype(np.uint8))
    for k in ["spaceranger", "he_mask", "expression"]:
        m = source == k; ax[0].scatter(pc[m] * hs, pr[m] * hs, s=6 if k == "spaceranger" else 12, c=col[k], lw=0, alpha=.55 if k == "spaceranger" else .9, label=f"{k} (n={m.sum()})")
    m = source == "excluded"; ax[0].scatter(pc[m] * hs, pr[m] * hs, s=14, marker="x", c="black", lw=.6, label=f"excluded (n={m.sum()})")
    ax[0].legend(loc="lower left", fontsize=8, frameon=True); ax[0].set_title(f"{sec}: tissue call by tier"); ax[0].axis("off")
    ax[1].imshow(img.astype(np.uint8)); s = ax[1].scatter(pc * hs, pr * hs, s=6, c=np.log10(umi + 1), cmap="magma", vmin=1.5, vmax=4.3, lw=0)
    m = source != "spaceranger"; ax[1].scatter(pc[m] * hs, pr[m] * hs, s=12, facecolors="none", edgecolors="cyan", lw=.35)
    plt.colorbar(s, ax=ax[1], fraction=.03, label="log10 UMI"); ax[1].set_title("UMI; cyan = not in the Space Ranger call"); ax[1].axis("off")
    b = np.linspace(1, 4.6, 37)
    for k in ["spaceranger", "he_mask", "expression", "excluded", "contaminant"]:
        m = source == k
        if m.any(): ax[2].hist(np.log10(umi[m] + 1), bins=b, color=col[k], alpha=.55, label=f"{k} (median {med(m):,})")
    ax[2].axvline(np.log10(EXPR_UMI), ls=":", c="grey"); ax[2].set_xlabel("log10 UMI"); ax[2].set_ylabel("barcodes"); ax[2].set_yscale("log"); ax[2].legend(frameon=False, fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(FIGD, f"{sec}_tissue_call_tiers.png"), dpi=130); plt.close()
    for k in n:
        m = source == k
        summary.append(dict(section=sec, source=k, n=int(m.sum()), median_umi=med(m), median_genes=int(np.median(ng[m])) if m.any() else 0,
                            median_pct_mt=round(float(np.median(pmt[m])), 1) if m.any() else 0, median_he_fraction=round(float(np.median(hef[m])), 3) if m.any() else 0))
json.dump(manifest, open(os.path.join(VIS, "VISIUM_OUTS_MANIFEST.json"), "w"), indent=2)
with open(os.path.join(PROJ, "tables", "visium_tissue_call_summary_FH.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
try: os.rmdir(STAGE)
except OSError: pass
print("manifest:", os.path.join(VIS, "VISIUM_OUTS_MANIFEST.json"), "\n=== DONE ===")
