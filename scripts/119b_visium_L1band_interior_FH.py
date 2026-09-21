#!/usr/bin/env python3
# 119b_visium_L1band_interior_FH.py — Band interior of the pathologists' layer-I annotation.
# Rasterises the digitised L1-line pixels (tables/annotation_dashes_20260918_mapunits.csv; component -> label rule
# identical to 119_visium_sulcus_L1_check_FH.R) at 1 map unit (~12 um), closes the band with a 10-unit disk
# (~120 um; fills the two-line band and the dash gaps), and marks every spot whose centre lies inside.
# OUT: tables/visium_L1band_membership_20260918.csv (spot_id, section, domain, inside_L1_band, WM_or_L6_inside_band) + a printed
# profile of the WM/L6-labelled spots inside the band against the section's deep WM and its L1/pia spots.
# Python environment: ctm_env (python_env_ctm_env.txt).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
import numpy as np, csv
from scipy.ndimage import binary_closing, binary_dilation
from collections import Counter
import os
if "NHD_PROJ" not in os.environ:
    raise SystemExit("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
P=os.environ["NHD_PROJ"]
ps=list(csv.DictReader(open(f"{P}/tables/visium_sulcus_check_perspot_FH.csv")))
dash=list(csv.DictReader(open(f"{P}/tables/annotation_dashes_20260918_mapunits.csv")))
secs=["CON_Frontal1","CON_Frontal2","NHD_Frontal1","NHD_Frontal2"]
import json
RULE=json.load(open(f"{P}/tables/annotation_component_rule_20260918.json"))   # single source (written by scripts/_visium_L1_annotation_rule_FH.R)
def is_L1(sec,fx,fy):
    expr=RULE["rule"][sec].replace("&"," and ").replace("|"," or "); return bool(eval(expr, {"fx":fx,"fy":fy}))
H_LINK=RULE["h_link"]
rows=[]
for sec in secs:
    S=[r for r in ps if r["sample_id"]==sec]; xs=np.array([float(r["x"]) for r in S]); ys=np.array([float(r["y"]) for r in S])
    D=[(float(d["x"]),float(d["y"])) for d in dash if d["section"]==sec]
    xr=(xs.min(),xs.max()); yr=(ys.min(),ys.max())
    # component label by centroid: single-linkage components as in 119 (approximate by labelling each pixel through its own position is
    # not equivalent for long lines, so reproduce the component step)
    from scipy.cluster.hierarchy import fcluster, linkage
    pts=np.array(D); thin=pts[np.unique((pts/2).round(),axis=0,return_index=True)[1]]
    Z=linkage(thin,"single"); comp=fcluster(Z,H_LINK,criterion="distance")
    lab={}
    for c in np.unique(comp):
        cx,cy=thin[comp==c].mean(0); lab[c]="L1" if is_L1(sec,(cx-xr[0])/(xr[1]-xr[0]),(cy-yr[0])/(yr[1]-yr[0])) else "WM"
    L1=thin[[lab[c]=="L1" for c in comp]]
    x0=int(min(xs.min(),pts[:,0].min()))-20; y0=int(min(ys.min(),pts[:,1].min()))-20
    x1=int(max(xs.max(),pts[:,0].max()))+20; y1=int(max(ys.max(),pts[:,1].max()))+20
    G=np.zeros((y1-y0+1,x1-x0+1),bool)
    for x,y in L1: G[int(round(y))-y0,int(round(x))-x0]=True
    G=binary_dilation(G,iterations=1)
    yy,xx=np.ogrid[-10:11,-10:11]; C=binary_closing(G,structure=(xx**2+yy**2)<=100)
    inside=np.array([C[int(round(y))-y0,int(round(x))-x0] for x,y in zip(xs,ys)])
    dom=np.array([r["domain"] for r in S]); dL1=np.array([float(r["d_L1"]) for r in S])
    conflict=inside&np.isin(dom,["WM","L6"]); far=(~inside)&(dom=="WM")&(dL1>300); l1=inside&(dom=="L1/pia")
    mv=lambda k,m: round(float(np.mean([float(r[k]) for r,f in zip(S,m) if f])),1) if m.any() else None
    print(f"\n{sec}: inside band {int(inside.sum())} {dict(Counter(dom[inside]))}")
    for name,m in [("WM/L6 inside band (to re-assign)",conflict),("deep WM (>300 um)",far),("L1/pia inside band",l1)]:
        print(f"   {name:34s} n={int(m.sum()):4d} myelin {mv('myelin',m)} glia {mv('glia_limitans',m)} lepto {mv('leptomeningeal',m)} neurons {mv('pan_neuronal',m)} immune {mv('immune',m)} HBB {mv('HBB',m)}")
    rows+=[[r["spot_id"],sec,r["domain"],int(f),int(c)] for r,f,c in zip(S,inside,conflict)]
with open(f"{P}/tables/visium_L1band_membership_20260918.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["spot_id","section","domain","inside_L1_band","WM_or_L6_inside_band"]); w.writerows(rows)
print("\nwrote tables/visium_L1band_membership_20260918.csv")
