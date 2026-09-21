#!/bin/bash
# 107_visium_HE_fetch_FH.sh — Fetch the four full-resolution Visium H&E TIFFs from RIKEN.
#   FrontalA = NHD_Frontal1, FrontalB = NHD_Frontal2, FrontalC = CON_Frontal1, FrontalD = CON_Frontal2
# (macOS bash 3.2: no associative arrays — plain "section file" pairs.)
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
set -euo pipefail
: "${NHD_PROJ:?Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)}"
DST="$(dirname "$NHD_PROJ")/Visium/_HE_fullres"
SRC="/osc-fs_home/julio.l/analysis/NHD_visium/Repeat_2025/Frontal/20241217_Frontal"
mkdir -p "$DST"; : > "$DST/SOURCE_MAP.txt"
while read -r s f; do
  rsync -a --partial -e "ssh -o BatchMode=yes" "d8-work:$SRC/$f" "$DST/${s}_HE_fullres.tif"
  echo "$s <- $f  $(md5 -q "$DST/${s}_HE_fullres.tif")" | tee -a "$DST/SOURCE_MAP.txt"
done <<PAIRS
NHD_Frontal1 20241217121656_FrontalA.tif
NHD_Frontal2 20241217122010_FrontalB.tif
CON_Frontal1 20241217122249_FrontalC.tif
CON_Frontal2 20241217122453_FrontalD.tif
PAIRS
echo "=== DONE ==="
