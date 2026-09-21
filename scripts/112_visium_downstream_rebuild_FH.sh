#!/bin/bash
# =============================================================================
# 112_visium_downstream_rebuild_FH.sh — Visium downstream rebuild driver: re-runs every step after integration + cell2location, in dependency order, with a gate on the deconvolution outputs.
# -----------------------------------------------------------------------------
# Upstream INPUTS (produced before this driver, not run here):
#   * Visium/VISIUM_OUTS_MANIFEST.json        117b — which outs folder each section is read from
#                                            (<sec>_masked/outs, three-tier tissue call; NHD_Frontal2's inherits the
#                                            reflected barcode map of NHD_Frontal2_registered, 108d/108e). Every step
#                                            resolves outs through _visium_outs_FH.R / the JSON, never a literal path.
#   * integrated_harmony.R                   -> NHD_frontal_integrated_harmony.rds
#   * integrated_cluster_identity.csv        curated cluster -> domain
#   * integrated_domain_levels.csv           domain order / palette / compartment / depth_rank (read by every step via _visium_domains_FH.R)
#   * cell2location 01f                      visium_{matrix,meta,coords} exported from the integrated object
#   * cell2location 02f/03e (run_c2l_NHD_spatial.sh, ~5 h)  -> c2l_MAIN/{q05_cell_abundance,
#                                            abundance_q05, proportions_q05}.csv + visium_c2l.h5ad
#
# Dependency ORDER (established from each script's inputs/outputs):
#   gate   cell2location finished (all three c2l tables newer than visium_meta.csv, row counts match,
#          no 02f process alive) — everything below reads proportions/abundance directly or transitively
#   01 108a  H&E crops + HE_frames.json            <- spot_coords.csv, manifest outs  (python ctm_env)
#   02 3M    myelin_scales_per_spot_visium_FH.csv  <- integrated rds (mtime-guarded; rebuilds because rds is new)
#   03 65    visium_dm_{marker_counts,module_spots,qc}_FH.csv  <- rds + proportions_q05 + spot_coords
#            (the cache KEY carries the input mtimes, so a new object/deconvolution rebuilds it)
#   04 66    6_depth_matched_markers               <- 65 marker table
#   05 67    F6e1/F6e2 spatial + 6_supplychain_matched + fig6_supplychain_FH.csv  <- 65 module table
#   06 71    6_tissue_depth_map, 6_depth_control, F6f_c2l_composition  <- 65 tables + abundance_q05
#   07 72    6_laminar_profile, 6_micro_vs_myelin, 6_crossmodality  <- 65 module table + 67's fig6_supplychain csv
#   08 73    F6d_supplychain_genes                 <- 65 module table + proportions_q05 (seurat_clusters) + identity csv
#   09 74    F6a/F6b/F6c                           <- visium_meta, spot_coords, identity csv, 65 module table, rds (UMAP cache mtime-guarded)
#   10 108b  F6b_visium_HE_domains                 <- 108a frames + visium_meta + spot_coords + identity csv
#   11 3L    4_myelin_scales_bridge                <- 3M spot table + proportions_q05
#   12 3P    4_visium_programme_reconcile          <- rds (score cache mtime-guarded) + proportions_q05
#   13 102   Data_deposition/Visium/spot_metadata.csv.gz + outs copy (+ tissue_call.csv per section)  <- rds + abundance_q05 + identity csv + manifest outs
#   14 111   SuppFig6 a proximity/severity           <- 65 module table + 102's spot_metadata.csv.gz + tissue_positions
#   15 114   SuppFig6 b/c composition per domain + ST21 <- 102's spot_metadata.csv.gz + tissue_positions + domain levels csv
#   16 suppfig4  Visium/figure_Visium/assemble_figure_Visium.R (VISIUM_SUPP_ONLY=1)  <- rds + identity csv;
#            then the four panels are copied into figures/Supplementary/SuppFig4_Visium_QC/panels/
#   17 108c  registration verdict — every section must print REGISTERED (exit != 0 otherwise)
#   18 115   QC summary tables quoted by Methods + Supp Fig 4 legend  <- visium_meta + integration log + tissue_call.csv + identity csv
#   not here: 90/100 (no supplementary table is Visium-derived — checked), 101/103/104/94.
#
# Usage
#   bash scripts/112_visium_downstream_rebuild_FH.sh                 run everything (stops at the first failure)
#   DRY=1 bash scripts/112_visium_downstream_rebuild_FH.sh           print the plan + gate status, run nothing
#   SKIP="108a 3M" bash scripts/112_visium_downstream_rebuild_FH.sh  skip steps by tag (space-separated)
#   only="65 66" bash scripts/112_visium_downstream_rebuild_FH.sh    run only these tags (in pipeline order)
#   FORCE_GATE=1 ...                                        bypass the cell2location gate (you know why)
# Per-step logs: logs/112_<nn>_<tag>.log ; master log: logs/112_visium_downstream_rebuild_FH.log
# macOS bash 3.2: no associative arrays / mapfile — plain indexed arrays and case statements.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
set -euo pipefail

# ---- roots (from NHD_PROJ) ---------------------------------------------------
: "${NHD_PROJ:?Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)}"
PROJ="$NHD_PROJ"
[ -n "$PROJ" ] || { echo "PROJ root not found"; exit 1; }
ROOT="$(dirname "$PROJ")"
SCR="$PROJ/scripts"; LOGD="$PROJ/logs"; mkdir -p "$LOGD"
C2L="$ROOT/Visium/cell2location/c2l_MAIN"
IH="$ROOT/Visium/integrated_harmony"
RSCRIPT="${RSCRIPT:-/usr/local/bin/Rscript}"
PY="${PY:-$HOME/miniforge3/envs/ctm_env/bin/python}"
DRY="${DRY:-0}"; SKIP="${SKIP:-}"; ONLY="${ONLY:-}"; FORCE_GATE="${FORCE_GATE:-0}"
MASTER="$LOGD/112_visium_downstream_rebuild_FH.log"

# gentle on a shared machine (cell2location may still be training beside this)
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 VECLIB_MAXIMUM_THREADS=2 KMP_DUPLICATE_LIB_OK=TRUE

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { printf '%s\n' "$*" | tee -a "$MASTER"; }
need() { [ -e "$1" ] || { log "MISSING $1 — $2"; exit 1; }; }
mtime() { stat -f %m "$1"; }
nrows() { awk 'END{print NR-1}' "$1"; }

# ---- the plan: "tag|kind|target|description" ------------------------------------------------------
STEPS=(
  "108a|py|$SCR/108a_visium_HE_crop_FH.py|H&E crops + HE_frames.json from the section outs + spot_coords"
  "3M|R|$SCR/3M_prep_opc_and_scales_FH.R|per-spot myelin share from the re-integrated object (mtime guard triggers)"
  "65|R|$SCR/65_visium_prep_depthmatched_FH.R|depth-matched spot tables (input-mtime key -> rebuild)"
  "66|R|$SCR/66_visium_markers_panel_FH.R|6_depth_matched_markers"
  "67|R|$SCR/67_visium_supplychain_spatial_FH.R|6_supplychain_matched (+ csv used by 72); F6e1/F6e2 depth-matched drafts"
  "71|R|$SCR/71_visium_fig6_overview_FH.R|6_tissue_depth_map / 6_depth_control / F6f_c2l_composition"
  "72|R|$SCR/72_visium_fig6_cascade_FH.R|6_laminar_profile / 6_micro_vs_myelin / 6_crossmodality"
  "73|R|$SCR/73_visium_gene_heatmap_FH.R|F6d_supplychain_genes (domains from identity csv)"
  "74|R|$SCR/74_visium_fig6_abc_FH.R|F6a / F6b / F6c (UMAP cache rebuilds from the new object)"
  "108b|R|$SCR/108b_visium_HE_panel_FH.R|F6b_visium_HE_domains (needs 108a frames)"
  "3L|R|$SCR/3L_myelin_scales_FH.R|4_myelin_scales_bridge (needs 3M spot table + c2l)"
  "3P|R|$SCR/3P_visium_programme_reconcile_FH.R|4_visium_programme_reconcile (score cache rebuilds)"
  "102|R|$SCR/102_export_deposition_bundle_FH.R|Data_deposition: spot_metadata.csv.gz + fixed outs copy (snRNA part reused)"
  "67b|R|$SCR/67b_visium_fig6e_alldepth_FH.R|F6e1/F6e2 at native depth on every spot + per-domain table (reads 102's spot_metadata for the domain labels -> must follow 102; moved 2026-09-18)"
  "111|R|$SCR/111_visium_proximity_severity_FH.R|proximity / severity (reads 102's spot_metadata)"
  "114|R|$SCR/114_visium_domain_composition_FH.R|S6b/S6c composition per domain + ST21 (reads 102's spot_metadata)"
  "suppfig4|assemble|$ROOT/Visium/figure_Visium/assemble_figure_Visium.R|SuppFig 4 panels a-d re-rendered + copied into figures/Supplementary/SuppFig4_Visium_QC/panels/"
  "108c|py|$SCR/108c_visium_HE_registration_check_FH.py|registration verdict: every section must be REGISTERED"
  "115|R|$SCR/115_visium_qc_summary_FH.R|QC numbers quoted in Methods / Supp Fig 4 legend (floor, medians, occupancy vs depth) -> tables/visium_qc_summary_FH.csv"
)

in_list() { # in_list needle "a b c"
  local n="$1" x; for x in $2; do [ "$x" = "$n" ] && return 0; done; return 1; }

# ---- gate: has cell2location finished on the current query? ------------------------------------------
gate() {
  local ok=1 f ref nmeta nab
  need "$IH/NHD_frontal_integrated_harmony.rds" "re-run integrated_harmony.R first"
  need "$IH/integrated_cluster_identity.csv"     "curate the identities first"
  need "$IH/integrated_domain_levels.csv"        "write the domain levels file first"
  need "$C2L/visium_meta.csv"                    "run cell2location 01f first"
  # The integrated object must postdate the tissue-call manifest (117b), else integrated_harmony.R was not re-run
  if [ "$(mtime "$IH/NHD_frontal_integrated_harmony.rds")" -lt "$(mtime "$ROOT/Visium/VISIUM_OUTS_MANIFEST.json")" ]; then
    log "  gate: NHD_frontal_integrated_harmony.rds is OLDER than Visium/VISIUM_OUTS_MANIFEST.json — re-run integrated_harmony.R"; ok=0; fi
  ref=$(mtime "$C2L/visium_meta.csv")
  for f in q05_cell_abundance.csv abundance_q05.csv proportions_q05.csv visium_c2l.h5ad; do
    if [ ! -s "$C2L/$f" ]; then log "  gate: $f missing"; ok=0
    elif [ "$(mtime "$C2L/$f")" -lt "$ref" ]; then
      log "  gate: $f ($(date -r "$(mtime "$C2L/$f")" '+%m-%d %H:%M')) is OLDER than visium_meta.csv ($(date -r "$ref" '+%m-%d %H:%M')) — cell2location has not re-exported"; ok=0
    fi
  done
  if pgrep -f "python.*02f_run_NHD_spatial\.py" >/dev/null 2>&1; then log "  gate: 02f_run_NHD_spatial.py is still running"; ok=0; fi
  if [ "$ok" = 1 ]; then
    nmeta=$(nrows "$C2L/visium_meta.csv"); nab=$(nrows "$C2L/abundance_q05.csv")
    if [ "$nmeta" != "$nab" ]; then log "  gate: abundance_q05.csv has $nab rows, visium_meta.csv has $nmeta — mismatch"; ok=0
    else log "  gate: cell2location tables are current ($nab spots, newer than visium_meta.csv)"; fi
  fi
  [ "$ok" = 1 ]
}

# ---- per-step runners ---------------------------------------------------------------------------------
run_R()  { "$RSCRIPT" "$1"; }
run_py() { "$PY" "$1"; }
run_assemble() {
  local src="$ROOT/Visium/figure_Visium" dst="$PROJ/figures/Supplementary/SuppFig4_Visium_QC/panels"
  local bak="$dst/_pre_2026-09-15" bak2="$dst/_pre_2026-09-18" p   # second snapshot folder, same once-only logic
  need "$dst" "SuppFig4 panel folder"
  ( cd "$src" && VISIUM_SUPP_ONLY=1 "$RSCRIPT" "$1" )
  mkdir -p "$bak" "$bak2"
  for p in panel_a_QC panel_b_UMAP_clusters panel_c_spatial_clusters panel_d_dotplot_clusters; do
    need "$src/$p.png" "assemble did not write $p.png"; need "$src/$p.pdf" "assemble did not write $p.pdf"
    [ -e "$dst/$p.png" ] && [ ! -e "$bak/$p.png" ] && cp -p "$dst/$p.png" "$dst/$p.pdf" "$bak/" 2>/dev/null || true
    [ -e "$dst/$p.png" ] && [ ! -e "$bak2/$p.png" ] && cp -p "$dst/$p.png" "$dst/$p.pdf" "$bak2/" 2>/dev/null || true   # pre-L1 copies, written once
    cp -p "$src/$p.png" "$src/$p.pdf" "$dst/"
    echo "  refreshed $p.{png,pdf} -> $dst  ($(sips -g pixelWidth -g pixelHeight "$dst/$p.png" 2>/dev/null | awk '/pixel/{printf "%s ", $2}'))"
  done
  local note="$PROJ/figures/Supplementary/SuppFig4_Visium_QC/PANELS_USED.md"
  if ! grep -q "2026-09-15 — NHD_Frontal2 fix" "$note"; then
    cat >> "$note" <<EOF

## 2026-09-15 — NHD_Frontal2 fix

All four panels re-rendered from the RE-INTEGRATED object (NHD_Frontal2 from \`NHD_Frontal2_registered/outs\`;
17,994 spots; TEN clusters g0–g9, identities from \`integrated_harmony/integrated_cluster_identity.csv\`,
see \`IDENTITY_CURATION_20260915.md\`). Previous copies kept in \`panels/_pre_2026-09-15/\`.
Cluster numbering CHANGED (e.g. WM is now cluster 2, L2/3 is cluster 1): the legend text
"nine spatial domains (clusters 0–8)" and the "b"/"d" rows above must be re-worded to ten clusters
(0–9), and \`Fig4_suppl.pdf\` must be re-linked/re-assembled by hand (script 104).
Canvas sizes are unchanged from the 2026-09-05 pass, so a plain relink is sufficient.
EOF
  fi
  if ! grep -q "2026-09-16 — Fig-6 frame" "$note"; then
    cat >> "$note" <<EOF

## 2026-09-16 — Fig-6 frame

Panel c is drawn as a ggplot in the Figure-6b frame (script 74's per-section normalised spot coordinates,
no H&E underlay) instead of Seurat's SpatialDimPlot; canvas 3.60 x 4.55 in unchanged. Domain order / palette
now come from \`integrated_harmony/integrated_domain_levels.csv\` (labels and colours unchanged).
EOF
  fi
  if ! grep -q "2026-09-18 — two-tier tissue call" "$note"; then   # L1 re-integration note, same once-only pattern
    cat >> "$note" <<EOF

## 2026-09-18 — two-tier tissue call + BANKSY-guided L1/pia
All four panels re-rendered from the object re-integrated on the final two-tier tissue call (every section read from the
outs named in \`Visium/VISIUM_OUTS_MANIFEST.json\`: Space Ranger call + relaxed H&E mask with the cortical-neuropil floor;
18,765 spots; 11 clusters g0-g10 after the script-120 refinement, identities from \`integrated_harmony/integrated_cluster_identity.csv\`,
see \`IDENTITY_CURATION_20260918.md\`; counts in \`tables/visium_qc_summary_FH.csv\`). (Panel e, the annotation overlay of script 124, was placed 18 Sep and retired 19 Sep.)
Previous copies kept in \`panels/_pre_2026-09-18/\`. Legend spot/cluster numbers re-derived from the new tables.
EOF
  fi
}

# ---- main -----------------------------------------------------------------------------------------------
log "=== 112_visium_downstream_rebuild_FH.sh  $(ts)  DRY=$DRY  SKIP='$SKIP'  ONLY='$ONLY' ==="
log "PROJ=$PROJ"
log "R=$RSCRIPT ($("$RSCRIPT" --version 2>&1 | head -1)) ; python=$PY"
[ -x "$PY" ] || { log "MISSING python env $PY (ctm_env) — 108a/108c need numpy + PIL"; exit 1; }

log "-- gate --"
if gate; then GATE=1; else GATE=0; fi
if [ "$GATE" = 0 ] && [ "$FORCE_GATE" != 1 ]; then
  if [ "$DRY" = 1 ]; then log "  (DRY) gate would BLOCK — wait for run_c2l_NHD_spatial.sh to print DONE"
  else log "GATE BLOCKED — cell2location has not finished on the fixed query. Wait, or FORCE_GATE=1."; exit 2; fi
fi

log "-- plan --"
i=0
for s in "${STEPS[@]}"; do
  i=$((i+1)); tag="${s%%|*}"; rest="${s#*|}"; kind="${rest%%|*}"; rest="${rest#*|}"; target="${rest%%|*}"; desc="${rest#*|}"
  state="run"
  if [ -n "$ONLY" ] && ! in_list "$tag" "$ONLY"; then state="skip(ONLY)"; fi
  if in_list "$tag" "$SKIP"; then state="skip(SKIP)"; fi
  [ -e "$target" ] || state="MISSING-SCRIPT"
  log "$(printf '  %02d  %-9s %-11s %-4s %s' "$i" "$tag" "$state" "$kind" "$desc")"
  [ "$state" = "MISSING-SCRIPT" ] && { log "  -> $target does not exist"; exit 1; }
done
[ "$DRY" = 1 ] && { log "=== DRY RUN — nothing executed ==="; exit 0; }

# provenance beside the logs
{
  echo "run: $(ts) host=$(hostname) user=$(whoami)"; echo "PROJ=$PROJ"; "$RSCRIPT" --version 2>&1 | head -1
  echo "identity csv:"; cat "$IH/integrated_cluster_identity.csv" | cut -d, -f1,2
  echo "domain levels csv:"; cat "$IH/integrated_domain_levels.csv"
  echo "inputs:"; for f in "$IH/NHD_frontal_integrated_harmony.rds" "$C2L/visium_meta.csv" "$C2L/spot_coords.csv" "$C2L/abundance_q05.csv" "$C2L/proportions_q05.csv"; do
    echo "  $(date -r "$(mtime "$f")" '+%Y-%m-%d %H:%M:%S')  $(stat -f %z "$f")  $f"; done
  echo "script md5:"; for s in "${STEPS[@]}"; do rest="${s#*|}"; rest="${rest#*|}"; t="${rest%%|*}"; md5 -q "$t" | sed "s|\$|  $(basename "$t")|"; done
} > "$LOGD/112_provenance_$(date +%Y%m%d_%H%M%S).txt"

i=0; t_all=$(date +%s)
for s in "${STEPS[@]}"; do
  i=$((i+1)); tag="${s%%|*}"; rest="${s#*|}"; kind="${rest%%|*}"; rest="${rest#*|}"; target="${rest%%|*}"; desc="${rest#*|}"
  if in_list "$tag" "$SKIP" || { [ -n "$ONLY" ] && ! in_list "$tag" "$ONLY"; }; then
    log "$(printf '[%02d] %-9s SKIPPED' "$i" "$tag")"; continue; fi
  lf="$LOGD/112_$(printf '%02d' "$i")_${tag}.log"
  t0=$(date +%s)
  log "$(printf '[%02d] %-9s START  %s  -> %s' "$i" "$tag" "$(ts)" "$(basename "$lf")")"
  set +e
  case "$kind" in
    R)        run_R "$target"        > "$lf" 2>&1 ;;
    py)       run_py "$target"       > "$lf" 2>&1 ;;
    assemble) run_assemble "$target" > "$lf" 2>&1 ;;
    *)        echo "unknown kind $kind" > "$lf"; false ;;
  esac
  rc=$?
  set -e
  dt=$(( $(date +%s) - t0 ))
  if [ "$rc" -ne 0 ]; then
    log "$(printf '[%02d] %-9s FAILED rc=%d after %ds — see %s' "$i" "$tag" "$rc" "$dt" "$lf")"
    log "---- last 25 lines of $(basename "$lf") ----"; tail -25 "$lf" | tee -a "$MASTER"
    log "=== STOPPED at step $i ($tag) $(ts) ==="; exit "$rc"
  fi
  # Every script in the plan prints "=== done ===" (65 prints "nothing to do" when its tables are
  # current). exit 0 without the sentinel = a silent quit() somewhere -> treated as a failure.
  if ! grep -q "=== DONE ===\|nothing to do" "$lf"; then
    log "$(printf '[%02d] %-9s FAILED exit 0 but no "=== DONE ===" sentinel in %s' "$i" "$tag" "$(basename "$lf")")"
    tail -15 "$lf" | tee -a "$MASTER"; log "=== STOPPED at step $i ($tag) $(ts) ==="; exit 3
  fi
  log "$(printf '[%02d] %-9s OK     %ds' "$i" "$tag" "$dt")"
done
log "=== ALL STEPS DONE in $(( $(date +%s) - t_all ))s  $(ts) ==="
