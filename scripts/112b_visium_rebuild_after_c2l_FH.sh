#!/bin/bash
# 112b_visium_rebuild_after_c2l_FH.sh — Poll until the local cell2location driver prints done, then run the Visium rebuild driver (112).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
: "${NHD_PROJ:?Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)}"
PROJ="$NHD_PROJ"
LOG="$(dirname "$NHD_PROJ")/Visium/cell2location/c2l_L1_driver_20260918.log"   # tee target of run_c2l_NHD_spatial.sh for the L1 run
until grep -q "] DONE" "$LOG" 2>/dev/null; do
  if grep -q "MISSING q05 . run failed\|MISSING inputs\|01f did not finish\|MISSING integrated object" "$LOG" 2>/dev/null; then echo "c2l FAILED — not starting 112"; exit 1; fi   # the driver prints a hyphen ("q05 - run failed"), the old em-dash pattern never matched; the 01f-step failures are caught too
  sleep 600
done
echo "[$(date +%F_%T)] c2l DONE — starting 112"
cd "$PROJ" && bash scripts/112_visium_downstream_rebuild_FH.sh
echo "[$(date +%F_%T)] 112 exit $?"
