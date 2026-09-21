# _visium_outs_FH.R — Single source for which Space Ranger outs folder each Visium section is read from.
# Reads Visium/VISIUM_OUTS_MANIFEST.json. Every R script that opens a section's outs
# sources this file and calls vis_outs("CON_Frontal1") instead of hard-coding "<sec>/outs" or
# "NHD_Frontal2_registered/outs". Python scripts read the same JSON directly.
# Exports: VIS_BASE (Visium folder), VIS_MANIFEST (parsed list), VIS_SECTIONS (chr), vis_outs(sec).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
.vis_candidates <- file.path(dirname(Sys.getenv("NHD_PROJ")), "Visium")
VIS_BASE <- .vis_candidates[dir.exists(.vis_candidates)][1]
stopifnot("Visium folder not found" = !is.na(VIS_BASE))
.mf <- file.path(VIS_BASE, "VISIUM_OUTS_MANIFEST.json")
stopifnot("VISIUM_OUTS_MANIFEST.json missing — run scripts/117b_visium_tissue_call_outs_FH.py" = file.exists(.mf))
VIS_MANIFEST <- jsonlite::fromJSON(.mf)
VIS_SECTIONS <- c("CON_Frontal1", "CON_Frontal2", "NHD_Frontal1", "NHD_Frontal2")
stopifnot("manifest lacks a section" = all(VIS_SECTIONS %in% names(VIS_MANIFEST$outs)))
vis_outs <- function(sec) {
  stopifnot(sec %in% VIS_SECTIONS)
  d <- file.path(VIS_BASE, VIS_MANIFEST$outs[[sec]], "outs")
  if (!dir.exists(d)) stop("outs folder missing for ", sec, ": ", d)
  d
}
