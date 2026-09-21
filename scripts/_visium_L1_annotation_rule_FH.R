# _visium_L1_annotation_rule_FH.R — Single source for how the digitised dashed lines of the pathologists' H&E annotation (tables/annotation_dashes_20260918_mapunits.csv; Fig6_260918_anatomy.pdf, 600 dpi) are grouped into components and
# labelled "L1" (layer-I lines) or "WM" (grey–white boundary lines). Sourced by 119, 121, 124; 119b (python) reads the
# JSON this file writes (tables/annotation_component_rule_20260918.json). Component = single-linkage cluster of the
# thinned pixels at H_LINK map units (~3 spot pitches); label = position of the component centroid as a fraction of the
# section's spot-cloud frame (x right, y up), by the rule below — the rule was set by eye against the labelled page.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
H_LINK <- 25
ANNOT_RULE <- list(
  CON_Frontal1 = "fx < 0.50 & fy > 0.50",                    # L1 = top-left arc; WM = lower-right arc
  CON_Frontal2 = "fx < 0.35 | fy < 0.25",                    # L1 = left edge + bottom; WM = top-centre loop
  NHD_Frontal1 = "fx < 0.50 & fy > 0.42",                    # L1 = the slit lines only; the U-band lines (fx ~0.59; fy ~0.33) are WM
  NHD_Frontal2 = "fy < 0.82 & fx > 0.15 & fx < 0.75")        # L1 = the diagonal band; WM = top-left / top-right corners
is_L1 <- function(sec, fx, fy) eval(parse(text = ANNOT_RULE[[sec]]))
if (exists("PROJ")) jsonlite::write_json(list(h_link = H_LINK, rule = ANNOT_RULE, frame = "centroid as fraction of the section spot-cloud range; x right, y up",
  source = "Fig6_260918_anatomy.pdf, domain-map row, 600 dpi; J.L. and S.Y., 2026-09-18"),
  file.path(PROJ, "tables", "annotation_component_rule_20260918.json"), auto_unbox = TRUE, pretty = TRUE)
