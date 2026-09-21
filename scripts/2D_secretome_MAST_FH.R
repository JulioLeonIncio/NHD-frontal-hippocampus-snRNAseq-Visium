# =============================================================================
# 2D_secretome_MAST_FH.R — Figure 2 panel d: microglial SECRETOME, MAST-forward, frontal + HIPPOCAMPUS (2 regions). Faithful port of _run_3I_secretome_MAST.R.
# -----------------------------------------------------------------------------
# Per-nucleus MAST contrast (NHD vs CON) of curated microglial secretome ligands
# across the 2 FH regions.  Atlas-free (reads the MAST discovery CSV only).
#
# Detection guardrail: a gene x region that is not in the region's MAST
#   background, or below the per-region detection floor (max(pct.1, pct.2) < 0.10),
#   is not a biological zero — it is drawn as a distinct grey open "x" at 0.
#
# gene-set re-validation:
#   Every curated secretome ligand was re-checked against the FH Micro-PVM MAST
#   tables.  The curated map below is unchanged; genes below the detection floor in
#   both FH regions are auto-dropped from the display (as in the reference).  On FH:
#     kept + moving (>=1 region): C1QA/B/C, F13A1, CTSB, PLAUR, GAS6, PSAP, PDGFC,
#       SEMA4D, HSP90AA1, HSPA1A, CD163, SPP1, FTL, FTH1, SLC11A1, TF.
#     Below-detection in both regions -> auto-dropped: CR1, LPL, FGF2, METRNL,
#       VEGFA, HSPB1, HMGB1, TREM1, COLEC12.
#     Detected but flat (kept if any region detected; grey "x" where undetected):
#       TGFB1, TNFSF13B, CXCL12, CH25H, CPM, ADAM17, SERPINE1, TIMP2, CHI3L1,
#       GPNMB, HMOX1.
#   Load-bearing biology confirmed: GAS6 + SPP1 co-lead present (GAS6 up both
#   regions; SPP1 up in Hippo, flat in Frontal); complement C1Q up; iron FTL/FTH1/
#   SLC11A1 up + TF down; PSAP up.  Consistent with DAM-1 arrest (no lipid-endpoint
#   secreted damp surge in cortex).
#
# Outputs:
#   figures/Figure_2/panels/F2f_secretome_MAST.{png,pdf}
#   diagnostics/10_microglia_core_signature/03I_secretome_MAST_values_FH.csv
# -----------------------------------------------------------------------------
# apparent-type lift 2D_secretome_MAST_FH.R: every text size in this script scaled
# by 1.14 so the panel reads at ~5.0 pt on the assembled page, matching the
# Figure-5 dot panels. Point sizes, line widths and unit() dimensions are not
# touched. Canvas size unchanged, so re-linking is a no-op.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  library(ggh4x); library(ragg)
})
set.seed(42)

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))

# microglia-only primary — see 40b/41b. Pooled table kept for sensitivity.
MAST_PATH <- file.path(PROJ, "tables", "mast_dual", "MAST_dual_discovery_all.csv")
THEME     <- file.path(PROJ, "scripts", "22_publication_theme_FH.R")
ARTIFACT  <- file.path(PROJ, "scripts", "_artifact_genes.R")
PANEL     <- file.path(PROJ, "figures", "Figure_2", "panels")
OUT_DIAG  <- file.path(PROJ, "diagnostics", "10_microglia_core_signature")
LOGS      <- file.path(PROJ, "logs")
for (d in c(PANEL, OUT_DIAG, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

for (f in c(MAST_PATH, THEME, ARTIFACT))
  if (!file.exists(f)) stop("MISSING ", f, " — upstream prep step did not run.")
source(THEME)       # theme_pub, PAL_DGE, PAL_REGION, REGION_FULL, REGION_ORDER
source(ARTIFACT)    # is_artifact()

CT      <- "Micro-PVM"   # POOLED microglia+PVM; PVM = 5.4% of compartment
REGIONS <- REGION_ORDER            # c("Frontal","Hippo")
PCT_FLOOR <- 0.10

# --- Curated secretome gene -> category map -
secretome_groups <- list(
  "Cytokine"      = c("TGFB1","TNFSF13B"),                    # baff
  "Chemokine"     = c("CXCL12"),                              # SDF-1
  "Complement"    = c("CR1","C1QB","C1QA","C1QC"),
  # F13A1 is an identity marker of the CD163+/F13A1+ myeloid population, not just a
  # secreted coagulation ligand, so it gets its own row group rather than sitting in
  # "Cytokine".
  # This row strip read "PVM/BAM (F13A1)", which asserts a
  # perivascular / border-associated macrophage origin for these cells on the face of a
  # main-figure panel. That origin is unresolved — CD163 can mark a microglial state —
  # and the house rule forbids asserting it. Script 2C (the Fig-2 volcano, same figure)
  # already calls this population "CD163+/F13A1+ state"; the two panels now agree
  # instead of giving one population two names.
  "CD163+/F13A1+" = c("F13A1"),
  "Lipid"         = c("CH25H","LPL","CST7","LGALS3"),   # LPL/CST7/LGALS3 = DAM-2 endpoint (kept visible even if below detection; see SHOW_AS_ABSENT)
  "Protease"      = c("CTSB","PLAUR","CPM","ADAM17"),         # ADAM17 sheddase
  "Prot. inhib."  = c("SERPINE1","TIMP2"),                    # PAI-1, MMP inhib
  "Trophic"       = c("GAS6","PSAP","PDGFC","FGF2","SEMA4D","METRNL","VEGFA"),  # GAS6 co-lead w/ SPP1
  "HSP / proteostasis" = c("HSP90AA1","HSPA1A","HSPB1"),
  "DAMP / activation"  = c("CHI3L1","CD163","GPNMB","SPP1","HMGB1"),
  "Soluble receptors"  = c("TREM1","COLEC12"),
  "Iron / metal"  = c("TF","FTL","FTH1","SLC11A1","HMOX1"))
sec_flat <- unique(unlist(secretome_groups, use.names = FALSE))
sec_cat  <- setNames(rep(names(secretome_groups), lengths(secretome_groups)),
                     unlist(secretome_groups, use.names = FALSE))
stopifnot(!any(duplicated(sec_flat)))
cat(sprintf("== curated secretome: %d genes / %d categories ==\n",
            length(sec_flat), length(secretome_groups)))

# --- Load MAST; restrict to Micro-PVM in the 2 FH regions ---------------------
mast <- read.csv(MAST_PATH, stringsAsFactors = FALSE)
stopifnot(all(c("cell_type","region","gene","avg_log2FC","p_val_adj",
                "pct.1","pct.2","cliffs_delta","discovery") %in% names(mast)))
mp <- mast %>% filter(cell_type == CT, region %in% REGIONS)
stopifnot(nrow(mp) > 0)

cat("== Microglia MAST background sizes per region ==\n")
for (r in REGIONS)
  cat(sprintf("  %-8s : %d genes tested\n", r,
              length(unique(mp$gene[mp$region == r]))))

mp <- mp %>% filter(!is_artifact(gene, CT))

grid <- expand_grid(gene = sec_flat, region = REGIONS)
df <- grid %>%
  left_join(mp %>% select(gene, region, avg_log2FC, p_val_adj,
                          pct.1, pct.2, cliffs_delta, discovery),
            by = c("gene", "region")) %>%
  mutate(
    min_pct  = pmin(pct.1, pct.2),
    max_pct  = pmax(pct.1, pct.2),
    detected = !is.na(avg_log2FC) & !is.na(max_pct) & max_pct >= PCT_FLOOR,
    dir = case_when(!detected            ~ NA_character_,
                    avg_log2FC > 0       ~ "NHD-up",
                    TRUE                 ~ "CON-up"),
    abs_cd  = abs(cliffs_delta),
    is_disc = detected & !is.na(discovery) & as.logical(discovery),  # padj<0.05 & |delta|>=0.15
    # Drop the ***/**/* stars
    # that were read straight off the per-nucleus MAST p_val_adj — that p is
    # pseudoreplicated (single NHD donor) and cannot carry a donor-level significance
    # claim.  Replace with the same effect-size encoding used on the volcano (script
    # 2C): the effect SIZE (Cliff's delta) sets a marker at the bar tip, gated
    # exactly like discovery.  "q*" = the house effect-gated mark (padj<0.05 and
    # |delta|<0.15, i.e. significant-but-small); "delta*" tier is drawn as a filled
    # effect-size dot whose SIZE = |Cliff's delta| (the strict/discovery tier).
    eff_tier = case_when(
      !detected | is.na(p_val_adj)          ~ "none",
      is_disc                               ~ "delta",   # padj<0.05 & |delta|>=0.15
      p_val_adj < 5e-2 & !is_disc           ~ "qstar",   # significant-but-small
      TRUE                                  ~ "none"),
    star = ifelse(eff_tier == "qstar", "q*", ""),        # only the house q* prints as text
    category = factor(sec_cat[as.character(gene)],
                      levels = names(secretome_groups)),
    gene     = factor(gene, levels = rev(sec_flat)),
    region   = factor(region, levels = REGIONS))

csv_path <- file.path(OUT_DIAG, "03I_secretome_MAST_values_FH.csv")
write.csv(df %>% arrange(category, gene, region) %>%
            mutate(gene = as.character(gene), region = as.character(region)),
          csv_path, row.names = FALSE)
cat(sprintf("Wrote values: %s (%d rows)\n", basename(csv_path), nrow(df)))

# Biology 3: show-as-absent.  LPL is the
# cardinal DAM-2-failure evidence — it (with CST7 / LGALS3, the DAM-2 lipid endpoint)
# must be visible as a greyed "below detection" row under Lipid, not silently dropped.
# These are the only genes exempted from the all-below-detection auto-drop; the other
# 8 non-DAM-2 drops (CR1/FGF2/METRNL/VEGFA/HSPB1/HMGB1/TREM1/COLEC12) still drop out.
SHOW_AS_ABSENT <- c("LPL","CST7","LGALS3")
# Show genes detected in >=1 region/condition; also always keep SHOW_AS_ABSENT rows.
keep_g <- df %>% group_by(gene) %>% summarise(a = any(detected), .groups = "drop") %>%
  filter(a | gene %in% SHOW_AS_ABSENT) %>% pull(gene) %>% as.character()
dropped <- setdiff(sec_flat, keep_g)
df <- df %>% filter(as.character(gene) %in% keep_g) %>%
  mutate(gene = droplevels(gene), category = droplevels(category))
cat(sprintf("kept %d/%d secretome genes detected in >=1 region; dropped: %s\n",
            length(keep_g), length(sec_flat),
            if (length(dropped)) paste(dropped, collapse = ", ") else "(none)"))

cat("\n== Per-region tally (curated secretome genes) ==\n")
tally <- df %>%
  mutate(class = case_when(!detected ~ "below-detection",
                           dir == "NHD-up" ~ "NHD-up",
                           TRUE ~ "CON-up")) %>%
  count(region, class) %>%
  pivot_wider(names_from = class, values_from = n, values_fill = 0)
print(as.data.frame(tally))

cat("\n== GAS6 + SPP1 co-lead + Iron/metal program (MAST) ==\n")
df %>% filter(gene %in% c("GAS6","SPP1","PSAP","FTL","FTH1","SLC11A1","TF")) %>%
  transmute(gene = as.character(gene), region = as.character(region),
            avg_log2FC = round(avg_log2FC, 2),
            max_pct = round(max_pct, 3),
            padj = formatC(p_val_adj, digits = 2, format = "g"),
            class = ifelse(!detected, "below-detection",
                           ifelse(dir == "NHD-up", "NHD-up", "CON-up"))) %>%
  arrange(gene, region) %>%
  write.table(file = stdout(), sep = "\t", row.names = FALSE, quote = FALSE)

# --- Render ------------------------------------------------------------------
L2FC_CAP <- 5
df_det <- df %>% filter(detected) %>%
  mutate(l2fc_cap = pmax(pmin(avg_log2FC, L2FC_CAP), -L2FC_CAP),
         # x anchor at the bar tip for the effect-size dot / text.  widen the tip gap 0.12 -> 0.28 so the value clears the
         # bar edge / zero-line (CHI3L1 Hippo -3.2 was reading into the zero line as
         x_tip  = l2fc_cap + ifelse(l2fc_cap > 0, 0.28, -0.28),
         hj_tip = ifelse(l2fc_cap > 0, 0, 1),
         # The grey house q* mark kept colliding
         # with the variable-width tip number (a fixed q* column could not know how wide
         # "%.1f" rendered).  Anchor q* just beyond the value text by estimating the
         # value-label width in x-units: n_chars * per-char width + a fixed gap, on the
         # outer side of the bar.  Keeps q* a distinct grey50 mark (not merged into the
         # value) with constant clearance from the number.
         val_lab = sprintf("%.1f", avg_log2FC),
         val_w   = nchar(val_lab) * 0.62 + 0.75,        # approx label width in x-units
         x_qstar = x_tip + ifelse(l2fc_cap > 0, val_w, -val_w),
         # x_qstar was unclamped and coord_cartesian(clip="off") drew it
         # outside the panel -- after the raw-RNA fix Hippocampus CHI3L1 moved out to
         # -3.2 and its q* was struck through by the facet's left border.  Simply
         # clamping it inward is not enough: on a long bar the value label already fills
         # the space between the tip and the frame, so the clamped mark lands on the
         # number ("q*-3.2").  Rule: keep the mark beyond the number when it fits, and
         # when it does not, flip it to the empty side of the zero axis on that row
         # (the bar occupies only one side, so the row stays unambiguous).
         qstar_outside = abs(x_qstar) > L2FC_CAP * 0.95,
         x_qstar = ifelse(qstar_outside,
                          ifelse(l2fc_cap > 0, -0.35, 0.35),
                          x_qstar),
         qstar_hjust = ifelse(qstar_outside,
                              ifelse(l2fc_cap > 0, 1, 0),
                              ifelse(l2fc_cap > 0, 0, 1)))
df_nd  <- df %>% filter(!detected)
df_eff <- df_det %>% filter(eff_tier == "delta")   # strict effect-size tier (discovery)

p <- ggplot() +
  geom_blank(data = df, aes(x = 0, y = gene)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(data = df_det,
           aes(x = l2fc_cap, y = gene, fill = dir),
           width = 0.72, colour = "grey30", linewidth = 0.15) +
  # numeric avg log2FC at every detected bar tip (no p-stars).
  geom_text(data = df_det,
            aes(x = x_tip, y = gene, hjust = hj_tip,
                label = sprintf("%.1f", avg_log2FC)),
            size = 1.7, colour = "grey15") +
  # effect-size encoding (replaces the p-value stars): a filled dot inside the bar,
  # near the zero axis, for the strict/discovery tier (padj<0.05 and |Cliff's delta|
  # >=0.15), SIZE = the effect magnitude |delta|.  Anchored at a fixed small offset
  # from 0 (into the bar, so it never collides with the tip number).  Mirrors the
  # volcano's two-tier effect encoding: a big dot = a large, gated effect.
  geom_point(data = df_eff,
             aes(x = ifelse(l2fc_cap > 0, 0.55, -0.55), y = gene,
                 size = abs_cd, fill = dir),
             shape = 21, stroke = 0.2, colour = "white",
             show.legend = c(fill = FALSE, size = TRUE)) +
  # house effect-gated "q*" mark for significant-but-small (padj<0.05 & |delta|<0.15),
  # anchored just beyond the tip number (x_qstar = value tip + estimated value width)
  # so it always clears the variable-width number with a constant grey gap.
  geom_text(data = df_det %>% filter(eff_tier == "qstar"),
            aes(x = x_qstar, y = gene, hjust = qstar_hjust, label = star),
            size = 1.7, colour = "grey50") +
  geom_point(data = df_nd, aes(x = 0, y = gene, shape = "below detection"),
             colour = "grey60", size = 1.7, stroke = 0.5) +
  facet_grid2(category ~ region, scales = "free_y", space = "free_y",
              switch = "y", labeller = labeller(region = as_labeller(REGION_FULL)),
              strip = strip_themed(
                background_y = elem_list_rect(fill = "grey97",
                                              colour = "grey80",
                                              linewidth = 0.25),
                text_y = elem_list_text(angle = 0, hjust = 1, size = 7.3,
                                        face = "plain", colour = "black"),
                # pale region banners + black labels (same lighten-0.55 blend as
                # the shared FH theme, so banners match the other Fig-2 panels).
                # The single-NHD-hippo-lane caveat is a legend note,
                # not an on-plot greyed banner.
                background_x = elem_list_rect(
                  fill = unname(PAL_REGION_PALE[REGIONS]), colour = NA),
                text_x = elem_list_text(colour = "black", size = 8.0,
                                        face = "plain"))) +
  scale_fill_manual(values = c(`NHD-up` = unname(PAL_COND[["NHD"]]), `CON-up` = unname(PAL_COND[["CON"]])),
                    name = NULL) +
  scale_shape_manual(values = c("below detection" = 4), name = NULL) +
  # effect-size dot size = |Cliff's delta|; fixed limits so the legend is stable
  # (0.15 = the discovery gate; only discovery genes carry a dot).
  # plotmath group("|", ., "|") renders TRUE vertical bars (the earlier literal
  # "|...|" string rendered as lowercase-l glyphs). ASCII apostrophe kept.
  scale_size_area(name = expression(group("|", "Cliff's " * delta, "|") ~ "(effect size)"),
                  max_size = 2.4, limits = c(0, 0.7),
                  breaks = c(0.2, 0.4, 0.6)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3),
                     expand = expansion(mult = c(0.25, 0.42))) +
  # clip="off": let the beyond-tip q* mark on the widest qstar bars (~2.7) draw into
  # the panel margin rather than being clipped at the +-5 viewport edge.
  coord_cartesian(xlim = c(-L2FC_CAP, L2FC_CAP), clip = "off") +
  guides(fill  = guide_legend(order = 1),
         size  = guide_legend(order = 2,
                   override.aes = list(fill = "grey55", colour = "white", shape = 21)),
         shape = guide_legend(order = 3, override.aes = list(size = 2.4))) +
  labs(x = expression(avg~log[2]~"FC (NHD/CON)"), y = NULL) +
  theme_pub(base_size = 9.1) +
  theme(
    plot.title         = element_blank(),
    axis.text.y        = element_text(size = 7.2, face = "italic", colour = "black"),
    axis.text.x        = element_text(size = 7.2, colour = "black"),
    axis.title.x       = element_text(size = 8.0, margin = margin(t = 3)),
    axis.line          = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA,
                                       linewidth = 0.25),
    panel.grid.major.x = element_line(colour = "grey95", linewidth = 0.2),
    strip.placement    = "outside",
    panel.spacing.x    = unit(0.25, "lines"),
    panel.spacing.y    = unit(0.15, "lines"),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    # The |Cliff's delta| size key ran off the
    # right edge (largest size dot + number clipped).  Wrap the bottom legend to 2
    # rows (legend.box="vertical" stacks fill / size / shape guides), widen the
    # canvas, and raise the right margin so the final size dot sits inside the frame.
    legend.box         = "vertical",
    legend.box.just    = "left",
    legend.key.size    = unit(0.22, "cm"),
    legend.spacing.x   = unit(0.2, "cm"),
    legend.spacing.y   = unit(0.02, "cm"),
    legend.margin      = margin(-2, 0, 0, 0),
    legend.text        = element_text(size = 7.1),
    plot.margin        = margin(6, 10, 4, 4))

# 2 region columns -> narrower.  Widened 3.35 -> 3.9 so the wrapped bottom
# legend / size key clears the right border.
W <- 3.9; H <- 6.0
png_path <- file.path(PANEL, "F2f_secretome_MAST.png")
pdf_path <- file.path(PANEL, "F2f_secretome_MAST.pdf")
ragg::agg_png(png_path, width = W, height = H, units = "in", res = 600)
print(p); dev.off()
ggsave(pdf_path, p, width = W, height = H, useDingbats = FALSE)
cat(sprintf("\nWrote panel: %s\n         and: %s\n", png_path, pdf_path))

# Provenance log next to outputs.
si_path <- file.path(OUT_DIAG, "03I_secretome_MAST_sessionInfo_FH.txt")
writeLines(capture.output({
  cat("Rendered:", format(Sys.time()), "\n")
  cat("PROJ:", PROJ, "\n")
  cat("MAST source:", MAST_PATH, "\n")
  cat("kept genes:", paste(keep_g, collapse = ", "), "\n")
  cat("dropped (below detection both regions):", paste(dropped, collapse = ", "), "\n\n")
  print(sessionInfo())
}), si_path)
cat(sprintf("Wrote provenance: %s\n", basename(si_path)))

cat("\n=== DONE ===\n", file = stderr())
