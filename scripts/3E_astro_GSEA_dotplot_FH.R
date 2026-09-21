#!/usr/bin/env Rscript
# =============================================================================
# 3E_astro_GSEA_dotplot_FH.R — Figure 3 panel e: astrocyte GSEA dot-matrix, frontal + hippocampus. Port of the astro-GSEA part of 45_fig3_astro_oligo_MAST.R.
# ATLAS-FREE (reads the MAST-cliffs GSEA table only).
# -----------------------------------------------------------------------------
# One row per program x 2 region columns (full names).  Encoding (APOE panel-g /
# 3e house style):
#   fill  = NES (steel CON-up -> white -> rose NHD-up), clamp +/-3
#   size  = -log10(p-value), underflow-guarded + finite-capped (both logged)
#   alpha = FDR<0.05 solid vs nominal faded
# Same reviewer-facing blocklists as the reference (cilia/axonemal = ependymal
# admixture; off-target cardiac/muscle/vascular; neuronal-ambient synaptic;
# clock-mislabel) — the same junk terms recur in FH; kept.
#
# FH re-validation of the reference 5-program set
# against diagnostics/05_GSEA/.../fGSEA_NHD_MAST_cliffs_all.csv (Astro rows):
#   * OXPHOS / respiration  down both regions (FDR-sig Frontal, Hippo)   -> common.
#   * Ubiquitin-proteasome  down Frontal (FDR-sig)                        -> region.
#   * HSF1 / heat-shock     up Hippocampus (strongly FDR-sig)             -> region.
#   * Glutamate clearance   down both regions but only NOMINAL (padj~0.5-0.7)
#       -> kept as a nominal (faded) row: SLC1A2/GLUL are the most reproducible
#          astro DEGs, so the programme is shown honestly at nominal significance.
#   * Metallothionein / metal-detox  -> dropped: no significant metal-detox GSEA
#       term exists in the FH astro pool (the MT loss is a gene-level finding,
#       carried by the volcano + violins + dot-matrix, not a GSEA pathway).
#
# Output: figures/Figure_3/panels/F3d_astro_GSEA_dotplot.{pdf,png} + stats CSV.
# -----------------------------------------------------------------------------
# Harmonised to figure 2. Fig 3 had a 1.5x internal spread in type
# (axis.text 5.6 in the LR/receptor dotplots vs 8.5 in the dotmatrix/gas6 panels) and none of
# it matched Fig 2's house spec. All Fig-3 panels now use the same tokens as Fig 2's
# 2D_secretome_MAST reference -- base_size 8 / axis.text 6.3 / axis.title 7 / strip 6.4-7 /
# legend 6.2 -- rendered at their placed width so 6.3 pt is 6.3 pt on the page.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(scales); library(stringr); library(forcats)
  library(ggh4x)   # facet_grid2 + strip_region_x (sibling-consistent region banner)
})
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "2", KMP_DUPLICATE_LIB_OK = "TRUE")

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))   # theme_pub, REGION_FULL, REGION_ORDER

GSEA_MAST <- file.path(PROJ, "diagnostics", "05_GSEA",
                       "fGSEA_NHD_MAST_cliffs", "fGSEA_NHD_MAST_cliffs_all.csv")
PANEL <- file.path(PROJ, "figures", "Figure_3", "panels")
TDIR  <- file.path(PROJ, "tables")
LOGS  <- file.path(PROJ, "logs")
for (d in c(PANEL, TDIR, LOGS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(GSEA_MAST)) stop(sprintf("MISSING %s — GSEA prep step failed", GSEA_MAST))

REGIONS <- REGION_ORDER            # c("Frontal","Hippo")

LOG <- file.path(LOGS, "3E_astro_GSEA_dotplot_FH.log")
logcon <- file(LOG, open = "wt"); sink(logcon, split = TRUE)
on.exit({ sink(); close(logcon) }, add = TRUE)
cat("== 3E_astro_GSEA_dotplot_FH.R ==\nrun:", format(Sys.time()), "| PROJ:", PROJ, "\n\n")

gsea <- read.csv(GSEA_MAST, stringsAsFactors = FALSE)

FADE_N <- 150L
if ("power_tier" %in% names(gsea)) {
  n_excl <- sum(gsea$power_tier == "excluded", na.rm = TRUE)
  gsea <- gsea[is.na(gsea$power_tier) | gsea$power_tier != "excluded", , drop = FALSE]
  cat(sprintf("power filter: dropped %d excluded-tier rows (n_nuclei_min < 50)\n", n_excl))
}
if ("dge_member" %in% names(gsea)) {
  n0 <- nrow(gsea)
  gsea <- gsea[is.na(gsea$dge_member) | as.logical(gsea$dge_member), , drop = FALSE]
  cat(sprintf("dge-member display filter: %d -> %d rows\n", n0, nrow(gsea)))
}

ast <- gsea %>%
  filter(grepl("^Astro_", subtype), !is.na(NES)) %>%
  mutate(Region = sub("^Astro_", "", subtype),
         pathway_short = pathway_clean %>%
           str_replace_all(regex("\\bEndoplasmic reticulum\\b", ignore_case = TRUE), "ER") %>%
           str_replace_all(regex("\\bMitochondrial\\b", ignore_case = TRUE), "Mito") %>%
           str_replace_all("\\s+", " ") %>% str_trim()) %>%
  filter(Region %in% REGIONS)

# ---- blocklists (same as reference; the same junk recurs in FH) -------------
CILIA_BL <- regex("cilium|cilia|axonem|dynein arm|axonemal dynein|flagell|microtubule[-_ ]based[-_ ]movement|ciliary",
                  ignore_case = TRUE)
ast$.bl <- str_detect(ast$pathway_clean, CILIA_BL)
cat(sprintf("\ncilium/axonemal blacklist (ependymal admixture): %d term(s) removed (%d sig)\n",
            sum(ast$.bl), sum(ast$.bl & ast$padj < 0.05, na.rm = TRUE)))
ast <- ast %>% filter(!.bl) %>% select(-.bl)

CLOCK_MISLABEL <- regex("thermogenesis|temperature homeostasis", ignore_case = TRUE)
ast$.cm <- str_detect(ast$pathway_clean, CLOCK_MISLABEL)
if (any(ast$.cm)) cat(sprintf("clock-mislabel drop: %d term(s)\n", sum(ast$.cm)))
ast <- ast %>% filter(!.cm) %>% select(-.cm)

OFFTARGET <- regex(paste(
  "(?<!neuro)immune|antigen|mhc class|interferon|interleukin|complement|neutrophil",
  "leukocyte|(?<!neuro)inflammat|bacteri|infecti|(?<!anti)viral|\\bvirus\\b|influenza|leishmania",
  "endothel|(?<!neuro)epithel|vegf|angiogen|blood vessel|vasculature|vascular|blood circulation",
  "hemostasis|platelet|coagulation|wound healing|interspecies|\\btube\\b",
  "cilium|cilia|ciliary|axoneme|flagell|dynein|motile",
  "osteoclast|muscle contraction|cardiac conduction|myofibril|myoblast|smooth muscle|myogenesis",
  "hair cell|cochlea|cochlear|auditory|inner ear|\\bsound\\b", sep = "|"), ignore_case = TRUE)
ast$.ot <- str_detect(ast$pathway_clean, OFFTARGET)
ot_sig <- ast %>% filter(.ot, padj < 0.05)
cat(sprintf("off-target NAME blocklist: %d term(s) dropped (%d sig): %s\n",
            sum(ast$.ot), nrow(ot_sig),
            if (nrow(ot_sig)) paste(unique(ot_sig$pathway_clean), collapse = "; ") else "(none sig)"))
ast <- ast %>% filter(!.ot) %>% select(-.ot)

SYN_AMBIENT <- regex(paste(
  "neurexin|neuroligin", "synaptic vesicle|vesicle mediated.*synap|synaptic vesicle cycle|exocytosis",
  "presynap|active zone|neurotransmitter secretion|neurotransmitter loading",
  "protein protein interactions at synapses", sep = "|"), ignore_case = TRUE)
ast$.sa <- str_detect(ast$pathway_clean, SYN_AMBIENT)
sa_sig <- ast %>% filter(.sa, padj < 0.05)
cat(sprintf("neuronal-ambient synaptic blocklist: %d term(s) dropped (%d sig): %s\n",
            sum(ast$.sa), nrow(sa_sig),
            if (nrow(sa_sig)) paste(unique(sa_sig$pathway_clean), collapse = "; ") else "(none sig)"))
ast <- ast %>% filter(!.sa) %>% select(-.sa)

cat("\nper-region significance (padj<0.05):\n")
print(as.data.frame(ast %>% group_by(Region) %>%
        summarise(n_tested = n(), n_sig = sum(padj < 0.05, na.rm = TRUE), .groups = "drop")))
cat("NOTE: Hippocampus has a single NHD donor lane — direction-consistent support, not replication.\n\n")

# ---- FH-revalidated program set (metal-detox dropped; glutamate dropped) --------
# visual fix: the Glutamate-clearance/EAAT row is dropped from
# this dotplot. It carried only two tiny near-white non-significant dots (padj~0.5-0.7,
# NES~0) in a full-height facet, reading as a broken/empty row against the dense
# OXPHOS + Proteostasis rows. The glutamate-clearance message is already carried at
# the gene level by SLC1A2/GLUL in the MAST volcano and the state violins, so nothing
# is lost — the dotplot is tightened to the two categories that actually show GSEA
# signal. not a data change: the underlying GSEA table/gene-sets are untouched; this
# only removes a visually weak row from one display panel.
# Row labels shortened to one line each. At the requested 1.79 in height the
# previous wrapped 3-line names ("Oxidative phosphorylation / respiration" etc.) overlapped
# their neighbours. The full pathway names belong in the legend anyway.
PROGRAMS <- tibble::tribble(
  ~label,                                     ~category,           ~type,    ~pattern,
  "OXPHOS / respiration",  "Energy / OXPHOS",   "common", "oxidative phosphorylation|cellular respiration|aerobic respiration|respiratory electron|electron transport chain|proton transmembrane",
  "Ubiquitin-proteasome",          "Proteostasis",      "region", "proteasom|ubiquitin dependent|protein catabolic process|proteolysis involved",
  "HSF1 / heat-shock",               "Proteostasis",      "region", "hsf1|heat shock|heat stress|response to heat|attenuation phase|cellular response to heat",
  # The 3-row panel was missing the
  # strongest astrocyte result in the corrected GSEA. In Astro_Frontal the significant sets
  # are all CON-up (lost in NHD) and include regulation of synapse maturation (NES -2.01,
  # padj 0.0042), synapse organization (-1.62) and vesicle-mediated transport in synapse
  # (-1.77) -- i.e. astrocytes losing their synaptic support programme, which is both the
  # canonical reactive-astrocyte phenotype and the natural bridge to the neuron figure.
  # Protein transport/localization is the other consistently lost block.
  "Synapse organization",        "Synaptic support",  "region", "synapse organization|synapse maturation|synapse assembly|vesicle mediated transport in synapse|synaptic vesicle cycle|regulation of synapse structure",
  "Protein transport",         "Trafficking",       "region", "establishment of protein localization|protein transport|intracellular transport|localization within membrane")
# Not added: GOBP_CILIUM_MOVEMENT / MICROTUBULE_BASED_MOVEMENT, which are the top NHD-up sets
# in Astro_Hippo (NES 1.97 / 1.70). Motile-cilium programmes in an astrocyte annotation most
# likely reflect ependymal cells captured with astrocytes, not astrocyte biology. Flagged in
# the manifest rather than featured.
CAT_LEVELS <- c("Energy / OXPHOS", "Proteostasis", "Synaptic support", "Trafficking")

eaat_from_le <- function(le) {
  gg <- unlist(strsplit(as.character(le), "[|,/;\t ]+"))
  hit <- intersect(c("SLC1A2","SLC1A3","SLC1A4","GLUL"), gg)
  if (!length(hit)) return(NA_character_); paste(hit, collapse = ", ")
}
pick <- function(pat, rg) {
  m <- ast %>% filter(Region == rg, str_detect(pathway_clean, regex(pat, ignore_case = TRUE)))
  if (!nrow(m)) return(NULL); m %>% arrange(padj, desc(abs(NES))) %>% slice(1)
}
rows <- list()
for (i in seq_len(nrow(PROGRAMS))) for (rg in REGIONS) {
  m <- pick(PROGRAMS$pattern[i], rg)
  rows[[length(rows) + 1]] <- tibble(
    label = PROGRAMS$label[i], category = PROGRAMS$category[i], type = PROGRAMS$type[i], Region = rg,
    NES  = if (is.null(m)) NA_real_ else m$NES,
    pval = if (is.null(m)) NA_real_ else m$pval,
    padj = if (is.null(m)) NA_real_ else m$padj,
    rep_term    = if (is.null(m)) NA_character_ else m$pathway_clean,
    leadingEdge = if (is.null(m)) NA_character_ else m$leadingEdge)
}
grid <- bind_rows(rows)

# EAAT annotation — no-op now the Glutamate row is dropped from this panel (kept as
# a guarded block so re-adding the row later would restore the SLC1A2/GLUL suffix).
eaat_rep <- grid %>% filter(category == "Glutamate / EAAT", !is.na(leadingEdge))
eaat_str <- if (nrow(eaat_rep)) eaat_from_le(paste(eaat_rep$leadingEdge, collapse = "|")) else NA_character_
if (!is.na(eaat_str) && nzchar(eaat_str)) {
  base_lab <- PROGRAMS$label[PROGRAMS$category == "Glutamate / EAAT"][1]
  grid$label[grid$category == "Glutamate / EAAT"] <- paste0(base_lab, " (", eaat_str, ")")
}
cat(sprintf("glutamate row dropped from panel; EAAT leading-edge (n.s.) not shown: %s\n",
            if (!is.na(eaat_str)) eaat_str else "(none)"))

coh <- grid %>% group_by(label, category, type) %>%
  summarise(Frontal = NES[Region == "Frontal"][1], Hippo = NES[Region == "Hippo"][1],
            n_region_sig = sum(padj < 0.05, na.rm = TRUE), .groups = "drop") %>%
  arrange(match(category, CAT_LEVELS), type != "common")
cat("\n== FINAL astro GSEA program set (per-region NES; type) ==\n")
print(as.data.frame(coh %>% mutate(across(c(Frontal, Hippo), ~round(.x, 2)))), row.names = FALSE)

# ---- degenerate-value guards (reviewer-facing) -----------------------------
XMIN <- .Machine$double.xmin; CAPL <- 5
n_under <- sum(!is.na(grid$pval) & grid$pval > 0 & grid$pval < XMIN)
raw_nlp <- ifelse(is.na(grid$pval), NA_real_, -log10(pmax(grid$pval, XMIN)))
n_cap   <- sum(!is.na(raw_nlp) & grid$pval < 0.05 & raw_nlp > CAPL, na.rm = TRUE)
cat(sprintf("\nsize guard: %d p-value(s) underflowed (floored to xmin); %d nominal cell(s) capped at -log10(p)=%d.\n",
            n_under, n_cap, CAPL))

grid <- grid %>% mutate(
  nom_sig = !is.na(pval) & pval < 0.05,
  fdr_sig = !is.na(padj) & padj < 0.05,
  nes_c   = ifelse(is.na(NES), 0, pmax(pmin(NES, 3), -3)),
  dotsize = ifelse(nom_sig, pmin(-log10(pmax(pval, XMIN)), CAPL), 0),
  category    = factor(category, levels = CAT_LEVELS),
  # region_col = short codes (Frontal/Hippo) so strip_region_x(REGIONS, clip = "off") matches the
  # pale banner fills by name; the labeller maps them to full Frontal/Hippocampus.
  region_col  = factor(Region, levels = REGIONS),
  Region_full = factor(unname(REGION_FULL[Region]), levels = unname(REGION_FULL[REGIONS])))

lab_order <- grid %>% filter(!is.na(padj)) %>% group_by(label, category, type) %>%
  summarise(best = min(padj, na.rm = TRUE), .groups = "drop") %>%
  arrange(match(category, CAT_LEVELS), type != "common", best) %>% pull(label)
lab_order <- c(lab_order, setdiff(unique(grid$label), lab_order))
grid$label <- factor(grid$label, levels = rev(lab_order))
ylab_map <- setNames(str_wrap(levels(grid$label), width = 24), levels(grid$label))

GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
# Was: facet_grid(rows = category) with Frontal/Hippocampus dropped to rotated 30-deg
# x-axis tick text and grey (grey90) vertical category strips on the left.
# Now: facet_grid2(category ~ Region) with the TOP region strip supplied via
# strip_region_x(REGIONS, clip = "off") so Frontal / Hippocampus get the same pale PAL_REGION_PALE
# banner + black plain text as the siblings; category names stay on the left but their
# strip is harmonised to the sibling look (plain black text, no grey fill). Region is
# taken from the banner, not from rotated axis text.
strip_reg <- strip_region_x(REGIONS, fontsize = 7, clip = "off")
p4b <- ggplot(grid, aes(x = Region, y = label)) +
  geom_point(aes(fill = nes_c, size = dotsize, alpha = fdr_sig), shape = 21, colour = "grey35", stroke = 0.25) +
  # round-4 visual fix: the two-category facet
  # (OXPHOS | Proteostasis) split the 3 programs across two vertically-stacked boxes,
  # and even with force_panelsizes(1,2) the single-program OXPHOS box read as a tall
  # rectangle with one dot floating at mid-height (whitespace-dominated vs the dense
  # sibling dot-matrices). Collapsed to one facet row per region (region-only facet,
  # like the 3d reactive dot-matrix sibling), so all 3 programs sit in a single
  # contiguous column at consistent row spacing and the OXPHOS dot no longer floats.
  # The category grouping is not lost: it is preserved as the row ORDER (OXPHOS on
  # top, then the two Proteostasis rows) and stays in the source CSV `category` column.
  # Not a data change — same programs, same NES/p, only the facet layout changed.
  # x is mapped to Region while
  # the panel also facets on region, and the x scale was not free -- so every facet retained
  # both region levels and its single dot sat at that level's slot (Frontal at slot 1 = left of
  # centre, Hippocampus at slot 2 = right) instead of in the middle of its own panel.
  # scales = "free" drops the unused level per facet so each dot centres.
  facet_grid2(. ~ region_col, scales = "free", space = "free_y",
              labeller = labeller(
                region_col = ggplot2::as_labeller(function(x)
                  ifelse(x %in% names(REGION_FULL), unname(REGION_FULL[x]), x))),
              strip = strip_reg) +
  # Single-line break labels (were 2-line "+3\n(NHD-up)") so the top break is compact,
  # and the "NES" title is lifted clear of the bar with title.vjust + a bottom margin
  # (theme legend.title) — fixes the "NES+3" title/top-break collision that read as a
  # garbled subscript.
  scale_fill_gradient2(low = GSEA_STEEL, mid = "white", high = GSEA_ROSE, midpoint = 0,
                       limits = c(-3, 3), breaks = c(-3, 0, 3),
                       # Break labels compacted to bare numbers. On the new
                       # horizontal bottom legend the words made the bar ~2x wider than the
                       # bar itself, so "+3 NHD-up" ran into the size key and pushed the
                       # alpha key off the right edge. The direction now lives in the title,
                       # which costs no extra width.
                       labels = c("-3", "0", "+3"),
                       # Title back to bare "NES". "(+ = NHD-up)" was the widest
                       # single element on the bottom legend row and was the last thing pushing the
                       # FDR key off the canvas -- and it is caption text, which the house rules put
                       # in the figure legend, not the panel. The legend must therefore state that
                       # positive NES = enriched in NHD.
                       name = "NES") +
  # 1e-5) and the ALPHA aesthetic separately encodes FDR<0.05. The review note asked
  # to rename this "BH q" to match the delta dot-matrices, but those size on
  # -log10(BH q) — sizing this on p and calling it "q" would be a factual mislabel a
  # reviewer would catch. Kept the encoding, made the label unambiguous ("nominal p")
  # GSEA should instead size on padj.
  # Wider size range so the significant dots carry more visual weight
  # against the sparse rows. Encoding unchanged.
  # Size floor lifted 0.9 -> 1.8: the smallest significant dot
  # (Hippo Ubiquitin-proteasome, -log10 p just above the 0.05 gate) rendered almost
  # invisible; a higher floor keeps every nominally-significant dot legible while the
  # ceiling stays 5.4 so the encoding contrast is preserved.
  scale_size_continuous(range = c(1.8, 5.4), limits = c(0, CAPL),
                        # Three size breaks were the widest block on the
                        # single-row bottom legend and pushed the FDR alpha key off the canvas.
                        # Two breaks still convey the scale (the range is what matters here).
                        breaks = c(1.3, 5), labels = c("0.05", "1e-5"), name = "nom. p") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4),
                     # Labels shortened with the 20% slim (4.15 -> 3.32 in). The
                     # bottom legend carries three blocks on one row and the alpha key is last,
                     # so it is the one that clips; "FDR<.05" / "p<.05" buys the two characters
                     # that keep it inside the canvas without a second legend row.
                     # Its keys therefore
                     # sat on the TITLES' baseline instead of the keys' baseline and the bottom
                     # legend read as three blocks scattered at two different heights. A single
                     # space reserves the row so all three key rows align.
                     labels = c(`TRUE` = "FDR<.05", `FALSE` = "p<.05"), name = " ") +
  # tighter y expansion (add=0.55) so the 3 rows sit at even spacing filling the
  # single facet box — no tall stranded cell now that the category split is gone.
  scale_y_discrete(labels = function(x) ylab_map[x], expand = expansion(add = 0.42)) +
  # These panels never
  # set scale_x_discrete, so ggplot's default discrete expansion (0.6 units either side) padded
  # every facet far wider than its dot column -- with a single receiver that is ~0.25 in of white
  # on each side of every dot. Tighten the padding to the dots.
  scale_x_discrete(expand = expansion(add = 0.30)) +
  # Legend moved to a horizontal bottom strip. Stacked on the right it needed
  # ~2.2 in of height for its three blocks (NES bar / size key / alpha key), so at the
  # requested 1.79 in the bottom block clipped off the canvas. Laid out horizontally it costs
  # ~0.4 in of height instead, which is what makes the shorter panel possible.
  # The size key's three breaks were the widest block on the single-row
  # bottom legend and were pushing the FDR alpha key off the canvas. Two breaks still convey the
  # scale, "nominal p" -> "nom. p", legend text 6.2 -> 5.8. This keeps the requested width
  # without adding a second legend row (which would have cost ~0.25 in of height).
  # Three guides laid side by side with
  # title.position="top" produced two staggered baselines -- each block sized its title to its
  # own key height, so "NES", "nom. p" and the alpha dots all sat at different y. Re-laid as a
  # vertical stack of three single rows with INLINE (left) titles: each guide is now one tidy
  # line, the block is narrow (which is what lets the panel go 20% slimmer), and there is only
  # one baseline per row instead of two per block.
  guides(fill  = guide_colourbar(order = 1, direction = "horizontal",
                                 barwidth = unit(1.2, "cm"), barheight = unit(0.24, "cm"),
                                 title.position = "left", title.vjust = 1),
         size  = guide_legend(order = 2, nrow = 1, title.position = "left", title.vjust = 0.5),
         alpha = guide_legend(order = 3, nrow = 1, title.position = "left", title.vjust = 0.5,
                              override.aes = list(fill = "grey40", size = 2.2))) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 8) +
  theme(
    # region now lives in the TOP pale banner (strip_region_x) -> no rotated x tick text
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6.3, color = "black"), axis.line = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
    panel.spacing.x = unit(3, "pt"),
    strip.placement = "outside", strip.clip = "off",
    # round-4: single region-only facet now — no left category strip to theme.
    # Legend moved to the right. Kept as a vertical stack with inline
    # (left) titles -- that is what removed the staggered baselines -- so it reads as three
    # tidy rows in a right-hand column rather than three blocks scattered under the plot.
    # Legend on the right. Vertical stack
    # of three single rows with INLINE (left) titles -- that is what removed the staggered
    # baselines. A bottom row also re-clipped "FDR<.05" at this width, so right is the safer
    # placement for this three-guide legend.
    legend.position = "right", legend.box = "vertical", legend.margin = margin(l = -2),
    legend.box.just = "left", legend.justification = "center",
    legend.spacing.x = unit(0.18, "cm"),
    legend.spacing.y = unit(0.02, "cm"), legend.key.height = unit(0.24, "cm"),
    legend.key.size = unit(0.4, "lines"), legend.text = element_text(size = 5.8),
    # margin r=3: with INLINE (left) titles the label sat flush against its key
    # ("NES" touching the colourbar). b=4 kept from the old top-title layout.
    legend.title = element_text(size = 6.2, margin = margin(b = 4, r = 3)),
    panel.grid.major.y = element_line(linewidth = 0.1, colour = "grey95"),
    panel.grid.major.x = element_blank(), plot.margin = margin(6, 4, 4, 4))

# H 2.0 -> 2.5.  The right-hand legend is a 3-block stack (NES colourbar +
# nominal-p size key + significance alpha key); at H=2.0 the last alpha row ("FDR<0.05")
# was sliced in half by the canvas edge and "p<0.05" sat on top of it.  The legend, not
# the plot, sets the minimum height here.  Also tighten the y expansion so the extra
# height goes to the legend rather than re-opening dead space between the 3 rows.
# Width FLOOR: the bottom legend carries three blocks on one row (NES bar, nominal-p size key,
# FDR alpha key). Below ~4.15 in the alpha key clips ("FDR<0.05" -> "FDR"). Panel area is still
# 7.4 in2 vs the 9.9 in2 it started at, and height is the requested 1.79.
# Moving the legend to the right costs ~1.2 in of width, which left each facet
# ~0.48 in -- under the ~0.68 in "Hippocampus" needs at 7 pt, so the banners collided
# ("FrontaHippocampus"). Widened to pay for the right-hand legend column; the dots keep the
# tightened 0.30 x-expansion.
# Right-hand legend column needs ~1.2 in; below ~3.70 the facets drop under the ~0.68 in that
# "Hippocampus" needs at 7 pt and the banners collide.
W4b <- 3.70; H4b <- 1.55
                           # the 3-row single facet was still whitespace-dominated — the 6
                           # dots swam in two tall boxes, reading far weaker than the dense
                           # sibling dot-matrices. Cut the height so the 3 program rows fill
                           # a shorter box and the dots sit at tight, even spacing (dots no
                           # longer floating). Width kept at 3.7 so the right legend stack
                           # (colourbar + 3 size breaks + p<0.05/FDR<0.05 alpha) still
                           # clears. not a data change — same 3 programs / NES / p, layout
                           # only; the add=0.55 y-expansion keeps the top/bottom rows off
                           # the panel border.
BN <- "F3d_astro_GSEA_dotplot"
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p4b, width = W4b, height = H4b, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p4b, width = W4b, height = H4b, dpi = 600, device = ragg::agg_png)
cat("\nWrote", BN, ".{pdf,png}\n")

write.csv(
  grid %>% transmute(program = as.character(label), category = as.character(category), type,
              Region, NES = round(NES, 3), pval = signif(pval, 3), padj = signif(padj, 3),
              fdr_sig, nom_sig, rep_term) %>%
    arrange(match(category, CAT_LEVELS), program, match(Region, REGIONS)),
  file.path(TDIR, "fig3e_astro_gsea_MAST_stats.csv"), row.names = FALSE)
cat("Wrote tables/fig3e_astro_gsea_MAST_stats.csv\n")
cat("=== DONE ===\n", file = stderr())
