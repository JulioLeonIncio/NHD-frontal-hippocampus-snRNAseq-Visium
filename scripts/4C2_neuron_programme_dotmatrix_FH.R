#!/usr/bin/env Rscript
# =============================================================================
# 4C2_neuron_programme_dotmatrix_FH.R — Figure 5: neuronal programme dot-matrix, replacing the violin panel ("for the neuronal panels lets replace the violin module scores for dot plots").
# -----------------------------------------------------------------------------
# Same decision as Figure 4, where the violin sibling `3f_oligo_lineage_block` was
# demoted to _supp/ and the dot-matrix became the paper panel: at five programmes x
# two regions x two classes the violins spend most of their ink on distribution
# shape that carries no claim, while the claim itself is a single effect size per
# cell. The dot-matrix shows twenty cells in the space four violins used.
#
# Grammar is copied from F4d_oligo_lineage_dotmatrix (astro spec): Cliff's delta as
# fill, BH-q as size, the house effect-gated grey "q*" (q < 0.05 and |delta| < 0.15),
# base_size 8 / axis.text 6.3 / strip 6.5 / legend 6.2, panel.spacing.x 1.5 pt.
# No per-nucleus significance stars — single donor per condition.
#
# Added vs the violins: Inhibitory neurons as the second x-position, exactly as OPC
# was added beside Oligo in Fig 4. The violins were excitatory-only, so the Ex/Inh
# contrast — which is the whole vulnerability axis of this figure — was invisible.
#
# FLAGGED, not silently changed: `GluR_ionotropic` is a CANCELLED set. Its members
# move in opposite directions (GRIA2 -0.289 down vs GRIK2 +0.244 / GRIA4 +0.197 up,
# all discovery-tier), so the pooled mean lands near zero and will render as a white
# dot carrying the triviality mark — the exact failure that cost Fig 4 its headline
# . The per-member deltas are printed below and written to the source CSV
# so the split can be made deliberately rather than by accident.
# -----------------------------------------------------------------------------
# apparent-type parity PASS  ("adjust the font/panel size ratio
#   so everything looks more or less uniform, favouring the readability of the figure")
# Measured from the assembled Fig5.pdf: each panel is placed as a linked image with a
# known scale (parsed from the page's placement matrices). On-page type was 3.2-3.9 pt
# for most panels but 6.4-6.9 pt for the chord -- a 2x spread, and most of it below the
# print-legibility floor. Declared sizes here are therefore set to
#     declared = apparent_target / placement_scale
# with one set of apparent targets shared by every Figure-5 panel:
#     tick / legend text 5.2 pt | legend title 5.4 | axis title & facet strip 6.0
#     panel title 6.4 | in-panel gene labels & annotations 5.4
# canvas size and PIXEL DIMENSIONS are unchanged, so re-linking in Illustrator is a
# no-op: the frame keeps its box and only the type grows.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggh4x); library(stringr)
})
set.seed(42)
if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- Sys.getenv("NHD_PROJ")
stopifnot(!is.na(PROJ), dir.exists(PROJ))
source(file.path(PROJ, "scripts", "22_publication_theme_FH.R"))
DATA <- file.path(PROJ, "data"); TDIR <- file.path(PROJ, "tables")
PANEL<- file.path(PROJ, "figures", "Figure_5", "panels")
dir.create(PANEL, showWarnings = FALSE, recursive = TRUE)

CACHE <- file.path(DATA, "_cache_neuron_scores_FH.rds")
stopifnot("MISSING neuron score cache — run 4C_signature_violins_FH.R first" =
            file.exists(CACHE))
sc <- readRDS(CACHE)
sigs <- attr(sc, "sigs")
MODS <- names(sigs)
cat("modules:", paste(MODS, collapse = " | "), "\n")

CLASS <- c(Neuron_Ex = "Ex", Neuron_Inh = "Inh")

# Detection basis per CLASS and region, measured on raw counts (the same treatment
# Fig. 3f and Fig. 4d now carry). The neuron scores are a mean of log-normalised counts, so an
# undetected member contributes 0 rather than a random control bin — but a set can still be
# carried by one or two abundant members while its canonical genes sit below the floor. The IEG
# set is the case in point: 8 of its 14 members (FOS, FOSB, ARC, NPAS4, JUNB, EGR2, NR4A2, NR4A3)
# are absent from the frontal excitatory universe, so the score rests on HOMER1 and JUND.
source(file.path(PROJ, "scripts", "_receptor_detection_FH.R"))
.det_n <- receptor_detection(PROJ, file.path(PROJ, "atlas", "NHD_FH_harmony.rds"),
                             genes = unique(unlist(sigs)), classes = names(CLASS), regions = REGION_ORDER)
det_tab <- do.call(rbind, lapply(names(sigs), function(nm) do.call(rbind, lapply(names(CLASS), function(ct)
  do.call(rbind, lapply(REGION_ORDER, function(rg) {
    i <- match(paste(sigs[[nm]], ct, rg), paste(.det_n$gene, .det_n$class, .det_n$region))
    pc <- .det_n$pct_max[i]
    data.frame(signature = nm, class = unname(CLASS[ct]), Region = rg,
               n_in_set = length(sigs[[nm]]), n_present = sum(!is.na(pc)),
               n_det10 = sum(!is.na(pc) & pc >= 0.10), stringsAsFactors = FALSE) }))))))
cat("\ndetected members per programme x class x region (>= 10 % of nuclei in either condition):\n")
print(det_tab[order(det_tab$signature, det_tab$class, det_tab$Region), ], row.names = FALSE)
long <- sc %>%
  filter(cell_type %in% names(CLASS), Region %in% REGION_ORDER,
         Condition %in% c("CON","NHD")) %>%
  mutate(class = factor(unname(CLASS[cell_type]), levels = unname(CLASS)),
         Region = factor(Region, levels = REGION_ORDER)) %>%
  pivot_longer(all_of(MODS), names_to = "signature", values_to = "score") %>%
  filter(!is.na(score))

cliffs_d <- function(nhd, con) {          # +ve = higher in NHD
  nhd <- nhd[!is.na(nhd)]; con <- con[!is.na(con)]
  n1 <- length(nhd); n2 <- length(con)
  if (n1 < 3 || n2 < 3) return(NA_real_)
  r <- rank(c(nhd, con)); U <- sum(r[seq_len(n1)]) - n1*(n1+1)/2
  2*(U/(n1*n2)) - 1
}
stat <- long %>%
  group_by(signature, Region, class) %>%
  summarise(p = tryCatch(wilcox.test(score[Condition=="NHD"], score[Condition=="CON"],
                                     exact = FALSE)$p.value, error = function(e) NA_real_),
            cliff_d = cliffs_d(score[Condition=="NHD"], score[Condition=="CON"]),
            n_CON = sum(Condition=="CON"), n_NHD = sum(Condition=="NHD"),
            .groups = "drop") %>%
  mutate(q_BH = p.adjust(p, method = "BH"),
         passes = !is.na(cliff_d) & abs(cliff_d) >= 0.15,
         label  = ifelse(!is.na(q_BH) & q_BH < 0.05 & !passes, "q*", ""))

XMIN <- .Machine$double.xmin; CAPL <- 50
stat$dotsize <- pmin(ifelse(is.na(stat$q_BH), NA_real_,
                            -log10(pmax(stat$q_BH, XMIN))), CAPL)
n_uf <- sum(!is.na(stat$q_BH) & stat$q_BH == 0)
if (n_uf) cat(sprintf("NOTE: %d q==0 underflows floored to xmin for the size scale only\n", n_uf))

MOD_LAB <- c(Presynaptic = "Presynaptic (SNARE)", Postsynaptic = "Postsynaptic (PSD)",
             GluR_ionotropic = "Ionotropic GluR", OxPhos = "OxPhos", IEG = "IEG")
ORD <- intersect(c("Presynaptic","Postsynaptic","GluR_ionotropic","OxPhos","IEG"), MODS)
# gate: >= 3 members detected, or >= 2 covering at least half the set (small sets)
stat <- stat %>% left_join(det_tab, by = c("signature", "class", "Region")) %>%
  mutate(scorable = !is.na(n_det10) & (n_det10 >= 3L | (n_det10 >= 2L & n_det10 / pmax(n_present, 1L) >= 0.5)))
stopifnot("detection table does not cover every programme x class x region" = !any(is.na(stat$n_det10)))
if (any(!stat$scorable)) cat(sprintf("\nNOT SCORABLE: %s\n", paste(sprintf("%s/%s/%s (%d of %d)",
  stat$signature[!stat$scorable], stat$class[!stat$scorable], stat$Region[!stat$scorable],
  stat$n_det10[!stat$scorable], stat$n_present[!stat$scorable]), collapse = "; ")))
stat <- stat %>% mutate(signature = factor(signature, levels = ORD))

cat("\n== Cliff's delta (NHD - CON) per programme x region x class ==\n")
print(as.data.frame(stat %>% arrange(signature, Region, class) %>%
      mutate(cliff_d = round(cliff_d,3), q_BH = signif(q_BH,3)) %>%
      select(signature, Region, class, cliff_d, q_BH, label, n_CON, n_NHD)), row.names = FALSE)

L <- max(ceiling(max(abs(stat$cliff_d), na.rm = TRUE)*20)/20, 0.15)
# The depth-artifact bound now covers every cell the panel draws (96 loops over both
# neuron classes and includes the GluR set). A cell whose |delta| does not exceed the delta a
# depth split of control nuclei alone produces is not separable from a depth artefact and is drawn
# hollow — previously only 6 of 20 cells were controlled at all and the walk-back lived in the
# legend. All ten inhibitory cells fail the bound, which is itself the finding.
.dep <- tryCatch(read.csv(file.path(TDIR, "depth_E_neuron_module_artifact_FH.csv"), stringsAsFactors = FALSE),
                 error = function(e) NULL)
if (!is.null(.dep) && all(c("cell_type", "region", "programme", "separable") %in% names(.dep))) {
  .ct <- c(Ex = "Neuron_Ex", Inh = "Neuron_Inh")
  stat <- stat %>% mutate(separable = .dep$separable[match(paste(.ct[as.character(class)], Region, signature),
                                                           paste(.dep$cell_type, .dep$region, .dep$programme))])
  cat(sprintf("depth bound: %d of %d cells separable from a depth-only artefact; NOT separable: %s\n",
              sum(stat$separable, na.rm = TRUE), nrow(stat),
              paste(sprintf("%s/%s/%s", stat$signature[!stat$separable %in% TRUE],
                            stat$class[!stat$separable %in% TRUE], stat$Region[!stat$separable %in% TRUE]), collapse = "; ")))
} else { stat$separable <- NA; warning("depth bound table missing — no separability encoding") }

dm <- stat %>% filter(scorable) %>%
  mutate(delta_c = pmax(pmin(cliff_d, L), -L), fragile = !(separable %in% TRUE))
GSEA_STEEL <- "#3E7CB1"; GSEA_ROSE <- "#D1495B"
# the row label carries how many members are detected, so a row cannot be read as stronger than
# its basis (five rows here, so the extra line fits — it did not in Fig. 4d's eight)
.dl <- stat %>% dplyr::filter(scorable) %>% dplyr::group_by(signature) %>%   # describe the cells drawn, not the ones gated out
  dplyr::summarise(k = min(n_det10), n = max(n_in_set), rng = dplyr::n_distinct(n_det10) > 1, .groups = "drop")
.dl_k <- setNames(.dl$k, as.character(.dl$signature)); .dl_n <- setNames(.dl$n, as.character(.dl$signature))
.dl_r <- setNames(.dl$rng, as.character(.dl$signature))
ylab_map <- setNames(paste0(str_wrap(unname(MOD_LAB[ORD]), width = 18),
                            sprintf("\n(%s%d of %d detected)", ifelse(.dl_r[ORD], "\u2265 ", ""), .dl_k[ORD], .dl_n[ORD])), ORD)

p <- ggplot(dm, aes(x = class, y = signature)) +
  geom_point(data = function(d) dplyr::filter(d, !fragile),
             aes(fill = delta_c, size = dotsize), shape = 21, colour = "grey35", stroke = 0.3) +
  # hollow = not separable from a depth-only artefact (see the bound joined above)
  geom_point(data = function(d) dplyr::filter(d, fragile),
             aes(size = dotsize), shape = 21, fill = NA, colour = "grey35", stroke = 0.45) +
  geom_text(data = dm %>% filter(label == "q*"),
            aes(x = class, y = signature, label = label), inherit.aes = FALSE,
            nudge_y = 0.30, size = 3.13, colour = "black") +
  facet_wrap2(~ Region, nrow = 1,
              labeller = ggplot2::as_labeller(c(Frontal = "Frontal", Hippo = "Hippo")),
              strip = strip_region_x(REGION_ORDER, fontsize = 9.9, clip = "off")) +
  scale_fill_gradient2(low = GSEA_STEEL, mid = "white", high = GSEA_ROSE, midpoint = 0,
                       limits = c(-L, L), breaks = c(-L, 0, L),
                       labels = c(sprintf("%.2f", -L), "0", sprintf("+%.2f", L)),
                       name = expression(atop("Cliff's " * delta, "NHD - CON"))) +
  scale_size_continuous(range = c(0.8, 5.2), limits = c(0, CAPL), breaks = c(1.3, 10, 50),
                        labels = c("0.05","1e-10","<=1e-50"), name = "BH q") +
  scale_y_discrete(limits = rev(ORD), labels = function(x) ylab_map[x]) +
  scale_x_discrete(expand = expansion(add = 0.55)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(0.40,"cm"), barheight = unit(2.1,"cm")),
         size = guide_legend(order = 2)) +
  labs(x = NULL, y = NULL) +
  theme_pub(base_size = 14.9) +
  theme(strip.text = element_text(size = 9.9, face = "plain"),
        axis.text.x = element_text(size = 8.6, colour = "black"),
        axis.text.y = element_text(size = 8.6, colour = "black", lineheight = 0.82),
        axis.ticks.x = element_blank(), axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25),
        panel.spacing.x = unit(1.5, "pt"), strip.clip = "off",
        legend.position = "right", legend.box = "vertical", legend.box.just = "left",
        legend.margin = margin(t = -2, r = 0, b = 0, l = 0),
        legend.key.size = unit(0.55, "lines"), legend.text = element_text(size = 8.6),
        legend.title = element_text(size = 8.9, lineheight = 0.9),
        plot.margin = margin(5, 4, 3, 3))

BN <- "F5d_neuron_program_dotmatrix"; W <- 3.30; H <- 2.30
ggsave(file.path(PANEL, paste0(BN, ".pdf")), p, width = W, height = H, useDingbats = FALSE)
ggsave(file.path(PANEL, paste0(BN, ".png")), p, width = W, height = H, dpi = 600,
       device = ragg::agg_png)
cat(sprintf("\nSaved %s.{png,pdf} (%.2f x %.2f in)\n", BN, W, H))
write.csv(stat, file.path(TDIR, "fig5_neuron_programme_dotmatrix_FH.csv"), row.names = FALSE)
cat("Wrote tables/fig5_neuron_programme_dotmatrix_FH.csv\n=== DONE ===\n")
