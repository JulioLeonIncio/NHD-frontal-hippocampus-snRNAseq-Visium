#!/usr/bin/env Rscript
# =============================================================================
# 41b_build_microglia_mast_tables_FH.R — Drop-in MAST tables for the split myeloid classes, in the canonical schema the Figure-2 panels already expect.
# -----------------------------------------------------------------------------
# 40b wrote MAST_{Microglia,PVM}_{Frontal,Hippo}.csv, but the Fig-2 panel scripts
# read the combined table `MAST_dual_discovery_all.csv`, whose schema carries a
# leading `comparison` column and orders n_NHD before n_CON. Rather than editing
# column handling in three panel scripts (three chances to introduce a bug), this
# emits tables that are byte-compatible with the canonical schema, so each panel
# needs one change: the cell-type string.
#
# Produces:
#   MAST_dual_discovery_all_MICROGLIA.csv — full drop-in replacement for the
#       combined table: every non-myeloid cell type copied through unchanged, and
#       the pooled "Micro-PVM" rows REPLACED by the split Microglia + PVM rows.
#   MAST_Microglia_{Frontal,Hippo}.csv / MAST_PVM_*.csv — rewritten in canonical
#       column order (2I reads per-region files).
# The pooled MAST_dual_discovery_all.csv and MAST_Micro-PVM_*.csv are left intact
# as the supplementary sensitivity tier.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr) })

if (!nzchar(Sys.getenv("NHD_PROJ"))) stop("Set the NHD_PROJ environment variable to the NHD_frontal_hippo_rebuild folder (see README.md)")
PROJ <- dirname(Sys.getenv("NHD_PROJ"))
OUT <- file.path(PROJ, "NHD_frontal_hippo_rebuild", "tables", "mast_dual")
POOLED <- file.path(OUT, "MAST_dual_discovery_all.csv")
stopifnot(file.exists(POOLED))

pool <- read.csv(POOLED, stringsAsFactors = FALSE)
CANON <- names(pool)
cat("canonical schema:", paste(CANON, collapse = ","), "\n")

split_files <- file.path(OUT, sprintf("MAST_%s_%s.csv",
                                      rep(c("Microglia","PVM"), each = 2),
                                      rep(c("Frontal","Hippo"), 2)))
missing <- split_files[!file.exists(split_files)]
if (length(missing)) stop("MISSING split MAST outputs — run 40b first:\n  ",
                          paste(basename(missing), collapse = "\n  "))

split <- bind_rows(lapply(split_files, read.csv, stringsAsFactors = FALSE)) %>%
  mutate(comparison = paste0(cell_type, "_", region))

# Hard schema check before writing anything
need <- setdiff(CANON, names(split))
if (length(need)) stop("split tables lack canonical column(s): ", paste(need, collapse = ", "))
split <- split[, CANON, drop = FALSE]
stopifnot(identical(names(split), CANON))

# per-region files in canonical order (2I reads these)
for (ct in c("Microglia","PVM")) for (rg in c("Frontal","Hippo")) {
  d <- split[split$cell_type == ct & split$region == rg, , drop = FALSE]
  write.csv(d, file.path(OUT, sprintf("MAST_%s_%s.csv", ct, rg)), row.names = FALSE)
}

# combined drop-in: non-myeloid cell types unchanged + split myeloid rows
other <- pool[pool$cell_type != "Micro-PVM", , drop = FALSE]
combined <- bind_rows(other, split)
stopifnot(identical(names(combined), CANON))
write.csv(combined, file.path(OUT, "MAST_dual_discovery_all_MICROGLIA.csv"), row.names = FALSE)

cat("\n=== rows by cell_type in the drop-in combined table ===\n")
print(table(combined$cell_type))
cat("\n=== discovery counts (split myeloid) ===\n")
print(split %>% group_by(cell_type, region) %>%
        summarise(tested = n(), discovery = sum(discovery),
                  up = sum(discovery & cliffs_delta > 0),
                  down = sum(discovery & cliffs_delta < 0), .groups = "drop") %>%
        as.data.frame())
cat("\npooled Micro-PVM rows replaced:", sum(pool$cell_type == "Micro-PVM"),
    "-> split rows added:", nrow(split), "\n")
cat("\n=== DONE: 41b_build_microglia_mast_tables_FH ===\n")
