# =============================================================================
# _visium_domains_FH.R — Single source for the Visium spatial-domain vocabulary.
# Sourced by 73, 74, 108b, 111, 114 (and mirrored by Visium/figure_Visium/
# assemble_figure_Visium.R, which reads the same two files) so that the cluster ->
# domain map, the domain order, the palette, the grey/white compartment and the
# laminar depth axis are read from disk and never typed into a script.
#
# Files (Visium/integrated_harmony/):
#   integrated_cluster_identity.csv   cluster (g0..g9) -> identity  (curated by hand)
#   integrated_domain_levels.csv      domain, order, colour, compartment, depth_rank
#
# Exports (all in superficial -> deep order as given by `order`):
#   DOM_INFO   data.frame  the levels file, ordered
#   DOM_MAP    named chr   cluster number ("0".."9") -> domain
#   DOM_LEV    chr         every domain, in `order`
#   DOM_PAL    named chr   domain -> hex colour
#   DOM_COMP   named chr   domain -> compartment (grey | white | other)
#   DOM_GREY   chr         grey-matter domains in depth_rank order (the laminar axis)
#   DOM_WHITE  chr         white-matter domain(s)
#   DOM_BANDS  chr         DOM_GREY then DOM_WHITE (the Fig-6d column order)
# Requires PROJ to be defined by the caller (the portable project-root resolver).
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================
stopifnot(exists("PROJ"), dir.exists(PROJ))
.IH <- file.path(dirname(PROJ), "Visium", "integrated_harmony")
.idf <- file.path(.IH, "integrated_cluster_identity.csv")
.lvf <- file.path(.IH, "integrated_domain_levels.csv")
stopifnot("MISSING integrated_cluster_identity.csv" = file.exists(.idf),
          "MISSING integrated_domain_levels.csv"    = file.exists(.lvf))

DOM_INFO <- read.csv(.lvf, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("domain","order","colour","compartment","depth_rank") %in% names(DOM_INFO)),
          !anyDuplicated(DOM_INFO$domain), !anyDuplicated(DOM_INFO$order),
          all(DOM_INFO$compartment %in% c("grey","white","other")),
          all(grepl("^#[0-9A-Fa-f]{6}$", DOM_INFO$colour)))
DOM_INFO <- DOM_INFO[order(DOM_INFO$order), ]
DOM_LEV  <- DOM_INFO$domain
DOM_PAL  <- setNames(DOM_INFO$colour, DOM_INFO$domain)
DOM_COMP <- setNames(DOM_INFO$compartment, DOM_INFO$domain)
.g <- DOM_INFO[DOM_INFO$compartment == "grey", ]
stopifnot("every grey-matter domain needs a depth_rank" = !anyNA(.g$depth_rank),
          "depth_rank must be unique within grey matter"  = !anyDuplicated(.g$depth_rank))
DOM_GREY  <- .g$domain[order(.g$depth_rank)]
DOM_WHITE <- DOM_INFO$domain[DOM_INFO$compartment == "white"]
DOM_BANDS <- c(DOM_GREY, DOM_WHITE)

.id <- read.csv(.idf, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("cluster","identity") %in% names(.id)))
DOM_MAP <- setNames(.id$identity, sub("^g", "", .id$cluster))
.unk <- setdiff(unique(DOM_MAP), DOM_LEV)
if (length(.unk)) stop("integrated_cluster_identity.csv names a domain absent from integrated_domain_levels.csv: ",
                       paste(.unk, collapse = ", "))
rm(.IH, .idf, .lvf, .g, .id, .unk)
