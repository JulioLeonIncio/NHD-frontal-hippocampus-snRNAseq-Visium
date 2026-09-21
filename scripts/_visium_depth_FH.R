# =============================================================================
# _visium_depth_FH.R — Single source for the derived cortical-depth axis used by the Figure-6 panels. Sourced by 71 and 72 so they cannot drift apart.
# -----------------------------------------------------------------------------
# Why this replaces the first attempt.
# The depth axis first defined white matter as each section's top-15% cell2location
# OLIGODENDROCYTE tail. Rendered, the "white matter" came out as single spots sprayed
# across the whole section, so almost every cortical spot had a WM spot next door and
# the depth axis collapsed to ~0 everywhere. It was measuring deconvolution noise.
#
# The diagnostic that settled it — spatial coherence, each spot's value against the mean
# of its six nearest neighbours, per section:
#     gene MBP            0.65 0.51 0.37 0.52   mean 0.51
#     Structural myelin   0.54 0.44 0.51 0.43   mean 0.48
#     gene PLP1           0.47 0.35 0.45 0.30   mean 0.39
#     c2l Oligo           0.47 0.27 0.48 0.33   mean 0.39   <- the one that was used
#     gene MOBP           0.14 0.08 0.07 0.08   mean 0.09   <- near-empty, do not use
# Myelin transcript forms a contiguous compartment on these sections; the deconvolved
# oligodendrocyte proportion does not. So white matter is defined from the structural
# myelin score, NEIGHBOUR-SMOOTHED first so that a single bright spot in cortex cannot
# open a false WM seed, and then taken as that section's own upper tail.
#
# Everything is per section: a difference in section size, shape or orientation cannot
# masquerade as a difference in depth.
# Leon et al., Nasu-Hakola disease frontal cortex and hippocampus.
# =============================================================================

add_cortical_depth <- function(d, k = 6L, wm_q = 0.80,
                               score_col = "Structural myelin",
                               verbose = TRUE) {
  stopifnot(requireNamespace("FNN", quietly = TRUE))
  stopifnot(all(c("x", "y", "sample_id") %in% names(d)),
            score_col %in% names(d))
  d$wm_score <- NA_real_; d$is_wm <- NA; d$depth <- NA_real_

  for (sm in unique(as.character(d$sample_id))) {
    ii <- which(as.character(d$sample_id) == sm)
    xy <- cbind(d$x[ii], d$y[ii])
    v  <- d[[score_col]][ii]

    # 1. neighbour-smooth: mean of the spot and its k nearest neighbours
    nn <- FNN::get.knn(xy, k = k)$nn.index
    sm_v <- (v + rowSums(matrix(v[nn], nrow = length(ii)))) / (k + 1)
    d$wm_score[ii] <- sm_v

    # 2. white matter = this section's own upper tail of the SMOOTHED score
    thr <- as.numeric(quantile(sm_v, wm_q, na.rm = TRUE))
    wm  <- ii[sm_v >= thr]
    gm  <- setdiff(ii, wm)
    d$is_wm[ii] <- sm_v >= thr

    if (!length(wm) || !length(gm)) next
    # 3. depth = distance to the nearest WM spot, scaled by this section's own
    #    95th percentile so the axis is comparable across sections without being
    #    stretched by a single outlying spot.
    nnd <- FNN::get.knnx(cbind(d$x[wm], d$y[wm]), cbind(d$x[gm], d$y[gm]), k = 1)$nn.dist[, 1]
    d$depth[gm] <- nnd / as.numeric(quantile(nnd, 0.95, na.rm = TRUE))
    d$depth[wm] <- 0
    if (verbose)
      cat(sprintf("  %-14s WM %4d spots (smoothed %s >= %.3f) | GM %4d | median depth %.2f\n",
                  sm, length(wm), score_col, thr, length(gm),
                  median(d$depth[gm], na.rm = TRUE)))
  }
  d$depth <- pmin(d$depth, 1)      # the >95th-pct tail is flattened onto the pial end
  d
}
