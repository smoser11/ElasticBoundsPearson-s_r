# 20261118_bounds_shape_examples.R
#
# When do the bounds of r move? Exact bounds for small three-category examples,
# to check the condition "for the floor/ceiling to change, the marginals must become
# more (less) asymmetric". Scores 1, 2, 3 for both items.
#
# Run from the repository root:
#   Rscript Notes/20261118_bounds_shape_examples.R > Notes/20261118_bounds_shape_examples_output.txt 2>&1

rmax_from_props <- function(sx, px, sy, py) {
  cx <- cumsum(px); cy <- cumsum(py); cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy))); mid <- (head(br, -1) + tail(br, -1)) / 2; w <- diff(br)
  qx <- sx[pmin(findInterval(mid, cx) + 1, length(sx))]; qy <- sy[pmin(findInterval(mid, cy) + 1, length(sy))]
  mx <- sum(px * sx); my <- sum(py * sy)
  (sum(w * qx * qy) - mx * my) / (sqrt(sum(px * sx^2) - mx^2) * sqrt(sum(py * sy^2) - my^2))
}
rmin_from_props <- function(sx, px, sy, py) -rmax_from_props(sx, px, -rev(sy), rev(py))
show <- function(label, px, py) {
  s <- 1:3
  cat(sprintf("%-58s X = (%s)  Y = (%s)  r_min %.3f  r_max %.3f\n", label,
              paste(px, collapse = ","), paste(py, collapse = ","),
              rmin_from_props(s, px, s, py), rmax_from_props(s, px, s, py)))
}
show("1. Both symmetric, same shape",                         c(.25, .5, .25), c(.25, .5, .25))
show("2. Both equally skewed, same direction",                c(.6, .3, .1),   c(.6, .3, .1))
show("3. Both more skewed, same direction",                   c(.8, .15, .05), c(.8, .15, .05))
show("4. One symmetric, one skewed",                          c(.25, .5, .25), c(.6, .3, .1))
show("5. Skewed in opposite directions",                      c(.6, .3, .1),   c(.1, .3, .6))
show("6. Both symmetric, different peakedness (no skew)",     c(1/3, 1/3, 1/3), c(.1, .8, .1))
