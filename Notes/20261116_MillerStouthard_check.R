# 20261116_MillerStouthard_check.R
#
# What can be computed from the tables PRINTED in Miller & Stouthard (1975,
# European Journal of Political Research 3:219-258), without their microdata?
# Their Table V (p. 231, 1970 provincial vote) and Table XI (p. 242, 1972
# parliamentary vote) give, for Catholics with a homogeneous religious background,
# the vote by six ordered levels of confessional attachment (Always, Usually,
# Sometimes, Seldom, Never attend; No religious preference), with the number of
# cases in each level. Joint counts are rebuilt as (printed %) x (printed N);
# the percentages are rounded to whole numbers, so the counts are approximate.
# Dutch Reformed attachment counts: Table V (p. 231) plus Table II (p. 223) for
# the 120 with no religious preference (77+47+61+70+90+120 = 465, the Table IV N).
#
# Run from the repository root:
#   Rscript Notes/20261116_MillerStouthard_check.R > Notes/20261116_MillerStouthard_check_output.txt 2>&1

rmax_from_props <- function(sx, px, sy, py) {
  cx <- cumsum(px); cy <- cumsum(py); cx[length(cx)] <- 1; cy[length(cy)] <- 1
  br <- sort(unique(c(0, cx, cy))); mid <- (head(br, -1) + tail(br, -1)) / 2; w <- diff(br)
  qx <- sx[pmin(findInterval(mid, cx) + 1, length(sx))]; qy <- sy[pmin(findInterval(mid, cy) + 1, length(sy))]
  mx <- sum(px * sx); my <- sum(py * sy)
  (sum(w * qx * qy) - mx * my) / (sqrt(sum(px * sx^2) - mx^2) * sqrt(sum(py * sy^2) - my^2))
}
rmin_from_props <- function(sx, px, sy, py) -rmax_from_props(sx, px, -rev(sy), rev(py))

# Attachment scored 1 (no religious preference) ... 6 (always attend); vote 1 = confessional, 0 = anything else
score <- 1:6
summarise <- function(label, n, pconf) {
  n <- rev(n); pconf <- rev(pconf)            # printed order is Always ... No pref; reverse to 1..6
  px <- n / sum(n); share <- sum(px * pconf)
  mx <- sum(px * score); sx <- sqrt(sum(px * (score - mx)^2))
  r <- (sum(px * score * pconf) - mx * share) / (sx * sqrt(share * (1 - share)))
  lo <- rmin_from_props(score, px, 0:1, c(1 - share, share)); hi <- rmax_from_props(score, px, 0:1, c(1 - share, share))
  cat(sprintf("%-55s confessional share %.3f | r %.3f | r_min %.3f | r_max %.3f | r/r_max %.3f\n",
              label, share, r, lo, hi, r / hi))
  invisible(list(r = r, hi = hi, rstar = r / hi, share = share))
}

cat("Printed tables (counts rebuilt from rounded percentages):\n")
# Catholics, 1970 (Table V): N per level and % confessional vote
n70 <- c(324, 122, 81, 49, 63, 96);  c70 <- c(.52, .41, .19, .14, .13, .00)
# Catholics, 1972 (Table XI)
n72 <- c(259, 92, 59, 39, 44, 60);   c72 <- c(.61, .47, .29, .23, .14, .03)
a <- summarise("Catholics 1970 (Table V, provincial)", n70, c70)
b <- summarise("Catholics 1972 (Table XI, parliamentary)", n72, c72)

cat("\nAttachment distribution, Always ... No pref (shares):\n")
cat("  1970:", round(n70 / sum(n70), 3), "\n  1972:", round(n72 / sum(n72), 3), "\n")
cat(sprintf("  total-variation distance between the two: %.3f\n", 0.5 * sum(abs(n70 / sum(n70) - n72 / sum(n72)))))

d_r <- b$r - a$r
mech <- (a$rstar + b$rstar) / 2 * (b$hi - a$hi)
pos  <- (a$hi + b$hi) / 2 * (b$rstar - a$rstar)
cat(sprintf("\n1970 -> 1972 change in r: %.3f = bounds moving %.3f + position within bounds %.3f\n", d_r, mech, pos))

cat("\nCounterfactuals: keep each level's confessional rate, change the attachment distribution\n")
cat("(Miller & Stouthard's second type of change, pp. 227-229):\n")
summarise("Catholic 1970 rates, 1972 Catholic attachment", n72, c70)
nDR <- c(77, 47, 61, 70, 90, 120)
summarise("Catholic 1970 rates, 1970 Dutch Reformed attachment", nDR, c70)
