# =============================================================================
# Uniform sampling of J x K cross-tab matrices with fixed total N
# (row/column margins are NOT fixed -- only the grand total N is fixed)
# =============================================================================
#
# THE KEY POINT:
# The set of J*K-cell tables whose entries sum to N is in 1-1 correspondence
# with compositions of N into m = J*K nonnegative integer parts. There are
# choose(N + m - 1, m - 1) such tables. "Sampling uniformly" means every one
# of these tables gets equal probability.
#
# This is NOT what you get by drawing cell counts i.i.d. multinomial(N, 1/m).
# Multinomial sampling overweights "balanced" tables (there are combinatorially
# more compositions near the mean), so it is biased relative to true uniform
# sampling over the discrete table space.
#
# THE METHOD (stars and bars):
# Lay out N "stars" (units) and m-1 "bars" (dividers) in a row of N+m-1 slots.
# Choosing which m-1 of the N+m-1 slots are bars, uniformly at random, induces
# a uniform distribution over compositions of N into m parts -- i.e. over all
# tables. The part sizes (= cell counts) are read off as the gaps between
# consecutive bars.
# =============================================================================
#
# See also: R/02b_mc_sampling_crosstab_given_r.R, which builds on this module
# to sample tables conditioned on a TARGET Pearson's r (not just uniform over
# all tables of a given size).
# =============================================================================

sample_uniform_crosstab <- function(N, J, K) {
  m <- J * K
  total_slots <- N + m - 1
  bars <- sort(sample.int(total_slots, m - 1))
  parts <- diff(c(0L, bars, total_slots + 1L)) - 1L
  matrix(parts, nrow = J, ncol = K)
}

# Draw many tables at once
sample_many_crosstabs <- function(n_draws, N, J, K) {
  replicate(n_draws, sample_uniform_crosstab(N, J, K), simplify = FALSE)
}

# =============================================================================
# Sanity check: the marginal distribution of any single cell count x_1 has a
# known closed form under true uniform sampling:
#
#   P(x_1 = k) = choose(N - k + m - 2, m - 2) / choose(N + m - 1, m - 1)
#
# Compare empirical frequencies from the sampler against this formula.
# =============================================================================

N <- 20; J <- 2; K <- 2; m <- J * K

draws <- sample_many_crosstabs(20000, N, J, K)
cell11 <- sapply(draws, function(tab) tab[1, 1])

emp  <- as.numeric(table(factor(cell11, levels = 0:N)) / length(cell11))
theo <- choose(N - (0:N) + m - 2, m - 2) / choose(N + m - 1, m - 1)

cat("k, empirical, theoretical\n")
print(round(cbind(k = 0:N, empirical = emp, theoretical = theo), 4))

# Example: draw one random 3x4 table with N = 50
print(sample_uniform_crosstab(N = 50, J = 3, K = 4))
