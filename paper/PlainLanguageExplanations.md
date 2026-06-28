# Plain-language explanations worth reusing later

A running collection of explanations from working sessions that landed well —
saved here so we don't have to reinvent them when we get to actually writing
the manuscript. Not manuscript prose yet; just a bank to draw from.

---

## The three "axes" of why two r_obs values aren't directly comparable

*(Scott's reaction, 2026-06-23: "makes perfect sense! ... This version of your
answer is really excellent.")*

I think the word "axes" was the problem — it sounds like one abstract
framework when it's really just three separate, very ordinary reasons why two
correlation numbers that look the same on paper might not mean the same
thing. Here's each one as a small story:

**Axis 1 — "how much room did this pair have to move?"** This is the
project's whole reason for existing. If 90% of people pick "strongly agree"
on item A, and only 10% pick the top category on item B, then even a
*perfect* underlying relationship between A and B can't produce r = 1,
because there just aren't enough "top-top" people to pair everyone up. The
ceiling is lower than 1 — and the floor might not be -1 either. So an r_obs
of .30 from a pair with a low ceiling (say, max possible is .40) is a much
bigger deal than an r_obs of .30 from a pair whose ceiling is .95. Same raw
number, very different amount of "stretch" used up. This is what the five
transforms (`U1`, `tight_linear`, `split_warrens`, `U3`) already fix.

**Axis 2 — "how many people did we measure?"** Even with identical marginals
and an identical true relationship, the r you compute from a sample jitters
around by pure chance, and that jitter is bigger with fewer respondents.
r = .30 from N = 50 is a much shakier estimate than r = .30 from N = 5,000.
This has nothing to do with bounds — it's just sampling noise. It's what
`07_permutation_zscore.R` already fixes via z-scores.

**Axis 3 — "how fine-grained was the measurement scale?"** A 2x2 table
(yes/no by yes/no) simply has less room to express a relationship than a 5x5
table does, independent of marginals or sample size. This one mostly matters
for sign-less measures like Cramer's V, which is part of why those are
disqualified as FA/SEM drop-ins — they throw away the sign entirely to solve
this problem, which is too high a price for our purposes.

The short version: same r_obs, three independent reasons it might not be
comparable to some other pair's r_obs — room to move, sample size, scale
coarseness. The project currently only has machinery for the first two.

---

## Percentile-in-the-feasible-polytope

*(Scott's reaction, 2026-06-23: "I still love your 'percentile-in-the-polytope'
idea! ... Please be sure to remember this for the future!")*

Instead of judging an observed r_obs only against the two FH-extreme corner
tables (r_min from the anti-comonotonic coupling, r_max from the comonotonic
one), build the *entire* distribution of r over every table consistent with
some conditioning scheme, and ask what percentile r_obs falls at within that
distribution. Two genuinely different versions of "every table" are in play,
and they ask different questions:

- **Margin-conditioned** ("Version B"): every table sharing this *specific*
  pair's actual observed row AND column marginals. This stays inside the
  project's core "elastic bounds depend on marginal shape" frame — the
  comparison population is exactly the set of tables this pair *could* have
  produced given its own marginals.
- **Size-conditioned, margin-free**: every table with the same J, K, and total
  N, with marginals left completely free to vary across the comparison
  tables. This is what `R/02a_mc_sampling_uniform_crosstab.R` already
  produces directly (stars-and-bars, exact, no MCMC needed) — Scott's stated
  preference as of 2026-06-23. Trade-off worth remembering: this version
  no longer corrects for *this pair's* marginal shape specifically, since most
  of the comparison tables won't share it — it answers "how unusual is r_obs
  among same-sized tables in general" rather than "...among tables this pair's
  marginals could actually produce."

Either version is a legitimate, well-defined statistic; they just answer
different questions, and which one is "the" percentile-in-the-polytope
depends on which comparison population you actually want.
