# Elastic Bounds: Novelty and Next Steps

**To:** Cees van der Eijk, Scott Moser
**From:** Scott Moser
**Date:** 15 October 2026 (revised 20 October 2026)

## Summary

Most of our theoretical results are already published. The paper's contribution therefore has to come from showing consequences in published applied work, and the most promising route is to re-analyse trend studies of issue alignment using exact bounds.

- **The core theory is in print.** Shojima (2024, *Behaviormetrika*) derives r_min and r_max for Likert items, notes that the two bounds are usually unequal in size, and proposes the same split rescaling we use. We don't cite it.
- **The bounds have an exact formula.** r_max and r_min follow from the Wasserstein distance between the two standardized marginals; we checked this on all 7,503 BES pairs. So the beta regression approximates a known formula rather than discovering what determines shrinkage.
- **Our Hypothesis Testing section needs a narrower claim.** Tight bounds don't create false positives. In small samples, though, they split the Type I error unevenly between the two tails and can make significance impossible in one direction. At BES sample sizes this doesn't matter.
- **The constraint argument is old.** That marginals distort comparisons of ideological constraint was argued in *AJPS* in 1978–79, and Fowler (2026, *Political Behavior*) recently restated it for two-category items.
- **What looks open:** applying exact bounds to the multi-category items in trend studies (Baldassarri & Gelman 2008; Kozlowski & Murphy 2021) to see how much of the measured change is mechanical. We propose a short ANES pilot before writing more.

## The beta-regression section

Cees's analysis is careful, and its descriptive numbers match our `bes_bounds.rds` exactly (C1 min .516, max .991, mean .885, SD .060). The catch is that C1 is an exact function of the two marginals.

For standardized variables, E[(Z_X − Z_Y)²] = 2 − 2r, and sorting both variables the same way makes that gap as small as possible. That smallest gap is the 2-Wasserstein distance W₂, which gives:

$$
r_{\max} = 1 - \tfrac{1}{2}W_2^2(Z_X, Z_Y), \qquad r_{\min} = -1 + \tfrac{1}{2}W_2^2(Z_X, -Z_Y), \qquad C_1 = 1 - \tfrac{1}{4}\left[W_2^2(Z_X, Z_Y) + W_2^2(Z_X, -Z_Y)\right]
$$

Here Z is the standardized marginal and −Z_Y is Y flipped end to end. On all 7,503 BES pairs these formulas match our computed bounds to about 10⁻¹³. Four consequences follow.

- **Mean and SD differences can't matter in themselves,** because r doesn't change when either variable is shifted or rescaled. That explains why the SD difference barely predicts C1 (r = .11) and why the mean-difference coefficient flips sign once skew and kurtosis are in the model.
- **C1 depends on each item's own skew, not only on differences.** r_min falls when both items lean the same way. Two identical items with p = (.6, .3, .1) have r_max = 1 but r_min = −.556, so C1 = .78. In BES, s₁² + s₂² correlates −.83 with C1, against −.63 for |s₁ − s₂| and −.67 for the kurtosis difference.
- **Which bound is tighter follows the product of the skews.** |r_min| − r_max correlates .94 with −s₁s₂: items skewed in opposite directions have |r_min| > r_max.
- **Significance tests don't apply.** The outcome is deterministic, and the pairs aren't independent: each of the 123 items appears in 122 pairs. The 7,503 pairs are also every pair from 123 items, not a sample of survey item pairs. So p-values, likelihood-ratio tests and BIC-based link selection aren't meaningful; dyad-robust standard errors (Aronow, Samii & Assenova 2015) would fix only the dependence. The section works better as description. Testing each pair's own r against zero is unaffected.

Smaller fixes to the draft:

- It says "BES 2024"; the data are BES 2019.
- It calls C1 "shrinkage s", but a higher C1 means *less* shrinkage; "range coverage" would be clearer.
- Footnote 1 treats r_min < 0 < r_max as a feature of the data. It holds for any two non-constant variables.
- Footnote 3 says C1 = 1 when the two items are identical. It also needs both distributions to be symmetric.
- φ is a precision parameter, Var(y|x) = μ(1 − μ)/(1 + φ), not "the variance of the conditional mean".
- In the memo table, kurtosis difference predicts C2 at .85, better than anything predicts C1. Is that a typo? The size and square/rectangular rows are still blank for C1 and C2.
- The June homework said "means and skew, not kurtosis"; the September memo says kurtosis is the strongest predictor.

## What in the current draft is already known

Every theoretical claim in `ElasticBounds-r_v6b.qmd` has a published precedent. Our own introduction already calls the problem "well-known in the theoretical statistical literature".

| Claim in our draft | Prior art | Status |
| --- | --- | --- |
| r_min and r_max from sorting the two variables the same or opposite way | Hoeffding 1940; Fréchet 1951; Whitt 1976; Shih & Huang 1992; Shojima 2024 for Likert items | Known |
| The two bounds usually differ in size | φ_max/φ_min for two-category items (Guilford 1965); Embrechts et al. 2002; Shojima 2024 | Known |
| A symmetric marginal gives \|r_min\| = r_max | Standard; follows from the φ formulas in the two-category case | Known |
| Rescaling r by its bounds (r\*, or ρ/ρ_max) | Cureton 1959 (φ/φ_max); Shih & Huang 1992; Shojima 2024 | Known. Goodman–Kruskal γ (1954) already reaches ±1 at the bounds |
| Pearson r on Likert items distorts factor analysis | Ferguson 1941 (difficulty factors); Flora & Curran 2004; Garrido et al. 2013 | Known |
| Uniform reference over all tables with the same margins (`02c`) | Diaconis & Efron 1985 | Known |
| A permutation test corrects the t-test | Not searched; Fisher's exact test has a similar small-sample issue | Partly right, in small samples only (below) |

**The Hypothesis Testing section needs a narrower claim.** The exact null distribution of r is not the t-distribution. With the marginals fixed, r can only fall in [r_min, r_max], and the null is skewed when the marginals are; its mean is exactly 0 and its SD exactly 1/√(n − 1).

We checked an extreme pair, p = (.85, .10, .04, .01) against its reverse (r_max ≈ .14–.15), using 40,000 permutations per row:

| n | Critical r (t-test, 5%) | Lower tail (nominal .025) | Upper tail (nominal .025) | Two-sided (nominal .05) |
| --- | --- | --- | --- | --- |
| 20 | .44 | .062 | .000 | .062 |
| 50 | .28 | .063 | .000 | .063 |
| 100 | .20 | .044 | .000 | .044 |
| 200 | .14 | .038 | .004 | .043 |
| 1000 | .06 | .029 | .019 | .048 |

Two-sided Type I error stays near .05, but two problems are real. When r_max is below the critical r, a positive correlation can never be significant, however strong the true association. Even at n = 1000 the tails are .029 and .019 rather than .025 each.

So the section's mechanism needs changing: tight bounds don't manufacture significance; they split the error unevenly and cap power in one direction. This matters when a bound is within a few standard errors of zero (roughly r_bound × √n below 2–4): small samples such as elite surveys, subgroups and pilots, with lopsided items. At BES sample sizes it doesn't, and Figure 15 remains an effect-size point.

**The Asymmetry section has two errors.** The Strict Asymmetry corollary attributes asymmetry to "a marginal sufficiently skewed". What matters is the relative *direction* of the two skews: two identically right-skewed items have r_max = 1 > |r_min|. And the claim that r_min "may not equal FH lower in discrete settings" is false: sorting one variable up and the other down always reaches the lower bound, for discrete and continuous marginals alike. The Fréchet–Hoeffding lower bound fails only with three or more variables.

**Asymmetry is not caused by discreteness.** Two identical lognormals have r_max = 1 but r_min = −.369, and two identical exponentials have r_min = 1 − π²/6 = −.645. Normal vs lognormal gives ±.764, equal because the normal is symmetric.

The driver is shape, mainly skew, and the Wasserstein formula holds for continuous variables too. Discreteness makes the problem nearly universal: with few categories most items are skewed or lumpy, and items with different numbers of categories can never share a shape.

## Literature check 1: shape and skew-direction results

None of these results is new. The closest paper, Shojima (2024), already covers most of our theory section for Likert items.

**[Shojima (2024), *Behaviormetrika*](https://link.springer.com/article/10.1007/s41237-024-00235-4)** treats finding r_max and r_min for R×C Likert items as a transportation problem. It:

- proves when the northwest-corner rule gives the optimum (equivalent to our sorting);
- notes that the bounds are "generally not symmetric";
- shows r_min = −1 exactly when one marginal is the reverse of the other;
- simulates how far the range shrinks for 5×5 and 5×7 tables;
- proposes ρ/ρ_max (ρ ≥ 0) or ρ/ρ_min (ρ < 0), our `split_warrens` adjustment.

It credits [Carroll (1945, 1961)](https://www.cambridge.org/core/journals/psychometrika/article/abs/nature-of-the-data-or-how-to-choose-a-correlation-coefficient/2AAE82C5FD725F889545B7B943465ACD) and Guilford (1965) with earlier maximal-correlation results for discrete items. What it lacks is a real-survey audit, the skew-direction rule, and any downstream application.

**The skew-direction rule is a textbook point.** McNemar's *Psychological Statistics* (1969, pp. 186–187) states it in words, according to a [secondary source](https://journal.sipsych.org/index.php/IJP/article/download/708/607). Marginals skewed in opposite directions compress the range near +1; marginals skewed the same way compress it near −1. [Wylie (1976, *EPM*)](https://consensus.app/papers/details/0a87af83e34451cea7e931f47589f781/) found r "severely limited" when marginals are "highly, oppositely skewed". For two-category items the rule follows exactly from the φ_max/φ_min formulas.

**The Wasserstein formula is a known decomposition.** The squared 2-Wasserstein distance splits into location, size and shape: (Δμ)² + (Δσ)² + 2σ₁σ₂(1 − ρ_QQ) ([waddR documentation](https://bioconductor.org/packages/release/bioc/vignettes/waddR/inst/doc/wasserstein_metric.html)). ρ_QQ, the correlation of the two quantile functions, is exactly r_max. This is usually credited to [Irpino & Romano and Verde & Irpino](https://arxiv.org/pdf/1110.1462), though we have that only secondhand.

**What we did not find:** the flipped-variable formula for the asymmetry; an explicit statement that one symmetric marginal guarantees |r_min| = r_max for *every* partner; and a check of the skew rule on real survey data (our .94). These are small corollaries and one empirical fact. They belong in the paper as a short proposition with these citations, not as its contribution.

## Literature check 2: constraint and issue alignment

The argument that marginals distort constraint comparisons is about 45 years old and was restated in 2026. What we found no one doing is applying exact bounds to the trend studies themselves.

**Political science had this debate in 1978–79.** After the 1964 ANES format change, Nie & Rabjohn (1979, *AJPS*) defended a rise in constraint. Their argument, as their critics summarise it, was that the pre-1964 items' marginals had "become skewed, thus attenuating the correlations among them" ([Sullivan et al. 1979](https://consensus.app/papers/details/dd3e080e2c5f578b89433842ad2c9e2c/)). Bishop et al. called this "extreme marginality… 'degenerate variance'" and disputed it ([1979 rejoinder](https://consensus.app/papers/details/024a891e0bda5cfe839afac2ebeb6bf0/)). [Sullivan, Piereson & Marcus (1978)](https://consensus.app/papers/details/6b27db40ae2f5631875f8f759dc2a44e/) and [Bishop et al. (1978)](https://consensus.app/papers/details/bdc34ff4aa1659428e561842b838f380/) traced the rise to the format change. All sides argued informally; nobody computed bounds.

**[Fowler (2026, *Political Behavior*, open access)](https://link.springer.com/article/10.1007/s11109-025-10111-w) is the modern version.** He shows that with two-category questions, even "maximally constrained" preferences can produce almost any r. Splitting off 1% on opposite sides gives about .01; 25% gives about .34. He also shows that elite-vs-mass comparisons can come from different distributions of ideology: .85 vs .27 under identical, maximal constraint.

He says the logic "would extend" to items with three or more options. His "maximal constraint" correlation is exactly our r_max. He cites neither the bounds literature nor the 1970s debate, and he analyses no trends.

**The trend studies haven't dealt with it.** [Baldassarri & Gelman (2008)](https://sites.stat.columbia.edu/gelman/surveys.course/BaldassarriGelman2008.pdf) model trends in Pearson r across 47 ANES items, 1972–2004. Their only comment on measurement is a footnote saying other correlation measures give similar results; they don't discuss marginals shifting over time.

We could not check how [Kozlowski & Murphy (2021)](https://consensus.app/papers/details/81b014eddec35a9bb4a3db7724386850/) measure correlation; they report rapid alignment growth from 2004 to 2016. The same goes for [DellaPosta (2020)](https://consensus.app/papers/details/7df727a81bc95b4fb23a76557fedd7c7/). [Hare (2021, *BJPS*)](https://consensus.app/papers/details/229b4b6b3a0c5e4e984e827e102a7dd9/) avoids the problem with a dynamic ordinal IRT model, so referees will ask what bounds add beyond IRT.

**The gap:** nobody we found computes exact bounds for the multi-category items in these studies. Nor has anyone asked how much of a reported rise in r simply follows rises in r_max over time.

## What might still be new

One empirical contribution looks open, and three smaller pieces could support it. Whether the paper fits Political Analysis or PSRM depends on what the first one finds.

1. **Re-analysing published findings with exact bounds (main candidate).** Trends in issue alignment are the strongest target. We would extend Fowler's argument from two-category to ordinal items, using r_max as the "maximal constraint" benchmark, and test whether published trends survive. If a meaningful share of a landmark trend is mechanical, that's a political-science finding; candidates are listed under Proposed next steps.
2. **A compact mechanism.** The Wasserstein formula and the skew-direction rule explain *why* bounds shrink and which side they shrink on. This goes in as a cited proposition, not as the contribution.
3. **An audit of a whole national survey.** None of the 7,503 BES pairs can reach the full [−1, 1] range; median C1 is about .90, and 5% of pairs fall below .76. This is real but descriptive.
4. **A narrower testing result.** In small samples with tight bounds, the t-test splits its Type I error unevenly between the tails and has no power in one direction. This is useful advice for applied work; whether it's new is unchecked.

Two referee questions to prepare for: why not polychoric correlations or ordinal IRT (Hare 2021)? And what "constraint" means here: observed r relative to its attainable range, or a latent correlation? If the pilot finds little, the fallback is a shorter methods note for a survey-methods journal, built around the R tool.

## Proposed next steps

We propose running a short pilot before writing anything more, then deciding together whether to go on.

Decisions for the two of us:

1. Do we agree that theoretical novelty can't carry the paper, and that Shojima (2024) should be cited prominently?
2. Should the beta-regression section be paused, or cut to a short descriptive appendix, until the pilot is done?
3. Should we leave the venue open (PA/PSRM vs a survey-methods journal) until the pilot is done?
4. What should the shrinkage section tell readers? If it's which item pairs to worry about, the formula answers that directly: each item's skew, whether the two lean the same way, and their numbers of categories. The R tool gives exact bounds anyway.

Replication candidates, from memory (data availability checked only for BES). The BES battery is the cheapest first test; Baldassarri & Gelman has the biggest potential payoff.

| Study | Data | What exact bounds would test |
| --- | --- | --- |
| Evans, Heath & Lalljee (1996): left-right and libertarian-authoritarian scales | Already in our BES pairs file (f01 items, v25–v35) | Whether skew direction predicts which dimension each item loads on |
| Baldassarri & Gelman (2008) | ANES cumulative file, 1972–2004 | Whether the trend in r survives the trend in r_max |
| Kozlowski & Murphy (2021) | ANES, 2004–2016 | The same, for the recent rise in alignment |
| DellaPosta (2020) | GSS belief networks | The same, for belief consolidation |
| Boutyline & Vaisey (2017) | ANES belief networks | Whether item centrality tracks item skew |
| Jennings (1992) | Convention delegates vs the public | Elite–mass constraint gaps when the two groups' marginals differ |
| Akkerman, Mudde & Zaslove (2014); Castanho Silva et al. (2020) | Populist-attitude scales | Dimensionality |
| Rosenberg self-esteem scale (Marsh 1996) | Various | Whether the positive/negative wording factors track skew direction |
| van der Eijk & Rose (2015), Risky Business | Simulated | Already planned |

Tasks:

- [ ] ANES pilot. Use the ANES cumulative file and the Baldassarri & Gelman item set. Compute r, r_min and r_max for each item pair in each year, and check how much of the trend in r tracks r_max (Scott).
- [ ] Cheapest pilot: factor-analyse the BES left-right/authoritarian battery (f01 items, v25–v35) and check whether skew direction predicts which dimension each item loads on.
- [ ] Read the methods of Kozlowski & Murphy (2021) and DellaPosta (2020): which correlation measure, and which items.
- [ ] Get full texts of Carroll (1961), McNemar (1969, pp. 186–187), Nie & Rabjohn (1979) and the Irpino & Romano decomposition.
- [ ] Narrow the Hypothesis Testing section to the small-sample claim, and fix the two errors in the Asymmetry section.
- [ ] Run a systematic literature search; this one was targeted.
- [ ] Meet after the pilot to decide whether to go ahead.

## References

We read Shojima, Fowler and the Baldassarri & Gelman methods section in full. The others below come from abstracts or secondary accounts. Classic works named only in the table (Hoeffding, Fréchet, Whitt, Embrechts et al., Cureton, Guilford, Goodman & Kruskal, Diaconis & Efron, Ferguson, Flora & Curran) are cited from memory or from our bib. So are Aronow, Samii & Assenova (2015) and the replication candidates under Proposed next steps.

- [Baldassarri & Gelman (2008), Partisans without Constraint, *AJS*](https://sites.stat.columbia.edu/gelman/surveys.course/BaldassarriGelman2008.pdf)
- [Bishop et al. (1978), Change in the Structure of American Political Attitudes, *AJPS*](https://consensus.app/papers/details/bdc34ff4aa1659428e561842b838f380/)
- [Bishop et al. (1979), Questions about Question Wording: A Rejoinder, *AJPS*](https://consensus.app/papers/details/024a891e0bda5cfe839afac2ebeb6bf0/)
- [Carroll (1961), The Nature of the Data, or How to Choose a Correlation Coefficient, *Psychometrika*](https://www.cambridge.org/core/journals/psychometrika/article/abs/nature-of-the-data-or-how-to-choose-a-correlation-coefficient/2AAE82C5FD725F889545B7B943465ACD)
- [DellaPosta (2020), Pluralistic Collapse, *ASR*](https://consensus.app/papers/details/7df727a81bc95b4fb23a76557fedd7c7/)
- [Fowler (2026), Reassessing Extremism, Polarization, and Constraint with Continuous Policy Questions, *Political Behavior*](https://link.springer.com/article/10.1007/s11109-025-10111-w)
- [Garrido et al. (2013), A New Look at Horn's Parallel Analysis with Ordinal Variables, *Psychological Methods*](https://consensus.app/papers/details/52feb0cef2775ddfa1cad827cfadd7e7/)
- [Hare (2021), Constrained Citizens?, *BJPS*](https://consensus.app/papers/details/229b4b6b3a0c5e4e984e827e102a7dd9/)
- [Irpino, Verde & de Carvalho (2011), on the Wasserstein location/size/shape decomposition](https://arxiv.org/pdf/1110.1462)
- [Kozlowski & Murphy (2021), Issue Alignment and Partisanship in the American Public, *Social Science Research*](https://consensus.app/papers/details/81b014eddec35a9bb4a3db7724386850/)
- [Leonov & Qaqish (2020), Correlated Endpoints: Simulation, Modeling, and Extreme Correlations, *Statistical Papers*](https://consensus.app/papers/details/9ef548020efd5c18a48535c62fbe4874/)
- [McNemar (1969), as summarised in a secondary source](https://journal.sipsych.org/index.php/IJP/article/download/708/607)
- [Shih & Huang (1992), Evaluating Correlation with Proper Bounds, *Biometrics*](https://consensus.app/papers/details/054f0e678f865a4d94fe45a9e0d9b668/)
- [Shojima (2024), Maximal Pearson's Correlation Coefficient between Likert Items, *Behaviormetrika*](https://link.springer.com/article/10.1007/s41237-024-00235-4)
- [Sullivan, Piereson & Marcus (1978), Ideological Constraint in the Mass Public, *AJPS*](https://consensus.app/papers/details/6b27db40ae2f5631875f8f759dc2a44e/)
- [Sullivan et al. (1979), The More Things Change, the More They Stay the Same, *AJPS*](https://consensus.app/papers/details/dd3e080e2c5f578b89433842ad2c9e2c/)
- [waddR documentation: 2-Wasserstein decomposition](https://bioconductor.org/packages/release/bioc/vignettes/waddR/inst/doc/wasserstein_metric.html)
- [Wylie (1976), Effects of Coarse Grouping and Skewed Marginal Distributions on the Pearson Correlation, *EPM*](https://consensus.app/papers/details/0a87af83e34451cea7e931f47589f781/)
