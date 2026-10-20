# Literature note: Shojima (2024/2025)

Shojima, K. (2025). Maximal Pearson's correlation coefficient between Likert items from the viewpoint of Hitchcock's transportation problem. *Behaviormetrika*, 52(2), 505–528. https://doi.org/10.1007/s41237-024-00235-4 (published online Aug 2024; open access, CC BY 4.0). Full text: `literature/Shojima - 2025 - …pdf`. Not yet in `paper/references.bib`.

Notes made 20 October 2026 (Scott, with Claude).

## What the paper does
- **Problem:** max/min Pearson r between two Likert items with fixed marginals, posed as a linear programme. This is Hitchcock's transportation problem. The objective is the cross-moment Σ w_r v_c p_rc, because the means and variances are fixed.
- **Method:** the northwest-corner rule (NWCR) gives r_max, and the southwest-corner rule gives r_min. He proves optimality using LP duality and complementary slackness.
- **Structure of the solution:** at most R + C − 1 non-empty cells, which form a spanning tree. He treats two degenerate cases: empty categories, and cumulative proportions of the two items that coincide. The NWCR reproduces Carroll's (1945, 1961) method.
- **Properties of the bounds:**
  - they are "generally not symmetric";
  - r_max = 1 only if R = C and the marginals are identical;
  - r_min = −1 only if one marginal is the reverse of the other.
- **Simulations:** 5×5 and 5×7 tables with binomial marginals. The range shrinks as the two means diverge. r_max has "ridges" (local maxima) where the two items' cut-points coincide.
- **Recommendations:**
  - Use polychorics. The polychoric of his maximal table is .990.
  - Or rescale with ρ* = ρ/ρ_max (ρ ≥ 0) or ρ/ρ_min (ρ < 0), following Cureton (1959). He leaves the properties of ρ* "for future study".

## What to include in our paper
1. **Cite him prominently as the closest prior work on the bounds themselves.** His NWCR and SWCR tables are exactly our same-order and opposite-order sort. Our code reproduces his examples exactly: Table 3 (3×5) gives r_max = .858 and r_min = −.845; §1.2 (2×3) gives .795 and −.767. Also cite Carroll (1945, 1961) and Guilford (1965) through him.
2. **Note the gap in his optimality proof (one sentence or a footnote).**
   - Appendix 5 establishes complementary slackness but never shows that the dual solution is feasible.
   - Dual feasibility is exactly where monotone category scores are needed. The NWCR is optimal for "Monge" cost matrices (Hoffman 1963).
   - Counterexample: scoring item A as 1, 3, 2 in his Table 3, the NWCR gives r = −.017, while the true maximum is .878.
   - Our rearrangement-inequality proof is complete and shorter.
3. **Position our contribution against his:**
   - Why the bounds shrink: shape, via the Wasserstein identity, not location or scale.
   - Which bound is tighter: the skew-direction rule, plus the condition that a symmetric marginal gives |r_min| = r_max.
   - Real data: the audit of all 7,503 BES pairs.
   - Simulations over random (Dirichlet) marginals rather than a one-parameter binomial family.
   - Inference and consequences.

   He has none of these.
4. **Reinterpret his "difference in means is a good indicator".** In a binomial family the mean fixes the skew, and r is unchanged by shifting a variable, so the mean only acts as a stand-in for shape. This is the same point as in our critique of the beta-regression section.
5. **Pick up his deferred question about ρ*.** His ρ* is exactly our `split_warrens` adjustment (`R/06_correlation_adjustment.R`). Our results on whether ρ* matrices stay valid correlation matrices (PSD):
   - 0% failures in real BES subsets of 8 and 12 items;
   - 0% in the mild synthetic setting;
   - about 10% in the adversarial synthetic setting (6–8 items, latent ρ .3–.85).
6. **Link his polychoric remark to Fowler (2026).** The maximal table corresponds to perfect latent association, so r_max is the Pearson r under perfect latent association. That is Fowler's "maximal constraint" correlation.
7. **Possible extension: the elasticity of the bounds.** See below.

## The dual: how the bounds move when marginals shift
**The dual problem:** minimise Σ_r α_r a_r + Σ_c β_c b_c subject to α_r + β_c ≥ w_r v_c, where a and b are the marginals. At the optimum it equals the maximal cross-moment. Only differences between the α's (or β's) are meaningful, because a constant can be shifted from α to β.

**Interpretation:** α_{r+1} − α_r is the change in the maximal cross-moment per unit of probability mass moved from category r to category r + 1 of item A.

**Closed form (derived 20 October 2026; checked against Shojima's dual values):**

    α_{r+1} − α_r = (w_{r+1} − w_r) · v*(A_r)

Here A_r is A's cumulative proportion up to category r, and v*(A_r) is item B's score at that quantile. In words, the shadow price of a cut-point is the other item's score at that cut-point. For r_min, use B reversed. Shojima's Table 3 duals are α = (−4, −3, 0) and β = (5, 7, 9, 12, 15). Their differences, (1, 3) and (2, 2, 3, 3), match the formula. These are the one-dimensional Kantorovich potentials from optimal transport, so the machinery is known. Applying it to bound elasticity would be new.

**Derivative of the bound** when mass moves from A_r to A_{r+1}, with Δw = w_{r+1} − w_r:

    d r_max = Δw · [ (v*(A_r) − m_B)/(σ_A σ_B) − r_max · (w_r + w_{r+1} − 2 m_A)/(2 σ_A²) ]

This matches finite differences. In Table 3, moving one percentage point of responses from A1 to A2 lowers r_max by .0073 and raises r_min by .0062.

**Kinks:** the bound surface is piecewise smooth. It has kinks where cut-points coincide, which is where the dual is non-unique. These kinks are local maxima of r_max: in an example, moving mass either way lowers r_max (by .0080 or .0026 per percentage point). This explains Shojima's "ridges".

**Use for us:** when marginals shift between survey waves or between groups, this decomposes the change in r_max (and r_min) into each category's contribution. That is the mechanical part of an apparent change in "constraint" or "alignment", which is the trend re-analysis in our novelty memo.
