From Contingency Tables to Correlation: Bounding Pearson’s \textit{r} for Ordinal Data with Fixed Marginals

# Introduction and Motivation

- Pearson’s correlation coefficient $r$ is widely used in social science
  research, including with ordinal variables such as Likert-type survey
  items.

- A common but problematic assumption is that ordinal variables can be
  treated as interval-scaled, allowing use of parametric statistics like
  Pearson’s $r$.

- We explore how marginal distributions constrain the possible values of
  $r$, leading to an *elastic* range bounded by a theoretical minimum
  and maximum.

- Our contributions are fourfold:

  1.  Derive analytic and algorithmic bounds for $r$ with fixed
      marginals.

  2.  Identify conditions under which bounds are symmetric or
      asymmetric.

  3.  Illustrate the practical consequences for applied statistical
      modeling.

  4.  Provide a software toolkit to support computation and testing.

# Background and Definitions

- **Ordinal Variables**: Discrete, ordered categories, often arising
  from rating scales. Their use in correlation poses challenges due to
  undefined distance between levels.

- **Common Correlation Alternatives**:

  - Spearman’s $\rho$: rank-based, but ignores marginal distributions.

  - Polychoric correlation: assumes underlying bivariate normality.

  - Polyserial correlation: mixes continuous and ordinal assumptions.

- **Coupling and Rearrangement**: A probabilistic method for generating
  joint distributions with fixed marginals; foundational for bounding
  statistics like $r$.

# Theoretical Framework

## Optimal Coupling and Correlation Bounds

We formalize the problem of identifying the joint distribution
(coupling) of two ordinal variables with fixed marginals that maximizes
or minimizes the Pearson correlation.

### Corollary: Analytic Forms for $r_{\min}$ and $r_{\max}$ [corollary-analytic-forms-for-r_min-and-r_max]

Closed-form expressions and bounds are derived from known results on
extremal expectations under fixed marginals. In some cases, extrema are
not strictly attainable due to discretization effects.

## Fréchet–Hoeffding and Boole–Fréchet Inequalities

These classical results define the upper and lower limits of joint
cumulative distributions given fixed marginals. When applied to ordinal
variables, they constrain the admissible space for $E[XY]$, and thus
$r$.

## Symmetry and the Role of Ties

- We explore when $r_{\min} = -r_{\max}$, typically under symmetric
  marginals.

- Ties in distributions (non-uniformity or repeated values) often lead
  to asymmetry in the feasible bounds of $r$.

## Special Cases and Boundary Conditions

- Our framework accommodates both equal (e.g., 4×4) and unequal (e.g.,
  4×5, 3×7) numbers of categories.

- We demonstrate how asymmetric marginals (e.g., skewed, unimodal)
  affect the feasible range of $r$.

## Relation to Copula Theory (Optional)

Though we work in discrete space, the problem of bounding correlation
given fixed marginals is conceptually related to copula theory, which
models dependence structures independent of marginal forms.

# Empirical Illustration and Practical Relevance

## Simulation Design

- We generate random contingency tables via permutation-based
  randomization, preserving marginal distributions.

- Scenarios include:

  - Category structures: 3×3, 4×4, 4×5, 3×7, 5×5, 7×10

  - Marginal types: uniform, skewed, bimodal

  - Distributions with varying entropy and variance

- For each scenario, we simulate the empirical distribution of $r$ and
  compare it to theoretical min/max values.

## Effects of Marginal Shapes

- Skewness and modality heavily influence the attainable range of $r$.

- Marginals with extreme concentration (e.g., most responses in one or
  two categories) sharply restrict the possible correlation range.

## Symmetry Breaking in Correlation Bounds

- In asymmetric tables, the min and max bounds of $r$ are not symmetric
  about zero.

- The shape of the null distribution of $r$ can be strongly skewed,
  especially in cases with unequal marginals.

## Applications and Implications

- Impacts standard methods:

  - **PCA and Factor Analysis**: Attenuated correlations may lead to
    over-factoring.

  - **SEM**: Biased estimation of structural paths if ordinal-level
    attenuation is ignored.

- We recommend:

  - Reporting observed $r$ in context of its attainable bounds.

  - Using our tools to test whether an observed $r$ is unusually strong
    or weak given marginals.

# Extensions and Open Questions

- How does marginal entropy relate to the flexibility of $r$’s bounds?

- Is there a principled way to adjust or normalize $r$ given its
  feasible interval?

- Can these results extend to multivariate correlation structures (e.g.,
  partial $r$, canonical correlation)?

- What is the potential for integrating with discrete copula frameworks?

# Conclusion

- Pearson’s $r$, when used with ordinal variables, is bounded in
  non-obvious ways by the marginal distributions.

- We provide theoretical derivations, empirical validation, and software
  tools to make these constraints visible and usable.

- Future work should explore latent modeling connections, correction
  strategies, and extensions to higher-dimensional data.

# Appendix A: Mathematical Proofs

- Formal derivations for optimal coupling.

- Closed-form expressions for max/min $E[XY]$ and $r$.

# Appendix B: Software Tools and Code Listings

- R functions to compute min/max correlation bounds.

- Permutation-based testing for fixed marginals.

- Usage examples and benchmarks.

# Appendix C: Additional Tables and Simulations

- Supplementary figures: empirical null distributions.

- Case studies: 4×5 and 3×7 examples.

[^1]: Corresponding author. <a href="scott.moser@nottingham.ac.uk"
    class="uri">scott.moser@nottingham.ac.uk</a>

