# Literature notes: categorised items, factor analysis, and φ/φ_max (19 November 2026)

Read for the question: is anything in our factor-analysis route new? Sources were read via NUsearch/LibKey (Scott's UoN login) or from Scott's Zotero.

## Bernstein & Teng (1989), *Psychological Bulletin* 105(3): 467–477

Read in full as Ovid HTML, which has no page numbers, so citations are to sections.

- **Design (Method).** 20 parallel items from one normal latent factor, n = 1,000, latent inter-item correlation .25, .50 or .75. Items categorised six ways:
  - continuous ("latent");
  - binary split at z = 0 (binary-.50);
  - binary split at z = +1 (binary-.16);
  - half split at +1 and half at −1 ("binary–H/L", i.e. opposite skew);
  - four categories with cut-offs at +.5, +1, +1.5 ("MC–equal");
  - half with those cut-offs and half mirrored ("MC–H/L").

  Dimensionality judged by principal-component eigenvalues (Kaiser), exploratory ML chi-square, and LISREL.
- **Abstract.** "All forms of categorization spuriously suggested multidimensionality." Varying item splits "further enhanced these spurious effects". Multicategory items "were more likely to yield artifacts than dichotomous items using inferential criteria".
- **Introduction.** The correlation between two items "will therefore be affected, in part, by the similarity of their distributions as well as by the similarity in their content". Difficulty factors (Carroll 1945; Ferguson 1941; Gorsuch 1983) "are essentially spurious".
- **Results, Principal Components.** "Differential categorization (binary–H/L and MC–H/L) produced spurious factors under all conditions."
- **Discussion.**
  - Unifactor data that have been differentially categorised "can be indistinguishable from" continuous two-factor data; "the difference would be noted in the univariate distributions".
  - Advice: "compare the univariate item statistics to the factor structure in order to determine the extent to which factors are defined statistically rather than substantively."
  - Remedies discussed: miniscales and extension analysis; confirmatory multiple-groups; IRT; Muthén's tetrachoric/polychoric approach. Collins et al. (1986) found frequent Heywood cases with tetrachorics.
- **Not discussed:** φ/φ_max or any rescaling by bounds.

**For us:** the spurious-factor mechanism, and even the informal diagnostic (compare item distributions with factors), are established. What Bernstein & Teng don't have is a formal diagnostic computed from the frequencies alone (the ceiling matrix R_max), or any evaluation of bound-based rescaling.

## Grønneberg & Foldnes (2024), *Psychological Methods* 29(1): 65–87

Zotero key `gron24-factor`; Scott's PDF. Online first 19 May 2022. Printed page = PDF page + 64.

- **p. 65 (abstract).** Ordinal factor analysis needs "substantive" knowledge of the *latent* item distributions, or of equally spaced thresholds; neither can be extracted from the data. An illustrative case turns a true two-factor model into a perfect one-factor fit. Proposes an adjusted polychoric estimator plus sensitivity analysis (R package `discnorm`, p. 67, n. 1).
- **p. 67.**
  - "It has long been known that Pearson correlations among ordinal variables tend to underestimate the underlying true response correlation."
  - Polychorics "may be severely biased" when the latent variables are not normal.
- **p. 74.** With seven symmetric categories, Pearson-based ML "mistakes the number of factors to be one instead of two". Paradoxically, it is approximately unbiased only with *asymmetric* observed distributions, when there are many equally spaced latent thresholds.
- **p. 76.** With a normal copula but non-normal latent marginals, the polychoric estimates "the maximum correlation", i.e. Gebelein–Rényi maximal correlation over all transformations. "As far as we know, it has not been studied in the psychometric literature." This is a different object from our r_max, the largest Pearson r compatible with fixed observed marginals.

**For us:** they show the mirror-image danger. Marginal effects can also *hide* real dimensions, and the polychoric remedy rests on latent marginals that can't be observed. Our simulation generated normal latent variables, which favours polychorics, so it should be rerun with non-normal latents before drawing conclusions about polychorics.

## Davenport & El-Sanhurry (1991), "Phi/Phimax: Review and synthesis", *Educational and Psychological Measurement* 51(4): 821–828

DOI 10.1177/001316449105100403. **UoN has no access** (SAGE "Restricted access"); abstract and references only.

- **Abstract:** φ/φ_max "is an asymmetrical, equal-interval step function with ties to probability" and "is shown to be nonrobust in situations where it should be most useful". The φ_min adjustment for inverse relations is "also considered".
- **References include:**
  - Carroll (1961);
  - Collins, Cliff, McCormick & Zatkin (1986), "Factor recovery in binary data sets: A simulation", *Multivariate Behavioral Research* 21: 377–391;
  - Cureton (1959);
  - Ferguson (1941);
  - Guilford (1965), "The minimal phi coefficient and the maximal phi", *EPM* 25: 3–8.
- **Still unknown:** whether they discuss matrices of φ/φ_max values, or factor analysis of them. That needs the full text (interlibrary loan).

## Carroll (1961), *Psychometrika* 26: 347–372

Scott's copy in ~/Documents.

- **pp. 356–357.** Constraints from marginals "operate to alter the rank of a correlation matrix". He credits Ferguson with showing that φ among items differing in difficulty precludes rank one, "also true whenever Pearsonian correlations are based on disparate marginal distributions".
- **p. 369.** Guilford's "difficulty factors" are "artifacts of the kinds of correlations used". Recommends tetrachoric correlations.

## Open question for novelty

No source found so far reports that matrices of bound-rescaled correlations (φ/φ_max and its ordinal extensions) are often not valid correlation matrices when items are skewed in opposite directions. A brief web search found nothing.

Next to read: Collins et al. (1986), Cureton (1959), Guilford (1965), and the full Davenport & El-Sanhurry.
