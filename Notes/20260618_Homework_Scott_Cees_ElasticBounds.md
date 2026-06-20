

### **Meeting Summary: Homework from the June 18, 2026 Working Session (Scott & Cees)**


---

#### **Agreed Manuscript Structure**

1. Introduction
2. Literature Review (brief "Problem 1" survey, pivoting to "Problem 2")
3. Algorithm / Theoretical Results — "metrics of constraint" (already in good shape)
4. Monte Carlo Validation — **Scott**
5. Empirical Illustration / BES 2019 — **Cees**
6. Implications for Applied Research — **Scott**, including two illustrative examples

---

#### **Scott's Homework**

1. **Finish/formalize the Monte Carlo simulation design ("types of randomness").** Vary the number of categories from 3–11, including both square (K×K) and rectangular (K×J) variable pairs; vary sample size; and broaden the marginal shapes tested beyond normal/uniform to include bipolar, tripolar, and embedded-zero distributions. Implement the bivariate K×J cell-allocation-with-replacement algorithm discussed in the meeting.
2. **Code and run the simulations comparing R-adjustment methods.** Compare: the stretch-to-[-1,1] "S-star" method, divide-by-range, and the two-step "2x" rescaling approach. Also explore a Kendall's-tau-based alternative. The chi-squared approach was explicitly discussed and rejected — do not pursue it.
3. **Build two illustrative examples** for the Implications section:
   - A 0.6-vs-0.5 rank-reversal example.
   - A factor-analysis dimensionality-change example, reusing/adapting the contrived item matrices from Cees & Jonathan Rose's "Buyer Beware: Risky Business" paper.
4. **Write the "Implications for Applied Research" manuscript section** — covering the hypothesis-testing bias toward non-significance, the permutation-test alternative, and proposed remedies.
5. **Lower priority / shelved:** keep percolating on the copula-based "C function" generalized-dependence idea, but do not actively pursue it right now.
6. **Finish verifying the remaining Consensus-AI-sourced literature references** (~20% still unchecked against original sources).

#### **Cees's Homework**

1. **Write up the BES 2019 empirical results into manuscript prose.** The analysis itself is done (range width is explained mainly by the difference in means and difference in skew between the two marginals, not by kurtosis) but currently exists only as bullet points/appendix material — needs to become full prose for the Empirical Illustration section.
2. **Write the Introduction and Literature Review**, including a concise (roughly half-to-two-thirds of a page) "Problem 1" treatment that tables the existing debate (citing O'Brien 1979, *American Sociological Review*, vol. 44 no. 5; Grether 1976, *American Sociological Review*, vol. 41 no. 5; and Norman 2010, *Advances in Health Sciences Education*, vol. 15) before pivoting to "Problem 2," the paper's actual contribution.
3. **Track down the citation** for the prior author (per Cees's own notes / a previously-sent email) who already does asymmetric positive/negative stretching of correlations against known bounds — this needs to be cited and distinguished from the current paper's approach.
4. **Draft the footnote/justification for restricting the paper's scope to ordered categorical (not continuous) variables**, using the argument that partial-identification intervals converge to a point as the number of categories and/or sample size grows. Cite Moss & Grønneberg (2023, *Psychometrika*, vol. 88) — specifically the figure on p. 245 and Theorem 1 on p. 244.

---

#### **Joint / Open / Shelved Items**

- **Copula-based "C function" generalized-dependence idea** — explicitly shelved for now ("time to percolate"); not assigned to either person as active work.
- **Whether to seek a real (vs. contrived) published empirical example** for the dimensionality-change illustration — discussed and decided against; defaulting to the contrived/reused item matrices from the Buyer Beware paper.
- **Whether the critique generalizes to all SEM-with-Likert-data work in political science** — raised by Scott, judged plausible by Cees, but left unresolved and unassigned. Worth revisiting once the core results are written up.

---


#### **Key Takeaways**

* Scott owns the Monte Carlo validation, the R-adjustment-method comparison, the two illustrative examples, and the Implications section.
* Cees owns the BES 2019 write-up, the Introduction/Literature Review (including the Problem 1 framing and the ordinal-scope footnote), and tracking down the asymmetric-stretching precedent citation.
* The copula "C function" idea and the broader political-science-SEM generalization question are both open but not currently active — park them for a future discussion.
