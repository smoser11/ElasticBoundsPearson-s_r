---
title: "Extremal Pearson Correlation for Ordinal Variables: Consolidated Analysis"
format: pdf
toc: true
number-sections: true
---

# 1. Problem Setting and Notation

We consider two ordinal variables:

- \( X \) with \(K_X\) ordered categories \(0, 1, \dots, K_X-1\)
- \( Y \) with \(K_Y\) ordered categories \(0, 1, \dots, K_Y-1\)

**Marginals** given by:

\[
p_X(i) = P(X = i), \quad p_Y(j) = P(Y = j)
\]
or counts:
\[
n_{X=i}, \quad n_{Y=j}, \quad N = \sum_i n_{X=i} = \sum_j n_{Y=j}
\]

**CDFs**:
\[
F_X(i) = \sum_{k \le i} p_X(k), \quad F_Y(j) = \sum_{k \le j} p_Y(k)
\]

**Pearson correlation**:
\[
r_{XY} = \frac{\mathrm{Cov}(X,Y)}{\sigma_X \sigma_Y}, \quad
E[XY] = \sum_{i=0}^{K_X-1} \sum_{j=0}^{K_Y-1} x_i y_j \, P(X=x_i, Y=y_j)
\]
For Spearman's correlation, replace \(X,Y\) by their ranks \(R_X,R_Y\).

*(Sources: Proofs-OrdinalVariables.md; Ordinal Variables and Correlation.md; BFbounds.md; Math Proof GPT – Whitt 1976 Paper.md)*

---

# 2. Fréchet–Hoeffding / Boole–Fréchet Bounds

Any joint CDF \(F_{XY}\) with marginals \(F_X, F_Y\) satisfies:

\[
\max\{F_X(i) + F_Y(j) - 1, 0\} \ \le \ F_{XY}(i,j) \ \le \ \min\{F_X(i), F_Y(j)\}
\]  
*(BFbounds.md)*

- **Upper bound (comonotonic)**:
\[
F^{\text{upper}}_{XY}(i,j) = \min\{F_X(i), F_Y(j)\}
\]
- **Lower bound (countermonotonic)**:
\[
F^{\text{lower}}_{XY}(i,j) = \max\{F_X(i) + F_Y(j) - 1, 0\}
\]

---

# 3. Maximizing and Minimizing \(r\)

## 3.1 Max \(r\) — Comonotonic Coupling

**Claim:** The maximum correlation \(r_{\max}\) is achieved by the comonotonic (FH upper) arrangement.

**Proofs / Arguments:**
- **FH upper pmf construction**: In discrete case, differences of CDF steps allocate mass along the diagonal, pairing high \(X\) with high \(Y\).  
  *(BFbounds.md; Proofs-OrdinalVariables.md)*
- **Rearrangement inequality (ranks)**: For ranks \(R_X, R_Y\), \(\sum R_X R_Y\) is maximized when sequences are ordered the same way.  
  *(Ordinal Variables and Correlation.md)*
- **Quantile coupling**: \(U \sim \mathrm{Unif}(0,1)\), \(X = F_X^{-1}(U), Y = F_Y^{-1}(U)\).  
  *(Math Proof GPT – Whitt 1976 Paper.md)*

---

## 3.2 Min \(r\) — Countermonotonic Coupling

**Claim:** The minimal correlation \(r_{\min}\) corresponds conceptually to the countermonotonic (FH lower) arrangement.

**Proofs / Arguments:**
- **FH lower pmf construction**: Differences of the lower-bound CDF.  
  *(BFbounds.md)*
- **Opposite-order pairing**: Sort \(X\) descending, \(Y\) ascending; achieves best feasible negative correlation.  
  *(Proofs-OrdinalVariables.md; Ordinal Variables and Correlation.md)*
- **Quantile coupling**: \(X = F_X^{-1}(U), Y = F_Y^{-1}(1-U)\) in ideal continuous case.  
  *(Math Proof GPT – Whitt 1976 Paper.md)*

**Note:** In discrete settings, FH lower may not be attainable; see Section 5.

---

# 4. Binary (2×2) Cases and Examples

*(From Cheemera v1.1 – Pearson Correlation Analysis Implications.md)*

## 4.1 Balanced marginals (70–30, 70–30)
- Perfect diagonal fill yields \(r_{\max} = 1\).

## 4.2 Skew vs balanced (90–10, 50–50)
- Perfect alignment impossible; \(r_{\max} < 1\) due to capacity limits.

---

# 5. Attainability and Symmetry (Updated)

## 5.1 Maximal correlation
- **Result:** In the discrete ordinal setting, \(r_{\max}\) is **always** attained by the FH upper (comonotonic) construction.  
  *(BFbounds.md; Proofs-OrdinalVariables.md; Ordinal Variables and Correlation.md)*

## 5.2 Minimal correlation
- **Result:** The FH lower (countermonotonic) construction may be **unattainable** in discrete settings; \(r_{\min}\) is then the best feasible opposite-order arrangement.  
  *(BFbounds.md; Proofs-OrdinalVariables.md; Ordinal Variables and Correlation.md)*

## 5.3 Magnitude relationships
- **Earlier view:** Often \(|r_{\min}| < r_{\max}\) for discrete asymmetrical marginals.  
  *(BFbounds.md; Ordinal Variables and Correlation.md)*
- **New evidence (BES data)**: Possible to have \(|r_{\min}| > r_{\max}|\).  
  Example: \(X\) with 4 categories, \(Y\) with 5 categories, skewed marginals;  
  \(r_{\max} = 0.85027\), \(r_{\min} = -0.89916\).  
  *(Cees – r_min MAGNITUDE Asymmetry in Pearson Correlation.md)*

**Reason:** Pearson correlation uses centered values; skewness can make negative deviations in antitone pairing larger in magnitude than positive deviations in comonotone pairing.

## 5.4 Symmetry condition
- If one marginal is **palindromically symmetric** about its center rank, then \(r_{\min} = -r_{\max}\).  
  *(Cees – r_min MAGNITUDE Asymmetry in Pearson Correlation.md)*

---

# 6. Simulation / Randomization

**Permutation distribution:**  
Generate full vectors from marginals, independently permute \(X\) and \(Y\), compute \(r\), repeat. Produces null distribution centered at 0.  
*(Proofs-OrdinalVariables.md)*

---

# 7. Key Tools Referenced

- **Fréchet–Hoeffding bounds** for extremal joint CDFs.  
  *(BFbounds.md)*
- **Hardy–Littlewood–Pólya rearrangement inequality** for rank-based maximization/minimization.  
  *(Ordinal Variables and Correlation.md)*
- **Comonotone/countermonotone quantile coupling**.  
  *(Math Proof GPT – Whitt 1976 Paper.md)*
- **Capacity limits** in discrete joint tables explain why FH lower may be unattainable.

---

# 8. Master Takeaways

1. **\(r_{\max}\) = FH upper (comonotonic), always attainable** in discrete ordinal settings.
2. **\(r_{\min}\)** may **not** equal FH lower in discrete settings; use best feasible opposite-order arrangement.
3. No universal bound on relative magnitudes: \(|r_{\min}|\) may be less than, equal to, or greater than \(r_{\max}\).
4. Palindromic symmetry of one marginal guarantees \(|r_{\min}| = r_{\max}|\).
5. Binary cases illustrate when perfect alignment is feasible and when skew prevents it.

---
