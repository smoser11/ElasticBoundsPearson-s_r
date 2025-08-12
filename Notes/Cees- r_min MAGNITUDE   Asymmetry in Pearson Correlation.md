Thank you, but I think I found a counter example.  I am using 'Lane A' for both R-min and r_man.  Here are the values of countsX: 641, 1648,  806,  306 (so, four categorys, with the sample size of 3401.  Using the co-monotonic ordering (using the values 0,1,2,3) I obtain r_max = 0.8502718.  Next I will show you a countsY that results in an r_min of -0.8991577.

* * *

**You:**

Continuing this example, Here are the values of countsY:  66,  820,  750, 1223,  542.  Using the co-monotonic ordering (using the values 0,1,2,3) I obtain a r_min of -0.8991577.  How can this be? Please note this comes from actual data (the BES), not simulated data.  

* * *

**ChatGPT:**

my best guess is the numbers you posted are **correct**, and they actually reveal the subtle point: even when you use “Lane A” (pure permutations of the same atoms), it is **possible** to have

$$|r_{\min}| \;>\; r_{\max}.$$

Let me show why, step-by-step, using your exact counts and then the general intuition.

* * *

1) Verify with your counts (what you’re doing is consistent)
------------------------------------------------------------

* $X$ has 4 categories scored $0,1,2,3$ with counts  
    $[641,\,1648,\,806,\,306]$ (total $N=3401$).
    
* $Y$ has 5 categories scored $0,1,2,3,4$ with counts  
    $[66,\,820,\,750,\,1223,\,542]$ (also $N=3401$).
    

Construct the two length-$N$ vectors by replication, then compute permutation extremes:

* $r_{\max}^{\text{perm}}=\mathrm{cor}(\text{sort}(X,\downarrow),\,\text{sort}(Y,\downarrow))=0.8502718$.
    
* $r_{\min}^{\text{perm}}=\mathrm{cor}(\text{sort}(X,\downarrow),\,\text{sort}(Y,\uparrow))=-0.8991577$.
    

So indeed $|r_{\min}^{\text{perm}}| > r_{\max}^{\text{perm}}$ for this **real** BES pair.

> Note: if you build the **Fréchet–Hoeffding** (FH) comonotone/antitone **joint tables** from these marginals and compute $r$, you get the **same** two values here. So this is not an artifact of the permutation vs FH choice.

* * *

2) Why the usual rearrangement logic doesn’t give $|r_{\min}| \le r_{\max}$
---------------------------------------------------------------------------

* The **rearrangement inequality** guarantees only that
    
    $$\sum x_i y_i^{\downarrow} \;\ge\; \sum x_i y_i^{\uparrow}$$
    
    (same-order inner product is **larger** than opposite-order).
    
* Pearson correlation, however, uses **centered** and **scaled** variables:
    
    $$r \;=\; \frac{\sum (x_i-\bar x)\,(y_i-\bar y)}{\sqrt{\sum (x_i-\bar x)^2}\,\sqrt{\sum (y_i-\bar y)^2}}.$$
* Subtracting means changes the contributions of each matched pair. Rearrangement still says
    
    $$S_{\max} \equiv \sum (x_i-\bar x)\,y_i^{\downarrow} \;\ge\; S_{\min} \equiv \sum (x_i-\bar x)\,y_i^{\uparrow},$$
    
    but **this does not imply** $S_{\max} \ge |S_{\min}|$ (equivalently $S_{\max}+S_{\min}\ge 0$). Your example is precisely a case where $S_{\max} < |S_{\min}|$.
    

* * *

3) Intuition from your marginals (why the negative side is “stronger”)
----------------------------------------------------------------------

Compute the marginal means (with the standard 0,1,2,3 and 0–4 scores):

* $$\bar X \approx 1.228$$
    
* $$\bar Y \approx 2.398$$
    

Two key facts:

1. Because $\bar Y$ is **well above** the mid-scale (2), the **low** Y categories (0 and 1) lie **far below** the mean:
    
    $$0 - \bar Y \approx -2.398,\quad 1 - \bar Y \approx -1.398,$$
    
    whereas the high categories are much closer:
    
    $$3 - \bar Y \approx +0.602,\quad 4 - \bar Y \approx +1.602.$$
2. In the **antitone** pairing, **large** $X-\bar X>0$ values get matched with those **very negative** $Y-\bar Y$ deviations (0s, then 1s). These products are **large in magnitude and negative**.  
    In the **comonotone** pairing, large $X-\bar X>0$ gets matched mostly with $Y-\bar Y$ that is **only moderately positive** (3s and 4s aren’t far above $\bar Y$), so the positive products are **smaller**.
    

Put differently: with a **skewed** $Y$ that has $\bar Y$ high, the **opposite-order** match creates **bigger deviations** (in absolute value) than the same-order match, making $|\mathrm{Cov}_{\min}|>\mathrm{Cov}_{\max}$ and hence $|r_{\min}|>r_{\max}$.

This effect is especially common when:

* the two variables have **different numbers of categories** (here 4 vs 5, so the score ranges and mean locations differ), and/or
    
* one marginal is **asymmetric**, pushing its mean away from the mid-scale.
    

* * *

4) When do we get symmetry $\;r_{\min} = -r_{\max}\;$?
------------------------------------------------------

A clean **sufficient** condition (no need for both sides):

* If **one** marginal (say $Y$) is **palindromically symmetric** about its center rank, i.e. $p_Y(k)=p_Y(K_Y-1-k)$, then defining $Y^\* = (K_Y-1) - Y$ leaves the marginal unchanged; comonotone coupling with $Y^\*$ is antitone with $Y$; and
    
    $$\mathrm{corr}(X, Y^\*) = -\,\mathrm{corr}(X, Y).$$
    
    Because $Y^\*$ has the **same marginal** as $Y$, the maximal correlation with $Y^\*$ equals $r_{\max}$ for $Y$, and therefore $r_{\min} = -r_{\max}$.
    

In your example $Y$ is **not** symmetric (66 ≠ 542, 820 ≠ 1223), so no magnitude symmetry is expected, and indeed we see $|r_{\min}|>r_{\max}$.

* * *

5) Practical takeaways
----------------------

* There is **no general law** forcing $|r_{\min}| \le r_{\max}$ for Pearson (or Spearman) extremes.  
    What always holds is the **ordering** $r_{\max} \ge r \ge r_{\min}$ for any feasible coupling—**not** symmetry of magnitudes.
    
* Your BES counterexample is **valid**. It arises because $Y$ is skewed (mean well above mid-scale) and has a different score range than $X$, so the antitone pairing creates larger negative leverage than the comonotone pairing creates positive leverage.
    