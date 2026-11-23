# Meeting Summary: 23 October 2026 (Scott & Cees)

About 110 minutes. It was our first working session since 7 July, and the main item was Scott's novelty memo of 15/20 October (`Notes/20261020_Memo_Scott_Cees_ElasticBounds_Novelty.md`). Transcripts: `Meetings/Voice 261023_133746_sd*.txt`.

---

## Summary

1. **We accept that most of the theory is already published.** Cees said the memo "reinforces a doubt that I expressed earlier": at a high level we don't have a large amount that is new. We both read the recent papers (Shojima 2024/25, Fowler 2026) as confirming our thinking rather than sinking the project. The literature is "not non-existent", but applied researchers don't notice it. We agreed the project is still worth doing, as long as it doesn't become "a 10-year project".

2. **New framing: mechanical versus substantive effects.** Instead of "the bounds of r", which are well described in the literature, the paper would ask how much of a difference in r between groups, or over time, comes from the marginals being different rather than the association being different. Scott put it as "are the betas changing, or are the Xs changing?" The mathematical results become smaller supporting contributions, along with comments on the recent literature. Precedents Cees named:
   - **Stouthardt** (Dutch political sociologist; *Acta Politica*, 1970s). He separated religion's effect on party choice from changes in religion's marginals, by asking what the result would be if the conditional relationship stayed fixed and only the margins changed.
   - ***Electoral Change*** (probably Franklin, Mackie & Valen 1992), which Cees contributed to. Social structure explains less and less of the vote, partly perhaps because the marginals of the social-structure variables changed.

3. **Cees prefers comparing groups within one survey to trend studies.** Trend studies need at least four roughly evenly spaced time points, plus a lot of finding and cleaning data. A comparison within one cross-section is simpler, and the corollary for trend studies can then be noted. Examples: high vs middle vs low political interest, sophistication or knowledge; men vs women; education. There is a necessary condition: the marginals must differ a lot, otherwise there is no mechanical effect to find. Candidate settings with large changes in the marginals:
   - **The Dutch Parliamentary Election Study 1986 panel around Chernobyl.** Cees was a principal investigator. The nuclear-energy items were repeated in wave 2, and about 80% of wave 1 was interviewed before Chernobyl. The 20% interviewed afterwards were reached before the parties had repositioned. Many publications use these data.
   - BES before and after the Brexit referendum.
   - Eurobarometer items on Russia and the EU, before and after the 2022 invasion of Ukraine (if such items exist).
   - The over-time data in Sullivan et al. (Vietnam war).

   Scott's quick ANES check of Baldassarri & Gelman found r_max roughly flat across years. We agreed their data change too slowly to make a good test bed.

4. **Option B: a paper for practitioners**, for a survey-methodology journal or *Political Behavior*. It would open with a few recent articles that use r as if there were no problem. It would then summarise what is "known but ignored", set out what survey researchers should and shouldn't do, and illustrate with one or two examples. As a remedy it would offer the split rescaling that Shojima also proposes (different for positive and negative r), plus evidence that rescaled correlation matrices remain usable. Fowler (2026) shows such venues exist. Cees summed up the two options as "(A) mechanical versus causal effects, for different times and different groups; (B) a paper for practitioners".

5. **Whether to combine the bounds into one measure (C1).** Cees re-ran his beta regressions separately for the upper and lower bound. The average marginal effects, and even which predictors are significant, differ between the two. So is it justified to combine them into C1? If we don't combine them, we must explain why the drivers differ. Scott linked this to skew and symmetry and suggested trying both. **Not resolved.**

6. **What happened to the 1978–79 AJPS debate.** Cees said it petered out, with big names on both sides. Everyone agreed that the ANES changes to question wording and format affected naive measures of constraint. Nobody could say how much real change remained after allowing for that.

7. **Invertibility and factor analysis.** Scott was disappointed that ignoring the problem, or rescaling, rarely produces an unusable correlation matrix, which leaves no fix to sell. He thought the SEM/factor-analysis angle "might be too much" for this paper. Cees argued for a narrower scope. Whether the corrected matrix stays invertible is "incredibly important": if it does, "there is no strong argument against" correcting, and that is a simple remedy we can offer. Scott: "here's the problem, here are the solutions", which gives practitioners reassurance.

8. **Points for footnotes.**
   - *Shojima's degenerate solutions:* he gets "degenerate" solutions when a category has zero frequency. Cees thinks the sorting approach avoids this.
   - *Shojima's simulations* use only binomial marginals.
   - *Shojima's concluding remark* that χ² and Cramér's V are also constrained by the marginals. This is interesting, but we won't pursue it now.
   - *Proportions vs counts and missing data* (Cees). Shojima works with proportions, which implicitly assumes that "don't know" responses would be distributed like the valid ones, whereas in our approach they simply drop out. Cees will check whether this matters.
   - *"Why not polychorics?"* (the referee question Scott expects). Polychorics are estimated rather than calculated, often give matrices that aren't positive definite, and their algorithms sometimes fail; cite van der Eijk & Rose (2015). Scott added that the usual fix, the nearest PSD matrix, is a hack.

---

## Homework

### Scott (SM)
1. **Write and send this summary and homework list.** Scott promised it at the end of the meeting.
2. **Literature on three topics:**
   - mechanical vs causal/substantive effects, starting from Stouthardt once Cees finds it;
   - decomposition methods. Scott recalled an older decomposition literature, probably the Kitagawa / Oaxaca–Blinder kind;
   - the concept and measurement of constraint, including the dissertation by Stefano Camettari's wife (ask Cees for her name and the reference).
3. **Re-read Shojima** on proportions and relative marginals, and check Cees's missing-data point.
4. **Re-read Fowler (2026).** It is more relevant than Scott first thought, and shows there is a venue for this kind of work.
5. **Look for candidate data sets for option A**, i.e. groups or periods with large differences in the marginals. Start with BES before and after the referendum.
6. **Not this project:** send Cees the 1½-page outline of the ESRC grant idea (with the Bayesian epidemiologist), for his comments.

### Cees (CvdE)
1. **Send his own summary of the meeting from his notes**, starting Monday 26 October. He also said he would "try to make these two things"; the transcript is unclear, but it most likely means sketching options A and B.
2. **Find the Stouthardt article** in his *Acta Politica* run, and send it even if it's in Dutch ("tables are tables").
3. **Check the missing-data point** about proportions vs counts in Shojima.
4. **Suggest cross-section comparisons and data**, such as levels of sophistication or interest, men vs women, the DPES 1986 Chernobyl panel (with fieldwork dates), and Eurobarometer. He said he had "a few things that I might think of".
5. **Send the name and reference of the dissertation on constraint** (by Stefano Camettari's wife).

### Joint / open
- **Choose between option A** (mechanical vs substantive; Political Analysis/PSRM if it works) **and option B** (practitioner paper; a survey-methodology journal or *Political Behavior*), or combine them. This needs to be settled before any more writing.
- **The role of the beta-regression section and of C1**, i.e. one combined measure vs the two bounds separately (item 5 of the summary).
- **Next meeting:** not scheduled.

### Earlier homework: status after this meeting
| Item | Owner | Status |
|---|---|---|
| Two illustrative examples (0.6-vs-0.5 rank reversal; factor analysis based on Risky Business) | SM | Not done. The factor-analysis example may be outside the narrower scope |
| Narrow the Hypothesis Testing section; fix the two errors in the Asymmetry section; cite Shojima prominently | SM | Not done. Still needed under either option |
| MC1 vs MC1A sampling comparison | SM | On hold; probably not needed under either option |
| ANES trend pilot (`R/08`, `R/09`) | SM | Done on branch `anes-trend-decomposition`, not merged. The moral-issue result wasn't discussed (see below) |
| Pilot simulation of spurious factors (`R/10`) | SM | Done on that branch; not discussed |
| Introduction and literature review; precedent for asymmetric stretching; footnote on the ordinal scope (Moss & Grønneberg) | CvdE | Not delivered. The introduction now depends on the A/B choice |
| Corrections to the September draft (BES 2024 → 2019; C1 naming; footnote 3) | CvdE | Not discussed |

---

## Notes added after the meeting (Claude, from the repo and transcripts)

- **We already have an option-A result that wasn't raised.** In the Kozlowski & Murphy re-analysis (`R/09`), moral-issue alignment looks flat or falling because r_max fell. Opinion on gay rights became lopsided (e.g. gay adoption 27% → 76% in favour), while alignment *within* the bounds held steady. That is a mechanical-versus-substantive result, with exactly the large change in the marginals that Cees says is necessary. It is over time, which Cees is wary of, but the data are already clean.
- **The bounds for upper and lower are *expected* to have different drivers.**
  - r_max = 1 − W₂²(Z_X, Z_Y)/2 depends on how different the two standardized shapes are.
  - r_min = −1 + W₂²(Z_X, −Z_Y)/2 depends on how different X is from Y *flipped end to end*.
  - So skew in opposite directions lowers r_max, and skew in the same direction pulls r_min towards zero. This answers Cees's question in summary item 5: modelling the two bounds separately follows from the theory.
- **The missing-data point is probably not an issue for the bounds.** They depend only on the proportions among the cases used for r, whether you write them as counts or proportions. What matters is that the marginals come from the same cases as r. In our BES file they are pairwise-complete, so they do. Worth a quick check with Cees.
- **State the invertibility evidence precisely.**
  - Raw matrices never failed.
  - Split rescaling (Shojima's ρ*, our `split_warrens`) failed PSD in about 10% of the adversarial synthetic sets and 0% of the BES subsets.
  - Also, PSD is not the same as invertible.
  - `R/06` currently can't show how close a matrix is to singular, because its minimum-eigenvalue tracker starts at 0.
- **Comparisons within one cross-section need respondent-level data.** The BES pairs file holds only summaries for each pair. Splitting by political interest, for example, needs the BES 2019 face-to-face microdata, which include interest and attention items (v1, v61).
