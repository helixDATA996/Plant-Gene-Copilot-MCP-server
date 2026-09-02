# PGCopilot: human expert evaluation

This folder contains everything for the **double-blind human expert panel**
(manuscript NPH-MET-2026-58029, Referee 1 Main point 4). It complements
`../workflow_evaluation/` (the automated 308/257-question multiple-choice
benchmark).

## Contents

| Path | Purpose |
|---|---|
| `system/` | The web-based scoring system used to collect ratings from the 15 expert assessors (Flask + MySQL). See its own `README.md`. |
| `data/expert_scores_anonymized.xlsx` | All **1,050 anonymized ratings** (15 raters × 10 questions × 7 metrics × 4 systems), plus ICC, significance-test, rubric, and benchmark-question sheets. Evaluator identities are replaced by `Rater01`–`Rater15`; roles are preserved. |
| `analysis/mixed_effects_bootstrap.py` | Mixed-effects linear model + clustered bootstrap (1,000 resamples) of the model contrasts. Re-runnable to reproduce Table S11. |
| `results/` | Script outputs written here at run time (Table S11 estimates + forest plot). |

## Running the analysis

```bash
# from analysis/
python mixed_effects_bootstrap.py
```

Input: `../data/expert_scores_anonymized.xlsx` (sheet "Reviewed Scores").
Output: `../results/Table S11 Mixed-effects estimates.xlsx` and
`../results/Table S11 forest plot.png`.

## Privacy note

Only the anonymized ratings are distributed. The raw workbook containing
real rater identities is **not** published, in line with the double-blind
design and the manuscript's ethics statement.
