#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mixed_effects_bootstrap.py
==========================
Supplementary statistics for the double-blind expert evaluation
(manuscript NPH-MET-2026-58029, response to Referee 1, Main point 4).

Complements the ICC analysis with:

1. A mixed-effects linear model
       score ~ C(model)                      (model = fixed effect)
       + random intercept for Evaluator      (crossed random effect)
       + random intercept for Question       (crossed random effect)
   fitted as variance components (statsmodels MixedLM, groups = observation).

2. A clustered bootstrap over questions (1,000 resamples, percentile CIs)
   for the paired mean differences against the reference system.

Outputs
-------
  Table S11 Mixed-effects estimates.xlsx  (Estimates / Coefficients / Notes)
  Table S11 forest plot.png               (forest plot, 300 dpi)

Input: anonymized expert-evaluation workbook with a sheet "Reviewed Scores"
containing 1,050 ratings (15 evaluators x 10 questions x 7 metrics x 4 systems).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

XLSX = r'C:\Users\a1835\Desktop\mcp\又又又重投\NP-R1\data\Supplementary Data1_anonymized.xlsx'
OUT_XLSX = r'C:\Users\a1835\Desktop\mcp\又又又重投\NP-R1\data\Table S11 Mixed-effects estimates.xlsx'
OUT_PNG = r'C:\Users\a1835\Desktop\mcp\又又又重投\NP-R1\data\Table S11 forest plot.png'

REF = 'Full PGC Agent'
OTHERS = ['PlantGPT', 'Ablated General LLM', 'PGC Agent without System Prompt']
MODELS = [REF] + OTHERS
N_BOOT = 1000
SEED = 20260831

# ---------- load ----------
df = pd.read_excel(XLSX, sheet_name='Reviewed Scores')
long = df.melt(id_vars=['Evaluator', 'Question', 'Metric'],
               value_vars=MODELS, var_name='model', value_name='score')
long['score'] = pd.to_numeric(long['score'])
long = long.dropna(subset=['score'])
print('ratings melted:', len(long), '(expect 4200)')

# ---------- 1) mixed-effects (crossed random intercepts) ----------
import statsmodels.formula.api as smf
ld = long.copy()
ld['row'] = range(len(ld))
FORMULA = 'score ~ C(model, Treatment(reference="%s"))' % REF
mixed_df = None
mixed_note = ''
attempts = [
    ('crossed variance components (lbfgs)',
     lambda: smf.mixedlm(FORMULA, data=ld, groups=ld['row'],
                         vc_formula={'evaluator': '0 + C(Evaluator)',
                                     'question': '0 + C(Question)'}).fit(method='lbfgs', disp=False)),
    ('crossed variance components (powell)',
     lambda: smf.mixedlm(FORMULA, data=ld, groups=ld['row'],
                         vc_formula={'evaluator': '0 + C(Evaluator)',
                                     'question': '0 + C(Question)'}).fit(method='powell', disp=False)),
    ('random intercept for evaluator, question as fixed strata',
     lambda: smf.mixedlm(FORMULA + ' + C(Question)', data=ld, groups=ld['Evaluator']).fit(method='lbfgs', disp=False)),
]
res = None
for note, fn in attempts:
    try:
        res = fn()
        mixed_note = note
        print('mixed-effects fitted via:', note)
        break
    except Exception as ex:
        print('fit failed (%s): %s' % (note, str(ex)[:80]))

if res is not None:
    coef = res.params
    ci = res.conf_int()
    names = {'C(model, Treatment(reference="%s"))[T.%s]' % (REF, m): m for m in OTHERS}
    mixed_rows = []
    for term, m in names.items():
        if term in coef.index:
            # statsmodels dummy = Other - Reference; flip sign so that
            # positive = Full PGC Agent higher (consistent with bootstrap)
            mixed_rows.append((m, -coef[term], -ci.loc[term, 1], -ci.loc[term, 0]))
    mixed_df = pd.DataFrame(mixed_rows, columns=['Contrast (vs %s)' % REF, 'Estimate', 'CI_low', 'CI_high'])
else:
    print('!! all mixed-model specifications failed -> bootstrap-only output')

# ---------- 2) clustered bootstrap over questions ----------
rng = np.random.default_rng(SEED)
questions = long['Question'].unique()
qidx = {q: i for i, q in enumerate(questions)}
arr = np.full((len(long), len(MODELS)), np.nan)
qi = np.zeros(len(long), dtype=int)
mpos = {m: i for i, m in enumerate(MODELS)}
for i, (_, r) in enumerate(long.iterrows()):
    arr[i, mpos[r['model']]] = r['score']
    qi[i] = qidx[r['Question']]
qmask = {q: (qi == i) for q, i in qidx.items()}

draws = {m: np.zeros(N_BOOT) for m in OTHERS}
for b in range(N_BOOT):
    sel = rng.choice(len(questions), size=len(questions), replace=True)
    mask = np.zeros(len(long), dtype=bool)
    for s in sel:
        mask |= qmask[questions[s]]
    means = np.nanmean(arr[mask], axis=0)
    for m in OTHERS:
        draws[m][b] = means[mpos[REF]] - means[mpos[m]]

boot_rows = []
for m in OTHERS:
    d = draws[m]
    boot_rows.append((m, d.mean(), np.percentile(d, 2.5), np.percentile(d, 97.5),
                      float(2 * min((d <= 0).mean(), (d >= 0).mean()))))
boot_df = pd.DataFrame(boot_rows, columns=['Contrast (vs %s)' % REF, 'Boot_mean',
                                           'Boot_CI_low', 'Boot_CI_high', 'Boot_p_two_sided'])
print(boot_df)

# ---------- 3) outputs ----------
if mixed_df is not None:
    est = pd.merge(mixed_df, boot_df, on='Contrast (vs %s)' % REF)
    est.columns = ['Contrast', 'Mixed_estimate', 'Mixed_CI_low', 'Mixed_CI_high',
                   'Bootstrap_mean', 'Bootstrap_CI_low', 'Bootstrap_CI_high', 'Bootstrap_p_two_sided']
else:
    est = boot_df.copy()
    est.columns = ['Contrast', 'Bootstrap_mean', 'Bootstrap_CI_low', 'Bootstrap_CI_high', 'Bootstrap_p_two_sided']

with pd.ExcelWriter(OUT_XLSX, engine='openpyxl') as xw:
    est.to_excel(xw, sheet_name='Estimates', index=False)
    if res is not None:
        coef_df = pd.DataFrame({'term': coef.index, 'estimate': coef.values,
                                'CI_low': ci[0].values, 'CI_high': ci[1].values})
        coef_df.to_excel(xw, sheet_name='Mixed_model_coefficients', index=False)
    notes = pd.DataFrame({'Item': [
        'Analysis', 'Model fixed effect', 'Random effects', 'Specification used',
        'Bootstrap', 'N ratings', 'Design', 'Significance', 'Seed'],
        'Description': [
        'Supplementary analysis of expert evaluation ratings (complements ICC, Table S6)',
        'score ~ C(model); reference = %s' % REF,
        'crossed random intercepts for Evaluator (15) and Question (10)',
        mixed_note or 'mixed model failed; clustered bootstrap only',
        'clustered bootstrap, %d resamples over questions, percentile 95%% CI' % N_BOOT,
        '%d (15 evaluators x 10 questions x 7 metrics x 4 systems)' % len(long),
        'all ratings of a resampled question enter together, preserving within-question dependence',
        'two-sided bootstrap p = 2 x min(P(diff<=0), P(diff>=0))',
        str(SEED)]})
    notes.to_excel(xw, sheet_name='Notes', index=False)
print('written', OUT_XLSX)

# ---------- 4) forest plot ----------
fig, ax = plt.subplots(figsize=(7.2, 3.4), dpi=300)
ypos = np.arange(len(est))[::-1]
for y, (_, r) in zip(ypos, est.iterrows()):
    ax.plot([r['Mixed_CI_low'], r['Mixed_CI_high']], [y + 0.13, y + 0.13],
            color='#1f4e79', lw=2, solid_capstyle='round')
    ax.plot(r['Mixed_estimate'], y + 0.13, 's', color='#1f4e79', ms=6, label='Mixed model (95% CI)' if y == ypos[0] else None)
    ax.plot([r['Bootstrap_CI_low'], r['Bootstrap_CI_high']], [y - 0.13, y - 0.13],
            color='#c0504d', lw=2, solid_capstyle='round')
    ax.plot(r['Bootstrap_mean'], y - 0.13, 'o', color='#c0504d', ms=6, label='Clustered bootstrap (95% CI)' if y == ypos[0] else None)
ax.axvline(0, color='grey', ls='--', lw=1)
ax.set_yticks(ypos + 0.0)
ax.set_yticklabels([c.replace(' (vs %s)' % REF, '') for c in est['Contrast']], fontsize=9)
ax.set_xlabel('Mean score difference vs %s (points, 0-10 scale)' % REF, fontsize=9)
ax.tick_params(axis='x', labelsize=8)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(fontsize=8, loc='lower right', frameon=False)
ax.set_title('Expert evaluation: adjusted model contrasts', fontsize=10)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
print('written', OUT_PNG)