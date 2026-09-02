#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_summary_tables.py
=========================
Recompute all summary tables from the per-question evaluation outputs.

Inputs : one CSV per model/dataset named  <Model>_<target|nontarget>.csv
         (as produced by run_dify_evaluation.py)
Outputs: Table1_Performance_Comparison.csv   accuracy + McNemar significance,
                                         LLM-judge flag AND objective choice-match
         Table_S_Species_Stratified.csv   per-species accuracy with Wilson 95% CI
         Table_S_Domain_SampleSizes.csv   questions per research domain

Usage:
    python compute_summary_tables.py --inputs model_outputs/ --outdir summary_tables/
"""
import argparse
import csv
import glob
import math
import os
import sys

MODELS_ORDER = ['DeepSeek-V3.2', 'Grok-4.1-Fast', 'Doubao-1.8',
                'MiMo-V2-Flash', 'Qwen3-235B', 'Gemini-3-Flash-High']
DATASETS = ['target', 'nontarget']


def load(p):
    for enc in ('utf-8-sig', 'gbk', 'latin1'):
        try:
            with open(p, encoding=enc, newline='') as f:
                return list(csv.DictReader(f))
        except UnicodeDecodeError:
            continue
    raise IOError('cannot decode ' + p)


def nb(v):
    s = str(v).strip().lower()
    if s in ('true', '1', '正确'):
        return 'True'
    if s in ('false', '0'):
        return 'False'
    return ''


def obj(r, side):
    ch = str(r.get(side + '_choice') or '').strip()
    ca = str(r.get('correct_answer') or '').strip()
    return ch in ('0', '1', '2') and ch == ca


def mcnemar_p(b, c):
    """Chi-square (df=1) p-value with continuity correction."""
    if b + c == 0:
        return 1.0
    chi2 = (abs(b - c) - 1) ** 2 / (b + c)
    return 0.5 * math.erfc(math.sqrt(chi2 / 2))


def sig(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    center = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, center - half) * 100, min(1.0, center + half) * 100)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', required=True, help='folder with <Model>_<dataset>.csv files')
    ap.add_argument('--outdir', required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    DATA = {}
    for p in glob.glob(os.path.join(args.inputs, '*.csv')):
        name = os.path.splitext(os.path.basename(p))[0]
        if name.endswith('_target'):
            m, ds = name[:-len('_target')], 'target'
        elif name.endswith('_nontarget'):
            m, ds = name[:-len('_nontarget')], 'nontarget'
        else:
            continue
        rows = load(p)
        for r in rows:
            r['with_mcp_correct'] = nb(r.get('with_mcp_correct'))
            r['no_mcp_correct'] = nb(r.get('no_mcp_correct'))
        DATA.setdefault(m, {})[ds] = rows
        print('loaded %s (%d rows)' % (os.path.basename(p), len(rows)))

    models = [m for m in MODELS_ORDER if m in DATA] + \
             [m for m in DATA if m not in MODELS_ORDER]

    # ---------- Table 1 ----------
    t1 = [['Model', 'Dataset', 'N', 'Failed', 'Empty_Agent_Answer',
           'With_Tool_flag_%', 'Without_Tool_flag_%', 'Delta_flag_%',
           'With_raw_flag', 'Without_raw_flag', 'Significance_flag',
           'With_Tool_obj_%', 'Without_Tool_obj_%', 'Delta_obj_%',
           'With_raw_obj', 'Without_raw_obj', 'Significance_obj']]
    for ds in DATASETS:
        for m in models:
            rows = DATA[m].get(ds, [])
            if not rows:
                continue
            n = len(rows)
            failed = sum(1 for r in rows if str(r['status']).strip().lower() != 'success')
            empty = sum(1 for r in rows if not str(r.get('with_mcp_answer') or '').strip())
            wf = sum(1 for r in rows if r['with_mcp_correct'] == 'True')
            lf = sum(1 for r in rows if r['no_mcp_correct'] == 'True')
            b = sum(1 for r in rows if r['no_mcp_correct'] == 'True' and r['with_mcp_correct'] != 'True')
            c = sum(1 for r in rows if r['with_mcp_correct'] == 'True' and r['no_mcp_correct'] != 'True')
            wo = sum(1 for r in rows if obj(r, 'with_mcp'))
            lo = sum(1 for r in rows if obj(r, 'no_mcp'))
            bo = sum(1 for r in rows if obj(r, 'no_mcp') and not obj(r, 'with_mcp'))
            co = sum(1 for r in rows if obj(r, 'with_mcp') and not obj(r, 'no_mcp'))
            t1.append([m, ds, n, failed, empty,
                       '%.2f' % (100 * wf / n), '%.2f' % (100 * lf / n), '%+.2f' % (100 * (wf - lf) / n),
                       '%d/%d' % (wf, n), '%d/%d' % (lf, n), sig(mcnemar_p(b, c)),
                       '%.2f' % (100 * wo / n), '%.2f' % (100 * lo / n), '%+.2f' % (100 * (wo - lo) / n),
                       '%d/%d' % (wo, n), '%d/%d' % (lo, n), sig(mcnemar_p(bo, co))])
    p1 = os.path.join(args.outdir, 'Table1_Performance_Comparison.csv')
    with open(p1, 'w', encoding='utf-8-sig', newline='') as f:
        csv.writer(f, lineterminator='\n').writerows(t1)
    print('written', p1)

    # ---------- Species stratified (Wilson CI, long format) ----------
    out = [['Dataset', 'Species', 'Model', 'n', 'Agent_correct', 'Agent_acc_%', 'Agent_CI95',
            'LLM_correct', 'LLM_acc_%', 'LLM_CI95']]
    for ds in DATASETS:
        ref = DATA[models[0]].get(ds, [])
        groups = {}
        for r in ref:
            sp = str(r.get('species') or '').strip() or '(unlabeled)'
            groups.setdefault(sp, []).append(str(r['index']).strip())
        for sp, idxs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            n = len(idxs)
            idxset = set(idxs)
            for m in models:
                sub = [r for r in DATA[m].get(ds, []) if str(r['index']).strip() in idxset]
                ac = sum(1 for r in sub if r['with_mcp_correct'] == 'True')
                lc = sum(1 for r in sub if r['no_mcp_correct'] == 'True')
                alo, ahi = wilson(ac, n)
                llo, lhi = wilson(lc, n)
                out.append([ds, sp, m, n, ac, '%.2f' % (100 * ac / n), '%.1f-%.1f' % (alo, ahi),
                            lc, '%.2f' % (100 * lc / n), '%.1f-%.1f' % (llo, lhi)])
    p2 = os.path.join(args.outdir, 'Table_S_Species_Stratified.csv')
    with open(p2, 'w', encoding='utf-8-sig', newline='') as f:
        csv.writer(f, lineterminator='\n').writerows(out)
    print('written', p2)

    # ---------- Domain sample sizes ----------
    out2 = [['Dataset', 'Research_Area', 'n']]
    for ds in DATASETS:
        cnt = {}
        for r in DATA[models[0]].get(ds, []):
            a = str(r.get('area') or '').strip()
            cnt[a] = cnt.get(a, 0) + 1
        for a, n in sorted(cnt.items(), key=lambda kv: -kv[1]):
            out2.append([ds, a, n])
    p3 = os.path.join(args.outdir, 'Table_S_Domain_SampleSizes.csv')
    with open(p3, 'w', encoding='utf-8-sig', newline='') as f:
        csv.writer(f, lineterminator='\n').writerows(out2)
    print('written', p3)


if __name__ == '__main__':
    main()