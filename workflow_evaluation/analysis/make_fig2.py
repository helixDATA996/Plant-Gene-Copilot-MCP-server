# -*- coding: utf-8 -*-
import csv, os, sys, io, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'Helvetica']
plt.rcParams['svg.fonttype'] = 'path'

FINAL = r'C:\Users\a1835\Desktop\mcp\评估\汇总结果\final_data'
OUT = r'C:\Users\a1835\Desktop\mcp\又又又重投\NP-R1\Figure 2. Performance comparison of AI models with and without MCP-based tool augmentation.svg'
MODELS = ['DeepSeek-V3.2', 'Grok-4.1-Fast', 'Doubao-1.8', 'MiMo-V2-Flash', 'Qwen3-235B', 'Gemini-3-Flash-High']
SHORT = ['DeepSeek-V3.2', 'Grok-4.1-Fast', 'Doubao-1.8', 'MiMo-V2-Flash', 'Qwen3-235B', 'Gemini-3-Flash-High']

def load(p):
    for enc in ('utf-8-sig', 'gbk', 'latin1'):
        try:
            with open(p, encoding=enc, newline='') as f:
                return list(csv.DictReader(f))
        except UnicodeDecodeError:
            continue
    raise IOError(p)

def nb(v):
    return str(v).strip().lower() == 'true'

def mcnemar(rows):
    b = sum(1 for r in rows if nb(r.get('no_mcp_correct')) and not nb(r.get('with_mcp_correct')))
    c = sum(1 for r in rows if nb(r.get('with_mcp_correct')) and not nb(r.get('no_mcp_correct')))
    if b + c == 0:
        return 1.0
    chi2 = (abs(b - c) - 1) ** 2 / (b + c)
    return 0.5 * math.erfc(math.sqrt(chi2 / 2))

def star(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else None

def acc(rows, side):
    return 100.0 * sum(1 for r in rows if nb(r.get(side + '_correct'))) / len(rows)

data = {}
for ds in ['target', 'nontarget']:
    for m in MODELS:
        rows = load('%s\\%s_%s.csv' % (FINAL, m, ds))
        data[(m, ds)] = (acc(rows, 'no_mcp'), acc(rows, 'with_mcp'), mcnemar(rows))

fig, axes = plt.subplots(1, 3, figsize=(174/25.4, 40/25.4), dpi=600,
                         gridspec_kw={'width_ratios': [1, 1, 0.95], 'wspace': 0.18})
c_with, c_without = '#dd8452', '#4c72b0'
x = np.arange(6)
w = 0.38

for ax, key, ylabel in [(axes[0], 'target', 'Target Plant Accuracy (%)'), (axes[1], 'nontarget', 'Non-Target Plant Accuracy (%)')]:
    no = [data[(m, key)][0] for m in MODELS]
    wi = [data[(m, key)][1] for m in MODELS]
    ax.bar(x - w/2, no, w, color=c_without, edgecolor='black', linewidth=0.5, label='Without MCP')
    ax.bar(x + w/2, wi, w, color=c_with, edgecolor='black', linewidth=0.5, label='With MCP')
    ax.set_xticks(x); ax.set_xticklabels(SHORT, fontsize=4.2)
    ax.set_ylim(0, 100); ax.set_ylabel(ylabel, fontsize=6)
    ax.tick_params(labelsize=5.5)
    for xi, m in zip(x, MODELS):
        s = star(data[(m, key)][2])
        if s:
            ax.text(xi, 5, s, ha='center', fontsize=6)

imp = [data[(m, 'target')][1] - data[(m, 'target')][0] for m in MODELS]
axes[2].bar(x, imp, 0.55, color='#dd8452', edgecolor='black', linewidth=0.5)
axes[2].axhline(0, color='black', linewidth=0.5)
axes[2].set_xticks(x); axes[2].set_xticklabels(SHORT, fontsize=4.2)
axes[2].set_ylim(-15, 25); axes[2].set_yticks([-15,-10,-5,0,5,10,15,20])
axes[2].set_ylabel('Accuracy Improvement (%)', fontsize=6)
axes[2].tick_params(labelsize=5.5)

for ax, lab in zip(axes, ['A', 'B', 'C']):
    ax.text(-0.07, 1.05, lab, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')
    ax.spines[['top', 'right']].set_visible(False)

h = [plt.Rectangle((0,0),1,1, color=c_without), plt.Rectangle((0,0),1,1, color=c_with)]
fig.legend(h, ['Without MCP', 'With MCP'], loc='upper center', ncol=2, fontsize=5, frameon=False,
           bbox_to_anchor=(0.5, 0.995), bbox_transform=fig.transFigure, borderpad=0.2)

plt.savefig(OUT, format='svg', bbox_inches=None, pad_inches=0)
print('已重新生成 Fig2 SVG:', OUT)
print('画布 174mm x 40mm（无 tight 裁剪）')