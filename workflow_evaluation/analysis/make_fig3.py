# -*- coding: utf-8 -*-
import csv, os, sys, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Liberation Sans', 'Helvetica']
plt.rcParams['svg.fonttype'] = 'path'

FINAL = r'C:\Users\a1835\Desktop\mcp\评估\汇总结果\final_data'
OUT = r'C:\Users\a1835\Desktop\mcp\又又又重投\NP-R1\Figure 3. Research area-specific performance of PGCopilot across six large language models.svg'
MODELS = ['DeepSeek-V3.2', 'Grok-4.1-Fast', 'Doubao-1.8', 'MiMo-V2-Flash', 'Qwen3-235B', 'Gemini-3-Flash-High']
MLAB = ['DeepSeek-V3.2', 'Grok-4.1-Fast', 'Doubao-1.8', 'MiMo-V2-Flash', 'Qwen3-235B', 'Gemini-3-Flash-High']

DISPLAY = {
 'CELL BIOLOGY AND CELL SIGNALING': 'Cell Biology',
 'ENVIRONMENT - ABIOTIC STRESS': 'Abiotic Stress',
 'ENVIRONMENT - BIOTIC STRESS': 'Biotic Stress',
 'ENVIRONMENT - LIGHT AND TEMPERATURE': 'Light & Temperature',
 'ENVIRONMENT - NUTRIENTS': 'Nutrients',
 'ENVIRONMENT - PLANT-SYMBIONTS': 'Plant-Symbionts',
 'EVOLUTION': 'Evolution',
 'GENE REGULATION - ALTERNATIVE SPLICING': 'Alternative Splicing',
 'GENE REGULATION - EPIGENETICS AND TGS': 'Epigenetics & TGS',
 'GENE REGULATION - EPITRANSCRIPTOMICS AND RNA STRUCTURE': 'Epitranscriptomics',
 'GENE REGULATION - POST-TRANSLATIONAL MODIFICATIONS': 'Post-translational Mod.',
 'GENE REGULATION - PTGS': 'PTGS',
 'GENE REGULATION - TRANSCRIPTION': 'Transcription',
 'GENE REGULATION - TRANSLATION': 'Translation',
 'GENOME AND GENOMICS': 'Genome & Genomics',
 'GROWTH AND DEVELOPMENT': 'Growth & Development',
 'HORMONES': 'Hormones',
 'PHYSIOLOGY AND METABOLISM': 'Physiology & Metabolism',
 'PLANT BIOTECHNOLOGY': 'Plant Biotechnology',
}
ORDER = ['Cell Biology', 'Abiotic Stress', 'Biotic Stress', 'Light & Temperature', 'Nutrients',
         'Plant-Symbionts', 'Evolution', 'Alternative Splicing', 'Epigenetics & TGS',
         'Post-translational Mod.', 'PTGS', 'Transcription', 'Translation', 'Genome & Genomics',
         'Growth & Development', 'Hormones', 'Physiology & Metabolism', 'Plant Biotechnology',
         'Epitranscriptomics']

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

def matrix(ds):
    rows = load('%s\\%s_%s.csv' % (FINAL, MODELS[0], ds))
    present = {}
    for r in rows:
        a = (r.get('area') or '').strip()
        present[DISPLAY.get(a, a)] = True
    order = [o for o in ORDER if o in present]
    avg = np.zeros((len(order), len(MODELS)))
    for j, m in enumerate(MODELS):
        for r in load('%s\\%s_%s.csv' % (FINAL, m, ds)):
            label = DISPLAY.get((r.get('area') or '').strip(), (r.get('area') or '').strip())
            pass
    for j, m in enumerate(MODELS):
        for k, label in enumerate(order):
            sub = [r for r in load('%s\\%s_%s.csv' % (FINAL, m, ds))
                   if DISPLAY.get((r.get('area') or '').strip(), (r.get('area') or '').strip()) == label]
            avg[k, j] = 100.0 * sum(1 for r in sub if nb(r.get('with_mcp_correct'))) / len(sub)
    return order, avg

torder, T = matrix('target')
norder, N = matrix('nontarget')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(130.5/25.4, 97.9/25.4), dpi=600,
                               gridspec_kw={'height_ratios': [1.05, 1.0], 'hspace': 0.05})
cmap = 'RdYlBu_r'
im1 = ax1.imshow(T, cmap=cmap, aspect='auto', vmin=0, vmax=100)
im2 = ax2.imshow(N, cmap=cmap, aspect='auto', vmin=0, vmax=100)
ax1.set_xticks(range(6)); ax1.set_xticklabels([])
ax1.set_yticks(range(len(torder))); ax1.set_yticklabels(torder, fontsize=4.5)
ax2.set_xticks(range(6)); ax2.set_xticklabels(MLAB, fontsize=4.5, rotation=0)
ax2.set_yticks(range(len(norder))); ax2.set_yticklabels(norder, fontsize=4.5)
ax1.set_title('A  Target Plant', fontsize=7, loc='left')
ax2.set_title('B  Non-Target Plant', fontsize=7, loc='left')
for ax in (ax1, ax2):
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

cbar = fig.colorbar(im1, ax=[ax1, ax2], fraction=0.02, pad=0.02)
cbar.set_label('Score', fontsize=6)
cbar.ax.tick_params(labelsize=5)

plt.savefig(OUT, format='svg', bbox_inches=None, pad_inches=0)
print('已生成 Fig3 SVG:', OUT)
print('target 领域数:', len(torder), ' non-target 领域数:', len(norder))
print('target 顺序:', torder)
print('non-target 顺序:', norder)