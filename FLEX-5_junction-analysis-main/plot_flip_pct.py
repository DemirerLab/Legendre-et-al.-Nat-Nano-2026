import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE    = os.path.dirname(os.path.abspath(__file__))
SVG_DIR = os.path.join(BASE, 'svg')
os.makedirs(SVG_DIR, exist_ok=True)

# lox2272-flip as % of all classified flip reads (from Table 1)
DATA = {
    'Neg':            [100*1566/1566,   100*92/92],
    'Agro\n(+ctrl)':  [100*4573/4574,   100*4352/4352],
    'PVC\n(normal)':  [100*3441/3441,   100*2601/2604],
    'PVC\npFLS2':     [100*4291/4291,   100*4364/4364],
    'PVC\npUBQ10':    [100*4143/4143,   100*4347/4347],
}

conditions = list(DATA.keys())
means = [np.mean(v) for v in DATA.values()]
x = np.arange(len(conditions))

fig, ax = plt.subplots(figsize=(5, 3.8))

ax.bar(x, means, color='#bcbddc', width=0.55, edgecolor='white', linewidth=0)
for i, vals in enumerate(DATA.values()):
    ax.scatter([i]*len(vals), vals, color='#6a51a3', s=22, zorder=3, linewidths=0)

ax.set_xticks(x)
ax.set_xticklabels(conditions, fontsize=9)
ax.set_ylabel('lox2272-flip (%)', fontsize=10)
ax.set_ylim(0, 100)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.yaxis.set_tick_params(labelsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
out = os.path.join(BASE, 'Figure 4 — lox2272-flip percentage per condition.png')
plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(SVG_DIR, os.path.basename(out).replace('.png', '.svg')), bbox_inches='tight', facecolor='white')
print('Saved:', out)
