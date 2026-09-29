import sys, io, os, random
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FASTQ_DIR = os.environ.get('FLEX_FASTQ_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fastq'))
OUT_DIR   = os.path.dirname(os.path.abspath(__file__))
SVG_DIR   = os.path.join(OUT_DIR, 'svg')
os.makedirs(SVG_DIR, exist_ok=True)
def _svg(png_path): return os.path.join(SVG_DIR, os.path.basename(png_path).replace('.png', '.svg'))

# ── cowplot theme ─────────────────────────────────────────────────────────────
def apply_cowplot(font_size=14):
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': font_size,
        'axes.labelsize': font_size,
        'axes.titlesize': font_size,
        'xtick.labelsize': font_size - 2,
        'ytick.labelsize': font_size - 2,
        'axes.linewidth': 1.2,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'pdf.fonttype': 42,
        'svg.fonttype': 'none',
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
    })

def theme_half_open(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ── Pipeline helpers (duplicated from analyze_flex.py) ───────────────────────
FWD    = 'CCACTGACGTAAGGGATGACGCACAATC'
ANCHOR = 'CCACTATCCTTCGCAAGACCCTTCCTCTATATAAGGAAGTTCATTTCATTTGGAGAGG'
PROBE_LOXP    = 'ACAACGGCCGCGTATTTTTACAACAATTACCAACAACAACAAACAACAAACAACATTACA'
PROBE_LOX2272 = 'ACACGCTGACGGCCGCGTATTTTTACAACAATTACCAACAACAACAAACAACAAACAACA'
FWD_MM = 4; ANCHOR_MM = 8; K = 12
MIN_FLIP, MAX_FLIP = 480, 680
FWD_LEN = len(FWD); ANCHOR_LEN = len(ANCHOR)

def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

def find_fuzzy(seq, primer, max_mm):
    plen = len(primer); best_pos, best_mm = -1, max_mm + 1
    for i in range(len(seq) - plen + 1):
        mm = sum(a != b for a, b in zip(seq[i:i+plen], primer))
        if mm < best_mm:
            best_mm, best_pos = mm, i
            if mm == 0: break
    return (best_pos, best_mm) if best_mm <= max_mm else (-1, -1)

def kmer_score(window, probe, k=12):
    pk = set(probe[i:i+k] for i in range(len(probe)-k+1))
    return sum(1 for i in range(len(window)-k+1) if window[i:i+k] in pk)

def nw_align(seq, ref):
    m, n = len(seq), len(ref)
    GAP, MATCH, MISM = -2, 2, -1
    dp = [[0]*(n+1) for _ in range(m+1)]
    tb = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1): dp[i][0] = i*GAP; tb[i][0] = 1
    for j in range(n+1): dp[0][j] = j*GAP; tb[0][j] = 2
    for i in range(1,m+1):
        for j in range(1,n+1):
            sc = MATCH if seq[i-1]==ref[j-1] else MISM
            diag,up,left = dp[i-1][j-1]+sc, dp[i-1][j]+GAP, dp[i][j-1]+GAP
            best = max(diag,up,left); dp[i][j] = best
            tb[i][j] = 0 if best==diag else (1 if best==up else 2)
    i,j = m,n; aseq,aref = [],[]
    while i>0 or j>0:
        if i>0 and j>0 and tb[i][j]==0: aseq.append(seq[i-1]); aref.append(ref[j-1]); i-=1; j-=1
        elif i>0 and (j==0 or tb[i][j]==1): aseq.append(seq[i-1]); aref.append('-'); i-=1
        else: aseq.append('-'); aref.append(ref[j-1]); j-=1
    return ''.join(reversed(aseq)), ''.join(reversed(aref))

def extract_junctions(fpath):
    with open(fpath) as f: lines = f.readlines()
    loxp, lox2272 = [], []
    for i in range(len(lines)//4):
        seq = lines[i*4+1].strip().upper()
        if not (MIN_FLIP <= len(seq) <= MAX_FLIP): continue
        r = rc(seq)
        pf,_ = find_fuzzy(seq,FWD,FWD_MM); pr,_ = find_fuzzy(r,FWD,FWD_MM)
        if pf==-1 and pr==-1: continue
        if pf!=-1 and (pr==-1 or pf<=pr): norm,fwd_pos = seq,pf
        else: norm,fwd_pos = r,pr
        a_start = fwd_pos+FWD_LEN
        search  = norm[max(0,a_start-5):a_start+ANCHOR_LEN+15]
        rel,_   = find_fuzzy(search,ANCHOR,ANCHOR_MM)
        if rel==-1: continue
        win_start = max(0,a_start-5)+rel+ANCHOR_LEN
        window    = norm[win_start:win_start+60]
        sl = kmer_score(window,PROBE_LOXP,K); s2 = kmer_score(window,PROBE_LOX2272,K)
        if sl>s2: loxp.append(window)
        elif s2>sl: lox2272.append(window)
    return loxp, lox2272

# ── Collect data ──────────────────────────────────────────────────────────────
FLIP_FILES = {
    ('Neg',          1): 'KQG442_2_c51-flip-1.fastq',
    ('Neg',          2): 'KQG442_4_c51-flip-2.fastq',
    ('Agro\n(+ctrl)',1): 'KQG442_6_c51c39-flip-1.fastq',
    ('Agro\n(+ctrl)',2): 'KQG442_8_c51c39-flip-2.fastq',
    ('PVC\n(normal)',1): 'KQG442_10_c30c41-flip-1.fastq',
    ('PVC\n(normal)',2): 'KQG442_12_c30c41-flip-2.fastq',
    ('PVC\npFLS2',   1): 'KQG442_14_pFLS2-flip-1.fastq',
    ('PVC\npFLS2',   3): 'KQG442_16_pFLS2-flip-3.fastq',
    ('PVC\npUBQ10',  2): 'KQG442_18_pUBQ10-flip-2.fastq',
    ('PVC\npUBQ10',  3): 'KQG442_20_pUBQ10-flip-3.fastq',
}

CONDITIONS = ['Neg', 'Agro\n(+ctrl)', 'PVC\n(normal)', 'PVC\npFLS2', 'PVC\npUBQ10']

print('Collecting junction data (this may take ~1 min)...')
counts    = {c: {'loxp':[], 'lox2272':[]} for c in CONDITIONS}
pool_lox2 = {c: [] for c in CONDITIONS if c != 'Neg'}

for (cond, rep), fname in FLIP_FILES.items():
    lp, l2 = extract_junctions(os.path.join(FASTQ_DIR, fname))
    counts[cond]['loxp'].append(len(lp))
    counts[cond]['lox2272'].append(len(l2))
    if cond != 'Neg':
        pool_lox2[cond].extend(l2)
    print(f'  done: {cond.replace(chr(10)," ")} rep{rep}  loxP={len(lp)}  lox2272={len(l2)}')

print('Computing per-position identity...')
rng = random.Random(42)
pos_id = {}   # condition → list of 60 identity values
for cond, windows in pool_lox2.items():
    sample = rng.sample(windows, min(2000, len(windows)))
    match  = [0]*60; total = [0]*60
    for w in sample:
        aw, ar = nw_align(w, PROBE_LOX2272)
        ri = 0
        for ac, rc_ in zip(aw, ar):
            if rc_ != '-' and ri < 60:
                total[ri] += 1
                if ac == rc_: match[ri] += 1
                ri += 1
    pos_id[cond] = [100*m/t if t else 0 for m,t in zip(match,total)]
    print(f'  {cond.replace(chr(10)," ")}: avg identity {sum(pos_id[cond])/60:.1f}%')

# ── FIGURE 1 — Total flipped read counts (lox2272 + loxP) ────────────────────
apply_cowplot(14)
fig1, ax = plt.subplots(figsize=(5.5, 4.5))
theme_half_open(ax)

DOT   = '#555555'
BAR   = '#888888'
NEG_D = '#AAAAAA'
jitter_w = 0.12
rng2 = np.random.default_rng(42)

x = np.arange(len(CONDITIONS))
for xi, cond in enumerate(CONDITIONS):
    vals = np.array([l + p for l, p in zip(counts[cond]['lox2272'], counts[cond]['loxp'])], dtype=float)
    mean = vals.mean()
    dc = NEG_D if cond == 'Neg' else DOT
    bc = NEG_D if cond == 'Neg' else BAR

    ax.bar(xi, mean, width=0.45, color=bc, alpha=0.20, zorder=1, linewidth=0)
    ax.plot([xi-0.225, xi+0.225], [mean, mean], color=bc, lw=2.0, zorder=3, solid_capstyle='butt')
    jitter = rng2.uniform(-jitter_w, jitter_w, size=len(vals))
    ax.scatter(xi+jitter, vals, color=dc, s=52, zorder=4, edgecolors='white', linewidths=0.6)

    if cond == 'Neg':
        ax.annotate('rep1*', xy=(xi+jitter[0], vals[0]), xytext=(xi+0.28, vals[0]),
                    fontsize=9, color='#888888', va='center')

ax.set_xticks(x)
ax.set_xticklabels(CONDITIONS, fontsize=11)
ax.set_xlim(-0.6, len(CONDITIONS)-0.4)
ax.set_ylabel('Total flipped reads', fontsize=14)
ax.set_ylim(0, 5200)
ax.spines['left'].set_bounds(0, 5000)
ax.set_yticks([0, 1000, 2000, 3000, 4000, 5000])
ax.yaxis.set_tick_params(left=True)
ax.xaxis.set_tick_params(bottom=True)

out1 = os.path.join(OUT_DIR, 'Figure 1 — Total flipped read counts per condition.png')
fig1.savefig(out1, dpi=300, bbox_inches='tight', facecolor='white')
fig1.savefig(_svg(out1), bbox_inches='tight', facecolor='white')
plt.close(fig1)
print(f'Saved: {out1}')

# ── FIGURE 2 — Per-position junction identity ─────────────────────────────────
COND_COLORS = {
    'Agro\n(+ctrl)': '#555555',
    'PVC\n(normal)': '#888888',
    'PVC\npFLS2':    '#333333',
    'PVC\npUBQ10':   '#AAAAAA',
}
COND_LABELS = {
    'Agro\n(+ctrl)': 'Agro (+ctrl)',
    'PVC\n(normal)': 'PVC (normal)',
    'PVC\npFLS2':    'PVC pFLS2',
    'PVC\npUBQ10':   'PVC pUBQ10',
}
COND_LS = {
    'Agro\n(+ctrl)': '-',
    'PVC\n(normal)': '--',
    'PVC\npFLS2':    '-.',
    'PVC\npUBQ10':   ':',
}

apply_cowplot(13)
fig2, ax2 = plt.subplots(figsize=(7.0, 3.5))
theme_half_open(ax2)

positions = np.arange(1, 61)
for cond in ['Agro\n(+ctrl)', 'PVC\n(normal)', 'PVC\npFLS2', 'PVC\npUBQ10']:
    ax2.plot(positions, pos_id[cond],
             color=COND_COLORS[cond], lw=1.6,
             linestyle=COND_LS[cond],
             label=COND_LABELS[cond])

ax2.set_xlabel('Junction position (bp)', fontsize=13)
ax2.set_ylabel('Identity to reference (%)', fontsize=13)
ax2.set_xlim(0.5, 60.5)
ax2.set_ylim(90, 101)
ax2.set_yticks([90, 92, 94, 96, 98, 100])
ax2.spines['left'].set_bounds(90, 100)
ax2.axhline(100, color='#CCCCCC', lw=0.8, zorder=0)

# reference sequence tick labels every 10 bp
xtick_pos = [1, 10, 20, 30, 40, 50, 60]
ax2.set_xticks(xtick_pos)
ax2.legend(frameon=False, fontsize=10, loc='lower left')

out2 = os.path.join(OUT_DIR, 'Figure 2 — Per-position junction identity to reference.png')
fig2.savefig(out2, dpi=300, bbox_inches='tight', facecolor='white')
fig2.savefig(_svg(out2), bbox_inches='tight', facecolor='white')
plt.close(fig2)
print(f'Saved: {out2}')

# ── FIGURE 3 — stacked bar composition (non-Neg, per replicate) ──────────────
NON_NEG = ['Agro\n(+ctrl)', 'PVC\n(normal)', 'PVC\npFLS2', 'PVC\npUBQ10']
REP_LABELS = {
    'Agro\n(+ctrl)': ['Agro r1', 'Agro r2'],
    'PVC\n(normal)': ['PVC r1',  'PVC r2'],
    'PVC\npFLS2':    ['pFLS2 r1','pFLS2 r3'],
    'PVC\npUBQ10':   ['pUBQ10 r2','pUBQ10 r3'],
}

apply_cowplot(12)
fig3, ax3 = plt.subplots(figsize=(6.0, 3.5))
theme_half_open(ax3)

all_labels, pct_l2, pct_lp = [], [], []
for cond in NON_NEG:
    for ri, (nl, n2) in enumerate(zip(counts[cond]['loxp'], counts[cond]['lox2272'])):
        tot = nl + n2
        if tot == 0: continue
        all_labels.append(REP_LABELS[cond][ri])
        pct_l2.append(100*n2/tot)
        pct_lp.append(100*nl/tot)

xi = np.arange(len(all_labels))
ax3.bar(xi, pct_l2, color='#5e3c99', label='lox2272-flip')
ax3.bar(xi, pct_lp, bottom=pct_l2, color='#e66101', label='loxP-flip')
ax3.set_xticks(xi)
ax3.set_xticklabels(all_labels, rotation=30, ha='right', fontsize=10)
ax3.set_ylabel('% of classified flip reads', fontsize=12)
ax3.set_ylim(0, 105)
ax3.set_yticks([0, 25, 50, 75, 100])
ax3.spines['left'].set_bounds(0, 100)
ax3.legend(frameon=False, fontsize=10, loc='center left', bbox_to_anchor=(1.01, 0.5))

out3 = os.path.join(OUT_DIR, 'Figure 3 — Flip read composition per replicate.png')
fig3.savefig(out3, dpi=300, bbox_inches='tight', facecolor='white')
fig3.savefig(_svg(out3), bbox_inches='tight', facecolor='white')
plt.close(fig3)
print(f'Saved: {out3}')

# build per-replicate percentage arrays for dot plots
pct_lox2 = {c: [] for c in NON_NEG}
pct_loxp = {c: [] for c in NON_NEG}
for cond in NON_NEG:
    for nl, n2 in zip(counts[cond]['loxp'], counts[cond]['lox2272']):
        tot = nl + n2
        pct_lox2[cond].append(100 * n2 / tot if tot else 0)
        pct_loxp[cond].append(100 * nl / tot if tot else 0)

# ── FIGURE 4 — lox2272-flip% dot + bar, y from 0–100% ───────────────────────
apply_cowplot(13)
fig4, ax4 = plt.subplots(figsize=(5.0, 4.0))
theme_half_open(ax4)
rng4 = np.random.default_rng(42)
x4 = np.arange(len(NON_NEG))
for xi, cond in enumerate(NON_NEG):
    vals = np.array(pct_lox2[cond])
    mean = vals.mean()
    ax4.bar(xi, mean, width=0.45, color='#5e3c99', alpha=0.20, zorder=1, linewidth=0)
    ax4.plot([xi-0.225, xi+0.225], [mean, mean],
             color='#5e3c99', lw=2.2, zorder=3, solid_capstyle='butt')
    jitter = rng4.uniform(-0.10, 0.10, size=len(vals))
    ax4.scatter(xi+jitter, vals, color='#5e3c99', s=52, zorder=4,
                edgecolors='white', linewidths=0.6)
ax4.set_xticks(x4)
ax4.set_xticklabels(NON_NEG, fontsize=11)
ax4.set_xlim(-0.6, len(NON_NEG)-0.4)
ax4.set_ylabel('lox2272-flip (%)', fontsize=13)
ax4.set_ylim(0, 105)
ax4.set_yticks([0, 25, 50, 75, 100])
ax4.spines['left'].set_bounds(0, 100)
ax4.yaxis.set_tick_params(left=True)
ax4.xaxis.set_tick_params(bottom=True)

out4 = os.path.join(OUT_DIR, 'Figure 4 — lox2272-flip percentage per condition.png')
fig4.savefig(out4, dpi=300, bbox_inches='tight', facecolor='white')
fig4.savefig(_svg(out4), bbox_inches='tight', facecolor='white')
plt.close(fig4)
print(f'Saved: {out4}')

# ── FIGURE 5 — loxP-flip% dot + bar (y from 0, honest) ──────────────────────
apply_cowplot(13)
fig5, ax5 = plt.subplots(figsize=(5.0, 4.0))
theme_half_open(ax5)
rng5 = np.random.default_rng(42)
x5 = np.arange(len(NON_NEG))
for xi, cond in enumerate(NON_NEG):
    vals = np.array(pct_loxp[cond])
    mean = vals.mean()
    ax5.bar(xi, mean, width=0.45, color='#e66101', alpha=0.20, zorder=1, linewidth=0)
    ax5.plot([xi-0.225, xi+0.225], [mean, mean],
             color='#e66101', lw=2.2, zorder=3, solid_capstyle='butt')
    jitter = rng5.uniform(-0.10, 0.10, size=len(vals))
    ax5.scatter(xi+jitter, vals, color='#e66101', s=52, zorder=4,
                edgecolors='white', linewidths=0.6)
ax5.set_xticks(x5)
ax5.set_xticklabels(NON_NEG, fontsize=11)
ax5.set_xlim(-0.6, len(NON_NEG)-0.4)
ax5.set_ylabel('loxP-flip (%)', fontsize=13)
ax5.set_ylim(0, 1.5)
ax5.set_yticks([0, 0.5, 1.0, 1.5])
ax5.spines['left'].set_bounds(0, 1.5)
ax5.yaxis.set_tick_params(left=True)
ax5.xaxis.set_tick_params(bottom=True)

out5 = os.path.join(OUT_DIR, 'Figure 5 — loxP-flip percentage per condition.png')
fig5.savefig(out5, dpi=300, bbox_inches='tight', facecolor='white')
fig5.savefig(_svg(out5), bbox_inches='tight', facecolor='white')
plt.close(fig5)
print(f'Saved: {out5}')

print('\nAll figures done.')
