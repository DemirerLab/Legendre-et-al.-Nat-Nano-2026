import sys, io, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter

FWD    = 'CCACTGACGTAAGGGATGACGCACAATC'
ANCHOR = 'CCACTATCCTTCGCAAGACCCTTCCTCTATATAAGGAAGTTCATTTCATTTGGAGAGG'
MOTIF_LOX2272 = 'CGCTGAC'; MOTIF_MM = 1; MOTIF_SLACK = 5
FWD_MM = 4; ANCHOR_MM = 8; MIN_FLIP, MAX_FLIP = 480, 900
ANCHOR_SHOW  = 20
JUNCTION_LEN = 164

REF_JUNCTION = ('ACACGCTGACGGCCGC'
                'GTATTTTTACAACAATTACCAACAACAACAAACAACAAACAACATTACAATTACTATTTACAATTACA'
                'GCGGCCGCCCCGG'
                'AATAACTTCGTATAGGATACTTTATACGAAGTTAT'
                'CCCGCCAAAAAA'
                'ATGGCTGCAGCCAAGGGCGA')
REF_SEQ = ANCHOR[-ANCHOR_SHOW:] + REF_JUNCTION
SEQ_LEN = len(REF_SEQ)

BASE    = os.environ.get('FLEX_FASTQ_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fastq'))
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
SVG_DIR = os.path.join(OUT_DIR, 'svg_no_ubq10')
os.makedirs(SVG_DIR, exist_ok=True)

CONDITIONS = [
    ('Agro (+ctrl)',['KQG442_6_c51c39-flip-1.fastq',  'KQG442_8_c51c39-flip-2.fastq']),
    ('PVC (normal)',['KQG442_10_c30c41-flip-1.fastq', 'KQG442_12_c30c41-flip-2.fastq']),
    ('PVC pFLS2',   ['KQG442_14_pFLS2-flip-1.fastq',  'KQG442_16_pFLS2-flip-3.fastq']),
]

NT_BG = {'A': '#8dd3c7', 'T': '#fb8072', 'G': '#ffffb3',
          'C': '#bebada', '-': '#D9D9D9'}
NT_FG = {'A': '#000000', 'T': '#000000', 'G': '#000000',
          'C': '#000000', '-': '#888888'}

ANNOTATIONS = [
    (0,   19,  '35S promoter', '#555555'),
    (20,  35,  '',             '#aaaaaa'),
    (36,  103, 'TMV UTR',      '#555555'),
    (104, 116, '',             '#aaaaaa'),
    (117, 151, 'lox2272',      '#555555'),
    (152, 163, '',             '#aaaaaa'),
    (164, 183, '2×YFP',        '#555555'),
]

def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

def find_fuzzy(seq, primer, max_mm):
    plen = len(primer); best_pos, best_mm = -1, max_mm + 1
    for i in range(len(seq) - plen + 1):
        mm = sum(a != b for a, b in zip(seq[i:i+plen], primer))
        if mm < best_mm:
            best_mm, best_pos = mm, i
            if mm == 0: break
    return (best_pos, best_mm) if best_mm <= max_mm else (-1, -1)

def get_junctions(fpath):
    with open(fpath) as f:
        lines = f.readlines()
    seqs = []
    for i in range(len(lines) // 4):
        seq = lines[i*4+1].strip().upper()
        if not (MIN_FLIP <= len(seq) <= MAX_FLIP): continue
        r = rc(seq)
        pf, _ = find_fuzzy(seq, FWD, FWD_MM)
        pr, _ = find_fuzzy(r,   FWD, FWD_MM)
        if pf == -1 and pr == -1: continue
        if pf != -1 and (pr == -1 or pf <= pr): norm, fwd_pos = seq, pf
        else:                                    norm, fwd_pos = r,   pr
        a_start = fwd_pos + len(FWD)
        search_anc = norm[max(0, a_start-5): a_start + len(ANCHOR) + 15]
        rel, _ = find_fuzzy(search_anc, ANCHOR, ANCHOR_MM)
        if rel == -1: continue
        win_start = max(0, a_start-5) + rel + len(ANCHOR)
        ext = norm[max(0, win_start-MOTIF_SLACK): win_start+MOTIF_SLACK+len(MOTIF_LOX2272)]
        if not any(sum(a!=b for a,b in zip(ext[j:], MOTIF_LOX2272)) <= MOTIF_MM
                   for j in range(len(ext)-len(MOTIF_LOX2272)+1)):
            continue
        anc_tail = norm[win_start - ANCHOR_SHOW: win_start]
        junc     = norm[win_start: win_start + JUNCTION_LEN]
        if len(anc_tail) == ANCHOR_SHOW and len(junc) == JUNCTION_LEN:
            seqs.append(anc_tail + junc)
    return seqs

CELL_W  = 0.165
CELL_H  = 0.28
FONT_SZ = 6.5
LABEL_W = 1.6
ANNOT_H = 0.55
GAP_H   = 0.08

rows = []
for cond_name, files in CONDITIONS:
    all_seqs = []
    for fname in files:
        fpath = os.path.join(BASE, fname)
        all_seqs.extend(get_junctions(fpath))
    if not all_seqs:
        print(f'WARNING: no reads for {cond_name}')
        continue
    counter = Counter(all_seqs)
    most_common_seq, count = counter.most_common(1)[0]
    rows.append((cond_name, most_common_seq))
    print(f'{cond_name}: most common = {most_common_seq}  ({count}/{len(all_seqs)})')

n_rows = 1 + len(rows)
fig_w  = LABEL_W + SEQ_LEN * CELL_W + 0.2
fig_h  = ANNOT_H + GAP_H + n_rows * CELL_H + 0.3

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
ax.set_xlim(0, fig_w)
ax.set_ylim(0, fig_h)
ax.axis('off')

seq_x0  = LABEL_W
seq_y0  = fig_h - ANNOT_H - GAP_H - CELL_H

annot_y = fig_h - ANNOT_H + 0.05
bar_h   = 0.13
for (a_start, a_end, a_label, a_col) in ANNOTATIONS:
    x0 = seq_x0 + a_start * CELL_W
    x1 = seq_x0 + (a_end + 1) * CELL_W
    rect = mpatches.FancyBboxPatch((x0 + 0.01, annot_y), x1 - x0 - 0.02, bar_h,
                                    boxstyle='round,pad=0.01',
                                    facecolor=a_col, edgecolor='none', alpha=0.85)
    ax.add_patch(rect)
    cx = (x0 + x1) / 2
    ax.text(cx, annot_y + bar_h/2, a_label,
            ha='center', va='center', fontsize=6.5, fontweight='bold',
            color='white', fontfamily='sans-serif')
    for xi in [x0 + 0.01, x1 - 0.01]:
        ax.plot([xi, xi], [annot_y, seq_y0 + CELL_H],
                color=a_col, lw=0.6, alpha=0.5, linestyle='--')

def draw_row(y_top, label, seq, is_ref=False):
    ax.text(LABEL_W - 0.06, y_top + CELL_H/2, label,
            ha='right', va='center', fontsize=6, fontfamily='monospace',
            color='#222222')
    for col, nt in enumerate(seq):
        x = seq_x0 + col * CELL_W
        bg = NT_BG.get(nt, '#FFFFFF')
        fg = NT_FG.get(nt, 'black')
        mismatch = (not is_ref) and (nt != REF_SEQ[col])
        edge_col = '#CC0000' if mismatch else 'none'
        edge_lw  = 1.2       if mismatch else 0
        rect = mpatches.FancyBboxPatch((x + 0.005, y_top + 0.01),
                                        CELL_W - 0.01, CELL_H - 0.02,
                                        boxstyle='round,pad=0.01',
                                        facecolor=bg, edgecolor=edge_col,
                                        linewidth=edge_lw)
        ax.add_patch(rect)
        ax.text(x + CELL_W/2, y_top + CELL_H/2, nt,
                ha='center', va='center', fontsize=FONT_SZ,
                fontfamily='monospace', fontweight='bold' if is_ref else 'normal',
                color=fg)

draw_row(seq_y0, 'Reference', REF_SEQ, is_ref=True)

for ri, (label, seq) in enumerate(rows):
    y = seq_y0 - (ri + 1) * CELL_H
    draw_row(y, label, seq, is_ref=False)

tick_y = seq_y0 + CELL_H + 0.01
for col in range(0, SEQ_LEN, 10):
    x = seq_x0 + col * CELL_W + CELL_W/2
    ax.text(x, tick_y, str(col), ha='center', va='bottom',
            fontsize=5, color='#666666', fontfamily='monospace')

plt.tight_layout(pad=0.1)
out = os.path.join(OUT_DIR, 'Figure 6 — lox2272-flip junction alignment.png')
plt.savefig(out, dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(SVG_DIR, os.path.basename(out).replace('.png', '.svg')), bbox_inches='tight', facecolor='white')
print('Saved:', out)
