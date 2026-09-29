import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Sequences ────────────────────────────────────────────────────────────────
FWD    = 'CCACTGACGTAAGGGATGACGCACAATC'   # 28 bp, shared by all amplicons

# 58 bp anchor immediately after FWD — identical in all three states
ANCHOR = 'CCACTATCCTTCGCAAGACCCTTCCTCTATATAAGGAAGTTCATTTCATTTGGAGAGG'

# Classification strategy (two-step):
#   Step 1 — lox2272 detection: search for motif CGCTGAC (≤1 mm) in a ±5 bp window
#            around the post-anchor junction start. This handles Nanopore errors in
#            the motif AND anchor-detection offsets of up to 5 bp. If motif found
#            → lox2272-flip.
#   Step 2 — loxP detection: if motif absent, sweep offsets −3…+3 and find the
#            minimum Hamming distance to the loxP 30-mer (≤5 mm → loxP-flip).
#            The sweep catches loxP reads where the first few junction bases carry
#            Nanopore errors, shifting the apparent best-match position by 1–3 bp.
#   Unclassified: anchor-confirmed flip reads where neither test passes.
#
# Why motif over Hamming-vs-lox2272:
#   The lox2272 and loxP 30-mers share 7/30 positions; an anchor offset of 5 bp
#   shifts the lox2272 head by 5 bp, making it look more like loxP than lox2272
#   (23 mm vs lox2272, 3 mm vs loxP). The motif CGCTGAC, searched in an extended
#   window, is immune to this offset and has no overlap with the loxP junction.
MOTIF_LOX2272 = 'CGCTGAC'   # 7-mer unique to lox2272-flip junction
MOTIF_MM      = 1           # allow 1 Nanopore error in the 7-mer
MOTIF_SLACK   = 5           # search win_start±5 bp to cover anchor offsets

HEAD_LOXP = 'ACAACGGCCGCGTATTTTTACAACAATTAC'   # 30 bp loxP-flip reference head
HEAD_LEN  = 30
LOXP_ED   = 3   # max edit distance (allows indels) to loxP head

# ── Thresholds ────────────────────────────────────────────────────────────────
FWD_MM    = 4   # Hamming mismatches for FWD primer (28 bp)
ANCHOR_MM = 8   # Hamming mismatches for anchor (58 bp, ~14% tolerance)

MIN_UNFLIP, MAX_UNFLIP = 800, 1050   # expected ~928 bp
MIN_FLIP,   MAX_FLIP   =  480,  680   # expected 579 or 584 bp

# ── Samples ───────────────────────────────────────────────────────────────────
import os
BASE = os.environ.get('FLEX_FASTQ_DIR', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fastq'))

# (condition_label, rep) → (unflip_file, flip_file)
PAIRS = {
    ('Neg',          1): ('KQG442_1_c51-unflip-1.fastq',       'KQG442_2_c51-flip-1.fastq'),
    ('Neg',          2): ('KQG442_3_c51-unflip-2.fastq',       'KQG442_4_c51-flip-2.fastq'),
    ('Agro(+ctrl)',   1): ('KQG442_5_c51c39-unflip-1.fastq',    'KQG442_6_c51c39-flip-1.fastq'),
    ('Agro(+ctrl)',   2): ('KQG442_7_c51c39-unflip-2.fastq',    'KQG442_8_c51c39-flip-2.fastq'),
    ('PVC(normal)',  1): ('KQG442_9_c30c41-unflip-1.fastq',    'KQG442_10_c30c41-flip-1.fastq'),
    ('PVC(normal)',  2): ('KQG442_11_c30c41-unflip-2.fastq',   'KQG442_12_c30c41-flip-2.fastq'),
    ('PVC-pFLS2',    1): ('KQG442_13_pFLS2-unflip-1.fastq',    'KQG442_14_pFLS2-flip-1.fastq'),
    ('PVC-pFLS2',    3): ('KQG442_15_pFLS2-unflip-3.fastq',    'KQG442_16_pFLS2-flip-3.fastq'),
    ('PVC-pUBQ10',   2): ('KQG442_17_pUBQ10-unflip-2.fastq',   'KQG442_18_pUBQ10-flip-2.fastq'),
    ('PVC-pUBQ10',   3): ('KQG442_19_pUBQ10-unflip-3.fastq',   'KQG442_20_pUBQ10-flip-3.fastq'),
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def edit_dist(a, b):
    n, m = len(a), len(b)
    dp = list(range(m + 1))
    for i in range(1, n + 1):
        prev = dp[0]; dp[0] = i
        for j in range(1, m + 1):
            tmp = dp[j]
            dp[j] = prev if a[i-1] == b[j-1] else 1 + min(prev, dp[j], dp[j-1])
            prev = tmp
    return dp[m]

def find_fuzzy(seq, primer, max_mm):
    plen = len(primer)
    best_pos, best_mm = -1, max_mm + 1
    for i in range(len(seq) - plen + 1):
        mm = 0
        for a, b in zip(seq[i:i+plen], primer):
            if a != b:
                mm += 1
                if mm >= best_mm:
                    break
        if mm < best_mm:
            best_mm, best_pos = mm, i
            if mm == 0:
                break
    return (best_pos, best_mm) if best_mm <= max_mm else (-1, -1)

FWD_LEN    = len(FWD)     # 28
ANCHOR_LEN = len(ANCHOR)  # 58

# ── Per-file analysis ─────────────────────────────────────────────────────────
def analyze(fpath, mode):
    """mode = 'unflip' or 'flip'. Returns count dict."""
    with open(fpath) as f:
        lines = f.readlines()

    total = pass_size = pass_fwd = pass_anchor = 0
    loxp = lox2272 = unclassified = 0

    lo = MIN_UNFLIP if mode == 'unflip' else MIN_FLIP
    hi = MAX_UNFLIP if mode == 'unflip' else MAX_FLIP

    for i in range(len(lines) // 4):
        seq = lines[i*4+1].strip().upper()
        total += 1

        if not (lo <= len(seq) <= hi):
            continue
        pass_size += 1

        # FWD primer on either strand — normalize to FWD strand
        r = rc(seq)
        pf, _ = find_fuzzy(seq, FWD, FWD_MM)
        pr, _ = find_fuzzy(r,   FWD, FWD_MM)
        if pf == -1 and pr == -1:
            continue
        pass_fwd += 1

        if pf != -1 and (pr == -1 or pf <= pr):
            norm, fwd_pos = seq, pf
        else:
            norm, fwd_pos = r, pr

        # anchor search in window right after FWD (allow ±15 bp slack)
        a_start = fwd_pos + FWD_LEN
        search   = norm[max(0, a_start - 5) : a_start + ANCHOR_LEN + 15]
        rel, _   = find_fuzzy(search, ANCHOR, ANCHOR_MM)
        if rel == -1:
            continue
        pass_anchor += 1

        if mode == 'unflip':
            continue  # just counting anchor-confirmed unflipped reads

        # FLIP mode: motif + edit-distance classification
        win_start = max(0, a_start - 5) + rel + ANCHOR_LEN

        # Step 1: lox2272 — CGCTGAC motif in ±MOTIF_SLACK window (≤1 mm)
        ext = norm[max(0, win_start - MOTIF_SLACK) : win_start + MOTIF_SLACK + len(MOTIF_LOX2272)]
        motif_found = any(
            sum(a != b for a, b in zip(ext[j:], MOTIF_LOX2272)) <= MOTIF_MM
            for j in range(len(ext) - len(MOTIF_LOX2272) + 1)
        )
        if motif_found:
            lox2272 += 1
            continue

        # Step 2: loxP — edit distance on 30 bp post-anchor (allows indels)
        head = norm[win_start : win_start + HEAD_LEN]
        if len(head) < HEAD_LEN:
            unclassified += 1
            continue
        if edit_dist(head, HEAD_LOXP) <= LOXP_ED:
            loxp += 1
        else:
            unclassified += 1

    return dict(total=total, pass_size=pass_size, pass_fwd=pass_fwd,
                pass_anchor=pass_anchor, loxp=loxp, lox2272=lox2272,
                unclassified=unclassified)

# ── Run & report ──────────────────────────────────────────────────────────────
print(f"{'Sample':<20} {'Rep':>3}  "
      f"{'Unflip(N)':>10}  "
      f"{'loxP-flip':>10}  {'lox2272-flip':>12}  "
      f"{'Unclass':>8}  "
      f"{'%loxP':>7}  {'%lox2272':>9}")
print('-' * 90)

for (cond, rep), (uf_file, fl_file) in PAIRS.items():
    uf = analyze(os.path.join(BASE, uf_file), 'unflip')
    fl = analyze(os.path.join(BASE, fl_file), 'flip')

    total_flip = fl['loxp'] + fl['lox2272']
    pct_loxp    = 100 * fl['loxp']   / total_flip if total_flip else 0
    pct_lox2272 = 100 * fl['lox2272']/ total_flip if total_flip else 0

    print(f"{cond:<20} {rep:>3}  "
          f"{uf['pass_anchor']:>10}  "
          f"{fl['loxp']:>10}  {fl['lox2272']:>12}  "
          f"{fl['unclassified']:>8}  "
          f"{pct_loxp:>6.1f}%  {pct_lox2272:>8.1f}%")

print()
print('Notes:')
print('  Unflip(N)    = reads with FWD + anchor confirmed (unflipped molecules)')
print('  loxP-flip    = flip reads classified as loxP-mediated inversion')
print('  lox2272-flip = flip reads classified as lox2272-mediated inversion')
print('  %loxP/%lox2272 = fraction of classified flip reads (excludes unclassified)')
print()
print('Detailed counts per file:')
print()

for (cond, rep), (uf_file, fl_file) in PAIRS.items():
    uf = analyze(os.path.join(BASE, uf_file), 'unflip')
    fl = analyze(os.path.join(BASE, fl_file), 'flip')
    label = f'{cond} rep{rep}'
    print(f'{label}:')
    print(f'  [UNFLIP] total={uf["total"]}  size={uf["pass_size"]}  FWD={uf["pass_fwd"]}  anchor={uf["pass_anchor"]}')
    print(f'  [FLIP]   total={fl["total"]}  size={fl["pass_size"]}  FWD={fl["pass_fwd"]}  anchor={fl["pass_anchor"]}  loxP={fl["loxp"]}  lox2272={fl["lox2272"]}  unclass={fl["unclassified"]}')
