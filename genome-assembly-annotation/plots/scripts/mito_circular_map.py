#!/usr/bin/env python3
"""Circular map of the final mitogenome (Liftoff annotation, GenomeMaker stage 05 / release).

usage: mito_circular_map.py style.json mitogenome.fasta mitogenome.liftoff.gff3 out_prefix
writes out_prefix.png (300 dpi), out_prefix.pdf and out_prefix.tsv (features drawn)
Tracks, outside to inside: gene names; genes (CDS exons as arrows in the direction of transcription, the
intron of a split gene as a shaded span, rRNA and tRNA as blocks); position ticks every 2 kb; GC content in
100 bp windows (20 bp step) as the deviation from the mean; GC skew (G - C) / (G + C) in the same windows.
Gene groups and the two GC tracks are keyed in the centre. Colours and font: setup.R via style.json.
env: conda env VTK2026_Annotate (pyCirclize 1.10.1, matplotlib 3.10.8)
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as fs  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from pycirclize import Circos  # noqa: E402

style, fasta, gff, out = sys.argv[1:5]
S = fs.load(style)
COL = S['col']['mito']
INK, MUTED, BAND = S['col']['ink'], S['col']['ink_muted'], S['col']['band']
GC_COL = SKEW_COL = S['col']['other']      # GC tracks in the neutral grey, skew lighter (alpha), apart from tRNA
SKEW_ALPHA = 0.5

seq = ''.join(l.strip() for l in open(fasta) if not l.startswith('>')).upper()
L = len(seq)
AA3 = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
       'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
       'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


def attr(s, key):
    m = re.search(key + r'=([^;]+)', s)
    return m.group(1) if m else ''


cds, rnas = {}, []
for line in open(gff):
    if line.startswith('#') or not line.strip():
        continue
    f = line.rstrip('\n').split('\t')
    s, e, strand = int(f[3]), int(f[4]), f[6]
    if f[2] == 'CDS':
        cds.setdefault(attr(f[8], 'gene'), []).append((s, e, strand))
    elif f[2] in ('tRNA', 'rRNA'):
        p = attr(f[8], 'product')
        name = 'trn' + AA3.get(p.replace('tRNA-', ''), p) if f[2] == 'tRNA' else \
            ('rnl' if p.startswith('16S') else 'rns')
        rnas.append((f[2], name, s, e, strand))


def group(g):
    for k in ('COX', 'ND', 'ATP'):
        if g.startswith(k):
            return k
    return 'CYTB'


rows = []
circos = Circos({'mito': L}, space=0, start=0, end=360)
sec = circos.sectors[0]
feat = sec.add_track((80, 92))
feat.axis(fc='none', ec='none')
feat.line([0, L], [0.5, 0.5], vmin=0, vmax=1, lw=0.6, color=S['col']['other'])     # backbone
lab_pc = dict(min_r=92, max_r=100, label_size=8, text_kws={'weight': 'bold', 'color': INK},
              line_kws={'color': MUTED, 'lw': 0.5})
for g, parts in cds.items():
    parts.sort()
    c = COL[group(g)]
    for (s, e, st) in parts:
        feat.arrow(s - 1, e, head_length=2, shaft_ratio=0.75, fc=c, ec='white', lw=0.4, r_lim=(81, 91))
        rows.append(('CDS', g, group(g), s, e, st))
    for (a, b) in zip(parts, parts[1:]):                                          # intron of a split gene
        feat.line([a[1], b[0] - 1], [0.5, 0.5], vmin=0, vmax=1, lw=2.2, color=c, alpha=0.45)
    # one label per exon block of a split gene, else one per gene
    for (s, e, _) in (parts if len(parts) > 1 else parts[:1]):
        feat.annotate((s + e) / 2, g, **lab_pc)
for typ, name, s, e, st in rnas:
    feat.rect(s - 1, e, r_lim=(82, 90) if typ == 'tRNA' else (81, 91), fc=COL[typ], ec='white', lw=0.3)
    rows.append((typ, name, typ, s, e, st))
    if typ == 'rRNA':
        feat.annotate((s + e) / 2, name, **lab_pc)
    else:
        feat.annotate((s + e) / 2, name, min_r=92, max_r=97, label_size=7, text_kws={'color': MUTED},
                      line_kws={'color': MUTED, 'lw': 0.4})

ticks = sec.add_track((78.5, 79))
ticks.axis(fc=S['col']['grid'], ec='none')
ticks.xticks_by_interval(2000, outer=False, label_formatter=lambda v: f'{v / 1000:.0f} kb',
                         label_size=7, tick_length=1.2, text_kws={'color': MUTED})

W, STEP = 100, 20
xs = list(range(0, L - W + 1, STEP))
gc = [(seq[x:x + W].count('G') + seq[x:x + W].count('C')) / W for x in xs]
sk = []
for x in xs:
    g_, c_ = seq[x:x + W].count('G'), seq[x:x + W].count('C')
    sk.append((g_ - c_) / (g_ + c_) if g_ + c_ else 0)
mean_gc = sum(gc) / len(gc)
pos = [x + W / 2 for x in xs]
dev = [v - mean_gc for v in gc]

gc_t = sec.add_track((58, 72))
gc_t.axis(fc=BAND, ec='none')
glim = max(abs(v) for v in dev)
gc_t.fill_between(pos, dev, vmin=-glim, vmax=glim, color=GC_COL, lw=0)
sk_t = sec.add_track((43, 57))
sk_t.axis(fc=BAND, ec='none')
lim = max(abs(v) for v in sk)
sk_t.fill_between(pos, sk, vmin=-lim, vmax=lim, color=SKEW_COL, alpha=SKEW_ALPHA, lw=0)

fig = circos.plotfig(figsize=(6, 6))
handles = [Patch(color=COL[k], label=l) for k, l in
           [('COX', 'Cytochrome c oxidase'), ('ND', 'NADH dehydrogenase'), ('CYTB', 'Cytochrome b'),
            ('ATP', 'ATP synthase'), ('rRNA', 'rRNA'), ('tRNA', 'tRNA')]]
handles += [Patch(color=GC_COL, label='GC content (deviation from mean)'), Patch(color=SKEW_COL, alpha=SKEW_ALPHA, label='GC skew')]
fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.5, 0.5), fontsize=7.5, frameon=False,
           handlelength=1.2, handleheight=0.9, labelspacing=0.45)
fs.save(fig, out, bbox_inches='tight', pad_inches=0.05)
with open(out + '.tsv', 'w') as f:
    f.write('type\tname\tgroup\tstart\tend\tstrand\n')
    for r in rows:
        f.write('\t'.join(map(str, r)) + '\n')
print(f'{L:,} bp, GC {100 * (seq.count("G") + seq.count("C")) / L:.1f}%, {len(cds)} protein-coding genes, '
      f'{sum(1 for r in rnas if r[0] == "rRNA")} rRNA, {sum(1 for r in rnas if r[0] == "tRNA")} tRNA')
