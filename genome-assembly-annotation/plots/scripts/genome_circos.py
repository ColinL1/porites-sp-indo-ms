#!/usr/bin/env python3
"""Circular overview of the assembly: gene, RNA, repeat and GC density in windows, near-identical gene copies.

usage: genome_circos.py style.json scaffolds.fa annotation.gff3 repeats.gff3 n_chrom out_prefix
                        [--window BP] [--unplaced] [--links possible_haplotigs.tsv]
writes out_prefix.png (300 dpi), out_prefix.pdf and out_prefix.tsv (values per window)

Sectors: scaffold_1..n_chrom; with --unplaced, all other scaffolds concatenated (in file order) into one
sector "U". Tracks, outside to inside: Mb ticks, protein-coding genes per window (heatmap), tRNA/ncRNA/rRNA
positions, repeat coverage (classified and unclassified families, each base counted once), GC content (N
excluded; windows with < 10% ACGT skipped). --links draws the gene pairs of the release's
possible_haplotigs.tsv (a gene and its near-identical copy on a chromosome-scale scaffold) in the centre.
Colours and font: setup.R via style.json.
env: conda env VTK2026_Annotate (pyCirclize 1.10.1, matplotlib 3.10.8)
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_style as fs  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
from pycirclize import Circos  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('style'); ap.add_argument('fasta'); ap.add_argument('gff'); ap.add_argument('rep_gff')
ap.add_argument('n_chrom', type=int); ap.add_argument('out')
ap.add_argument('--window', type=int, default=1_000_000)
ap.add_argument('--unplaced', action='store_true')
ap.add_argument('--links')
a = ap.parse_args()
S = fs.load(a.style)
C = S['col']
W = a.window
chroms = [f'scaffold_{i}' for i in range(1, a.n_chrom + 1)]

# sequences; with --unplaced, every other scaffold is placed at an offset inside sector "U"
seqs, offset, ublocks = {}, {}, []
name, buf, upos = None, [], 0


def flush():
    global upos
    if name is None:
        return
    s = ''.join(buf).upper()
    if name in chroms:
        seqs[name], offset[name] = s, (name, 0)
    elif a.unplaced:
        offset[name] = ('U', upos)
        ublocks.append(s)
        upos += len(s)


for line in open(a.fasta):
    if line.startswith('>'):
        flush()
        name, buf = line[1:].split()[0], []
    else:
        buf.append(line.strip())
flush()
sectors = list(chroms)
if a.unplaced:
    seqs['U'] = ''.join(ublocks)
    sectors.append('U')
size = {c: len(seqs[c]) for c in sectors}
n_unplaced = len(ublocks)


def where(scaf, pos):
    """(sector, coordinate) of a position on any scaffold, or None"""
    if scaf not in offset:
        return None
    sec, off = offset[scaf]
    return sec, off + pos


def bins(c):
    return np.arange(0, size[c], W)


gc = {}
for c in sectors:
    s, v = seqs[c], []
    for x in bins(c):
        w = s[x:x + W]
        acgt = len(w) - w.count('N')
        v.append((w.count('G') + w.count('C')) / acgt if acgt > W / 10 else np.nan)
    gc[c] = np.array(v)

gene_pos, mrna_parent, rna_pts = {}, set(), []
for line in open(a.gff):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9:
        continue
    p = where(f[0], (int(f[3]) + int(f[4])) // 2)
    if p is None:
        continue
    if f[2] == 'gene':
        gene_pos[re.search(r'ID=([^;]+)', f[8]).group(1)] = p
    elif f[2] == 'mRNA':
        mrna_parent.add(re.search(r'Parent=([^;]+)', f[8]).group(1))
    elif f[2] in ('tRNA', 'ncRNA', 'rRNA'):
        rna_pts.append((p[0], p[1], f[2]))
genes = {c: np.zeros(len(bins(c))) for c in sectors}
for g in mrna_parent:
    if g in gene_pos:
        c, x = gene_pos[g]
        genes[c][min(x // W, len(genes[c]) - 1)] += 1

known = {c: np.zeros(size[c] + 1, dtype=np.int8) for c in sectors}
unk = {c: np.zeros(size[c] + 1, dtype=np.int8) for c in sectors}
for line in open(a.rep_gff):
    if line.startswith('#'):
        continue
    f = line.split('\t')
    if len(f) < 9:
        continue
    p = where(f[0], int(f[3]) - 1)
    if p is None:
        continue
    arr = unk if 'Class=Unknown' in f[8] else known
    arr[p[0]][p[1]:p[1] + int(f[4]) - int(f[3]) + 1] = 1
cov_k, cov_u = {}, {}
for c in sectors:
    k, u = known[c], unk[c] & (1 - known[c])      # each base once, classified first
    cov_k[c] = np.array([k[x:x + W].mean() for x in bins(c)])
    cov_u[c] = np.array([u[x:x + W].mean() for x in bins(c)])

links = []
if a.links:
    for line in open(a.links):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene' or len(f) < 2 or not f[1]:
            continue
        if f[0] in gene_pos and f[1] in gene_pos:
            links.append((gene_pos[f[0]], gene_pos[f[1]]))

gmax = max(v.max() for v in genes.values())
gc_lo = np.nanmin([np.nanmin(v) for v in gc.values()])
gc_hi = np.nanmax([np.nanmax(v) for v in gc.values()])
LINK_COL = {'U': C['link']['unplaced'], 'other': C['link']['between'], 'same': C['link']['within']}
RNA_COL = C['genome_rna']
cmap = LinearSegmentedColormap.from_list('genes', C['seq_ramp'])

circos = Circos(size, space=2.0, start=0, end=360)
for sec in circos.sectors:
    c = sec.name
    x = np.minimum(bins(c) + W / 2, (bins(c) + size[c]) / 2)
    sec.text('U' if c == 'U' else c.replace('scaffold_', ''), r=108, size=9, weight='bold', color=C['ink'])
    ideo = sec.add_track((97, 100))
    ideo.axis(fc=C['other'] if c == 'U' else C['ink_muted'], ec='none')
    ideo.xticks_by_interval(10_000_000 if c != 'U' else 20_000_000, label_size=6,
                            label_orientation='vertical', label_formatter=lambda v: f'{v / 1e6:.0f}',
                            text_kws={'color': C['ink_muted']}, line_kws={'color': C['ink_muted']})
    gd = sec.add_track((88, 95))
    gd.axis(ec='none')
    gd.heatmap(genes[c].reshape(1, -1), vmin=0, vmax=gmax, cmap=cmap, rect_kws={'ec': 'none'})
    rna = sec.add_track((84, 87))
    rna.axis(ec='none')
    for (cc, p, t) in rna_pts:
        if cc == c:
            rna.line([p, p], [0, 1], vmin=0, vmax=1, lw=0.3, color=RNA_COL[t])
    rp = sec.add_track((64, 82))
    rp.axis(fc=C['band'], ec='none')
    rp.fill_between(x, cov_k[c] + cov_u[c], vmin=0, vmax=1, color=C['repeat']['unclassified'], lw=0)
    rp.fill_between(x, cov_k[c], vmin=0, vmax=1, color=C['repeat']['classified'], lw=0)
    g = sec.add_track((48, 62))
    g.axis(fc=C['band'], ec='none')
    ok = ~np.isnan(gc[c])
    g.line(x[ok], gc[c][ok], vmin=gc_lo, vmax=gc_hi, color=C['gc'], lw=0.6)


def link_kind(s1, s2):
    return 'U' if 'U' in (s1, s2) else ('same' if s1 == s2 else 'other')


for (s1, p1), (s2, p2) in links:
    kind = link_kind(s1, s2)
    circos.link((s1, p1, p1 + 1), (s2, p2, p2 + 1), r1=46, r2=46, color=LINK_COL[kind],
                alpha=0.18 if kind == 'U' else 0.45, lw=0.35)

fig = circos.plotfig(figsize=(7, 7))
handles = [Line2D([], [], color=RNA_COL['tRNA'], lw=1.2, label='tRNA'),
           Line2D([], [], color=RNA_COL['ncRNA'], lw=1.2, label='ncRNA'),
           Line2D([], [], color=RNA_COL['rRNA'], lw=1.2, label='rRNA'),
           Patch(color=C['repeat']['classified'], label='Repeats, classified families'),
           Patch(color=C['repeat']['unclassified'], label='Repeats, unclassified families'),
           Line2D([], [], color=C['gc'], lw=1.5, label='GC content')]
n_link = {k: sum(1 for (s1, _), (s2, _) in links if link_kind(s1, s2) == k) for k in LINK_COL}
lab = {'U': 'Gene copies, unplaced scaffold and chromosome', 'other': 'Gene copies, two chromosomes',
       'same': 'Gene copies, same chromosome'}
link_handles = [Line2D([], [], color=LINK_COL[k], lw=1.5, label=lab[k]) for k in LINK_COL if n_link[k]]
leg = fig.legend(handles=handles + link_handles, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=3,
                 fontsize=7.5, frameon=False, handlelength=1.6, columnspacing=1.8, labelspacing=0.45)
# colour bar of the gene heatmap, left-aligned under the legend
fig.canvas.draw()
bb = leg.get_window_extent().transformed(fig.transFigure.inverted())
cax = fig.add_axes([bb.x0 + 0.012, bb.y0 - 0.055, 0.24, 0.014])
cb = fig.colorbar(ScalarMappable(norm=Normalize(0, gmax), cmap=cmap), cax=cax, orientation='horizontal')
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=7, length=2, colors=C['ink_muted'])
cb.set_label(f'Protein-coding genes per {W / 1e6:g} Mb', size=7.5, color=C['ink'], labelpad=3)
cb.ax.xaxis.set_label_position('top')
fs.save(fig, a.out, bbox_inches='tight', pad_inches=0.05)

with open(a.out + '.tsv', 'w') as f:
    f.write('sector\tstart\tprotein_coding_genes\trepeat_classified\trepeat_unclassified\tgc\n')
    for c in sectors:
        for i, x0 in enumerate(bins(c)):
            f.write(f'{c}\t{x0}\t{int(genes[c][i])}\t{cov_k[c][i]:.4f}\t{cov_u[c][i]:.4f}\t{gc[c][i]:.4f}\n')
print(f'{len(sectors)} sectors, {sum(size.values()):,} bp; unplaced scaffolds {n_unplaced:,}; '
      f'max genes/window {int(gmax)}; GC {100 * gc_lo:.1f}-{100 * gc_hi:.1f}%; links {n_link}')
