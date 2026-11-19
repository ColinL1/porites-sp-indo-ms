"""fig_style.py: the setup.R customisation for the Python (matplotlib / pyCirclize) figure scripts.

style.json is written by style_json.R from shared/plots/setup.R; this module only reads it:
  S = load(path)        the JSON as a dict (S['col'], S['new_genome_label'], ...)
  italic(text)          mathtext with the italic_terms of setup.R in italics ("Porites sp. Indonesia" ->
                        "$\\mathit{Porites}$ sp. Indonesia")
  save(fig, prefix)     PDF (fonts embedded as TrueType) and 300 dpi PNG
"""
import json
import re

import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt

S = {}


def load(path):
    S.update(json.load(open(path)))
    files = S['font']['files']
    for f in files.values():
        font_manager.fontManager.addfont(f)
    fam = font_manager.FontProperties(fname=files['plain']).get_name()
    plt.rcParams.update({
        'font.family': fam, 'font.size': 8, 'text.color': S['col']['ink'],
        'mathtext.fontset': 'custom', 'mathtext.rm': fam, 'mathtext.it': f'{fam}:italic',
        'mathtext.bf': f'{fam}:bold', 'pdf.fonttype': 42, 'svg.fonttype': 'none',
    })
    S['_regex'] = re.compile('(' + '|'.join(S['italic_terms']) + ')')
    return S


def italic(text):
    def mt(m):
        return r'$\mathit{' + m.group(1).replace(' ', r'\ ') + '}$'
    return S['_regex'].sub(mt, text)


def save(fig, prefix, **kw):
    fig.savefig(prefix + '.png', dpi=300, facecolor='white', **kw)
    fig.savefig(prefix + '.pdf', facecolor='white', **kw)
