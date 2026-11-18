#!/usr/bin/env python3
"""gbk_to_proteinortho.py <outdir> <gbk>...

Writes the protein FASTA files `funannotate compare` hands to ProteinOrtho (<name>.faa, <name>.gff,
<name>.transcripts.fa) and prints the names, one per line, in input order. Same calls and the same
naming as funannotate 1.8.17 compare.py: <organism>[_<isolate>] of the GenBank file, a "-1" suffix on
a repeated name. Run inside the funannotate image.
"""
import os
import sys

import funannotate.library as lib

outdir, gbks = sys.argv[1], sys.argv[2:]
lib.setupLogging("gbk_to_proteinortho.log")   # library functions log through lib.log
os.makedirs(outdir, exist_ok=True)
seen = []
for gbk in gbks:
    stats = lib.genomeStats(gbk)
    name = stats[0].replace(" ", "_") + ("_" + stats[1] if stats[1] else "")
    if name in seen:
        name = name + "-1"
    seen.append(name)
    lib.gb2proteinortho(gbk, outdir, name)
    print(name)
