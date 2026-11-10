#!/usr/bin/env python3
"""Write one protein FASTA per BUSCO locus that is complete and single-copy in every genome.

usage: shared_loci.py lineage outdir name=busco_dir [name=busco_dir ...]
Headers are the genome names. outdir/shared_genes.txt lists the loci.
"""
import os
import sys

lineage, outdir = sys.argv[1], sys.argv[2]
genomes = dict(a.split('=', 1) for a in sys.argv[3:])
os.makedirs(outdir, exist_ok=True)

sc = {}
for name, d in genomes.items():
    path = os.path.join(d, f'run_{lineage}', 'busco_sequences', 'single_copy_busco_sequences')
    if not os.path.isdir(path):
        sys.exit(f'ERROR: {path} not found for {name}')
    sc[name] = path
    genes = {f[:-4] for f in os.listdir(path) if f.endswith('.faa')}
    shared = genes if name == next(iter(genomes)) else shared & genes

for gene in sorted(shared):
    with open(os.path.join(outdir, f'{gene}.fasta'), 'w') as o:
        for name, path in sc.items():
            seq = ''.join(l.strip() for l in open(os.path.join(path, f'{gene}.faa')) if not l.startswith('>'))
            o.write(f'>{name}\n{seq}\n')
with open(os.path.join(outdir, 'shared_genes.txt'), 'w') as o:
    o.write('\n'.join(sorted(shared)) + '\n')
print(f'{len(shared)} loci complete and single-copy in all {len(genomes)} genomes')
