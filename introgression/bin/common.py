"""Shared helpers for the introgression scripts: the references table, paftools callable
regions and SNV projections (all in target-assembly coordinates)."""
import collections
import csv
import re

import pysam

ROLES = ("parent1", "parent2", "outgroup", "extra")
NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")   # safe in VCF, Newick, Dsuite and HyDe files


def read_references(path):
    """rows of references.tsv (name, fasta, preset, role, species, mito), in file order"""
    with open(path) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["name"] and not r["name"].startswith("#")]
    for r in rows:
        if r["role"] not in ROLES:
            raise SystemExit(f"{path}: role '{r['role']}' of {r['name']} not in {ROLES}")
        if not NAME_RE.match(r["name"]):
            raise SystemExit(f"{path}: name '{r['name']}' must match {NAME_RE.pattern}")
        r["mito"] = (r.get("mito") or "").strip()
    for role in ("parent1", "parent2"):
        if sum(r["role"] == role for r in rows) != 1:
            raise SystemExit(f"{path}: need exactly one {role}")
    if not any(r["role"] == "outgroup" for r in rows):
        raise SystemExit(f"{path}: need at least one outgroup")
    return rows


def by_role(rows, role):
    return [r["name"] for r in rows if r["role"] == role]


def load_callable(bed):
    """chrom -> sorted, merged list of [start, end) intervals"""
    d = collections.defaultdict(list)
    with open(bed) as fh:
        for line in fh:
            c, s, e = line.split()[:3]
            d[c].append((int(s), int(e)))
    for c in d:
        d[c].sort()
        m = []
        for s, e in d[c]:
            if m and s <= m[-1][1]:
                m[-1] = (m[-1][0], max(m[-1][1], e))
            else:
                m.append((s, e))
        d[c] = m
    return d


def contains(iv, p):
    """p 0-based; iv sorted list of [s, e)"""
    lo, hi = 0, len(iv)
    while lo < hi:
        m = (lo + hi) // 2
        s, e = iv[m]
        if p < s:
            hi = m
        elif p >= e:
            lo = m + 1
        else:
            return True
    return False


def load_snv(vcf):
    """chrom -> pos (1-based) -> alt base, biallelic SNVs only"""
    d = collections.defaultdict(dict)
    for r in pysam.VariantFile(vcf):
        if r.alts is None or len(r.alts) != 1:
            continue
        if len(r.ref) != 1 or len(r.alts[0]) != 1:
            continue
        if r.ref.upper() not in "ACGT" or r.alts[0].upper() not in "ACGT":
            continue
        d[r.chrom][r.pos] = r.alts[0].upper()
    return d


def read_bed_names(bed):
    with open(bed) as fh:
        return [line.split()[0] for line in fh if line.strip()]


def twisst_weights(path, p1, p2, target):
    """Read a twisst weights file (topologies in its '#topoN <newick>' header lines) and return
    an array (windows x 3) of normalised weights in the order T1 (target with parent1),
    T2 (target with parent2), T3 (parent1 with parent2: the species tree). Windows whose
    weights sum to 0 are NaN."""
    import gzip
    import numpy as np
    opener = gzip.open if path.endswith(".gz") else open
    topos, rows = [], []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                topos.append(line.split(None, 1)[1].strip())
            elif line.strip() and not line.startswith("topo"):
                rows.append([float(x) for x in line.split()])
    col = [topo_class(t, p1, p2, target) for t in topos]
    if sorted(col) != [0, 1, 2]:
        raise SystemExit(f"{path}: topologies {topos} do not map one-to-one onto T1/T2/T3")
    W = np.array(rows, float)
    out = np.full_like(W, np.nan)
    s = W.sum(1)
    ok = s > 0
    for i, k in enumerate(col):
        out[ok, k] = W[ok, i] / s[ok]
    return out


def topo_class(newick, p1, p2, target):
    """0 = (target, parent1), 1 = (target, parent2), 2 = (parent1, parent2) for an unrooted
    four-group topology such as '((A,B),C,D);'"""
    m = re.search(r"\(([^(),;]+),([^(),;]+)\)", newick.replace(" ", ""))
    if not m:
        raise SystemExit(f"cannot parse topology {newick}")
    pair = {m.group(1), m.group(2)}
    names = set(re.findall(r"[^(),;]+", newick.replace(" ", "")))
    if target in pair:
        partner = (pair - {target}).pop()
    else:
        partner = (names - pair - {target}).pop()
    return 0 if partner == p1 else 1 if partner == p2 else 2


def nj_all(taxa, dist):
    """Neighbour-joining that follows every tie. `dist` maps frozenset({a, b}) -> distance
    (use fractions.Fraction so that ties are exact). Returns one Newick string (no branch
    lengths) per distinct unrooted topology reachable by some tie resolution, in a fixed order."""
    from itertools import combinations

    def canon(splits):
        return tuple(sorted(tuple(sorted(s)) for s in splits))

    found = {}

    def leaves(node):
        return node if isinstance(node, frozenset) else frozenset([node])

    def newick(node, tree):
        if node in tree:
            return "(" + ",".join(sorted(newick(c, tree) for c in tree[node])) + ")"
        return node

    def rec(nodes, d, tree, splits):
        n = len(nodes)
        if n <= 3:
            all_leaves = frozenset(taxa)
            key = canon({s if taxa[0] not in s else all_leaves - s for s in splits})
            if key not in found:
                found[key] = "(" + ",".join(sorted(newick(x, tree) for x in nodes)) + ");"
            return
        r = {i: sum(d[frozenset((i, k))] for k in nodes if k != i) for i in nodes}
        q = {(i, j): (n - 2) * d[frozenset((i, j))] - r[i] - r[j] for i, j in combinations(nodes, 2)}
        qmin = min(q.values())
        for (i, j), v in sorted(q.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
            if v != qmin:
                continue
            u = leaves(i) | leaves(j)
            rest = [k for k in nodes if k not in (i, j)]
            d2 = dict(d)
            for k in rest:
                d2[frozenset((u, k))] = (d[frozenset((i, k))] + d[frozenset((j, k))] - d[frozenset((i, j))]) / 2
            t2 = dict(tree); t2[u] = (i, j)
            s2 = splits | ({u} if 1 < len(u) < len(taxa) - 1 else set())
            rec(rest + [u], d2, t2, s2)

    rec(list(taxa), dist, {}, frozenset())
    return [found[k] for k in sorted(found)]


def p_distances(seqs, taxa):
    """exact p-distances (Fraction) between aligned sequences, sites with a non-ACGT base skipped"""
    from fractions import Fraction
    from itertools import combinations
    d = {}
    for x, y in combinations(taxa, 2):
        pairs = [(u, v) for u, v in zip(seqs[x], seqs[y]) if u in "ACGT" and v in "ACGT"]
        d[frozenset((x, y))] = Fraction(sum(u != v for u, v in pairs), len(pairs)) if pairs else Fraction(0)
    return d
