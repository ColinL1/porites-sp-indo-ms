#!/usr/bin/env python3
"""Root the mitochondrial ML tree on the given sample and report the target's sister group."""
import argparse
import dendropy

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--tree", required=True)
ap.add_argument("--root", required=True)
ap.add_argument("--target_name", required=True)
a = ap.parse_args()
t = dendropy.Tree.get(path=a.tree, schema="newick", preserve_underscores=True)
t.reroot_at_edge(t.find_node_with_taxon_label(a.root).edge, update_bipartitions=True)
print(t.as_string(schema="newick").strip())
print(t.as_ascii_plot(plot_metric="length"))
node = t.find_node_with_taxon_label(a.target_name)
parent = node.parent_node
sister = sorted(l.taxon.label for l in parent.leaf_iter() if l.taxon.label != a.target_name)
support = parent.label
if parent.edge.length is not None and parent.edge.length <= 1e-5 or support is None:
    print(f"{a.target_name} sister group: unresolved (closest: {', '.join(sister)}, joined by a zero-length "
          f"branch without bootstrap support)")
else:
    print(f"{a.target_name} sister group: {', '.join(sister)} (UFBoot {support})")
