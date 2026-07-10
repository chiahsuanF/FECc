#!/usr/bin/env python3

import os, logging, csv
from bisect import bisect_left, bisect_right

# -----------------------
# 6. Cycle Detection
# -----------------------
def dfs_visit(node, graph, path, visited, cycles):
    if node in path:
        idx=path.index(node)
        cycles.append(path[idx:]+[node])
        return
    if node in visited:
        return
    visited.add(node)
    for nbr in graph.get(node,{}):
        dfs_visit(nbr, graph, path+[node], visited, cycles)

def find_cycles(graph, outdir, sample):
    t0 = tic("Detecting cycles")
    visited, cycles = set(), []
    for n in graph:
        if n not in visited:
            dfs_visit(n, graph, [], visited, cycles)
    fn = os.path.join(outdir, f"06_{sample}_raw_cycles.txt")
    with open(fn,'w') as f:
        f.write("cycle_nodes\n")
        for c in cycles:
            f.write(",".join(map(str,c))+"\n")
    logging.info(f"Wrote {len(cycles)} raw cycles → {fn}")
    toc(t0, "Detecting cycles")
    return cycles

# -----------------------
# Write cycles with real coords
# -----------------------
def cycle_size(cycle, regions):
    """Return total bp size of a cycle."""
    size = 0
    # ignore last node because it repeats the first node
    for idx in cycle[:-1]:
        chrom, start, end = regions[idx]
        size += (end - start + 1)
    return size


def filter_maximal_cycles(cycles):
    """
    Remove cycles that are strict subsets of another cycle.
    """
    keep = []

    node_sets = [set(c[:-1]) for c in cycles]

    for i, s1 in enumerate(node_sets):
        is_subset = False

        for j, s2 in enumerate(node_sets):
            if i == j:
                continue

            if s1 < s2:      # strict subset
                is_subset = True
                break

        if not is_subset:
            keep.append(cycles[i])

    return keep

def write_cycles_with_coords(cycles, regions, outdir, sample, label):
    cycles = filter_maximal_cycles(cycles) #
    fn = os.path.join( outdir, f"07_{sample}_{label}_cycles_coords.txt")

    with open(fn, 'w') as f:
        f.write("cycle_id\tsize\tpath\n")
        for cid, cycle in enumerate(cycles):
            size = cycle_size(cycle, regions)
            path = ";".join(
                f"{regions[i][0]}:{regions[i][1]}-{regions[i][2]}"
                for i in cycle
            )
            f.write( f"{cid}\t{size}\t{path}\n")

    logging.info(
        f"Wrote {len(cycles)} {label} cycles with coordinates → {fn}"
    )