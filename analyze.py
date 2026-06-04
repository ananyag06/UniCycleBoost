#!/usr/bin/env python3
"""
analyze.py — parse and summarize UniCycleBoost output files.

Usage:
    python analyze.py <output_file.txt> <K>
    python analyze.py results/search_output_1.txt 4

Prints a summary table and saves a plot of the score distribution.
"""

import sys
import os
import math
from pathlib import Path
from collections import Counter

try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False

N = 33  # must match constants.jl

# ── known bounds ──────────────────────────────────────────────────────────────

def known_bound(k, n=N):
    if k == 3:
        return n * n / 4
    if k == 4:
        return 96 if n == 33 else 0.5 * n ** 1.5
    if k == 5:
        return 0.5 * n ** 1.5
    if k == 6:
        return 0.5 * n ** (4 / 3)
    if k == 7:
        return 0.5 * n ** (4 / 3)
    t = k // 2
    return 0.5 * n ** (1 + 1 / t)

# ── graph string utilities ────────────────────────────────────────────────────

def edge_count(s):
    return s.count('1')

def parse_adj(s, n=N):
    adj = [[0] * (n + 1) for _ in range(n + 1)]
    chars = [c for c in s if c in ('0', '1')]
    idx = 0
    for i in range(1, n):
        for j in range(i + 1, n + 1):
            if idx >= len(chars):
                return adj
            if chars[idx] == '1':
                adj[i][j] = 1
                adj[j][i] = 1
            idx += 1
    return adj

def degree_sequence(s, n=N):
    adj = parse_adj(s, n)
    return sorted([sum(adj[i]) for i in range(1, n + 1)], reverse=True)

def avg_clustering(s, n=N):
    """Local clustering coefficient averaged over vertices with degree >= 2."""
    adj = parse_adj(s, n)
    coeffs = []
    for v in range(1, n + 1):
        nbrs = [u for u in range(1, n + 1) if adj[v][u]]
        d = len(nbrs)
        if d < 2:
            continue
        links = sum(adj[nbrs[i]][nbrs[j]]
                    for i in range(d) for j in range(i + 1, d))
        possible = d * (d - 1) / 2
        coeffs.append(links / possible)
    return sum(coeffs) / len(coeffs) if coeffs else 0.0

# ── load output file ──────────────────────────────────────────────────────────

def load_file(path):
    lines = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line)
    return lines

# ── summary table ─────────────────────────────────────────────────────────────

def print_summary(graphs, k):
    bound = known_bound(k)
    counts = sorted([edge_count(g) for g in graphs], reverse=True)
    freq   = Counter(counts)

    best   = counts[0]
    worst  = counts[-1]
    mean   = sum(counts) / len(counts)
    median = counts[len(counts) // 2]

    print()
    print(f"  ┌{'─'*52}┐")
    print(f"  │  UniCycleBoost — Results Summary                   │")
    print(f"  │  C{k}-free graphs on N={N} vertices" + " " * (34 - len(str(k))) + "│")
    print(f"  ├{'─'*52}┤")
    print(f"  │  Graphs analyzed      : {len(graphs):<28}│")
    print(f"  │  Known bound ex(n,C{k}) : {bound:<28.1f}│")
    print(f"  ├{'─'*52}┤")
    print(f"  │  Best  (edges)        : {best:<28}│")
    print(f"  │  % of known bound     : {100*best/bound:<28.1f}│")
    print(f"  │  Mean  (edges)        : {mean:<28.2f}│")
    print(f"  │  Median(edges)        : {median:<28}│")
    print(f"  │  Worst (edges)        : {worst:<28}│")
    print(f"  ├{'─'*52}┤")

    # top 5 unique scores
    print(f"  │  Score distribution (top 10 unique values):        │")
    for score, cnt in sorted(freq.items(), reverse=True)[:10]:
        bar = '█' * min(cnt, 20)
        print(f"  │    {score:>4} edges  x{cnt:<4}  {bar:<20}  │")

    print(f"  └{'─'*52}┘")
    print()

    # best graph detail
    best_graph = next(g for g in graphs if edge_count(g) == best)
    deg = degree_sequence(best_graph)
    cc  = avg_clustering(best_graph)

    print(f"  Best graph details:")
    print(f"    Edges          : {best}")
    print(f"    Max degree     : {max(deg)}")
    print(f"    Avg degree     : {sum(deg)/len(deg):.2f}")
    print(f"    Min degree     : {min(deg)}")
    print(f"    Avg clustering : {cc:.4f}")
    print(f"    Degree seq     : {deg[:12]}{'...' if len(deg)>12 else ''}")
    print()


# ── plot ──────────────────────────────────────────────────────────────────────

def plot_distribution(graphs, k, outpath=None):
    if not HAS_PLOT:
        print("  (matplotlib not installed — skipping plot)")
        return

    bound  = known_bound(k)
    counts = [edge_count(g) for g in graphs]
    freq   = Counter(counts)

    scores = sorted(freq.keys())
    cnts   = [freq[s] for s in scores]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(f"UniCycleBoost — C{k}-free on N={N}", fontsize=13, fontweight='bold')

    # left: bar chart of score distribution
    ax = axes[0]
    ax.bar(scores, cnts, color='steelblue', edgecolor='white', linewidth=0.5)
    ax.axvline(bound, color='tomato', linestyle='--', linewidth=1.5,
               label=f'Known bound ({bound:.0f})')
    ax.set_xlabel('Edge count')
    ax.set_ylabel('Number of graphs')
    ax.set_title('Score distribution')
    ax.legend(fontsize=9)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # right: cumulative % of bound
    ax2 = axes[1]
    pcts = sorted([100 * c / bound for c in counts])
    ax2.plot(range(len(pcts)), pcts, color='steelblue', linewidth=1.2)
    ax2.axhline(100, color='tomato', linestyle='--', linewidth=1.2,
                label='Known bound (100%)')
    ax2.set_xlabel('Graph index (sorted)')
    ax2.set_ylabel('% of known bound')
    ax2.set_title('Cumulative % of bound')
    ax2.legend(fontsize=9)

    plt.tight_layout()

    if outpath is None:
        outpath = f"distribution_C{k}.png"
    plt.savefig(outpath, dpi=150)
    print(f"  Plot saved to: {outpath}")
    plt.close()


# ── export CSV ────────────────────────────────────────────────────────────────

def export_csv(graphs, k, outpath=None):
    if outpath is None:
        outpath = f"results_C{k}.csv"
    bound = known_bound(k)
    with open(outpath, 'w') as f:
        f.write("rank,edges,pct_of_bound,max_degree,avg_degree,graph_string\n")
        sorted_graphs = sorted(graphs, key=edge_count, reverse=True)
        for rank, g in enumerate(sorted_graphs, 1):
            ec   = edge_count(g)
            pct  = 100 * ec / bound
            deg  = degree_sequence(g)
            maxd = max(deg)
            avgd = sum(deg) / len(deg)
            f.write(f"{rank},{ec},{pct:.2f},{maxd},{avgd:.2f},{g}\n")
    print(f"  CSV saved to: {outpath}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 3:
        print("\nUsage: python analyze.py <output_file.txt> <K>\n")
        print("  Example: python analyze.py results/search_output_1.txt 4\n")
        sys.exit(1)

    filepath = sys.argv[1]
    k        = int(sys.argv[2])

    if not Path(filepath).exists():
        print(f"\n  File not found: {filepath}\n")
        sys.exit(1)

    graphs = load_file(filepath)
    if not graphs:
        print("\n  No graphs found in file.\n")
        sys.exit(1)

    print(f"\n  Loaded {len(graphs)} graphs from {filepath}")

    print_summary(graphs, k)
    plot_distribution(graphs, k)
    export_csv(graphs, k)


if __name__ == "__main__":
    main()
