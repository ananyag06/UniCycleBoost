#!/usr/bin/env python3
"""
bounds.py — visualize known extremal bounds ex(n, Ck) and overlay results.

Usage:
    python bounds.py                          # plot bounds only
    python bounds.py <results.txt> <K>        # overlay your results too

Shows the Bondy–Simonovits / Turán bounds for k = 3..7 across a range of n,
and marks where UniCycleBoost results sit relative to the bound.
"""

import sys
import math
from pathlib import Path

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False
    print("  matplotlib / numpy not found. Install with: pip install matplotlib numpy")
    print("  Running in text-only mode.\n")


# ── bound formulas ────────────────────────────────────────────────────────────

BOUND_LABELS = {
    3: "ex(n,C₃) ~ n²/4  [Turán]",
    4: "ex(n,C₄) ~ ½n^{3/2}  [Reiman]",
    5: "ex(n,C₅) ~ ½n^{3/2}  [Bondy–Simonovits]",
    6: "ex(n,C₆) ~ ½n^{4/3}  [Bondy–Simonovits]",
    7: "ex(n,C₇) ~ ½n^{4/3}  [Bondy–Simonovits]",
}

COLORS = {3: '#e74c3c', 4: '#3498db', 5: '#2ecc71', 6: '#f39c12', 7: '#9b59b6'}

def bound(k, n):
    if k == 3:
        return n * n / 4
    if k in (4, 5):
        return 0.5 * n ** 1.5
    if k in (6, 7):
        return 0.5 * n ** (4 / 3)
    t = k // 2
    return 0.5 * n ** (1 + 1 / t)


# ── text table ────────────────────────────────────────────────────────────────

def print_bounds_table(n_values=None):
    if n_values is None:
        n_values = [10, 20, 33, 50, 100, 200, 500]

    ks = [3, 4, 5, 6, 7]

    header = f"  {'n':>6}" + "".join(f"  {'C'+str(k)+'-free':>12}" for k in ks)
    print()
    print("  Known bounds ex(n, Ck)  [upper bound, rounded]")
    print()
    print(header)
    print("  " + "─" * (7 + 14 * len(ks)))

    for n in n_values:
        row = f"  {n:>6}"
        for k in ks:
            b = bound(k, n)
            row += f"  {b:>12.1f}"
        print(row)

    print()
    print("  Formulas:")
    for k, label in BOUND_LABELS.items():
        print(f"    C{k}: {label}")
    print()


# ── ratio table (what % of bound is each result) ─────────────────────────────

def print_ratio_table(results_by_k):
    print()
    print("  Results vs known bounds:")
    print()
    print(f"  {'K':>4}  {'N':>4}  {'Edges':>8}  {'Bound':>8}  {'%':>8}  {'Gap':>8}")
    print("  " + "─" * 48)
    for k, (edges, n) in sorted(results_by_k.items()):
        b    = bound(k, n)
        pct  = 100 * edges / b
        gap  = b - edges
        print(f"  {k:>4}  {n:>4}  {edges:>8}  {b:>8.1f}  {pct:>7.1f}%  {gap:>8.1f}")
    print()


# ── load output file ──────────────────────────────────────────────────────────

def best_from_file(path):
    best = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                best = max(best, line.count('1'))
    return best


# ── plot ──────────────────────────────────────────────────────────────────────

def plot_bounds(result_points=None):
    """
    result_points: list of (k, n, edges) tuples to overlay on the bound curves.
    """
    if not HAS_PLOT:
        return

    n_range = np.linspace(5, 200, 300)
    ks      = [3, 4, 5, 6, 7]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Extremal bounds ex(n, Cₖ) — UniCycleBoost reference", fontsize=13)

    # left: absolute edge counts
    ax = axes[0]
    for k in ks:
        y = [bound(k, n) for n in n_range]
        ax.plot(n_range, y, color=COLORS[k], linewidth=1.8, label=f"C{k}")

    if result_points:
        for k, n_val, edges in result_points:
            ax.scatter([n_val], [edges], color=COLORS.get(k, 'black'),
                       s=80, zorder=5, marker='*',
                       label=f"C{k} result ({edges} edges)")

    ax.set_xlabel("n (vertices)")
    ax.set_ylabel("ex(n, Cₖ)  [edges]")
    ax.set_title("Absolute bounds")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # right: as % of C4 bound (normalized view)
    ax2 = axes[1]
    base = np.array([bound(4, n) for n in n_range])
    for k in ks:
        y    = np.array([bound(k, n) for n in n_range])
        frac = y / base
        ax2.plot(n_range, frac, color=COLORS[k], linewidth=1.8, label=f"C{k}")

    ax2.axhline(1.0, color='gray', linestyle=':', linewidth=1)
    ax2.set_xlabel("n (vertices)")
    ax2.set_ylabel("Relative to ex(n, C₄)")
    ax2.set_title("Normalized (C₄ = 1.0)")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = "bounds_overview.png"
    plt.savefig(outpath, dpi=150)
    print(f"  Plot saved to: {outpath}")
    plt.close()


def plot_single_k(k, result_edges=None, n_max=200):
    """Detailed plot for a single k, showing UniCycleBoost result if provided."""
    if not HAS_PLOT:
        return

    n_range = np.linspace(5, n_max, 400)
    y       = [bound(k, n) for n in n_range]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(n_range, y, color=COLORS.get(k, 'steelblue'), linewidth=2,
            label=BOUND_LABELS.get(k, f"ex(n, C{k})"))

    # mark n=33 specifically
    b33 = bound(k, 33)
    ax.axvline(33, color='gray', linestyle=':', linewidth=1, alpha=0.7)
    ax.scatter([33], [b33], color='gray', s=60, zorder=4,
               label=f"n=33, bound={b33:.1f}")

    if result_edges is not None:
        pct = 100 * result_edges / b33
        ax.scatter([33], [result_edges], color='tomato', s=100, zorder=5,
                   marker='*', label=f"UniCycleBoost: {result_edges} edges ({pct:.1f}%)")
        ax.vlines(33, result_edges, b33, colors='tomato', linestyles='--',
                  linewidth=1, alpha=0.6)

    ax.set_xlabel("n (vertices)")
    ax.set_ylabel("edges")
    ax.set_title(f"C{k}-free graphs: known bound vs UniCycleBoost result")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    outpath = f"bounds_C{k}.png"
    plt.savefig(outpath, dpi=150)
    print(f"  Plot saved to: {outpath}")
    plt.close()


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    result_file = sys.argv[1] if len(sys.argv) > 1 else None
    k           = int(sys.argv[2]) if len(sys.argv) > 2 else None

    print_bounds_table()

    result_points = None
    results_by_k  = {}

    if result_file and k:
        path = Path(result_file)
        if not path.exists():
            print(f"  File not found: {result_file}")
        else:
            edges = best_from_file(path)
            print(f"  Best result from {result_file}: {edges} edges (C{k}-free, N={N})")
            b   = bound(k, N)
            pct = 100 * edges / b
            print(f"  Known bound: {b:.1f}  →  {pct:.1f}% of bound\n")
            result_points = [(k, N, edges)]
            results_by_k[k] = (edges, N)
            print_ratio_table(results_by_k)
            plot_single_k(k, result_edges=edges)

    plot_bounds(result_points=result_points)


N = 33

if __name__ == "__main__":
    main()
