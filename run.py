#!/usr/bin/env python3
"""
UniCycleBoost — interactive runner
Asks for K, runs the Julia local search, and prints results to screen.
"""

import subprocess
import sys
import os
import math
import tempfile
import shutil
from pathlib import Path

# ── known upper bounds for ex(n, Ck) at n=33 ──────────────────────────────────
# These are the best known values / tight asymptotic bounds we use for comparison.
# C4 at n=33: exact value from OEIS A006855 is 96
KNOWN_BOUNDS = {
    3: lambda n: n**2 / 4,           # Turán
    4: lambda n: 96 if n == 33 else 0.5 * n**1.5,
    5: lambda n: 0.5 * n**1.5,       # same order as C4
    6: lambda n: 0.5 * n**(4/3),
    7: lambda n: 0.5 * n**(4/3),     # same order as C6
    8: lambda n: 0.5 * n**(5/4),
}

N = 33  # number of vertices (must match constants.jl)


def bound_for(k, n=N):
    if k in KNOWN_BOUNDS:
        return KNOWN_BOUNDS[k](n)
    # Bondy–Simonovits general: ex(n, C_{2t}) ~ c * n^{1+1/t}
    # rough fallback
    t = k // 2
    return 0.5 * n ** (1 + 1/t)


def parse_output_file(path):
    """
    Read the Julia output file, return list of (edge_count, graph_string) sorted descending.
    Graph strings are 0/1/2 encoded; count of '1' chars = edge count.
    """
    results = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                edges = line.count('1')
                results.append((edges, line))
    results.sort(reverse=True)
    return results


def decode_graph(s, n=N):
    """
    Reconstruct adjacency list from the flat string encoding.
    Upper triangle entries separated by '2' as row delimiter.
    Returns list of edges as (i, j) 1-indexed.
    """
    edges = []
    chars = [c for c in s if c in ('0', '1')]
    idx = 0
    for i in range(1, n):
        for j in range(i + 1, n + 1):
            if idx >= len(chars):
                break
            if chars[idx] == '1':
                edges.append((i, j))
            idx += 1
    return edges


def degree_sequence(edges, n=N):
    deg = [0] * (n + 1)
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    return sorted(deg[1:], reverse=True)


def print_result(k, edge_count, graph_string, rank=1):
    bound = bound_for(k)
    pct = 100 * edge_count / bound if bound > 0 else 0
    edges = decode_graph(graph_string)
    deg_seq = degree_sequence(edges)
    max_deg = max(deg_seq)
    avg_deg = sum(deg_seq) / len(deg_seq)

    print()
    print(f"  {'─'*50}")
    print(f"  Result #{rank}  (C{k}-free, N={N})")
    print(f"  {'─'*50}")
    print(f"  Edges found    : {edge_count}")
    print(f"  Known bound    : {bound:.1f}  (ex({N}, C{k}))")
    print(f"  % of bound     : {pct:.1f}%")
    print(f"  Max degree     : {max_deg}")
    print(f"  Avg degree     : {avg_deg:.2f}")
    print(f"  Degree sequence: {deg_seq[:10]}{'...' if len(deg_seq) > 10 else ''}")
    print(f"  {'─'*50}")
    print(f"  Graph string   : {graph_string[:80]}{'...' if len(graph_string) > 80 else ''}")
    print()


def check_julia():
    try:
        subprocess.run(["julia", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def get_k_from_user():
    print()
    print("  ╔══════════════════════════════════════╗")
    print("  ║        UniCycleBoost Runner          ║")
    print("  ║  Ck-free graph construction engine   ║")
    print("  ╚══════════════════════════════════════╝")
    print()
    print("  Finds the densest Ck-free graph on 33 vertices")
    print("  using DFS-based cycle detection (PatternBoost extension).")
    print()

    while True:
        try:
            raw = input("  Enter K (forbidden cycle length, e.g. 3, 4, 5, 6): ").strip()
            k = int(raw)
            if k < 3:
                print("  K must be at least 3.")
                continue
            if k > 10:
                print(f"  Warning: K={k} will be slow (DFS depth {k-1}). Continue? [y/n] ", end="")
                if input().strip().lower() != 'y':
                    continue
            return k
        except ValueError:
            print("  Please enter a whole number.")


def get_search_params():
    print()
    print("  Search parameters (press Enter to use defaults):")

    def ask(prompt, default):
        val = input(f"  {prompt} [{default}]: ").strip()
        return int(val) if val else default

    local_searches   = ask("Local searches per round", 50)
    initial_pool     = ask("Initial pool size      ", 200)
    final_db_size    = ask("Final DB size          ", 500)
    target_db_size   = ask("Working DB size        ", 1000)
    return local_searches, initial_pool, final_db_size, target_db_size


def run(k, local_searches, initial_pool, final_db_size, target_db_size):
    script_dir = Path(__file__).parent
    output_dir = tempfile.mkdtemp(prefix="unicycleboost_")

    cmd = [
        "julia",
        str(script_dir / "search_fc.jl"),
        output_dir + "/",
        str(local_searches),
        str(initial_pool),
        str(final_db_size),
        str(target_db_size),
        str(k),
    ]

    print()
    print(f"  Running C{k}-free search on N={N} vertices...")
    print(f"  (this may take a minute — watch Julia output below)\n")
    print("  " + "─" * 50)

    try:
        proc = subprocess.run(cmd, cwd=script_dir)
    except KeyboardInterrupt:
        print("\n  Interrupted.")
        shutil.rmtree(output_dir, ignore_errors=True)
        sys.exit(0)

    # find output file(s)
    output_files = sorted(Path(output_dir).glob("search_output_*.txt"))
    if not output_files:
        print("\n  No output file found. Julia may have errored above.")
        shutil.rmtree(output_dir, ignore_errors=True)
        return

    print("\n  " + "─" * 50)
    print(f"  Search complete.")

    results = parse_output_file(output_files[-1])
    if not results:
        print("  Output file was empty.")
    else:
        top_n = min(3, len(results))
        print(f"\n  Top {top_n} constructions found:\n")
        for i, (edge_count, graph_string) in enumerate(results[:top_n], 1):
            print_result(k, edge_count, graph_string, rank=i)

    shutil.rmtree(output_dir, ignore_errors=True)


def main():
    if not check_julia():
        print("\n  ERROR: Julia not found. Please install Julia and make sure it's on your PATH.")
        print("  https://julialang.org/downloads/\n")
        sys.exit(1)

    k = get_k_from_user()
    local_searches, initial_pool, final_db_size, target_db_size = get_search_params()
    run(k, local_searches, initial_pool, final_db_size, target_db_size)


if __name__ == "__main__":
    main()
