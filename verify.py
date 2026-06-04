#!/usr/bin/env python3
"""
verify.py — check whether a graph string is truly Ck-free.

Usage:
    python verify.py                        # interactive
    python verify.py <k> <graph_string>     # direct
    python verify.py <k> <output_file.txt>  # check all graphs in a file

Useful for sanity-checking Julia output independently of the Julia code.
"""

import sys
from pathlib import Path


N = 33  # must match constants.jl


def parse_string(s, n=N):
    """Decode flat graph string → symmetric adjacency matrix (1-indexed, size n+1)."""
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


def has_simple_path(adj, u, v, length, n=N):
    """DFS: does a simple path of exactly `length` steps exist from u to v?"""
    if length == 0:
        return u == v
    visited = [False] * (n + 1)
    visited[u] = True

    def dfs(cur, remaining):
        for nb in range(1, n + 1):
            if adj[cur][nb] and not visited[nb]:
                if remaining == 1:
                    if nb == v:
                        return True
                else:
                    visited[nb] = True
                    if dfs(nb, remaining - 1):
                        visited[nb] = False
                        return True
                    visited[nb] = False
        return False

    return dfs(u, length)


def find_cycles(adj, k, n=N):
    """Find all simple Ck cycles. Returns list of cycles (each a list of vertices)."""
    cycles = []
    visited = [False] * (n + 1)

    def dfs(start, cur, depth, path):
        if depth == k:
            if adj[cur][start]:
                cycles.append(path[:])
            return
        for nb in range(1, n + 1):
            if adj[cur][nb] and not visited[nb]:
                if depth == 1 and nb <= start:
                    continue  # anchor at min vertex to avoid duplicates
                visited[nb] = True
                path.append(nb)
                dfs(start, nb, depth + 1, path)
                path.pop()
                visited[nb] = False

    for s in range(1, n + 1):
        visited[s] = True
        dfs(s, s, 1, [s])
        visited[s] = False

    return cycles


def verify(graph_string, k):
    adj = parse_string(graph_string)
    edge_count = graph_string.count('1')
    cycles = find_cycles(adj, k)
    return edge_count, cycles


def report(graph_string, k, idx=None):
    label = f"Graph {idx}" if idx is not None else "Graph"
    edge_count, cycles = verify(graph_string, k)
    status = "✓  C{}-free".format(k) if not cycles else "✗  Contains {} C{} cycle(s)".format(len(cycles), k)
    print(f"  {label}: {edge_count} edges — {status}")
    if cycles:
        print(f"    First bad cycle: {cycles[0]}")
    return len(cycles) == 0


def interactive():
    print()
    print("  UniCycleBoost — Graph Verifier")
    print()

    while True:
        try:
            k = int(input("  K (forbidden cycle length): ").strip())
            break
        except ValueError:
            print("  Enter a whole number.")

    print("  Paste graph string (or path to output .txt file):")
    raw = input("  > ").strip()

    path = Path(raw)
    if path.exists() and path.suffix == '.txt':
        lines = [l.strip() for l in path.read_text().splitlines() if l.strip()]
        print(f"\n  Checking {len(lines)} graphs for C{k}-freeness...\n")
        ok = sum(report(line, k, i + 1) for i, line in enumerate(lines))
        print(f"\n  {ok}/{len(lines)} graphs are C{k}-free.")
    else:
        print()
        report(raw, k)
    print()


def main():
    if len(sys.argv) == 1:
        interactive()
    elif len(sys.argv) == 3:
        k = int(sys.argv[1])
        raw = sys.argv[2]
        path = Path(raw)
        if path.exists():
            lines = [l.strip() for l in path.read_text().splitlines() if l.strip()]
            print(f"\nChecking {len(lines)} graphs for C{k}-freeness...\n")
            ok = sum(report(line, k, i + 1) for i, line in enumerate(lines))
            print(f"\n{ok}/{len(lines)} graphs are C{k}-free.")
        else:
            report(raw, k)
    else:
        print("Usage: python verify.py [k] [graph_string_or_file]")
        sys.exit(1)


if __name__ == "__main__":
    main()
