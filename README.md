# UniCycleBoost

> **A unified generalized Cₖ-free graph construction engine — one score function for all cycle lengths, extending PatternBoost**

---

## Motivation

[PatternBoost](https://arxiv.org/abs/2411.00566) (Wagner et al.) is a powerful framework that alternates between local search and transformer-based global learning to find extremal graph constructions. However, its original implementation requires a **separate hand-crafted file for each problem**:

| Problem | File in original repo |
|---|---|
| Triangle-free (C₃) | `problem_triangle_free.jl` |
| Square-free (C₄) | `problem_4_cycle_free.jl` |
| Any new Cₖ | ❌ write a new file from scratch |

Each file has its own hardcoded detection logic, its own safe-edge heuristic, and its own string encoding. A mathematician who wants to explore C₅-free or C₇-free graphs must re-engineer everything manually.

**UniCycleBoost fixes this with a single unified file: `problem_ck_free.jl`.**

---

## Novelty

### 1. One function for all cycle lengths

Pass `K` as a command-line argument. Everything — detection, repair, extension, reward — adapts automatically.

```bash
# Triangle-free (C3)
julia search_fc.jl output/ 100 500 1000 2000 3

# Square-free (C4)  
julia search_fc.jl output/ 100 500 1000 2000 4

# C5-free
julia search_fc.jl output/ 100 500 1000 2000 5

# C6-free
julia search_fc.jl output/ 100 500 1000 2000 6
```

### 2. DFS-based simple-path cycle detection

The original PatternBoost files use **matrix power checks** (`A^(k-1)[u,v]`) to decide whether adding edge `(u,v)` would create a forbidden cycle. This counts **all walks** of length `k-1`, including non-simple ones that revisit vertices. The result: some safe edges are wrongly rejected.

**The key insight:**

> A Cₖ is a closed **simple** path — exactly k **distinct** vertices where vertex k+1 = vertex 1. Therefore, only simple paths (no vertex revisits) should be checked.

UniCycleBoost uses **DFS with a visited-vertex tracker**:

```
For C4: walk 4 steps from vertex u
        track every vertex visited
        never step on a visited vertex
        if step 5 = step 1 → true C4 found
```

This is the **same logic for all k** — only the depth changes. It is provably more accurate than matrix powers for k ≥ 4.

### 3. Tighter safe-edge check → denser graphs → higher rewards

Because the safe-edge check only considers simple paths, **fewer edges are wrongly rejected** during the extension phase. This allows the greedy search to build denser Cₖ-free graphs, producing higher rewards and potentially closing the gap to known upper bounds.

---

## How it works

```
Input: graph string + K (forbidden cycle length)
         │
         ▼
┌─────────────────────┐
│   PHASE 1: REPAIR   │  ← DFS finds all simple Ck cycles
│                     │    Remove most-frequent edge
│                     │    Repeat until Ck-free
└─────────────────────┘
         │
         ▼
┌─────────────────────┐
│   PHASE 2: EXTEND   │  ← For each candidate edge (u,v):
│                     │    DFS checks if simple path of
│                     │    length K-1 exists u→v
│                     │    If no → safe to add
│                     │    Repeat until no safe edges left
└─────────────────────┘
         │
         ▼
┌─────────────────────┐
│   PHASE 3: RETURN   │  ← Original + 3 permuted copies
│                     │    (diversity for transformer)
└─────────────────────┘
         │
         ▼
Output: denser Ck-free graph string
```

---

## Known bounds (for comparison with results)

| Problem | ex(n, Cₖ) lower bound | ex(n, Cₖ) upper bound |
|---|---|---|
| C₃-free | ~ n²/4 (Turán) | ~ n²/4 |
| C₄-free | ~ ½ n^(3/2) | ~ ½ n^(3/2) |
| C₅-free | ~ ½ n^(3/2) | ~ ½ n^(3/2) |
| C₆-free | ~ ½ n^(4/3) | ~ ½ n^(4/3) |

See [OEIS A006855](https://oeis.org/A006855) for exact values for C₄-free.

---

## Installation

### Prerequisites
- Julia 1.8+
- Python 3.10+

### Setup
```bash
git clone https://github.com/ananyag06/UniCycleBoost.git
cd UniCycleBoost
```

Install Julia dependencies:
```julia
using Pkg
Pkg.add(["Dictionaries", "StatsBase", "Plots", "Combinatorics"])
```

---

## Usage

```bash
julia search_fc.jl <output_dir> <nb_local_searches> <num_initial_objects> \
                   <final_db_size> <target_db_size> <K>
```

| Argument | Description |
|---|---|
| `output_dir` | Where to write results |
| `nb_local_searches` | Local searches per iteration |
| `num_initial_objects` | Starting pool size |
| `final_db_size` | Max constructions to save |
| `target_db_size` | Working database size |
| `K` | **Forbidden cycle length (3, 4, 5, 6, ...)** |

---

## Repository structure

```
UniCycleBoost/
├── problem_ck_free.jl   ← core novelty: unified Ck-free engine
├── search_fc.jl         ← patched main loop (reads K from args)
├── constants.jl         ← shared type definitions
├── README.md
└── results/
    └── analysis.md      ← bounds comparison and expected gains
```

---

## Comparison with original PatternBoost

| | Original PatternBoost | UniCycleBoost |
|---|---|---|
| Supported problems | C₃, C₄, permanent | Any Cₖ |
| Files needed per problem | 1 new file | 0 (just change K) |
| Cycle detection | Matrix powers (walks) | DFS (simple paths) |
| Safe-edge check | A^(k-1)[u,v] == 0 | No simple path of length k-1 |
| False positives in check | Yes (non-simple walks) | No |
| Diversity in output | Varies per file | Always 4 permutations |

---

## Citation

If you use this work, please cite the original PatternBoost paper:

```
@article{wagner2024patternboost,
  title={PatternBoost: Constructions in Mathematics with a Little Help from AI},
  author={Wagner, Adam Zsolt and others},
  year={2024}
}
```

---

*Built at IISc as an extension of PatternBoost.*
