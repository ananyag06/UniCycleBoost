# Results & Analysis — UniCycleBoost

## Why DFS beats matrix powers: the core argument

The original PatternBoost uses `A^(k-1)[u,v]` to check if adding edge `(u,v)` 
is safe. This counts the number of walks of length `k-1` from `u` to `v`.

A **walk** allows vertex revisits. A **simple path** does not.

A Cₖ requires k **distinct** vertices — so only simple paths matter.

### Example: C₄ on N=6

Suppose the current graph has edges: `1-2, 2-3, 3-1` (a triangle).

Now check if edge `(1, 4)` is safe to add (no C₄ would form).

**Matrix power check (A³[1,4]):**
Counts walks of length 3 from 1 to 4. One such walk: `1 → 2 → 1 → 2`... 
wait, 4 is not reachable here. But consider a graph where `1-2, 2-3, 3-2` 
creates a non-simple walk `1 → 2 → 3 → 2` which makes A³[1,2] > 0 even 
though no simple path of length 3 exists from 1 to 2.

**DFS check:**
Walks from 1, marks visited. Finds only true simple paths. 
Correctly identifies whether a real C₄ would form.

---

## Expected gains by cycle length

| K | Matrix power sensitivity | DFS improvement | Why |
|---|---|---|---|
| 3 | Exact (A²[u,v] correct for simple graphs) | Minimal | Short paths, few non-simple walks |
| 4 | Slightly over-conservative | Small-moderate | A³ can count non-simple walks |
| 5 | Increasingly over-conservative | Moderate | A⁴ has more non-simple walk noise |
| 6+ | Significantly over-conservative | Large | Non-simple walks dominate at large k |

**Takeaway:** UniCycleBoost's advantage grows with K. For C₅, C₆, C₇ the 
DFS approach is expected to find measurably denser constructions.

---

## Known upper bounds (OEIS / literature)

### C₄-free, N=33
Best known: **102 edges** (from OEIS A006855)

### C₅-free
ex(n, C₅) ~ ½ n^(3/2) — same asymptotic as C₄-free but exact small values 
are less well-studied, leaving more room to improve lower bounds.

### C₆-free  
ex(n, C₆) ~ ½ n^(4/3) — fewer edges allowed, constructions harder to find.

---

## How to reproduce results

```bash
# C4-free baseline (compare with OEIS A006855)
julia -t 4 search_fc.jl results/c4/ 100 500 1000 2000 4

# C5-free  
julia -t 4 search_fc.jl results/c5/ 100 500 1000 2000 5

# C6-free
julia -t 4 search_fc.jl results/c6/ 100 500 1000 2000 6
```

Results will be written to `results/c4/search_output_1.txt` etc.
The best construction and its edge count (reward) will print to stdout.

---

## What to look for

After running, compare the best reward (edge count) against:
- The known upper bound for that (N, K) pair
- The result from the original `problem_4_cycle_free.jl` at the same N

Any improvement in edge count for the same N is a valid experimental result 
demonstrating the tighter simple-path check finds better constructions.
