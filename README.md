# UniCycleBoost

An extension of [PatternBoost](https://arxiv.org/abs/2411.00566) (Wagner et al.) for constructing dense Cₖ-free graphs, where k is any cycle length you want.

---

## Background / why this exists

We were working with PatternBoost at IISc and wanted to try it on C₅-free and C₆-free graphs, not just C₃ and C₄. The problem: each forbidden pattern in the original repo is its own hardcoded file with its own detection logic. To do C₅ you'd basically copy `problem_4_cycle_free.jl`, change some numbers, and hope you got it right. That seemed fragile, so we unified everything into a single parameterized file.

While doing that, we noticed the original matrix-power check (`A^(k-1)[u,v]`) is counting *all walks* of length k-1, not just simple paths. For C₄ this barely matters. For C₅ and above, non-simple walks become more common and the check starts rejecting edges that are actually safe — meaning you end up with sparser graphs than necessary. The fix is straightforward: use DFS with a visited-vertex set instead.

---

## What changed

**One file for all cycle lengths.** Pass `K` as the last argument; the detection depth, safe-edge check, and repair phase all adapt to it automatically. No new file needed for any Cₖ.

**DFS-based cycle detection instead of matrix powers.** A Cₖ requires k *distinct* vertices — it's a simple cycle by definition. So we only check simple paths (DFS + visited set), not all walks. This means the safe-edge check has no false positives: an edge is only rejected if it would genuinely close a forbidden cycle.

**Degree-aware edge selection in the extend phase.** Instead of picking randomly from safe candidate edges, we prefer edges whose endpoints have lower combined degree. This tends to spread degree more evenly across vertices, which matches how known extremal constructions (e.g. polarity graphs) are structured.

---

## Running it

```bash
# C4-free on 33 vertices
julia search_fc.jl output/ 100 500 1000 2000 4

# C5-free
julia search_fc.jl output/ 100 500 1000 2000 5

# C6-free
julia search_fc.jl output/ 100 500 1000 2000 6
```

Arguments: `<output_dir> <local_searches> <initial_pool> <final_db_size> <target_db_size> <K>`

Julia dependencies:
```julia
using Pkg
Pkg.add(["Dictionaries", "StatsBase", "Plots", "Combinatorics"])
```

---

## Files

```
problem_ck_free.jl   # unified engine — the main change from original PatternBoost
search_fc.jl         # patched main loop (reads K from args, includes above)
constants.jl         # shared type definitions (unchanged)
results/             # outputs go here
```

---

## Known bounds

For context when reading results:

| k | ex(n, Cₖ) | reference |
|---|-----------|-----------|
| 3 | ~ n²/4 | Turán (exact) |
| 4 | ~ ½ n^{3/2} | [OEIS A006855](https://oeis.org/A006855) |
| 5 | ~ ½ n^{3/2} | Bondy–Simonovits |
| 6 | ~ ½ n^{4/3} | Bondy–Simonovits |

For n=33, C₄-free: the best known construction has 96 edges. We're trying to match or approach that.

---

## What's still missing

- Actual benchmark numbers comparing DFS vs matrix-power baseline (running now)
- Transformer integration (this is just the local search component; the full PatternBoost loop needs a trained model)
- Testing on larger n

---

## Citation

If you use this, please also cite the original PatternBoost paper:

```
@article{wagner2024patternboost,
  title={PatternBoost: Constructions in Mathematics with a Little Help from AI},
  author={Wagner, Adam Zsolt and others},
  year={2024},
  url={https://arxiv.org/abs/2411.00566}
}
```

*IISc Bangalore*
