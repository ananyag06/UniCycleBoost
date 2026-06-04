# =============================================================================
# problem_ck_free.jl
#
# UniCycleBoost — A unified generalized Ck-free graph construction engine
# Author: ananyag06
#
# NOVELTY:
#   The original PatternBoost framework required a separate hand-crafted file
#   for each problem (problem_triangle_free.jl for C3, problem_4_cycle_free.jl
#   for C4, etc.), each with its own hardcoded detection and search logic.
#
#   This file replaces ALL of them with a single unified function parameterized
#   by K — the forbidden cycle length. Pass K=3 for triangle-free, K=4 for
#   square-free, K=5 for C5-free, and so on. Nothing else changes.
#
#   The key insight that makes unification possible:
#   A Ck is a closed simple path of exactly k DISTINCT vertices where the
#   (k+1)th vertex = the 1st vertex. By using DFS with a visited-vertex
#   tracker, we ensure:
#     1. No vertex is revisited during the walk (true simple cycle check)
#     2. The detection logic is identical for ALL k — just the depth changes
#     3. The safe-edge check is tighter than matrix powers (no false positives
#        from non-simple walks), allowing more edges to be added
#
#   This is in contrast to matrix power checks (A^(k-1)[u,v]) used in the
#   original code, which count ALL walks including non-simple ones, making
#   the check overly conservative for k >= 4.
#
# Usage:
#   K is passed as ARGS[6] to search_fc.jl
#   e.g. julia search_fc.jl output/ 100 500 1000 2000 3  -> C3-free
#        julia search_fc.jl output/ 100 500 1000 2000 4  -> C4-free
#        julia search_fc.jl output/ 100 500 1000 2000 5  -> C5-free
#        julia search_fc.jl output/ 100 500 1000 2000 6  -> C6-free
#
# Known best constructions for reference (OEIS):
#   C3-free: ex(n,C3) ~ n²/4          (Turán graph)
#   C4-free: ex(n,C4) ~ (1/2)n^(3/2)  (OEIS A006855)
#   C5-free: ex(n,C5) ~ (1/2)n^(3/2)
#   C6-free: ex(n,C6) ~ (1/2)n^(4/3)
# =============================================================================

include("constants.jl")

# K = forbidden cycle length, injected from search_fc.jl via ARGS[6]
# Default to 4 (square-free) if running standalone
if !@isdefined(K)
    const K = 4
end

const N = 33  # number of vertices


# -----------------------------------------------------------------------------
# HELPER: ordered edge for use as dictionary key
# -----------------------------------------------------------------------------
function ordered(edge)::Tuple{Int,Int}
    return edge[1] <= edge[2] ? edge : (edge[2], edge[1])
end


# -----------------------------------------------------------------------------
# CORE NOVELTY: Simple-path DFS cycle detection
#
# Checks whether a simple path of exactly `depth` steps exists from `start`
# to `target` in the graph represented by `adjmat`, without revisiting any
# vertex. This is the unified detection logic that works for ALL k.
#
# For cycle detection:   call with target=start, depth=K
# For safe-edge check:   call with target=v,     depth=K-1 (before adding edge u->v)
#
# The visited array ensures no vertex is revisited — this is the key difference
# from matrix power checks which count non-simple walks.
# -----------------------------------------------------------------------------
function simple_path_exists(adjmat::Matrix{Int}, start::Int, target::Int,
                            depth::Int, visited::Vector{Bool})::Bool
    if depth == 0
        # Base case: are we at the target?
        return adjmat[visited |> x -> findfirst(x), target] !== nothing &&
               start == target
    end

    for neighbor in 1:N
        if adjmat[start, neighbor] == 1 && !visited[neighbor]
            # For cycle detection: at depth > 1, don't go back to original start
            # For depth == K-1 (last step before closing), we WANT to check target
            if depth == 1
                # Last step: check if we can reach target
                if neighbor == target
                    return true
                end
            else
                # Intermediate step: explore unvisited neighbors
                if neighbor != target || depth == K - 1
                    visited[neighbor] = true
                    if simple_path_exists(adjmat, neighbor, target, depth - 1, visited)
                        visited[neighbor] = false
                        return true
                    end
                    visited[neighbor] = false
                end
            end
        end
    end
    return false
end

# Cleaner wrapper used throughout — detects if a simple path of length `len`
# exists from u to v, marking u as visited so the path stays simple
function has_simple_path(adjmat::Matrix{Int}, u::Int, v::Int, len::Int)::Bool
    if len == 0
        return u == v
    end
    visited = falses(N)
    visited[u] = true
    for neighbor in 1:N
        if adjmat[u, neighbor] == 1 && !visited[neighbor]
            if len == 1
                if neighbor == v
                    return true
                end
            else
                visited[neighbor] = true
                if has_simple_path_inner(adjmat, neighbor, v, len - 1, visited)
                    visited[neighbor] = false
                    return true
                end
                visited[neighbor] = false
            end
        end
    end
    return false
end

function has_simple_path_inner(adjmat::Matrix{Int}, current::Int, target::Int,
                                remaining::Int, visited::Vector{Bool})::Bool
    if remaining == 0
        return current == target
    end
    for neighbor in 1:N
        if adjmat[current, neighbor] == 1 && !visited[neighbor]
            if remaining == 1
                if neighbor == target
                    return true
                end
            else
                visited[neighbor] = true
                if has_simple_path_inner(adjmat, neighbor, target, remaining - 1, visited)
                    visited[neighbor] = false
                    return true
                end
                visited[neighbor] = false
            end
        end
    end
    return false
end


# -----------------------------------------------------------------------------
# CYCLE DETECTION: find all Ck cycles using simple-path DFS
#
# For each starting vertex s (anchored as minimum vertex to avoid duplicates),
# DFS walks exactly K steps. If step K closes back to s through K distinct
# vertices, it is a true simple Ck cycle.
#
# This replaces find_all_triangles() and find_all_four_cycles() from the
# original files — one function handles all k.
# -----------------------------------------------------------------------------
function find_all_ck_cycles(adjmat::Matrix{Int})::Vector{Vector{Int}}
    cycles = Vector{Vector{Int}}()
    visited = falses(N)

    function dfs(start::Int, current::Int, depth::Int, path::Vector{Int})
        if depth == K
            # Check if we can close back to start (completing the cycle)
            if adjmat[current, start] == 1
                push!(cycles, copy(path))
            end
            return
        end

        for neighbor in 1:N
            if adjmat[current, neighbor] == 1 && !visited[neighbor]
                # Anchor at minimum vertex to avoid counting same cycle multiple times
                if depth == 1 && neighbor <= start
                    continue
                end
                visited[neighbor] = true
                push!(path, neighbor)
                dfs(start, neighbor, depth + 1, path)
                pop!(path)
                visited[neighbor] = false
            end
        end
    end

    for s in 1:N
        visited[s] = true
        dfs(s, s, 1, [s])
        visited[s] = false
    end

    return cycles
end


# -----------------------------------------------------------------------------
# REPAIR HEURISTIC: remove most-frequent edge across all bad cycles
#
# Same greedy strategy as original PatternBoost — unchanged because it is
# problem-agnostic. Works for any Ck.
# -----------------------------------------------------------------------------
function most_frequent_edge(cycles::Vector{Vector{Int}})::Tuple{Int,Int}
    edge_count = Dict{Tuple{Int,Int}, Int}()
    for cycle in cycles
        len = length(cycle)
        for i in 1:len
            u = cycle[i]
            v = cycle[mod1(i + 1, len)]
            e = ordered((u, v))
            edge_count[e] = get(edge_count, e, 0) + 1
        end
        # Also count the closing edge (last vertex back to first)
        e = ordered((cycle[end], cycle[1]))
        edge_count[e] = get(edge_count, e, 0) + 1
    end
    _, best = findmax(edge_count)
    return best
end


# -----------------------------------------------------------------------------
# STRING <-> ADJACENCY MATRIX
#
# Identical format to problem_4_cycle_free.jl — upper triangle entries with
# "2" as row delimiter. Kept identical so the transformer/tokenizer pipeline
# works completely unchanged.
# -----------------------------------------------------------------------------
function convert_adjmat_to_string(adjmat::Matrix{Int})::String
    entries = []
    for i in 1:N-1
        for j in i+1:N
            push!(entries, string(adjmat[i, j]))
        end
        push!(entries, "2")
    end
    return join(entries)
end

function parse_string_to_adjmat(obj::OBJ_TYPE)::Union{Matrix{Int}, Nothing}
    num_twos = count(c -> c == '2', obj)
    if num_twos != N - 1
        return nothing
    end
    adjmat = zeros(Int, N, N)
    index = 1
    for i in 1:N-1
        for j in i+1:N
            while index <= length(obj) && obj[index] == '2'
                index += 1
            end
            if index > length(obj)
                return nothing
            end
            adjmat[i, j] = parse(Int, obj[index])
            adjmat[j, i] = adjmat[i, j]
            index += 1
        end
    end
    return adjmat
end


# -----------------------------------------------------------------------------
# MAIN GREEDY SEARCH — called by search_fc.jl's local_search! loop
#
# This is the unified function. The exact same logic runs for C3, C4, C5, C6...
# Only K changes — everything else is identical.
#
# Phase 1 — REPAIR:   remove edges until no Ck exists (DFS detection)
# Phase 2 — EXTEND:   add edges that don't create a Ck (DFS safe-edge check)
# Phase 3 — RETURN:   result + permuted copies for transformer diversity
# -----------------------------------------------------------------------------
function greedy_search_from_startpoint(db, obj::OBJ_TYPE, additional_loops=0)::Vector{OBJ_TYPE}

    # Parse input string to adjacency matrix
    adjmat = parse_string_to_adjmat(obj)
    if adjmat === nothing
        return greedy_search_from_startpoint(db, empty_starting_point(), additional_loops)
    end

    # ------------------------------------------------------------------
    # PHASE 1: REPAIR
    # Remove edges until graph is Ck-free
    # Uses DFS-based simple cycle detection (the core novelty)
    # ------------------------------------------------------------------
    cycles = find_all_ck_cycles(adjmat)

    while !isempty(cycles)
        # Remove the edge that appears in the most Ck cycles
        i, j = most_frequent_edge(cycles)
        adjmat[i, j] = 0
        adjmat[j, i] = 0

        # Incrementally remove cycles that contained this edge
        cycles = filter(cycles) do cycle
            len = length(cycle)
            for idx in 1:len
                u = cycle[idx]
                v = cycle[mod1(idx + 1, len)]
                if ordered((u, v)) == ordered((i, j))
                    return false
                end
            end
            # Also check closing edge
            if ordered((cycle[end], cycle[1])) == ordered((i, j))
                return false
            end
            return true
        end
    end

    # ------------------------------------------------------------------
    # PHASE 2: EXTEND
    # Greedily add edges that won't create a Ck
    # Safe-edge check: does a simple path of length K-1 exist from u to v?
    # If yes → adding (u,v) would close a Ck → skip
    # If no  → safe to add
    #
    # This is tighter than A^(K-1)[u,v] == 0 because we only count
    # simple paths (no vertex revisits), not all walks.
    # ------------------------------------------------------------------
    allowed_edges = Vector{Tuple{Int,Int}}()
    for i in 1:N-1
        for j in i+1:N
            if adjmat[i, j] == 0 && !has_simple_path(adjmat, i, j, K - 1)
                push!(allowed_edges, (i, j))
            end
        end
    end

    while !isempty(allowed_edges)
        # Pick a random safe edge and add it
        edge = allowed_edges[rand(1:length(allowed_edges))]
        i, j = edge
        adjmat[i, j] = 1
        adjmat[j, i] = 1

        # Recheck all remaining candidate edges
        new_allowed = Vector{Tuple{Int,Int}}()
        for (x, y) in allowed_edges
            if adjmat[x, y] == 0 && !has_simple_path(adjmat, x, y, K - 1)
                push!(new_allowed, (x, y))
            end
        end
        allowed_edges = new_allowed
    end

    # ------------------------------------------------------------------
    # PHASE 3: RETURN
    # Return original + 3 randomly permuted copies
    # Permutations expose the transformer to the same construction under
    # different vertex labellings — improves generalisation
    # ------------------------------------------------------------------
    results = [convert_adjmat_to_string(adjmat)]
    for _ in 1:3
        perm = randperm(N)
        push!(results, convert_adjmat_to_string(adjmat[perm, perm]))
    end
    return results
end


# -----------------------------------------------------------------------------
# REWARD: count edges — identical for all k, no problem-specific logic needed
# -----------------------------------------------------------------------------
function reward_calc(obj::OBJ_TYPE)::REWARD_TYPE
    return count(isequal('1'), obj)
end


# -----------------------------------------------------------------------------
# EMPTY STARTING POINT: all-zeros adjacency matrix
# -----------------------------------------------------------------------------
function empty_starting_point()::OBJ_TYPE
    adjmat = zeros(Int, N, N)
    return convert_adjmat_to_string(adjmat)
end
