using Base.Threads
using Random
using LinearAlgebra
using Statistics
import StatsBase: countmap
using Dictionaries
using Printf
using Plots
using Combinatorics
using Dates

# =============================================================================
# UniCycleBoost — Patched search_fc.jl
#
# Changes from original PatternBoost search_fc.jl:
#   1. Reads K (forbidden cycle length) from ARGS[6]
#   2. Includes problem_ck_free.jl instead of a hardcoded problem file
#   3. K is injected as a global before the problem file is loaded
#
# Usage:
#   julia search_fc.jl <write_path> <nb_local_searches> <num_initial_empty>
#                      <final_database_size> <target_db_size> <K>
#
#   Example (C4-free, N=33):
#   julia search_fc.jl output/ 100 500 1000 2000 4
# =============================================================================

# Parse K first so it is available when problem_ck_free.jl is included
const K = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 4
println("Running UniCycleBoost for C", K, "-free graphs on N=", 33, " vertices")

include("problem_ck_free.jl")
include("constants.jl")

if Threads.nthreads() > 1
    BLAS.set_num_threads(1)
end
println("Using ", nthreads(), " thread(s)")

function find_next_available_filename(base::String, extension::String)
    i = 1
    while true
        filename = @sprintf("%s/%s_%d.%s", write_path, base, i, extension)
        if !isfile(filename)
            return filename
        end
        i += 1
    end
end

function write_output_to_file(db)
    rewards = [ rew for rew in keys(db.rewards) ]
    sort!(rewards, rev=true)
    base_name = "search_output"
    extension = "txt"
    filename = find_next_available_filename(base_name, extension)
    curr_rew_index = 1
    lines_written::Int = 0
    open(filename, "w") do file
        while lines_written < final_database_size && curr_rew_index <= length(rewards)
            curr_rew = rewards[curr_rew_index]
            for obj in db.rewards[curr_rew][1:min(final_database_size - lines_written, length(db.rewards[curr_rew]))]
                write(file, obj * "\n")
            end
            lines_written += length(db.rewards[curr_rew])
            curr_rew_index += 1
        end
    end
    println("Data written to $(filename)")
    println("Best construction (reward=", rewards[1], "):")
    println(db.rewards[rewards[1]][1])
end

function write_plot_to_file(db)
    rewards = [ rew for rew in keys(db.rewards) ]
    sort!(rewards, rev=true)
    reward_counts = [ length(db.rewards[rew]) for rew in rewards ]
    bar(rewards, reward_counts, xlabel="Scores", ylabel="Count",
        title="Score Distribution (C$(K)-free, N=$(N))", legend=false)
    base_name = "plot"
    extension = "png"
    filename = find_next_available_filename(base_name, extension)
    savefig(filename)
    println("Plot saved to $(filename)")
    txt_filename = @sprintf("%s/%s.%s", write_path, "distribution", "txt")
    open(txt_filename, "w") do f
        for (rew, count) in zip(rewards, reward_counts)
            println(f, "Score: $rew, Count: $count")
        end
    end
    println("Score distribution saved to $(txt_filename)")
end

function new_db()
    return Database(
        Dictionary{OBJ_TYPE, REWARD_TYPE}(),
        Dictionary{REWARD_TYPE, Vector{OBJ_TYPE}}(),
        Dictionary{REWARD_TYPE, UInt}()
    )
end

function initial_lines()
    input_file = ""
    for arg in ARGS
        if arg == "-i" || arg == "--input"
            input_file_index = findfirst(==(arg), ARGS) + 1
            if input_file_index <= length(ARGS)
                input_file = ARGS[input_file_index]
            end
            break
        end
    end
    println("Input file: ", input_file)
    lines = String[]
    if input_file != ""
        println("Using input file")
        open(input_file, "r") do file
            for line in eachline(file)
                if length(line) == length(empty_starting_point())
                    push!(lines, line)
                end
            end
        end
    else
        println("No input file provided, starting from empty graphs")
        for _ in 1:num_initial_empty_objects
            push!(lines, empty_starting_point())
        end
    end
    return lines
end

function reward(obj)
    return reward_calc(obj)
end

function reward(db, obj)
    if haskey(db.objects, obj)
        return db.objects[obj], false
    end
    return reward(obj), true
end

function local_search_on_object(db, obj)
    objects = Vector{OBJ_TYPE}(undef, 0)
    rewards = Vector{REWARD_TYPE}(undef, 0)
    greedily_expanded_objs = greedy_search_from_startpoint(db, obj)
    for greedily_expanded_obj in greedily_expanded_objs
        rew, new = reward(db, greedily_expanded_obj)
        if new
            push!(objects, greedily_expanded_obj)
            push!(rewards, rew)
        end
    end
    return objects, rewards
end

function print_db(db)
    rewards = [ rew for rew in keys(db.rewards) ]
    sort!(rewards, rev=true)
    db_size = sum(length(db.rewards[r]) for r in rewards)
    if db_size > 2 * target_db_size
        println(" - Shrinking database to $target_db_size best objects")
        shrink!(db)
        rewards = [ rew for rew in keys(db.rewards) ]
        sort!(rewards, rev=true)
    end
    if !isempty(rewards)
        println("  Best reward so far: ", rewards[1],
                " | DB size: ", db_size,
                " | C$(K)-free on N=$(N)")
    end
end

function local_search!(db, lines, start_ind, nb=nb_local_searches)
    local_search_results_threads = [ [[], []] for _ in 1:nthreads() ]
    pool = OBJ_TYPE[]
    append!(pool, lines[start_ind:min(start_ind + nb - 1, length(lines))])
    @threads for obj in pool
        list_obj, list_rew = local_search_on_object(db, obj)
        append!(local_search_results_threads[threadid()][1], list_obj)
        append!(local_search_results_threads[threadid()][2], list_rew)
    end
    for j in 1:nthreads()
        add_db!(db, local_search_results_threads[j][1], local_search_results_threads[j][2])
    end
    return nothing
end

struct Database
    objects::Dictionary{OBJ_TYPE, REWARD_TYPE}
    rewards::Dictionary{REWARD_TYPE, Vector{OBJ_TYPE}}
    local_search_indices::Dictionary{REWARD_TYPE, UInt}
end

function add_db!(db, list_obj, list_rew=nothing)
    rewards_new_objects = []
    if list_rew !== nothing
        for i in 1:length(list_obj)
            obj = list_obj[i]
            if !haskey(db.objects, obj)
                rew = list_rew[i]
                push!(rewards_new_objects, rew)
                set!(db.objects, obj, rew)
                if !haskey(db.rewards, rew)
                    insert!(db.rewards, rew, [obj])
                    insert!(db.local_search_indices, rew, 0)
                else
                    push!(db.rewards[rew], obj)
                end
            end
        end
    else
        list_indices = Int[]
        for i in 1:length(list_obj)
            if !haskey(db.objects, list_obj[i])
                push!(list_indices, i)
            end
        end
        list_rew = zeros(Float32, length(list_obj))
        @threads for i in list_indices
            list_rew[i] = reward(list_obj[i])
        end
        for i in list_indices
            obj = list_obj[i]
            rew = list_rew[i]
            push!(rewards_new_objects, rew)
            set!(db.objects, obj, rew)
            if !haskey(db.rewards, rew)
                insert!(db.rewards, rew, [obj])
                insert!(db.local_search_indices, rew, 0)
            else
                push!(db.rewards[rew], obj)
            end
        end
    end
    return rewards_new_objects
end

function shrink!(db)
    count = 0
    rewards = [ rew for rew in keys(db.rewards) ]
    sort!(rewards, rev=true)
    for rew in rewards
        if count < target_db_size
            lg = length(db.rewards[rew])
            count += lg
            if count > target_db_size
                k = count - target_db_size
                for obj in db.rewards[rew][lg-k+1:end]
                    try delete!(db.objects, obj) catch e end
                end
                db.rewards[rew] = db.rewards[rew][1:lg-k]
                db.local_search_indices[rew] = min(db.local_search_indices[rew], lg-k)
            end
        else
            for obj in db.rewards[rew]
                unset!(db.objects, obj)
            end
            delete!(db.rewards, rew)
            delete!(db.local_search_indices, rew)
        end
    end
    return nothing
end

function main()
    db = new_db()
    lines = initial_lines()
    println("Total starting objects: ", length(lines))
    println("Unique starting objects: ", length(Set(lines)))
    println("Using ", nthreads(), " thread(s)")

    start_idx = 1
    while start_idx < length(lines)
        time_local_search = @elapsed local_search!(db, lines, start_idx)
        start_idx += nb_local_searches
        print_db(db)
    end

    print_db(db)
    write_output_to_file(db)
    write_plot_to_file(db)
end

# Parse command-line arguments
write_path           = ARGS[1]
nb_local_searches    = parse(Int, ARGS[2])
num_initial_empty_objects = parse(Int, ARGS[3])
final_database_size  = parse(Int, ARGS[4])
target_db_size       = parse(Int, ARGS[5])
# K already parsed above (ARGS[6])

main()
