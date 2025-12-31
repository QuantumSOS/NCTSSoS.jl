# NCTSSOS Oracle: Basis generation with constraint reduction
# 
# This script generates expected unique word counts for projector/unipotent algebras.
# Run against original NCTSSOS to regenerate oracle values.
#
# Usage:
#   cd /Users/yushengzhao/projects/NCTSSOS && julia --project test/oracles/basis_counts.jl
#
# Or copy-paste the functions below into a Julia REPL with NCTSSOS loaded.

function get_ncbasis(n, d; ind=Vector{UInt16}(1:n))
    basis = [UInt16[]]
    for i = 1:d
        append!(basis, _get_ncbasis_deg(n, i, ind=ind))
    end
    return basis
end

function _get_ncbasis_deg(n, d; ind=Vector{UInt16}(1:n))
    if d > 0
        basis = Vector{UInt16}[]
        for i = 1:n
            temp = _get_ncbasis_deg(n, d-1, ind=ind)
            push!.(temp, ind[i])
            append!(basis, temp)
        end
        return basis
    else
        return [UInt16[]]
    end
end

# Projector constraint: P² = P → delete one adjacent duplicate
function constraint_reduce_projector!(word)
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, i)
        else
            i += 1
        end
    end
    return word
end

# Unipotent constraint: U² = I → delete both adjacent duplicates
function constraint_reduce_unipotent!(word)
    i = 1
    while i < length(word)
        if word[i] == word[i+1]
            deleteat!(word, i)
            deleteat!(word, i)
            i > 1 && (i -= 1)
        else
            i += 1
        end
    end
    return word
end

# Site-based commutation: _comm from NCTSSOS src/utils.jl
# - partition: first `partition` vars commute with remaining vars (placed first)
# - comm_var: first `comm_var` vars commute among themselves (sorted)
#
# KEY DIFFERENCE: NCTSSOS _comm only sorts the FIRST partition (w1).
# NCTSSoS sorts ALL vars within EACH site.
#
# Example with partition=2, comm_var=2:
#   NCTSSOS: [4,3,2,1] → [1,2,4,3] (only vars 1,2 sorted, not 3,4)
#   NCTSSoS: [4,3,2,1] → [1,2,3,4] (both sites sorted)
#
# To match NCTSSoS, use comm_both_sites below instead.
function _comm(word::Vector{UInt16}, partition, comm_var)
    if partition > 0
        w1 = copy(word[word .<= partition])
        w2 = word[word .> partition]
    else
        w1 = copy(word)
        w2 = UInt16[]
    end
    if comm_var > 0
        i = 1
        while i < length(w1)
            if w1[i] <= comm_var && w1[i+1] <= comm_var && w1[i] > w1[i+1]
                w1[i], w1[i+1] = w1[i+1], w1[i]
                if i > 1
                    i -= 1
                else
                    i = 2
                end
            else
                i += 1
            end
        end
    end
    return [w1; w2]
end

# NCTSSoS-equivalent: sort BOTH partitions (sites)
# This matches NCTSSoS behavior where all vars within each site commute.
function comm_both_sites(word::Vector{UInt16}, partition)
    w1 = sort(word[word .<= partition])
    w2 = sort(word[word .> partition])
    return [w1; w2]
end

# Generate oracle data
function generate_oracle()
    println("# NCTSSOS Oracle: Constraint Basis Counts")
    println("# Generated: $(Dates.now())")
    println("# (n, d, projector_count, unipotent_count)")
    println("const CONSTRAINT_BASIS_COUNTS_ORACLE = [")
    for n in 1:4
        entries = String[]
        for d in 0:3
            basis = get_ncbasis(n, d)
            proj = Set(constraint_reduce_projector!(copy(w)) for w in basis)
            unip = Set(constraint_reduce_unipotent!(copy(w)) for w in basis)
            push!(entries, "($n, $d, $(length(proj)), $(length(unip)))")
        end
        println("    ", join(entries, ", "), ",")
    end
    println("]")
end

# Print detailed breakdown for specific case
function print_detailed(n, d)
    basis = get_ncbasis(n, d)
    println("\nDetailed breakdown for n=$n, d=$d:")
    println("  Total raw words: $(length(basis))")
    
    proj_words = Set(constraint_reduce_projector!(copy(w)) for w in basis)
    unip_words = Set(constraint_reduce_unipotent!(copy(w)) for w in basis)
    
    println("  Unique projector-reduced: $(length(proj_words))")
    println("  Unique unipotent-reduced: $(length(unip_words))")
    
    println("\n  Projector words (sorted by length):")
    for (i, w) in enumerate(sort(collect(proj_words), by=length))
        println("    $i: $w")
    end
end

# Demonstrate _comm site-based commutation
function demo_comm()
    println("\nSite-based commutation (_comm):")
    println("  partition=2, comm_var=2:")
    println("    [3,1,2] → $(_comm(UInt16[3,1,2], 2, 2))")
    println("    [3,2,1] → $(_comm(UInt16[3,2,1], 2, 2))")
    println("    [2,3,1] → $(_comm(UInt16[2,3,1], 2, 2))")
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    using Dates
    generate_oracle()
    print_detailed(3, 3)  # The specific case: n=3, d=3 → 22
    demo_comm()
end
