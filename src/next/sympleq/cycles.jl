# =============================================================================
# SympleQ Step 4: linear relations (checked algebraically during search)
# =============================================================================

function _gf2_rref(A::AbstractMatrix{<:Integer})
    B = UInt8.(mod.(A, 2))
    nrows, ncols = size(B)
    pivots = Int[]
    row = 1

    for col in 1:ncols
        pivot = 0
        for r in row:nrows
            if B[r, col] == 1
                pivot = r
                break
            end
        end
        pivot == 0 && continue

        if pivot != row
            B[row, :], B[pivot, :] = copy(B[pivot, :]), copy(B[row, :])
        end

        for r in 1:nrows
            if r != row && B[r, col] == 1
                @inbounds for c in col:ncols
                    B[r, c] ⊻= B[row, c]
                end
            end
        end

        push!(pivots, col)
        row += 1
        row > nrows && break
    end

    return B, pivots
end

function _gf2_rank(A::AbstractMatrix{<:Integer})
    _, pivots = _gf2_rref(A)
    return length(pivots)
end

function _gf2_nullspace_basis(A::AbstractMatrix{<:Integer})
    rref, pivots = _gf2_rref(A)
    ncols = size(rref, 2)
    pivot_set = Set(pivots)
    free_cols = [col for col in 1:ncols if !(col in pivot_set)]
    basis = Vector{UInt8}[]

    for free_col in free_cols
        v = zeros(UInt8, ncols)
        v[free_col] = 1
        for (r, pivot_col) in enumerate(pivots)
            rref[r, free_col] == 1 && (v[pivot_col] = 1)
        end
        push!(basis, v)
    end

    return basis
end

function _gf2_iszero(v::AbstractVector{<:Integer})
    return all(x -> !isodd(x), v)
end

function _gf2_row_in_span(row::AbstractVector{<:Integer}, rows::AbstractMatrix{<:Integer})
    isempty(row) && return true
    size(rows, 2) == length(row) || throw(DimensionMismatch("GF(2) span check has incompatible row width."))
    isempty(rows) && return _gf2_iszero(row)
    return _gf2_rank(vcat(UInt8.(mod.(rows, 2)), reshape(UInt8.(mod.(row, 2)), 1, :))) == _gf2_rank(rows)
end

function _gf2_independent_row_indices(A::AbstractMatrix{<:Integer})
    selected = Int[]
    current = zeros(UInt8, 0, size(A, 2))
    current_rank = 0

    for i in axes(A, 1)
        candidate = vcat(current, reshape(UInt8.(mod.(A[i, :], 2)), 1, :))
        candidate_rank = _gf2_rank(candidate)
        if candidate_rank > current_rank
            push!(selected, i)
            current = candidate
            current_rank = candidate_rank
        end
    end

    return selected
end

"""
    pauli_cycle_basis(tab::SymplecticTableau; cycle_strategy=:cycle_basis)

Return cycles as vectors of Pauli-term vertex indices. The default and currently
only implemented strategy is a GF(2) cycle basis of `pᵀ`, i.e. a basis for
linear relations among tableau rows.

This algebraic basis is not invariant under term permutations and must not be
used as a set of graph vertices. The detection search checks its span instead.
"""
function pauli_cycle_basis(tab::SymplecticTableau; cycle_strategy::Symbol=:cycle_basis)
    cycle_strategy == :cycle_basis || throw(ArgumentError(
        "Unsupported SympleQ cycle strategy $cycle_strategy. Implemented: :cycle_basis."
    ))

    isempty(tab.paulis) && return Vector{Int}[]
    null_basis = _gf2_nullspace_basis(transpose(tab.paulis))
    cycles = Vector{Int}[]
    for relation in null_basis
        cycle = findall(==(UInt8(1)), relation)
        isempty(cycle) || push!(cycles, cycle)
    end
    return cycles
end

"""
    cycle_augmented_graph(tab; cycle_strategy=:cycle_basis)

Return the anticommutation graph without auxiliary relation vertices.
The full relation space is checked algebraically during automorphism search:
an arbitrary basis must not be required to permute into itself, and enumerating
all circuits can be exponential. The empty auxiliary set cannot exclude any
symmetry. The legacy name and `:cycle_basis` keyword remain for compatibility.
"""
function cycle_augmented_graph(tab::SymplecticTableau; cycle_strategy::Symbol=:cycle_basis)
    cycle_strategy == :cycle_basis || throw(ArgumentError("Unsupported SympleQ cycle strategy $cycle_strategy."))
    return anticommutation_graph(tab)
end
