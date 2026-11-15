# Generator-only decomposition of a closed Pauli row representation over C.
# Pair orbits give the exact dimension of the commutant. Two Hermitian samples
# split irreducible copies and align equivalent copies (RepLAB, arXiv:1911.09154).
# We verify the generator intertwining equations and the commutant dimension
# before discarding repeated blocks. These are floating-point checks, not an
# interval-arithmetic certificate.
struct _CliffordGenerators{T<:Integer}
    generators::Vector{CliffordSymmetry{T}}
end

_sw_action_elements(group::_CliffordGenerators) = group.generators

function _signed_row_actions(basis, generators)
    lookup = Dict(w => i for (i, w) in enumerate(basis))
    length(lookup) == length(basis) || throw(ArgumentError("Moment rows must be distinct."))
    actions = Tuple{Vector{Int},Vector{Int}}[]
    for g in generators
        perm, signs = Int[], Int[]
        for w in basis
            sign, image = _act_monomial(g, w)
            haskey(lookup, image) || throw(ArgumentError(
                "Generator-only reduction requires a closed row basis; $w maps to missing row $image."
            ))
            push!(perm, lookup[image])
            push!(signs, sign)
        end
        push!(actions, (perm, signs))
    end
    return actions
end

function _signed_pair_orbits(n, actions)
    labels = zeros(Int, n, n)
    signs = zeros(Int8, n, n)
    count = 0
    queue = CartesianIndex{2}[]
    for seed in CartesianIndices(labels)
        signs[seed] == 0 || continue
        empty!(queue)
        push!(queue, seed)
        signs[seed] = 1
        zero_orbit = false
        pos = 1
        while pos <= length(queue)
            i, j = Tuple(queue[pos])
            pos += 1
            for (p, s) in actions
                target = CartesianIndex(p[i], p[j])
                sign = signs[i,j] * s[i] * s[j]
                if signs[target] == 0
                    signs[target] = sign
                    push!(queue, target)
                elseif signs[target] != sign
                    zero_orbit = true
                end
            end
        end
        if !zero_orbit
            count += 1
            for index in queue
                labels[index] = count
            end
        end
    end
    return labels, signs, count
end

function _commutant_sample(labels, signs, count, seed)
    weights = ComplexF64[
        complex(_pseudorandom_weight(i, seed), _pseudorandom_weight(i, seed + UInt64(1)))
        for i in 1:count
    ]
    C = zeros(ComplexF64, size(labels))
    for index in eachindex(labels)
        labels[index] == 0 && continue
        C[index] = signs[index] * weights[labels[index]]
    end
    return Hermitian((C + C') / 2)
end

function _generator_blocks(n, actions, labels, signs, count, trial; tol=1e-8)
    seed = UInt64(0x636c6966666f7264) + UInt64(4trial)
    C = _commutant_sample(labels, signs, count, seed)
    values, vectors = eigen(C)
    threshold = tol * max(1.0, maximum(abs, values))
    ranges = UnitRange{Int}[]
    start = 1
    for i in 2:n
        if abs(values[i] - values[start]) > threshold
            push!(ranges, start:i-1)
            start = i
        end
    end
    push!(ranges, start:n)
    copies = [vectors[:, r] for r in ranges]
    second = _commutant_sample(labels, signs, count, seed + UInt64(2))
    second_scale = tol * max(1.0, norm(second))
    # Transform once; every copy-to-copy intertwiner is a small submatrix.
    second_in_copies = vectors' * (second * vectors)
    groups = Vector{Matrix{ComplexF64}}[]
    used = falses(length(copies))
    for i in eachindex(copies)
        used[i] && continue
        Q = copies[i]
        aligned = [Q]
        used[i] = true
        for j in i+1:length(copies)
            used[j] && continue
            F = second_in_copies[ranges[i], ranges[j]]
            norm(F) <= second_scale && continue
            size(Q, 2) == size(copies[j], 2) || return nothing
            decomp = svd(F)
            minimum(decomp.S) > tol || return nothing
            maximum(decomp.S) - minimum(decomp.S) <= tol * max(1.0, maximum(decomp.S)) || return nothing
            push!(aligned, copies[j] * (decomp.U * decomp.Vt)')
            used[j] = true
        end
        push!(groups, aligned)
    end

    # If all grouped copies intertwine and sum(multiplicity^2) is the EXACT
    # commutant dimension, no irreducible type was missed or split incorrectly.
    sum(length(g)^2 for g in groups) == count || return nothing
    for (p, s) in actions, group in groups
        Q = first(group)
        acted = similar(Q)
        acted[p, :] = s .* Q
        D = Q' * acted
        for V in group
            acted[p, :] = s .* V
            norm(acted - V * D) <= tol * max(1.0, sqrt(size(V, 2))) || return nothing
        end
    end
    # Select one coordinate of each irreducible type. Schur's lemma makes
    # the other coordinates identical PSD blocks, so they need not be built.
    return [reduce(vcat, [Matrix(Q[:, 1:1]') for Q in group]) for group in groups]
end

function _sw_decompose_half_basis(
    basis::Vector{NormalMonomial{PauliAlgebra,T}}, group::_CliffordGenerators{T},
) where {T<:Integer}
    actions = _signed_row_actions(basis, group.generators)
    labels, signs, count = _signed_pair_orbits(length(basis), actions)
    for trial in 1:3
        blocks = _generator_blocks(length(basis), actions, labels, signs, count, trial)
        isnothing(blocks) || return blocks
    end
    throw(ArgumentError(
        "Generator-only block decomposition failed its numerical intertwining/dimension checks; no reduced problem was returned."
    ))
end
