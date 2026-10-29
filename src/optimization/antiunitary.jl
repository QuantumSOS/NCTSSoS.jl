"""
    PauliAntiunitary(x, z)

The antiunitary operator `Θ = P(x,z) K`, where `K` is computational-basis
conjugation and `P` is a Hermitian Pauli word. Binary vectors use site order.
Its square is `(-1)^sum(x .* z) I`; both signs are supported.
"""
struct PauliAntiunitary
    x::BitVector
    z::BitVector
    function PauliAntiunitary(x::AbstractVector, z::AbstractVector)
        length(x) == length(z) || throw(DimensionMismatch("x and z must have equal length"))
        all(a -> a == 0 || a == 1, [x; z]) || throw(ArgumentError("Expected binary vectors"))
        new(BitVector(x), BitVector(z))
    end
end

"""Return whether a Hermitian Pauli word changes sign under `theta`."""
function antiunitary_parity(theta::PauliAntiunitary, mono::NormalMonomial{PauliAlgebra})
    parity = false
    for idx in mono.word
        site, axis = _pauli_site(idx), _pauli_type(idx)
        site <= length(theta.x) || throw(ArgumentError("Pauli word exceeds antiunitary site count"))
        # X=(1,0), Y=(1,1), Z=(0,1). Evaluate q₀ + ω without a dense vector.
        parity ⊻= (axis == 1) ⊻ ((axis != 0) && theta.x[site]) ⊻ ((axis != 2) && theta.z[site])
    end
    return parity
end

function _antiunitary_polynomials(pop::PolyOpt{PauliAlgebra,T,P}) where {T,P<:Polynomial}
    polys = [pop.objective; pop.eq_constraints; pop.ineq_constraints; pop.moment_eq_constraints]
    all(_is_hermitian_polynomial, polys) || throw(ArgumentError(
        "Antiunitary reduction requires every constraint to be Hermitian and individually fixed."
    ))
    return polys
end

"""
    detect_pauli_antiunitary(pop) -> Union{Nothing,PauliAntiunitary}

Solve `ω(r,v) = q₀(v)` over the binary field for every nonzero word in the
objective and all constraints. Here `ω(r,v)=rₓ⋅v_z+r_z⋅vₓ` and `q₀(v)=vₓ⋅v_z`.
Return `nothing` if inconsistent. Free variables are set to zero deterministically;
this finds one symmetry, not all symmetries or a prescribed square sign.
"""
function detect_pauli_antiunitary(pop::PolyOpt{PauliAlgebra,T,P}) where {T,P<:Polynomial}
    polys = _antiunitary_polynomials(pop)
    n = maximum(_pauli_nqubits, polys)
    rows = BitVector[]
    for poly in polys, (c, mono) in terms(poly)
        iszero(c) && continue
        v, _ = _pauli_word_to_symplectic(mono, n)
        push!(rows, BitVector([v[n+1:end]; v[1:n]; isodd(sum(v[i] * v[n+i] for i in 1:n; init=0))]))
    end
    pivot = 1
    columns = Int[]
    for col in 1:2n
        found = findnext(row -> row[col], rows, pivot)
        isnothing(found) && continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        for k in eachindex(rows)
            k != pivot && rows[k][col] && (rows[k] .⊻= rows[pivot])
        end
        push!(columns, col)
        pivot += 1
    end
    any(row -> !any(row[1:2n]) && row[end], rows) && return nothing
    r = falses(2n)
    for (i, col) in enumerate(columns)
        r[col] = rows[i][end]
    end
    return PauliAntiunitary(r[1:n], r[n+1:end])
end

function _antiunitary_project(poly::Polynomial{PauliAlgebra,T,C}, theta, phase=1) where {T,C}
    R = typeof(real(zero(C)))
    kept = Tuple{R,NormalMonomial{PauliAlgebra,T}}[]
    for (c, mono) in terms(poly)
        antiunitary_parity(theta, mono) && continue
        value = phase * c
        iszero(imag(value)) || throw(ArgumentError("Transformed Pauli block is not real"))
        push!(kept, (real(value), mono))
    end
    return Polynomial{PauliAlgebra,T,R}(kept)
end

"""
    antiunitary_moment_relaxation(pop, sparsity; theta=detect_pauli_antiunitary(pop))

Build an invariant Pauli relaxation with odd expectations eliminated and each
Hermitian block changed directly to a real symmetric block of the same size.
Returns `(moment_problem, original, theta, phases)`. `phases[b]` is the diagonal
of `U` for PSD block `b`; rows and sparse block provenance are preserved.
No parity splitting or complex-to-real doubling is performed. Supply `theta`
explicitly to select another solution, including one with square `-I`.
"""
function antiunitary_moment_relaxation(pop::PolyOpt{PauliAlgebra,T,P}, sparsity::SparsityResult;
    theta=detect_pauli_antiunitary(pop),
) where {T,P<:Polynomial}
    isnothing(theta) && throw(ArgumentError("No individually invariant Pauli antiunitary exists"))
    polys = _antiunitary_polynomials(pop)
    for poly in polys, (c, mono) in terms(poly)
        !iszero(c) && antiunitary_parity(theta, mono) && throw(ArgumentError(
            "The supplied antiunitary does not fix the objective and every constraint"
        ))
    end
    # Include unused sites present only in the relaxation basis with trivial action.
    n = maximum((_pauli_site(i) for i in indices(pop.registry)); init=0)
    theta = PauliAntiunitary(vcat(theta.x, falses(max(0, n-length(theta.x)))),
                            vcat(theta.z, falses(max(0, n-length(theta.z)))))
    original = moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    objective = _antiunitary_project(original.objective, theta)
    RP = typeof(objective)
    M = eltype(original.total_basis)
    constraints = Tuple{Symbol,Matrix{RP}}[]
    metadata = Dict{Int,BlockMeta{M}}()
    phases = Vector{Complex{Int}}[]
    block_index = Dict(c => b for (b,c) in enumerate(original.linear.psd_block_constraint_idx))
    for (idx, (cone, mat)) in enumerate(original.constraints)
        if cone == :HPSD
            block = original.linear.psd_blocks_lin[block_index[idx]]
            u = Complex{Int}[antiunitary_parity(theta, m) ? im : 1 for m in block.meta.row_labels]
            transformed = [_antiunitary_project(mat[i,j], theta, u[i]*conj(u[j])) for i in axes(mat,1), j in axes(mat,2)]
            push!(constraints, (:PSD, transformed))
            metadata[length(constraints)] = BlockMeta(:PSD, block.meta.origin, block.meta.row_labels)
            push!(phases, u)
        elseif cone == :Zero
            # Equality matrices need not carry row metadata. Their real and
            # imaginary scalar equations are equivalent after odd elimination.
            for poly in mat
                even = RP([(real(c), m) for (c,m) in terms(poly) if !antiunitary_parity(theta,m)])
                imagpart = RP([(imag(c), m) for (c,m) in terms(poly) if !antiunitary_parity(theta,m)])
                iszero(even) || push!(constraints, (:Zero, reshape([even], 1, 1)))
                iszero(imagpart) || push!(constraints, (:Zero, reshape([imagpart], 1, 1)))
            end
        else
            throw(ArgumentError("Unsupported antiunitary input cone $cone"))
        end
    end
    basis = filter(m -> !antiunitary_parity(theta, m), original.total_basis)
    moment_keys = Set(key for block in original.linear.psd_blocks_lin
        if block.meta.origin isa MomentMatrixOrigin for form in block.entries for (key, _) in form
        if !antiunitary_parity(theta, original.linear.key_to_monomial[key]))
    mp = MomentProblem{PauliAlgebra,T,M,RP}(objective, constraints, basis, length(moment_keys);
        block_meta_by_constraint=metadata, real_moments=true)
    return (moment_problem=mp, original=original, theta=theta, phases=phases)
end

"""
    recover_antiunitary_gram(reduction, block, Z)

Recover the original complex positive semidefinite multiplier `Q = U† Z U`
from a real symmetric multiplier `Z`. The certificate convention is `tr(Q M)`.
Because the localizer is individually invariant, odd coefficients cancel in
this full polynomial expression, not merely after setting moments to zero.
Equality multipliers remain scalar real/imaginary equations in the returned
moment problem and must also be included when checking a constrained certificate.
"""
function recover_antiunitary_gram(reduction, block::Integer, Z::AbstractMatrix)
    u = reduction.phases[block]
    size(Z) == (length(u), length(u)) || throw(DimensionMismatch("Gram matrix size does not match block"))
    return [conj(u[i]) * Z[i,j] * u[j] for i in eachindex(u), j in eachindex(u)]
end
