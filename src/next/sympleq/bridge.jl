# SympleQ's binary search is shared with the historical implementation. Its
# boundary here uses only packed words and checked CliffordMap actions.
function _sympleq_clifford(S::SymplecticMatrix,phase::PhaseVector,::Type{W}) where {W<:PauliWord}
    n = size(S.data,1) ÷ 2
    length(phase) == 2n || throw(DimensionMismatch("Phase vector does not match symplectic matrix."))
    columns = W[]
    negative = Bool[]
    for i in 1:2n
        word,eta = _symplectic_row_to_pauli_word(view(S.data,i,:),W)
        isone(word) && throw(ArgumentError("SympleQ generator image became identity."))
        k = mod(Int(phase[i])-eta,4)
        iseven(k) || throw(ArgumentError("SympleQ phase lift must give Hermitian Pauli images."))
        push!(columns,word)
        push!(negative,k==2)
    end
    return CliffordMap(columns,negative)
end

function _sympleq_phases(tab::SymplecticTableau{W},generator::SympleQGenerator) where {W}
    base = PhaseVector([last(_symplectic_row_to_pauli_word(view(generator.S.data,i,:),W))
        for i in 1:2tab.nqubits],true)
    clifford = _sympleq_clifford(generator.S,base,W)
    rhs = UInt8[]
    for i in 1:length(tab)
        j = generator.permutation.images[i]
        phase,image = act(clifford,tab.monomials[i])
        image == tab.monomials[j] || throw(ArgumentError("SympleQ phase solve mapped to a different term."))
        source,target = tab.coeffs[i],tab.coeffs[j]
        target == source || target == -source || throw(ArgumentError("Incompatible SympleQ coefficients."))
        push!(rhs, UInt8(xor(phase.k==2,target == -source)))
    end
    solution = _gf2_affine_solution_space(tab.paulis,rhs)
    isnothing(solution) && throw(ArgumentError("Inconsistent SympleQ phase equations."))
    particular,kernel = solution
    lift(bits) = PhaseVector(base.values .+ 2 .* bits,true)
    phases = [lift(particular)]
    active = sort!(unique!(Int[s for w in tab.monomials for s in sites(w)]))
    for direction in kernel
        any(s->direction[s]==1 || direction[tab.nqubits+s]==1,active) || continue
        push!(phases,lift(particular .⊻ direction))
    end
    return phases
end

function _sympleq_preserves(g::CliffordMap{W},h::Polynomial) where {W}
    mapped = PolynomialBuffer{W,ComplexF64}(length(h.data))
    for (c,w) in terms(h)
        phase,image = act(g,w)
        push_term!(mapped,scale(phase,c),image)
    end
    return _polynomial(h.registry,mapped) == h
end

# Preserve the historical signed-permutation generator reduction and ordering.
# If the small domain is not closed, keep the full candidate list.
function _sympleq_reduce_generators(candidates::Vector{CliffordMap{W}},n::Int) where {W}
    words = [pauli_letter(W,i,a) for i in 1:n for a in 0:2]
    lookup = Dict(w=>i for (i,w) in enumerate(words))
    for candidate in candidates, i in 1:n, a in 0:2
        _,word = act(candidate,pauli_letter(W,i,a))
        haskey(lookup,word) && continue
        length(words) < 12n || return candidates
        push!(words,word)
        lookup[word] = length(words)
    end
    d = length(words)
    permutations = PermutationGroups.Perm{Int}[]
    generators = CliffordMap{W}[]
    group = PermutationGroups.PermGroup(PermutationGroups.Perm{Int}(collect(1:2d)))
    for candidate in candidates
        images = Vector{Int}(undef,2d)
        for (i,word) in enumerate(words)
            phase,image = act(candidate,word)
            j = get(lookup,image,0)
            iszero(j) && return candidates
            images[i] = phase.k==0 ? j : j+d
            images[i+d] = phase.k==0 ? j+d : j
        end
        permutation = PermutationGroups.Perm{Int}(images)
        permutation in group && continue
        push!(generators,candidate)
        push!(permutations,permutation)
        group = PermutationGroups.PermGroup(permutations)
    end
    return generators
end

function _sympleq_clifford_candidates(h::Polynomial{W};kwargs...) where {W<:PauliWord}
    tab = SymplecticTableau(h)
    tab.nqubits > 0 || throw(ArgumentError("SympleQ needs a non-constant Pauli Hamiltonian."))
    generators = sympleq_generators(tab;kwargs...)
    push!(generators,SympleQGenerator(TermPermutation(collect(1:length(tab))),
        SymplecticMatrix(_gf2_identity(2tab.nqubits)),
        PhaseVector(zeros(Int8,2tab.nqubits),true),true))
    candidates = CliffordMap{W}[]
    seen = Set{CliffordMap{W}}()
    for generator in generators
        phases = try
            _sympleq_phases(tab,generator)
        catch err
            @warn "Skipping SympleQ generator whose phase lift is not representable as a CliffordMap." exception=(err,catch_backtrace())
            continue
        end
        for phase in phases
            clifford = try
                _sympleq_clifford(generator.S,phase,W)
            catch err
                @warn "Skipping SympleQ phase solution that does not define a valid CliffordMap." exception=(err,catch_backtrace())
                continue
            end
            isone(clifford) && continue
            _sympleq_preserves(clifford,h) || continue
            clifford in seen && continue
            push!(seen,clifford)
            push!(candidates,clifford)
        end
    end
    return _sympleq_reduce_generators(candidates,tab.nqubits)
end
