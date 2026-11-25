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

function _sympleq_phases(tab::SymplecticTableau{W},perm,S; kernel_lifts=true) where {W}
    base = PhaseVector([last(_symplectic_row_to_pauli_word(view(S.data,i,:),W))
        for i in 1:2tab.nqubits],true)
    clifford = _sympleq_clifford(S,base,W)
    rhs = UInt8[]
    for i in 1:length(tab)
        j = perm.images[i]
        phase,image = act(clifford,tab.monomials[i])
        image == tab.monomials[j] || throw(ArgumentError("SympleQ phase solve mapped to a different term."))
        source,target = tab.coeffs[i],tab.coeffs[j]
        target == source || target == -source || throw(ArgumentError("Incompatible SympleQ coefficients."))
        push!(rhs, UInt8(xor(phase.k==2,target == -source)))
    end
    !kernel_lifts && all(iszero,rhs) && return [base]
    solution = _gf2_affine_solution_space(tab.paulis,rhs)
    isnothing(solution) && throw(ArgumentError("Inconsistent SympleQ phase equations."))
    particular,kernel = solution
    lift(bits) = PhaseVector(base.values .+ 2 .* bits,true)
    phases = [lift(particular)]
    kernel_lifts || return phases
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

function _sympleq_clifford_candidates(h::Polynomial{W};ga_backend=:backtracking,known_generators=(),kwargs...) where {W<:PauliWord}
    tab = SymplecticTableau(h)
    tab.nqubits > 0 || throw(ArgumentError("SympleQ needs a non-constant Pauli Hamiltonian."))
    ga_backend==:seeded && return _seeded_clifford_candidates(h,tab,known_generators;kwargs...)
    isempty(known_generators) || throw(ArgumentError("Use ga_backend=:seeded with known_generators."))
    generators = sympleq_generators(tab;ga_backend,kwargs...)
    push!(generators,SympleQGenerator(TermPermutation(collect(1:length(tab))),
        SymplecticMatrix(_gf2_identity(2tab.nqubits)),
        PhaseVector(zeros(Int8,2tab.nqubits),true),true))
    candidates = CliffordMap{W}[]
    seen = Set{CliffordMap{W}}()
    for generator in generators
        phases = try
            _sympleq_phases(tab,generator.permutation,generator.S)
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

function _seeded_clifford_candidates(h,tab::SymplecticTableau{W},known;
    cycle_strategy=:cycle_basis,max_automorphisms=100_000) where {W}
    max_automorphisms>0 || throw(ArgumentError("Supply a positive automorphism budget."))
    graph=cycle_augmented_graph(tab;cycle_strategy)
    n=tab.nqubits
    active=unique([i for w in tab.monomials for i in sites(w)])
    verified(g)=_sympleq_preserves(g,h)
    cs=filter(verified,_site_symmetry_candidates(W,n,(Dict(w=>c for (c,w) in terms(h)),),verified))
    for action in known
        g=_group_action(W,action)
        g isa CliffordMap{W} && !is_antiunitary(g) && length(g.columns)≤2h.registry.nsites && verified(g) ||
            throw(ArgumentError("A known generator must be a unitary Clifford symmetry of this objective."))
        push!(cs,g)
    end
    # Sign coordinates of the identity lift are [sign(X_i) | sign(Z_i)].
    # Keep only the active-site kernel, as in the unpruned detector.
    kernel=_gf2_nullspace_basis(tab.paulis)
    for v in kernel
        any(i->!iszero(v[i]|v[n+i]),active) || continue
        push!(cs,CliffordMap(_identity_columns(W,n),Bool.(v)))
    end
    lookup=Dict(w=>i for (i,w) in enumerate(tab.monomials))
    P=PermutationGroups.Perm{Int}
    pruners=Vector{Int}[]
    perms=[P(collect(1:length(tab)))]
    image=PermutationGroups.PermGroup(copy(perms))
    function remember!(p,permutation)
        push!(pruners,p); push!(perms,permutation); image=PermutationGroups.PermGroup(copy(perms))
    end
    for g in cs
        p=[lookup[last(act(g,w))] for w in tab.monomials]
        # Perm owns and trims its input; the search indexes the full vertex domain.
        permutation=P(copy(p))
        permutation in image || remember!(p,permutation)
    end

    function lift!(images)
        perm=TermPermutation(images)
        try
            S=symplectic_matrix_from_permutation(tab,perm)
            phase=only(_sympleq_phases(tab,perm,S;kernel_lifts=false))
            g=_sympleq_clifford(S,phase,W)
            verified(g) || return false
            push!(cs,g)
            return true
        catch err
            err isa ArgumentError || rethrow()
            @warn "Skipping term permutation without a verified Clifford lift." exception=err
            return false
        end
    end
    # A supplied action can also move unused trailing registry sites, which
    # the tableau omits. Its term image then does not specify the full map.
    full_rank=isempty(kernel) && all(g -> length(g.columns)≤2n,cs)
    function witness(p)
        permutation=P(copy(p))
        permutation in image && return
        # Full rank: one binary lift, and only VERIFIED lifts may prune the
        # phase search. Otherwise grow the term group and replay its lifts below.
        (!full_rank || lift!(p)) && remember!(p,permutation)
    end
    leaves=_enumerate_color_preserving_automorphisms(graph;paulis=tab.paulis,
        max_automorphisms,pruners,on_leaf=witness)
    @debug "Seeded SympleQ search" leaves=length(leaves) term_order=GroupsCore.order(BigInt,image) full_rank
    if !full_rank
        # A product of chosen lifts need not be the historical Witt extension
        # of the product. Replay term permutations, NOT graph branches, so we
        # retain differences acting trivially on terms (including the ring kernel).
        for (i,p) in enumerate(image)
            if i>max_automorphisms
                @warn "Seeded Clifford lift replay reached max_automorphisms=$max_automorphisms; returning a verified subgroup."
                break
            end
            isone(p) || lift!([j^p for j in 1:length(tab)])
        end
    end
    _sympleq_reduce_generators(unique(cs),h.registry.nsites)
end
