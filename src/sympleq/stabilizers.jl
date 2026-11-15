# Exact stabilizer-chain search. Each level finds one transporter per reachable
# image of a base point, then fixes that point. Transversals at all levels
# generate the group; we never form their Cartesian product (the group list).
# Search can still be exponential. Exhausting the budget is an error, not a
# partially detected group advertised as complete.
function _sympleq_search_tick!(budget)
    budget[] -= 1
    budget[] >= 0 || throw(ArgumentError(
        "Complete Clifford detection exceeded max_search_nodes; no partial symmetry spec was returned."
    ))
end

function _sympleq_lift(tab::SymplecticTableau{T}, perm, S) where {T}
    candidate = SympleQGenerator(perm, S, PhaseVector(zeros(Int8, 2tab.nqubits), false), false)
    phases = _sympleq_phase_vectors_for_hamiltonian(tab, candidate; include_inactive=true)
    return sympleq_clifford_symmetry(S, first(phases), T, tab.nqubits)
end

function _sympleq_term_stabilizers(tab, budget)
    graph = anticommutation_graph(tab)
    order = _automorphism_order(graph)
    candidates = _automorphism_candidates(graph)
    base = _gf2_independent_row_indices(tab.paulis)
    assigned = zeros(Int, length(tab))
    used = falses(length(tab))
    result = CliffordSymmetry[]

    function search()
        v_index = findfirst(v -> assigned[v] == 0, order)
        if isnothing(v_index)
            perm = TermPermutation(assigned)
            S = symplectic_matrix_from_permutation(tab, perm)
            # A valid linear isometry can still fail the signed coefficient
            # equations. Test INSIDE the transporter search, not by filtering
            # generators of the larger unsigned group afterwards.
            try
                return _sympleq_lift(tab, perm, S)
            catch err
                if err isa ArgumentError && occursin("phase equations are inconsistent", err.msg)
                    return nothing
                end
                rethrow()
            end
        end
        v = order[v_index]
        for w in candidates[v]
            used[w] && continue
            _sympleq_search_tick!(budget)
            (tab.coeffs[w] == tab.coeffs[v] || tab.coeffs[w] == -tab.coeffs[v]) || continue
            _automorphism_consistent(graph.graph, assigned, v, w) || continue
            assigned[v] = w
            used[w] = true
            selected = findall(!iszero, assigned)
            P = tab.paulis[selected, :]
            Q = tab.paulis[assigned[selected], :]
            rank = _gf2_rank(P)
            # Equality of relation spaces, independent of any chosen basis.
            valid = _gf2_rank(Q) == rank && _gf2_rank(hcat(P, Q)) == rank
            found = valid ? search() : nothing
            assigned[v] = 0
            used[w] = false
            isnothing(found) || return found
        end
        return nothing
    end

    fixed = Int[]
    for b in base
        for w in candidates[b]
            (w == b || w in fixed) && continue
            Base.fill!(assigned, 0)
            Base.fill!(used, false)
            for v in fixed
                assigned[v] = v
                used[v] = true
            end
            _automorphism_consistent(graph.graph, assigned, b, w) || continue
            (tab.coeffs[w] == tab.coeffs[b] || tab.coeffs[w] == -tab.coeffs[b]) || continue
            assigned[b] = w
            used[w] = true
            selected = [fixed; b]
            Q = tab.paulis[assigned[selected], :]
            _gf2_rank(Q) == length(selected) || continue
            found = search()
            isnothing(found) || push!(result, found)
        end
        push!(fixed, b)
    end
    return result
end

function _sympleq_kernel_stabilizers(tab::SymplecticTableau{T}, budget) where {T}
    dim = 2tab.nqubits
    D = tab.paulis[_gf2_independent_row_indices(tab.paulis), :]
    perm = TermPermutation(collect(1:length(tab)))
    result = CliffordSymmetry{T}[]
    while size(D, 1) < dim
        d = _choose_domain_extension(D)
        A = _gf2_matmul(D, symplectic_form(tab.nqubits))
        rhs = UInt8[_gf2_pairing(view(D, i, :), d) for i in axes(D, 1)]
        particular, directions = _gf2_affine_solution_space(A, rhs)
        # Enumerate possible images of ONE extension vector, not symplectic
        # matrices or group elements. Witt extension supplies a transporter.
        function visit(r, depth)
            if depth <= length(directions)
                visit(r, depth + 1)
                visit(r .⊻ directions[depth], depth + 1)
                return
            end
            _sympleq_search_tick!(budget)
            (r == d || _gf2_row_in_span(r, D)) && return
            domain = _gf2_append_row(D, d)
            target = _gf2_append_row(D, r)
            U, V = _extend_symplectic_isometry(domain, target)
            S = SymplecticMatrix(_gf2_matmul(_gf2_inverse(U), V))
            push!(result, _sympleq_lift(tab, perm, S))
        end
        visit(particular, 1)
        D = _gf2_append_row(D, d)
    end
    # Kernel of Clifford -> symplectic maps: Pauli sign characters which
    # vanish on all Hamiltonian labels. Include unused sites within 1:n.
    S = SymplecticMatrix(_gf2_identity(dim))
    candidate = SympleQGenerator(perm, S, PhaseVector(zeros(Int8, dim), true), true)
    for phase in _sympleq_phase_vectors_for_hamiltonian(tab, candidate; include_inactive=true)
        g = sympleq_clifford_symmetry(S, phase, T, tab.nqubits)
        _sympleq_clifford_is_identity(g, tab.nqubits) || push!(result, g)
    end
    return result
end

function _sympleq_complete_candidates(H::Polynomial{PauliAlgebra,T}; max_search_nodes) where {T}
    max_search_nodes > 0 || throw(ArgumentError("max_search_nodes must be positive."))
    tab = SymplecticTableau(H)
    tab.nqubits > 0 || throw(ArgumentError("SympleQ needs a non-constant Pauli Hamiltonian."))
    budget = Ref(max_search_nodes)
    candidates = [_sympleq_term_stabilizers(tab, budget); _sympleq_kernel_stabilizers(tab, budget)]
    result = CliffordSymmetry{T}[]
    seen = Set{Any}()
    for g in candidates
        _act_polynomial(g, H) == H || error("Internal Clifford lift failed Hamiltonian invariance.")
        key = _sympleq_clifford_key(g, tab.nqubits)
        key in seen && continue
        push!(seen, key)
        push!(result, g)
    end
    return result
end
