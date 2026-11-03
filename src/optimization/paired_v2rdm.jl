"""
    build_paired_moment_data(J; pairs, conditions=:PQG) -> MomentLinearData

Build the real seniority-zero paired-fermion relaxation of
`sum(J[i,j] * K⁺[i] * K⁻[j])`, where each level is empty or doubly occupied,
`K⁺[i] = c†[2i-1]c†[2i]`, and `pairs` fixes the number of occupied levels.
`J` must be real symmetric. This is not unrestricted fermionic PQG.

Implements Eqs. (15)–(35) of Rubio-García et al., arXiv:1805.01322v1.
`:PQGT` adds both T1 and T2; number and contraction constraints are exact
equalities, unlike the paper's small inequality tolerance. Real moments suffice
because conjugation preserves this feasible set and the real objective.
"""
function build_paired_moment_data(J::AbstractMatrix{<:Real}; pairs::Integer, conditions::Symbol=:PQG)
    L = size(J, 1)
    L >= 1 && size(J, 2) == L || throw(ArgumentError("J must be nonempty and square"))
    issymmetric(J) && all(isfinite, J) || throw(ArgumentError("J must be finite and symmetric"))
    0 <= pairs <= L || throw(ArgumentError("pairs must lie between zero and the number of levels"))
    conditions in (:PQG, :PQGT) || throw(ArgumentError("conditions must be :PQG or :PQGT"))

    T = _v2rdm_index_type(2L)
    K = Vector{T}
    M = NormalMonomial{FermionicAlgebra,T}
    F = LinearMomentForm{K,Float64}
    monomials = Dict{K,M}()
    # Canonical real representatives retain actual fermionic words, rather than
    # treating pair operators as ordinary bosons.
    function raw_form(word)
        terms = Pair{K,Float64}[]
        for (key, coefficient) in normal_order_terms(T[word...])
            adj = _v2rdm_adjoint_word(key)
            representative = key_lt(adj, key) ? adj : key
            monomials[representative] = M(representative)
            push!(terms, representative => Float64(coefficient))
        end
        return F(terms)
    end
    combine(terms...) = F(key => coefficient * value
        for (coefficient, form) in terms for (key, value) in form)
    one_form = raw_form(())
    P = [raw_form((-(2i-1), -2i, 2j, 2j-1)) for i in 1:L, j in 1:L]
    D = [i == j ? P[i,i] : raw_form((-(2i-1), 2i-1, -(2j-1), 2j-1))
         for i in 1:L, j in 1:L]
    blocks = PSDBlockLin{K,Float64,M}[]
    function block!(family, indices, entries)
        n = size(entries, 1)
        n == 0 && return
        # These include sums of Gram matrices (T1/T2), not a single monomial
        # basis. Identity placeholders carry no basis semantics; origin and
        # the equation's row ordering identify each structured block.
        meta = BlockMeta{M}(:PSD, V2RDMBlockOrigin(family, indices), fill(one(M), n))
        push!(blocks, PSDBlockLin{K,Float64,M}(n, entries, meta))
    end
    scalar!(family, indices, form) = block!(family, indices, reshape([form], 1, 1))
    block!(:P, (), P)
    block!(:Q, (), [i == j ? combine(1=>one_form, -1=>P[i,i]) : P[i,j]
                    for i in 1:L, j in 1:L])
    block!(:g, (), D)
    # The i=j scalar conditions duplicate P/Q diagonal positivity.
    for i in 1:L, j in 1:i-1
        scalar!(:D, (i,j), D[i,j])
        scalar!(:q, (i,j), combine(1=>D[i,j], 1=>one_form, -1=>P[i,i], -1=>P[j,j]))
        block!(:G, (i,j), [combine(1=>P[i,i], -1=>D[i,j]) combine(-1=>P[i,j]);
                          combine(-1=>P[j,i]) combine(1=>P[j,j], -1=>D[i,j])])
    end
    if conditions == :PQGT
        for i in 1:L
            others = [j for j in 1:L if j != i]
            block!(:T1, (i,), [j == k ? combine(1=>one_form, -1=>P[j,j], -1=>P[i,i], 2=>D[i,j]) : P[j,k]
                              for j in others, k in others])
            n = L - 1
            entries = fill(F(Pair{K,Float64}[]), 2n+1, 2n+1)
            for (a,j) in enumerate(others), (b,k) in enumerate(others)
                entries[a,b] = D[j,k]
                entries[n+a,n+b] = j == k ? combine(1=>P[i,i], -2=>D[i,k], 1=>P[j,k]) : P[j,k]
                entries[a,n+b] = j == k ? combine(-1=>P[i,k]) : F(Pair{K,Float64}[])
                entries[n+b,a] = entries[a,n+b]
            end
            for (a,j) in enumerate(others)
                entries[a,end] = entries[end,a] = D[i,j]
                entries[n+a,end] = entries[end,n+a] = P[i,j]
            end
            entries[end,end] = P[i,i]
            block!(:T2, (i,), entries)
        end
        for i in 1:L, j in 1:i-1, k in 1:j-1
            scalar!(:t1, (i,j,k), combine(1=>one_form, -1=>P[i,i], -1=>P[j,j], -1=>P[k,k],
                                        1=>D[i,j], 1=>D[j,k], 1=>D[k,i]))
            indices = (i,j,k)
            entries = [P[a,b] for a in indices, b in indices]
            for (a,(u,v,w)) in enumerate(((i,j,k),(j,i,k),(k,i,j)))
                entries[a,a] = combine(1=>P[u,u], -1=>D[u,v], -1=>D[u,w], 1=>D[v,w])
            end
            block!(:T2triple, indices, entries)
        end
    end
    equalities = ScalarLinearConstraint{K,Float64}[]
    number = combine((-pairs=>one_form), (1=>P[i,i] for i in 1:L)...)
    push!(equalities, ScalarLinearConstraint(number, :zero, MomentEqOrigin(1,1)))
    for i in 1:L
        contraction = combine((-pairs=>P[i,i]), (1=>D[i,j] for j in 1:L)...)
        push!(equalities, ScalarLinearConstraint(contraction, :zero, MomentEqOrigin(2,i)))
    end
    objective = combine((J[i,j]=>P[i,j] for i in 1:L for j in 1:L)...)
    moments = sort!(collect(keys(monomials)); lt=key_lt)
    adjoints = Dict(key=>key for key in moments)
    pivots = _discover_linear_pivots(blocks, adjoints)
    return MomentLinearData{K,Float64,M}(moments,
        Dict(key=>i for (i,key) in enumerate(moments)), T[], monomials, adjoints,
        blocks, collect(eachindex(blocks)), equalities, objective, pivots,
        _build_pivot_at(pivots), K[key for key in moments if !haskey(pivots,key)])
end
