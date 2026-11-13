# Full normalized subsystem density matrix in the computational basis.
# Shared canonical Pauli moments supply overlap consistency automatically.
function _pauli_subsystem_density(::Type{T}, sites) where {T}
    C = _ActionCoefficient
    M = NormalMonomial{PauliAlgebra,T}
    k = length(sites)
    dim = 2^k
    entries = [Tuple{C,M}[] for _ in 1:dim, _ in 1:dim]
    for code0 in 0:(4^k-1)
        code, flip, phase_mask, ny = code0, 0, 0, 0
        word = T[]
        for (position, site) in enumerate(sites)
            axis = code % 4
            code ÷= 4
            axis == 0 && continue
            push!(word, T(_pauli_index(site, axis-1)))
            bit = 1 << (k-position)
            axis in (1,2) && (flip |= bit)
            axis in (2,3) && (phase_mask |= bit)
            ny += axis == 2
        end
        monomial = M(word)
        for col in 0:(dim-1)
            coefficient = C(im)^ny * (isodd(count_ones(col & phase_mask)) ? -1 : 1) / dim
            push!(entries[xor(col,flip)+1,col+1], (coefficient,monomial))
        end
    end
    return Polynomial.(entries)
end

# Exact physical constraints, before choosing the shared moment domain.
function _pauli_strengthening_constraints(pop::PolyOpt{PauliAlgebra,T}, original;
    subsystems, perturbations, stationarity, max_block_size,
    stationarity_operators=perturbations, magnetization::Bool=false) where {T}
    max_block_size >= 1 || throw(ArgumentError("max_block_size must be positive."))
    M = NormalMonomial{PauliAlgebra,T}
    P = typeof(original.objective)
    patches = Vector{Int}[]
    for patch in subsystems
        sites = sort!(Int[patch...])
        !isempty(sites) && allunique(sites) || throw(ArgumentError("Subsystems must contain distinct sites."))
        side = magnetization ? binomial(big(length(sites)),fld(length(sites),2)) : big(2)^length(sites)
        side <= max_block_size || throw(ArgumentError("Subsystem matrix exceeds max_block_size."))
        for site in sites, axis in 0:2
            _pauli_index(site,axis) in indices(pop.registry) || throw(ArgumentError("Subsystem leaves the Pauli registry."))
        end
        sites in patches || push!(patches,sites)
    end
    all(q isa M for q in perturbations) || throw(ArgumentError("Perturbations must be Pauli monomials from this registry."))
    qs = sort!(unique!(M[perturbations...]))
    length(qs) <= max_block_size || throw(ArgumentError("Perturbative matrix exceeds max_block_size."))
    all(i in indices(pop.registry) for q in qs for i in q.word) || throw(ArgumentError("Perturbation leaves the registry."))
    all(q isa M for q in stationarity_operators) || throw(ArgumentError("Stationarity operators must be Pauli monomials from this registry."))
    ss = stationarity ? sort!(unique!(M[stationarity_operators...])) : M[]
    all(i in indices(pop.registry) for q in ss for i in q.word) || throw(ArgumentError("Stationarity operator leaves the registry."))
    constraints = Tuple{Symbol,Matrix{P}}[]
    labels = Any[]
    for sites in patches
        if magnetization
            for (down_spins,mat) in _pauli_rdm_blocks(length(sites),nothing,P;sign_symmetry=false,sites)
                push!(constraints, (:HPSD,mat / big(2)^length(sites)))
                push!(labels, (kind=:subsystem,sites=sites,down_spins=down_spins))
            end
        else
            push!(constraints, (:HPSD,_pauli_subsystem_density(T,sites)))
            push!(labels, (kind=:subsystem,sites=sites))
        end
    end
    h = original.objective
    if !isempty(qs)
        q = _action_polynomial.(qs)
        mat = Matrix{P}(undef,length(q),length(q))
        for i in eachindex(q), j in eachindex(q)
            right = _action_adjoint(q[j])
            product = _action_product(q[i],right)
            mat[i,j] = _action_product(_action_product(q[i],h),right) -
                (_action_product(h,product)+_action_product(product,h))/2
        end
        push!(constraints,(:HPSD,mat))
        push!(labels,(kind=:perturbative,operators=qs))
    end
    equalities = P[]
    if stationarity
        for q in _action_polynomial.(ss)
            # Q and H are Hermitian: -i[H,Q] has real Pauli coefficients.
            p = -im*(_action_product(h,q)-_action_product(q,h))
            iszero(p) || push!(equalities,p)
        end
    end
    return (;constraints,labels,equalities,patches,qs,stationarity_operators=ss)
end

# Added variables do not add rows to the main positive matrix. Close their
# shared domain before forming an invariant-coordinate map, never zero-fill.
function _pauli_strengthening_domain(words, extra, policy; actions=[], clifford_generators=[])
    policy in (:unchanged,:extend) || throw(ArgumentError("moment_domain must be :unchanged or :extend."))
    domain = Set(words)
    for p in Iterators.flatten(((p for (_,mat) in extra.constraints for p in mat),extra.equalities))
        for (_,word) in terms(p)
            if policy == :unchanged
                word in domain || throw(ArgumentError(
                    "Strengthening requires missing moment $word; use moment_domain=:extend or enlarge the main basis."))
            else
                push!(domain,word)
            end
        end
    end
    if policy == :extend
        pending = collect(domain)
        head = 1
        while head <= length(pending)
            word = pending[head]
            head += 1
            images = [differentiate(a,word) for a in actions]
            append!(images,[_act_polynomial(g,_action_polynomial(word)) for g in clifford_generators])
            for image in images, (_,m) in terms(image)
                m in domain && continue
                push!(domain,m)
                push!(pending,m)
            end
        end
    end
    return sort!(collect(domain))
end

# Public assembly for unconstrained finite Pauli Hamiltonians. This uses the
# existing basis/sparsity selection and, optionally, the checked continuous /
# Clifford composition. Extra cones retain their full matrices; only the main
# moment cones use the continuous representation compression. No magnetization
# sector, translated patch, or sign partner is dropped on an assumed symmetry.
function pauli_strengthened_relaxation(pop::PolyOpt{PauliAlgebra,T}, config::SolverConfig;
    subsystems=Vector{Int}[], perturbations=NormalMonomial{PauliAlgebra,T}[],
    stationarity::Bool=true, max_block_size::Integer=256,
    stationarity_operators=perturbations, moment_domain::Symbol=:unchanged,
    actions::AbstractVector{<:SuppliedDerivation}=SuppliedDerivation[],
    clifford_generators::AbstractVector{<:CliffordSymmetry}=CliffordSymmetry[],
    compact_metric=Matrix{Int}(I,length(pop.registry),length(pop.registry)),
    compress::Bool=true) where {T}
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) && isempty(pop.moment_eq_constraints) ||
        throw(ArgumentError("Pauli strengthening currently requires an unconstrained Hamiltonian."))
    isnothing(config.symmetry) || throw(ArgumentError("Supply actions/clifford_generators for checked strengthening composition, not config.symmetry."))
    h = _action_polynomial(pop.objective)
    h == _action_adjoint(h) || throw(ArgumentError("The Hamiltonian must be Hermitian."))
    sparsity = compute_sparsity(pop,config)
    original = _continuous_exact_problem(pop,sparsity)
    extra = _pauli_strengthening_constraints(pop,original;
        subsystems,perturbations,stationarity,max_block_size,stationarity_operators)
    words = _pauli_strengthening_domain(original.words,extra,moment_domain;actions,clifford_generators)
    original = merge(original,(words=words,))
    # Check preservation of the *added* constraint family, not only the main
    # moment basis. Local algebra preservation makes density positivity covariant.
    for action in actions
        for sites in extra.patches, site in sites, axis in 0:2
            letter = NormalMonomial{PauliAlgebra,T}(T[T(_pauli_index(site,axis))])
            all(_pauli_site(i) in sites for (_,m) in terms(differentiate(action,letter)) for i in m.word) ||
                throw(ArgumentError("Continuous action does not preserve a subsystem algebra."))
        end
        isempty(extra.qs) || _action_coordinates(_action_polynomial.(extra.qs),
            differentiate.(Ref(action),extra.qs); context="Perturbation row span")
        _action_coordinates(extra.equalities,differentiate.(Ref(action),extra.equalities);
            context="Stationarity equality span")
    end
    for g in clifford_generators
        for sites in extra.patches
            image = sort!(unique!(Int[_pauli_site(i) for site in sites for axis in 0:2
                for i in last(_act_monomial(g,NormalMonomial{PauliAlgebra,T}(T[T(_pauli_index(site,axis))]))).word]))
            image in extra.patches || throw(ArgumentError("Clifford action does not preserve the subsystem family; include the entire orbit."))
        end
        all(last(_act_monomial(g,q)) in extra.qs for q in extra.qs) ||
            throw(ArgumentError("Clifford action does not preserve the perturbation basis."))
        _action_coordinates(extra.equalities,[_act_polynomial(g,p) for p in extra.equalities];
            context="Stationarity equality span")
    end
    reduction = nothing
    if !isempty(actions) || !isempty(clifford_generators)
        supplied = isempty(actions) ? [SuppliedDerivation(pop.registry,
            Dict(i=>zero(pop.objective) for i in indices(pop.registry)))] : actions
        reduction = continuous_moment_relaxation(pop,config,supplied;
            compact_metric,compress,clifford_generators,additional_moments=words)
    end
    if isnothing(reduction)
        FP = Polynomial{PauliAlgebra,T,ComplexF64}
        constraints = [(cone,map(p->convert(FP,p),mat)) for (cone,mat) in original.constraints]
        objective, words = convert(FP,original.objective), original.words
    else
        FP = typeof(reduction.moment_problem.objective)
        constraints = copy(reduction.moment_problem.constraints)
        objective, words = reduction.moment_problem.objective, reduction.coordinate_labels
    end
    extra_embedded = Bool[]
    for (cone,mat) in extra.constraints
        if isnothing(reduction)
            push!(constraints,(cone,map(p->convert(FP,p),mat)))
            push!(extra_embedded,false)
        else
            projected = _continuous_symbolic(reduction.moment_map,reduction.section,reduction.coordinate_labels,mat)
            realmat, embedded = _continuous_real_cone(projected)
            push!(constraints,(:PSD,map(p->convert(FP,p),realmat)))
            push!(extra_embedded,embedded)
        end
    end
    for p in extra.equalities
        q = isnothing(reduction) ? p : _continuous_real_polynomial(
            _continuous_symbolic(reduction.moment_map,reduction.section,reduction.coordinate_labels,p))
        push!(constraints,(:Zero,reshape([convert(FP,q)],1,1)))
    end
    mp = MomentProblem(objective,constraints,copy(words),length(words);real_moments=!isnothing(reduction))
    return (moment_problem=mp,original=original,extra=extra,reduction=reduction,
        extra_embedded=extra_embedded,ground_state_only=!isempty(extra.qs) || !isempty(extra.equalities),
        moment_domain=moment_domain,extra_cones_compressed=false)
end

# Verify supplied multipliers against the exact original constraint matrices.
# `factors=true` interprets each input as B and forms B*B' in exact rationals,
# avoiding an unjustified numerical PSD claim. A coefficient l1 residual bound
# corrects the energy because every canonical Pauli word has operator norm one.
# With perturbative/stationarity constraints this certifies ground energy only,
# not positivity of H-bound on every state via these particular multipliers.
function verify_pauli_strengthening_certificate(result,bound,grams::AbstractVector,
    multipliers::AbstractVector=Rational{BigInt}[]; factors::Bool=false)
    original, extra = result.original, result.extra
    length(multipliers) == length(extra.equalities) || throw(DimensionMismatch("Wrong stationarity multiplier count."))
    main_count = isnothing(result.reduction) ? length(original.positive_indices) : length(result.reduction.blocks)
    length(grams) == main_count+length(extra.constraints) || throw(DimensionMismatch("Wrong positive multiplier count."))
    values = [factors ? (b=map(_action_coefficient,g); b*b') : map(_action_coefficient,g) for g in grams]
    lifted = isnothing(result.reduction) ? Dict(idx=>q for (idx,q) in zip(original.positive_indices,values[1:main_count])) :
        recover_continuous_grams(result.reduction,values[1:main_count])
    lambda = _action_coefficient(bound)
    isreal(lambda) || throw(ArgumentError("The bound must be real."))
    residual = original.objective-lambda*one(original.objective)
    positive = all(_continuous_is_psd,values)
    for (idx,q) in lifted
        size(q) == size(original.constraints[idx][2]) || throw(DimensionMismatch("Wrong main multiplier size."))
        residual -= _continuous_trace(original.constraints[idx][2],q)
    end
    for (i,(_,mat)) in enumerate(extra.constraints)
        q = values[main_count+i]
        if result.extra_embedded[i]
            d = size(mat,1)
            size(q) == (2d,2d) || throw(DimensionMismatch("Wrong realified multiplier size."))
            q = q[1:d,1:d]+q[d+1:end,d+1:end]+im*(q[d+1:end,1:d]-q[1:d,d+1:end])
        end
        size(q) == size(mat) || throw(DimensionMismatch("Wrong subsystem/perturbative multiplier size."))
        residual -= _continuous_trace(mat,q)
    end
    for (c,p) in zip(multipliers,extra.equalities)
        isreal(c) || throw(ArgumentError("Stationarity multipliers must be real."))
        residual -= _action_coefficient(c)*p
    end
    if !isnothing(result.reduction)
        invariant_map = result.reduction.moment_map
        coefficients = transpose(invariant_map.reconstruction)*_continuous_coefficients(invariant_map.index,residual)
        residual = Polynomial([(_ActionCoefficient(c),invariant_map.words[invariant_map.free_indices[i]]) for (i,c) in enumerate(coefficients)])
    end
    correction = sum(abs(real(c))+abs(imag(c)) for (c,_) in terms(residual);init=big(0)//1)
    return (certified=positive, certified_bound=positive ? real(lambda)-correction : nothing,
        residual_bound=correction, identity_exact=iszero(residual), ground_state_only=result.ground_state_only,
        symmetry_averaged=!isnothing(result.reduction))
end

# The rotation Casimir C=sum_a D_a^2 commutes with lattice permutations and
# orthogonal spin transformations (including conjugation). On a compact spin
# representation, ker(C) is exactly the common kernel of the three D_a.
# Thus its equations on the discrete-invariant moment coordinates impose full
# spin invariance without assembling derivations on the much larger word space.
function _pauli_chain_spin_section(pop, ops, words, reducer; extend::Bool=false)
    vectors = [[axis[i] for axis in ops] for i in eachindex(first(ops))]
    actions = [SuppliedDerivation(pop.registry,Dict(v[i].word[1] =>
        (i == a ? Polynomial(v[b]) : i == b ? -Polynomial(v[a]) : zero(Polynomial(v[i])))
        for v in vectors for i in 1:3)) for a in 1:3 for b in a+1:3]
    all(iszero(differentiate(a,pop.objective)) for a in actions) ||
        throw(ArgumentError("Full spin reduction requires a rotation-invariant Hamiltonian."))
    function discrete(p)
        even = Polynomial([(c,m) for (c,m) in terms(p) if iszero(_pauli_sign_signature(m))])
        return _reduce_moment_polynomial(even,reducer)
    end
    if extend
        words = collect(words)
        seen = Set(words)
        head = 1
        while head <= length(words)
            p = discrete(sum(differentiate(a,differentiate(a,words[head])) for a in actions))
            head += 1
            for (_,m) in terms(p)
                m in seen && continue
                push!(seen,m)
                push!(words,m)
            end
        end
        sort!(words)
    end
    index = Dict(m=>i for (i,m) in enumerate(words))
    equations = [discrete(sum(differentiate(a,differentiate(a,m)) for a in actions)) for m in words]
    kernel,free = _continuous_nullspace(sparse(transpose(_continuous_coefficients(index,equations))))
    return (;actions,discrete,index,kernel,free,words=words[free],equations)
end

# Numerical cone splitting after exact invariant-moment substitution. Every
# coefficient matrix is checked, not merely a random evaluation or a solution.
# The orthogonal transforms are retained for reconstruction. This is a floating
# point reduction, like the chain Fourier transform, not an exact certificate.
function _pauli_coefficient_blocks(mat; atol=1e-9)
    n = size(mat,1)
    n <= 64 || throw(ArgumentError("Coefficient block splitting is limited to matrices of size 64."))
    P = eltype(mat)
    index = Dict(m=>i for (i,m) in enumerate(sort!(unique([m for p in mat for (_,m) in terms(p)]))))
    coefficients = _continuous_coefficients(index,vec(mat))
    slices = [reshape(Float64.(real.(Array(coefficients[i,:]))),n,n) for i in 1:length(index)]
    all(isreal,coefficients) || throw(ArgumentError("Realify a complex matrix before coefficient splitting."))
    all(norm(a-a',Inf) <= atol for a in slices) || error("Non-Hermitian coefficient matrix")
    eye = Matrix{Float64}(I,n,n)
    if all(norm(a-a[1,1]*eye,Inf) <= atol*max(1.0,norm(a,Inf)) for a in slices)
        splitter = Matrix(Diagonal(collect(1.0:n)))
    else
        # A common symmetric commutant must preserve the eigenspaces of ANY
        # linear combination of the coefficient matrices. Start there, then
        # intersect the much smaller candidate space with each commutant. This
        # avoids forming and diagonalizing a dense n²-by-n² Laplacian.
        seed = zeros(n,n)
        for (j,a) in enumerate(slices)
            scale = norm(a)
            iszero(scale) || (seed .+= (sin(j*sqrt(2.0))/scale).*a)
        end
        spectrum = eigen(Symmetric(seed))
        eigenspaces = Vector{Int}[]
        for i in eachindex(spectrum.values)
            if isempty(eigenspaces) || abs(spectrum.values[i]-spectrum.values[first(last(eigenspaces))]) > 1e-8*max(1.0,norm(seed))
                push!(eigenspaces,[i])
            else
                push!(last(eigenspaces),i)
            end
        end
        candidates = Matrix{Float64}[]
        for group in eigenspaces, (j,q) in enumerate(group), p in group[1:j]
            u,v = spectrum.vectors[:,p],spectrum.vectors[:,q]
            push!(candidates,p == q ? u*u' : (u*v'+v*u')/sqrt(2.0))
        end
        for a in slices
            scale = norm(a)
            iszero(scale) && continue
            a = a/scale
            commutators = [a*x-x*a for x in candidates]
            all(norm(c,Inf) <= atol for c in commutators) && continue
            decomposition = svd(reduce(hcat,vec.(commutators));full=false)
            null = findall(s -> s <= 1e-10*max(1.0,maximum(decomposition.S)),decomposition.S)
            isempty(null) && error("Common commutant lost the identity")
            basis = reduce(hcat,vec.(candidates))*decomposition.V[:,null]
            candidates = [reshape(basis[:,j],n,n) for j in axes(basis,2)]
        end
        splitter = zeros(n,n)
        for (j,x) in enumerate(candidates)
            x = (x+x')/2
            x -= LinearAlgebra.tr(x)/n*eye
            all(norm(a*x-x*a,Inf) <= atol*max(1.0,norm(a,Inf)) for a in slices) || continue
            splitter += sin(j*sqrt(2.0))*x
        end
    end
    split = eigen(Symmetric(splitter))
    groups = Vector{Int}[]
    for i in eachindex(split.values)
        if isempty(groups) || abs(split.values[i]-split.values[first(last(groups))]) > 1e-7*max(1.0,norm(splitter))
            push!(groups,[i])
        else
            push!(last(groups),i)
        end
    end
    # An unsplit cone needs no change of coordinates. Rotating it anyway
    # densifies its polynomial entries without reducing the solver block.
    length(groups) == 1 && return (;blocks=[mat],rows=[eye],error_bound=0.0)
    rows = [split.vectors[:,g] for g in groups]
    norm(split.vectors*split.vectors'-eye,Inf) <= atol || error("Incomplete orthogonal cone transform")
    C = coeff_type(P)
    term_type = Tuple{C,keytype(index)}
    entries = [[term_type[] for _ in 1:size(u,2), _ in 1:size(u,2)] for u in rows]
    error_bound = 0.0
    for (word,i) in index
        a = slices[i]
        reconstructed = zeros(n,n)
        for (u,b) in zip(rows,entries)
            part = u'*a*u
            part = (part+part')/2
            reconstructed += u*part*u'
            for q in axes(part,2), p in axes(part,1)
                iszero(part[p,q]) || push!(b[p,q],(C(part[p,q]),word))
            end
        end
        residual = norm(a-reconstructed,Inf)/max(1.0,norm(a,Inf))
        error_bound = max(error_bound,residual)
    end
    error_bound <= atol || error("Coefficientwise cone reconstruction failed: $error_bound")
    # Canonicalize each entry once, rather than re-sorting its growing term list
    # after every coefficient contribution.
    blocks = [P.(b) for b in entries]
    return (;blocks,rows,error_bound)
end

# Keep the original exact matrix inspectable without storing one polynomial per
# entry. Canonical Pauli words are Hermitian, so the entry is simply rows[i]*rows[j].
struct _PauliMomentMatrix{T<:Integer} <: AbstractMatrix{Polynomial{PauliAlgebra,T,_ActionCoefficient}}
    rows::Vector{NormalMonomial{PauliAlgebra,T}}
end
Base.size(mat::_PauliMomentMatrix) = (length(mat.rows),length(mat.rows))
function Base.getindex(mat::_PauliMomentMatrix, i::Int, j::Int)
    @boundscheck checkbounds(mat,i,j)
    c,m = _mul_term(mat.rows[i],mat.rows[j])
    return Polynomial((_action_coefficient(c),m))
end

# Full physical symmetry pipeline for rotation-invariant periodic Pauli chains.
"""
    pauli_chain_stationarity_operators(pop, ops; basis, subsystems=Vector{Int}[],
        five_site_window=nsites ÷ 2, seven_site_window=7)

Return QMBCertify's historical periodic-chain stationarity candidates whose
commutator equations are supported by the reduced main-basis products and the
supplied subsystem density moments. Zero equations are discarded; equations
that agree up to sign are kept once. Unlike the historical random fingerprint,
duplicate detection is exact. The result also includes `moments`, the exact
sorted supported domain before equation cancellation, and `equations`, the
reduced commutators in the same order as `operators`. Reuse these rather than
rebuilding commutators or inferring support from nonzero matrix entries.
"""
function pauli_chain_stationarity_operators(pop::PolyOpt{PauliAlgebra,T}, ops;
    basis, subsystems=Vector{Int}[], five_site_window::Integer=0,
    seven_site_window::Integer=7) where {T}
    _,_,_,n = _validate_pauli_chain_ops(ops)
    _check_pauli_chain_support(pop.objective,n;context="objective")
    _check_translation_invariance(pop.objective,n)
    _check_reflection_invariance(pop.objective,n)
    _check_pauli_sign_invariance(pop.objective)
    _check_conjugate_invariant_objective(pop.objective)
    _check_pauli_axis_permutation_invariance(pop.objective)
    five = five_site_window == 0 ? fld(n,2) : Int(five_site_window)
    seven = Int(seven_site_window)
    2 <= five <= n && 2 <= seven <= n || throw(ArgumentError("Stationarity windows must lie in 2:nsites."))
    M = NormalMonomial{PauliAlgebra,T}
    reducer = _PauliMomentReducer{M}(n,true,true,true)
    reduce(m) = iszero(_pauli_sign_signature(m)) ? _reduce_moment(reducer,m) : nothing
    domain = Set{M}()
    rows = collect(basis)
    _check_pauli_chain_support(rows,n;context="basis")
    row_set = Set(rows)
    translation_closed = all(_translate_pauli_monomial(row,1,n) in row_set for row in rows)
    # Translate both factors so the first is an orbit representative. This
    # preserves the reduced product word if, and only if, all translated rows
    # remain available. Pair each representative with ALL rows, not a triangle.
    left_rows = translation_closed ? _translation_orbit_representatives(rows,n) : rows
    for (i,row) in enumerate(left_rows), j in (translation_closed ? eachindex(rows) : i:length(rows))
        rep = reduce(last(_mul_term(row,rows[j])))
        isnothing(rep) || push!(domain,rep)
    end
    for patch in subsystems
        sites = sort!(Int[patch...])
        !isempty(sites) && allunique(sites) && all(in(1:n),sites) ||
            throw(ArgumentError("Subsystems must contain distinct chain sites."))
        for code0 in 0:(4^length(sites)-1)
            code, word = code0, T[]
            for site in sites
                axis = code % 4
                code ÷= 4
                axis == 0 || push!(word,T(_pauli_index(site,axis-1)))
            end
            rep = reduce(M(word))
            isnothing(rep) || push!(domain,rep)
        end
    end

    letter(site,axis) = 3*(site-1)+axis
    candidates = Vector{Int}[]
    for i in 2:n-1, j in i+1:n
        push!(candidates,[1,letter(i,2),letter(j,3)])
    end
    for i in 2:five-1, j in i+1:five, k in 2:five-1, l in k+1:five
        allunique((i,j,k,l)) || continue
        for t in 1:3
            push!(candidates,sort([1,letter(i,2),letter(j,3),letter(k,t),letter(l,t)]))
        end
    end
    for i in 2:seven-1, j in i+1:seven, k in 2:seven-3, l in k+1:seven-2, u in l+1:seven-1, v in u+1:seven
        allunique((i,j,k,l,u,v)) || continue
        for t in 1:3
            push!(candidates,sort([1,letter(i,2),letter(j,3),letter(k,t),letter(l,t),letter(u,t),letter(v,t)]))
        end
    end
    for i in 2:seven-1, j in i+1:seven, k in 2:seven-1, l in k+1:seven, u in 2:seven-1, v in u+1:seven
        allunique((i,j,k,l,u,v)) || continue
        for s in 1:2, t in s+1:3
            push!(candidates,sort([1,letter(i,2),letter(j,3),letter(k,s),letter(l,s),letter(u,t),letter(v,t)]))
        end
    end
    unique!(candidates)

    h = _action_polynomial(pop.objective)
    selected = M[]
    equations = typeof(h)[]
    seen = Set{Vector{Tuple{M,_ActionCoefficient}}}()
    unsupported = zero_equations = duplicates = 0
    by_site = _pauli_hamiltonian_terms_by_site(h,n)
    relevant,marks = Int[],zeros(Int,length(h.terms))
    generation = 0
    for raw in candidates
        q = M(T.(raw))
        generation = _pauli_relevant_hamiltonian_terms!(relevant,marks,generation,q,by_site)
        p = _pauli_linear_state_polynomial(q,h,reducer,by_site,domain,true,typeof(h);
            atol=0,relevant_terms=relevant,support_before_reduction=true)
        if isnothing(p)
            unsupported += 1
            continue
        end
        coefficients = Dict{M,_ActionCoefficient}(m=>c for (c,m) in terms(p))
        if isempty(coefficients)
            zero_equations += 1
            continue
        end
        key = sort!([(m,c) for (m,c) in coefficients]; by=first)
        real(last(first(key))) < 0 && (key = [(m,-c) for (m,c) in key])
        if key in seen
            duplicates += 1
            continue
        end
        push!(seen,key)
        push!(selected,M(T.(raw)))
        push!(equations,Polynomial([(c,m) for (m,c) in coefficients]))
    end
    return (operators=selected,equations=equations,moments=sort!(collect(domain)),
        candidate_count=length(candidates),unsupported=unsupported,
        zero_equations=zero_equations,duplicates=duplicates,domain_size=length(domain))
end

# Full physical symmetry pipeline for rotation-invariant periodic Pauli chains.
# The caller supplies the SAME basis, patches and perturbations as the control.
# No singlet-sector constraint is imposed. This does not claim exhaustive
# discovery of Hamiltonian-preserving Clifford transformations.
function pauli_chain_strengthened_relaxation(pop::PolyOpt{PauliAlgebra,T}, ops;
    basis, subsystems=Vector{Int}[], perturbations=eltype(basis)[],
    stationarity::Bool=true, compress::Bool=true, max_block_size::Integer=256,
    stationarity_operators=perturbations, moment_domain::Symbol=:unchanged,
    subsystem_symmetry::Symbol=:none) where {T}
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) && isempty(pop.moment_eq_constraints) ||
        throw(ArgumentError("Chain strengthening requires an unconstrained Hamiltonian."))
    subsystem_symmetry in (:none,:magnetization) || throw(ArgumentError(
        "subsystem_symmetry must be :none or :magnetization."))
    _,_,_,n = _validate_pauli_chain_ops(ops)
    subsystem_symmetry == :magnetization && _check_pauli_u1_invariance(pop.objective,ops[3])
    matrix = _PauliMomentMatrix(collect(basis))
    # Reversing a Pauli product changes only its phase, not its canonical word.
    # Traverse one triangle to retain the SAME original moment domain, without
    # materializing the full exact matrix before translation reduction.
    words = Set([one(eltype(basis))])
    for i in eachindex(basis), j in i:length(basis)
        push!(words,last(_mul_term(basis[i],basis[j])))
    end
    h = _action_polynomial(pop.objective)
    original = (objective=h,constraints=[(:HPSD,matrix)],words=sort!(collect(words)),
        positive_rows=[matrix.rows],positive_indices=[1],equalities=typeof(h)[])
    _continuous_coefficients(Dict(m=>i for (i,m) in enumerate(original.words)),h)
    extra = _pauli_strengthening_constraints(pop,original;
        subsystems,perturbations,stationarity,max_block_size,stationarity_operators,
        magnetization=subsystem_symmetry == :magnetization)
    original = merge(original,(words=_pauli_strengthening_domain(original.words,extra,moment_domain),))
    main,report = pauli_translation_invariant_moment_relaxation(pop,ops,maximum(degree,basis);
        basis,axis_permutation_symmetry=true)
    reducer = _PauliMomentReducer{eltype(basis)}(n,true,true,true)
    spin_words = Set(main.total_basis)
    for word in original.words
        iszero(_pauli_sign_signature(word)) || continue
        image = _reduce_moment(reducer,word)
        isnothing(image) || push!(spin_words,image)
    end
    spin = _pauli_chain_spin_section(pop,ops,sort!(collect(spin_words)),reducer;extend=moment_domain == :extend)
    row_set = Set(basis)
    all(m in row_set for a in spin.actions for b in basis for (_,m) in terms(differentiate(a,b))) ||
        throw(ArgumentError("The selected row basis is not closed under spin rotations."))
    domain = Set(keys(spin.index))
    all(m in domain for word in original.words for (_,m) in terms(spin.discrete(_action_polynomial(word)))) ||
        error("The chain reduction omitted a discrete-invariant moment")
    labels = spin.words
    FP = typeof(main.objective)
    section = sparse(transpose(spin.kernel))
    function project(p)
        c = section*_continuous_coefficients(spin.index,spin.discrete(p))
        return convert(FP,Polynomial([(_ActionCoefficient(v),labels[i]) for (i,v) in enumerate(c) if !iszero(v)]))
    end
    # The site ordering of a density matrix is part of its tensor basis.
    # Keep translated patches unless their projected matrices actually agree;
    # retaining equivalent constraints is safe and does not drop a whole orbit.
    for sites in extra.patches
        sort!(mod1.(sites .+ 1,n)) in extra.patches ||
            throw(ArgumentError("Subsystem patches must include every translated window."))
        sort!(n .+ 1 .- sites) in extra.patches ||
            throw(ArgumentError("Subsystem patches must be reflection closed."))
    end
    constraints = Tuple{Symbol,Matrix{FP}}[]
    provenance = Symbol[]
    transforms = []
    max_error = 0.0
    function append_matrix(mat,kind)
        projected = map(project,mat)
        if compress && size(projected,1) <= 64
            split = _pauli_coefficient_blocks(projected)
            max_error = max(max_error,split.error_bound)
            for b in split.blocks
                push!(constraints,(:PSD,b))
                push!(provenance,kind)
            end
            push!(transforms,(kind=kind,projected=projected,rows=split.rows))
        else
            push!(constraints,(:PSD,projected))
            push!(provenance,kind)
            push!(transforms,(kind=kind,projected=projected,rows=[Matrix{Float64}(I,size(mat,1),size(mat,1))]))
        end
    end
    for (_,mat) in main.constraints
        append_matrix(mat,:moment)
    end
    density_seen = Matrix{FP}[]
    for (label,(_,mat)) in zip(extra.labels,extra.constraints)
        label.kind == :subsystem || continue
        realmat,_ = _continuous_real_cone(map(spin.discrete,mat))
        projected = map(project,realmat)
        projected in density_seen && continue
        push!(density_seen,projected)
        append_matrix(realmat,:subsystem)
    end
    if !isempty(extra.qs)
        _check_translation_basis_closure(extra.qs,n)
        _check_reflection_basis_closure(extra.qs,n)
        _check_pauli_axis_permutation_basis_closure(extra.qs)
        qset = Set(extra.qs)
        all(m in qset for a in spin.actions for q in extra.qs for (_,m) in terms(differentiate(a,q))) ||
            throw(ArgumentError("Perturbation rows are not closed under spin rotations."))
        reps = _translation_orbit_representatives(extra.qs,n)
        translated = Dict(q=>[_translate_pauli_monomial(q,r,n)
            for r in 0:(_translation_orbit_length(q,n)-1)] for q in reps)
        pairing = _mirror_pairing(reps,n)
        by_site = _pauli_hamiltonian_terms_by_site(pop.objective,n)
        CP = Polynomial{PauliAlgebra,T,ComplexF64}
        for k in report.momentum_sectors, (_,rows) in _pauli_signature_blocks(
            [q for q in reps if mod(k*length(translated[q]),n) == 0];axis_permutation_symmetry=true)
            isempty(rows) && continue
            mat = _translation_state_optimality_block(rows,k,n,translated,pop.objective,by_site,reducer,CP)
            for (_,block) in _mirror_real_blocks(mat,rows,pairing,k,n,FP;atol=1e-10)
                append_matrix(block,:perturbative)
            end
        end
    end
    for p in extra.equalities
        q = project(p)
        if !iszero(q)
            push!(constraints,(:Zero,reshape([q],1,1)))
            push!(provenance,:stationarity)
        end
    end
    mp = MomentProblem(project(main.objective),constraints,labels,length(labels);real_moments=true)
    return (moment_problem=mp,spin=spin,extra=extra,original=original,transforms=transforms,
        provenance=provenance,coefficient_error=max_error,main_report=report,
        continuous_spin=true,ground_state_only=!isempty(extra.qs) || !isempty(extra.equalities),compressed=compress,
        moment_domain=moment_domain,subsystem_symmetry=subsystem_symmetry,
        arithmetic=:exact_spin_section_numerical_cone_transforms)
end
