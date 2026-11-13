# Exact linear algebra shared by invariant moments and representation splitting.
function _continuous_rref(input)
    a = map(_action_coefficient, input)
    pivots = Int[]
    r = 1
    for col in axes(a, 2)
        k = findfirst(i -> !iszero(a[i, col]), r:size(a, 1))
        isnothing(k) && continue
        k += r - 1
        a[r, :], a[k, :] = copy(a[k, :]), copy(a[r, :])
        a[r, :] ./= a[r, col]
        for i in axes(a, 1)
            (i == r || iszero(a[i, col])) && continue
            a[i, :] -= a[i, col] * a[r, :]
        end
        push!(pivots, col)
        r += 1
        r > size(a, 1) && break
    end
    return a, pivots
end

function _continuous_nullspace(a)
    issparse(a) && return _continuous_nullspace(sparse(a))
    r, pivots = _continuous_rref(a)
    free = setdiff(collect(axes(a, 2)), pivots)
    n = zeros(_ActionCoefficient, size(a, 2), length(free))
    for (j, col) in enumerate(free)
        n[col, j] = 1
        n[pivots, j] = -r[1:length(pivots), col]
    end
    all(iszero, sparse(a) * n) || error("Exact nullspace verification failed")
    return n, free
end

# Sparse echelon elimination keeps only nonzero row coefficients. Back-substitute
# into free coordinates instead of constructing a dense reduced echelon matrix.
function _continuous_nullspace(input::SparseMatrixCSC)
    a = map(_action_coefficient, input)
    rows = sparse(transpose(a))
    pivots = Dict{Int,Dict{Int,_ActionCoefficient}}()
    z = zero(_ActionCoefficient)
    for i in axes(rows,2)
        isempty(nzrange(rows,i)) && continue
        row = Dict(rowvals(rows)[k]=>nonzeros(rows)[k] for k in nzrange(rows,i) if !iszero(nonzeros(rows)[k]))
        while !isempty(row)
            col = minimum(keys(row))
            scale = row[col]
            if !haskey(pivots,col)
                if !isone(scale)
                    for (j,c) in row
                        row[j] = c/scale
                    end
                end
                pivots[col] = row
                break
            end
            for (j,c) in pivots[col]
                value = get(row,j,z)-scale*c
                iszero(value) ? delete!(row,j) : (row[j] = value)
            end
        end
    end
    free = [j for j in axes(a,2) if !haskey(pivots,j)]
    directions = Dict(j=>Dict(i=>one(_ActionCoefficient)) for (i,j) in enumerate(free))
    for col in sort!(collect(keys(pivots));rev=true)
        direction = Dict{Int,_ActionCoefficient}()
        for (j,c) in pivots[col]
            j == col && continue
            for (k,v) in directions[j]
                value = get(direction,k,z)-c*v
                iszero(value) ? delete!(direction,k) : (direction[k] = value)
            end
        end
        directions[col] = direction
    end
    ii,jj,vv = Int[],Int[],_ActionCoefficient[]
    for (i,direction) in directions, (j,c) in direction
        push!(ii,i); push!(jj,j); push!(vv,c)
    end
    n = sparse(ii,jj,vv,size(a,2),length(free))
    all(iszero,a*n) || error("Exact nullspace verification failed")
    return n,free
end

function _continuous_columns(a)
    _, pivots = _continuous_rref(a)
    return a[:, pivots]
end

# Reynolds projection: projection onto fixed vectors along generator images.
# For compact representations these are complementary subspaces. This also
# works in a nonorthogonal frame; no Euclidean symmetrization is used.
function _continuous_projector(generators)
    fixed, _ = _continuous_nullspace(reduce(vcat, generators))
    dual, _ = _continuous_nullspace(reduce(vcat, transpose.(generators)))
    pairing = transpose(dual) * fixed
    size(pairing, 1) == size(pairing, 2) || throw(ArgumentError(
        "Generator images do not complement the fixed space; compact projection unavailable."))
    projector = sparse(fixed) * sparse(inv(Matrix(pairing))) * sparse(transpose(dual))
    sparse(projector) * projector == projector || error("Exact projector verification failed")
    all(all(iszero, sparse(k) * projector) && all(iszero, sparse(projector) * k) for k in generators) ||
        error("Exact invariant projector verification failed")
    return projector
end

struct InvariantMomentMap{M}
    words::Vector{M}
    index::Dict{M,Int}
    equations::SparseMatrixCSC{_ActionCoefficient,Int}
    reconstruction::SparseMatrixCSC{_ActionCoefficient,Int}
    free_indices::Vector{Int}
    generators::Vector{SparseMatrixCSC{_ActionCoefficient,Int}}
end

function _continuous_coefficients(index, p)
    v = spzeros(_ActionCoefficient, length(index))
    for (c, m) in terms(_action_polynomial(p))
        haskey(index, m) || throw(ArgumentError("The supplied moment space is not closed: missing $m."))
        v[index[m]] += c
    end
    return v
end

function _continuous_coefficients(index, polynomials::AbstractVector{<:Polynomial})
    ii,jj,values = Int[],Int[],_ActionCoefficient[]
    for (j,p) in enumerate(polynomials), (c,m) in terms(p)
        haskey(index,m) || throw(ArgumentError("The supplied moment space is not closed: missing $m."))
        push!(ii,index[m]); push!(jj,j); push!(values,_action_coefficient(c))
    end
    return sparse(ii,jj,values,length(index),length(polynomials))
end

# The map is shared by the objective and every cone. It keeps ordered words
# distinct: adjoint/real-section equations are a separate, subsequent step.
function invariant_moment_map(actions::AbstractVector{<:SuppliedDerivation}, words::AbstractVector{M};
    clifford_generators::AbstractVector{<:CliffordSymmetry}=CliffordSymmetry[]) where {M<:NormalMonomial}
    isempty(actions) && throw(ArgumentError("Supply at least one derivation."))
    allunique(words) || throw(ArgumentError("The moment word basis must be independent."))
    isempty(words) && throw(ArgumentError("The moment word basis is empty."))
    basis = collect(words)
    index = Dict(m => i for (i, m) in enumerate(basis))
    haskey(index, one(first(words))) || throw(ArgumentError("The moment space must contain the identity."))
    for action in actions
        _validate_derivation(action)
    end
    generators = [_continuous_coefficients(index, differentiate.(Ref(a),basis)) for a in actions]
    for g in clifford_generators
        images = [_act_polynomial(g, _action_polynomial(m)) - _action_polynomial(m) for m in basis]
        push!(generators, _continuous_coefficients(index, images))
    end
    equations = sparse(reduce(vcat, transpose.(generators)))
    n, free = _continuous_nullspace(equations)
    return InvariantMomentMap(basis, index, equations, sparse(n), free, generators)
end

function reconstruct_continuous_moments(map::InvariantMomentMap, values::AbstractVector)
    length(values) == size(map.reconstruction, 2) || throw(DimensionMismatch("Wrong invariant-coordinate count."))
    y = map.reconstruction * values
    return Dict(m => y[i] for (i, m) in enumerate(map.words))
end

function continuous_moment_coefficients(map::InvariantMomentMap, p)
    return transpose(map.reconstruction) * _continuous_coefficients(map.index, p)
end

# Exact real coordinates for the chosen relaxation. Complex algebras impose
# y(p†)=conj(y(p)); the existing real-algebra relaxation additionally imposes
# Im(y)=0. This is not an extra physical sector restriction.
function _continuous_real_section(map::InvariantMomentMap{M}, real_section::Bool) where {A,T,M<:NormalMonomial{A,T}}
    n = size(map.reconstruction,2)
    eye = sparse(Matrix{_ActionCoefficient}(I, n, n))
    adj = _continuous_coefficients(map.index, [_action_adjoint(_action_polynomial(m)) for m in map.words])
    # Adjoint compatibility acts on the already invariant coordinates. Their
    # free-word rows are the identity, so no second word-sized nullspace is needed.
    star = (transpose(adj)*map.reconstruction)[map.free_indices,:]
    equations = [real.(star)-real.(eye) -imag.(star); imag.(star) real.(star)+real.(eye)]
    real_section && (equations = vcat(equations, hcat(imag.(map.reconstruction),real.(map.reconstruction))))
    w, _ = _continuous_nullspace(equations)
    all(isreal, w) || error("Real-section elimination produced complex coordinates")
    b = map.reconstruction*(w[1:n, :] + im * w[n+1:end, :])
    transpose(adj)*b == conj.(b) || error("Exact adjoint-section reconstruction failed")
    real_section && !all(isreal,b) && error("Exact real-section reconstruction failed")
    identity = map.index[one(first(map.words))]
    pivot = findfirst(!iszero, b[identity, :])
    isnothing(pivot) && throw(ArgumentError("Invariant section cannot normalize the identity."))
    # Normalize by a small sparse coordinate change, avoiding a many-argument
    # concatenation whose type promotion grows with the invariant dimension.
    normalizer = spzeros(_ActionCoefficient,size(b,2),size(b,2))
    normalizer[pivot,1] = inv(b[identity,pivot])
    column = 2
    for j in axes(b,2)
        j == pivot && continue
        normalizer[j,column] = 1
        normalizer[pivot,column] = -b[identity,j]/b[identity,pivot]
        column += 1
    end
    b = b*normalizer
    # The real slice must itself be invariant before averaging can be used.
    if real_section
        all(all(isreal, k) for k in map.generators) || throw(ArgumentError(
            "The supplied action does not preserve the chosen real moment section."))
    end
    return b
end

function _continuous_row_metric(generators)
    n = size(first(generators), 1)
    eye = Matrix{_ActionCoefficient}(I, n, n)
    # An orthonormal row representation already has a positive invariant
    # metric. Verify it exactly before forming the n²-dimensional metric action.
    all(k' == -k for k in generators) && return eye
    metric_actions = [kron(eye, k') + kron(transpose(k), eye) for k in generators]
    h = reshape(_continuous_projector(metric_actions) * vec(eye), n, n)
    _action_metric(h, n)
    all(all(iszero, k' * h + h * k) for k in generators) || error("Invariant row metric failed verification")
    return h
end

function _continuous_intertwiners(left, right)
    m, n = size(first(left), 1), size(first(right), 1)
    imat = sparse(Matrix{_ActionCoefficient}(I, m, m))
    jmat = sparse(Matrix{_ActionCoefficient}(I, n, n))
    equations = reduce(vcat, [kron(jmat, sparse(l)) - kron(sparse(transpose(r)), imat) for (l,r) in zip(left,right)])
    return first(_continuous_nullspace(equations))
end

function _continuous_restrict(generators, basis, metric)
    basis, metric = sparse(basis), sparse(metric)
    left = sparse(inv(Matrix(basis' * metric * basis))) * basis' * metric
    restricted = [Matrix(left * sparse(k) * basis) for k in generators]
    all(sparse(k) * basis == basis * sparse(r) for (k,r) in zip(generators,restricted)) ||
        error("Invariant-subspace reconstruction failed")
    return restricted
end

# Split a cyclic module. A one-dimensional complex commutant certifies absolute
# irreducibility. Otherwise try self-adjoint commutant elements. Floating point
# proposes rational eigenvalues only; exact eigenspaces must span the whole
# module before any split is accepted. Algebraic-number splits are unsupported.
function _continuous_irreducibles(generators, metric)
    n = size(metric, 1)
    commutant = _continuous_intertwiners(generators, generators)
    size(commutant, 2) == 1 && return [Matrix{_ActionCoefficient}(I, n, n)]
    eye = Matrix{_ActionCoefficient}(I, n, n)
    for j in axes(commutant, 2), phase in (1, im)
        c = phase * reshape(commutant[:, j], n, n)
        s = c + inv(metric) * c' * metric
        # A scalar element does not split the module.
        all(iszero, s - s[1,1] * eye) && continue
        candidates = unique(rationalize.(BigInt, real.(eigvals(ComplexF64.(s))); tol=1e-8))
        spaces = [first(_continuous_nullspace(s - value * eye)) for value in candidates]
        filter!(v -> size(v, 2) > 0, spaces)
        length(spaces) > 1 && sum(v -> size(v, 2), spaces) == n || continue
        pieces = Matrix{_ActionCoefficient}[]
        for v in spaces
            restricted = _continuous_restrict(generators, v, metric)
            append!(pieces, [v * w for w in _continuous_irreducibles(restricted, v' * metric * v)])
        end
        return pieces
    end
    throw(ArgumentError("Exact representation splitting over complex rationals is unsupported for this module; no numerical eigenspace is accepted as exact."))
end

function _continuous_coordinate_pieces(generators, h)
    n = size(h, 1)
    eye = sparse(Matrix{_ActionCoefficient}(I, n, n))
    used = zeros(_ActionCoefficient, n, 0)
    remaining = copy(eye)
    pieces = Matrix{_ActionCoefficient}[]
    # Closing coordinate seeds avoids constructing a full commutant for simple
    # repeated modules, without assuming a group name or a spin representation.
    for j in 1:n
        seed = Matrix(remaining[:, j:j])
        all(iszero, seed) && continue
        v = seed
        while true
            expanded = _continuous_columns(hcat(v, [k * v for k in generators]...))
            size(expanded, 2) == size(v, 2) && break
            v = expanded
        end
        restricted = _continuous_restrict(generators, v, h)
        append!(pieces, [v * w for w in _continuous_irreducibles(restricted, v' * h * v)])
        used = hcat(used, v)
        # The new invariant module is H-orthogonal to all previous modules.
        # Subtract its projector once, rather than invert the growing used space
        # again for every coordinate seed (including seeds already exhausted).
        metric_rows = sparse(v)' * h
        remaining -= sparse(v) * sparse(inv(Matrix(metric_rows * v))) * metric_rows
    end
    size(used, 2) == n || error("Representation decomposition is incomplete")
    all(iszero, remaining) || error("Invariant-subspace projectors do not resolve the identity")
    return pieces
end

function _continuous_decomposition(generators)
    generators = sparse.(generators)
    n = size(first(generators), 1)
    graph = SimpleGraph(n)
    for k in generators
        rows, cols, values = findnz(k)
        for (i, j, value) in zip(rows, cols, values)
            iszero(value) || add_edge!(graph, i, j)
        end
    end
    components = connected_components(graph)
    if length(components) > 1
        # Disconnected coordinates are invariant under every generator. Solve
        # each distinct representation once; identical blocks need no new split.
        h = spzeros(_ActionCoefficient, n, n)
        pieces = Matrix{_ActionCoefficient}[]
        cache = Dict{Vector{Matrix{_ActionCoefficient}},Any}()
        for indices in components
            local_generators = [Matrix(k[indices, indices]) for k in generators]
            local_decomposition = get!(cache, local_generators) do
                _continuous_decomposition(local_generators)
            end
            h[indices, indices] = local_decomposition.metric
            for group in local_decomposition.groups, v in group.copies
                embedded = zeros(_ActionCoefficient, n, size(v, 2))
                embedded[indices, :] = v
                push!(pieces, embedded)
            end
        end
    else
        h = sparse(_continuous_row_metric(generators))
        pieces = _continuous_coordinate_pieces(generators, h)
    end
    eye = sparse(Matrix{_ActionCoefficient}(I, n, n))
    groups = []
    for v in pieces
        current = _continuous_restrict(generators, v, h)
        found = false
        for group in groups
            size(v, 2) == size(first(group.copies), 2) || continue
            intertwiners = _continuous_intertwiners(current, group.generators)
            size(intertwiners, 2) == 0 && continue
            size(intertwiners, 2) == 1 || error("Module was not absolutely irreducible")
            t = reshape(intertwiners[:, 1], size(v, 2), size(v, 2))
            size(_continuous_columns(t), 2) == size(t, 2) || error("Singular irreducible intertwiner")
            push!(group.copies, v * t)
            found = true
            break
        end
        found || push!(groups, (generators=current, copies=[v], metric=v' * h * v))
    end
    basis = hcat([hcat(g.copies...) for g in groups]...)
    # The irreducible copies are mutually H-orthogonal. Invert their small
    # metrics separately, then verify the full inverse exactly as before.
    inverse_rows = SparseMatrixCSC{_ActionCoefficient,Int}[]
    for group in groups, v in group.copies
        w = sparse(v)
        push!(inverse_rows, sparse(inv(Matrix(w' * h * w))) * w' * h)
    end
    inverse = reduce(vcat,inverse_rows)
    inverse * sparse(basis) == eye || error("Representation basis inverse failed")
    return (groups=groups, basis=basis, inverse=inverse, metric=h)
end

function _continuous_congruence(mat::AbstractMatrix{P}, s) where {P<:Polynomial}
    # Apply the two factors separately, visiting only nonzero basis entries.
    # The former four-index loop scanned n⁴ coefficients even for sparse bases.
    columns = [findall(!iszero, s[:,j]) for j in axes(s,2)]
    right = [zero(P) for _ in axes(s,1), _ in axes(s,2)]
    for j in axes(s,2), b in columns[j], a in axes(s,1)
        right[a,j] += s[b,j] * mat[a,b]
    end
    out = [zero(P) for _ in axes(s, 2), _ in axes(s, 2)]
    for i in axes(s,2), a in columns[i], j in axes(s,2)
        out[i,j] += conj(s[a,i]) * right[a,j]
    end
    return out
end

function _continuous_trace(mat::AbstractMatrix{P}, z) where {P<:Polynomial}
    size(mat) == size(z) || throw(DimensionMismatch("Gram matrix size mismatch."))
    coefficients = eltype(terms(zero(P)))[]
    for j in axes(mat,2), i in axes(mat,1)
        iszero(z[j,i]) && continue
        for (c,m) in terms(mat[i,j])
            push!(coefficients,(z[j,i]*c,m))
        end
    end
    return Polynomial(coefficients)
end

function _continuous_exact_problem(pop::PolyOpt{A,T}, sparsity) where {A,T}
    P = Polynomial{A,T,_ActionCoefficient}
    constraints = Tuple{Symbol,Matrix{P}}[]
    positive_rows = Vector{NormalMonomial{A,T}}[]
    positive_indices = Int[]
    corr = sparsity.corr_sparsity
    function append_block(poly, rows)
        q = _action_polynomial(poly)
        basis = _action_polynomial.(rows)
        mat = [_action_product(_action_product(_action_adjoint(a),q),b) for a in basis, b in basis]
        cone = poly in pop.eq_constraints ? :Zero : :HPSD
        push!(constraints, (cone, mat))
        if cone != :Zero
            push!(positive_rows, rows)
            push!(positive_indices, length(constraints))
        end
    end
    for (term_spaces, cons) in zip(sparsity.cliques_term_sparsities, corr.clq_cons)
        for (ts, poly) in zip(term_spaces, [one(pop.objective); corr.cons[cons]])
            for rows in ts.block_bases
                append_block(poly, rows)
            end
        end
    end
    for i in corr.global_cons
        append_block(corr.cons[i], [one(NormalMonomial{A,T})])
    end
    _append_parity_constraints!(constraints, A, P)
    if !isempty(pop.moment_eq_constraints)
        rows, degrees = _collect_moment_eq_row_bases(sparsity.cliques_term_sparsities)
        for g in pop.moment_eq_constraints
            for row in _truncate_moment_eq_row_bases(rows, degrees, g)
                q = _action_product(_action_adjoint(_action_polynomial(row)), _action_polynomial(g))
                push!(constraints, (:Zero, reshape([q], 1, 1)))
            end
        end
    end
    objective = _action_polynomial(pop.objective)
    equalities = P[]
    for (cone, mat) in constraints
        cone == :Zero || continue
        for p in mat
            for q in ((one(_ActionCoefficient)/2)*(p + _action_adjoint(p)),
                (one(_ActionCoefficient)/(2im))*(p - _action_adjoint(p)))
                iszero(q) || push!(equalities, q)
            end
        end
    end
    # Like the existing complex lowering, equality tests include both Hermitian
    # components. Their adjoint words can lie outside the positivity products.
    words = sort!(unique([one(NormalMonomial{A,T});
        [m for (_, mat) in constraints for p in mat for (_,m) in terms(p)];
        [m for p in equalities for (_,m) in terms(p)]]))
    # Do not expand an incomplete relaxation to accommodate the objective.
    _continuous_coefficients(Dict(m=>i for (i,m) in enumerate(words)), objective)
    return (; objective, constraints, words, positive_rows, positive_indices, equalities)
end

function _continuous_symbolic(map, section, labels, p)
    c = transpose(section) * _continuous_coefficients(map.index, p)
    return Polynomial([(_ActionCoefficient(c[i]), labels[i]) for i in findall(!iszero,c)])
end

function _continuous_symbolic(map, section, labels, mat::AbstractMatrix{<:Polynomial})
    # A sparse matrix product visits only moment/coordinate incidences. Repeated
    # transposed matrix-vector products scan every invariant coordinate per entry.
    coefficients = sparse(transpose(section))*_continuous_coefficients(map.index,vec(mat))
    P = typeof(Polynomial((one(_ActionCoefficient),first(labels))))
    projected = Matrix{P}(undef,size(mat))
    for j in eachindex(mat)
        projected[j] = Polynomial([(nonzeros(coefficients)[k],labels[rowvals(coefficients)[k]])
            for k in nzrange(coefficients,j)])
    end
    return projected
end

function _continuous_word_polynomial(map, c)
    return Polynomial([(_ActionCoefficient(c[i]), map.words[i]) for i in eachindex(map.words)])
end

function _continuous_real_polynomial(p::Polynomial{A,T}) where {A,T}
    all(isreal(c) for (c,_) in terms(p)) || throw(ArgumentError("Expected real coefficients after selecting real moment coordinates."))
    return Polynomial(Tuple{Rational{BigInt},NormalMonomial{A,T}}[(real(c),m) for (c,m) in terms(p)])
end

function _continuous_real_cone(mat)
    realmat = map(p -> Polynomial([(real(c),m) for (c,m) in terms(p)]), mat)
    imagmat = map(p -> Polynomial([(imag(c),m) for (c,m) in terms(p)]), mat)
    if all(iszero, imagmat)
        return realmat, false
    end
    return [realmat -imagmat; imagmat realmat], true
end

# Finite averaging must preserve the continuous invariant space and permute
# whole positivity blocks. Initially this composition supports unconstrained
# Pauli problems only; localizer/equality orbit lifting is not implicit.
function _continuous_clifford_actions(pop::PolyOpt{A,T}, original, actions, generators;
    max_group_order::Int=typemax(Int)) where {A,T}
    isempty(generators) && return []
    A === PauliAlgebra || throw(ArgumentError("Clifford composition requires Pauli algebra."))
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) &&
        isempty(pop.moment_eq_constraints) || throw(ArgumentError(
            "Finite/continuous composition currently supports unconstrained Hamiltonians only."))
    domain = indices(pop.registry)
    letters = [_action_polynomial(NormalMonomial{A,T}(T[i])) for i in domain]
    for g in generators
        _validate_clifford_symmetry(_convert_clifford_symmetry(T, g), collect(domain))
        _act_polynomial(g, original.objective) == original.objective || throw(ArgumentError(
            "Clifford generator does not preserve the objective exactly."))
        left = [[_act_polynomial(g, a.images[i]) for i in domain] for a in actions]
        right = [[differentiate(a, _act_polynomial(g, x)) for x in letters] for a in actions]
        words = unique([m for family in (left, right) for ps in family for p in ps for (_,m) in terms(p)])
        index = Dict(m => i for (i,m) in enumerate(words))
        l = hcat([reduce(vcat, [_continuous_coefficients(index, p) for p in ps]) for ps in left]...)
        r = hcat([reduce(vcat, [_continuous_coefficients(index, p) for p in ps]) for ps in right]...)
        l, r = vcat(real.(l), imag.(l)), vcat(real.(r), imag.(r))
        # δg and gδ must span the same derivations, with one coefficient vector
        # shared by all letters. Per-letter span checks would be insufficient.
        all(iszero, transpose(first(_continuous_nullspace(transpose(l)))) * r) &&
            all(iszero, transpose(first(_continuous_nullspace(transpose(r)))) * l) ||
            throw(ArgumentError("Clifford generator does not normalize the supplied continuous generator span."))
    end
    group = _enumerate_symmetry_group([_convert_clifford_symmetry(T, g) for g in generators], collect(domain);
        max_group_order)
    lookup = Dict(Tuple(sort(rows)) => idx for (idx, rows) in zip(original.positive_indices, original.positive_rows))
    result = []
    for g in group
        transforms = Dict{Int,Tuple{Int,SparseMatrixCSC{_ActionCoefficient,Int}}}()
        for (idx, rows) in zip(original.positive_indices, original.positive_rows)
            images = [_act_monomial(g, m) for m in rows]
            target = get(lookup, Tuple(sort(last.(images))), 0)
            target != 0 || throw(ArgumentError("Clifford action does not preserve the positivity block family."))
            target_rows = original.positive_rows[findfirst(==(target), original.positive_indices)]
            positions = Dict(m => i for (i,m) in enumerate(target_rows))
            s = sparse([positions[m] for (_,m) in images], collect(eachindex(rows)),
                _ActionCoefficient[first(image) for image in images], length(rows), length(rows))
            transforms[idx] = (target, s)
        end
        push!(result, transforms)
    end
    return result
end

# Integrates with the existing basis selection and JuMP lowering. The synthetic
# self-adjoint letters label real scalar coordinates only, never new physical
# operators. All algebra and reconstruction data retain their original words.
function continuous_moment_relaxation(pop::PolyOpt{A,T,P}, config::SolverConfig,
    actions::AbstractVector{<:SuppliedDerivation}; compact_metric,
    compress::Bool=true, physics::Union{Nothing,PhysicsSparsity}=nothing,
    additional_moments::AbstractVector=NormalMonomial{A,T}[],
    clifford_generators::AbstractVector{<:CliffordSymmetry}=CliffordSymmetry[],
    max_group_order::Int=typemax(Int)) where {A,T,P<:Polynomial{A,T}}
    isnothing(config.symmetry) || throw(ArgumentError("Use clifford_generators for checked finite/continuous composition; other SymmetrySpec options are unsupported."))
    if !_is_complex_problem(A)
        all(isreal(c) for p in [pop.objective; pop.eq_constraints; pop.ineq_constraints; pop.moment_eq_constraints]
            for (c,_) in terms(p)) || throw(ArgumentError("The real-algebra relaxation requires real input coefficients."))
    end
    sparsity = isnothing(physics) ? compute_sparsity(pop, config) :
        compute_sparsity(pop, config, physics; actions, clifford_generators)
    original = _continuous_exact_problem(pop, sparsity)
    original = merge(original,(words=sort!(unique([original.words; additional_moments])),))
    validation = validate_continuous_symmetry(pop, actions;
        positivity_bases=original.positive_rows, compact_metric)
    validation.compactness == :verified || throw(ArgumentError("Exact reduction requires a verified compact action."))
    # Check the actual finite equality tests, including parity and one-sided
    # localizers, not just the defining constraint-polynomial span.
    for action in actions
        _action_coordinates(original.equalities, differentiate.(Ref(action),original.equalities);
            context="Finite relaxation equality-test span")
    end
    finite_actions = _continuous_clifford_actions(pop, original, actions, clifford_generators; max_group_order)
    moment_map = invariant_moment_map(actions, original.words; clifford_generators)
    section = _continuous_real_section(moment_map, !_is_complex_problem(A))
    projector = isempty(original.equalities) ? nothing : _continuous_projector(moment_map.generators)
    reg, (z,) = create_noncommutative_variables([("continuous_coordinate",1:max(1,size(section,2)-1))])
    labels = [one(first(z)); z[1:size(section,2)-1]]
    objective = _continuous_real_polynomial(_continuous_symbolic(moment_map, section, labels, original.objective))
    RP = typeof(objective)
    constraints = Tuple{Symbol,Matrix{RP}}[]
    blocks = []
    original_decompositions = []
    row_generators = Dict{Int,Vector{Matrix{_ActionCoefficient}}}()
    for (idx, rows) in zip(original.positive_indices, original.positive_rows)
        mat = original.constraints[idx][2]
        basis = _action_polynomial.(rows)
        generators = [_action_coordinates(basis, differentiate.(Ref(a),basis); independent=true) for a in actions]
        row_generators[idx] = generators
        projected = _continuous_symbolic(moment_map, section, labels, mat)
        if !compress
            realmat, embedded = _continuous_real_cone(projected)
            push!(constraints, (:PSD, realmat))
            push!(blocks, (original_index=idx, columns=Matrix{_ActionCoefficient}(I,length(rows),length(rows)),
                metric=ones(_ActionCoefficient,1,1), dimension=1, multiplicity=length(rows),
                embedded=embedded, projected=projected, constraint_index=length(constraints)))
            push!(original_decompositions, nothing)
            continue
        end
        decomposition = _continuous_decomposition(generators)
        transformed = _continuous_congruence(projected, decomposition.basis)
        reconstructed = [zero(eltype(transformed)) for _ in axes(transformed,1), _ in axes(transformed,2)]
        offset = 0
        for group in decomposition.groups
            d, m = size(group.metric,1), length(group.copies)
            metric_inverse = inv(group.metric)
            reduced = [(one(_ActionCoefficient)/d) * _continuous_trace(transformed[offset+(i-1)*d+1:offset+i*d,
                offset+(j-1)*d+1:offset+j*d], metric_inverse) for i in 1:m, j in 1:m]
            for i in 1:m, j in 1:m, a in 1:d, b in 1:d
                reconstructed[offset+(i-1)*d+a,offset+(j-1)*d+b] = group.metric[a,b] * reduced[i,j]
            end
            realmat, embedded = _continuous_real_cone(reduced)
            push!(constraints, (:PSD, realmat))
            push!(blocks, (original_index=idx, columns=hcat(group.copies...), metric=group.metric,
                dimension=d, multiplicity=m, embedded=embedded, projected=reduced,
                constraint_index=length(constraints)))
            offset += m*d
        end
        all(iszero, transformed - reconstructed) || error("Coefficientwise multiplicity reconstruction failed")
        push!(original_decompositions, decomposition)
    end
    equality_lifts = []
    for p in original.equalities
        averaged = _continuous_word_polynomial(moment_map, projector * _continuous_coefficients(moment_map.index,p))
        reduced = _continuous_real_polynomial(_continuous_symbolic(moment_map, section, labels, averaged))
        iszero(reduced) && continue
        lift = vec(_action_coordinates(original.equalities,[averaged]; context="Averaged equality-test span"))
        push!(constraints, (:Zero, reshape([reduced],1,1)))
        push!(equality_lifts, lift)
    end
    # Only the solver boundary converts coefficients to Float64. Exact matrices,
    # maps, and original-algebra certificate data remain in the reduction.
    FP = Polynomial{NonCommutativeAlgebra,eltype(first(labels).word),Float64}
    numerical = [(cone, map(p -> convert(FP,p),mat)) for (cone,mat) in constraints]
    mp = MomentProblem(convert(FP,objective), numerical, labels, length(labels); real_moments=true)
    return (moment_problem=mp, original=original, moment_map=moment_map, section=section,
        projector=projector, coordinate_labels=labels, blocks=blocks, equalities=equality_lifts,
        decompositions=original_decompositions, exact_constraints=constraints, exact_objective=objective,
        row_generators=row_generators, compressed=compress, finite_actions=finite_actions,
        sparsity=sparsity, relation_to_input=:equivalent, physical_target=:unchanged,
        compactness=:verified, arithmetic=:complex_rationals, solver_arithmetic=:float64)
end

function continuous_moment_relaxation(pop, config, actions; kwargs...)
    throw(ArgumentError("Continuous reduction supports ordinary PolyOpt polynomials and supplied derivations in one built-in algebra; state/mixed-algebra actions are unsupported."))
end

function reconstruct_continuous_moments(reduction, coordinates::AbstractVector)
    length(coordinates) == size(reduction.section,2) || throw(DimensionMismatch("Wrong real-coordinate count."))
    y = reduction.section * coordinates
    return Dict(m => y[i] for (i,m) in enumerate(reduction.moment_map.words))
end

function _continuous_evaluate(p, index, values)
    return sum(c * values[index[m]] for (c,m) in terms(p); init=zero(eltype(values)))
end

function reconstruct_continuous_matrices(reduction, coordinates::AbstractVector)
    index = Dict(m=>i for (i,m) in enumerate(reduction.coordinate_labels))
    values = map(_action_coefficient, coordinates)
    matrices = Dict{Int,Matrix{_ActionCoefficient}}()
    for (position,idx) in enumerate(reduction.original.positive_indices)
        selected = [b for b in reduction.blocks if b.original_index == idx]
        s = hcat([b.columns for b in selected]...)
        transformed = zeros(_ActionCoefficient,size(s,2),size(s,2))
        offset = 0
        for block in selected
            b = map(p -> _continuous_evaluate(p,index,values),block.projected)
            lifted = kron(b,block.metric)
            range = offset+1:offset+size(lifted,1)
            transformed[range,range] = lifted
            offset += size(lifted,1)
        end
        if reduction.compressed
            inverse = reduction.decompositions[position].inverse
            matrices[idx] = inverse' * sparse(transformed) * inverse
        else
            matrices[idx] = transformed
        end
    end
    return matrices
end

function _continuous_effective_gram(block, gram)
    z = map(_action_coefficient, gram)
    m = block.multiplicity
    if block.embedded
        size(z) == (2m,2m) || throw(DimensionMismatch("Wrong real-embedded Gram size."))
        return z[1:m,1:m] + z[m+1:end,m+1:end] + im*(z[m+1:end,1:m]-z[1:m,m+1:end])
    end
    size(z) == (m,m) || throw(DimensionMismatch("Wrong reduced Gram size."))
    return z
end

# Schur averaging of a Gram matrix uses small multiplicity traces. The checked
# irreducible decomposition supplies the same compact-group projector without
# constructing its n²-by-n² action on all matrix entries.
function _continuous_average_gram(q, decomposition)
    inverse = decomposition.inverse
    transformed = inverse*sparse(q)*inverse'
    basis = sparse(decomposition.basis)
    averaged = spzeros(_ActionCoefficient,size(q)...)
    offset = 0
    for group in decomposition.groups
        d,m = size(group.metric,1),length(group.copies)
        z = [sum(group.metric[b,a]*transformed[offset+(i-1)*d+a,offset+(j-1)*d+b]
            for a in 1:d for b in 1:d if !iszero(group.metric[b,a]); init=zero(_ActionCoefficient))
            for i in 1:m,j in 1:m]
        columns = basis[:,offset+1:offset+d*m]
        averaged += columns*kron(sparse(z),sparse(inv(group.metric)/d))*columns'
        offset += d*m
    end
    return averaged
end

function recover_continuous_grams(reduction, grams::AbstractVector)
    length(grams) == length(reduction.blocks) || throw(DimensionMismatch("Wrong reduced Gram count."))
    lifted = Dict(idx => zeros(_ActionCoefficient,size(reduction.original.constraints[idx][2]))
        for idx in reduction.original.positive_indices)
    for (block, gram) in zip(reduction.blocks, grams)
        z = _continuous_effective_gram(block,gram)
        columns = sparse(block.columns)
        q = columns * kron(sparse(z),sparse(inv(block.metric)/block.dimension)) * columns'
        if !reduction.compressed
            position = findfirst(==(block.original_index),reduction.original.positive_indices)
            decomposition = reduction.decompositions[position]
            if isnothing(decomposition)
                decomposition = try
                    _continuous_decomposition(reduction.row_generators[block.original_index])
                catch err
                    # Moment-only reduction does not require rational irreducible
                    # splitting. Retain the general exact projector when unavailable.
                    err isa ArgumentError || rethrow()
                    nothing
                end
                reduction.decompositions[position] = decomposition
            end
            if isnothing(decomposition)
                n = size(q,1)
                eye = sparse(Matrix{_ActionCoefficient}(I,n,n))
                generators = [kron(eye,sparse(k))+kron(sparse(conj(k)),eye)
                    for k in reduction.row_generators[block.original_index]]
                q = reshape(_continuous_projector(generators)*vec(q),n,n)
            else
                q = _continuous_average_gram(q,decomposition)
            end
        end
        # A real moment section uses the real symmetric part of a Hermitian
        # multiplier. This preserves positivity and removes imaginary identities
        # which hold only on that section, not in the original operator algebra.
        A = typeof(reduction.original.objective).parameters[1]
        _is_complex_problem(A) || (q = _ActionCoefficient.(real.(q)))
        lifted[block.original_index] += q
    end
    if !isempty(reduction.finite_actions)
        averaged = Dict(idx => zero(q) for (idx,q) in lifted)
        for transforms in reduction.finite_actions, (idx,q) in lifted
            target, s = transforms[idx]
            averaged[target] += s * q * s' / length(reduction.finite_actions)
        end
        return averaged
    end
    return lifted
end

function _continuous_is_psd(input)
    a = map(_action_coefficient,input)
    a == a' || return false
    for j in axes(a,1)
        pivot = real(a[j,j])
        pivot < 0 && return false
        if iszero(pivot)
            all(iszero,a[j,j+1:end]) || return false
            continue
        end
        for k in j+1:size(a,1), i in j+1:size(a,1)
            a[i,k] -= a[i,j]*a[j,k]/pivot
        end
    end
    return true
end

# Checks the identity in the original quotient algebra, including averaged
# equality multipliers. No numerical residual or negative eigenvalue is chopped.
# For bounded-letter algebras, the exact coefficient l1 bound supplies a rigorous
# residual correction. For free/bosonic algebras only a zero residual certifies.
function verify_continuous_certificate(reduction, bound, grams::AbstractVector,
    multipliers::AbstractVector=Rational{BigInt}[])
    length(multipliers) == length(reduction.equalities) || throw(DimensionMismatch("Wrong equality-multiplier count."))
    lambda = _action_coefficient(bound)
    isreal(lambda) || throw(ArgumentError("The proposed lower bound must be real."))
    lifted = recover_continuous_grams(reduction,grams)
    objective = reduction.original.objective
    residual = objective - lambda * one(objective)
    for (idx,q) in lifted
        residual -= _continuous_trace(reduction.original.constraints[idx][2],q)
    end
    equality_coefficients = zeros(_ActionCoefficient,length(reduction.original.equalities))
    for (c,lift) in zip(multipliers,reduction.equalities)
        equality_coefficients += _action_coefficient(c)*lift
    end
    for (c,p) in zip(equality_coefficients,reduction.original.equalities)
        residual -= c*p
    end
    # Positivity survives the checked metric congruences, compact averaging and
    # real section. Test the smaller effective Grams first; keep the original
    # test as a fallback when cancellation makes their sum positive.
    psd = all(_continuous_is_psd(_continuous_effective_gram(b,g)) for (b,g) in zip(reduction.blocks,grams)) ||
        all(_continuous_is_psd, values(lifted))
    exact = iszero(residual)
    coefficient_bound = sum(abs(real(c))+abs(imag(c)) for (c,_) in terms(residual); init=big(0)//1)
    A = typeof(objective).parameters[1]
    bounded = A in (PauliAlgebra,ProjectorAlgebra,UnipotentAlgebra)
    correction = exact ? big(0)//1 : bounded ? coefficient_bound : nothing
    certified = psd && !isnothing(correction)
    return (identity_exact=exact, grams_psd_exact=psd, residual=residual,
        residual_coefficient_bound=coefficient_bound, residual_operator_bound=correction,
        certified=certified, certified_bound=certified ? real(lambda)-correction : nothing,
        grams=lifted, equality_coefficients=equality_coefficients)
end

function continuous_nctssos(pop, config::SolverConfig, actions; compact_metric,
    compress::Bool=true, dualize::Bool=false, physics::Union{Nothing,PhysicsSparsity}=nothing,
    clifford_generators::AbstractVector{<:CliffordSymmetry}=CliffordSymmetry[])
    reduction = continuous_moment_relaxation(pop,config,actions;
        compact_metric,compress,physics,clifford_generators)
    if dualize
        sos = sos_dualize(reduction.moment_problem)
        model = sos.model
        set_optimizer(model,config.optimizer)
        set_silent(model)
        optimize!(model)
        status = _check_solver_status(model)
        grams = [value.(q) for q in model[:psd_duals]]
        multipliers = value.(model[:zero_duals])
        bound = objective_value(model)
        certificate = verify_continuous_certificate(reduction,bound,grams,multipliers)
        return (; objective=bound, model, status, reduction, grams, multipliers, certificate)
    end
    model, extract = build_jump_model(reduction.moment_problem)
    set_optimizer(model,config.optimizer)
    set_silent(model)
    optimize!(model)
    status = _check_solver_status(model)
    values = extract()
    K = typeof(reduction.moment_problem.linear.identity)
    # Some invariant coordinates occur only in equality tests annihilated by
    # averaging. They affect no remaining form; choose zero when reconstructing
    # a representative, without adding constraints to the relaxation.
    coordinates = [get(values,_moment_key(K,m),0.0) for m in reduction.coordinate_labels]
    moments = reconstruct_continuous_moments(reduction,coordinates)
    matrices = reconstruct_continuous_matrices(reduction,coordinates)
    return (; objective=objective_value(model), model, status, reduction, coordinates, moments, matrices)
end

# Search only the real linear span of supplied compact candidate derivations.
# This is an exact linear search, not a claim to discover all automorphisms.
function detect_continuous_symmetries(pop::PolyOpt{A,T}, candidates::AbstractVector{<:SuppliedDerivation};
    positivity_bases, compact_metric) where {A,T}
    isempty(candidates) && throw(ArgumentError("Supply a nonempty candidate search space."))
    all(a isa SuppliedDerivation{A,T} for a in candidates) || throw(ArgumentError(
        "Every candidate must use the problem's algebra and index type."))
    domain = indices(pop.registry)
    letters = [_action_polynomial(NormalMonomial{A,T}(T[i])) for i in domain]
    metric = _action_metric(compact_metric,length(letters))
    for action in candidates
        Set(keys(action.images)) == Set(domain) || throw(ArgumentError("Candidate registry mismatch."))
        _validate_derivation(action)
        k = _action_coordinates(letters,differentiate.(Ref(action),letters); independent=true,
            context="Compact candidate letter span")
        all(iszero,k'*metric+metric*k) || throw(ArgumentError("Candidate is outside the supplied compact family."))
    end
    equations = zeros(_ActionCoefficient,0,length(candidates))
    # Build equations for derivatives modulo a preserved polynomial span.
    function residuals(polynomials, span)
        isempty(polynomials) && return zeros(_ActionCoefficient,0,length(candidates))
        result = zeros(_ActionCoefficient,0,length(candidates))
        for p in polynomials
            derivatives = differentiate.(candidates,Ref(p))
            words = unique([m for q in vcat(span,derivatives) for (_,m) in terms(q)])
            index = Dict(m=>i for (i,m) in enumerate(words))
            d = hcat([_continuous_coefficients(index,q) for q in derivatives]...)
            if !isempty(span)
                basis = hcat([_continuous_coefficients(index,q) for q in span]...)
                annihilator = first(_continuous_nullspace(transpose(basis)))
                d = transpose(annihilator)*d
            end
            result = vcat(result,d)
        end
        return result
    end
    equations = vcat(equations,residuals(_action_polynomial.([pop.objective;pop.ineq_constraints]),[]))
    for polynomials in (pop.eq_constraints,pop.moment_eq_constraints)
        span = _action_polynomial.(polynomials)
        equations = vcat(equations,residuals(span,span))
    end
    for rows in positivity_bases
        span = _action_polynomial.(rows)
        equations = vcat(equations,residuals(span,span))
    end
    weights = first(_continuous_nullspace(vcat(real.(equations),imag.(equations))))
    # The dimension counts derivations, not presentations: duplicate or zero
    # supplied candidates must not create spurious symmetry directions.
    image_words = unique([m for a in candidates for p in values(a.images) for (_,m) in terms(p)])
    image_index = Dict(m=>i for (i,m) in enumerate(image_words))
    images = hcat([reduce(vcat, [_continuous_coefficients(image_index,a.images[i]) for i in domain]) for a in candidates]...)
    _, independent = _continuous_rref(vcat(real.(images*weights),imag.(images*weights)))
    weights = weights[:,independent]
    actions = [SuppliedDerivation(pop.registry, Dict(i =>
        sum(weights[a,j]*candidates[a].images[i] for a in eachindex(candidates)) for i in domain))
        for j in axes(weights,2)]
    isempty(actions) || validate_continuous_symmetry(pop,actions;positivity_bases,compact_metric)
    return (generators=actions, coefficients=real.(weights), dimension=size(weights,2),
        search_class=:real_span_of_supplied_compact_generators, arithmetic=:exact)
end

# Complete within Q = sum_i q_i σ_i over registered on-site Pauli letters.
# δ(p) = i[Q,p]/2. Exact binary floating-point inputs are not rounded to a
# nearby symmetric Hamiltonian. Constraints and finite row closure are checked
# by the reduction, not inferred from Hamiltonian invariance here.
function detect_onsite_spin_symmetries(pop::PolyOpt{PauliAlgebra,T}) where {T}
    domain = indices(pop.registry)
    positions = Dict(i => j for (j,i) in enumerate(domain))
    all(T(_pauli_index(_pauli_site(i), a)) in domain for i in domain for a in 0:2) ||
        throw(ArgumentError("On-site spin discovery requires all three Pauli axes at each registered site."))
    h = _action_polynomial(pop.objective)
    _action_zero(h - _action_adjoint(h), "Objective must be Hermitian")
    all(i in domain for (_,m) in terms(h) for i in m.word) ||
        throw(ArgumentError("Hamiltonian contains an unregistered Pauli letter."))
    M = NormalMonomial{PauliAlgebra,T}
    # Only letters on the support of a term contribute. Avoid constructing 3n
    # dense candidate derivations and differentiating every term with each one.
    equations = Dict{M,Dict{Int,_ActionCoefficient}}()
    for (c,m) in terms(h), (k,i) in enumerate(m.word), a in 0:2
        axis = _pauli_type(i)
        a == axis && continue
        word = copy(m.word)
        word[k] = T(_pauli_index(_pauli_site(i), 3-a-axis))
        image = M(word)
        col = positions[T(_pauli_index(_pauli_site(i), a))]
        scale = mod(axis-a,3) == 1 ? -c : c
        row = get!(equations,image,Dict{Int,_ActionCoefficient}())
        row[col] = get(row,col,zero(_ActionCoefficient)) + scale
    end
    ii,jj,vv = Int[],Int[],_ActionCoefficient[]
    for (r,m) in enumerate(sort!(collect(keys(equations)))), (j,c) in equations[m]
        iszero(c) && continue
        push!(ii,r); push!(jj,j); push!(vv,c)
    end
    system = sparse(ii,jj,vv,length(equations),length(domain))
    weights = real.(first(_continuous_nullspace(system)))
    actions = SuppliedDerivation{PauliAlgebra,T}[]
    for j in axes(weights,2)
        images = Dict{T,Polynomial{PauliAlgebra,T,_ActionCoefficient}}()
        for i in domain
            axis = _pauli_type(i)
            images[i] = Polynomial(Tuple{_ActionCoefficient,M}[
                ((mod(axis-a,3) == 1 ? -1 : 1)*weights[positions[T(_pauli_index(_pauli_site(i),a))],j],
                 M(T[_pauli_index(_pauli_site(i),3-a-axis)])) for a in 0:2 if a != axis])
        end
        push!(actions,SuppliedDerivation(pop.registry,images))
    end
    return (; generators=actions, coefficients=weights, letters=domain,
        dimension=size(weights,2), search_class=:onsite_spin, arithmetic=:exact)
end

# Automatic entry point for nonconstant Hamiltonians. SympleQ is a bounded
# candidate search, not an exhaustive detector. Full validation is retained;
# incompatible candidates fail explicitly rather than silently changing rows.
function continuous_moment_relaxation(pop::PolyOpt{PauliAlgebra,T}, config::SolverConfig;
    cycle_strategy::Symbol=:cycle_basis, ga_backend::Symbol=:backtracking,
    max_automorphisms::Int=100_000, max_group_order::Int=4096, kwargs...) where {T}
    max_automorphisms >= 1 || throw(ArgumentError("max_automorphisms must be positive."))
    max_group_order >= 1 || throw(ArgumentError("max_group_order must be positive."))
    discovery = detect_onsite_spin_symmetries(pop)
    finite = _sympleq_clifford_candidates(pop.objective; cycle_strategy, ga_backend, max_automorphisms)
    actions = discovery.generators
    if isempty(actions)
        actions = [SuppliedDerivation(pop.registry,
            Dict(i=>zero(pop.objective) for i in indices(pop.registry)))]
    end
    reduction = continuous_moment_relaxation(pop,config,actions;
        compact_metric=Matrix{Int}(I,length(indices(pop.registry)),length(indices(pop.registry))),
        clifford_generators=finite, max_group_order, kwargs...)
    return merge(reduction,(; discovery, clifford_generators=finite,
        clifford_search=:nonexhaustive, max_automorphisms, max_group_order))
end
