# Physics selects principal submatrices of a fixed dense moment matrix.
# No completion/equality-of-optima claim is needed for intermediate stages.
struct PhysicsSparsity
    max_region_sites::Int
    max_block_size::Int
    strategy::Symbol
    function PhysicsSparsity(; max_region_sites::Int=2, max_block_size::Int=32, strategy::Symbol=:cut)
        max_region_sites >= 1 || throw(ArgumentError("max_region_sites must be positive"))
        max_block_size >= 3 || throw(ArgumentError("max_block_size must be at least 3 for objective coverage"))
        strategy in (:cut, :radius, :connected) || throw(ArgumentError("strategy must be :cut, :radius, or :connected"))
        new(max_region_sites, max_block_size, strategy)
    end
end

# Each canonical Pauli string has norm one. Keep many-site supports intact;
# multiple strings on the same support remain separate weighted interactions.
function hamiltonian_interactions(pop::PolyOpt{PauliAlgebra})
    supports = Vector{Int}[]
    weights = Float64[]
    for (coefficient, word) in terms(pop.objective)
        isempty(word.word) && continue
        weight = Float64(abs(coefficient))
        isfinite(weight) || throw(ArgumentError("interaction weights must be finite"))
        push!(supports, sort!(unique!(Int[_pauli_site(i) for i in word.word])))
        push!(weights, weight)
    end
    sites = sort!(unique!(Int[_pauli_site(i) for i in indices(pop.registry)]))
    return (; sites, supports, weights)
end

# Exact row count for a complete Pauli registry; a cost proxy, not solve time.
_physics_row_count(k, order) = sum(binomial(big(k), j) * big(3)^j for j in 0:min(k, order))

# Geometry-only joint basis. Graph vertices correspond to sorted registry sites;
# edges are unweighted geometric bonds, independent of Hamiltonian coefficients.
# Include every Pauli word whose support lies in at least one radius ball,
# including disconnected supports inside a ball. This is NOT a clique basis of
# the graph power, nor a decomposition into separate regional moment matrices.
# Reject oversized balls/unions without silently truncating words or symmetry.
function physics_moment_basis(reg::VariableRegistry{PauliAlgebra,T}, geometry::AbstractGraph;
    radius::Integer=1, max_rows::Integer=100_000) where {T}
    radius >= 0 || throw(ArgumentError("Geometry radius must be nonnegative."))
    max_rows >= 1 || throw(ArgumentError("max_rows must be positive."))
    is_directed(geometry) && throw(ArgumentError("Lattice geometry must be undirected."))
    sites = sort!(unique!(Int[_pauli_site(i) for i in indices(reg)]))
    nv(geometry) == length(sites) || throw(ArgumentError(
        "Geometry vertices must correspond to the sorted registry sites."))
    for site in sites, axis in 0:2
        _pauli_index(site, axis) in indices(reg) || throw(ArgumentError(
            "Geometry selection requires all three Pauli axes at each site."))
    end
    M = NormalMonomial{PauliAlgebra,T}
    basis = Set{M}([one(M)])
    seen_balls = Set{Tuple{Vararg{Int}}}()
    for center in vertices(geometry)
        ball, frontier = Set([Int(center)]), [Int(center)]
        for _ in 1:min(radius, length(sites) - 1)
            next = Set{Int}()
            for site in frontier, adjacent in neighbors(geometry, site)
                adjacent in ball || push!(next, adjacent)
            end
            isempty(next) && break
            union!(ball, next)
            frontier = collect(next)
        end
        region = sites[sort!(collect(ball))]
        key = Tuple(region)
        key in seen_balls && continue
        push!(seen_balls, key)
        big(4)^length(region) <= max_rows || throw(ArgumentError(
            "A complete geometry ball exceeds max_rows=$max_rows; reduce radius or increase the budget."))
        local_indices = T[i for i in indices(reg) if _pauli_site(i) in region]
        local_registry = subregistry(reg, local_indices)
        for word in get_ncbasis(local_registry, length(region))
            push!(basis, word)
            length(basis) <= max_rows || throw(ArgumentError(
                "The joint geometry basis exceeds max_rows=$max_rows; no truncated basis was returned."))
        end
    end
    return sort!(collect(basis))
end

# Sparse rows of ONE joint moment matrix, rather than separate regional cones.
# On a chain, connected supports recover contiguous words. On higher-valence
# graphs, retain connected triples and induced four-cycles (plaquettes on a
# square lattice). Graph-distance pairs add separated supports without a list
# of lattice-specific coordinates. Every selected support gets all Pauli axes.
function physics_moment_basis(pop::PolyOpt{PauliAlgebra,T}, order::Integer;
    pair_distance::Integer=1, max_rows::Integer=100_000) where {T}
    1 <= order <= 4 || throw(ArgumentError("Physics moment basis supports orders 1 through 4."))
    pair_distance >= 1 || throw(ArgumentError("pair_distance must be positive."))
    max_rows >= 1 || throw(ArgumentError("max_rows must be positive."))
    interactions = hamiltonian_interactions(pop)
    all(length(s) <= 2 for s in interactions.supports) || throw(ArgumentError(
        "Graph-based moment selection currently requires one- or two-site Hamiltonian terms."))
    neighbors = Dict(s => Set{Int}() for s in interactions.sites)
    for support in interactions.supports
        length(support) == 2 || continue
        a, b = support
        push!(neighbors[a], b)
        push!(neighbors[b], a)
    end
    for site in interactions.sites, axis in 0:2
        _pauli_index(site,axis) in indices(pop.registry) || throw(ArgumentError(
            "Physics moment basis requires all three Pauli axes at each site."))
    end
    chain_graph = all(length(adjacent) <= 2 for adjacent in values(neighbors))
    supports = Set{Tuple{Vararg{Int}}}()
    for region in _physics_connected_regions(interactions, Int(order))
        if length(region) <= 3 || chain_graph ||
            all(count(in(region), neighbors[s]) == 2 for s in region)
            push!(supports, Tuple(region))
        end
    end
    if order >= 2
        for source in interactions.sites
            seen = Set([source])
            frontier = [source]
            for _ in 1:min(pair_distance,length(interactions.sites)-1)
                frontier = sort!(collect(setdiff(union((neighbors[s] for s in frontier)...),seen)))
                isempty(frontier) && break
                union!(seen,frontier)
                for target in frontier
                    push!(supports, minmax(source,target))
                end
            end
        end
    end
    count_rows = 1 + sum(big(3)^length(s) for s in supports; init=big(0))
    count_rows <= max_rows || throw(ArgumentError(
        "Selected joint basis needs $count_rows rows, exceeding max_rows=$max_rows; no rows were truncated."))
    M = NormalMonomial{PauliAlgebra,T}
    basis = M[one(M)]
    for support in sort!(collect(supports); by=s -> (length(s),s))
        for code in 0:(3^length(support)-1)
            push!(basis, _pauli_chain_word(T,collect(support),code))
        end
    end
    return sort!(basis)
end

function _physics_regions(interactions, order, max_sites)
    clusters = [[s] for s in interactions.sites]
    regions = copy(clusters)
    while length(clusters) > 1
        owner = Dict(s => i for (i, cluster) in enumerate(clusters) for s in cluster)
        crossing = Dict{Tuple{Int,Int},Float64}()
        for (support, weight) in zip(interactions.supports, interactions.weights)
            touched = sort!(unique!([owner[s] for s in support]))
            for a in 1:length(touched), b in (a + 1):length(touched)
                pair = (touched[a], touched[b])
                crossing[pair] = get(crossing, pair, 0.0) + weight
            end
        end
        # Deterministic zero-weight merges make the endpoint dense even for
        # disconnected Hamiltonians. A hyperedge counts once per candidate cut.
        candidates = isempty(crossing) ? [(1, 2)] : sort!(collect(keys(crossing)))
        best = first(candidates)
        best_score = -Inf
        for (a, b) in candidates
            na, nb = length(clusters[a]), length(clusters[b])
            cost = _physics_row_count(na + nb, order)^3 -
                _physics_row_count(na, order)^3 - _physics_row_count(nb, order)^3
            score = 2get(crossing, (a, b), 0.0) / Float64(cost)
            if score > best_score
                best, best_score = (a, b), score
            end
        end
        a, b = best
        region = sort!(vcat(clusters[a], clusters[b]))
        # Stop, rather than skip an expensive proposal: increasing the cap
        # must reproduce a prefix of exactly the same merge sequence.
        length(region) > max_sites && break
        push!(regions, region)
        clusters[a] = region
        deleteat!(clusters, b)
    end
    return regions
end

# Baseline: weighted shortest-path balls centered at every site. Pair expansion
# is used only for distances, never interpreted as an interaction-norm bound.
# Retain whole distance ties and every smaller ball to keep the family nested.
function _physics_radius_regions(interactions, max_sites)
    sites = interactions.sites
    n = length(sites)
    position = Dict(s => i for (i, s) in enumerate(sites))
    strength = zeros(n, n)
    for (support, weight) in zip(interactions.supports, interactions.weights)
        for a in 1:length(support), b in (a + 1):length(support)
            i, j = position[support[a]], position[support[b]]
            strength[i, j] += weight
            strength[j, i] += weight
        end
    end
    distance = [i == j ? 0.0 : (strength[i, j] > 0 ? inv(strength[i, j]) : Inf)
        for i in 1:n, j in 1:n]
    for k in 1:n, j in 1:n, i in 1:n
        distance[i, j] = min(distance[i, j], distance[i, k] + distance[k, j])
    end
    regions = Vector{Int}[]
    for center in 1:n
        for radius in sort!(unique(distance[center, :]))
            region = sites[findall(d -> d <= radius, distance[center, :])]
            length(region) > max_sites && break
            push!(regions, region)
        end
    end
    return unique!(regions)
end

# Enumerate connected supports from the interaction graph, not a lattice-specific
# list of paths or plaquettes. Complete regional rows make the degree meaningful.
function _physics_connected_regions(interactions, max_sites)
    neighbors = Dict(s => Set{Int}() for s in interactions.sites)
    for support in interactions.supports, a in support, b in support
        a == b || push!(neighbors[a], b)
    end
    regions = [[s] for s in interactions.sites]
    seen = Set{Tuple{Vararg{Int}}}(Tuple(r) for r in regions)
    head = 1
    while head <= length(regions)
        region = regions[head]
        head += 1
        length(region) >= max_sites && continue
        boundary = sort!(collect(setdiff(union((neighbors[s] for s in region)...), region)))
        for site in boundary
            candidate = sort!([region; site])
            key = Tuple(candidate)
            key in seen && continue
            push!(seen, key)
            push!(regions, candidate)
        end
    end
    return regions
end

function _physics_close_rows(rows, actions, max_rows)
    closed = copy(rows)
    seen = Set(closed)
    head = 1
    while head <= length(closed)
        word = closed[head]
        head += 1
        for action in actions, (_, image) in terms(differentiate(action, word))
            image in seen && continue
            length(closed) < max_rows || throw(ArgumentError(
                "Continuous row closure exceeds max_block_size=$max_rows; increase the budget instead of truncating the action."))
            push!(seen, image)
            push!(closed, image)
        end
    end
    return sort!(closed)
end

function compute_sparsity(
    pop::PolyOpt{PauliAlgebra,T,P}, config::SolverConfig, physics::PhysicsSparsity;
    actions::AbstractVector{<:SuppliedDerivation}=SuppliedDerivation[],
    clifford_generators::AbstractVector{<:CliffordSymmetry}=CliffordSymmetry[],
) where {T,P<:Polynomial{PauliAlgebra,T}}
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) &&
        isempty(pop.moment_eq_constraints) || throw(ArgumentError(
            "Physics sparsity currently supports unconstrained Pauli Hamiltonians only."))
    isnothing(config.symmetry) && isnothing(config.moment_basis) &&
        config.cs_algo isa NoElimination && config.ts_algo isa NoElimination ||
        throw(ArgumentError("Physics sparsity requires an order target, no symmetry, and no other sparsity algorithm."))
    order = compute_relaxation_order(pop, config.order)
    order >= 1 && degree(pop.objective) <= 2order ||
        throw(ArgumentError("The target order must be positive and cover the Hamiltonian degree."))
    interactions = hamiltonian_interactions(pop)
    M = NormalMonomial{PauliAlgebra,T}
    identity = one(M)
    blocks = Vector{M}[[identity]]
    # Seed each energy moment with a 2- or 3-row block, even when its support
    # exceeds the region cap. Splitting uses distinct sites, hence no phase.
    for (_, word) in terms(pop.objective)
        split = min(order, length(word.word))
        left = M(word.word[1:split])
        right = M(word.word[(split + 1):end])
        push!(blocks, sort!(unique!([identity, left, right])))
    end
    regions = physics.strategy == :cut ?
        _physics_regions(interactions, order, physics.max_region_sites) :
        physics.strategy == :radius ? _physics_radius_regions(interactions, physics.max_region_sites) :
        _physics_connected_regions(interactions, physics.max_region_sites)
    for region in regions
        region_indices = T[i for i in indices(pop.registry) if _pauli_site(i) in region]
        registry = subregistry(pop.registry, region_indices)
        if physics.strategy == :connected
            _physics_row_count(length(region), order) <= physics.max_block_size || throw(ArgumentError(
                "Complete connected-region rows exceed max_block_size=$(physics.max_block_size); increase the budget or reduce the declared region/order."))
            push!(blocks, get_ncbasis(registry, min(order, length(region))))
            continue
        end
        basis = M[]
        for local_order in 0:min(order, length(region))
            basis = get_ncbasis(registry, local_order)
            length(basis) >= physics.max_block_size && break
        end
        # Prioritize low-degree words, then words supported by strong terms.
        # Ordering is independent of both caps, preserving nested prefixes.
        function priority(word)
            support = Int[_pauli_site(i) for i in word.word]
            weight = sum((w for (s, w) in zip(interactions.supports, interactions.weights)
                if issubset(support, s)); init=0.0)
            return (degree(word), -weight, word)
        end
        sort!(basis; by=priority)
        resize!(basis, min(length(basis), physics.max_block_size))
        push!(blocks, sort!(basis))
    end
    # Close row spaces under continuous actions and the family under finite
    # actions. This is explicit selection, before comparing reduced/plain SDPs.
    blocks = unique!([_physics_close_rows(b, actions, physics.max_block_size) for b in blocks])
    all(degree(m) <= order for b in blocks for m in b) || throw(ArgumentError(
        "Continuous row closure exceeds the declared order."))
    seen = Set{Tuple{Vararg{M}}}(Tuple(b) for b in blocks)
    head = 1
    while head <= length(blocks)
        block = blocks[head]
        head += 1
        for g in clifford_generators
            rows = sort!([last(_act_monomial(g, m)) for m in block])
            rows = _physics_close_rows(rows, actions, physics.max_block_size)
            all(degree(m) <= order for m in rows) || throw(ArgumentError(
                "Clifford row closure exceeds the declared order."))
            all(i in indices(pop.registry) for m in rows for i in m.word) || throw(ArgumentError(
                "Clifford row closure leaves the problem registry."))
            key = Tuple(rows)
            key in seen && continue
            push!(seen, key)
            push!(blocks, rows)
        end
    end
    # A larger principal matrix implies every contained block. Removing only
    # contained blocks saves cones without weakening any earlier constraint.
    unique!(blocks)
    sets = Set.(blocks)
    keep = [i for i in eachindex(blocks) if !any(j ->
        length(blocks[j]) > length(blocks[i]) && issubset(sets[i], sets[j]), eachindex(blocks))]
    blocks = blocks[keep]
    cliques = [sort!(unique!(T[i for word in block for i in word.word])) for block in blocks]
    supports = [init_activated_supp(project_to_clique(pop.objective, clique), P[], block)
        for (clique, block) in zip(cliques, blocks)]
    corr = CorrelativeSparsity{PauliAlgebra,T,P,M,Nothing}(
        cliques, pop.registry, P[], [Int[] for _ in blocks], Int[], blocks,
        [Vector{M}[] for _ in blocks])
    ts = [[TermSparsity(support, [block])] for (support, block) in zip(supports, blocks)]
    return SparsityResult(corr, supports, ts)
end
