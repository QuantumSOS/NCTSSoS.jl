"Adapt the result of `solve_ground_state` to the `price_square` callback contract."
function square_pricing_state(result; domain=nothing)
    plan = result.plan
    return (; rows=only(plan.reports.original_rows), plan.table, result.moments,
        plan.constraint_ids, result.bound, result.status, domain)
end

_pricing_existing_terms(t, w, allowed) = allowed ? NC.moment_terms(t,w) : nothing

function _pricing_existing_terms(t::NC.OrbitMoments, w, allowed)
    if haskey(t.cache,w)
        pairs = NC.moment_terms(t,w)  # Guaranteed cache hit; never extend the orbit.
        return allowed || isempty(pairs) ? pairs : nothing
    end
    # An ambient Hermitian Pauli fixed with sign -1 has an exact zero
    # moment even when it is not in the registered domain.
    for g in t.annihilators
        phase, image = NC.act(g,w)
        image == w && phase.k == 0x02 && return ()
    end
    # Missing from the coordinate map: only an explicitly registered zero
    # orbit counts, not an arbitrary zero outside the pricing domain.
    return NC.moment_terms(t,w;registered_only=true)
end

function _pricing_existing_terms(t::Union{NC.ContinuousMoments,NC.SectorMoments}, w, allowed)
    pairs = _pricing_existing_terms(t.base,w,allowed)
    isnothing(pairs) && return nothing
    return ((id,c*factor) for (old,c) in pairs for (id,factor) in t.coordinates[old])
end

_pricing_native_table(t) = t isa Union{NC.OrbitMoments,NC.MomentTable}
_pricing_native_table(t::Union{NC.ContinuousMoments,NC.SectorMoments}) = _pricing_native_table(t.base)

function _pricing_moment(state, w, partial, domain)
    NC.isone(w) && return 1.0 + 0.0im  # The normalized identity needs no lookup.
    allowed = isnothing(domain) || w in domain
    pairs = partial ? _pricing_existing_terms(state.table,w,allowed) : NC.moment_terms(state.table,w)
    isnothing(pairs) && return nothing
    # Continuous elimination can return several contributions to one coordinate.
    coordinates = Dict{Int,Number}()
    for (id,c) in pairs
        coordinates[id] = get(coordinates,id,0) + c
    end
    value = 0.0 + 0.0im
    for (id, c) in coordinates
        iszero(c) && continue
        id == 1 || (id in state.constraint_ids && 1 <= id-1 <= length(state.moments) &&
            isfinite(state.moments[id-1])) || return nothing
        value += ComplexF64(c) * (id == 1 ? 1.0 : state.moments[id-1])
    end
    return isfinite(value) ? value : nothing
end

function _pricing_images(points, L)
    return ([(mod(sx*(swap ? y : x), L), mod(sy*(swap ? x : y), L)) for (x,y) in points]
        for swap in (false,true) for sx in (-1,1) for sy in (-1,1))
end

function _pricing_motif(points, L)
    return minimum(Tuple(sort!([(mod(x-a,L), mod(y-b,L)) for (x,y) in image]))
        for image in _pricing_images(points,L) for (a,b) in image)
end

function _pricing_motif_supports(motif, L)
    return Set(Tuple(sort!([1+mod(x+a,L)+L*mod(y+b,L) for (x,y) in image]))
        for image in _pricing_images(motif,L) for a in 0:L-1 for b in 0:L-1)
end

function _pricing_word_rep(w::W, L) where {W}
    letters = [((mod(i,L),i÷L), a) for (i,a) in (divrem(k-1,3) for k in NC.legacy_letters(w))]
    best = nothing
    for swap in (false,true), sx in (-1,1), sy in (-1,1), ((x0,y0),_) in letters
        image = [(mod(sx*(swap ? y-y0 : x-x0),L),mod(sy*(swap ? x-x0 : y-y0),L),a)
            for ((x,y),a) in letters]
        sort!(image)
        colors = fill(-1,3)
        next = 0
        for (i,(x,y,a)) in enumerate(image)
            if colors[a+1] == -1
                colors[a+1] = next
                next += 1
            end
            image[i] = (x,y,colors[a+1])
        end
        key = Tuple(image)
        (isnothing(best) || key < best) && (best = key)
    end
    p = one(W)
    for (x,y,a) in best
        p = last(NC.mul(p,NC.pauli_letter(W,1+x+L*y,a)))
    end
    return p
end

function _pricing_shell(H, rows, graph, L, max_diameter, partial)
    distances = [Graphs.gdistances(graph,i) for i in Graphs.vertices(graph)]
    # For Pauli words each commutator term is zero or a scalar multiple of
    # the same boundary product, so the product shell already contains its support.
    words = [w for (c,w) in NC.terms(H) if !iszero(c) && !isone(w)]
    W = eltype(rows)
    candidates = Set{W}()
    products = 0
    if partial
        by_site = [W[] for _ in Graphs.vertices(graph)]
        for h in words, i in NC.sites(h)
            push!(by_site[i],h)
        end
        nearby = [findall(<=(1),d) for d in distances]
        for b in rows
            # Every incumbent orbit has its canonical word in the row family.
            # Its nonempty support touches site 1, so test that cheap filter first.
            # The Hamiltonian support is invariant under the same actions.
            # Disjoint, nonadjacent factors cannot leave a connected product.
            sites = collect(NC.sites(b))
            isempty(sites) || (1 in sites && b == _pricing_word_rep(b,L)) || continue
            hs = isempty(sites) ? by_site[1] : Set(h for i in sites for j in nearby[i] for h in by_site[j])
            products += length(hs)
            union!(candidates,(last(NC.mul(b,h)) for h in hs))
        end
    else
        products = length(rows)*length(words)
        union!(candidates,(last(NC.mul(b,h)) for b in rows for h in words))
    end
    setdiff!(candidates,rows)
    motifs = Dict{Tuple{Vararg{Int}},Union{Nothing,Tuple}}()
    result = NamedTuple{(:word,:motif),Tuple{eltype(rows),Tuple}}[]
    seen = Set{W}()
    for p in candidates
        motif = get!(motifs,Tuple(NC.sites(p))) do
            sites = collect(NC.sites(p))
            isempty(sites) && return nothing
            all(distances[i][j] <= max_diameter for i in sites for j in sites) || return nothing
            local_graph, _ = Graphs.induced_subgraph(graph, sites)
            Graphs.is_connected(local_graph) || return nothing
            _pricing_motif([(mod(i-1,L),(i-1)÷L) for i in sites],L)
        end
        isnothing(motif) && continue
        rep = partial ? _pricing_word_rep(p,L) : p
        rep in seen && continue
        push!(seen,rep)
        push!(result,(; word=rep,motif))
    end
    # Group equal supports so the local moment matrix and masks can be reused.
    sort!(result; by=c -> partial ? (Tuple(NC.sites(c.word)),c.word) : ((),c.word))
    return (; candidates=result,products,distances)
end

function _pricing_minors(indices, known)
    isempty(indices) && return Vector{Int}[]
    all(known[i,j] for i in indices for j in indices) && return [indices]
    # A deterministic cover by maximal known cliques, not enumeration of
    # every maximal clique (which can be exponential even in a dense graph).
    uncovered = BitSet(indices)
    minors = Vector{Int}[]
    while !isempty(uncovered)
        clique = [first(uncovered)]
        for j in indices
            j == first(clique) && continue
            all(i -> known[i,j],clique) && push!(clique,j)
        end
        setdiff!(uncovered,clique)
        push!(minors,sort!(clique))
    end
    return minors
end

_pricing_rank(cut) = (cut.kind === :range, cut.kind === :range ? cut.range_residual : cut.schur)

function _pricing_local_pool(sites, by_site, row_sites, distances, radius, limit)
    distance(i) = minimum(distances[s][i] for s in sites)
    nearby = [distance(i) for i in eachindex(by_site)]
    indices = BitSet(i for site in eachindex(nearby) if nearby[site] <= radius for i in by_site[site])
    options = [(0,0,findfirst(isempty,row_sites))]
    for i in indices
        d = maximum(nearby[s] for s in row_sites[i])
        d <= radius && push!(options,(d,length(row_sites[i]),i))
    end
    sort!(options)
    # Restore the incumbent row order for deterministic masks and witnesses.
    return sort!([options[i][3] for i in 1:min(limit,length(options))])
end

function _price_square_state(H, state, L, max_diameter; rank_tol, range_tol, violation_tol, witness_tol,
    partial, local_radius, local_row_limit)
    rows = state.rows
    graph, _ = NC._interaction_graph(NC.Problem(H))
    shell = _pricing_shell(H, rows, graph, L, max_diameter, partial)
    W = eltype(rows)
    moments = Dict{W,Union{Nothing,ComplexF64}}()
    domain = partial && hasproperty(state,:domain) ? state.domain : nothing
    if partial && isnothing(domain) && !_pricing_native_table(state.table)
        throw(ArgumentError("Orbit pricing needs a registered-word domain for a custom moment table."))
    end
    function entry(b, p)
        # Packed Pauli strings are Hermitian. The phase of b†p is not optional.
        phase, w = NC.mul(b,p)
        if partial && length(moments) >= 2^18 && !haskey(moments,w)
            empty!(moments)  # Bound both positive and negative lookup caching.
        end
        v = get!(() -> _pricing_moment(state,w,partial,domain),moments,w)
        return isnothing(v) ? nothing : NC.scale(phase,v)
    end
    function context(pool)
        local_rows = rows[pool]
        M = Matrix{ComplexF64}(undef,length(pool),length(pool))
        known = falses(size(M))
        for j in eachindex(local_rows), i in 1:j
            value = entry(local_rows[i],local_rows[j])
            isnothing(value) && continue
            known[i,j] = known[j,i] = true
            M[i,j] = value
            M[j,i] = conj(value)
        end
        return (; pool, local_rows, M, known, masks=Dict{Vector{Int},Vector{Vector{Int}}}())
    end
    factors = Dict{Vector{Int},LinearAlgebra.Eigen}()
    cached_entries = 0
    function price(candidate, base)
        (; pool, local_rows, M, known, masks) = base
        cross = entry.(local_rows,Ref(candidate.word))
        I = findall(!isnothing,cross)
        invalid, unconfirmed, largest_range, largest_schur = 0, 0, 0.0, -Inf
        best = nothing
        complete = get!(() -> _pricing_minors(I,known),masks,I)
        for J in complete
            indices = pool[J]
            F = get(factors,indices,nothing)
            if isnothing(F)
                F = eigen(Hermitian(M[J,J]))
                # Cache recurring patterns without retaining a dense factor
                # for every distinct mask in a large shell.
                if cached_entries + length(J)^2 <= 2^21
                    factors[indices] = F
                    cached_entries += length(J)^2
                end
            end
            if first(F.values) < -witness_tol
                invalid += 1  # The old principal minor is already indefinite.
                continue
            end
            m = ComplexF64[cross[j] for j in J]
            coeffs = F.vectors' * m
            retained = F.values .> rank_tol * max(1.0,last(F.values))
            residual = LinearAlgebra.norm(coeffs[.!retained])
            largest_range = max(largest_range,residual)
            schur = if residual > range_tol
                Inf  # Check range first; do not use a pseudoinverse on this branch.
            else
                sum(abs2.(coeffs[retained]) ./ F.values[retained]) - 1
            end
            isfinite(schur) && (largest_schur = max(largest_schur,schur))
            schur > violation_tol || continue
            # A rank threshold can classify a small, positive eigenvalue as null.
            # Only report a cut if the original, untruncated major minor fails PSD.
            lambda_min = only(LinearAlgebra.eigvals(Hermitian([M[J,J] m; m' 1.0]),1:1))
            if lambda_min >= -witness_tol
                unconfirmed += 1
                continue
            end
            cut = (; kind=isinf(schur) ? :range : :schur, indices, range_residual=residual,
                schur, lambda_min)
            (isnothing(best) || _pricing_rank(cut) > _pricing_rank(best)) && (best = cut)
        end
        return (; candidate..., local_rows=length(pool), known_cross=length(I), minors=length(complete), invalid,
            unconfirmed, largest_range, largest_schur, cut=best)
    end
    scores = if isempty(shell.candidates)
        NamedTuple[]
    elseif !partial
        base = context(collect(eachindex(rows)))
        map(c -> price(c,base),shell.candidates)
    else
        row_sites = [Tuple(NC.sites(w)) for w in rows]
        by_site = [Int[] for _ in Graphs.vertices(graph)]
        for (i,sites) in enumerate(row_sites), site in sites
            push!(by_site[site],i)
        end
        prior_sites, base = (), nothing
        map(shell.candidates) do candidate
            sites = Tuple(NC.sites(candidate.word))
            if sites != prior_sites
                pool = _pricing_local_pool(sites,by_site,row_sites,shell.distances,local_radius,local_row_limit)
                (isnothing(base) || pool != base.pool) && (base = context(pool))
                prior_sites = sites
            end
            price(candidate,base)
        end
    end
    return (; scores, products=shell.products, largest_pool=maximum(s -> s.local_rows,scores; init=0))
end

"""
    price_square(H, L, solve_motifs; initial_motifs=(), batch_size=2, max_rounds=2, ...)

Run adaptive pricing on the periodic `L × L` square (site `x + L*(y-1)`).
`solve_motifs(motifs)` receives the cumulative tuple of support motifs, each a
tuple of `(dx,dy)` offsets modulo `L`. Return `(; rows, table, moments,
constraint_ids, bound, status)` or use `square_pricing_state` on a workflow
solve. A failed solve must return `bound=nothing`. Use the same exact symmetry
and cone policy and retain all old constraints at each rebuild; the driver
checks row nesting and complete translation/D4/spin support orbits, not cone
equivalence. The seed builder owns its initial rows and any extra partners.

Only complete known principal minors are tested. An incomplete base uses a
deterministic clique cover, so no cut found is not an exhaustive PSD test.
Missing coordinates are not assigned zero. Rank is relative to the largest
incumbent eigenvalue (floored at 1); other pricing tolerances are absolute.
An accepted cut excludes the old solved point, not necessarily its objective.

Use `pricing_mode=:orbit` only if the Hamiltonian, actual row family and solved
moments are invariant under translations, D4 and global spin permutations.
This caller assertion enables anchored shell generation and one representative
per Pauli orbit. Local pools include identity and at most `local_row_limit`
rows (identity included) whose entire support lies within interaction-graph
distance `local_radius` of the candidate. `history[i].pricing.partial` marks
this scope; lack of a cut returns `:partial_no_violation`. A callback may
optionally return `domain` (fast ambient Pauli membership) to restrict existing
coordinate lookups further; custom tables need this field in orbit mode.
`history` records each numerical bound, candidate scores and selected batch.
`max_rounds` counts re-solves, not the initial solve. An optional physical
`energy_upper_bound` enables relative-gap stopping; no bound is certified here.
"""
function price_square(H, L::Integer, solve_motifs; initial_motifs=(), batch_size::Int=2,
    max_rounds::Int=2, max_diameter::Int=2, max_rows::Int=typemax(Int),
    pricing_mode::Symbol=:full, local_radius::Int=1, local_row_limit::Int=64,
    rank_tol::Real=1e-8, range_tol::Real=1e-7, violation_tol::Real=1e-6,
    witness_tol::Real=1e-8, bound_tol::Real=1e-6, energy_upper_bound=nothing, gap_tol::Real=1e-4)
    _check_hamiltonian(H)
    L >= 2 && L^2 == H.registry.nsites || throw(ArgumentError("Expected a square Pauli registry."))
    batch_size > 0 && max_rounds >= 0 && max_diameter >= 0 && max_rows > 0 ||
        throw(ArgumentError("Invalid pricing batch, round, diameter, or row limit."))
    pricing_mode in (:full,:orbit) && local_radius >= 0 && local_row_limit >= 1 ||
        throw(ArgumentError("Invalid pricing mode, local radius, or local row limit."))
    partial = pricing_mode === :orbit
    all(t -> isfinite(t) && t >= 0, (rank_tol,range_tol,violation_tol,witness_tol,bound_tol,gap_tol)) ||
        throw(ArgumentError("Pricing tolerances must be finite and nonnegative."))
    isnothing(energy_upper_bound) || isfinite(energy_upper_bound) || throw(ArgumentError("Nonfinite energy upper bound."))
    all(!isempty,initial_motifs) || throw(ArgumentError("A pricing motif must contain a site."))
    motifs = Tuple(unique(_pricing_motif(unique([(mod(x,L),mod(y,L)) for (x,y) in m]),L) for m in initial_motifs))
    history = NamedTuple[]
    previous = nothing
    for round in 0:max_rounds
        elapsed = @elapsed state = solve_motifs(motifs)
        rows = Set(state.rows)
        if !isnothing(previous)
            previous.rows ⊆ rows || throw(ArgumentError("Pricing callback did not retain the old rows."))
        end
        W = eltype(state.rows)
        one(W) in rows && length(rows) == length(state.rows) &&
            all(i <= L^2 for w in rows for i in NC.sites(w)) || throw(ArgumentError("Invalid pricing rows."))
        counts = Dict{Tuple{Vararg{Int}},Int}()
        for w in rows
            sites = Tuple(NC.sites(w))
            counts[sites] = get(counts,sites,0) + 1
        end
        active = isnothing(previous) ? motifs : previous.selected
        for motif in active, sites in _pricing_motif_supports(motif,L)
            get(counts,sites,0) == big(3)^length(sites) ||
                throw(ArgumentError("Pricing callback omitted a motif orbit member."))
        end
        bound = state.bound
        delta = isnothing(previous) || isnothing(bound) ? nothing : bound - previous.bound
        gap = isnothing(bound) || isnothing(energy_upper_bound) ? nothing :
            (energy_upper_bound - bound) / (iszero(energy_upper_bound) ? 1.0 : abs(energy_upper_bound))
        scores, selected = (), ()
        pricing_seconds = 0.0
        pricing = (; mode=pricing_mode, partial, products=0, largest_pool=0)
        stop_reason = if isnothing(bound) || !isfinite(bound)
            :solve_failed
        elseif !isnothing(delta) && delta < -bound_tol
            :bound_decreased
        elseif !isnothing(gap) && bound > energy_upper_bound + bound_tol
            :upper_bound_inconsistent
        elseif !isnothing(gap) && gap <= gap_tol
            :gap_reached
        elseif length(rows) > max_rows
            :row_limit
        elseif length(state.moments) < maximum(state.constraint_ids; init=1)-1 ||
            !all(id -> id == 1 || isfinite(state.moments[id-1]), state.constraint_ids)
            :moments_unavailable
        else
            pricing_seconds = @elapsed scored = _price_square_state(H,state,L,max_diameter;
                rank_tol,range_tol,violation_tol,witness_tol,partial,local_radius,local_row_limit)
            scores = scored.scores
            pricing = (; mode=pricing_mode, partial, scored.products, scored.largest_pool)
            cuts = filter(s -> !isnothing(s.cut),scores)
            sort!(cuts; by=s->_pricing_rank(s.cut), rev=true)  # Equal scores keep shell order.
            if isempty(cuts)
                if any(s -> s.invalid > 0 || s.unconfirmed > 0,scores)
                    :unreliable_moments
                else
                    partial ? :partial_no_violation : :no_violation
                end
            elseif round == max_rounds
                :round_limit
            else
                batch = Tuple[]
                size_after = length(rows)
                for motif in unique(s.motif for s in cuts)
                    motif in motifs && continue
                    added = sum(big(3)^length(sites) - get(counts,sites,0) for sites in _pricing_motif_supports(motif,L))
                    size_after + added <= max_rows || continue
                    push!(batch,motif)
                    size_after += added
                    length(batch) == batch_size && break
                end
                selected = Tuple(batch)
                isempty(selected) ? :row_limit : :continue
            end
        end
        push!(history,(; round, motifs, rows=copy(state.rows), row_count=length(rows), bound, state.status, delta, gap,
            scores, selected, pricing, solve_seconds=elapsed, pricing_seconds))
        stop_reason === :continue || return (; history, motifs, stop_reason)
        previous = (; rows, bound, selected)
        motifs = (motifs...,selected...)
    end
end
