# Local Kramers-invariant second-order cones, using existing scalar moments only.
# Encoding: v = x | (z << k), with Y represented by x=z=1.
_soc_y_parity(v, k) = isodd(count_ones((v & ((1 << k) - 1)) & (v >> k)))
_soc_swap(v, k) = (v >> k) | ((v & ((1 << k) - 1)) << k)
_soc_invariant(v, r, k) = !xor(_soc_y_parity(v, k), isodd(count_ones(v & _soc_swap(r, k))))

function _soc_label(v, k)
    return String(["IXZY"[1 + ((v >> j) & 1) + 2 * ((v >> (k + j)) & 1)] for j in 0:(k - 1)])
end

function _soc_code(label::AbstractString, k)
    length(label) == k || throw(ArgumentError("Pauli labels must have $k letters."))
    code = 0
    for (j, letter) in enumerate(label)
        letter in "IXYZ" || throw(ArgumentError("Expected a label over I, X, Y, Z."))
        letter in "XY" && (code |= 1 << (j - 1))
        letter in "ZY" && (code |= 1 << (k + j - 1))
    end
    return code
end

function _soc_cluster(cluster)
    sites = Int.(collect(cluster))
    length(sites) in (3, 4) || throw(ArgumentError("Invariant SOC clusters must have three or four sites."))
    issorted(sites) && allunique(sites) && all(>(0), sites) ||
        throw(ArgumentError("Cluster sites must be positive, distinct, and sorted."))
    return sites
end

function _soc_available(mp::MomentProblem{PauliAlgebra}, cluster)
    k = length(cluster)
    available = Dict{Int,Int}()
    for (idx, key) in enumerate(mp.linear.moments)
        mono = mp.linear.key_to_monomial[key]
        isempty(mono.word) && continue
        code = 0
        local_word = true
        for letter in mono.word
            j = findfirst(==(_pauli_site(letter)), cluster)
            if isnothing(j)
                local_word = false
                break
            end
            axis = _pauli_type(letter)
            axis != 2 && (code |= 1 << (j - 1))
            axis != 0 && (code |= 1 << (k + j - 1))
        end
        local_word && (available[code] = idx)
    end
    return available
end

# In-place Walsh–Hadamard transform. Swapping the x/z halves of its output
# index makes the ordinary binary dot product the symplectic pairing.
function _soc_walsh!(a)
    stride = 1
    while stride < length(a)
        for start in 1:(2 * stride):length(a), j in start:(start + stride - 1)
            u, v = a[j], a[j + stride]
            a[j], a[j + stride] = u + v, u - v
        end
        stride *= 2
    end
    return a
end

function _soc_scores(available, values, k)
    a = zeros(Float64, 1 << (2 * k))
    total = 0.0
    for (v, idx) in available
        weight = abs2(values[idx])
        total += weight
        a[v + 1] = _soc_y_parity(v, k) ? -weight : weight
    end
    _soc_walsh!(a)
    return [(total + a[_soc_swap(r, k) + 1]) / 2 for r in 0:(length(a) - 1)]
end

"""
    add_pauli_invariant_soc!(model, mp, cluster, r; labels=nothing)

Add a second-order cone (SOC), a constraint `norm(y[labels]) ≤ constant`, to
the direct moment-variable model returned by `build_jump_model(mp)`. `mp` must
be an unreduced Pauli moment problem. `cluster` contains three or four sorted
site numbers; strings `r` and `labels` use that site order and the alphabet
`IXYZ`. `r` must contain an odd number of Y letters. By default select every
available, nonidentity Pauli label fixed by `P_r * conjugation`.

No scalar moments or basis words are added. Returns an audit record including
the selected labels, radius, squared radius, and JuMP constraint reference.
The caller must pass the matching model and moment problem, without symmetry
identifications. This pilot supports the default `build_jump_model` lowering.
"""
function add_pauli_invariant_soc!(model, mp::MomentProblem{PauliAlgebra}, cluster, r::AbstractString; labels=nothing)
    sites = _soc_cluster(cluster)
    k = length(sites)
    code = _soc_code(r, k)
    _soc_y_parity(code, k) || throw(ArgumentError("r must have odd Y parity (Θ² = −I)."))
    available = _soc_available(mp, sites)
    selected = isnothing(labels) ? sort!([v for v in keys(available) if _soc_invariant(v, code, k)]) :
        sort!([_soc_code(label, k) for label in labels])
    !isempty(selected) && allunique(selected) || throw(ArgumentError("Select a nonempty set of distinct labels."))
    all(v -> v != 0 && haskey(available, v) && _soc_invariant(v, code, k), selected) ||
        throw(ArgumentError("Selected labels must be available, nonidentity, and invariant under r."))
    y = model[:y_re]
    radius_squared = 2^(k - 1) - 1
    radius = sqrt(Float64(radius_squared))
    constraint = @constraint(model, [radius; [y[available[v]] for v in selected]] in SecondOrderCone())
    return (; cluster=sites, r=String(r), labels=[_soc_label(v, k) for v in selected],
        constant=radius, constant_squared=radius_squared, constraint)
end

# Bounded connected-subset growth on the Hamiltonian's site-interaction graph.
# Restrict growth to correlative sparsity cliques. This is a deterministic
# heuristic, not exhaustive separation over all clusters of a large system.
function _soc_clusters(pop, sparsity, max_clusters)
    adjacency = Dict{Int,Set{Int}}()
    for mono in monomials(pop.objective)
        sites = [_pauli_site(letter) for letter in mono.word]
        for site in sites
            union!(get!(adjacency, site, Set{Int}()), setdiff(sites, [site]))
        end
    end
    cliques = [Set(_pauli_site.(clique)) for clique in sparsity.corr_sparsity.cliques]
    frontier = [[site] for site in sort!(collect(keys(adjacency)))]
    result = Vector{Int}[]
    for size in 2:4
        next = Set{Tuple{Vararg{Int}}}()
        for cluster in frontier
            neighbors = sort!(collect(union((adjacency[site] for site in cluster)...)))
            for site in neighbors
                site in cluster && continue
                candidate = sort!([cluster; site])
                any(clique -> issubset(candidate, clique), cliques) || continue
                push!(next, Tuple(candidate))
                length(next) >= max_clusters && break
            end
            length(next) >= max_clusters && break
        end
        frontier = [collect(c) for c in sort!(collect(next))]
        size >= 3 && append!(result, frontier)
    end
    return result
end

"""
    pauli_invariant_soc_nctssos(pop, config; clusters=nothing,
        max_clusters=64, max_rounds=3, batch_size=4, violation_tolerance=1e-6)

Solve a Pauli moment relaxation, score local invariant SOCs with a symplectic
Walsh–Hadamard transform, add the most violated constraints in limited batches,
and re-solve. Neither the moment order nor the scalar-moment universe changes.
Hamiltonian antiunitary symmetry is not required. Symmetry-reduced models are
not supported. This function accepts ordinary Pauli polynomials only.

Automatic clusters are connected three-/four-site subsets of the Hamiltonian
support graph inside sparsity cliques, with at most `max_clusters` candidates
per size. Supply `clusters` to choose another local search domain. The stopping
test uses squared norm minus squared radius, not the norm itself. Duplicate
selected sets are skipped. `max_rounds=0` solves only the bare relaxation.

Returns the model, numerical primal objective and dual objective (`bound`),
moment map, sparsity, symbolic problem, audit records, per-solve objectives,
and accumulated solve/separation wall times. Solver dual values are numerical
lower-bound estimates, not interval-certified bounds. Inspect solver accuracy
before interpreting small improvements.
"""
function pauli_invariant_soc_nctssos(
    pop::PolyOpt{PauliAlgebra,T,P}, config::SolverConfig;
    clusters=nothing, max_clusters::Int=64, max_rounds::Int=3,
    batch_size::Int=4, violation_tolerance::Real=1e-6,
) where {T,P<:Polynomial{PauliAlgebra}}
    isnothing(config.symmetry) || throw(ArgumentError("Invariant SOC separation requires unreduced scalar moments."))
    max_clusters > 0 || throw(ArgumentError("max_clusters must be positive."))
    sparsity = compute_sparsity(pop, config)
    mp = moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    model, extract = build_jump_model(mp)
    set_optimizer(model, config.optimizer)
    set_silent(model)
    sites = isnothing(clusters) ? _soc_clusters(pop, sparsity, max_clusters) : clusters
    result = separate_pauli_invariant_socs!(model, mp, sites; max_rounds, batch_size, violation_tolerance)
    return merge(result, (; monomap=extract(), sparsity, moment_problem=mp))
end

"""
    separate_pauli_invariant_socs!(model, mp, clusters;
        max_rounds=3, batch_size=4, violation_tolerance=1e-6)

Solve, score, and add local invariant cones to an existing unreduced Pauli
moment-variable model with an attached optimizer. This allows composition with
other valid constraints, such as anticommutation inequalities. Model and moment
problem must match, as for `add_pauli_invariant_soc!`. Existing constraints are
preserved. Each call tracks duplicates among its own newly added constraints.
Returns numerical primal/dual objectives, audit records, objectives after each
solve, and solve/separation wall times. No interval certification is performed.
"""
function separate_pauli_invariant_socs!(model, mp::MomentProblem{PauliAlgebra}, clusters;
    max_rounds::Int=3, batch_size::Int=4, violation_tolerance::Real=1e-6)
    max_rounds >= 0 && batch_size > 0 || throw(ArgumentError("Invalid separation limits."))
    isfinite(violation_tolerance) && violation_tolerance >= 0 || throw(ArgumentError("Invalid violation tolerance."))
    sites = _soc_cluster.(clusters)
    available = [_soc_available(mp, cluster) for cluster in sites]
    records = NamedTuple[]
    objectives = Float64[]
    seen = Set{Tuple}()
    solve_seconds = 0.0
    separation_seconds = 0.0
    for round in 0:max_rounds
        solve_seconds += @elapsed optimize!(model)
        _check_solver_status(model)
        has_values(model) && has_duals(model) || error("Separation requires primal and dual solver results.")
        push!(objectives, objective_value(model))
        round == max_rounds && break
        added = 0
        separation_seconds += @elapsed begin
            values = value.(model[:y_re])
            candidates = Tuple{Float64,Int,Int}[]
            for (i, cluster) in enumerate(sites)
                k = length(cluster)
                scores = _soc_scores(available[i], values, k)
                for r in 0:(length(scores) - 1)
                    _soc_y_parity(r, k) || continue
                    violation = scores[r + 1] - (2^(k - 1) - 1)
                    violation > violation_tolerance && push!(candidates, (violation, i, r))
                end
            end
            sort!(candidates; by=c -> (-c[1], c[2], c[3]))
            for (violation, i, r) in candidates
                k = length(sites[i])
                selected = sort!([v for v in keys(available[i]) if _soc_invariant(v, r, k)])
                signature = (Tuple(sites[i]), Tuple(selected))
                signature in seen && continue
                record = add_pauli_invariant_soc!(model, mp, sites[i], _soc_label(r, k))
                push!(records, merge(record, (; round=round + 1, violation)))
                push!(seen, signature)
                added += 1
                added == batch_size && break
            end
        end
        added == 0 && break
    end
    return (; objective=objective_value(model), bound=dual_objective_value(model), model,
        records, objectives, solve_seconds, separation_seconds)
end
