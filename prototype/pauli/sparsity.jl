# Included by the prototype module, which imports NCTSSoS as NC.
const Graphs = NC.Graphs

_workflow_basis_count(n, d) = sum(big(3)^k * binomial(big(n), k) for k in 0:min(n, d))

function _workflow_power(graph, radius)
    powered = copy(graph)
    for u in Graphs.vertices(graph), v in Graphs.neighborhood(graph, u, radius)
        u < v && Graphs.add_edge!(powered, u, v)
    end
    return powered
end

# Merge overlapping cliques up to a union-size cap (0 disables the cap).
# With just a minimum size, merge only when one input is below the floor.
function _workflow_merge_cliques(cliques; min_size=0, max_union=0)
    cliques = sort!(unique(sort!.(cliques)))
    min_size == max_union == 0 && return cliques
    while true
        pair, overlap = nothing, 0
        for i in eachindex(cliques), j in i+1:length(cliques)
            a, b = cliques[i], cliques[j]
            min_size > 0 && min(length(a), length(b)) >= min_size && continue
            shared = length(intersect(a, b))
            shared > overlap || continue
            max_union > 0 && length(a) + length(b) - shared > max_union && continue
            pair, overlap = (i, j), shared
        end
        isnothing(pair) && return sort!(cliques)
        i, j = pair
        merged = sort!(union(cliques[i], cliques[j]))
        cliques[i] = merged
        deleteat!(cliques, j)
        filter!(c -> c === merged || !issubset(c, merged), cliques)
    end
end

function _workflow_lattice_regions(shape, d, range_r)
    n = prod(shape)
    d == 0 && return [Int[]]
    regions = Set{Vector{Int}}([Int[i] for i in 1:n])
    d == 1 && return sort!(collect(regions))
    nx, ny = first(shape), length(shape) == 1 ? 1 : last(shape)
    site(x, y) = mod1(x, nx) + (mod1(y, ny) - 1) * nx
    if length(shape) == 1
        for start in 1:n, width in 2:d
            push!(regions, sort!([mod1(start + offset, n) for offset in 0:width-1]))
        end
    else
        for x in 1:nx, y in 1:ny
            square = [site(x,y), site(x+1,y), site(x,y+1), site(x+1,y+1)]
            for pair in ((1,2), (1,3))
                push!(regions, sort!(unique(square[collect(pair)])))
            end
            if d >= 3
                for k in 1:4
                    push!(regions, sort!(unique([square[j] for j in 1:4 if j != k])))
                end
                d >= 4 && push!(regions, sort!(unique(square)))
            end
        end
    end
    # Graph-distance pairs on the periodic nearest-neighbor lattice.
    for x in 1:nx, y in 1:ny
        u = site(x, y)
        dxmax = min(range_r, nx ÷ 2)
        for dx in -dxmax:dxmax
            dymax = min(range_r - abs(dx), ny ÷ 2)
            for dy in -dymax:dymax
                v = site(x + dx, y + dy)
                u < v && push!(regions, [u, v])
            end
        end
    end
    return sort!(collect(regions))
end

# Use one visitation array, not a new length-n array for each local BFS.
function _workflow_distance_pairs!(add, ::Type{W}, graph, radius) where {W}
    visited = zeros(Int, Graphs.nv(graph))
    queue = Tuple{Int,Int}[]
    for root in Graphs.vertices(graph)
        empty!(queue)
        push!(queue, (root, 0))
        visited[root] = root
        i = 1
        while i <= length(queue)
            u, depth = queue[i]
            i += 1
            depth == radius && continue
            for v in Graphs.neighbors(graph, u)
                visited[v] == root && continue
                visited[v] = root
                push!(queue, (v, depth + 1))
                root < v || continue
                for a in 0:2, b in 0:2
                    add(last(NC.mul(NC.pauli_letter(W, root, a), NC.pauli_letter(W, v, b))))
                end
            end
        end
    end
end

# The footprint is the union of factor supports, not the support after Pauli
# cancellations. Two prefixes with the same footprint and word have identical
# continuations; the shallower one already explores every eligible suffix.
function _workflow_local_products!(add, words::Vector{W}, n, d, depth) where {K,W<:NC.PauliWord{K}}
    masks = [w.x .| w.z for w in words]
    incident = [Int[] for _ in 1:n]
    for (i, mask) in enumerate(masks), s in NC.MaskSites(mask)
        push!(incident[s], i)
    end
    frontier = collect(zip(masks, words))
    seen = Set(frontier)
    neighbors = Set{Int}()
    products = 0
    for level in 1:depth
        isempty(frontier) && break
        next = empty(frontier)
        for (mask, word) in frontier
            NC.degree(word) <= d && add(word)
            level == depth && continue
            empty!(neighbors)
            for site in NC.MaskSites(mask)
                union!(neighbors, incident[site])
            end
            for i in neighbors
                products += 1
                state = (mask .| masks[i], last(NC.mul(word, words[i])))
                state in seen && continue
                push!(seen, state)
                push!(next, state)
            end
        end
        frontier = next
    end
    return products
end

# Balanced factorizations of a requested word only. Include both factors of
# every cut so the coverage augmentation commutes with relabelings. No global
# degree-d enumeration, and no query against a materialized Gram matrix.
function _workflow_coverage!(add, w::W, d) where {W}
    degree = NC.degree(w)
    degree <= d && return add(w)
    degree <= 2d || throw(ArgumentError("Requested moment $w exceeds row-product degree $(2d); increase order."))
    letters = W[l for l in NC.word_letters(w)]
    function cuts(left, start, remaining)
        if remaining == 0
            add(left)
            add(last(NC.mul(left, w)))
            return
        end
        for i in start:length(letters)-remaining+1
            cuts(last(NC.mul(left, letters[i])), i + 1, remaining - 1)
        end
    end
    cuts(one(W), 1, degree ÷ 2)
end

function _workflow_physics_support(H, observable, ::Type{W}, mode, d, radius, depth, generators, max_rows) where {W}
    registry = H.registry
    n = registry.nsites
    rows, seen = W[], Set{W}()
    considered = 0
    function add(w)
        considered += 1
        w in seen && return
        length(rows) < max_rows || throw(ArgumentError("Pauli basis exceeds max_rows=$max_rows."))
        push!(seen, w)
        push!(rows, w)
        return
    end
    add(one(W))
    if d >= 1
        for i in 1:n, a in 0:2
            add(NC.pauli_letter(W, i, a))
        end
    end
    products = 0
    if mode === :graph_distance && d >= 2
        graph, _ = NC._interaction_graph(NC.Problem(H))
        _workflow_distance_pairs!(add, W, graph, radius)
    elseif mode === :local_products && d >= 1
        words = W[w for (c, w) in NC.terms(H) if !iszero(c) && !isone(w)]
        products = _workflow_local_products!(add, words, n, d, depth)
    end
    # Protect the original physics stencils from unrelated observable factors;
    # e.g. a distant query must not implicitly select a distant RDM region.
    regions = sort!(collect(Set(collect(Int, NC.sites(w)) for w in rows if !isone(w))))
    isempty(regions) && push!(regions, Int[])
    for p in (isnothing(observable) ? (H,) : (H, observable)), (c, w) in NC.terms(p)
        iszero(c) || _workflow_coverage!(add, w, d)
    end
    i = 1
    while i <= length(rows)
        for g in generators
            g isa NC.AbstractAction && !NC.is_antiunitary(g) || continue
            w = last(NC.act(g, rows[i]))
            NC.degree(w) <= d && all(s -> s <= n, NC.sites(w)) ||
                throw(ArgumentError("A declared block action leaves the degree-$d Pauli basis on 1:$n."))
            add(w)
        end
        i += 1
    end
    dense_rows = _workflow_basis_count(n, d)
    retained = length(rows)
    return (family=[map(w -> NC._polynomial(registry, w), sort!(rows))], mode,
        relation=retained == dense_rows ? :exact : :relaxation, regions, row_counts=[retained],
        generation=(; considered, retained, dense_rows, products))
end

"""
    workflow_support(H, order; support_mode=:full, ...)

Build coefficient-one Pauli rows on the original registry for `Relaxation(family=...)`.
`:auto` conservatively selects `:full`. The correlative route completes the powered
interaction graph, then closes the selected family under declared block actions;
the redesign builder folds block orbits and uses their stabilizers. The lattice
route makes one union of rows on periodic intervals (1D) or edges, L's and
plaquettes (2D), plus pairs within periodic graph distance `two_body_range_r`.
These routes retain all three local axes. `regions` lists the initial site stencils
or post-merge cliques, and `row_counts` describes the returned, closed family.

`:graph_distance` selects all axes on each single site and each pair within
`power_graph_radius` in the interaction graph of `H`. `:local_products` selects
all single-site axes and degree-eligible words from products of at most
`local_product_depth` overlapping Hamiltonian factors, ignoring scalar phases.
Neither mode uses a lattice hint or an elimination ordering. Both add factors
to cover `H` and `observable` and close their single row family under declared
unitary actions. `regions` excludes these coverage/symmetry additions.
`generation` counts row candidates offered (including duplicates), unique rows,
the full-basis size (without generating it), and local product attempts.
A proper subset can weaken the dense lower bound. Term splitting is separate.
"""
function workflow_support(H, order::Integer; support_mode=:full, power_graph_radius::Integer=2,
    cs_algo=NC.MF(), min_cluster_size::Integer=0, merge_threshold::Integer=0,
    two_body_range_r::Integer=4, lattice_shape=(H.registry.nsites,), symmetry_generators=(),
    max_rows::Integer=1_000_000, local_product_depth::Integer=2, observable=nothing)
    registry = H.registry
    registry.algebra isa NC.Pauli || throw(ArgumentError("Expected a Pauli polynomial."))
    order >= 0 && power_graph_radius >= 1 && max_rows >= 1 &&
        min_cluster_size >= 0 && merge_threshold >= 0 && two_body_range_r >= 0 && local_product_depth >= 1 ||
        throw(ArgumentError("Invalid order, range, radius, product depth, clique size or row budget."))
    isnothing(observable) || observable.registry === registry ||
        throw(ArgumentError("Observable and Hamiltonian must use the same registry."))
    n, mode = registry.nsites, support_mode === :auto ? :full : support_mode
    d = min(order, n)
    W = NC.word_type(registry, d)
    if mode in (:graph_distance, :local_products)
        return _workflow_physics_support(H, observable, W, mode, d, power_graph_radius,
            local_product_depth, symmetry_generators, max_rows)
    end
    if mode === :full
        regions = [collect(1:n)]
    elseif mode === :correlative_cliques
        graph, _ = NC._interaction_graph(NC.Problem(H))
        alg = cs_algo === :MF ? NC.MF() : cs_algo === :NoElimination ? NC.NoElimination() : cs_algo
        cliques = NC.clique_decomp(_workflow_power(graph, power_graph_radius), alg)
        regions = _workflow_merge_cliques(isempty(cliques) ? [Int[]] : cliques;
            min_size=min_cluster_size, max_union=merge_threshold)
    elseif mode === :contiguous_lattice
        shape = Tuple(lattice_shape)
        length(shape) in (1, 2) && all(x -> x isa Integer && x > 0, shape) &&
            prod(big.(shape)) == n || throw(ArgumentError("Expected a positive periodic 1D or 2D lattice shape covering $n sites."))
        regions = _workflow_lattice_regions(shape, d, two_body_range_r)
    else
        throw(ArgumentError("Unknown Pauli support mode $mode."))
    end
    if mode === :contiguous_lattice
        rows = Set{W}((one(W),))
        sum(big(3)^length(region) for region in regions if !isempty(region); init=big(1)) <= max_rows ||
            throw(ArgumentError("Pauli basis exceeds max_rows=$max_rows."))
        for region in regions
            # Generate words with exactly this occupied-site set. In particular,
            # a punctured interval is not a contiguous word unless selected separately.
            for word in NC.basis_rows(W, registry.algebra, region, d)
                NC.degree(word) == length(region) && push!(rows, word)
            end
        end
        family = [sort!(collect(rows))]
    else
        sum(_workflow_basis_count(length(c), d) for c in (d == 0 ? [Int[]] : regions)) <= max_rows ||
            throw(ArgumentError("Pauli basis exceeds max_rows=$max_rows."))
        family = unique([NC.basis_rows(W, registry.algebra, region, d) for region in regions])
    end
    seen = Set(family)
    total = sum(length, family)
    i = 1
    while i <= length(family)
        for g in symmetry_generators
            g isa NC.AbstractAction && !NC.is_antiunitary(g) || continue
            image = sort!([NC.act(g, w)[2] for w in family[i]])
            image in seen && continue
            mode === :contiguous_lattice && throw(ArgumentError("A declared block action does not preserve the lattice row set."))
            all(w -> NC.degree(w) <= d && all(s -> s <= n, NC.sites(w)), image) ||
                throw(ArgumentError("A declared block action leaves the degree-$d Pauli basis on 1:$n."))
            total += length(image)
            total <= max_rows || throw(ArgumentError("Symmetry closure exceeds max_rows=$max_rows."))
            push!(seen, image)
            push!(family, image)
        end
        i += 1
    end
    if !isnothing(observable)
        sets = mode === :full ? () : Set.(family)
        for (c, w) in NC.terms(observable)
            iszero(c) && continue
            (mode === :full ? NC.degree(w) <= 2d :
                any(any(v -> last(NC.mul(v, w)) in block, block) for block in sets)) ||
                throw(ArgumentError("Observable moment $w is not covered by $mode; increase order or change supports."))
        end
    end
    return (family=[map(w -> NC._polynomial(registry, w), rows) for rows in family], mode,
        relation=mode === :full ? :exact : :relaxation, regions, row_counts=length.(family), generation=nothing)
end

_workflow_entry(b, i, j) = b.first_form + NC.tri_offset(b.side, i, j)
_workflow_psd(b) = b.kind in (NC.PSD, NC.HPSD)

# Compute support after exact cancellation of repeated IDs, with no tolerance.
# Most redesign forms have no repeats; use exact arithmetic only on repeats.
function _workflow_form_support(forms, f)
    ids = Int32[]
    repeated = Set{Int32}()
    for t in NC.term_range(forms, f)
        iszero(forms.coefs[t]) && continue
        id = forms.ids[t]
        id in ids ? push!(repeated, id) : push!(ids, id)
    end
    for id in repeated
        total = sum(Complex{Rational{BigInt}}(forms.coefs[t]) for t in NC.term_range(forms, f) if forms.ids[t] == id)
        iszero(total) && filter!(!=(id), ids)
    end
    return ids
end

function _workflow_pattern(block, support; active=nothing)
    graph = Graphs.SimpleGraph(block.side)
    for j in 1:block.side, i in j+1:block.side
        ids = support[_workflow_entry(block, i, j)]
        edge = isnothing(active) ? !isempty(ids) : any(in(active), ids)
        edge && Graphs.add_edge!(graph, i, j)
    end
    return graph
end

function _workflow_copy_form!(out, program, support, f; scale=1)
    for t in NC.term_range(program.forms, f)
        id = program.forms.ids[t]
        id in support[f] && NC.push_term!(out.forms, scale * program.forms.coefs[t], id)
    end
end

function _workflow_clique_parts!(out, program, support, block, cliques)
    for rows in cliques
        for j in eachindex(rows), i in j:length(rows)
            f = _workflow_entry(block, rows[i], rows[j])
            _workflow_copy_form!(out, program, support, f)
            NC.close_form!(out.forms)
        end
        NC.push_block!(out, block.kind, length(rows), block.label)
    end
end

# Enforce A = sum(E_C' Z_C E_C), Z_C >= 0. Matching *all* clique entries,
# including structural zeros at fill edges, is essential: principal-clique PSD
# alone would test for a PSD completion, not PSD of the supplied matrix.
function _workflow_clique_sum!(out, program, support, block, cliques)
    entries = Dict{Int,Vector{Tuple{Int,ComplexF64}}}()
    for rows in cliques
        for j in eachindex(rows), i in j:length(rows)
            f = _workflow_entry(block, rows[i], rows[j])
            terms = get!(() -> Tuple{Int,ComplexF64}[], entries, f)
            out.nvars += 1
            push!(terms, (out.nvars, 1))
            NC.push_term!(out.forms, 1, out.nvars)
            if block.kind == NC.HPSD && i != j
                out.nvars += 1
                push!(terms, (out.nvars, im))
                NC.push_term!(out.forms, im, out.nvars)
            end
            NC.close_form!(out.forms)
        end
        NC.push_block!(out, block.kind, length(rows), block.label)
    end
    fields = sort!(collect(keys(entries)))
    for f in fields
        for (id, c) in entries[f]
            NC.push_term!(out.forms, isreal(c) ? real(c) : c, id)
        end
        _workflow_copy_form!(out, program, support, f; scale=-1)
        NC.close_form!(out.forms)
    end
    NC.push_block!(out, NC.ZERO, length(fields), block.label)
end

"""
    workflow_term_sparsity(program; ts_algo=:Maximal, ts_step=0, ...)

Transform a symmetry-adapted `NC.MomentProgram` without changing original IDs.
By default only label-zero (ordinary moment) PSD/HPSD blocks are split; all
other cones seed closure and are retained. Use `split_blocks` to narrow that
selection. `protected_ids` pins additional observables in the zero-extension audit.
`:Maximal` has a zero lift outside `report.active_ids`. With
`prune_constant_singletons=true`, it removes only isolated constraints whose
diagonal is a nonnegative constant at normalized identity. `:AggregateDualChordal`
adds exact sparse PSD clique sums, and `:MMD` keeps only principal cliques and is
a relaxation. Each MMD step completes the current graph and extends support
using the clique entries, as in redesign `TermSplit`. `:None` does nothing.
"""
function workflow_term_sparsity(program::NC.MomentProgram{C}; ts_algo=:Maximal, ts_step::Integer=0,
    merge_threshold::Integer=0, split_blocks=findall(b -> b.label == 0 && _workflow_psd(b), program.blocks),
    protected_ids=(), prune_constant_singletons=false) where {C}
    ts_algo in (:None, :Maximal, :AggregateDualChordal, :MMD) ||
        throw(ArgumentError("Unknown term sparsity algorithm $ts_algo."))
    (ts_algo === :MMD ? ts_step >= 1 : ts_step == 0) && merge_threshold >= 0 ||
        throw(ArgumentError("Use ts_step >= 1 for :MMD and ts_step=0 otherwise; merge_threshold must be nonnegative."))
    selected = Set(Int.(collect(split_blocks)))
    all(i -> 1 <= i <= length(program.blocks) && program.blocks[i].label == 0 && _workflow_psd(program.blocks[i]), selected) ||
        throw(ArgumentError("split_blocks must index ordinary (label-zero) PSD or HPSD cones."))
    all(id -> id isa Integer && 1 <= id <= program.nvars, protected_ids) ||
        throw(ArgumentError("protected_ids must be original program IDs."))
    active = BitSet([1])
    union!(active, program.objective_ids, protected_ids)
    NC.validate(program)
    ts_algo === :None && return (program=program, report=(algorithm=ts_algo, relation=:exact, steps=0,
        zero_lift_audit=nothing, active_ids=collect(1:program.nvars), auxiliary_ids=program.nvars+1:program.nvars,
        cliques=Dict{Int,Vector{Vector{Int}}}(), pruned_singletons=Tuple{Int,Int}[]))
    support = [_workflow_form_support(program.forms, f) for f in 1:length(program.forms)]
    cliques = Dict{Int,Vector{Vector{Int}}}()
    steps = 0
    audit = nothing
    pruned_singletons = Tuple{Int,Int}[]
    if ts_algo in (:Maximal, :MMD)
        for (i, block) in enumerate(program.blocks)
            seed = i in selected ? (_workflow_entry(block, j, j) for j in 1:block.side) : NC.form_range(block)
            for f in seed
                union!(active, support[f])
            end
        end
        algo = ts_algo === :Maximal ? NC.MaximalElimination() : NC.MMD()
        while true
            extended = copy(active)
            for i in selected
                block = program.blocks[i]
                rows = sort!.(NC.clique_decomp(_workflow_pattern(block, support; active), algo))
                cliques[i] = rows
                for clique in rows, j in eachindex(clique), k in j:length(clique)
                    union!(extended, support[_workflow_entry(block, clique[k], clique[j])])
                end
            end
            steps += 1
            converged = extended == active
            active = extended
            (ts_algo === :Maximal ? converged : steps == ts_step) && break
        end
        if ts_algo === :Maximal
            # Audit the persisted partition, not just the graph used to build it.
            audit = all(in(active), program.objective_ids) && all(in(active), protected_ids)
            for (i, block) in enumerate(program.blocks)
                if i in selected
                    component = zeros(Int, block.side)
                    for (k, rows) in enumerate(cliques[i])
                        component[rows] .= k
                    end
                    audit &= all(!iszero, component)
                    for j in 1:block.side, k in j:block.side
                        ids = support[_workflow_entry(block, k, j)]
                        audit &= component[k] == component[j] ? all(in(active), ids) : isdisjoint(active, ids)
                    end
                else
                    audit &= all(f -> all(in(active), support[f]), NC.form_range(block))
                end
            end
            audit || error("Term-sparsity zero-extension audit failed.")
        end
    else
        union!(active, 1:program.nvars)
        if ts_algo === :AggregateDualChordal
            for i in selected
                graph = _workflow_pattern(program.blocks[i], support)
                cliques[i] = _workflow_merge_cliques(NC.clique_decomp(graph, NC.MF()); max_union=merge_threshold)
            end
        end
    end

    out = NC.MomentProgram{C}(program.nvars)
    append!(out.objective_ids, program.objective_ids)
    append!(out.objective_coefs, program.objective_coefs)
    for (i, block) in enumerate(program.blocks)
        if ts_algo === :Maximal && prune_constant_singletons && haskey(cliques,i)
            rows = filter(cliques[i]) do clique
                length(clique)==1 || return true
                j = only(clique)
                f = _workflow_entry(block,j,j)
                all(==(1),support[f]) || return true
                constant = sum((Complex{Rational{BigInt}}(program.forms.coefs[t]) for t in
                    NC.term_range(program.forms,f) if program.forms.ids[t]==1);init=0)
                isreal(constant) && real(constant)>=0 || return true
                push!(pruned_singletons,(i,j))
                false
            end
            _workflow_clique_parts!(out,program,support,block,rows)
        elseif !haskey(cliques, i) || length(cliques[i]) == 1
            NC._copy_cone!(out, program, block)
        elseif ts_algo === :AggregateDualChordal
            _workflow_clique_sum!(out, program, support, block, cliques[i])
        else
            _workflow_clique_parts!(out, program, support, block, cliques[i])
        end
    end
    NC.validate(out)
    return (program=out, report=(algorithm=ts_algo, relation=ts_algo === :MMD ? :relaxation : :exact,
        steps, zero_lift_audit=audit, active_ids=collect(active), auxiliary_ids=program.nvars+1:out.nvars,
        cliques, pruned_singletons))
end
