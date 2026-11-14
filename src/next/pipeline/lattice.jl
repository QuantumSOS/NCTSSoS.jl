"Induced signed/phase action on the cyclic row representatives at momentum k."
function _momentum_action(orbits, index, h, k; before=nothing)
    images = Vector{Tuple{ComplexF64,Int}}(undef, length(orbits.words))
    for (a, ws) in enumerate(orbits.words)
        p, w = _action_image(h, first(ws), before)
        haskey(index, w) || throw(ArgumentError("The point group does not preserve the row family."))
        b, t = index[w]
        q = p * conj(orbits.phases[b][t])
        images[a] = (scale(q, cispi(2 * k * (t - 1) / orbits.order)), b)
    end
    images
end

function _momentum_projectors(active, actions, projections)
    sectors = [Vector{Tuple{ComplexF64,Int}}[] for _ in projections]
    domain, seen = Set(active), Set{Int}()
    for a in active
        a in seen && continue
        orbit = sort!(unique([action[a][2] for action in actions]))
        all(in(domain), orbit) || error("Little group changed the momentum sector.")
        union!(seen, orbit)
        index = Dict(a => i for (i, a) in enumerate(orbit))
        for (s, coefficients) in enumerate(projections)
            P = zeros(ComplexF64, length(orbit), length(orbit))
            for (action, c) in zip(actions, coefficients), (j, a) in enumerate(orbit)
                phase, b = action[a]
                P[index[b],j] += c * phase
            end
            all(iszero, P) && continue
            F = LinearAlgebra.svd(P)
            for j in findall(>(1e-9), F.S)
                push!(sectors[s], [(F.U[i,j], a) for (i,a) in enumerate(orbit) if abs(F.U[i,j]) > 1e-12])
            end
        end
    end
    sectors
end

function _real_momentum_sectors(sectors, orbits, index, k, anti, point, conjugations)
    candidates=[h for (h,r) in zip(point.elements,conjugations) if mod((r+1)*k,orbits.order)==0]
    _real_character_sectors(sectors,orbits,index,k,anti,candidates)
end

function _real_character_sectors(sectors,orbits,index,k,anti,candidates)
    result = [Vector{Tuple{ComplexF64,Int}}[] for _ in sectors]
    n = length(orbits.words)
    for (s, rows) in enumerate(sectors)
        isempty(rows) && continue
        Q = zeros(ComplexF64, n, length(rows))
        for (j,row) in enumerate(rows), (c,a) in row
            Q[a,j] = c
        end
        found = false
        for h in candidates
            action = _momentum_action(orbits, index, h, k; before=anti)
            image = zeros(ComplexF64, size(Q))
            for a in 1:n, j in axes(Q,2)
                phase,b = action[a]
                image[b,j] += phase * conj(Q[a,j])
            end
            fixed = _fixed_real_columns(Q,image)
            isnothing(fixed) && continue
            result[s] = [[(c,a) for (a,c) in enumerate(v) if abs(c)>1e-12] for v in fixed]
            found = true
            break
        end
        found || throw(ArgumentError("No involutive antiunitary preserves this reduced sector; a real block without doubling is not available from the supplied point group."))
    end
    result
end

function _orbit_pair_kernel!(cache, table, orbits, a, b)
    get!(cache,(a,b)) do
        forms = LinearForms{ComplexF64}()
        sink = OrbitFormSink(table,Dict{Int32,ComplexF64}())
        wa,wb = orbits.words[a],orbits.words[b]
        p,sa = star(first(wa))
        for t in eachindex(wb)
            empty!(sink.terms)
            mul!(sink,sa,wb[t],scale(p*orbits.phases[b][t],ComplexF64(sqrt(length(wa)/length(wb)))))
            _close_orbit_form!(forms,sink.terms)
        end
        forms
    end
end

function _fourier_kernel!(fourier, cache, table, orbits, a, b, characters)
    get!(fourier,(a,b)) do
        raw = _orbit_pair_kernel!(cache,table,orbits,a,b)
        terms = Dict{Int32,ComplexF64}()
        for t in 1:length(raw), j in term_range(raw,t)
            id = raw.ids[j]
            terms[id] = get(terms,id,0im) + characters[t]*raw.coefs[j]
        end
        filter!(p->abs(last(p))>1e-12,terms)
        terms
    end
end

function _assemble_momentum!(program, sink, orbits, k, sectors, sector_labels, cache; real_blocks=false)
    characters = [cispi(-2k*t/orbits.order) for t in 0:(orbits.order-1)]
    _assemble_character!(program,sink,orbits,characters,k,sectors,sector_labels,cache;real_blocks)
end

function _assemble_character!(program,sink,orbits,characters,k,sectors,sector_labels,cache;real_blocks=false)
    fourier = Dict{Tuple{Int,Int},Dict{Int32,ComplexF64}}()
    for (s, reduced_rows) in enumerate(sectors)
        isempty(reduced_rows) && continue
        for j in eachindex(reduced_rows), i in j:length(reduced_rows)
            empty!(sink.terms)
            for (ca,a) in reduced_rows[i], (cb,b) in reduced_rows[j]
                kernel = _fourier_kernel!(fourier,cache,sink.table,orbits,max(a,b),min(a,b),characters)
                c = conj(ca)*cb
                for (id,coefficient) in kernel
                    sink.terms[id] = get(sink.terms,id,0im) + c*(a < b ? conj(coefficient) : coefficient)
                end
            end
            _close_orbit_form!(program.forms, sink.terms; diagonal=real_blocks || i == j)
        end
        push_block!(program, real_blocks ? PSD : HPSD, length(reduced_rows), length(sector_labels) + 1)
        push!(sector_labels, (k,s))
    end
end

# The selected elements already form a subgroup. Keep their order, but give
# character calculations a generating set instead of every subgroup element.
function _little_group(elements)
    generators = [first(elements)]
    generated = OrbitGroup(generators)
    for g in elements
        haskey(generated.indices,g) && continue
        push!(generators,g)
        generated = OrbitGroup(generators)
    end
    indices = Dict(g=>i for (i,g) in enumerate(elements))
    OrbitGroup(elements,indices,[indices[g] for g in generators])
end

"""
    lattice_program(rows, objective, translation, order, point_generators; moment_generators=())

One-dimensional symmorphic lattice spike. Name the translation and its order;
the finite point group must conjugate it to itself or its inverse. Fourier
blocks are built from cyclic representatives. Shared little-group projectors
are applied to one orbit of momentum rows at a time, not to all site rows.
The full lattice group is never enumerated.
"""
function lattice_program(rows::AbstractVector{W}, objective, translation, order,
    point_generators::AbstractVector{<:AbstractAction}; moment_generators=()) where {W<:AbstractWord}
    translation=_group_action(W,translation)
    point_generators=map(g->_group_action(W,g),point_generators)
    moment_generators=map(g->_group_action(W,g),moment_generators)
    all_generators = (translation, point_generators..., moment_generators...)
    _check_objective(objective, all_generators)
    orbits = CyclicRows(rows, translation, order)
    point = OrbitGroup(point_generators)
    conjugations = map(point.elements) do h
        t = _conjugate_action(h,translation)
        t == translation && return 1
        t == inv(translation) && return -1
        throw(ArgumentError("Name a normal translation subgroup: each point map must conjugate the translation to itself or its inverse."))
    end
    index = Dict(w => (a,t) for (a,ws) in enumerate(orbits.words) for (t,w) in enumerate(ws))
    # Diagonal elements expose zero character sectors before moment-orbit
    # enumeration. This avoids storing all translations of forced zeros.
    diagonal = filter([point.elements; collect(moment_generators)]) do h
        !isone(h) && all(w -> act(h,w)[2] == w, rows)
    end
    table = OrbitMoments(W, all_generators; annihilators=diagonal)
    program = MomentProgram{ComplexF64}()
    sink = OrbitFormSink(table, Dict{Int32,ComplexF64}())
    # The little group changes only at k=0 and, for even order, k=order/2.
    cache = Dict{Bool,Tuple{typeof(point),Vector{Vector{ComplexF64}}}}()
    kernels = Dict{Tuple{Int,Int},LinearForms{ComplexF64}}()
    sector_labels = Tuple{Int,Int}[]
    little_orders = Int[]
    anti_index = findfirst(is_antiunitary, moment_generators)
    if !isnothing(anti_index)
        anti = moment_generators[anti_index]
        _conjugate_action(anti,translation) == translation ||
            throw(ArgumentError("The M1 lattice antiunitary must commute with the named translation."))
    end
    for k in 0:(order - 1)
        paired = any(==(-1), conjugations)
        paired && k > mod(-k, order) && continue
        full = mod(2k, order) == 0 || !paired
        little, projections = get!(cache, full) do
            elements = point.elements[findall(i -> full || conjugations[i] == 1, eachindex(point.elements))]
            group = full ? point : _little_group(elements)
            ps, _ = shared_projections(group)
            (group, ps)
        end
        push!(little_orders, length(little))
        active = [a for a in eachindex(orbits.words) if
            mod(4k * length(orbits.words[a]) - Int(orbits.stabilizers[a].k) * order, 4order) == 0]
        actions = [_momentum_action(orbits, index, h, k) for h in little.elements]
        sectors = _momentum_projectors(active, actions, projections)
        if !isnothing(anti_index)
            sectors = _real_momentum_sectors(sectors, orbits, index, k,
                moment_generators[anti_index], point, conjugations)
        end
        _assemble_momentum!(program, sink, orbits, k, sectors, sector_labels, kernels; real_blocks=!isnothing(anti_index))
    end
    p, table = _orbit_objective!(program, table, objective)
    return _anti_real_program(p, all_generators), table, (; sectors=sector_labels, little_orders, point_order=length(point))
end
