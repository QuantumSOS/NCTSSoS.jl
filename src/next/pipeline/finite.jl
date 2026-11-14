import GroupsCore
import SymbolicWedderburn

# This adapter stores O(|G|) group data, not a multiplication table.
struct OrbitGroup{A} <: GroupsCore.Group
    elements::Vector{A}
    indices::Dict{A,Int}
    generators::Vector{Int}
end
struct OrbitGroupElement{A} <: GroupsCore.GroupElement
    group::OrbitGroup{A}
    index::Int
end

function OrbitGroup(generators::AbstractVector{A}; limit=65535) where {A<:AbstractAction}
    isempty(generators) && throw(ArgumentError("Supply at least one finite-group generator."))
    any(is_antiunitary, generators) && throw(ArgumentError("Finite block generators must be unitary; pass antiunitary maps as moment generators."))
    1 <= limit <= 65535 || throw(ArgumentError("The generic group limit must be at most 65,535."))
    elements = A[one(first(generators))]
    indices = Dict{A,Int}(first(elements) => 1)
    i = 1
    while i <= length(elements)
        for g in generators
            h = compose(g, elements[i])
            haskey(indices, h) && continue
            length(elements) < limit || throw(ArgumentError(
                "Generic finite group exceeds the limit $limit (SymbolicWedderburn 0.4.4 uses 16-bit indices; hard limit 65,535). Name translations for the lattice engine or use the abelian engine."))
            push!(elements, h); indices[h] = length(elements)
        end
        i += 1
    end
    OrbitGroup(elements, indices, [indices[g] for g in generators])
end
Base.eltype(::Type{OrbitGroup{A}}) where {A} = OrbitGroupElement{A}
Base.IteratorSize(::Type{<:OrbitGroup}) = Base.HasLength()
Base.length(g::OrbitGroup) = length(g.elements)
Base.iterate(g::OrbitGroup, i=1) = i > length(g) ? nothing : (OrbitGroupElement(g, i), i + 1)
Base.one(g::OrbitGroup) = OrbitGroupElement(g, 1)
Base.parent(g::OrbitGroupElement) = g.group
Base.copy(g::OrbitGroupElement) = g
Base.isone(g::OrbitGroupElement) = g.index == 1
Base.:(==)(a::OrbitGroupElement, b::OrbitGroupElement) = a.group === b.group && a.index == b.index
Base.hash(g::OrbitGroupElement, h::UInt) = hash(g.index, h)
Base.inv(g::OrbitGroupElement) = OrbitGroupElement(g.group, g.group.indices[inv(g.group.elements[g.index])])
function Base.:*(g::OrbitGroupElement, h::OrbitGroupElement)
    g.group === h.group || throw(ArgumentError("Elements belong to different groups."))
    OrbitGroupElement(g.group, g.group.indices[compose(g.group.elements[g.index], h.group.elements[h.index])])
end
GroupsCore.order(::Type{I}, g::OrbitGroup) where {I<:Integer} = I(length(g))
function GroupsCore.order(::Type{I}, g::OrbitGroupElement) where {I<:Integer}
    h, n = g, 1
    while !isone(h)
        h *= g; n += 1
    end
    I(n)
end
GroupsCore.gens(g::OrbitGroup) = [OrbitGroupElement(g, i) for i in g.generators]

"Compute the group-algebra projectors once, shared by all row orbits."
function shared_projections(group::OrbitGroup)
    sw = SymbolicWedderburn
    tbl = sw.CharacterTable(Rational{Int}, group)
    irr = sw.irreducible_characters(tbl)
    real_irr, _ = sw.affordable_real(irr, ones(Int, length(irr)))
    algebra = sw._group_algebra(group)
    projectors, ranks = sw.minimal_projection_system(real_irr, algebra)
    coefficients = [[ComplexF64(p[g]) for g in group] for p in projectors]
    return coefficients, ranks
end

function _row_orbits(rows::AbstractVector{W}, generators) where {W}
    domain, seen = Set(rows), Set{W}()
    length(domain) == length(rows) || throw(ArgumentError("The row list contains duplicates."))
    result = Vector{W}[]
    for r in sort!(collect(domain))
        r in seen && continue
        orbit = W[r]; push!(seen, r)
        i = 1
        while i <= length(orbit)
            for g in generators
                _, w = act(g, orbit[i])
                w in domain || throw(ArgumentError("A block generator does not preserve the rows; use it only for moment closure or supply a closed row family."))
                w in seen && continue
                push!(seen, w); push!(orbit, w)
            end
            i += 1
        end
        push!(result, orbit)
    end
    result
end

# Each reduced row is a short list of coefficients and words. Different row
# orbits are joined only after applying the SAME group-algebra projector.
function _project_orbits(rows::AbstractVector{W}, generators, group, projections) where {W}
    sectors = [Vector{Tuple{ComplexF64,W}}[] for _ in projections]
    for orbit in _row_orbits(rows, generators)
        index = Dict(w => i for (i, w) in enumerate(orbit))
        for (s, coefficients) in enumerate(projections)
            P = zeros(ComplexF64, length(orbit), length(orbit))
            for (g, c) in zip(group.elements, coefficients)
                iszero(c) && continue
                for (j, w) in enumerate(orbit)
                    p, u = act(g, w)
                    P[index[u], j] += scale(p, c)
                end
            end
            F = LinearAlgebra.svd(P)
            for j in findall(>(1e-9), F.S)
                push!(sectors[s], [(F.U[i,j], w) for (i, w) in enumerate(orbit) if abs(F.U[i,j]) > 1e-12])
            end
        end
    end
    sectors
end

function _projected_program(::Type{W}, sectors, objective, generators;real_blocks=false) where {W<:AbstractWord}
    table = OrbitMoments(W, generators)
    program = MomentProgram{ComplexF64}()
    sink = OrbitFormSink(table, Dict{Int32,ComplexF64}())
    for (s, rows) in enumerate(sectors)
        isempty(rows) && continue
        for j in eachindex(rows), i in j:length(rows)
            empty!(sink.terms)
            for (ca, a) in rows[i], (cb, b) in rows[j]
                p, sa = star(a)
                mul!(sink, sa, b, scale(p, conj(ca) * cb))
            end
            _close_orbit_form!(program.forms, sink.terms; diagonal=real_blocks || i == j)
        end
        push_block!(program, real_blocks ? PSD : HPSD, length(rows), s)
    end
    _orbit_objective!(program, table, objective)
end

"Direct finite-group reduction, using shared projectors on one row orbit at a time."
function finite_program(rows::AbstractVector{W}, objective, generators::AbstractVector{<:AbstractAction};
    moment_generators=(), limit=65535) where {W<:AbstractWord}
    generators=map(g->_group_action(W,g),generators)
    all_generators = (generators..., moment_generators...)
    _check_objective(objective, all_generators)
    group = OrbitGroup(generators; limit)
    projections, ranks = shared_projections(group)
    sectors = _project_orbits(rows, generators, group, projections)
    anti_index=findfirst(is_antiunitary,moment_generators)
    if !isnothing(anti_index)
        sectors=_real_finite_sectors(rows,sectors,group,moment_generators[anti_index])
    end
    p, table = _projected_program(W, sectors, objective, all_generators;real_blocks=!isnothing(anti_index))
    return _anti_real_program(p, all_generators), table, (; group_order=length(group), projection_ranks=ranks, sectors)
end
