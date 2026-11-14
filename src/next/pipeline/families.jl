"One principal moment/localizing block; source 0 is the moment matrix."
struct RowBlock{W,C}
    rows::Vector{W}
    polynomial::Vector{Tuple{C,W}}
    kind::ConeKind
    source::Int
end

_block_key(b::RowBlock) = (b.source,Tuple(b.rows))
function _mapped_block(b::RowBlock,g)
    rows=sort!([act(g,w)[2] for w in b.rows])
    RowBlock(rows,b.polynomial,b.kind,b.source)
end

"Close algorithm families, or check custom families, using generators only."
function close_family(family,generators;custom=false)
    isempty(family) && throw(ArgumentError("The block family is empty; supply at least one nonempty block."))
    blocks=copy(family)
    all(b->!isempty(b.rows) && issorted(b.rows) && allunique(b.rows),blocks) ||
        throw(ArgumentError("Block rows must be nonempty, sorted and distinct."))
    index=Dict(_block_key(b)=>i for (i,b) in enumerate(blocks))
    length(index)==length(blocks) || throw(ArgumentError("The custom family has duplicate blocks; remove duplicates explicitly."))
    i=1
    while i<=length(blocks)
        for g in generators
            image=_mapped_block(blocks[i],g)
            key=_block_key(image)
            haskey(index,key) && continue
            custom && throw(ArgumentError("The custom block family is not closed under a block generator: block $i has a missing image. Supply the complete family or select an automatic grouping/split algorithm."))
            push!(blocks,image); index[key]=length(blocks)
        end
        i+=1
    end
    return blocks,index
end

# Schreier generators preserve the representative block. This lists block
# images, not all elements of the full lattice group.
function _family_orbits(blocks,index,generators)
    isempty(generators) && return [(i,[i],AbstractAction[],nothing) for i in eachindex(blocks)]
    length(blocks)==1 && return [(1,[1],generators,nothing)]
    identity=one(first(generators))
    visited=Set{Int}()
    result=Tuple{Int,Vector{Int},typeof(generators),Dict{Int,typeof(identity)}}[]
    for root in eachindex(blocks)
        root in visited && continue
        orbit=Int[root]; transport=Dict(root=>identity)
        stabilizers=eltype(generators)[]
        seen=Set{eltype(generators)}()
        cursor=1
        while cursor<=length(orbit)
            i=orbit[cursor]
            for g in generators
                j=index[_block_key(_mapped_block(blocks[i],g))]
                edge=compose(g,transport[i])
                if !haskey(transport,j)
                    transport[j]=edge; push!(orbit,j)
                else
                    h=compose(inv(transport[j]),edge)
                    if !isone(h) && !(h in seen)
                        push!(seen,h); push!(stabilizers,h)
                    end
                end
            end
            cursor+=1
        end
        union!(visited,orbit)
        push!(result,(root,orbit,stabilizers,transport))
    end
    result
end

struct BlockProductSink{S,W}
    sink::S
    right::W
end
push_term!(s::BlockProductSink,c::Number,w) = mul!(s.sink,w,s.right,c)

function _block_product!(sink,a,b,polynomial,coefficient=1)
    phase,left=star(a)
    right=BlockProductSink(sink,b)
    for (c,w) in polynomial
        mul!(right,left,w,scale(phase,coefficient*c))
    end
end

function _append_family_sectors!(program,table,block,sectors;real_blocks=false)
    sink=OrbitFormSink(table,Dict{Int32,ComplexF64}())
    for rows in sectors
        isempty(rows) && continue
        for j in eachindex(rows),i in j:length(rows)
            empty!(sink.terms)
            for (ca,a) in rows[i],(cb,b) in rows[j]
                _block_product!(sink,a,b,block.polynomial,conj(ca)*cb)
            end
            _close_orbit_form!(program.forms,sink.terms;diagonal=real_blocks || i==j)
        end
        n=length(rows)
        kind=block.kind==ZERO ? ZERO : real_blocks ? PSD : HPSD
        push_block!(program,kind,kind==ZERO ? n*(n+1)÷2 : n,block.source)
    end
end

function _whole_lattice_family(block,objective,generators,moment_generators,site_letters)
    (isnothing(site_letters) || length(site_letters)<2 || isempty(generators) || block.source!=0) && return nothing
    structure=discover_lattice(generators,site_letters;required=false)
    isnothing(structure) && return nothing
    translations=structure.translations
    point=structure.point_generators
    moments=Tuple(unique([collect(moment_generators);structure.moment_generators]))
    if length(translations.orders)==1 && all(point) do g
        conjugate=_conjugate_action(g,only(translations.generators))
        conjugate in (only(translations.generators),inv(only(translations.generators)))
    end
        return lattice_program(block.rows,objective,only(translations.generators),only(translations.orders),point;
            moment_generators=moments)
    end
    abelian_lattice_program(block.rows,objective,translations,point;moment_generators=moments)
end

"Build one cone per block orbit, sharing orbit moments across all blocks."
function family_program(family::AbstractVector{<:RowBlock{W}},objective;
    generators=AbstractAction[],moment_generators=(),custom=false,site_letters=nothing) where {W<:AbstractWord}
    generators=map(g->_group_action(W,g),generators)
    moment_generators=map(g->_group_action(W,g),moment_generators)
    any(is_antiunitary,generators) && throw(ArgumentError("Put antiunitary maps in moment_only; block generators must be unitary."))
    all_generators=(generators...,moment_generators...)
    _check_objective(objective,all_generators)
    for block in family
        _check_objective(block.polynomial,all_generators;name="constraint $(block.source)")
    end
    blocks,index=close_family(family,generators;custom)
    orbits=_family_orbits(blocks,index,generators)
    report=(initial_blocks=length(family),closed_blocks=length(blocks),
        added_blocks=length(blocks)-length(family),block_orbits=[orbit for (_,orbit,_,_) in orbits],
        representatives=[root for (root,_,_,_) in orbits],family=blocks,
        block_generators=generators,moment_only=moment_generators)
    if length(blocks)==1
        lattice=_whole_lattice_family(only(blocks),objective,generators,moment_generators,site_letters)
        if !isnothing(lattice)
            p,table,route=lattice
            covered=Set(p.forms.ids)
            all(id==1 || id in covered for id in p.objective_ids) ||
                throw(ArgumentError("The block family does not cover the objective; increase the row order."))
            return p,table,merge(report,(route=:lattice,reduction=route))
        end
    end
    table=OrbitMoments(W,all_generators)
    program=MomentProgram{ComplexF64}()
    anti_index=findfirst(is_antiunitary,moment_generators)
    for (root,_,stabilizers,transport) in orbits
        block=blocks[root]
        if isempty(stabilizers)
            sectors=[[[(1.0+0im,w)] for w in block.rows]]
        else
            group=OrbitGroup(stabilizers)
            projections,_=shared_projections(group)
            sectors=_project_orbits(block.rows,stabilizers,group,projections)
        end
        if !isnothing(anti_index)
            anti=moment_generators[anti_index]
            # Bring the antiunitary image back to this representative using
            # the unitary block orbit, before testing its real structure.
            target=get(index,_block_key(_mapped_block(block,anti)),0)
            if target!=root && !isnothing(transport) && haskey(transport,target)
                anti=compose(inv(transport[target]),anti)
            end
            isempty(stabilizers) && (group=OrbitGroup([one(anti)]))
            sectors=_real_finite_sectors(block.rows,sectors,group,anti)
        end
        _append_family_sectors!(program,table,block,sectors;real_blocks=!isnothing(anti_index))
    end
    # A representative found while following an orbit is not necessarily
    # covered by an entry. Check the actual cone coordinates, not table size.
    covered=Set(program.forms.ids)
    sink=OrbitFormSink(table,Dict{Int32,ComplexF64}())
    for (c,w) in objective; push_term!(sink,c,w); end
    all(id==1 || iszero(c) || id in covered for (id,c) in sink.terms) ||
        throw(ArgumentError("The block family does not cover the objective; increase the row order or supply blocks containing its moments."))
    program,table=_orbit_objective!(program,table,objective)
    return program,table,merge(report,(route=:block_stabilizers,))
end
