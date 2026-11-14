"Block generators close the family; moment-only generators add moment relations only."
struct FamilySymmetry{G,M}
    generators::G
    moment_only::M
    function FamilySymmetry(generators=AbstractAction[];moment_only=())
        gs=collect(generators); ms=Tuple(moment_only)
        all(g->g isa AbstractAction,gs) && all(g->g isa AbstractAction,ms) ||
            throw(ArgumentError("FamilySymmetry needs multiplicative actions; continuous derivations are a separate reduction."))
        any(is_antiunitary,gs) && throw(ArgumentError("Put antiunitary actions in moment_only."))
        new{typeof(gs),typeof(ms)}(gs,ms)
    end
end

"A row rule per site group, a split algorithm, or an unchanged custom family."
struct BlockRelaxation{R,G,S,F,Y,K} <: AbstractRelaxation
    order::Union{Int,Nothing}
    row_rule::R
    grouping::G
    split::S
    family::F
    symmetry::Y
    strengthen::K
    function BlockRelaxation(;order=nothing,row_rule=nothing,grouping=nothing,split=nothing,
        family=nothing,symmetry=FamilySymmetry(),strengthen=())
        xor(isnothing(order),isnothing(family)) || throw(ArgumentError("Supply exactly one of order or a custom family."))
        isnothing(order) || order isa Integer && order>=0 || throw(ArgumentError("The row order must be a nonnegative integer."))
        symmetry isa FamilySymmetry || throw(ArgumentError("Use FamilySymmetry to name block generators and moment-only maps."))
        strengthen=Tuple(strengthen)
        all(k->k isa StrengtheningKind,strengthen) || throw(ArgumentError("Use strengthening descriptions such as DensityMatrix(sites) or InvariantCone(sites, r)."))
        if !isnothing(family)
            all(isnothing,(row_rule,grouping,split)) || throw(ArgumentError("A custom family is used as supplied; remove the row rule, grouping and split options."))
            family=collect.(family)
            !isempty(family) && all(!isempty,family) || throw(ArgumentError("A custom family and each of its blocks must be nonempty."))
        end
        new{typeof(row_rule),typeof(grouping),typeof(split),typeof(family),typeof(symmetry),typeof(strengthen)}(
            isnothing(order) ? nothing : Int(order),row_rule,grouping,split,family,symmetry,strengthen)
    end
end

"Sizes and closure are available through both inspect and SolveResult.sizes."
struct FamilySizes{R}
    program::ProgramSizes
    closure::R
end
function Base.show(io::IO,s::FamilySizes)
    show(io,s.program)
    r=s.closure
    print(io,"; family ",r.initial_blocks," → ",r.closed_blocks," blocks, ",length(r.block_orbits),
        " orbits; block generators: ",[nameof(typeof(g)) for g in r.block_generators],
        "; moment-only: ",[nameof(typeof(g)) for g in r.moment_only])
    if !isnothing(r.strengthening)
        print(io,"; strengthening: ",r.strengthening.added_blocks," blocks, ",r.strengthening.added_equalities," equalities; ")
        show(io,r.valid_for)
    end
end
Base.show(io::IO,::MIME"text/plain",s::FamilySizes) = show(io,s)

function _custom_family_rows(problem,description,W)
    registry=problem.objective.registry
    map(description.family) do block
        all(p->p isa Polynomial && p.registry===registry,block) ||
            throw(ArgumentError("Every custom row must be a polynomial from the problem registry."))
        all(p->length(p.data)==1 && p.data.coefficients[1]==1,block) ||
            throw(ArgumentError("Custom rows must be single normal words with coefficient 1."))
        rows=W[_pack_word(W,p.data.words[1]) for p in block]
        allunique(rows) || throw(ArgumentError("A custom block has duplicate rows; remove them explicitly."))
        sort!(rows)
    end
end

function build(problem::Problem,description::BlockRelaxation)
    Problem(problem.objective;constraints=problem.constraints)
    _check_hamiltonian_strengthening(problem,description.strengthen)
    custom=!isnothing(description.family)
    if custom
        all(p->p isa Polynomial,Iterators.flatten(description.family)) ||
            throw(ArgumentError("Custom family rows must be polynomials from the problem registry."))
    end
    order=custom ? maximum(degree,Iterators.flatten(description.family)) : description.order
    max_letters=max(2order,degree(problem.objective),maximum((degree(c.polynomial) for c in problem.constraints);init=0))
    max_letters=max(max_letters,maximum((_strengthening_width(k,problem,order) for k in description.strengthen);init=0))
    registry=problem.objective.registry
    W=word_type(registry,max_letters)
    objective=[(ComplexF64(c),_pack_word(W,w)) for (c,w) in terms(problem.objective)]
    identity=[(1.0+0im,one(W))]
    if custom
        rows=_custom_family_rows(problem,description,W)
        groups=[sort!(unique!(reduce(vcat,collect.(sites.(block));init=Int[]))) for block in rows]
    else
        groups=_site_groups(problem,description.grouping)
        rows=[_group_rows(description.row_rule,W,registry,group,order) for group in groups]
    end
    family=[RowBlock(block,identity,HPSD,0) for block in rows]
    custom && length(unique(_block_key,family))!=length(family) &&
        throw(ArgumentError("The custom family has duplicate blocks; remove duplicates explicitly."))
    any(b->one(W) in b.rows,family) || throw(ArgumentError("The family must contain the identity row in at least one block."))
    local_sizes=Int[]
    for (i,constraint) in enumerate(problem.constraints)
        local_order=order-cld(degree(constraint.polynomial),2)
        local_order>=0 || throw(ArgumentError("Constraint $i exceeds the row order; increase order to at least $(cld(degree(constraint.polynomial),2))."))
        polynomial=[(ComplexF64(c),_pack_word(W,w)) for (c,w) in terms(constraint.polynomial)]
        support=sort!(unique!(reduce(vcat,[collect(sites(w)) for (_,w) in polynomial];init=Int[])))
        assigned=findall(group->issubset(support,group),groups)
        isempty(assigned) && throw(ArgumentError("No site group covers constraint $i; increase the grouping level or supply a group containing its complete support."))
        for j in assigned
            local_rows=_group_rows(description.row_rule,W,registry,groups[j],local_order)
            push!(family,RowBlock(local_rows,polynomial,constraint isa Zero ? ZERO : HPSD,i))
            push!(local_sizes,length(local_rows))
        end
    end
    family=unique(_block_key,custom ? family : _split_family(family,objective,description.split))
    letters=[filter(w->!isone(w),basis_rows(W,registry.algebra,[site],1)) for site in 1:registry.nsites]
    program,table,closure=family_program(family,objective;generators=description.symmetry.generators,
        moment_generators=description.symmetry.moment_only,custom,site_letters=letters,real_form=_real_form(problem))
    extra=_strengthen!(program,table,problem,description.strengthen,closure)
    valid_for=isempty(description.strengthen) ? ValidStates(problem.constraints) :
        StrengthenedStates(ValidStates(problem.constraints),description.strengthen)
    closure=merge(closure,(strengthening=extra,valid_for=valid_for))
    sizes=ProgramSizes(sum(length(b.rows) for b in closure.family),length(table.representatives),program.nvars-1,
        length(program.blocks),maximum((b.side for b in program.blocks if b.kind in (PSD,HPSD));init=0),local_sizes)
    return (program=program,table=table,sizes=FamilySizes(sizes,closure),valid_for=valid_for)
end
