# A discarded free coordinate retains an id outside the solved program. The
# coverage check rejects it, after polynomial contributions have cancelled.
struct SectorMoments{T,W,R}
    base::T
    labels::Vector{MomentLabel{W}}
    coordinates::Vector{Vector{Pair{Int32,R}}}
    covered::BitSet
end
function moment_terms(t::SectorMoments,w::AbstractWord)
    pairs=moment_terms(t.base,w)
    isnothing(pairs) && return nothing
    ((id,c*factor) for (old,c) in pairs for (id,factor) in t.coordinates[old])
end
moment_covered(t::SectorMoments,id::Integer) = id in t.covered
function refresh_coverage!(t::SectorMoments,program::MomentProgram)
    empty!(t.covered)
    union!(t.covered,(1,),program.forms.ids,program.objective_ids)
    t
end

function _sector_add!(row,id,c)
    value=get(row,id,zero(c))+c
    iszero(value) ? delete!(row,id) : (row[id]=value)
    nothing
end

function _sector_directions(equations,table,::Type{R}) where R
    labels=table.labels
    order=[1;sort!(collect(2:length(labels));by=id->(-degree(labels[id].word),id))]
    inverse=invperm(order)
    pivots=Dict{Int,Dict{Int,R}}()
    for equation in equations
        row=Dict{Int,R}()
        for (old,c) in equation
            if table isa ContinuousMoments
                for (id,factor) in table.coordinates[old]
                    _sector_add!(row,inverse[id],R(c)*R(factor))
                end
            else
                _sector_add!(row,inverse[old],R(c))
            end
        end
        _insert_generator_equation!(pivots,row;label="Sector")
    end
    free=[id for id in eachindex(labels) if !haskey(pivots,id)]
    directions=[Pair{Int32,R}[] for _ in labels]
    for (id,old) in enumerate(free)
        push!(directions[old],Int32(id)=>one(R))
    end
    coefficients=Dict{Int32,R}()
    for pivot in sort!(collect(keys(pivots));rev=true)
        empty!(coefficients)
        for (old,c) in pivots[pivot]
            old==pivot && continue
            for (id,factor) in directions[old]
                _sector_add!(coefficients,id,-c*factor)
            end
        end
        directions[pivot]=sort!(collect(coefficients);by=first)
    end
    (coordinates=directions[inverse],free=order[free],rank=length(pivots))
end

function _sector_rewrite!(row,ids,coefs,directions)
    empty!(row)
    C=valtype(typeof(row))
    for (old,c) in zip(ids,coefs)
        coefficient=C(c)
        for (id,factor) in directions[old]
            _sector_add!(row,id,coefficient*factor)
        end
    end
    row
end

function _sector_rewrite_exact!(row,wide,ids,coefs,directions)
    try
        return _sector_rewrite!(row,ids,coefs,directions)
    catch error
        error isa Union{OverflowError,InexactError} || rethrow()
        # A tiny spectral coefficient can require a large binary denominator.
        # Retry this form, not every other form and the exact elimination.
        return _sector_rewrite!(wide,ids,coefs,directions)
    end
end

function _eliminate_sectors(program::MomentProgram,table,equations,sector_blocks)
    isnothing(equations) && return program,table,nothing
    try
        return _eliminate_sectors(program,table,equations,sector_blocks,Rational{Int})
    catch error
        error isa Union{OverflowError,InexactError} || rethrow()
        return _eliminate_sectors(program,table,equations,sector_blocks,Rational{BigInt})
    end
end

function _eliminate_sectors(program::MomentProgram{C},table,equations,sector_blocks,::Type{R}) where {C,R}
    reduction=_sector_directions(equations,table,R)
    directions=reduction.coordinates
    out=MomentProgram{C}(length(reduction.free))
    row=Dict{Int32,Complex{R}}()
    wide=Dict{Int32,Complex{Rational{BigInt}}}()
    arithmetic=R
    before=0
    for (i,block) in enumerate(program.blocks)
        i in sector_blocks && continue
        for f in form_range(block)
            range=term_range(program.forms,f)
            before+=length(range)
            exact=_sector_rewrite_exact!(row,wide,view(program.forms.ids,range),view(program.forms.coefs,range),directions)
            exact===wide && (arithmetic=Rational{BigInt})
            for (id,c) in sort!(collect(exact);by=first)
                # Exact arithmetic ends at emission into the conic contract.
                push_term!(out.forms,C(c),id)
            end
            close_form!(out.forms)
        end
        push_block!(out,block.kind,block.side,block.label)
    end
    after=length(out.forms.ids)
    report=(accepted=after<=before,reason=after<=before ? :accepted : :coefficient_growth,
        rank=reduction.rank,coefficients_before=before,
        coefficients_after=after,arithmetic=arithmetic,initial_variables=program.nvars-1)
    after<=before || return program,table,merge(report,(variables=program.nvars-1,))
    objective=_sector_rewrite_exact!(row,wide,program.objective_ids,program.objective_coefs,directions)
    objective===wide && (report=merge(report,(arithmetic=Rational{BigInt},)))
    used=BitSet([1]); union!(used,out.forms.ids,keys(objective))
    active=collect(used)
    # Live ids lead the numbering. Uncovered ids remain distinct so polynomial
    # combinations can cancel them before expectation checks coverage.
    all_free=[active;[id for id in eachindex(reduction.free) if !(id in used)]]
    renumber=invperm(all_free)
    out.nvars=length(active)
    for i in eachindex(out.forms.ids)
        out.forms.ids[i]=Int32(renumber[out.forms.ids[i]])
    end
    for (id,c) in sort!(collect(objective);by=first)
        push!(out.objective_ids,Int32(renumber[id]))
        push!(out.objective_coefs,Float64(real(c)))
    end
    for pairs in directions
        for i in eachindex(pairs)
            id,c=pairs[i]; pairs[i]=Int32(renumber[id])=>c
        end
        sort!(pairs;by=first)
    end
    reduced=SectorMoments(table,table.labels[reduction.free[active]],directions,BitSet([1]))
    out,reduced,merge(report,(variables=out.nvars-1,))
end
