# Sparse exact elimination, with coordinate 1 reserved for the normalized identity.
# A stored pivot is the least nonconstant coordinate in its row. Consequently
# back substitution visits pivots in descending order, without a dense matrix.
function _insert_generator_equation!(pivots,row;constant=true)
    while !isempty(row)
        pivot=minimum((id for id in keys(row) if !constant || id!=1);init=typemax(Int))
        pivot==typemax(Int) && throw(ArgumentError("Continuous equations contradict the normalized identity."))
        coefficient=row[pivot]
        if haskey(pivots,pivot)
            for (id,c) in pivots[pivot]
                value=get(row,id,zero(c))-coefficient*c
                iszero(value) ? delete!(row,id) : (row[id]=value)
            end
        else
            for id in keys(row)
                row[id]/=coefficient
            end
            pivots[pivot]=row
            return false
        end
    end
    return true
end

function _derivative(d,w::W) where {W}
    sink=SectorTerms(W)
    derive!(sink,d,w,one(Complex{Rational{BigInt}}))
    return sink.terms
end

function _check_derivation_polynomial(p,W,derivations,label)
    for (i,d) in enumerate(derivations)
        sink=SectorTerms(W)
        applicable(derive!,sink,d,one(W),one(Complex{Rational{BigInt}})) ||
            throw(ArgumentError("Continuous derivation $i has an incompatible algebra or word width for $label ($p)."))
        for (c,w) in terms(p)
            derive!(sink,d,_pack_word(W,w),Complex{Rational{BigInt}}(c))
        end
        isempty(sink.terms) || throw(ArgumentError(
            "Continuous derivation $i does not annihilate $label ($p) exactly; remove it or supply an invariant polynomial."))
    end
end

function _check_continuous(problem,symmetry,W)
    ds=symmetry.continuous
    isempty(ds) && return
    problem.objective.registry.algebra isa Union{Pauli,Fermionic} || throw(ArgumentError(
        "Continuous reduction currently requires Pauli or fermionic variables."))
    _check_derivation_polynomial(problem.objective,W,ds,"objective")
    for (i,c) in enumerate(problem.constraints)
        _check_derivation_polynomial(c.polynomial,W,ds,"constraint $i")
    end
    letters=filter(!isone,basis_rows(W,problem.objective.registry.algebra,1:problem.objective.registry.nsites,1))
    # A derivation is determined on the algebra's generators. Vectorize its
    # real and imaginary coefficients on those generators, not on moment rows.
    index=Dict{Tuple{Int,W,Bool},Int}()
    function encode(d,g=nothing)
        row=Dict{Int,Rational{BigInt}}()
        for (j,w) in enumerate(letters)
            image=if isnothing(g)
                _derivative(d,w)
            else
                p,u=act(inv(g),w)
                derivative=SectorTerms(W)
                derive!(derivative,d,u,scale(p,one(Complex{Rational{BigInt}})))
                mapped=SectorTerms(W)
                for (v,c) in derivative.terms
                    act!(mapped,g,v,is_antiunitary(g) ? conj(c) : c)
                end
                mapped.terms
            end
            for (u,c) in sort!(collect(image);by=first),imaginary in (false,true)
                value=imaginary ? imag(c) : real(c)
                iszero(value) && continue
                id=get!(index,(j,u,imaginary)) do
                    length(index)+1
                end
                row[id]=value
            end
        end
        row
    end
    pivots=Dict{Int,Dict{Int,Rational{BigInt}}}()
    for d in ds
        _insert_generator_equation!(pivots,encode(d);constant=false)
    end
    for (i,g) in enumerate((symmetry.generators...,symmetry.moment_only...)),(j,d) in enumerate(ds)
        # Use a copy: a failed membership check must not enlarge the declared span.
        _insert_generator_equation!(copy(pivots),encode(d,g);constant=false) || throw(ArgumentError(
            "Finite generator $i does not preserve the continuous span (derivation $j); add the missing derivations or remove the generator."))
    end
end

function _continuous_family(family,derivations;custom=false)
    isempty(derivations) && return family
    blocks=copy(family)
    i=1
    while i<=length(blocks)
        block=blocks[i]
        rows=Set(block.rows)
        queue=copy(block.rows)
        cursor=1
        while cursor<=length(queue)
            for (j,d) in enumerate(derivations),w in keys(_derivative(d,queue[cursor]))
                w in rows && continue
                custom && throw(ArgumentError(
                    "Custom block $i is not closed under continuous derivation $j; supply its complete linear closure or select automatic rows/grouping/splitting."))
                push!(rows,w); push!(queue,w)
                # Merge a block only when a derivative reaches one of its rows,
                # not merely because both blocks already contain the identity.
                for k in reverse(eachindex(blocks))
                    k==i && continue
                    other=blocks[k]
                    other.source==block.source && w in other.rows || continue
                    for v in other.rows
                        if !(v in rows)
                            push!(rows,v); push!(queue,v)
                        end
                    end
                    deleteat!(blocks,k)
                    k<i && (i-=1)
                end
            end
            cursor+=1
        end
        blocks[i]=RowBlock(sort!(collect(rows)),block.polynomial,block.kind,block.source)
        i+=1
    end
    unique!(_block_key,blocks)
    return blocks
end

"An exact reparameterization of orbit moments; new words require rebuilding it."
struct ContinuousMoments{T,W,D}
    base::T
    labels::Vector{MomentLabel{W}}
    coordinates::Vector{Vector{Pair{Int32,Rational{BigInt}}}}
    derivations::D
    covered::BitSet
end
function moment_terms(t::ContinuousMoments,w::AbstractWord)
    pairs=moment_terms(t.base,w)
    if isnothing(pairs)
        # A nonzero eigenvalue forces zero even outside the built word domain.
        for d in t.derivations
            image=_derivative(d,w)
            length(image)==1 && haskey(image,w) && return ()
        end
        return nothing
    end
    return ((id,c*factor) for (old,c) in pairs for (id,factor) in t.coordinates[old])
end
moment_covered(t::ContinuousMoments,id::Integer) = id in t.covered
function refresh_coverage!(t::ContinuousMoments,program::MomentProgram)
    empty!(t.covered)
    union!(t.covered,(1,),program.forms.ids,program.objective_ids)
    return t
end

function _continuous_program(program::MomentProgram{C},table::OrbitMoments{W},derivations;diagonal=true) where {C,W}
    initial=program.nvars
    support=sort!(unique!(reduce(vcat,[collect(sites(label.word)) for label in table.labels];init=Int[])))
    letters=basis_rows(W,W<:PauliWord ? Pauli() : Fermionic(),support,1)
    diagonal &= all(d->all(w->all(==(w),keys(_derivative(d,w))),letters),derivations)
    pivots=Dict{Int,Dict{Int,Rational{BigInt}}}()
    # Close words, not only nonzero orbit coordinates: D(w) can be nonzero
    # when finite invariance forces w itself to zero (a normalizer can rotate
    # the derivations). Such equations must not disappear with that orbit.
    queue=sort!(unique!([collect(keys(table.cache));[label.word for label in table.labels]]))
    seen=Set(queue)
    cursor=1
    labels_seen=length(table.labels)
    while cursor<=length(queue)
        w=queue[cursor]
        cursor+=1
        for d in derivations
            if diagonal
                isempty(_derivative(d,w)) && continue
                m=orbit_coordinates!(table,w)
                for id in (m.re,m.im)
                    id==0 && continue
                    id==1 && throw(ArgumentError("A derivation must annihilate the identity."))
                    pivots[id]=Dict(Int(id)=>one(Rational{BigInt}))
                end
            else
                image=_derivative(d,w)
                sink=OrbitFormSink(table,Dict{Int32,Complex{Rational{BigInt}}}())
                for (v,c) in sort!(collect(image);by=first)
                    if !(v in seen)
                        push!(seen,v); push!(queue,v)
                    end
                    push_term!(sink,c,v)
                end
                equation=sink.terms
                for part in (real,imag)
                    row=Dict(Int(id)=>part(c) for (id,c) in equation if !iszero(part(c)))
                    _insert_generator_equation!(pivots,row)
                end
            end
        end
        while labels_seen<length(table.labels)
            labels_seen+=1
            v=table.labels[labels_seen].word
            if !(v in seen)
                push!(seen,v); push!(queue,v)
            end
        end
    end
    free=[id for id in eachindex(table.labels) if !haskey(pivots,id)]
    directions=[Pair{Int32,Rational{BigInt}}[] for _ in table.labels]
    for (id,old) in enumerate(free)
        push!(directions[old],Int32(id)=>one(Rational{BigInt}))
    end
    for pivot in sort!(collect(keys(pivots));rev=true)
        coefficients=Dict{Int32,Rational{BigInt}}()
        for (old,c) in pivots[pivot]
            old==pivot && continue
            for (id,factor) in directions[old]
                coefficients[id]=get(coefficients,id,zero(c))-c*factor
            end
        end
        directions[pivot]=sort!([id=>c for (id,c) in coefficients if !iszero(c)];by=first)
    end
    out=MomentProgram{C}(length(free))
    function rewrite(ids,coefs)
        result=Dict{Int32,C}()
        for (old,c) in zip(ids,coefs),(id,factor) in directions[old]
            result[id]=get(result,id,zero(C))+c*Float64(factor)
        end
        sort!([id=>c for (id,c) in result if !iszero(c)];by=first)
    end
    for f in 1:length(program.forms)
        range=term_range(program.forms,f)
        for (id,c) in rewrite(view(program.forms.ids,range),view(program.forms.coefs,range))
            push_term!(out.forms,c,id)
        end
        close_form!(out.forms)
    end
    append!(out.blocks,program.blocks)
    for (id,c) in rewrite(program.objective_ids,program.objective_coefs)
        push!(out.objective_ids,id); push!(out.objective_coefs,real(c))
    end
    reduced=ContinuousMoments(table,table.labels[free],directions,derivations,BitSet([1]))
    refresh_coverage!(reduced,out)
    report=(method=diagonal ? :diagonal : :exact_elimination,initial_variables=initial-1,
        closed_variables=length(table.labels)-1,variables=length(free)-1,equations=length(pivots),arithmetic=:complex_rationals)
    return validate(out),reduced,report
end
