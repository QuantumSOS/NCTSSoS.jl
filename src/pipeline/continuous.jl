# Sparse exact elimination, with coordinate 1 reserved for the normalized identity.
# A stored pivot is the least nonconstant coordinate in its row. Consequently
# back substitution visits pivots in descending order, without a dense matrix.
function _insert_generator_equation!(pivots,row;constant=true,label="Continuous")
    while !isempty(row)
        pivot=minimum((id for id in keys(row) if !constant || id!=1);init=typemax(Int))
        pivot==typemax(Int) && throw(ArgumentError("$label equations contradict the normalized identity."))
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
    problem.objective.registry.algebra isa Union{Pauli,Fermionic,Bosonic} || throw(ArgumentError(
        "Continuous reduction currently requires Pauli, fermionic or bosonic variables."))
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
    for (name,generators) in (("Finite",symmetry.generators),("moment-only",symmetry.moment_only)),
        (i,g) in enumerate(generators),(j,d) in enumerate(ds)
        # Use a copy: a failed membership check must not enlarge the declared span.
        _insert_generator_equation!(copy(pivots),encode(d,g);constant=false) || throw(ArgumentError(
            "$name generator $i does not preserve the continuous span (derivation $j); add the missing derivations or remove the generator."))
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

# Floating coefficients can promote differently with machine and BigInt
# rationals. Preserve the original exact conversion path for such derivations.
function push_term!(s::SectorTerms{W,Complex{Rational{Int}}},c::Union{AbstractFloat,Complex{<:AbstractFloat}},w::W) where W
    throw(InexactError(:continuous_coefficient,Complex{Rational{Int}},c))
end

function _continuous_program(program::MomentProgram,table,derivations;diagonal=true,algebra=nothing)
    if isnothing(algebra)
        W=eltype(table.representatives.keys)
        algebra=W<:PauliWord ? Pauli() : W<:FermiWord ? Fermionic() : Bosonic()
    end
    # Closing the word domain may extend the table. A failed machine attempt
    # must not change the BigInt retry's coordinate or equation order.
    trial=deepcopy(table)
    try
        return _continuous_program_rational(program,trial,derivations,Rational{Int};diagonal,algebra)
    catch error
        error isa Union{OverflowError,InexactError} || rethrow()
        return _continuous_program_rational(program,table,derivations,Rational{BigInt};diagonal,algebra)
    end
end

# Translating [S,w] gives [S,T(w)] only if T preserves S. Verify the
# signed one-site terms of every spin sum before using this identity.
function _ward_translation(table::OrbitMoments{W},derivations) where W
    isnothing(table.cyclic) && return nothing
    t=table.cyclic.generator
    is_antiunitary(t) && return nothing
    for d in derivations
        d isa SpinRotation || return nothing
        terms=Set(pauli_letter(W,i,d.axis) for i in MaskSites(d.support))
        for w in terms
            p,v=act(t,w)
            isone(p) && v in terms || return nothing
        end
    end
    return t
end

# A representative table may opt in only when the finite action preserves
# the span of the supplied derivations under conjugation. Return one word
# per registered orbit, including zero orbits: D(w) need not vanish there.
# These seeds must already be canonical under _ward_key.
_ward_seed_words(table,derivations) = nothing
_ward_key(table,w) = w

function _continuous_program_rational(program::MomentProgram{C},table,derivations,::Type{R};diagonal=true,algebra) where {C,R}
    W=eltype(table.representatives.keys)
    initial=program.nvars
    support=sort!(unique!([site for label in table.labels for site in sites(label.word)]))
    letters=basis_rows(W,algebra,support,1)
    diagonal &= all(d->all(w->all(==(w),keys(_derivative(d,w))),letters),derivations)
    pivots=Dict{Int,Dict{Int,R}}()
    derivative=SectorTerms(Dict{W,Complex{R}}())
    image=Pair{W,Complex{R}}[]
    sink=OrbitFormSink(table,Dict{Int32,Complex{R}}())
    row=Dict{Int,R}()
    # Close words, not only nonzero orbit coordinates: D(w) can be nonzero
    # when finite invariance forces w itself to zero (a normalizer can rotate
    # the derivations). Such equations must not disappear with that orbit.
    representatives=_ward_seed_words(table,derivations)
    if isnothing(representatives)
        queue=sort!(unique!([collect(keys(table.cache));[label.word for label in table.labels]]))
        seen=Set(queue)
        translation=_ward_translation(table,derivations)
    else
        queue=sort!(representatives)
        seen=Set(queue)
        translation=nothing
    end
    translated=Set{W}()
    cursor=1
    labels_seen=length(table.labels)
    while cursor<=length(queue)
        w=queue[cursor]
        cursor+=1
        if !isnothing(translation)
            w in translated && continue
            # Only aliases under T, not the whole finite orbit: a word
            # forced to zero can still generate an independent equation.
            v=w
            while true
                push!(translated,v)
                _,v=act(translation,v)
                v==w && break
            end
        end
        for d in derivations
            empty!(derivative.terms)
            derive!(derivative,d,w,one(Complex{R}))
            if diagonal
                isempty(derivative.terms) && continue
                m=orbit_coordinates!(table,w)
                for id in (m.re,m.im)
                    id==0 && continue
                    id==1 && throw(ArgumentError("A derivation must annihilate the identity."))
                    pivots[id]=Dict(Int(id)=>one(R))
                end
            else
                empty!(image); append!(image,derivative.terms)
                empty!(sink.terms)
                for (v,c) in sort!(image;by=first)
                    key=_ward_key(table,v)
                    if !(key in seen)
                        # Keep the first encountered image as the seed, even
                        # when its membership key is an orbit representative.
                        push!(seen,key); push!(queue,v)
                    end
                    push_term!(sink,c,v)
                end
                equation=sink.terms
                for part in (real,imag)
                    empty!(row)
                    for (id,c) in equation
                        iszero(part(c)) || (row[Int(id)]=part(c))
                    end
                    # A new pivot owns its row; dependent rows reuse the buffer.
                    _insert_generator_equation!(pivots,row) || (row=Dict{Int,R}())
                end
            end
        end
        while labels_seen<length(table.labels)
            labels_seen+=1
            v=table.labels[labels_seen].word
            key=_ward_key(table,v)
            if !(key in seen)
                push!(seen,key); push!(queue,v)
            end
        end
    end
    free=[id for id in eachindex(table.labels) if !haskey(pivots,id)]
    directions=[Pair{Int32,R}[] for _ in table.labels]
    for (id,old) in enumerate(free)
        push!(directions[old],Int32(id)=>one(R))
    end
    for pivot in sort!(collect(keys(pivots));rev=true)
        coefficients=Dict{Int32,R}()
        for (old,c) in pivots[pivot]
            old==pivot && continue
            for (id,factor) in directions[old]
                coefficients[id]=get(coefficients,id,zero(c))-c*factor
            end
        end
        directions[pivot]=sort!([id=>c for (id,c) in coefficients if !iszero(c)];by=first)
    end
    # Keep the same rational-to-Float64 rounding as the original implementation.
    exact_directions=convert(Vector{Vector{Pair{Int32,Rational{BigInt}}}},directions)
    float_directions=[Pair{Int32,Float64}[id=>Float64(factor) for (id,factor) in row]
        for row in exact_directions]
    out=MomentProgram{C}(length(free))
    function rewrite(ids,coefs)
        result=Dict{Int32,C}()
        for (old,c) in zip(ids,coefs),(id,factor) in float_directions[old]
            result[id]=get(result,id,zero(C))+c*factor
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
    reduced=ContinuousMoments(table,table.labels[free],exact_directions,derivations,BitSet([1]))
    refresh_coverage!(reduced,out)
    report=(method=diagonal ? :diagonal : :exact_elimination,initial_variables=initial-1,
        closed_variables=length(table.labels)-1,variables=length(free)-1,equations=length(pivots),
        arithmetic=R==Rational{Int} ? :checked_machine_rationals : :complex_rationals)
    return validate(out),reduced,report
end
