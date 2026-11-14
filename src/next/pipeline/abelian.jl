"A finite product of cyclic translations; coordinates are checked, not inferred."
struct TranslationGroup{A}
    generators::Vector{A}
    orders::Vector{Int}
    elements::Vector{A}
    index::Dict{A,Int}
    coordinates::Vector{Vector{Int}}
    modulus::Int
    characters::Matrix{Int} # exact phase numerators modulo 4*lcm(orders)
end
function TranslationGroup(generators::AbstractVector{A}, orders) where {A<:AbstractAction}
    length(generators)==length(orders)>0 || throw(ArgumentError("Give one positive order per translation."))
    all(>(0),orders) || throw(ArgumentError("Translation orders must be positive."))
    any(is_antiunitary,generators) && throw(ArgumentError("Translations must be unitary."))
    for g in generators,h in generators
        compose(g,h)==compose(h,g) || throw(ArgumentError("The named translations do not commute."))
    end
    powers=map(zip(generators,orders)) do (g,n)
        ps=A[one(g)]
        for _ in 1:n
            push!(ps,compose(g,last(ps)))
        end
        isone(pop!(ps)) || throw(ArgumentError("A declared translation order does not close."))
        ps
    end
    coords=[collect(t) for t in vec(collect(Iterators.product((0:(n-1) for n in orders)...)))]
    elements=Vector{A}(undef,length(coords))
    for (j,t) in enumerate(coords)
        g=one(first(generators))
        for i in eachindex(t)
            g=compose(powers[i][t[i]+1],g)
        end
        elements[j]=g
    end
    index=Dict{A,Int}(g=>i for (i,g) in enumerate(elements))
    length(index)==length(elements) || throw(ArgumentError("Translation generators are dependent; give independent cyclic factors."))
    L=lcm(orders...)
    chars=[mod(4sum(k[i]*t[i]*div(L,orders[i]) for i in eachindex(k)),4L) for k in coords,t in coords]
    TranslationGroup(collect(generators),Int.(orders),elements,index,coords,4L,chars)
end

struct AbelianRows{W,T}
    translations::T
    words::Vector{Vector{W}}
    phases::Vector{Vector{Phase}}
    positions::Vector{Vector{Int}}
    admissible::BitMatrix
end
function AbelianRows(rows::AbstractVector{W}, translations::TranslationGroup) where {W<:AbstractWord}
    domain=Set(rows); seen=Set{W}()
    length(domain)==length(rows) || throw(ArgumentError("The row list contains duplicates."))
    words=Vector{W}[]; phases=Vector{Phase}[]; positions=Vector{Int}[]; admissible=BitVector[]
    for rep in sort!(collect(domain))
        rep in seen && continue
        ws=W[]; ps=Phase[]; ts=Int[]; index=Dict{W,Int}()
        allowed=trues(length(translations.elements))
        for (t,g) in enumerate(translations.elements)
            p,w=act(g,rep)
            w in domain || throw(ArgumentError("Translations do not preserve the supplied row family."))
            if haskey(index,w)
                j=index[w]; q=p*conj(ps[j])
                for k in eachindex(allowed)
                    allowed[k] &= mod(translations.characters[k,t]-translations.characters[k,ts[j]]-
                        Int(q.k)*div(translations.modulus,4),translations.modulus)==0
                end
            else
                push!(ws,w); push!(ps,p); push!(ts,t); index[w]=length(ws)
            end
        end
        union!(seen,ws); push!(words,ws); push!(phases,ps); push!(positions,ts); push!(admissible,allowed)
    end
    AbelianRows(translations,words,phases,positions,BitMatrix(hcat(admissible...)))
end

function _character_permutations(translations, point)
    generator_indices=[translations.index[g] for g in translations.generators]
    signatures=Dict(Tuple(translations.characters[k,generator_indices])=>k for k in axes(translations.characters,1))
    permutations=Vector{Int}[]
    for h in point.elements
        conjugates=Int[]
        for g in translations.generators
            c=_conjugate_action(inv(h),g)
            haskey(translations.index,c) || throw(ArgumentError("The translations are not a normal subgroup of the supplied group."))
            push!(conjugates,translations.index[c])
        end
        push!(permutations,[signatures[Tuple(translations.characters[k,conjugates])] for k in axes(translations.characters,1)])
    end
    permutations
end

function _momentum_action(orbits::AbelianRows,index,h,k;before=nothing)
    tr=orbits.translations
    map(orbits.words) do ws
        p,w=_action_image(h,first(ws),before)
        haskey(index,w) || throw(ArgumentError("The point map does not preserve the row family."))
        b,t=index[w]; position=orbits.positions[b][t]
        q=p*conj(orbits.phases[b][t])
        (scale(q,cispi(2tr.characters[k,position]/tr.modulus)),b)
    end
end

function _orbit_pair_kernel!(cache,table,orbits::AbelianRows,a,b)
    get!(cache,(a,b)) do
        forms=LinearForms{ComplexF64}(); sink=OrbitFormSink(table,Dict{Int32,ComplexF64}())
        wa,wb=orbits.words[a],orbits.words[b]; p,sa=star(first(wa))
        for t in eachindex(wb)
            empty!(sink.terms)
            mul!(sink,sa,wb[t],scale(p*orbits.phases[b][t],ComplexF64(sqrt(length(wa)/length(wb)))))
            _close_orbit_form!(forms,sink.terms)
        end
        forms
    end
end

function _fourier_kernel!(fourier,cache,table,orbits::AbelianRows,a,b,characters)
    get!(fourier,(a,b)) do
        raw=_orbit_pair_kernel!(cache,table,orbits,a,b)
        terms=Dict{Int32,ComplexF64}()
        for t in 1:length(raw),j in term_range(raw,t)
            id=raw.ids[j]
            terms[id]=get(terms,id,0im)+characters[orbits.positions[b][t]]*raw.coefs[j]
        end
        filter!(p->abs(last(p))>1e-12,terms)
        terms
    end
end

"""
Direct Fourier/little-group construction for a finite product of cyclic groups.
`selected_characters` selects PSD row blocks, never a physical state-charge
sector. A selection must contain whole point-group character orbits, unless
`union_characters=true` requests their closure. Omitting blocks weakens the
relaxation; it does not add particle-number or other state-sector equations.
"""
function abelian_lattice_program(rows::AbstractVector{W},objective,translations::TranslationGroup,
    point_generators::AbstractVector{<:AbstractAction};moment_generators=(),
    selected_characters=nothing,union_characters=false) where {W<:AbstractWord}
    canonical=map(g->_group_action(W,g),translations.generators)
    if any(canonical[i] !== translations.generators[i] for i in eachindex(canonical))
        translations=TranslationGroup(canonical,translations.orders)
    end
    point_generators=map(g->_group_action(W,g),point_generators)
    all_generators=(translations.generators...,point_generators...,moment_generators...)
    _check_objective(objective,all_generators)
    orbits=AbelianRows(rows,translations)
    point=OrbitGroup(point_generators)
    permutations=_character_permutations(translations,point)
    selected=isnothing(selected_characters) ? collect(axes(translations.characters,1)) :
        character_sector_closure(translations,point_generators,selected_characters;union=union_characters)
    index=Dict(w=>(a,t) for (a,ws) in enumerate(orbits.words) for (t,w) in enumerate(ws))
    table=OrbitMoments(W,all_generators)
    program=MomentProgram{ComplexF64}(); sink=OrbitFormSink(table,Dict{Int32,ComplexF64}())
    cache=Dict{Tuple{Int,Int},LinearForms{ComplexF64}}()
    projection_cache=Dict{Tuple{Vararg{Int}},Tuple{typeof(point),Vector{Vector{ComplexF64}}}}()
    sector_labels=Tuple{Int,Int}[]; seen=Set{Int}(); little_orders=Int[]
    anti_index=findfirst(is_antiunitary,moment_generators)
    if !isnothing(anti_index)
        anti=moment_generators[anti_index]
        all(_conjugate_action(anti,t)==t for t in translations.generators) ||
            throw(ArgumentError("The lattice antiunitary must commute with each named translation."))
    end
    for k in selected
        k in seen && continue
        union!(seen,(p[k] for p in permutations))
        stabilizer=Tuple(i for (i,p) in enumerate(permutations) if p[k]==k)
        little,projections=get!(projection_cache,stabilizer) do
            group=OrbitGroup(point.elements[collect(stabilizer)])
            ps,_=shared_projections(group)
            (group,ps)
        end
        active=findall(orbits.admissible[k,:])
        actions=[_momentum_action(orbits,index,h,k) for h in little.elements]
        sectors=_momentum_projectors(active,actions,projections)
        if !isnothing(anti_index)
            negative=findfirst(j->all(iszero(mod(translations.characters[k,t]+translations.characters[j,t],translations.modulus))
                for t in axes(translations.characters,2)),axes(translations.characters,1))
            candidates=point.elements[[i for (i,p) in enumerate(permutations) if p[negative]==k]]
            sectors=_real_character_sectors(sectors,orbits,index,k,moment_generators[anti_index],candidates)
        end
        characters=cispi.(-2 .* translations.characters[k,:] ./ translations.modulus)
        _assemble_character!(program,sink,orbits,characters,k,sectors,sector_labels,cache;real_blocks=!isnothing(anti_index))
        push!(little_orders,length(little))
    end
    p,table=_orbit_objective!(program,table,objective)
    return _anti_real_program(p,all_generators),table,(;sectors=sector_labels,little_orders,
        translation_order=length(translations.elements),point_order=length(point),selected_characters=selected)
end
