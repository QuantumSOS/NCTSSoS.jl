function _site_action(g,site_letters,index;pure=false)
    permutation=Int[]
    for letters in site_letters
        target=0
        for (label,w) in enumerate(letters)
            p,u=act(g,w)
            location=get(index,u,(0,0))
            location[1]==0 && return nothing
            pure && (!isone(p) || location[2]!=label) && return nothing
            target!=0 && location[1]!=target && return nothing
            target=location[1]
        end
        push!(permutation,target)
    end
    length(unique(permutation))==length(permutation) ? permutation : nothing
end

"""
    discover_lattice(generators, site_letters; candidate_limit=256)

Find a regular abelian site-permutation subgroup of rank at most two using
generators and products of two generators. `site_letters[site]` lists the same
local letter labels at every site; this witness comes from the registry or
detection helper, not from inspecting an action's storage. Check commutation,
independence, transitivity, and normality before using the lattice route.

Maps outside the discovered lattice normalizer are returned explicitly as
moment-only generators. Their degree-changing orbit relations remain active;
they are not incorrectly advertised as block symmetries. The caller can inspect
this partition. A failed bounded search asks for named translations.
"""
function discover_lattice(generators::AbstractVector{<:AbstractAction},site_letters;candidate_limit=256)
    isempty(generators) && throw(ArgumentError("Supply symmetry generators."))
    n=length(site_letters)
    n>1 && all(length(ws)==length(first(site_letters))>0 for ws in site_letters) ||
        throw(ArgumentError("Supply a nonempty, equal local letter list at each site."))
    W=typeof(first(first(site_letters)))
    generators=map(g->_group_action(W,g),generators)
    A=eltype(generators)
    index=Dict(w=>(site,label) for (site,ws) in enumerate(site_letters) for (label,w) in enumerate(ws))
    length(index)==sum(length,site_letters) || throw(ArgumentError("Site-letter witnesses must be distinct."))
    candidate_limit>0 || throw(ArgumentError("The candidate limit must be positive."))
    pure=[g for g in generators if !is_antiunitary(g) && !isnothing(_site_action(g,site_letters,index;pure=true))]
    candidates=collect(Iterators.take(unique(pure),candidate_limit))
    for g in pure,h in pure
        length(candidates)>=candidate_limit && break
        p=compose(g,h)
        !isone(p) && !(p in candidates) && push!(candidates,p)
    end
    filter!(g->!isone(g),candidates)
    orders=map(candidates) do g
        p=one(g)
        for order in 1:n
            p=compose(g,p)
            isone(p) && return order
        end
        0
    end
    choices=[Int[i] for i in eachindex(candidates) if orders[i]==n]
    append!(choices,[[i,j] for i in eachindex(candidates) for j in (i+1):length(candidates)
        if orders[i]>1 && orders[j]>1 && orders[i]*orders[j]==n && compose(candidates[i],candidates[j])==compose(candidates[j],candidates[i])])
    for choice in choices
        translations=try
            TranslationGroup(candidates[choice],orders[choice])
        catch error
            error isa ArgumentError || rethrow()
            continue
        end
        locations=[_site_action(g,site_letters,index;pure=true)[1] for g in translations.elements]
        length(unique(locations))==n || continue
        to_site=Dict(locations[i]=>translations.elements[i] for i in eachindex(locations))
        point=A[]; remainder=A[]
        for g in generators
            permutation=_site_action(g,site_letters,index)
            normalizes=!is_antiunitary(g) && !isnothing(permutation) && all(
                haskey(translations.index,_conjugate_action(g,t)) for t in translations.generators)
            if normalizes
                h=compose(inv(to_site[permutation[1]]),g)
                !isone(h) && !(h in point) && push!(point,h)
            else
                push!(remainder,g)
            end
        end
        isempty(point) && push!(point,one(first(candidates)))
        return (;translations,point_generators=point,moment_generators=remainder,
            searched_candidates=length(candidates),source=:verified_generator_search)
    end
    throw(ArgumentError("No regular abelian translation subgroup was found in the bounded generator/product search. Supply named translations and their orders, or use the generic finite route."))
end

function generators_program(rows,objective,generators,site_letters)
    structure=discover_lattice(generators,site_letters)
    p,table,report=abelian_lattice_program(rows,objective,structure.translations,structure.point_generators;
        moment_generators=Tuple(structure.moment_generators))
    p,table,(;report...,discovery=structure)
end

"Check selected character labels under the normalizer; this does not impose a physical state charge."
function character_sector_closure(group::TranslationGroup,point_generators,selected;union=false)
    all(k->1<=k<=length(group.elements),selected) || throw(ArgumentError("A selected character index is outside the group."))
    point=OrbitGroup(point_generators)
    permutations=_character_permutations(group,point)
    closure=sort!(unique([p[k] for k in selected for p in permutations]))
    union || closure==sort!(unique(collect(selected))) || throw(ArgumentError(
        "The point group moves a selected character sector. Select its full orbit explicitly or request union=true; selecting row sectors does not impose a physical state charge."))
    closure
end
