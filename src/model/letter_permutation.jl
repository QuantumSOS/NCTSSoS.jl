LetterPermutation(pairs::Pair...) = LetterPermutation(pairs)
function LetterPermutation(pairs)
    entries=collect(pairs)
    !isempty(entries) && all(p->p isa Pair,entries) || throw(ArgumentError(
        "LetterPermutation needs operator-letter pairs; use one(g) for an identity action."))
    first(entries).first isa Polynomial || throw(ArgumentError(
        "LetterPermutation pair $(first(entries)) needs operator-letter polynomials, not integer codes."))
    r=first(entries).first.registry
    r.algebra isa Union{SiteAlgebra,FunctionalAlgebra} || throw(ArgumentError(
        "LetterPermutation is for Free, Projector or Unipotent letters; use SitePermutation for Pauli, fermionic or bosonic sites, or CliffordMap for Pauli settings."))
    base=r.algebra isa FunctionalAlgebra ? r.algebra.base : r.algebra
    n,m=r.nsites,base.nletters
    images=[(s,l) for s in 1:n,l in 1:m]
    sources=Set{Tuple{Int,Int}}()
    negative=Set{Tuple{Int,Int}}()
    name(location)=_letter_name(r.names,location)
    pairname(location)="$(name(location)) => $(name(images[location...]))"
    for pair in entries
        for (i,endpoint) in enumerate(pair)
            endpoint isa Polynomial && endpoint.registry===r || throw(ArgumentError(
                "LetterPermutation pair $pair must use operator letters from one registry."))
            length(endpoint.data)==1 && (endpoint.data.coefficients[1]==1 || (i==2 && endpoint.data.coefficients[1]==-1)) &&
                degree(endpoint.data.words[1])==1 &&
                (!(endpoint.data.words[1] isa FunctionalWord) || _ftail(endpoint.data.words[1])==1) ||
                throw(ArgumentError("LetterPermutation pair $pair needs coefficient-one source operator letters and target coefficients +1 or -1, not symbols, sums or other scaled letters."))
        end
        source,target=_letter_location(pair.first.data.words[1]),_letter_location(pair.second.data.words[1])
        source in sources && throw(ArgumentError("LetterPermutation pair $pair repeats its source."))
        if pair.second.data.coefficients[1]==-1
            base isa SiteAlgebra{Projector} && throw(ArgumentError(
                "LetterPermutation pair $pair has a negative Projector target: -P is not a projector; remove the sign or use Free or Unipotent variables."))
            push!(negative,source)
        end
        push!(sources,source); images[source...]=target
    end
    seen=Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for s in 1:n,l in 1:m
        source=(s,l); target=images[s,l]
        if haskey(seen,target)
            bad=source in sources ? source : seen[target]
            throw(ArgumentError("LetterPermutation pair $(pairname(bad)) is not bijective when omitted letters are fixed; another source also maps to $(name(target))."))
        end
        seen[target]=source
    end
    for s in 1:n,l in 2:m
        if images[s,l][1]!=images[s,1][1]
            bad=(s,l) in sources ? (s,l) : (s,1)
            throw(ArgumentError("LetterPermutation pair $(pairname(bad)) splits the commuting-site partition: all letters of site $s must map to one site."))
        end
    end
    site_map=SitePermutation([images[s,1][1] for s in 1:n])
    labels=[SitePermutation([images[s,l][2] for l in 1:m]) for s in 1:n]
    LetterPermutation{_letter_rule(eltype(first(entries).first.data.words))}(site_map,labels,negative,r.names,Val(:normal))
end
