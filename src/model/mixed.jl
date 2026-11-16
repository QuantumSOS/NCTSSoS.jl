"""
    NCTSSoS.mixed_variables(; name=(algebra, sites), ...)

Declare named commuting tensor components in one registry. For example,
`r, p = NCTSSoS.mixed_variables(spin=(Pauli(),1), matter=(Fermionic(),2))` gives
`X,Y,Z = p.spin` and `a,ad = p.matter`. A site component uses
`(SiteAlgebra(Unipotent(),2),2)` and returns vectors of letters per site.
Names in the registry are qualified, for example `Symbol("spin.X1")`.

Different components commute, but moments of their products remain joint
moments and do not factor. At most one fermionic component is allowed: combine
all mutually anticommuting modes there. Functional state/trace components,
symmetry, block routing, strengthening and GNS are not yet supported; use
plain `Relaxation(order=...)` or custom rows and `expectation`.
"""
function mixed_variables(;components...)
    isempty(components) && throw(ArgumentError("Name at least one component, for example spin=(Pauli(),1)."))
    specs=Tuple(values(components))
    for (name,spec) in pairs(components)
        spec isa Tuple && length(spec)==2 && spec[1] isa Union{Pauli,Fermionic,Bosonic,SiteAlgebra} &&
            spec[2] isa Integer && spec[2]>0 || throw(ArgumentError(
                "Component $name needs (Pauli(), n), (Fermionic(), n), (Bosonic(), n), or (SiteAlgebra(rule, settings), n) with n > 0; omit state/trace functionals."))
    end
    count(s->s[1] isa Fermionic,specs)<=1 || throw(ArgumentError(
        "Combine mutually anticommuting modes in one Fermionic component; separate tensor components commute."))
    algebras=map(first,specs); sizes=map(s->Int(last(s)),specs)
    algebra=MixedAlgebra{typeof(algebras),sizes}(algebras)
    offset=0
    entries=Pair{Symbol,Tuple{Int,Int}}[]
    groups=map(Tuple(keys(components)),specs) do name,(a,n)
        axes=a isa Pauli ? ((:X,0),(:Y,1),(:Z,2)) :
            a isa Fermionic ? ((:a,1),(:ad,-1)) : ((:b,1),(:bd,-1))
        names=if a isa SiteAlgebra
            [[Symbol(name,".s",i,"_",j) for j in 1:a.nletters] for i in 1:n]
        else
            Tuple([Symbol(name,".",axis,i) for i in 1:n] for (axis,_) in axes)
        end
        if a isa SiteAlgebra
            for i in 1:n,j in 1:a.nletters; push!(entries,names[i][j]=>(offset+i,j)); end
        else
            for (k,(_,label)) in enumerate(axes),i in 1:n
                push!(entries,names[k][i]=>(offset+i,label))
            end
        end
        offset+=n
        names
    end
    r=Registry(algebra,sum(sizes),entries)
    vars=variables(r)
    parts=map(group->map(names->[vars[name] for name in names],group),groups)
    return r,NamedTuple{Tuple(keys(components))}(parts)
end
function _model_word_type(r::Registry{MixedAlgebra{A,N}}) where {A,N}
    types=map((a,n)->_model_word_type(Registry(a,n,())),r.algebra.components,N)
    return MixedWord{Tuple{types...},_mixed_offsets(N)}
end
_locations(w::MixedWord{T,O}) where {T,O} = Iterators.flatten(
    map((part,offset)->((s+offset,l) for (s,l) in _locations(part)),w.components,Base.front(O)))
