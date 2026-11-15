"""
    Sector(operator, eigenvalue)

Restrict states to an eigenspace of a Hermitian Pauli, fermionic or bosonic operator.
For particle number use `Sector(N, n)`; for singlets use `Sector(S², 0)`.
The operator must use the problem's registry. This is a state restriction,
not the operator identity `Zero(operator - eigenvalue)`.
"""
struct Sector{P<:Polynomial,T<:Real} <: StrengtheningKind
    operator::P
    eigenvalue::T
    function Sector(operator::Polynomial,eigenvalue::Real)
        _checked_hermitian(operator,"Sector operator")
        isfinite(eigenvalue) || throw(ArgumentError("The sector eigenvalue must be finite."))
        new{typeof(operator),typeof(eigenvalue)}(operator,eigenvalue)
    end
end
_validity(k::Sector) = "valid for states supported on the eigenspace ($(k.operator)) = $(k.eigenvalue)"
_strengthening_width(k::Sector,problem,order) = max(2order+degree(k.operator),2degree(k.operator))

# Exact checks use rational versions of the supplied coefficients, including
# the exact binary value of a floating-point input. No commutator tolerance.
struct SectorTerms{W,C}
    terms::Dict{W,C}
end
SectorTerms(::Type{W}) where W = SectorTerms(Dict{W,Complex{Rational{BigInt}}}())
function push_term!(s::SectorTerms{W,C},c::Number,w::W) where {W,C}
    value=get(s.terms,w,zero(C))+c
    iszero(value) ? delete!(s.terms,w) : (s.terms[w]=value)
    s
end
function _sector_terms(p,W)
    s=SectorTerms(W)
    for (c,w) in terms(p)
        push_term!(s,Complex{Rational{BigInt}}(c),_pack_word(W,w))
    end
    s
end
function _sector_product(a::SectorTerms{W},b::SectorTerms{W}) where W
    s=SectorTerms(W)
    for (u,c) in a.terms,(v,d) in b.terms
        mul!(s,u,v,c*d)
    end
    s
end

function _checked_sectors(problem,kinds,symmetry,W)
    result=Tuple{Int,SectorTerms{W,Complex{Rational{BigInt}}}}[]
    for (source,k) in enumerate(kinds)
        k isa Sector || continue
        problem.objective.registry.algebra isa Union{Pauli,Fermionic,Bosonic} ||
            throw(ArgumentError("Sector currently requires Pauli, fermionic or bosonic variables."))
        k.operator.registry===problem.objective.registry ||
            throw(ArgumentError("The sector operator must use the problem registry."))
        _check_derivation_polynomial(k.operator,W,symmetry.continuous,"Sector operator")
        any(x->x isa StateOptimality,kinds) && throw(ArgumentError(
            "StateOptimality with a sector is not supported in this delivery; remove StateOptimality or the sector restriction."))
        s=_sector_terms(k.operator,W)
        push_term!(s,-Rational{BigInt}(k.eigenvalue),one(W))
        for g in (symmetry.generators...,symmetry.moment_only...)
            image=SectorTerms(W)
            for (w,c) in s.terms
                act!(image,g,w,is_antiunitary(g) ? conj(c) : c)
            end
            same=keys(image.terms)==keys(s.terms)
            if same && !isempty(s.terms)
                w=first(keys(s.terms)); ratio=image.terms[w]/s.terms[w]
                same=!iszero(ratio) && all(image.terms[v]==ratio*c for (v,c) in s.terms)
            end
            same || throw(ArgumentError(
                "A supplied symmetry moves the selected sector; remove the generator or choose an invariant sector."))
        end
        if any(x->x isa Stationarity,kinds)
            h=_sector_terms(problem.objective,W)
            _sector_product(h,s).terms==_sector_product(s,h).terms || throw(ArgumentError(
                "Stationarity with a sector requires an exactly zero commutator of the objective and sector operator; remove Stationarity or choose a commuting sector operator."))
        end
        push!(result,(source,s))
    end
    result
end

function _sector_equalities!(program,table,closure,sectors)
    isempty(sectors) && return
    W=eltype(table.representatives.keys)
    domain=SectorTerms(Dict{W,Complex{Int}}())
    words=Set{W}([one(W)])
    for block in closure.family
        if block.source==0
            for a in block.rows,b in block.rows
                empty!(domain.terms)
                phase,left=star(a)
                mul!(domain,left,b,scale(phase,one(Complex{Int})))
                union!(words,keys(domain.terms))
            end
        else
            local_terms=SectorTerms(W)
            polynomial=[(Complex{Rational{BigInt}}(c),w) for (c,w) in block.polynomial]
            for a in block.rows,b in block.rows
                empty!(local_terms.terms)
                _block_product!(local_terms,a,b,polynomial)
                union!(words,keys(local_terms.terms))
            end
        end
    end
    ordered=sort!(collect(words))
    nterms,nforms,nblocks=length(program.forms.ids),length(program.forms.offsets),length(program.blocks)
    try
        _sector_equations!(program,table,ordered,sectors,Rational{Int})
    catch error
        error isa Union{OverflowError,InexactError} || rethrow()
        resize!(program.forms.ids,nterms); resize!(program.forms.coefs,nterms)
        resize!(program.forms.offsets,nforms); resize!(program.blocks,nblocks)
        # The failed traversal registered only a prefix of the same words.
        # Keep those coordinates: replay uses the same order and ids, without
        # copying the entire moment table on the ordinary machine path.
        _sector_equations!(program,table,ordered,sectors,Rational{BigInt})
    end
    nothing
end

function _sector_equations!(program,table,ordered,sectors,::Type{R}) where R
    sink=OrbitFormSink(table,Dict{Int32,Complex{R}}())
    seen=Set{Vector{Pair{Int32,R}}}()
    key=Pair{Int32,R}[]
    # Moment and Hermitian localizer entries occur with their adjoints.
    # Every right equation <u(Q-q)>=0 is therefore already the conjugate
    # (up to a word phase) of a left equation <(Q-q)u†>=0.
    for (source,s) in sectors
        count=0
        # The variance enforces support in the eigenspace even at order zero.
        # Left multiplication is essential: a†(Q-q)b need not vanish when
        # b takes a physical state out of the selected sector.
        variance=(_sector_product(s,s),)
        for u in Iterators.flatten((variance,ordered))
            empty!(sink.terms)
            if u isa SectorTerms
                for (w,c) in u.terms; push_term!(sink,Complex{R}(c),w); end
            else
                for (w,c) in s.terms
                    mul!(sink,w,u,Complex{R}(c))
                end
            end
            for part in (real,imag)
                empty!(key)
                for (id,c) in sink.terms
                    iszero(part(c)) || push!(key,id=>part(c))
                end
                isempty(key) && continue
                sort!(key;by=first)
                pivot=last(first(key))
                for i in eachindex(key)
                    id,c=key[i]; key[i]=id=>c/pivot
                end
                key in seen && continue
                push!(seen,copy(key))
                # Keep the reference rounding at the conic boundary. The
                # sink, normalization and duplicate checks remain exact.
                for (id,c) in key; push_term!(program.forms,Rational{BigInt}(c),id); end
                close_form!(program.forms)
                count+=1
            end
        end
        count==0 || push_block!(program,ZERO,count,source)
    end
end
