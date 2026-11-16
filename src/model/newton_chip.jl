"""
    NCTSSoS.newton_chip_basis(problem, order=cld(degree(problem.objective), 2); solver=nothing)

Select custom rows for an unconstrained Free problem with one noncommuting
group. Pass the result to `Relaxation(rows=rows)`; this is not a second solve
interface. The default order covers the objective; a larger order adds no rows.
Ordinary objectives retain right suffixes of Hermitian-square
support words. Scalar single-trace objectives in `NCTSSoS.RealTrace()` retain
operator words whose exponent vectors are in half the support polytope,
including the origin for the unknown bound.

The trace route requires an explicit solver backend. An inconclusive solve
keeps the row. A numerical separating vector removes a row only after its
strict inequalities are checked with exact rational arithmetic. Products of
traces, residual operator tails, state symbols, constraints, declared symmetry
and additional algebraic relations are unsupported and give `ArgumentError`.
"""
function newton_chip_basis(problem::Problem,order::Integer=cld(degree(problem.objective),2);solver=nothing)
    p=problem.objective; r=p.registry
    order>=0 && order>=cld(degree(p),2) || throw(ArgumentError(
        "Newton-chip order must be nonnegative and at least $(cld(degree(p),2)) to cover the objective."))
    isempty(problem.constraints) && !_has_symmetry(problem.symmetry) || throw(ArgumentError(
        "Newton-chip rows require an unconstrained problem without declared symmetry; use complete rows instead."))
    trace=r.algebra isa FunctionalAlgebra{RealTrace,Free}
    trace || r.algebra isa SiteAlgebra{Free} || throw(ArgumentError(
        "Newton-chip rows support ordinary Free or scalar single-trace RealTrace/Free objectives only; use complete rows for this algebra."))
    trace && !(solver isa AbstractSolver) && throw(ArgumentError(
        "Trace Newton-chip selection needs solver=SCSBackend(), MosekBackend(), or JuMPBackend(factory); no solver is chosen implicitly."))
    W=eltype(p.data.words)
    words=W[]; used_sites=Set{Int}()
    for (_,w) in terms(p)
        if trace && !isone(w)
            _ftail(w)==0 && _fstop(w,1)==degree(w) || throw(ArgumentError(
                "Trace Newton-chip rows require a scalar single-trace objective: remove products of traces and residual operator tails."))
        end
        union!(used_sites,sites(w))
        push!(words,foldl((a,b)->mul(a,b)[2],support_letters(w);init=one(W)))
    end
    length(used_sites)<=1 || throw(ArgumentError(
        "Newton-chip rows require one noncommuting group; cross-site commuting objectives need complete rows."))
    selected=trace ? _cyclic_chip_words(words,Int(min(order,degree(p)÷2)),solver) : _right_chip_words(words,order)
    return map(selected) do w
        buffer=PolynomialBuffer{W,ComplexF64}(1)
        push_term!(buffer,1,w)
        _polynomial(r,buffer)
    end
end

function _right_chip_words(words::Vector{W},order) where {W}
    rows=Set{W}([one(W)])
    for w in words
        letters=collect(word_letters(w)); n=length(letters)
        iseven(n) || continue
        all(letters[i]==letters[n-i+1] for i in 1:n÷2) || continue
        chip=one(W)
        for i in 1:min(n÷2,order)
            chip=mul(letters[n-i+1],chip)[2]
            push!(rows,chip)
        end
    end
    sort!(collect(rows))
end

# Separate 2*target from conv(support). Free variables are the normal h and
# margin t: maximize t with (2*target-v)⋅h >= t for every vertex, |h_i| <= 1.
function _chip_separated(target,support,solver)
    doubled=2 .* target
    doubled in support && return false
    n=length(target)
    any(doubled[i]>maximum(v[i] for v in support) for i in 1:n) && return true
    m=length(support)
    A=zeros(m+2n,n+1); b=zeros(m+2n); c=zeros(n+1); c[end]=-1
    for (j,v) in enumerate(support)
        A[j,1:n].=v.-doubled; A[j,end]=1
    end
    for i in 1:n
        A[m+i,i]=1; A[m+n+i,i]=-1
        b[m+i]=b[m+n+i]=1
    end
    form=MomentStandardForm(sparse(A),b,c,0.,[DataCone(:nonnegative,m+2n,1,m+2n,Int32(1))])
    result=solve_standard(form,solver)
    result.status==JuMP.MOI.OPTIMAL && !isnothing(result.bound) &&
        length(result.values)==n+1 && all(isfinite,result.values) || return false
    # Float64 -> Rational is exact, not a tolerance-based reconstruction.
    normal=Rational{BigInt}.(result.values[1:n])
    all(sum(normal[i]*(doubled[i]-v[i]) for i in 1:n)>0 for v in support)
end

function _cyclic_chip_words(words::Vector{W},order,solver) where {W}
    letters=sort!(unique!(W[l for w in words for l in word_letters(w)]))
    isempty(letters) && return W[one(W)]
    positions=Dict(w=>i for (i,w) in enumerate(letters))
    function exponents(w)
        counts=zeros(Int,length(letters))
        for l in word_letters(w); counts[positions[l]]+=1; end
        Tuple(counts)
    end
    support=[exponents(w) for w in words]
    push!(support,exponents(one(W))); unique!(support)
    limits=[maximum(v[i] for v in support)÷2 for i in eachindex(letters)]
    cache=Dict{typeof(first(support)),Bool}()
    rows=W[one(W)]; level=W[one(W)]
    for _ in 1:order
        next=W[]
        for w in level,l in letters
            v=mul(w,l)[2]
            exponent=exponents(v)
            any(exponent[i]>limits[i] for i in eachindex(limits)) && continue
            push!(next,v)
            keep=get!(cache,exponent) do
                !_chip_separated(exponent,support,solver)
            end
            keep && push!(rows,v)
        end
        level=next
    end
    sort!(rows)
end
