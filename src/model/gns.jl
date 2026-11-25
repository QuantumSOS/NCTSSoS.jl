"""
    gns(result::SolveResult; atol=nothing, rtol=0.0, check_atol=1e-6, psd_rtol=1e-8)

Reconstruct a finite-dimensional model from solved moments. The first covered,
positive, flat adjacent degree pair is selected automatically. A nonflat
candidate is also returned if its matrices satisfy the problem and attain the
solved bound. These are numerical checks at the reported tolerances, not a
rigorous certificate. `verdict` distinguishes `:moments_reproduced`,
`:optimizer_certified`, and `:withheld`.
With `atol=nothing`, a resolved singular-value gap above the numerical tail
sets the rank cutoff; otherwise `atol` sets an explicit absolute cutoff.
The automatic gap must exceed 1000, and its discarded singular value must be
at most `check_atol * max(1, σmax)`. `rtol * σmax` is a lower bound on the
cutoff. `check_atol` bounds the absolute matrix, moment and objective errors;
`psd_rtol * max(1, σmax)` bounds negativity of the full Hankel matrix.
Use `g[x]`, `evaluate(g, p)` and `expectation(g, p)` to read an accepted model.
"""
function gns end

"""
Numerical GNS reconstruction and its diagnostics; not a rigorous certificate.

`operators` and `state` are `nothing` for `verdict == :withheld`; otherwise
index the result with a variable (`g[x]`) to obtain its matrix. `rank` and
`singular_values` describe the principal Hankel matrix; `full_rank` and
`full_singular_values` describe its adjacent extension. `principal_degree`,
`full_degree` and `relaxation_order` record the selected pair and solved order.

`flatness` records the verdict, rank cutoffs and gaps, PSD tolerance, minimum
eigenvalue and extension residuals. `checks` records both Hankel and complete
program moment errors, objective error, constraint minimum eigenvalues, Zero
constraint norms, sector residuals, algebra/adjoint and state-norm errors,
distance to flatness and conditioning. `advice` explains a withheld result.
"""
struct GNSResult{R,W,F,C}
    registry::R
    operators::Union{Nothing,Dict{W,Matrix{ComplexF64}}}
    state::Union{Nothing,Vector{ComplexF64}}
    basis::Vector{W}
    full_basis::Vector{W}
    rank::Int
    full_rank::Int
    singular_values::Vector{Float64}
    full_singular_values::Vector{Float64}
    principal_degree::Int
    full_degree::Int
    relaxation_order::Int
    verdict::Symbol
    status::Symbol
    flatness::F
    checks::C
    advice::String
end
public GNSResult

function Base.show(io::IO,g::GNSResult)
    print(io,"GNS: ")
    g.verdict===:optimizer_certified && print(io,"optimizer certified within tolerance (objective equals bound, constraints hold); ")
    print(io,"moments ",g.verdict===:moments_reproduced ? "reproduced" : "NOT reproduced",
        "; max error ",round(g.checks.program_moment_error;sigdigits=3))
    g.status===:not_positive && print(io,"; Hankel matrix not positive")
    g.status===:not_flat && print(io,"; not flat")
    print(io,"; ranks ",g.rank," → ",g.full_rank,
        "; degrees (",g.full_degree,", ",g.principal_degree,")")
    g.full_degree>g.relaxation_order && print(io," above relaxation order ",g.relaxation_order)
    isnothing(g.operators) && print(io,"; operators withheld. ",g.advice)
end
Base.show(io::IO,::MIME"text/plain",g::GNSResult)=show(io,g)

function _gns_accepted(g::GNSResult)
    isnothing(g.operators) && throw(ArgumentError("GNS operators are withheld ($(g.status)). $(g.advice)"))
    return nothing
end
function Base.getindex(g::GNSResult,x::Polynomial)
    _gns_accepted(g)
    x.registry===g.registry || throw(ArgumentError("GNS variable uses a different registry; use a variable of the solved problem."))
    length(x.data)==1 && x.data.coefficients[1]==1 && degree(x.data.words[1])==1 ||
        throw(ArgumentError("Index GNS with a single variable; use evaluate(g, p) for a polynomial."))
    return g.operators[_pack_word(eltype(g.basis),x.data.words[1])]
end

function _gns_matrix(operators,w::AbstractWord,rank)
    out=Matrix{ComplexF64}(LinearAlgebra.I,rank,rank)
    for l in word_letters(w)
        out=out*operators[l]
    end
    return out
end

# Advance the bra in normal-form order: (ψ† A B)† = B† A† ψ.
# Two caller-owned vectors suffice, even when a fermionic word has many letters.
function _gns_observe!(a,b,operators,w,state)
    copyto!(a,state)
    for l in word_letters(w)
        LinearAlgebra.mul!(b,adjoint(operators[l]),a)
        a,b=b,a
    end
    return LinearAlgebra.dot(a,state)
end

function _gns_evaluate(operators,p::Polynomial,rank,::Type{W}) where W
    out=zeros(ComplexF64,rank,rank)
    for (c,w) in terms(p)
        out .+= c .* _gns_matrix(operators,_pack_word(W,w),rank)
    end
    return out
end
"Evaluate a user polynomial in an accepted GNS model, returning its matrix."
function evaluate(g::GNSResult,p::Polynomial)
    _gns_accepted(g)
    p.registry===g.registry || throw(ArgumentError("GNS polynomial uses a different registry; use the solved problem's variables."))
    return _gns_evaluate(g.operators,p,g.rank,eltype(g.basis))
end
"Evaluate a polynomial in the reconstructed state. For an optimizer-only verdict this can differ from expectation(solved_result, p)."
function expectation(g::GNSResult,p::Polynomial)
    _gns_accepted(g)
    p.registry===g.registry || throw(ArgumentError("GNS polynomial uses a different registry; use the solved problem's variables."))
    a,b=similar(g.state),similar(g.state)
    value=0.0im
    for (c,w) in terms(p)
        value+=c*_gns_observe!(a,b,g.operators,_pack_word(eltype(g.basis),w),g.state)
    end
    return LinearAlgebra.ishermitian(p) ? real(value) : value
end

# Enumerate the finite stored domain, but read every value exclusively through
# the coordinate hooks. Orbit caches include aliases; continuous tables retain
# the base domain even after eliminating coordinates.
_gns_domain(t::MomentTable)=[t.keys[i] for i in eachindex(t.ids) if t.ids[i]!=(0,0)]
_gns_domain(t::OrbitMoments)=_orbit_domain(t)
_gns_domain(t::ContinuousMoments)=_gns_domain(t.base)
_gns_domain(t::SectorMoments)=_gns_domain(t.base)

function _gns_coordinates(table,words,coefficients)
    coordinates=Dict{Int,Number}()
    for (w,c) in zip(words,coefficients)
        iszero(c) && continue
        pairs=moment_terms(table,w)
        isnothing(pairs) && return nothing
        for (id,factor) in pairs
            coordinates[id]=get(coordinates,id,0)+c*factor
        end
    end
    return coordinates
end
function _gns_value(result,coordinates;check_coverage=true)
    isnothing(coordinates) && return nothing
    value=0.0+0.0im
    for (id,c) in coordinates
        iszero(c) && continue
        if !moment_covered(result.table,id)
            check_coverage && return nothing
            continue
        end
        value+=ComplexF64(c)*(id==1 ? 1.0 : result.solution.values[id-1])
    end
    return value
end

mutable struct GNSProductValue{W}
    values::Dict{W,ComplexF64}
    value::ComplexF64
end
@inline function push_term!(s::GNSProductValue,c::Number,w::AbstractWord)
    s.value+=c*s.values[w]
    return s
end
@inline function _gns_entry!(s::GNSProductValue,a,b)
    s.value=0.0im
    p,u=star(a)
    mul!(s,u,b,scale(p,1.0+0.0im))
    return s.value
end
function _gns_fill!(H,rows,s::GNSProductValue)
    for j in eachindex(rows), i in eachindex(rows)
        H[i,j]=_gns_entry!(s,rows[i],rows[j])
    end
    return H
end

function _gns_required!(values,result,rows,work,partial=Set{eltype(rows)}())
    for b in rows, a in rows
        empty!(work)
        p,u=star(a)
        mul!(work,u,b,scale(p,1.0+0.0im))
        canonicalize!(work)
        needs_check=false
        for i in 1:length(work)
            w=work.words[i]
            needs_check |= w in partial
            haskey(values,w) && continue
            coordinates=_gns_coordinates(result.table,(w,),(1,))
            isnothing(coordinates) && return w
            value=_gns_value(result,coordinates)
            if isnothing(value)
                push!(partial,w)
                needs_check=true
                # An arbitrary zero for a free coordinate is safe only after
                # checking that it cancels from the complete product below.
                value=_gns_value(result,coordinates;check_coverage=false)
            end
            values[w]=value
        end
        if needs_check
            coordinates=_gns_coordinates(result.table,view(work.words,1:length(work)),view(work.coefficients,1:length(work)))
            isnothing(_gns_value(result,coordinates)) && return first(w for w in view(work.words,1:length(work)) if w in partial)
        end
    end
    return nothing
end

function _gns_algebra_error(operators,letters,rank,work)
    rank==0 && return Inf
    error=0.0
    for a in letters
        p,u=star(a)
        error=max(error,LinearAlgebra.opnorm(operators[a]' - scale(p,1.0+0im)*operators[u]))
        for b in letters
            empty!(work); mul!(work,a,b,1.0+0im)
            expected=zeros(ComplexF64,rank,rank)
            for i in 1:length(work)
                expected .+= work.coefficients[i] .* _gns_matrix(operators,work.words[i],rank)
            end
            error=max(error,LinearAlgebra.opnorm(operators[a]*operators[b]-expected))
        end
    end
    return error
end

function _gns_rank_advice(principal,full)
    gaps(s)=length(s)<2 ? 0 : argmax(s[1:end-1]./max.(s[2:end],eps(Float64)*first(s)))
    a,b=gaps(principal),gaps(full)
    if a==b && a>0 && min(principal[a],full[a])>1000max(principal[a+1],full[a+1])
        suggested=10max(principal[a+1],full[a+1],eps(Float64)*max(first(principal),first(full)))
        return "The common singular-value gap permits atol=$(round(suggested,RoundUp;sigdigits=3)) with rtol=0; inspect the residuals before accepting that cutoff."
    end
    return "No common rank gap is resolved; solve at higher order or with more accurate moments rather than flattening the matrix."
end

function _gns_cutoff(s,atol,rtol,check_atol)
    gap=length(s)<2 ? 0 : argmax(s[1:end-1]./max.(s[2:end],eps(Float64)*first(s)))
    retained=gap==0 ? first(s) : s[gap]
    discarded=gap==0 ? 0.0 : s[gap+1]
    resolved=gap>0 && retained>1000max(discarded,eps(Float64)*first(s)) && discarded<=check_atol*max(1.0,first(s))
    automatic=resolved ? min(retained/2,10max(discarded,eps(Float64)*first(s))) : 1e-8
    cutoff=max(isnothing(atol) ? automatic : atol,rtol*first(s))
    return cutoff,(retained=retained,discarded=discarded,resolved=resolved,automatic=isnothing(atol))
end

function _gns_candidate(result,domain,full_basis,basis,H,values,work,d,order;atol,rtol,check_atol,psd_rtol)
    W=eltype(basis)
    n=length(basis)
    A=H[1:n,1:n]
    F=LinearAlgebra.svd(A)
    # A valid Hankel is Hermitian. One eigensolve gives both its singular
    # values and PSD check; the Frobenius skew bound dominates the spectral
    # skew without another factorization.
    eigenvalues=LinearAlgebra.eigvals(LinearAlgebra.Hermitian((H+H')/2))
    full_s=sort!(abs.(eigenvalues);rev=true)
    cutoff,gap=_gns_cutoff(F.S,atol,rtol,check_atol)
    full_cutoff,full_gap=_gns_cutoff(full_s,atol,rtol,check_atol)
    r=count(>(cutoff),F.S); rf=count(>(full_cutoff),full_s)
    mineig=first(eigenvalues)
    psd_tolerance=psd_rtol*max(1.0,first(full_s))
    skew=LinearAlgebra.norm(H-H')
    U=F.U[:,1:r]; s=F.S[1:r]
    B=H[1:n,n+1:end]; C=H[n+1:end,n+1:end]
    Z=U*(LinearAlgebra.Diagonal(inv.(s))*(U'*B))
    projected=A*Z
    Cflat=Z'*projected
    range_error=LinearAlgebra.norm(B-projected)/(1+LinearAlgebra.norm(B)+LinearAlgebra.norm(projected))
    extension_error=LinearAlgebra.norm(C-Cflat)/(1+LinearAlgebra.norm(C)+LinearAlgebra.norm(Cflat))
    distance=isempty(C) ? 0.0 : LinearAlgebra.opnorm(C-Cflat)
    flat=r>0 && r==rf && range_error<=check_atol && extension_error<=check_atol
    psd=mineig>=-psd_tolerance && skew<=psd_tolerance
    letters=basis_rows(W,result.problem.objective.registry.algebra,1:result.problem.objective.registry.nsites,1)[2:end]
    operators=Dict{W,Matrix{ComplexF64}}()
    transform=U*LinearAlgebra.Diagonal(inv.(sqrt.(s)))
    full_index=Dict(w=>i for (i,w) in enumerate(full_basis))
    for l in letters
        K=zeros(ComplexF64,n,n)
        for j in 1:n
            empty!(work); mul!(work,l,basis[j],1.0+0im)
            for k in 1:length(work)
                col=full_index[work.words[k]]
                for i in 1:n
                    K[i,j]+=work.coefficients[k]*H[i,col]
                end
            end
        end
        operators[l]=transform'*K*transform
    end
    state=sqrt.(s).*conj.(U[1,:])
    left,right=similar(state),similar(state)
    observed(w)=_gns_observe!(left,right,operators,w,state)
    hankel_error=0.0
    for w in keys(values)
        target=_gns_value(result,_gns_coordinates(result.table,(w,),(1,)))
        isnothing(target) && continue
        hankel_error=max(hankel_error,abs(observed(w)-target))
    end
    program_error=0.0; program_count=0
    # Check complex moments and partially covered real/imaginary components.
    for w in domain
        target=_gns_value(result,_gns_coordinates(result.table,(w,),(1,)))
        if !isnothing(target)
            program_error=max(program_error,abs(observed(w)-target))
            program_count+=1
            continue
        end
        p,u=star(w)
        for (a,b) in ((1//2,scale(p,1//2+0im)),(-im//2,scale(p,im//2)))
            target=_gns_value(result,_gns_coordinates(result.table,(w,u),(a,b)))
            isnothing(target) && continue
            program_error=max(program_error,abs(a*observed(w)+b*observed(u)-target))
            program_count+=1
        end
    end
    objective=_gns_evaluate(operators,result.problem.objective,r,W)
    objective_error=abs(LinearAlgebra.dot(state,objective*state)-result.bound)
    constraint_min=Float64[]; equality_errors=Float64[]
    for c in result.problem.constraints
        M=_gns_evaluate(operators,c.polynomial,r,W)
        if c isa ExpectationConstraint
            value=real(LinearAlgebra.dot(state,M*state))
            push!(constraint_min,r==0 ? -Inf : value)
            push!(equality_errors,c isa ExpectationZero ? (r==0 ? Inf : abs(value)) : 0.0)
        else
            push!(constraint_min,r==0 ? -Inf : LinearAlgebra.eigmin(LinearAlgebra.Hermitian((M+M')/2)))
            push!(equality_errors,c isa Zero ? (r==0 ? Inf : LinearAlgebra.opnorm(M)) : 0.0)
        end
    end
    sector_errors=Float64[]
    for condition in result.relaxation.strengthen
        condition isa Sector || continue
        Q=_gns_evaluate(operators,condition.operator,r,W)
        push!(sector_errors,LinearAlgebra.norm(Q*state-condition.eigenvalue*state))
    end
    algebra_error=_gns_algebra_error(operators,letters,r,work)
    normalization_error=abs(LinearAlgebra.dot(state,state)-1)
    reproduced=max(hankel_error,program_error,objective_error,algebra_error,normalization_error,
        maximum(equality_errors;init=0.0),maximum(sector_errors;init=0.0),-minimum(constraint_min;init=0.0))<=check_atol
    optimizer=max(objective_error,algebra_error,normalization_error,
        maximum(equality_errors;init=0.0),maximum(sector_errors;init=0.0),-minimum(constraint_min;init=0.0))<=check_atol
    status=!psd ? :not_positive : !flat ? :not_flat : !reproduced ? :not_reproduced : :reproduced
    verdict=status===:reproduced ? :moments_reproduced : optimizer && (d<=order || psd) ? :optimizer_certified : :withheld
    advice=r==0 ? "No principal singular value exceeds the cutoff; lower atol below $(round(first(F.S),RoundDown;sigdigits=3)) and rtol below 1." :
        -mineig>1e-3*max(1.0,first(full_s)) ? "The full Hankel matrix has substantial negativity (minimum eigenvalue $(round(mineig;sigdigits=3))); solve at a higher relaxation order." :
        !psd ? "Full Hankel minimum eigenvalue is $(round(mineig;sigdigits=3)); inspect solver accuracy before increasing psd_rtol to $(round(max(-mineig,skew)/max(1.0,first(full_s)),RoundUp;sigdigits=3))." :
        !flat ? _gns_rank_advice(F.S,full_s)*" Measured extension residual=$(round(max(range_error,extension_error);sigdigits=3)); check_atol=$(round(check_atol;sigdigits=3))." :
        !reproduced ? "Measured residual requires check_atol ≥ $(round(max(hankel_error,program_error,objective_error,algebra_error,normalization_error,maximum(equality_errors;init=0.0),maximum(sector_errors;init=0.0),-minimum(constraint_min;init=0.0)),RoundUp;sigdigits=3)); inspect solver accuracy before accepting it." : ""
    flatness=(is_flat=flat,is_positive=psd,cutoff=cutoff,full_cutoff=full_cutoff,gap=gap,full_gap=full_gap,
        range_error=range_error,extension_error=extension_error,minimum_eigenvalue=mineig,
        psd_tolerance=psd_tolerance,check_atol=check_atol)
    checks=(hankel_moment_error=hankel_error,program_moment_error=program_error,program_components=program_count,
        objective_error=objective_error,constraint_min_eigenvalues=constraint_min,equality_errors=equality_errors,sector_errors=sector_errors,
        algebra_error=algebra_error,normalization_error=normalization_error,distance_to_flat=distance,
        sigma_min=isempty(s) ? 0.0 : last(s),sigma_max=first(F.S),
        condition_number=isempty(s) ? Inf : first(s)/last(s),distance_over_sigma_min=isempty(s) ? Inf : distance/last(s))
    return GNSResult(result.problem.objective.registry,verdict!==:withheld ? operators : nothing,
        verdict!==:withheld ? state : nothing,basis,full_basis,r,rf,collect(F.S),full_s,d-1,d,order,verdict,status,flatness,checks,advice)
end

function gns(result::SolveResult;atol=nothing,rtol=0.0,check_atol=1e-6,psd_rtol=1e-8)
    result.problem.objective.registry.algebra isa MixedAlgebra && throw(ArgumentError("Mixed GNS needs component-relation checks before reconstruction; use expectation(result, p)."))
    result.problem.objective.registry.algebra isa FunctionalAlgebra && throw(ArgumentError("State/trace-polynomial moments need scalar-factor decomposition before operator reconstruction; use expectation(result, p)."))
    result.problem.objective.registry.algebra isa Bosonic && throw(ArgumentError("No finite matrices satisfy the bosonic commutation relation; use expectation(result, p)."))
    all(x->isfinite(x)&&x>=0,(isnothing(atol) ? 0.0 : atol,rtol,check_atol,psd_rtol)) || throw(ArgumentError("GNS tolerances must be finite and nonnegative; supply nonnegative atol, rtol, check_atol and psd_rtol."))
    isnothing(result.bound) && throw(ArgumentError("GNS needs an accepted solve; resolve the problem before reconstructing."))
    relaxation=result.relaxation
    isnothing(relaxation.grouping) || throw(ArgumentError("GNS sparse clique reconstruction is not supported; solve with grouping=nothing and complete rows."))
    registry=result.problem.objective.registry
    domain=_gns_domain(result.table)
    W=eltype(domain)
    !isnothing(relaxation.family) && throw(ArgumentError("GNS sparse clique reconstruction is not supported for a custom family; solve with complete rows instead of family."))
    order=isnothing(relaxation.order) ? maximum(degree,relaxation.rows;init=0) : relaxation.order
    limit=registry.algebra isa Pauli ? registry.nsites+1 : registry.algebra isa Fermionic ? 2registry.nsites+1 : order
    values=Dict{W,ComplexF64}()
    partial=Set{W}()
    work=PolynomialBuffer{W,ComplexF64}()
    previous=nothing
    for d in 1:limit
        rows=basis_rows(W,registry.algebra,1:registry.nsites,d)
        basis=filter(w->degree(w)<d,rows)
        missing=_gns_required!(values,result,rows,work,partial)
        if !isnothing(missing)
            (isnothing(previous) || d<=order) && throw(ArgumentError("GNS needs moment $(_word_string(registry,missing)); $relaxation does not cover a full Hankel matrix. Remove split, custom rows or row_rule; solve with complete rows and a higher order."))
            return previous
        end
        H=_gns_fill!(zeros(ComplexF64,length(rows),length(rows)),rows,GNSProductValue(values,0.0im))
        previous=_gns_candidate(result,domain,rows,basis,H,values,work,d,order;atol,rtol,check_atol,psd_rtol)
        previous.flatness.is_flat && previous.flatness.is_positive && return previous
        d>order && !previous.flatness.is_positive && return previous
        length(rows)==length(basis) && return previous
    end
    isnothing(previous) && throw(ArgumentError("GNS needs positive-degree rows; solve with Relaxation(order=1) or higher."))
    return previous
end
