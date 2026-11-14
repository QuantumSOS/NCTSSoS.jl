"Declare block generators and optional moment-only actions on `Problem`."
Symmetry(generators...;moment_only=()) = FamilySymmetry(collect(generators);moment_only)

_has_symmetry(::Nothing) = false
_has_symmetry(s::FamilySymmetry) = !isempty(s.generators) || !isempty(s.moment_only)
_uses_blocks(p::Problem,r::Relaxation) = _has_symmetry(p.symmetry) ||
    !all(isnothing,(r.row_rule,r.grouping,r.split,r.family)) || !isempty(r.strengthen)

function _check_symmetry_domain(p::Polynomial,s)
    s isa FamilySymmetry || throw(ArgumentError("Use Problem(H; symmetry=Symmetry(generators...; moment_only=()))."))
    W=eltype(p.data.words)
    n=p.registry.nsites
    for (i,g) in enumerate(Iterators.flatten((s.generators,s.moment_only)))
        applicable(act,g,one(W)) || throw(ArgumentError("Symmetry generator $i has an incompatible algebra or word type; construct its images with this problem's variables and word type $W."))
        lastsite = g isa SitePermutation ? length(g.images) :
            g isa CliffordMap ? length(g.columns)÷2 :
            g isa PauliConjugation ? maximum(sites(g.word);init=0) : 0
        lastsite<=n || throw(ArgumentError("Symmetry generator $i acts outside sites 1:$n; rebuild it for this problem's sites."))
    end
    return s
end

function _check_named_invariance(p::Polynomial,generators,label)
    objective=collect(terms(p))
    try
        _check_objective(objective,generators)
    catch error
        error isa ArgumentError || rethrow()
        for (i,g) in enumerate(generators)
            try
                _check_objective(objective,(g,))
            catch failure
                failure isa ArgumentError || rethrow()
                throw(ArgumentError("Symmetry generator $i does not preserve $label ($p); remove this generator or supply a polynomial it preserves."))
            end
        end
        rethrow()
    end
end

function _build_blocks(p::Problem,r::Relaxation)
    family=isnothing(r.rows) ? r.family : [r.rows]
    description=BlockRelaxation(;order=r.order,row_rule=r.row_rule,grouping=r.grouping,
        split=r.split,family,strengthen=r.strengthen)
    return build(p,description)
end

# The pipeline method remains unchanged. This more specific model method merges
# declarations, checks named polynomials, then invokes that implementation.
build(p::Problem{P,C,Nothing},r::BlockRelaxation) where {P<:Polynomial,C<:Tuple} = _build_declared(p,r)
build(p::Problem{P,C,S},r::BlockRelaxation) where {P<:Polynomial,C<:Tuple,S<:FamilySymmetry} = _build_declared(p,r)
function _build_declared(p::Problem,r::BlockRelaxation)
    a,b=p.symmetry,r.symmetry
    _has_symmetry(a) && _has_symmetry(b) && throw(ArgumentError(
        "Symmetry is declared on both Problem and BlockRelaxation; keep it on Problem and remove the BlockRelaxation symmetry."))
    symmetry=_has_symmetry(a) ? a : b
    _check_symmetry_domain(p.objective,symmetry)
    if _has_symmetry(symmetry)
        generators=(symmetry.generators...,symmetry.moment_only...)
        _check_named_invariance(p.objective,generators,"objective")
        for (i,c) in enumerate(p.constraints)
            _check_named_invariance(c.polynomial,generators,"constraint $i")
        end
    end
    description = symmetry === r.symmetry ? r : BlockRelaxation(;order=r.order,
        row_rule=r.row_rule,grouping=r.grouping,split=r.split,family=r.family,
        symmetry,strengthen=r.strengthen)
    return invoke(build,Tuple{Problem,BlockRelaxation},p,description)
end

function Base.show(io::IO,p::Problem)
    print(io,"Problem: minimize ",p.objective,"; ",length(p.constraints)," operator constraints")
    _has_symmetry(p.symmetry) && print(io,"; ",p.symmetry)
end
Base.show(io::IO,::MIME"text/plain",p::Problem) = show(io,p)

function Base.show(io::IO,r::Union{Relaxation,BlockRelaxation})
    print(io,"Relaxation(")
    separator=""
    for key in (:order,:rows,:row_rule,:grouping,:split,:family,:strengthen,:symmetry)
        hasfield(typeof(r),key) || continue
        value=getfield(r,key)
        (isnothing(value) || key===:strengthen && isempty(value)) && continue
        key===:symmetry && !_has_symmetry(value) && continue
        print(io,separator,key,"=",value)
        separator=", "
    end
    print(io,")")
end
Base.show(io::IO,::MIME"text/plain",r::Union{Relaxation,BlockRelaxation}) = show(io,r)

"Inspection of the selected options and the resulting program sizes."
struct Inspection{P,R,S}
    problem::P
    relaxation::R
    sizes::S
end
inspect(p::Problem,r::AbstractRelaxation) = Inspection(p,r,build(p,r).sizes)
function Base.show(io::IO,r::Inspection)
    print(io,r.problem,"\n",r.relaxation,"\nSizes: ",r.sizes)
end
Base.show(io::IO,::MIME"text/plain",r::Inspection) = show(io,r)
