"Declare block generators and optional moment-only actions on `Problem`."
Symmetry(generators...;moment_only=(),continuous=()) = FamilySymmetry(collect(generators);moment_only,continuous)

_has_symmetry(::Nothing) = false
_has_symmetry(s::FamilySymmetry) = !isempty(s.generators) || !isempty(s.moment_only) || !isempty(s.continuous)
_uses_blocks(p::Problem,r::Relaxation) = _has_symmetry(p.symmetry) ||
    !all(isnothing,(r.row_rule,r.grouping,r.split,r.family)) || !isempty(r.strengthen)

function _check_symmetry_domain(p::Polynomial,s)
    if s isa Union{Tuple,AbstractVector}
        items=Tuple(s)
        for (i,g) in enumerate(items)
            g isa Union{AbstractAction,AbstractDerivation} || throw(ArgumentError(
                "Symmetry item $i must be an action or derivation, such as SitePermutation(...) or SpinRotation(...)."))
        end
        block=filter(g->g isa AbstractAction && !is_antiunitary(g),items)
        moments=filter(g->g isa AbstractAction && is_antiunitary(g),items)
        continuous=filter(g->g isa AbstractDerivation,items)
        s=Symmetry(block...;moment_only=moments,continuous)
    end
    s isa FamilySymmetry || throw(ArgumentError(
        "Use Problem(H; symmetry=(generators...,)) or Next.Symmetry(...; moment_only=(...), continuous=(...))."))
    p.registry.algebra isa MixedAlgebra && _has_symmetry(s) && throw(ArgumentError(
        "Mixed symmetry needs component-aware actions; remove symmetry and use the plain relaxation."))
    p.registry.algebra isa FunctionalAlgebra && _has_symmetry(s) && throw(ArgumentError(
        "State/trace polynomials currently support the plain dense route only; remove symmetry until functional block routing is available."))
    W=eltype(p.data.words)
    n=p.registry.nsites
    for (i,g) in enumerate(Iterators.flatten((s.generators,s.moment_only)))
        applicable(act,g,one(W)) || throw(ArgumentError("Symmetry generator $i has an incompatible algebra or word type; construct its images with this problem's variables and word type $W. Use SitePermutation for whole-site maps, LetterPermutation for Free/Projector/Unipotent letters, or CliffordMap/AntiunitaryPauli for Pauli variables."))
        lastsite = g isa SitePermutation ? length(g.images) :
            g isa LetterPermutation ? _action_sites(g) :
            g isa CliffordMap ? length(g.columns)÷2 :
            g isa PauliConjugation ? maximum(sites(g.word);init=0) : 0
        lastsite<=n || throw(ArgumentError("Symmetry generator $i acts outside sites 1:$n; rebuild it for this problem's sites."))
        if g isa LetterPermutation
            algebra=p.registry.algebra isa FunctionalAlgebra ? p.registry.algebra.base : p.registry.algebra
            all(s->_label_count(g,s)<=algebra.nletters,1:n) || throw(ArgumentError(
                "LetterPermutation generator $i uses settings outside this problem's registry; rebuild it from this problem's variables."))
            (isempty(g.names) || g.names===p.registry.names) || throw(ArgumentError(
                "LetterPermutation generator $i ($g) uses a different registry from the objective; rebuild it from this problem's variables."))
        end
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
    rows=r.split isa TermSplit ? r.rows : nothing
    family=isnothing(r.rows) || !isnothing(rows) ? r.family : [r.rows]
    description=BlockRelaxation(;order=r.order,row_rule=r.row_rule,grouping=r.grouping,
        split=r.split,rows,family,strengthen=r.strengthen)
    return _build_declared(p,description;relaxation=r)
end

# The pipeline method remains unchanged. This more specific model method merges
# declarations, checks named polynomials, then invokes that implementation.
build(p::Problem{P,C,Nothing},r::BlockRelaxation) where {P<:Polynomial,C<:Tuple} = _build_declared(p,r)
build(p::Problem{P,C,S},r::BlockRelaxation) where {P<:Polynomial,C<:Tuple,S<:FamilySymmetry} = _build_declared(p,r)
function _build_declared(p::Problem,r::BlockRelaxation;relaxation=r)
    p.objective.registry.algebra isa FunctionalAlgebra && throw(ArgumentError("State/trace polynomials do not yet support the block pipeline; use Relaxation(order=...) without symmetry or strengthening."))
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
    any(c->c isa ExpectationConstraint,p.constraints) && throw(ArgumentError(
        "Scalar expectation constraints do not yet support the block pipeline; use a plain Relaxation(order=...) or rows=... without symmetry, grouping, split or strengthening until scalar graph seeding and block integration are available."))
    description = symmetry === r.symmetry ? r : BlockRelaxation(;order=r.order,
        row_rule=r.row_rule,grouping=r.grouping,split=r.split,rows=r.rows,family=r.family,
        symmetry,strengthen=r.strengthen)
    built=invoke(build,Tuple{Problem,BlockRelaxation},p,description)
    return BuiltRelaxation(p,relaxation,built.program,built.table,built.sizes,built.valid_for)
end

function Base.show(io::IO,p::Problem)
    print(io,"Problem: minimize ",p.objective,"; ",count(c->c isa OperatorConstraint,p.constraints)," operator constraints")
    scalar=count(c->c isa ExpectationConstraint,p.constraints)
    scalar>0 && print(io,", ",scalar," scalar expectation constraints")
    if p.objective.registry.algebra isa FunctionalAlgebra
        print(io,"; domain: ",p.objective.registry.algebra isa FunctionalAlgebra{RealTrace} ? "normalized real traces" : "normalized real states")
    end
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
"""
    report = inspect(problem, relaxation)

Build the relaxation and return a report without calling a solver.
The report has `problem`, `relaxation` and `sizes` fields.
Print `report` to see the selected options and program sizes.
See [`solve`](@ref) for the fields of `sizes`.
The build uses the algebra relations of the problem's variables.
It does not choose a solver or compute a bound.

Example: `using NCTSSoS.Next; registry, (x, y, z) = pauli_variables(1); report = inspect(Problem(z[1]), Relaxation(order=1))`.
"""
inspect(p::Problem,r::AbstractRelaxation) = Inspection(p,r,build(p,r).sizes)
function Base.show(io::IO,r::Inspection)
    print(io,r.problem,"\n",r.relaxation,"\nSizes: ",r.sizes)
end
Base.show(io::IO,::MIME"text/plain",r::Inspection) = show(io,r)
