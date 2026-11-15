abstract type OperatorConstraint end
"Operator inequality p ⪰ 0, enforced by a localizing matrix (not just ⟨p⟩ ≥ 0)."
struct Positive{P<:Polynomial} <: OperatorConstraint
    polynomial::P
end
"Operator equality p = 0, enforced by zero localizing entries."
struct Zero{P<:Polynomial} <: OperatorConstraint
    polynomial::P
end

abstract type ExpectationConstraint end
"""
    ExpectationPositive(p)

The single scalar condition ⟨p⟩ ≥ 0 on a Hermitian polynomial, without a
localizing matrix. For example, `Problem(A; constraints=(ExpectationPositive(E-H),))`
bounds A subject to an energy expectation at most E. The selected rows must
already cover p. Functional state/trace polynomials are not supported here.
"""
struct ExpectationPositive{P<:Polynomial} <: ExpectationConstraint
    polynomial::P
    ExpectationPositive(p::Polynomial) = new{typeof(p)}(_checked_expectation(p,"ExpectationPositive"))
end
"""
    ExpectationZero(p)

The single scalar condition ⟨p⟩ = 0 on a Hermitian polynomial. Write
`ExpectationZero(p-c)` for ⟨p⟩ = c. Unlike Zero(p), it does not constrain
all localized products. The selected rows must already cover p.
"""
struct ExpectationZero{P<:Polynomial} <: ExpectationConstraint
    polynomial::P
    ExpectationZero(p::Polynomial) = new{typeof(p)}(_checked_expectation(p,"ExpectationZero"))
end
Base.show(io::IO,c::ExpectationConstraint) =
    print(io,"⟨",c.polynomial,c isa ExpectationPositive ? "⟩ ≥ 0 (scalar)" : "⟩ = 0 (scalar)")
Base.show(io::IO,::MIME"text/plain",c::ExpectationConstraint) = show(io,c)

"Normalized states satisfying the input operator and scalar expectation constraints."
struct ValidStates{C}
    constraints::C
end

function Base.show(io::IO,v::ValidStates)
    print(io,"Normalized states")
    isempty(v.constraints) && return print(io,"; no additional state conditions")
    print(io," satisfying: ")
    for (i,c) in enumerate(v.constraints)
        i>1 && print(io,"; ")
        if c isa ExpectationConstraint
            show(io,c)
        else
            print(io,c.polynomial,c isa Positive ? " ⪰ 0 (operator)" : " = 0 (operator)")
        end
    end
end
Base.show(io::IO,::MIME"text/plain",v::ValidStates) = show(io,v)

function _checked_hermitian(p::Polynomial,name)
    all(isfinite(c) for (c,_) in terms(p)) || throw(ArgumentError("$name contains a nonfinite coefficient; replace NaN or Inf before constructing Problem."))
    LinearAlgebra.ishermitian(p) || throw(ArgumentError("$name ($p) is not Hermitian; supply a self-adjoint polynomial. If the Hermitian part is intended, explicitly use (p + p')/2 via 0.5*(p + p')."))
    return p
end
function _checked_expectation(p::Polynomial,name)
    p.registry.algebra isa FunctionalAlgebra && throw(ArgumentError(
        "$name is a single outer moment constraint for ordinary algebras. For state/trace variables, use Positive(state_symbol(p)) / Zero(state_symbol(p)) (or Positive(p) / Zero(p) if p is already scalar); outer averaged functional conditions are not supported."))
    return _checked_hermitian(p,name)
end
"""
Checked minimization problem with operator or scalar expectation constraints and optional symmetry.
Pass a tuple or vector: derivations constrain continuous symmetry, antiunitary
actions constrain moments only, and other actions close the block family.
Use `Next.Symmetry` to explicitly place a unitary action in `moment_only`.
Declare symmetry here, not on `Relaxation`. Empty symmetries mean none declared.
For a legacy `BlockRelaxation`, use whichever declaration is nonempty; if both
are nonempty, remove the legacy declaration and keep the one on the problem.
Positive/Zero describe operator constraints; ExpectationPositive/ExpectationZero
describe single scalar conditions. Neither is a state-sector restriction.
"""
struct Problem{P<:Polynomial,C<:Tuple,S}
    objective::P
    constraints::C
    symmetry::S
    function Problem(objective::Polynomial;constraints=(),symmetry=nothing)
        _checked_hermitian(objective,"Objective")
        cs = Tuple(constraints)
        for (i,c) in enumerate(cs)
            c isa Union{OperatorConstraint,ExpectationConstraint} || throw(ArgumentError("Constraint $i must be Positive(p), Zero(p), ExpectationPositive(p) or ExpectationZero(p); do not pass a Boolean comparison."))
            c.polynomial.registry === objective.registry || throw(ArgumentError("Constraint $i uses a different registry from the objective; rebuild it with the objective's variables."))
            _checked_hermitian(c.polynomial,"Constraint $i")
        end
        isnothing(symmetry) || (symmetry=_check_symmetry_domain(objective,symmetry))
        new{typeof(objective),typeof(cs),typeof(symmetry)}(objective,cs,symmetry)
    end
end
add_constraint(p::Problem,c::Union{OperatorConstraint,ExpectationConstraint}) = Problem(p.objective;constraints=(p.constraints...,c),symmetry=p.symmetry)

abstract type AbstractRelaxation end
"""
    Relaxation(; order, rows, row_rule, grouping, split, family, strengthen=())

Choose exactly one of `order`, custom `rows`, or a custom block `family`.
Declare symmetry on `Problem`. Without symmetry or block options, use complete
moment/localizing matrices. Otherwise use the selected block family.
Custom rows can be split with `split=TermSplit(...)`. Before coverage checks
and splitting, these starting rows are closed under declared block generators
and continuous derivations. Localizing order uses the largest closed row degree.
Moment-only maps do not add rows. The selected
block family is then closed as well. Without a split, custom rows are unchanged.
Custom families are used as supplied. Remove grouping and
row_rule with custom rows, and all three options with custom families.
A single strengthening kind is accepted without a tuple.
"""
struct Relaxation{R,U,G,S,F,K} <: AbstractRelaxation
    order::Union{Nothing,Int}
    rows::R
    row_rule::U
    grouping::G
    split::S
    family::F
    strengthen::K
    function Relaxation(;order=nothing,rows=nothing,row_rule=nothing,grouping=nothing,
        split=nothing,family=nothing,strengthen=(),symmetry=nothing)
        isnothing(symmetry) || throw(ArgumentError("Declare symmetry on Problem(H; symmetry=(generators...,)), not on Relaxation; remove the symmetry keyword here."))
        count(!isnothing,(order,rows,family))==1 || throw(ArgumentError("Relaxation needs exactly one of order, rows, or family; remove the other options."))
        if !isnothing(rows)
            all(isnothing,(row_rule,grouping)) && (isnothing(split) || split isa TermSplit) || throw(ArgumentError("Custom rows cannot be combined with grouping, row_rule or PhysicsInformed; use rows with split=TermSplit(...) or remove these options."))
        elseif !isnothing(family)
            all(isnothing,(row_rule,grouping,split)) || throw(ArgumentError("A custom family is used as supplied; remove row_rule, grouping and split, or use rows with split=TermSplit(...)."))
        end
        isnothing(row_rule) || row_rule isa Union{Contiguous,Function} || throw(ArgumentError("Use row_rule=Contiguous(...) or a row-rule function, or remove row_rule."))
        isnothing(grouping) || grouping isa Union{Correlative,PhysicsInformed} || throw(ArgumentError("Use grouping=Correlative() or PhysicsInformed(...), or remove grouping."))
        isnothing(split) || split isa Union{TermSplit,PhysicsInformed} || throw(ArgumentError("Use split=TermSplit() or PhysicsInformed(...), or remove split."))
        strengthen isa StrengtheningKind && (strengthen=(strengthen,))
        strengthen isa Union{Tuple,AbstractVector} || throw(ArgumentError("Use strengthen=(DensityMatrix(...), Stationarity(), ...) or a single strengthening kind; use () for none."))
        kinds=Tuple(strengthen)
        all(k->k isa StrengtheningKind,kinds) || throw(ArgumentError("Use strengthening kinds such as DensityMatrix(...), Stationarity(), or Sector(Q,q), not types or arbitrary values."))
        if !isnothing(order)
            order isa Integer && order>=0 || throw(ArgumentError("Relaxation order must be a nonnegative integer; use order=1 or a larger integer."))
            order=Int(order)
        end
        rs=nothing
        if !isnothing(rows)
            rs = collect(rows)
            !isempty(rs) || throw(ArgumentError("Custom rows are empty; include the identity and rows that cover the objective."))
            all(r->r isa Polynomial,rs) || throw(ArgumentError("Custom rows must be polynomials from the problem registry; write one(X[1]) for the identity."))
            all(r->length(r.data)==1 && r.data.coefficients[1]==1,rs) || throw(ArgumentError("Custom rows must each be a single normal word with coefficient 1; supply such rows explicitly."))
            allunique(rs[i].data.words[1] for i in eachindex(rs)) || throw(ArgumentError("Custom rows contain duplicates; remove duplicate rows explicitly."))
            any(r->isone(r.data.words[1]),rs) || throw(ArgumentError("Custom rows lack the identity; include one(variable)."))
        end
        if !isnothing(family)
            family=collect.(family)
            !isempty(family) && all(!isempty,family) || throw(ArgumentError("The custom family and each block must be nonempty; supply rows in each block."))
        end
        new{typeof(rs),typeof(row_rule),typeof(grouping),typeof(split),typeof(family),typeof(kinds)}(
            order,rs,row_rule,grouping,split,family,kinds)
    end
end

function _rows(p::Problem,r::Relaxation,::Type{W}) where {W}
    registry = p.objective.registry
    if isnothing(r.rows)
        return basis_rows(W,registry.algebra,1:registry.nsites,r.order)
    end
    for (i,row) in enumerate(r.rows)
        row.registry === registry || throw(ArgumentError("Custom row $i uses a different registry; reconstruct it with the objective's variables."))
    end
    return W[_pack_word(W,row.data.words[1]) for row in r.rows]
end

struct CoveredMomentSink{S,R}
    sink::S
    name::String
    registry::R
end
function push_term!(s::CoveredMomentSink,c::Number,w)
    iszero(c) && return s
    key,_,_ = moment_key(w)
    moment_ids(s.sink.table,key)!=(0,0) || throw(ArgumentError("$(s.name) needs word $(_word_string(s.registry,w)) which the custom rows do not cover; use complete rows through degree $(cld(degree(w),2)) or add rows covering this word."))
    push_term!(s.sink,c,w)
    return s
end
struct RightProductSink{S,W}
    sink::S
    right::W
end
push_term!(s::RightProductSink,c::Number,w) = mul!(s.sink,w,s.right,c)

_localizing_order(order,g::Polynomial) = order-cld(degree(g),2)

function _localize!(program::MomentProgram{C},table,c::OperatorConstraint,order,index) where {C}
    g = c.polynomial
    local_order = _localizing_order(order,g)
    local_order>=0 || throw(ArgumentError("Constraint $index ($g) has degree $(degree(g)), above the relaxation capacity; use order ≥ $(cld(degree(g),2))."))
    W = eltype(table.keys)
    r = g.registry
    rows = basis_rows(W,r.algebra,1:r.nsites,local_order)
    sink = CoveredMomentSink(MomentSink(table,program.forms),"Constraint $index",r)
    n = length(rows)
    for j in 1:n,i in j:n
        phase,left = star(rows[i])
        right = RightProductSink(sink,rows[j])
        for (coef,w) in terms(g); mul!(right,left,_pack_word(W,w),scale(phase,coef)); end
        if i==j && C <: Complex
            # g is checked Hermitian: r'*g*r is real. Its separate terms can
            # have cancelling imaginary parts; store the real diagonal directly.
            for t in program.forms.offsets[end]:length(program.forms.coefs)
                program.forms.coefs[t] = real(program.forms.coefs[t])
            end
        end
        close_form!(program.forms)
    end
    kind = c isa Positive ? (C <: Real ? PSD : HPSD) : ZERO
    push_block!(program,kind,c isa Positive ? n : n*(n+1)÷2,index)
    return n
end

function _expectation_constraint!(program,table,constraint::ExpectationConstraint,index)
    p=constraint.polynomial
    W=eltype(p.data.words)
    coordinates=Dict{Int,Tuple{Number,W}}()
    missing(w)=throw(ArgumentError("$(nameof(typeof(constraint))) constraint $index needs uncovered word $(_word_string(p.registry,w)); add rows covering this word or increase the order."))
    for (c,w) in terms(p)
        pairs=moment_terms(table,w)
        isnothing(pairs) && missing(w)
        for (id,factor) in pairs
            coefficient=c*factor
            coordinates[id]=(get(coordinates,id,(0,w))[1]+coefficient,w)
        end
    end
    # Combine before coverage: adjoint terms can cancel an unused imaginary part.
    for id in sort!(collect(keys(coordinates)))
        coefficient,w=coordinates[id]
        iszero(coefficient) && continue
        moment_covered(table,id) || missing(w)
        push_term!(program.forms,real(coefficient),id)
    end
    close_form!(program.forms)
    push_block!(program,constraint isa ExpectationPositive ? NONNEGATIVE : ZERO,1,index)
    return program
end

"Flat program counts on result.sizes; result.sizes.closure is nothing on the plain route."
struct ProgramSizes
    rows::Int
    moments::Int
    real_variables::Int
    blocks::Int
    largest_block::Int
    localizing_rows::Vector{Int}
end
Base.getproperty(s::ProgramSizes,name::Symbol) = name===:closure ? nothing : getfield(s,name)
Base.propertynames(::ProgramSizes,private::Bool=false) = (fieldnames(ProgramSizes)...,:closure)
Base.:(==)(a::ProgramSizes,b::ProgramSizes) = all(getfield(a,k)==getfield(b,k) for k in fieldnames(ProgramSizes))
Base.hash(s::ProgramSizes,h::UInt) = hash((s.rows,s.moments,s.real_variables,s.blocks,s.largest_block,s.localizing_rows),h)
function Base.show(io::IO,s::ProgramSizes)
    print(io,s.rows," rows, ",s.moments," moments, ",s.real_variables,
        " real variables, ",s.blocks," blocks (largest ",s.largest_block,
        "); localizing rows: ",s.localizing_rows)
end
Base.show(io::IO,::MIME"text/plain",s::ProgramSizes) = show(io,s)

# These algebras have real structure constants. With real objective and
# constraint coefficients, averaging a functional with its conjugate preserves
# positivity, normalization and the objective, so a real optimum suffices.
# FunctionalAlgebra instead declares real-valued symbols explicitly; no
# averaging argument for products of arbitrary complex expectations is used.
_real_form(p::Problem) = _real_structure(p.objective.registry.algebra) &&
    all(isreal(c) for (c,_) in terms(p.objective)) &&
    all(c->all(isreal(a) for (a,_) in terms(c.polynomial)),p.constraints)

"A built plain or block relaxation, with lazily cached solver data. Treat its program and inputs as read-only; rebuild after changing them."
struct BuiltRelaxation{P,R,G,S,V}
    program::G
    table::Any
    sizes::S
    valid_for::V
    problem::P
    relaxation::R
    cache::Dict{Any,Any}
end
BuiltRelaxation(p,r,program,table,sizes,valid_for) =
    BuiltRelaxation(program,table,sizes,valid_for,p,r,Dict{Any,Any}())
function Base.show(io::IO,b::BuiltRelaxation)
    print(io,"Built relaxation\n",b.problem,"\n",b.relaxation,"\nSizes: ",b.sizes)
end
Base.show(io::IO,::MIME"text/plain",b::BuiltRelaxation) = show(io,b)

"""
    Next.build(problem, relaxation) -> Next.BuiltRelaxation

Formulate the plain or block relaxation without lowering or loading a solver.
Use `Next.prepare(built; solver)` to lower and prepare, then `solve(built; solver)`
to optimize, or omit prepare to do both on the first solve. The program, table,
sizes and valid_for fields remain available for inspection.
"""
function build(p::Problem,r::Relaxation)
    if p.objective.registry.algebra isa MixedAlgebra && _uses_blocks(p,r)
        throw(ArgumentError("Mixed algebras support the plain route only; remove symmetry, grouping, split, row_rule, family and strengthening, and use order or custom rows."))
    end
    if p.objective.registry.algebra isa FunctionalAlgebra && _uses_blocks(p,r)
        throw(ArgumentError("State/trace polynomials currently support the plain dense route only; remove symmetry, grouping, split, row_rule, family and strengthening, and use Relaxation(order=...)."))
    end
    _uses_blocks(p,r) && return _build_blocks(p,r)
    # Recheck mutable polynomial storage rather than trusting stale construction checks.
    Problem(p.objective;constraints=p.constraints)
    order = isnothing(r.order) ? maximum(degree,r.rows) : r.order
    max_letters = max(2order,degree(p.objective),maximum((degree(c.polynomial) for c in p.constraints);init=0))
    W = word_type(p.objective.registry,max_letters)
    rows = _rows(p,r,W)
    program,table = plain_program(rows,((c,_pack_word(W,w)) for (c,w) in terms(p.objective));
        C=_real_form(p) ? Float64 : ComplexF64,
        word_string=w->_word_string(p.objective.registry,w))
    local_sizes = Int[]
    for (i,c) in enumerate(p.constraints)
        if c isa ExpectationConstraint
            _expectation_constraint!(program,table,c,i)
        else
            push!(local_sizes,_localize!(program,table,c,order,i))
        end
    end
    program.nvars = nvars(table)
    validate(program)
    sizes = ProgramSizes(length(rows),length(table),nvars(program)-1,
        length(program.blocks),maximum(b.side for b in program.blocks if b.kind in (PSD,HPSD)),
        local_sizes)
    return BuiltRelaxation(p,r,program,table,sizes,ValidStates(p.constraints))
end
