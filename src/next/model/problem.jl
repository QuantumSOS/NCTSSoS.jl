abstract type OperatorConstraint end
"Operator inequality p ⪰ 0, enforced by a localizing matrix (not just ⟨p⟩ ≥ 0)."
struct Positive{P<:Polynomial} <: OperatorConstraint
    polynomial::P
end
"Operator equality p = 0, enforced by zero localizing entries."
struct Zero{P<:Polynomial} <: OperatorConstraint
    polynomial::P
end

"Normalized states satisfying the input operator constraints; the plain route adds no state restriction."
struct ValidStates{C}
    constraints::C
end

function _checked_hermitian(p::Polynomial,name)
    all(isfinite(c) for (c,_) in terms(p)) || throw(ArgumentError("$name contains a nonfinite coefficient; replace NaN or Inf before constructing Problem."))
    LinearAlgebra.ishermitian(p) || throw(ArgumentError("$name is not Hermitian; supply a self-adjoint polynomial. If the Hermitian part is intended, explicitly use (p + p')/2 via 0.5*(p + p')."))
    return p
end
"Checked minimization problem. Constraints describe operator inequalities/equalities, not a state sector."
struct Problem{P<:Polynomial,C<:Tuple}
    objective::P
    constraints::C
    function Problem(objective::Polynomial;constraints=())
        _checked_hermitian(objective,"Objective")
        cs = Tuple(constraints)
        for (i,c) in enumerate(cs)
            c isa OperatorConstraint || throw(ArgumentError("Constraint $i must be Positive(p) or Zero(p); do not pass a Boolean comparison."))
            c.polynomial.registry === objective.registry || throw(ArgumentError("Constraint $i uses a different registry from the objective; rebuild it with the objective's variables."))
            _checked_hermitian(c.polynomial,"Constraint $i")
        end
        new{typeof(objective),typeof(cs)}(objective,cs)
    end
end
add_constraint(p::Problem,c::OperatorConstraint) = Problem(p.objective;constraints=(p.constraints...,c))

abstract type AbstractRelaxation end
"Plain relaxation. Specify exactly one of order or custom polynomial rows; custom rows are never filtered."
struct Relaxation{R} <: AbstractRelaxation
    order::Union{Nothing,Int}
    rows::R
    function Relaxation(;order=nothing,rows=nothing)
        xor(isnothing(order),isnothing(rows)) || throw(ArgumentError("Relaxation needs exactly one of order or rows; remove the other option."))
        if !isnothing(order)
            order isa Integer && order>=0 || throw(ArgumentError("Relaxation order must be a nonnegative integer; use order=1 or a larger integer."))
            return new{Nothing}(Int(order),nothing)
        end
        rs = collect(rows)
        !isempty(rs) || throw(ArgumentError("Custom rows are empty; include the identity and rows that cover the objective."))
        all(r->r isa Polynomial,rs) || throw(ArgumentError("Custom rows must be polynomials from the problem registry; write one(X[1]) for the identity."))
        all(r->length(r.data)==1 && r.data.coefficients[1]==1,rs) || throw(ArgumentError("Custom rows must each be a single normal word with coefficient 1; supply such rows explicitly."))
        allunique(rs[i].data.words[1] for i in eachindex(rs)) || throw(ArgumentError("Custom rows contain duplicates; remove duplicate rows explicitly."))
        any(r->isone(r.data.words[1]),rs) || throw(ArgumentError("Custom rows lack the identity; include one(variable)."))
        new{typeof(rs)}(nothing,rs)
    end
end

_pack_word(::Type{W},w::AbstractWord) where {W<:AbstractWord} = convert(W,w)
_pack_word(::Type{SiteWord{A,V,B,S}},w::SiteWord{A,U,B,S}) where {A,U,V,B,S} =
    SiteWord{A,V,B,S}(V(w.bits),Val(:normal))

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

struct CoveredMomentSink{S}
    sink::S
    name::String
end
function push_term!(s::CoveredMomentSink,c::Number,w)
    iszero(c) && return s
    key,_,_ = moment_key(w)
    moment_ids(s.sink.table,key)!=(0,0) || throw(ArgumentError("$(s.name) needs word $w which the custom rows do not cover; use complete rows through degree $(cld(degree(w),2)) or add rows covering this word."))
    push_term!(s.sink,c,w)
    return s
end
struct RightProductSink{S,W}
    sink::S
    right::W
end
push_term!(s::RightProductSink,c::Number,w) = mul!(s.sink,w,s.right,c)

function _localize!(program,table,c::OperatorConstraint,order,index)
    g = c.polynomial
    local_order = order-cld(degree(g),2)
    local_order>=0 || throw(ArgumentError("Constraint $index has degree $(degree(g)), above the relaxation capacity; use order ≥ $(cld(degree(g),2))."))
    W = eltype(table.keys)
    r = g.registry
    rows = basis_rows(W,r.algebra,1:r.nsites,local_order)
    sink = CoveredMomentSink(MomentSink(table,program.forms),"Constraint $index")
    n = length(rows)
    for j in 1:n,i in j:n
        phase,left = star(rows[i])
        right = RightProductSink(sink,rows[j])
        for (coef,w) in terms(g); mul!(right,left,_pack_word(W,w),scale(phase,coef)); end
        close_form!(program.forms)
    end
    push_block!(program,c isa Positive ? HPSD : ZERO,c isa Positive ? n : n*(n+1)÷2,index)
    return n
end

"Build without a solver. The returned program and table are suitable for pipeline/conic inspection."
function build(p::Problem,r::Relaxation)
    # Recheck mutable polynomial storage rather than trusting stale construction checks.
    Problem(p.objective;constraints=p.constraints)
    order = isnothing(r.order) ? maximum(degree,r.rows) : r.order
    max_letters = max(2order,degree(p.objective),maximum((degree(c.polynomial) for c in p.constraints);init=0))
    W = word_type(p.objective.registry,max_letters)
    rows = _rows(p,r,W)
    program,table = plain_program(rows,((c,_pack_word(W,w)) for (c,w) in terms(p.objective)))
    local_sizes = Int[_localize!(program,table,c,order,i) for (i,c) in enumerate(p.constraints)]
    program.nvars = nvars(table)
    validate(program)
    sizes = (rows=length(rows),moments=length(table),real_variables=nvars(program)-1,
        blocks=length(program.blocks),largest_block=maximum(b.side for b in program.blocks if b.kind in (PSD,HPSD)),
        localizing_rows=local_sizes)
    return (program=program,table=table,sizes=sizes,valid_for=ValidStates(p.constraints))
end
inspect(p::Problem,r::AbstractRelaxation) = build(p,r).sizes
