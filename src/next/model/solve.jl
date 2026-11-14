struct SolveResult{P,R,S,B,V,Q}
    bound::Union{Nothing,Float64}
    status::S
    sizes::B
    valid_for::V
    problem::P
    relaxation::R
    solution::Q
    # Width is chosen at build time. Do not make solve's return type depend on it.
    table::Any
end

# Reuse the redundant tag's storage for the already-built moment table.
@inline function Base.getproperty(r::SolveResult,name::Symbol)
    name===:bound_kind && return isnothing(getfield(r,:bound)) ? :unavailable : :numerical_optimum_not_certificate
    return getfield(r,name)
end
Base.propertynames(::SolveResult,private::Bool=false) =
    (:bound,:bound_kind,:status,:sizes,:valid_for,:problem,:relaxation,:solution,:table)

function Base.show(io::IO,r::SolveResult)
    print(io,"Bound: ",r.bound," (",r.bound_kind,"); status: ",r.status,
        "\nSizes: ",r.sizes,"\nValid for: ",r.valid_for)
end
Base.show(io::IO,::MIME"text/plain",r::SolveResult) = show(io,r)

"One checked build and solve. The backend receives only standard data, not model or algebra objects."
function solve(problem::Problem,relaxation::AbstractRelaxation;solver=nothing)
    solver isa AbstractSolver || throw(ArgumentError("Missing or invalid solver description; load SCS with `using SCS` and pass solver=SCSBackend(), or pass solver=JuMPBackend(optimizer_factory), with attributes on the factory if needed."))
    built = build(problem,relaxation)
    solution = solve_standard(moment_standard_form(built.program),solver)
    return SolveResult(solution.bound,solution.status,built.sizes,built.valid_for,problem,relaxation,solution,built.table)
end
"Add an operator constraint and solve again on the same open moment-form handle."
function solve(previous::SolveResult,constraint::OperatorConstraint)
    problem = add_constraint(previous.problem,constraint)
    built = build(problem,previous.relaxation)
    solution = solve_standard(moment_standard_form(built.program),previous.solution.handle)
    return SolveResult(solution.bound,solution.status,built.sizes,built.valid_for,problem,previous.relaxation,solution,built.table)
end

public expectation
"""
    expectation(result::SolveResult, p::Polynomial)

Evaluate ⟨p⟩ in the numerical relaxation solution, not a reconstructed state or
certificate. `p` must use the problem's registry and only represented moments.
Return a real number for Hermitian `p`, and a complex number otherwise.
"""
function expectation(result::SolveResult,p::Polynomial)
    isnothing(result.bound) && throw(ArgumentError("Cannot read expectations: solve status $(result.status) has no accepted solution."))
    p.registry === result.problem.objective.registry ||
        throw(ArgumentError("Expectation polynomial uses a different registry; construct it from the problem's variables."))
    return _expectation(result.table,result,p)
end

# Specialize the evaluation loop on the table's packed word type only on demand.
function _expectation(table,result::SolveResult,p::Polynomial)
    value = 0.0+0.0im
    for (c,w) in terms(p)
        iszero(c) && continue
        moment = moment_value(table,w,result.solution.values)
        if isnothing(moment)
            required = cld(max(degree(w),degree(result.problem.objective),
                maximum((degree(c.polynomial) for c in result.problem.constraints);init=0)),2)
            throw(ArgumentError("Moment $(_word_string(p.registry,w)) is not in this program. Use Relaxation(order=$required) with complete rows, or custom rows whose adjoint products contain this word, and solve again."))
        end
        value += c*moment
    end
    return LinearAlgebra.ishermitian(p) ? real(value) : value
end
