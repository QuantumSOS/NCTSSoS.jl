struct SolveResult{P,R,S,B,V,Q}
    bound::Union{Nothing,Float64}
    bound_kind::Symbol
    status::S
    sizes::B
    valid_for::V
    problem::P
    relaxation::R
    solution::Q
end

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
    kind = isnothing(solution.bound) ? :unavailable : :numerical_optimum_not_certificate
    return SolveResult(solution.bound,kind,solution.status,built.sizes,built.valid_for,problem,relaxation,solution)
end
"Add an operator constraint and solve again on the same open moment-form handle."
function solve(previous::SolveResult,constraint::OperatorConstraint)
    problem = add_constraint(previous.problem,constraint)
    built = build(problem,previous.relaxation)
    solution = solve_standard(moment_standard_form(built.program),previous.solution.handle)
    kind = isnothing(solution.bound) ? :unavailable : :numerical_optimum_not_certificate
    return SolveResult(solution.bound,kind,solution.status,built.sizes,built.valid_for,problem,previous.relaxation,solution)
end
