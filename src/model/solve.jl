struct SolveResult{P,R,S,B,V,Q,D}
    bound::Union{Nothing,Float64}
    status::S
    sizes::B
    valid_for::V
    problem::P
    relaxation::R
    solution::Q
    # Width is chosen at build time. Do not make solve's return type depend on it.
    table::Any
    solver::D
    # Appending constraints consumes a handle; evict it from its original build.
    build_cache::Union{Nothing,Dict{Any,Any}}
end
SolveResult(bound,status,sizes,valid_for,problem,relaxation,solution,table,solver) =
    SolveResult(bound,status,sizes,valid_for,problem,relaxation,solution,table,solver,nothing)

# Reuse the redundant tag's storage for the already-built moment table.
@inline function Base.getproperty(r::SolveResult,name::Symbol)
    name===:bound_kind && return isnothing(getfield(r,:bound)) ? :unavailable : :numerical_optimum_not_certificate
    return getfield(r,name)
end
Base.propertynames(::SolveResult,private::Bool=false) =
    (:bound,:bound_kind,:status,:sizes,:valid_for,:problem,:relaxation,:solution,:table,:solver)

function Base.show(io::IO,r::SolveResult)
    print(io,"Bound: ",r.bound," (",r.bound_kind,"); status: ",r.status,
        "\n",r.problem,"\n",r.relaxation,
        "\nSizes: ",r.sizes,"\nValid for: ",r.valid_for)
end
Base.show(io::IO,::MIME"text/plain",r::SolveResult) = show(io,r)

function _check_solver(solver)
    solver isa AbstractSolver || throw(ArgumentError("Missing or invalid solver description; load SCS with `using SCS` and pass solver=SCSBackend(), load Mosek with `using Mosek` and pass solver=MosekBackend(), or pass solver=JuMPBackend(optimizer_factory), with attributes on the factory if needed."))
    return solver
end
function _standard_form!(built::BuiltRelaxation)
    get!(built.cache,:moment) do
        moment_standard_form(built.program)
    end::MomentStandardForm
end
_prepare(built::BuiltRelaxation,solver) = prepare_standard(_standard_form!(built),solver)

"""
    NCTSSoS.prepare(built::NCTSSoS.BuiltRelaxation; solver)

Lower solver data and prepare its task without optimizing; cache and return the
prepared handle. Repeating this call with the same solver description reuses it.
Forms are created only when requested, and backend options have separate tasks.
Repeated solve(built; solver) calls through native SCSBackend and MosekBackend
start cold; only solve(result, constraint) may warm-start an appended constraint.
JuMPBackend reuses the optimizer's workspace, so repeated results depend on the
optimizer (COSMO starts from retained state); use a new NCTSSoS.build for a cold solve.
For SCS, this caches cone ordering and sparse arrays, but solve still includes
the vendor's one-shot setup, reported separately as solution.diagnostics.setup_seconds.
Use NCTSSoS.build for formulation, this call for canonicalization/preparation, and
solve(built; solver) for optimization. Do not mutate built.program or its inputs;
call build again after an edit. A cached native task is not safe for concurrent solves.
"""
function prepare(built::BuiltRelaxation;solver=nothing)
    _check_solver(solver)
    get!(built.cache,solver) do
        _prepare(built,solver)
    end
end

"""
    result = solve(built; solver)
    result = solve(problem, relaxation; solver)
    result = solve(previous, constraint)

Solve a relaxation and return its numerical result.
The first form prepares solver data when needed.
The second form builds the relaxation first.
The third form adds a constraint and keeps the previous solver settings.

# Result fields

- `bound`: the numerical relaxation optimum, or `nothing` when no solution is accepted.
  This value is not a rigorous certificate.
- `bound_kind`: `:numerical_optimum_not_certificate` when `bound` is available;
  `:unavailable` otherwise.
- `valid_for`: a printable description of the states covered by the relaxation.
  It includes the input constraints and any conditions from strengthening.
  For example, `Stationarity()` restricts validity to eigenstates of the objective.
  `StateOptimality()` restricts validity to its ground space.
- `status`: the solver termination status, such as `OPTIMAL`, `TIME_LIMIT` or `INFEASIBLE`.
  Check `bound !== nothing` before reading expectations.
- `sizes`: the program counts listed below.
- `solution.diagnostics`: the backend reports listed below.

# Program sizes

Read these fields from `result.sizes` or `inspect(problem, relaxation).sizes`.
They describe the relaxation before the solver changes its representation or runs presolve.

- `rows`: the number of moment-matrix rows on a plain relaxation.
  With a block family, it is the sum of row counts in the closed family, before symmetry reduction.
- `moments`: the number of retained moment words or symmetry representatives, including the constant.
- `real_variables`: the number of real optimization variables, excluding the fixed constant.
- `blocks`: the total number of cone blocks, including equality blocks.
- `largest_block`: the largest PSD or Hermitian PSD matrix side.
  A Hermitian side is its complex dimension, before conversion to a real matrix.
- `localizing_rows`: a vector of row counts for the operator-constraint localizing matrices.
- `closure`: `nothing` for a plain relaxation.
  For a block relaxation, it reports family closure, symmetry reduction and strengthening.
  Print `result.sizes` for a summary of these changes.

# Backend diagnostics

Read `d = result.solution.diagnostics`.
These are solver reports, not certificates.

- `JuMPBackend`: `d` is `nothing`.
- `SCSBackend`: `d` has `primal_residual`, `dual_residual`, `gap`, `iterations`,
  `setup_seconds`, `solve_seconds` and `warm_started`.
  The times are SCS's reported setup and solve times in seconds.
  They exclude package-side building and preparation.
- `MosekBackend`: `d` has `representation`, `termination`, `warm_started`,
  `solve_seconds`, `iterations`, `primal_residual` and `dual_residual`.
  `representation` is `:moment` or `:sos`.
  `termination` is the raw optimizer termination code as text.
  `warm_started` is false.
  `solve_seconds` is `MSK_DINF_OPTIMIZER_TIME`, including presolve for this optimizer call.
  It excludes package-side building, task preparation and solution extraction.
  Residuals describe the native formulation, not the original problem.
  A residual of `-1` means that Mosek did not report that measure.
  See [`MosekBackend`](@ref) for the information items and interpretation.

Example: `using NCTSSoS, SCS; registry, (x, y, z) = pauli_variables(1); result = solve(Problem(z[1]), Relaxation(order=1); solver=SCSBackend())`.
"""
function solve(built::BuiltRelaxation;solver=nothing)
    solution = solve_standard(prepare(built;solver))
    return SolveResult(solution.bound,solution.status,built.sizes,built.valid_for,built.problem,
        built.relaxation,solution,built.table,solver,built.cache)
end

"One checked build and solve, equivalent to solve(NCTSSoS.build(problem, relaxation); solver)."
function solve(problem::Problem,relaxation::AbstractRelaxation;solver=nothing)
    _check_solver(solver)
    return solve(build(problem,relaxation);solver)
end
"Add an operator or scalar expectation constraint, retaining symmetry and solver settings. Complete matrices use the backend's append path (Mosek rebuilds its task); selected blocks rebuild because closure can change their coordinates."
function solve(previous::SolveResult,constraint::Union{OperatorConstraint,ExpectationConstraint})
    problem = add_constraint(previous.problem,constraint)
    built = build(problem,previous.relaxation)
    plain = previous.relaxation isa Relaxation && !_uses_blocks(problem,previous.relaxation)
    # A complex Hermitian constraint can invalidate averaging over real moments.
    (!plain || previous.solver isa MosekBackend || _real_form(problem)!=_real_form(previous.problem)) &&
        return solve(built;solver=previous.solver)
    isnothing(previous.build_cache) || delete!(previous.build_cache,previous.solver)
    solution = solve_standard(_standard_form!(built),previous.solution.handle)
    built.cache[previous.solver]=solution.handle
    return SolveResult(solution.bound,solution.status,built.sizes,built.valid_for,problem,
        previous.relaxation,solution,built.table,previous.solver,built.cache)
end

"""
    expectation(result::SolveResult, p::Polynomial)

Evaluate ⟨p⟩ in the numerical relaxation solution, not a reconstructed state or
certificate. `p` must use the problem's registry and only represented moments.
Return a real number for Hermitian `p`, and a complex number otherwise.
For state/trace polynomials this reads the outer moment of the whole expression.
In particular, `expectation(result, state_symbol(A)*state_symbol(B))` need not
equal `expectation(result,state_symbol(A))*expectation(result,state_symbol(B))`.
"""
function expectation(result::SolveResult,p::Polynomial)
    isnothing(result.bound) && throw(ArgumentError("Cannot read expectations: solve status $(result.status) has no accepted solution."))
    p.registry === result.problem.objective.registry ||
        throw(ArgumentError("Expectation polynomial uses a different registry; construct it from the problem's variables."))
    return _expectation(result.table,result,p)
end

function _missing_moment(result::SolveResult,p::Polynomial,w)
    if result.relaxation isa Relaxation && !_uses_blocks(result.problem,result.relaxation)
        required = cld(max(degree(w),degree(result.problem.objective),
            maximum((degree(c.polynomial) for c in result.problem.constraints);init=0)),2)
        throw(ArgumentError("Moment $(_word_string(p.registry,w)) is not in this program. Use Relaxation(order=$required) with complete rows, or custom rows whose adjoint products contain this word, and solve again."))
    end
    throw(ArgumentError("Moment $(_word_string(p.registry,w)) is not among the moments of the selected blocks. Try a higher order, no split (split=nothing), or a site group containing the word's sites ($(join(sites(w),", "))), and solve again."))
end

# Specialize the evaluation loop on the table's packed word type only on demand.
function _expectation(table,result::SolveResult,p::Polynomial{W}) where {W}
    value = 0.0+0.0im
    if applicable(moment_terms,table,one(W))
        # Keep coefficient precision; narrowing a tiny coefficient could falsely
        # remove its dependence on an unconstrained coordinate.
        coordinates=Dict{Int,Tuple{Number,W}}()
        for (c,w) in terms(p)
            iszero(c) && continue
            pairs=moment_terms(table,w)
            isnothing(pairs) && _missing_moment(result,p,w)
            for (id,factor) in pairs
                coefficient=c*factor
                iszero(coefficient) && continue
                previous=get(coordinates,id,(0,w))[1]
                coordinates[id]=(previous+coefficient,w)
            end
        end
        for (id,(c,w)) in coordinates
            # Zero tolerance, including general floating coefficients: any
            # nonzero dependence on an unconstrained coordinate is undetermined.
            # Opposite contributions cancel before this coverage check.
            iszero(c) && continue
            moment_covered(table,id) || _missing_moment(result,p,w)
            value+=c*(id==1 ? 1.0 : result.solution.values[id-1])
        end
    else
        # Compatibility until a table provider implements both coordinate hooks.
        for (c,w) in terms(p)
            iszero(c) && continue
            moment=moment_value(table,w,result.solution.values)
            isnothing(moment) && _missing_moment(result,p,w)
            value+=c*moment
        end
    end
    return LinearAlgebra.ishermitian(p) ? Float64(real(value)) : ComplexF64(value)
end
