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
end

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

"One checked build and solve. The backend receives only standard data, not model or algebra objects."
function solve(problem::Problem,relaxation::AbstractRelaxation;solver=nothing)
    solver isa AbstractSolver || throw(ArgumentError("Missing or invalid solver description; load SCS with `using SCS` and pass solver=SCSBackend(), load Mosek with `using Mosek` and pass solver=MosekBackend(), or pass solver=JuMPBackend(optimizer_factory), with attributes on the factory if needed."))
    built = build(problem,relaxation)
    solution = solve_standard(moment_standard_form(built.program),solver)
    return SolveResult(solution.bound,solution.status,built.sizes,built.valid_for,problem,relaxation,solution,built.table,solver)
end
"Add an operator constraint, retaining symmetry and solver settings. Complete matrices use the backend's append path (Mosek rebuilds its task); selected blocks rebuild because closure can change their coordinates."
function solve(previous::SolveResult,constraint::OperatorConstraint)
    problem = add_constraint(previous.problem,constraint)
    built = build(problem,previous.relaxation)
    plain = previous.relaxation isa Relaxation && !_uses_blocks(problem,previous.relaxation)
    solution = solve_standard(moment_standard_form(built.program),plain ? previous.solution.handle : previous.solver)
    return SolveResult(solution.bound,solution.status,built.sizes,built.valid_for,problem,previous.relaxation,solution,built.table,previous.solver)
end

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
