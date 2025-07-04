abstract type OptimizationProblem end

"""
    PolyOpt{P} <: OptimizationProblem

A polynomial optimization problem structure.

# Type Parameters
- `P`: Type of polynomial, either `Polynomial{T}` or `NCStatePolynomial{T}`

# Fields
- `objective::P`: The polynomial objective function to be optimized
- `eq_constraints::Vector{P}`: Vector of equality constraints (assumed to equal 0)
- `ineq_constraints::Vector{P}`: Vector of inequality constraints (assumed to be >= 0)
- `variables::Vector{Variable}`: All variables appearing in the problem
- `comm_gps::Vector{Vector{Variable}}`: Commutative groups - vectors of variables that commute with variables not in the same group
- `is_unipotent::Bool`: Whether variables square to 1 (e.g., Pauli operators, SWAP operators)
- `is_projective::Bool`: Whether variables are projective (X² = X)

# Notes
- All constraints are assumed to be simplified using `comm_gp`, `is_unipotent`, and `is_projective`
- The problem cannot be both unipotent and projective simultaneously
- Commutative groups must be disjoint sets
"""
struct PolyOpt{P} <: OptimizationProblem
    objective::P
    eq_constraints::Vector{P} # NOTE: assuming constraints are all simplified using comm_gp, is_unipotent, and is_projective
    ineq_constraints::Vector{P}
    variables::Vector{Variable}
    comm_gps::Vector{Vector{Variable}} # Vectors of Set of variables that commutes with variables not in the same set
    is_unipotent::Bool # square to 1. Examples: Pauli Operators, SWAP Operators
    is_projective::Bool # X^2 = X. Is projective.
end

"""
    polyopt(objective::P; eq_constraints=Any[], ineq_constraints=Any[], comm_gps=Vector{Variable}[], is_unipotent::Bool=false, is_projective::Bool=false) where {T,P<:AbstractPolynomial{T}}

Create a polynomial optimization problem.

# Arguments
- `objective::P`: The polynomial objective function to optimize.
- `eq_constraints=Any[]`: Equality constraints as polynomials (p = 0).
- `ineq_constraints=Any[]`: Inequality constraints as polynomials (p ≥ 0).
- `comm_gps=Vector{Variable}[]`: Groups of variables that commute. If empty, all variables are assumed to commute.
- `is_unipotent::Bool=false`: Flag indicating if the problem is unipotent.
- `is_projective::Bool=false`: Flag indicating if the problem is projective.

# Returns
A `PolyOpt{P}` structure representing the polynomial optimization problem.

# Notes
- The polynomial coefficients cannot be integers as they are not supported by JuMP solvers.
- Commutative groups must be disjoint, and all commutative variables must be a subset of all variables.
- A problem cannot be both unipotent and projective simultaneously.
"""
function polyopt(objective::P; eq_constraints=Any[], ineq_constraints=Any[], comm_gps=Vector{Variable}[], is_unipotent::Bool=false, is_projective::Bool=false) where {T,P<:AbstractPolynomial{T}}
    @assert !(T <: Integer) "The polynomial coefficients can not be integers (not supported by JuMP solvers)."
    eq_cons = unique!(collect(P, eq_constraints))
    ineq_cons = unique!(collect(P, ineq_constraints))
    vars = sorted_union(variables(objective), variables.(eq_cons)..., variables.(ineq_cons)...)
    if !isempty(comm_gps)
        @assert all([isempty(intersect(gp_a, gp_b)) for gp_a in comm_gps, gp_b in comm_gps if gp_a != gp_b]) "The commutative groups must be disjoint."
        @assert issubset(union(comm_gps...), vars) "The commutative variables must be a subset of the variables."
    else
        push!(comm_gps, vars)
    end
    @assert !(is_unipotent && is_projective) "The problem cannot be both unipotent and projective."
    return PolyOpt{P}(objective, eq_cons, ineq_cons, vars, comm_gps, is_unipotent, is_projective)
end

function Base.show(io::IO, pop::PolyOpt)
    cons_str(cons::Vector{P}, iseq::Bool) where {P} =
        join(["$(string(c)) " * (iseq ? "= 0" : ">= 0 \n") for c in cons], " \t")
    res_str = """
        obj: \n
            $(string(pop.objective)) \n
        constraints: \n
            $(cons_str(pop.eq_constraints,true))
            $(cons_str(pop.ineq_constraints,false))
        variables:
            $(join(string.(pop.variables)," ")) \n
        is_unipotent:
            $(pop.is_unipotent) \n
        is_projective:
            $(pop.is_projective) \n
    """
    print(io, res_str)
end
