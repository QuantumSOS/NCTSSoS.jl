# Temporary fix for the new interface
function cpolyopt(objective::P, algebra::NamedTuple; eq_constraints=Any[], ineq_constraints=Any[]) where {T,P<:AbstractPolynomial{T}}
    # Extract properties from algebra
    comm_gps = algebra.comm_gps
    is_unipotent = algebra.simplify_algo.is_unipotent
    is_projective = algebra.simplify_algo.is_projective
    
    # Merge algebra constraints with user-provided constraints
    # Convert user constraints to match the coefficient type of the objective
    user_eq = isempty(eq_constraints) ? typeof(algebra.equality_constraints)() : [T(1) * poly for poly in eq_constraints]
    merged_eq_constraints = vcat(algebra.equality_constraints, user_eq)
    
    # Call the original cpolyopt with merged constraints
    return cpolyopt(objective;
        eq_constraints=merged_eq_constraints,
        ineq_constraints=ineq_constraints,
        comm_gps=comm_gps,
        is_unipotent=is_unipotent,
        is_projective=is_projective)
end
