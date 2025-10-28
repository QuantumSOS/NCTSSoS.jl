"""
    pauli_algebra(N::Int)

Create a Pauli algebra system for N spin-1/2 sites.
"""
function pauli_algebra(N::Int)
    @assert N >= 1 "Number of sites N must be at least 1"
    
    # Declare non-commuting variables for Pauli operators
    @ncpolyvar x[1:N] y[1:N] z[1:N]
    
    # Create commutation groups (operators at same site commute with each other)
    comm_gps = [[x[i], y[i], z[i]] for i in 1:N]
    
    # Build equality constraints encoding Pauli commutation relations
    equality_constraints = reduce(vcat, [
        [x[i]*y[i] - im*z[i], 
         y[i]*x[i] + im*z[i],
         y[i]*z[i] - im*x[i], 
         z[i]*y[i] + im*x[i],
         z[i]*x[i] - im*y[i], 
         x[i]*z[i] + im*y[i]]
        for i in 1:N
    ])
    
    # Create SimplifyAlgorithm with unipotent property
    simplify_algo = SimplifyAlgorithm(comm_gps=comm_gps, is_unipotent=true, is_projective=false)
    
    # No inequality constraints for Pauli algebra
    inequality_constraints = empty(equality_constraints)
    
    return (
        variables = (x, y, z),
        simplify_algo = simplify_algo,
        equality_constraints = equality_constraints,
        inequality_constraints = inequality_constraints,
        comm_gps = comm_gps
    )
end


"""
    bosonic_algebra(N::Int)

Create a bosonic algebra system for N harmonic oscillator modes.
"""
function bosonic_algebra(N::Int)
    @assert N >= 1 "Number of modes N must be at least 1"
    
    # Declare non-commuting variables for position and momentum operators
    @ncpolyvar q[1:N] p[1:N]
    
    # Create commutation groups (each mode's q and p form a group)
    comm_gps = [[q[i], p[i]] for i in 1:N]
    
    # Build equality constraints encoding canonical commutation relations
    equality_constraints = [q[i]*p[i] - p[i]*q[i] - 1.0im for i in 1:N]
    
    # Create SimplifyAlgorithm with no special properties
    simplify_algo = SimplifyAlgorithm(comm_gps=comm_gps, is_unipotent=false, is_projective=false)
    
    # No inequality constraints for bosonic algebra
    inequality_constraints = empty(equality_constraints)
    
    return (
        variables = (q, p),
        simplify_algo = simplify_algo,
        equality_constraints = equality_constraints,
        inequality_constraints = inequality_constraints,
        comm_gps = comm_gps
    )
end
