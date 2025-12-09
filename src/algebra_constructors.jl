"""
    pauli_algebra(N::Int)

Create a Pauli algebra system for N spin-1/2 sites.

# Returns
A NamedTuple with fields:
- `variables`: Tuple of (x, y, z) Pauli operator arrays
- `is_unipotent`: true (Pauli operators square to identity)
- `is_projective`: false
- `equality_constraints`: Vector of Pauli commutation relations
- `inequality_constraints`: Empty vector
- `comm_gps`: Commutation groups (operators at same site)
"""
function pauli_algebra(N::Int)
    @assert N >= 1 "Number of sites N must be at least 1"

    # Declare non-commuting variables for Pauli operators
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    # Create commutation groups: operators at different sites commute, while operators at the same site form a commutation group for tracking purposes (Pauli operators at the same site anti-commute).
    comm_gps = [[x[i], y[i], z[i]] for i in 1:N]

    # Build equality constraints encoding Pauli commutation relations
    equality_constraints = reduce(vcat, [
        [x[i] * y[i] - im * z[i],
            y[i] * x[i] + im * z[i],
            y[i] * z[i] - im * x[i],
            z[i] * y[i] + im * x[i],
            z[i] * x[i] - im * y[i],
            x[i] * z[i] + im * y[i]]
        for i in 1:N
    ])

    # No inequality constraints for Pauli algebra
    inequality_constraints = empty(equality_constraints)

    return (
        variables=(x, y, z),
        is_unipotent=true,
        is_projective=false,
        equality_constraints=equality_constraints,
        inequality_constraints=inequality_constraints,
        comm_gps=comm_gps
    )
end
