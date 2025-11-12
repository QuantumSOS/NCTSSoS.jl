"""
    AbstractAlgebra

Abstract base type for all algebra systems in NCTSSoS.
Concrete subtypes include `PauliAlgebra` and `FermionicAlgebra`.
"""
abstract type AbstractAlgebra end

"""
    PauliAlgebra <: AbstractAlgebra

Pauli algebra system for spin-1/2 sites.

# Fields
- `N::Int`: Number of spin-1/2 sites
- `variables::Tuple`: Tuple of (x, y, z) variable arrays
- `simplify_algo::SimplifyAlgorithm`: Simplification algorithm configuration
- `equality_constraints::Vector`: Polynomial equality constraints
- `inequality_constraints::Vector`: Polynomial inequality constraints
- `comm_gps::Vector`: Commutation groups
"""
struct PauliAlgebra <: AbstractAlgebra
    N::Int
    variables::Tuple{Vector{Variable},Vector{Variable},Vector{Variable}}
    simplify_algo::SimplifyAlgorithm
    equality_constraints::Vector{Polynomial{ComplexF64}}
    inequality_constraints::Vector{Polynomial{ComplexF64}}
    comm_gps::Vector{Vector{Variable}}
end

"""
    FermionicAlgebra <: AbstractAlgebra

Fermionic algebra system for fermionic modes.

# Fields
- `N::Int`: Number of fermionic modes
- `variables::Tuple`: Tuple of (c, c_dag) variable arrays
- `simplify_algo::SimplifyAlgorithm`: Simplification algorithm configuration
- `equality_constraints::Vector`: Polynomial equality constraints
- `inequality_constraints::Vector`: Polynomial inequality constraints
- `comm_gps::Vector`: Commutation groups
"""
struct FermionicAlgebra <: AbstractAlgebra
    N::Int
    variables::Tuple{Vector{Variable},Vector{Variable}}
    simplify_algo::SimplifyAlgorithm
    equality_constraints::Vector{Polynomial{ComplexF64}}
    inequality_constraints::Vector{Polynomial{ComplexF64}}
    comm_gps::Vector{Vector{Variable}}
end

"""
    pauli_algebra(N::Int)

Create a Pauli algebra system for N spin-1/2 sites.
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

    # Create SimplifyAlgorithm with unipotent property
    simplify_algo = SimplifyAlgorithm(comm_gps=comm_gps, is_unipotent=true, is_projective=false)

    # No inequality constraints for Pauli algebra
    inequality_constraints = empty(equality_constraints)

    return PauliAlgebra(
        N,
        (x, y, z),
        simplify_algo,
        equality_constraints,
        inequality_constraints,
        comm_gps
    )
end
