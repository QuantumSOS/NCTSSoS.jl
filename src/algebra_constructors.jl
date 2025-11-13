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

"""
    fermionic_algebra(N::Int)

Create a fermionic algebra system for N fermionic modes.

Fermionic operators satisfy anti-commutation relations and nilpotency:
- Anti-commutation: {cᵢ, cⱼ} = 0, {cᵢ†, cⱼ†} = 0
- Canonical anti-commutation: {cᵢ, cⱼ†} = δᵢⱼ
- Nilpotency: cᵢ² = 0, (cᵢ†)² = 0

All properties are enforced through polynomial equality constraints.

# Arguments
- `N::Int`: Number of fermionic modes

# Returns
A `FermionicAlgebra` instance with fields:
- `N::Int`: Number of fermionic modes
- `variables::Tuple`: (c, c_dag) where c[i] and c_dag[i] are annihilation/creation operators
- `simplify_algo::SimplifyAlgorithm`: Simplification algorithm for commutation groups
- `equality_constraints::Vector`: All fermionic algebra constraints
- `inequality_constraints::Vector`: Empty (no inequality constraints)
- `comm_gps::Vector{Vector{Variable}}`: Commutation group structure

# Example
```julia
using NCTSSoS

# Create 2-mode fermionic system
sys = fermionic_algebra(2)
c, c_dag = sys.variables

# Number operator for mode 1
n1 = c_dag[1] * c[1]

# Create optimization problem
pop = cpolyopt(n1, sys)
```

# Mathematical Details
For N modes, the system includes:
- N(N+1)/2 constraints: {cᵢ, cⱼ} = 0
- N(N+1)/2 constraints: {cᵢ†, cⱼ†} = 0
- N² constraints: {cᵢ, cⱼ†} = δᵢⱼ
- 2N constraints: cᵢ² = 0, (cᵢ†)² = 0

Total: 2N² + 3N equality constraints
"""
function fermionic_algebra(N::Int)
    @assert N >= 1 "Number of modes N must be at least 1"

    # Declare variables for fermionic operators
    @ncpolyvar c[1:N] c_dag[1:N]

    # All operators in single group (no automatic commutation)
    all_ops = vcat(collect(c), collect(c_dag))
    comm_gps = [all_ops]

    # Build equality constraints
    equality_constraints = Polynomial{ComplexF64}[]

    # Anti-commutation: {cᵢ, cⱼ} = 0
    for i in 1:N, j in i:N
        push!(equality_constraints, ComplexF64(1.0) * (c[i] * c[j] + c[j] * c[i]))
    end

    # Anti-commutation: {cᵢ†, cⱼ†} = 0
    for i in 1:N, j in i:N
        push!(equality_constraints, ComplexF64(1.0) * (c_dag[i] * c_dag[j] + c_dag[j] * c_dag[i]))
    end

    # Canonical anti-commutation: {cᵢ, cⱼ†} = δᵢⱼ
    for i in 1:N, j in 1:N
        if i == j
            push!(equality_constraints, ComplexF64(1.0) * (c[i] * c_dag[i] + c_dag[i] * c[i]) - 1)
        else
            push!(equality_constraints, ComplexF64(1.0) * (c[i] * c_dag[j] + c_dag[j] * c[i]))
        end
    end

    # Nilpotent constraints: cᵢ² = 0, (cᵢ†)² = 0
    for i in 1:N
        push!(equality_constraints, ComplexF64(1.0) * c[i] * c[i])
        push!(equality_constraints, ComplexF64(1.0) * c_dag[i] * c_dag[i])
    end

    # Create SimplifyAlgorithm (no special flags needed)
    simplify_algo = SimplifyAlgorithm(
        comm_gps=comm_gps,
        is_unipotent=false,
        is_projective=false
    )

    # No inequality constraints
    inequality_constraints = empty(equality_constraints)

    # Return fermionic algebra system
    return FermionicAlgebra(
        N,
        (c, c_dag),
        simplify_algo,
        equality_constraints,
        inequality_constraints,
        comm_gps
    )
end
