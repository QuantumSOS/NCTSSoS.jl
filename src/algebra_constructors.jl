"""
    pauli_algebra(N::Int)

Create a Pauli algebra system for N spin-1/2 sites using FastPolynomials registry.

Pauli operators satisfy:
- Each Pauli matrix squares to identity
- Anticommutation relation: {sig_i, sig_j} = 2 delta_ij

These rules are automatically handled by FastPolynomials simplification.

# Arguments
- `N::Int`: Number of spin-1/2 sites (must be >= 1)

# Returns
A NamedTuple with fields:
- `registry::VariableRegistry{PauliAlgebra}`: Variable registry for the algebra
- `sx::Vector{Monomial}`: Pauli X operators for each site
- `sy::Vector{Monomial}`: Pauli Y operators for each site
- `sz::Vector{Monomial}`: Pauli Z operators for each site

# Examples
```julia
julia> sys = pauli_algebra(2);

julia> sys.sx[1] * sys.sy[1]  # Simplifies via Pauli rules
...

julia> ham = sys.sx[1]*sys.sx[2] + sys.sy[1]*sys.sy[2] + sys.sz[1]*sys.sz[2];

julia> pop = polyopt(ham, sys.registry)
...
```
"""
function pauli_algebra(N::Int)
    @assert N >= 1 "Number of sites N must be at least 1"
    registry, (sx, sy, sz) = create_pauli_variables(1:N)
    return (registry = registry, sx = sx, sy = sy, sz = sz)
end

"""
    fermionic_algebra(N::Int)

Create a fermionic algebra system for N modes using FastPolynomials registry.

Fermionic operators satisfy anticommutation relations:
- {a_i, a_j^dag} = delta_ij (creation-annihilation)
- {a_i, a_j} = 0 (annihilation-annihilation)
- {a_i^dag, a_j^dag} = 0 (creation-creation)

These rules are automatically handled by FastPolynomials simplification.

# Arguments
- `N::Int`: Number of fermionic modes (must be >= 1)

# Returns
A NamedTuple with fields:
- `registry::VariableRegistry{FermionicAlgebra}`: Variable registry for the algebra
- `a::Vector{Monomial}`: Annihilation operators for each mode
- `a_dag::Vector{Monomial}`: Creation operators for each mode

# Examples
```julia
julia> sys = fermionic_algebra(2);

julia> sys.a[1] * sys.a_dag[1]  # Simplifies via anticommutation
...
```
"""
function fermionic_algebra(N::Int)
    @assert N >= 1 "Number of modes N must be at least 1"
    registry, (a, a_dag) = create_fermionic_variables(1:N)
    return (registry = registry, a = a, a_dag = a_dag)
end

"""
    bosonic_algebra(N::Int)

Create a bosonic algebra system for N modes using FastPolynomials registry.

Bosonic operators satisfy commutation relations:
- [c_i, c_j^dag] = delta_ij (creation-annihilation)
- [c_i, c_j] = 0 (annihilation-annihilation)
- [c_i^dag, c_j^dag] = 0 (creation-creation)

These rules are automatically handled by FastPolynomials simplification.

# Arguments
- `N::Int`: Number of bosonic modes (must be >= 1)

# Returns
A NamedTuple with fields:
- `registry::VariableRegistry{BosonicAlgebra}`: Variable registry for the algebra
- `c::Vector{Monomial}`: Annihilation operators for each mode
- `c_dag::Vector{Monomial}`: Creation operators for each mode

# Examples
```julia
julia> sys = bosonic_algebra(2);

julia> sys.c[1] * sys.c_dag[1]  # Simplifies via commutation
...
```
"""
function bosonic_algebra(N::Int)
    @assert N >= 1 "Number of modes N must be at least 1"
    registry, (c, c_dag) = create_bosonic_variables(1:N)
    return (registry = registry, c = c, c_dag = c_dag)
end

"""
    projector_algebra(prefix::String, N::Int)
    projector_algebra(prefixes::Vector{String}, N::Int)

Create a projector algebra system using FastPolynomials registry.

Projector operators satisfy:
- P_i^2 = P_i (idempotency: projectors square to themselves)

This rule is automatically handled by FastPolynomials simplification.

# Arguments
- `prefix::String`: Symbol prefix for projector variables (e.g., "P")
- `prefixes::Vector{String}`: Multiple prefixes for grouped projectors
- `N::Int`: Number of projector indices per prefix (must be >= 1)

# Returns
A NamedTuple with fields:
- `registry::VariableRegistry{ProjectorAlgebra}`: Variable registry for the algebra
- `projectors::Tuple{Vector{Monomial},...}`: Tuple of projector vectors, one per prefix

# Examples
```julia
julia> sys = projector_algebra("P", 3);

julia> sys.projectors[1][1]  # P_1
...

julia> sys = projector_algebra(["P", "Q"], 2);

julia> length(sys.projectors)  # Two groups: P and Q
2
```
"""
function projector_algebra(prefix::String, N::Int)
    @assert N >= 1 "Number of indices N must be at least 1"
    registry, projectors = create_projector_variables([(prefix, 1:N)])
    return (registry = registry, projectors = projectors)
end

function projector_algebra(prefixes::Vector{String}, N::Int)
    @assert N >= 1 "Number of indices N must be at least 1"
    @assert !isempty(prefixes) "At least one prefix is required"
    prefix_subscripts = [(p, 1:N) for p in prefixes]
    registry, projectors = create_projector_variables(prefix_subscripts)
    return (registry = registry, projectors = projectors)
end

"""
    unipotent_algebra(prefix::String, N::Int)
    unipotent_algebra(prefixes::Vector{String}, N::Int)

Create a unipotent algebra system using FastPolynomials registry.

Unipotent operators satisfy:
- U^2 = I (squares to identity)

This is simpler than Pauli algebra - no cyclic products or cross-operator interactions.
This rule is automatically handled by FastPolynomials simplification.

# Arguments
- `prefix::String`: Symbol prefix for unipotent variables (e.g., "U")
- `prefixes::Vector{String}`: Multiple prefixes for grouped variables
- `N::Int`: Number of indices per prefix (must be >= 1)

# Returns
A NamedTuple with fields:
- `registry::VariableRegistry{UnipotentAlgebra}`: Variable registry for the algebra
- `variables::Tuple{Vector{Monomial},...}`: Tuple of variable vectors, one per prefix

# Examples
```julia
julia> sys = unipotent_algebra("U", 3);

julia> sys.variables[1][1] * sys.variables[1][1]  # U_1^2 = I
...

julia> sys = unipotent_algebra(["U", "V"], 2);

julia> length(sys.variables)  # Two groups: U and V
2
```
"""
function unipotent_algebra(prefix::String, N::Int)
    @assert N >= 1 "Number of indices N must be at least 1"
    registry, variables = create_unipotent_variables([(prefix, 1:N)])
    return (registry = registry, variables = variables)
end

function unipotent_algebra(prefixes::Vector{String}, N::Int)
    @assert N >= 1 "Number of indices N must be at least 1"
    @assert !isempty(prefixes) "At least one prefix is required"
    prefix_subscripts = [(p, 1:N) for p in prefixes]
    registry, variables = create_unipotent_variables(prefix_subscripts)
    return (registry = registry, variables = variables)
end

"""
    noncommutative_algebra(prefix::String, N::Int)
    noncommutative_algebra(prefixes::Vector{String}, N::Int)

Create a generic non-commutative algebra system using FastPolynomials registry.

Non-commutative variables have no simplification rules - word order is preserved
exactly as given.

# Arguments
- `prefix::String`: Symbol prefix for variables (e.g., "x")
- `prefixes::Vector{String}`: Multiple prefixes for grouped variables
- `N::Int`: Number of indices per prefix (must be >= 1)

# Returns
A NamedTuple with fields:
- `registry::VariableRegistry{NonCommutativeAlgebra}`: Variable registry for the algebra
- `variables::Tuple{Vector{Monomial},...}`: Tuple of variable vectors, one per prefix

# Examples
```julia
julia> sys = noncommutative_algebra("x", 3);

julia> sys.variables[1][1] * sys.variables[1][2]  # x_1 * x_2 (order preserved)
...

julia> sys = noncommutative_algebra(["x", "y"], 2);

julia> length(sys.variables)  # Two groups: x and y
2
```
"""
function noncommutative_algebra(prefix::String, N::Int)
    @assert N >= 1 "Number of indices N must be at least 1"
    registry, variables = create_noncommutative_variables([(prefix, 1:N)])
    return (registry = registry, variables = variables)
end

function noncommutative_algebra(prefixes::Vector{String}, N::Int)
    @assert N >= 1 "Number of indices N must be at least 1"
    @assert !isempty(prefixes) "At least one prefix is required"
    prefix_subscripts = [(p, 1:N) for p in prefixes]
    registry, variables = create_noncommutative_variables(prefix_subscripts)
    return (registry = registry, variables = variables)
end
