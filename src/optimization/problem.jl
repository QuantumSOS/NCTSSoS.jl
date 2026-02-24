"""
    OptimizationProblem{A<:AlgebraType, P}

Abstract type for polynomial optimization problems with algebra type tracking.

# Type Parameters
- `A<:AlgebraType`: The algebra type (PauliAlgebra, FermionicAlgebra, etc.)
- `P`: Type of polynomial (Polynomial{A,T,C})

The algebra type parameter enables dispatch on algebraic structure:
- `PauliAlgebra`: Pauli spin matrices with σ² = I and cyclic products
- `FermionicAlgebra`: Fermionic operators with anticommutation
- `BosonicAlgebra`: Bosonic operators with commutation
- `ProjectorAlgebra`: Projectors with P² = P (idempotent)
- `UnipotentAlgebra`: Operators with U² = I (involutory)
- `NonCommutativeAlgebra`: Generic non-commutative (no simplification)
"""
abstract type OptimizationProblem{A<:AlgebraType, P} end

"""
    PolyOpt{A<:AlgebraType, T<:Integer, P<:AbstractPolynomial} <: OptimizationProblem{A, P}

A polynomial optimization problem structure with algebra type tracking.
Handles both standard polynomials and state polynomials via the P parameter.

# Type Parameters
- `A<:AlgebraType`: The algebra type governing simplification rules
- `T<:Integer`: The integer type for monomial word indices
- `P<:AbstractPolynomial`: Type of polynomial (Polynomial{A,T,C} or NCStatePolynomial{C,ST,A,T})

# Fields
- `objective::P`: The polynomial objective function to be optimized
- `eq_constraints::Vector{P}`: Vector of equality constraints (p = 0)
- `ineq_constraints::Vector{P}`: Vector of inequality constraints (p >= 0)
- `registry::VariableRegistry{A,T}`: Variable registry mapping symbols to indices

# Notes
- Algebra type determines simplification rules (no manual comm_gps, is_unipotent, is_projective)
- Registry provides bidirectional symbol <-> index mapping for variable access
- For NCStatePolynomial objectives, the state type (Arbitrary, MaxEntangled) is embedded in P

# Examples
```julia
# Create Pauli algebra
reg, (σx, σy, σz) = create_pauli_variables(1:2)

# Build Hamiltonian
ham = 0.5 * (σx[1]*σx[2] + σy[1]*σy[2] + σz[1]*σz[2])

# Create optimization problem
pop = polyopt(ham, reg)

# State polynomial optimization (Bell inequalities)
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
sp = -1.0 * ς(x[1]*y[1]) - 1.0 * ς(x[1]*y[2]) - 1.0 * ς(x[2]*y[1]) + 1.0 * ς(x[2]*y[2])
pop = polyopt(sp * one(NormalMonomial), reg)
```

See also: [`polyopt`](@ref), [`VariableRegistry`](@ref), [`AlgebraType`](@ref)
"""
struct PolyOpt{A<:AlgebraType,T<:Integer,P<:AbstractPolynomial} <: OptimizationProblem{A,P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    registry::VariableRegistry{A,T}
end


"""
    polyopt(objective::P, registry::VariableRegistry{A,T};
            eq_constraints=P[], ineq_constraints=P[]) where {P<:AbstractPolynomial, A, T}

Create a polynomial optimization problem from objective, registry, and optional constraints.

Works with any `AbstractPolynomial` subtype including `Polynomial{A,T,C}` and
`NCStatePolynomial{C,ST,A,T}`.

# Arguments
- `objective::P`: The polynomial objective function to optimize
- `registry::VariableRegistry{A,T}`: Variable registry for the algebra

# Keyword Arguments
- `eq_constraints`: Equality constraints as polynomials (p = 0). Default: empty
- `ineq_constraints`: Inequality constraints as polynomials (p >= 0). Default: empty

# Returns
A `PolyOpt{A,T,P}` structure representing the optimization problem.

# Notes
- Algebra type `A` is inferred from the registry
- Coefficient type cannot be an integer subtype (JuMP solver requirement)
- Simplification rules are determined by the algebra type, not manual flags
- For `FermionicAlgebra`: objectives should have even parity (parity superselection rule).
  Odd-parity operators have zero expectation value. Validation is done during moment
  relaxation via `_add_parity_constraints!`.

# Examples
```julia
# Pauli algebra optimization
reg, (σx, σy, σz) = create_pauli_variables(1:3)
ham = 0.5 * (σx[1]*σx[2] + σy[1]*σy[2])
pop = polyopt(ham, reg)

# With equality constraints
constraint = σx[1]*σx[1] - one(typeof(ham))  # σx² = I (auto-simplified anyway)
pop = polyopt(ham, reg; eq_constraints=[constraint])

# State polynomial optimization (Bell inequalities)
reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
sp = -1.0 * ς(x[1]*y[1]) - 1.0 * ς(x[1]*y[2]) - 1.0 * ς(x[2]*y[1]) + 1.0 * ς(x[2]*y[2])
pop = polyopt(sp * one(NormalMonomial), reg)
```

See also: [`PolyOpt`](@ref), [`VariableRegistry`](@ref), [`NCStatePolynomial`](@ref)
"""
function polyopt(
    objective::P,
    registry::VariableRegistry{A,T};
    eq_constraints::Vector{P}=P[],
    ineq_constraints::Vector{P}=P[]
) where {P<:AbstractPolynomial, A<:AlgebraType, T<:Integer}
    C = coeff_type(P)
    if C <: Integer
        throw(ArgumentError("Polynomial coefficients cannot be integers (not supported by JuMP solvers). Use Float64 or other floating-point types."))
    end

    # Deduplicate constraints
    eq_cons = unique!(copy(eq_constraints))
    ineq_cons = unique!(copy(ineq_constraints))

    return PolyOpt{A,T,P}(objective, eq_cons, ineq_cons, registry)
end


"""
    Base.show(io::IO, pop::OptimizationProblem{A,P})

Display an optimization problem showing objective, constraints, and algebra type.
Uses registry-aware display for polynomials (shows symbolic variable names).
"""
function Base.show(io::IO, pop::OptimizationProblem{A,P}) where {A,P}
    # Helper to format a polynomial as string with registry for symbolic display
    function poly_str(p)
        buf = IOBuffer()
        show(IOContext(buf, :registry => pop.registry), p)
        String(take!(buf))
    end

    function cons_str(cons::Vector, iseq::Bool)
        join([poly_str(c) * (iseq ? " = 0" : " >= 0") for c in cons], "\n            ")
    end

    nvars = length(pop.registry)
    var_syms = symbols(pop.registry)
    var_str = nvars <= 10 ? join(string.(var_syms), ", ") :
              join(string.(var_syms[1:5]), ", ") * ", ..., " * join(string.(var_syms[end-2:end]), ", ")

    res_str = """
        Optimization Problem ($(nameof(A)))
        ────────────────────────────────────
        Objective:
            $(poly_str(pop.objective))

        Equality constraints ($(length(pop.eq_constraints))):
            $(isempty(pop.eq_constraints) ? "(none)" : cons_str(pop.eq_constraints, true))

        Inequality constraints ($(length(pop.ineq_constraints))):
            $(isempty(pop.ineq_constraints) ? "(none)" : cons_str(pop.ineq_constraints, false))

        Variables ($nvars):
            $var_str
        """
    print(io, res_str)
end


# =============================================================================
# Algebra Trait Functions
# =============================================================================

"""
    _is_complex_problem(::Type{A}) where {A<:AlgebraType} -> Bool

Determine if an algebra type requires complex moment relaxation (Hermitian PSD).

Returns true for algebras that produce complex phases during simplification:
- `PauliAlgebra`: Pauli products produce i phases (sigma_x * sigma_y = i*sigma_z)
- `FermionicAlgebra`: Anticommutation can produce phases
- `BosonicAlgebra`: Normal ordering can produce complex terms

Returns false for "real" algebras:
- `NonCommutativeAlgebra`: No simplification rules, no complex phases
- `ProjectorAlgebra`: P^2 = P, no phases
- `UnipotentAlgebra`: U^2 = I, no phases

This trait is used by interface.jl to dispatch between real and Hermitian
constraint handling in moment_relax.
"""
# @noinline: one-liner trait methods are inlined away by the compiler, making them
# invisible to Julia's code-coverage instrumentation. The runtime cost is negligible
# because these only gate SDP-solver calls.
@noinline function _is_complex_problem(::Type{PauliAlgebra})
    return true
end

@noinline function _is_complex_problem(::Type{FermionicAlgebra})
    return true
end

@noinline function _is_complex_problem(::Type{BosonicAlgebra})
    return true
end

@noinline function _is_complex_problem(::Type{NonCommutativeAlgebra})
    return false
end

@noinline function _is_complex_problem(::Type{ProjectorAlgebra})
    return false
end

@noinline function _is_complex_problem(::Type{UnipotentAlgebra})
    return false
end

