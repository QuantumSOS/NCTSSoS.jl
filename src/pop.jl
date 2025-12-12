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
    PolyOpt{A<:AlgebraType, P<:Polynomial{A}} <: OptimizationProblem{A, P}

A polynomial optimization problem structure with algebra type tracking.

# Type Parameters
- `A<:AlgebraType`: The algebra type governing simplification rules
- `P<:Polynomial{A}`: Type of polynomial matching the algebra

# Fields
- `objective::P`: The polynomial objective function to be optimized
- `eq_constraints::Vector{P}`: Vector of equality constraints (p = 0)
- `ineq_constraints::Vector{P}`: Vector of inequality constraints (p >= 0)
- `registry::VariableRegistry{A}`: Variable registry mapping symbols to indices

# Notes
- Algebra type determines simplification rules (no manual comm_gps, is_unipotent, is_projective)
- Registry provides bidirectional symbol <-> index mapping for variable access

# Examples
```julia
# Create Pauli algebra
reg, (σx, σy, σz) = create_pauli_variables(1:2)

# Build Hamiltonian
ham = 0.5 * (σx[1]*σx[2] + σy[1]*σy[2] + σz[1]*σz[2])

# Create optimization problem
pop = polyopt(ham, reg)
```

See also: [`polyopt`](@ref), [`VariableRegistry`](@ref), [`AlgebraType`](@ref)
"""
struct PolyOpt{A<:AlgebraType, P<:Polynomial{A}} <: OptimizationProblem{A, P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    registry::VariableRegistry{A}
end


"""
    polyopt(objective::Polynomial{A,T,C}, registry::VariableRegistry{A,T};
            eq_constraints=Polynomial{A,T,C}[],
            ineq_constraints=Polynomial{A,T,C}[]) where {A<:AlgebraType, T<:Integer, C<:Number}

Create a polynomial optimization problem from objective, registry, and optional constraints.

# Arguments
- `objective::Polynomial{A,T,C}`: The polynomial objective function to optimize
- `registry::VariableRegistry{A,T}`: Variable registry for the algebra

# Keyword Arguments
- `eq_constraints`: Equality constraints as polynomials (p = 0). Default: empty
- `ineq_constraints`: Inequality constraints as polynomials (p >= 0). Default: empty

# Returns
A `PolyOpt{A, Polynomial{A,T,C}}` structure representing the optimization problem.

# Notes
- Algebra type `A` is inferred from the polynomial and registry (must match)
- Coefficient type `C` cannot be an integer subtype (JuMP solver requirement)
- Simplification rules are determined by the algebra type, not manual flags

# Examples
```julia
# Pauli algebra optimization
reg, (σx, σy, σz) = create_pauli_variables(1:3)
ham = 0.5 * (σx[1]*σx[2] + σy[1]*σy[2])
pop = polyopt(ham, reg)

# With equality constraints
constraint = σx[1]*σx[1] - one(typeof(ham))  # σx² = I (auto-simplified anyway)
pop = polyopt(ham, reg; eq_constraints=[constraint])
```

See also: [`PolyOpt`](@ref), [`VariableRegistry`](@ref)
"""
function polyopt(
    objective::Polynomial{A,T,C},
    registry::VariableRegistry{A,T};
    eq_constraints::Vector{Polynomial{A,T,C}}=Polynomial{A,T,C}[],
    ineq_constraints::Vector{Polynomial{A,T,C}}=Polynomial{A,T,C}[]
) where {A<:AlgebraType, T<:Integer, C<:Number}
    @assert !(C <: Integer) "The polynomial coefficients cannot be integers (not supported by JuMP solvers)."

    # Deduplicate constraints
    eq_cons = unique!(copy(eq_constraints))
    ineq_cons = unique!(copy(ineq_constraints))

    P = Polynomial{A,T,C}
    return PolyOpt{A, P}(objective, eq_cons, ineq_cons, registry)
end


"""
    Base.show(io::IO, pop::OptimizationProblem{A,P})

Display an optimization problem showing objective, constraints, and algebra type.
"""
function Base.show(io::IO, pop::OptimizationProblem{A,P}) where {A,P}
    function cons_str(cons::Vector, iseq::Bool)
        join(["$(string(c)) " * (iseq ? "= 0" : ">= 0") for c in cons], "\n            ")
    end

    nvars = length(pop.registry)
    var_syms = symbols(pop.registry)
    var_str = nvars <= 10 ? join(string.(var_syms), ", ") :
              join(string.(var_syms[1:5]), ", ") * ", ..., " * join(string.(var_syms[end-2:end]), ", ")

    res_str = """
        Optimization Problem ($(nameof(A)))
        ────────────────────────────────────
        Objective:
            $(string(pop.objective))

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

This trait is used by interface.jl to dispatch between moment_relax (real)
and complex_moment_relax (Hermitian).
"""
_is_complex_problem(::Type{PauliAlgebra}) = true
_is_complex_problem(::Type{FermionicAlgebra}) = true
_is_complex_problem(::Type{BosonicAlgebra}) = true
_is_complex_problem(::Type{NonCommutativeAlgebra}) = false
_is_complex_problem(::Type{ProjectorAlgebra}) = false
_is_complex_problem(::Type{UnipotentAlgebra}) = false
_is_complex_problem(::Type{A}) where {A<:AlgebraType} = false  # Default fallback


# =============================================================================
# Accessor Functions
# =============================================================================

"""
    Base.getproperty(pop::PolyOpt{A,P}, sym::Symbol) where {A,P}

Custom property access for backward compatibility.

Provides access to:
- `variables`: Returns Vector{Variable} derived from registry (for backward compat)
- `is_unipotent`: Returns true if algebra type is UnipotentAlgebra
- `is_projective`: Returns true if algebra type is ProjectorAlgebra
- Standard fields: objective, eq_constraints, ineq_constraints, registry
"""
function Base.getproperty(pop::PolyOpt{A,P}, sym::Symbol) where {A<:AlgebraType, P}
    if sym === :variables
        # Backward compatibility: derive Variable vector from registry
        reg = getfield(pop, :registry)
        syms = symbols(reg)
        T = index_type(reg)
        return [Variable(s, false, Int(reg[s])) for s in syms]
    elseif sym === :is_unipotent
        return A === UnipotentAlgebra
    elseif sym === :is_projective
        return A === ProjectorAlgebra
    elseif sym === :comm_gps
        # Backward compatibility: return all variables as a single commutative group
        # (actual commutation is handled by algebra type)
        reg = getfield(pop, :registry)
        syms = symbols(reg)
        T = index_type(reg)
        vars = [Variable(s, false, Int(reg[s])) for s in syms]
        return [vars]
    else
        return getfield(pop, sym)
    end
end

# Prevent setproperty for computed properties
function Base.setproperty!(pop::PolyOpt, sym::Symbol, value)
    if sym in (:variables, :is_unipotent, :is_projective, :comm_gps)
        error("Cannot set computed property :$sym")
    end
    setfield!(pop, sym, value)
end

# Define propertynames for tab-completion and introspection
function Base.propertynames(::PolyOpt, private::Bool=false)
    return (:objective, :eq_constraints, :ineq_constraints, :registry,
            :variables, :is_unipotent, :is_projective, :comm_gps)
end


# =============================================================================
# Backward Compatibility (temporary - will be removed in Phase 4)
# =============================================================================

"""
    ComplexPolyOpt{A, P}

Type alias for `PolyOpt{A, P}`.

In the legacy API, `ComplexPolyOpt` was a separate type for optimization problems
with complex coefficients or non-Hermitian products. In the new design, `PolyOpt{A, P}`
handles both cases, with the algebra type `A` determining the appropriate behavior.

# Deprecated
Use `PolyOpt` directly for new code. This alias will be removed in Phase 4.

See also: [`PolyOpt`](@ref)
"""
const ComplexPolyOpt = PolyOpt

"""
    cpolyopt(objective, registry; kwargs...)

Backward compatibility alias for `polyopt`.

This function exists to ease migration from the legacy API where `cpolyopt` was
used for complex polynomial optimization problems. In the new design, `polyopt`
handles both real and complex problems based on the algebra type and coefficient type.

# Deprecated
Use `polyopt` directly for new code. This alias will be removed in Phase 4.

See also: [`polyopt`](@ref)
"""
const cpolyopt = polyopt
