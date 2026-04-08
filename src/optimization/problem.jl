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
pop = polyopt(sp, reg)  # promoted internally to NCStatePolynomial
```

See also: [`polyopt`](@ref), [`VariableRegistry`](@ref), [`AlgebraType`](@ref)
"""
struct PolyOpt{A<:AlgebraType,T<:Integer,P<:AbstractPolynomial} <: OptimizationProblem{A,P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    moment_eq_constraints::Vector{P}
    registry::VariableRegistry{A,T}
end

@inline _is_hermitian_polynomial(poly::Polynomial) = iszero(poly - adjoint(poly))
@inline _validate_polyopt_problem_data(::AbstractPolynomial, _eq_constraints, _ineq_constraints, _moment_eq_constraints) = nothing

function _validate_polyopt_problem_data(
    objective::Polynomial{A,T,C},
    _eq_constraints,
    ineq_constraints,
    _moment_eq_constraints,
) where {A<:AlgebraType,T<:Integer,C<:Number}
    _is_complex_problem(A) || return nothing

    _is_hermitian_polynomial(objective) || throw(ArgumentError(
        "Complex-algebra objectives must be Hermitian. Got a non-Hermitian objective for $(nameof(A))."
    ))

    for (i, poly) in pairs(ineq_constraints)
        _is_hermitian_polynomial(poly) || throw(ArgumentError(
            "Inequality constraint $i must be Hermitian for $(nameof(A)) problems. " *
            "Non-Hermitian PSD constraints are not currently supported."
        ))
    end

    return nothing
end

"""
    polyopt(objective::P, registry::VariableRegistry{A,T};
            eq_constraints=P[], ineq_constraints=P[]) where {P<:AbstractPolynomial, A, T}

Create a polynomial optimization problem from objective, registry, and optional constraints.

Works with any `AbstractPolynomial` subtype including `Polynomial{A,T,C}` and
`NCStatePolynomial{C,ST,A,T}`. For state/trace problems over `MonoidAlgebra`,
a `StatePolynomial{C,ST,A,T}` objective is accepted as a convenience and is
promoted to the equivalent `NCStatePolynomial{C,ST,A,T}` with identity
non-commutative words.

# Arguments
- `objective::P`: The polynomial objective function to optimize
- `registry::VariableRegistry{A,T}`: Variable registry for the algebra

# Keyword Arguments
- `eq_constraints`: Equality constraints as polynomials (p = 0). Default: empty
- `ineq_constraints`: Inequality constraints as polynomials (p >= 0). Default: empty
- `moment_eq_constraints`: One-sided moment equality constraints encoding
  state conditions `g|ψ⟩ = 0` via `⟨b† g⟩ = 0` inside the relaxation. Default: empty.
  Currently supported for ordinary polynomial optimization, not state/trace problems.

# Returns
A `PolyOpt{A,T,P}` structure representing the optimization problem.
For state-polynomial inputs, the stored polynomial type is the promoted
`NCStatePolynomial` form.

# Notes
- Algebra type `A` is inferred from the registry
- Coefficient type cannot be an integer subtype (JuMP solver requirement)
- Simplification rules are determined by the algebra type, not manual flags
- For complex-algebra polynomial problems (`PauliAlgebra`, `FermionicAlgebra`, `BosonicAlgebra`),
  the objective and inequality constraints must be Hermitian. Non-Hermitian PSD constraints are
  rejected.
- For `FermionicAlgebra`: objectives should have even parity (parity superselection rule).
  Odd-parity operators have zero expectation value. Validation is done during moment
  relaxation via `_add_parity_constraints!`.
- For `NCStatePolynomial` objectives: all `NCStateWord` terms in the objective must have
  identity NC words (i.e., pure expectation-value expressions). Non-identity NC words in
  the objective cause an `ArgumentError`. Constraints may have non-identity NC words.
  Use [`expect`](@ref) or [`tr`](@ref) to construct valid state-polynomial objectives.

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
pop = polyopt(sp, reg)  # promoted internally to NCStatePolynomial
```

See also: [`PolyOpt`](@ref), [`VariableRegistry`](@ref), [`NCStatePolynomial`](@ref)
"""
function polyopt(
    objective::P,
    registry::VariableRegistry{A,T};
    eq_constraints::Vector{P}=P[],
    ineq_constraints::Vector{P}=P[],
    moment_eq_constraints::Vector{P}=P[]
) where {P<:AbstractPolynomial, A<:AlgebraType, T<:Integer}
    C = coeff_type(P)
    if C <: Integer
        throw(ArgumentError("Polynomial coefficients cannot be integers (not supported by JuMP solvers). Use Float64 or other floating-point types."))
    end

    _validate_polyopt_problem_data(objective, eq_constraints, ineq_constraints, moment_eq_constraints)

    # Deduplicate constraints
    eq_cons = unique!(copy(eq_constraints))
    ineq_cons = unique!(copy(ineq_constraints))
    meq_cons = unique!(copy(moment_eq_constraints))

    return PolyOpt{A,T,P}(objective, eq_cons, ineq_cons, meq_cons, registry)
end

@inline function _lift_to_nc_state_poly(
    poly::StatePolynomial{C,ST,A,T}
) where {C<:Number,ST<:StateType,A<:MonoidAlgebra,T<:Integer}
    nc_words = NCStateWord{ST,A,T}[
        NCStateWord(sw, one(NormalMonomial{A,T})) for sw in monomials(poly)
    ]
    return NCStatePolynomial(copy(coefficients(poly)), nc_words)
end

@inline _lift_to_nc_state_constraint(::Type{NCStatePolynomial{C,ST,A,T}}, poly::StatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:MonoidAlgebra,T<:Integer} =
    _lift_to_nc_state_poly(poly)
@inline _lift_to_nc_state_constraint(::Type{NCStatePolynomial{C,ST,A,T}}, poly::NCStatePolynomial{C,ST,A,T}) where {C<:Number,ST<:StateType,A<:MonoidAlgebra,T<:Integer} =
    poly

function _lift_to_nc_state_constraint(::Type{NCStatePolynomial{C,ST,A,T}}, poly) where {C<:Number,ST<:StateType,A<:MonoidAlgebra,T<:Integer}
    throw(ArgumentError("State-polynomial constraints must be `StatePolynomial{$C,$ST,$A,$T}` or `NCStatePolynomial{$C,$ST,$A,$T}`; got $(typeof(poly))."))
end

function polyopt(
    objective::StatePolynomial{C,ST,A,T},
    registry::VariableRegistry{A,T};
    eq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[],
    ineq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[],
    moment_eq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[]
) where {C<:Number,ST<:StateType,A<:MonoidAlgebra,T<:Integer}
    NCSP = NCStatePolynomial{C,ST,A,T}
    return polyopt(
        _lift_to_nc_state_poly(objective),
        registry;
        eq_constraints=NCSP[
            _lift_to_nc_state_constraint(NCSP, con) for con in eq_constraints
        ],
        ineq_constraints=NCSP[
            _lift_to_nc_state_constraint(NCSP, con) for con in ineq_constraints
        ],
        moment_eq_constraints=NCSP[
            _lift_to_nc_state_constraint(NCSP, con) for con in moment_eq_constraints
        ]
    )
end

function polyopt(
    objective::NCStatePolynomial{C,ST,A,T},
    registry::VariableRegistry{A,T};
    eq_constraints::Vector{NCStatePolynomial{C,ST,A,T}}=NCStatePolynomial{C,ST,A,T}[],
    ineq_constraints::Vector{NCStatePolynomial{C,ST,A,T}}=NCStatePolynomial{C,ST,A,T}[],
    moment_eq_constraints::Vector{NCStatePolynomial{C,ST,A,T}}=NCStatePolynomial{C,ST,A,T}[]
) where {C<:Number,ST<:StateType,A<:MonoidAlgebra,T<:Integer}
    _validate_identity_nc_words(objective, "objective")
    isempty(moment_eq_constraints) || throw(ArgumentError(
        "moment_eq_constraints are not yet supported for state polynomial optimization."
    ))
    return invoke(polyopt, Tuple{AbstractPolynomial, VariableRegistry{A,T}},
        objective, registry; eq_constraints=eq_constraints, ineq_constraints=ineq_constraints, moment_eq_constraints=moment_eq_constraints)
end

"""
    _validate_identity_nc_words(poly::NCStatePolynomial, label::AbstractString)

Validate that all `NCStateWord` terms in `poly` have identity NC words.

State/trace polynomial optimization objectives must be pure expectation-value
expressions (i.e., each term is of the form `c * ⟨state_word⟩`). Non-identity
NC word factors indicate an ill-formed objective. Use [`expect`](@ref) or
[`tr`](@ref) to construct valid state-polynomial objectives.
"""
function _validate_identity_nc_words(poly::NCStatePolynomial, label::AbstractString)
    for ncsw in monomials(poly)
        isone(ncsw.nc_word) || throw(ArgumentError(
            "NCStatePolynomial $label contains a term with non-identity NC word " *
            "`$(ncsw.nc_word)`. State/trace polynomial objectives must be pure " *
            "expectation-value expressions where all NC words are identity. " *
            "Use `expect(...)` or `tr(...)` to construct valid objectives."
        ))
    end
end

function polyopt(
    objective::StatePolynomial{C,ST,A,T},
    registry::VariableRegistry{A,T};
    eq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[],
    ineq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[],
    moment_eq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[]
) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    throw(ArgumentError("`polyopt(::StatePolynomial, ...)` is only supported for state/trace problems over `MonoidAlgebra`; got `$(nameof(A))`. This algebra is not supported by the state-polynomial optimization pipeline."))
end

function polyopt(
    objective::StatePolynomial{C,ST,A,T},
    registry::VariableRegistry;
    eq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[],
    ineq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[],
    moment_eq_constraints::AbstractVector=StatePolynomial{C,ST,A,T}[]
) where {C<:Number,ST<:StateType,A<:AlgebraType,T<:Integer}
    throw(ArgumentError("`polyopt(::StatePolynomial, ...)` requires the objective and registry to use the same algebra and index types; got objective type `$(typeof(objective))` with registry type `$(typeof(registry))`."))
end

@inline function _newton_chip_single_site_objective(
    poly::Polynomial{NonCommutativeAlgebra,T,C}
) where {T<:Unsigned,C<:Number}
    return length(unique(decode_site.(collect(variable_indices(poly))))) <= 1
end

@inline function _newton_chip_single_site_objective(
    poly::NCStatePolynomial{C,ST,NonCommutativeAlgebra,T}
) where {C<:Number,ST<:StateType,T<:Unsigned}
    return length(unique(decode_site.(collect(variable_indices(poly))))) <= 1
end

@inline function _newton_cyclic_trace_word(
    ncsw::NCStateWord{MaxEntangled,NonCommutativeAlgebra,T}
) where {T<:Unsigned}
    isone(ncsw.nc_word) ||
        throw(ArgumentError("`newton_chip_basis` for tracial problems only supports scalar trace objectives of the form `tr(f)`; found a basis term with a residual nc-word factor. Pass a scalar trace objective and inject the resulting basis through `moment_basis`."))

    sw = ncsw.sw
    isone(sw) && return one(NormalMonomial{NonCommutativeAlgebra,T})
    length(sw.state_syms) == 1 ||
        throw(ArgumentError("`newton_chip_basis` for tracial problems only supports single-trace objectives `tr(f)`. Products of traces and general state-polynomial terms are outside the Newton cyclic chip theorem implemented here."))

    return NormalMonomial{NonCommutativeAlgebra,T}(copy(only(sw.state_syms).mono))
end

function _newton_cyclic_support_exponents(
    poly::NCStatePolynomial{C,MaxEntangled,NonCommutativeAlgebra,T}
) where {C<:Number,T<:Unsigned}
    used_indices = sort!(collect(variable_indices(poly)))
    isempty(used_indices) && return used_indices, Tuple{Vararg{Int}}[()]

    index_to_pos = Dict(idx => i for (i, idx) in enumerate(used_indices))
    support_exponents = Set{Tuple{Vararg{Int}}}()
    push!(support_exponents, Tuple(fill(0, length(used_indices))))

    for ncsw in monomials(poly)
        trace_word = _newton_cyclic_trace_word(ncsw)
        exponents = zeros(Int, length(used_indices))
        for idx in trace_word.word
            exponents[index_to_pos[idx]] += 1
        end
        push!(support_exponents, Tuple(exponents))
    end

    return used_indices, sort!(collect(support_exponents))
end

function _newton_cyclic_half_polytope_oracle(support_exponents::Vector{<:Tuple})
    nvars = length(first(support_exponents))
    nvars == 0 && return (_target_exp::Tuple) -> true

    support_set = Set(support_exponents)
    m = length(support_exponents)
    support_matrix = Matrix{Float64}(undef, nvars, m)
    coordwise_max = zeros(Int, nvars)

    for (j, exponent) in pairs(support_exponents)
        for i in 1:nvars
            support_matrix[i, j] = exponent[i]
            coordwise_max[i] = max(coordwise_max[i], exponent[i])
        end
    end

    model = JuMP.Model(Clarabel.Optimizer)
    JuMP.set_silent(model)
    @variable(model, λ[1:m] >= 0)
    @constraint(model, sum(λ) == 1)
    coord_constraints = [
        @constraint(model, sum(support_matrix[i, j] * λ[j] for j in 1:m) == 0.0)
        for i in 1:nvars
    ]
    @objective(model, Min, 0.0)

    function oracle(target_exp::Tuple{Vararg{Int}})
        all(iszero, target_exp) && return true

        doubled_target = ntuple(i -> 2 * target_exp[i], nvars)
        doubled_target in support_set && return true
        any(doubled_target[i] > coordwise_max[i] for i in 1:nvars) && return false

        for i in 1:nvars
            JuMP.set_normalized_rhs(coord_constraints[i], float(doubled_target[i]))
        end

        JuMP.optimize!(model)
        status = JuMP.termination_status(model)
        status in (JuMP.MOI.OPTIMAL, JuMP.MOI.ALMOST_OPTIMAL) && return true
        status in (JuMP.MOI.INFEASIBLE, JuMP.MOI.INFEASIBLE_OR_UNBOUNDED) && return false

        throw(ErrorException("Newton cyclic chip polytope membership solve failed with status $(status)."))
    end

    return oracle
end

function _newton_cyclic_chip_basis(
    poly::NCStatePolynomial{C,MaxEntangled,NonCommutativeAlgebra,T},
    registry::VariableRegistry{NonCommutativeAlgebra,T},
    d::Int
) where {C<:Number,T<:Unsigned}
    d < 0 && throw(DomainError(d, "`d` must be non-negative."))

    used_indices, support_exponents = _newton_cyclic_support_exponents(poly)
    M = NCStateWord{MaxEntangled,NonCommutativeAlgebra,T}
    isempty(used_indices) && return [one(M)]

    sub_reg = subregistry(registry, used_indices)
    dense_words = get_ncbasis(sub_reg, d)
    index_to_pos = Dict(idx => i for (i, idx) in enumerate(used_indices))
    in_half_polytope = _newton_cyclic_half_polytope_oracle(support_exponents)
    exponent_cache = Dict{Tuple{Vararg{Int}},Bool}()
    sw_identity = one(StateWord{MaxEntangled,NonCommutativeAlgebra,T})

    basis = M[]
    sizehint!(basis, length(dense_words))

    for mono in dense_words
        exponents = zeros(Int, length(used_indices))
        for idx in mono.word
            exponents[index_to_pos[idx]] += 1
        end
        exponent_key = Tuple(exponents)
        keep = get!(exponent_cache, exponent_key) do
            in_half_polytope(exponent_key)
        end
        keep && push!(basis, NCStateWord(sw_identity, mono))
    end

    return basis
end

"""
    newton_chip_basis(pop::PolyOpt, d::Int)

Construct a strict Newton-chip basis for an unconstrained ordinary polynomial
optimization problem over `NonCommutativeAlgebra`.

The returned basis contains the identity and every right chip (suffix) of each
word `w` such that `w' * w` appears in the support of `pop.objective`, with
chips truncated to degree at most `d`. The intended integration path is
`SolverConfig(moment_basis=newton_chip_basis(pop, d))`.

Warning: the classical Newton-chip theorem is stated for the free
noncommutative `*`-algebra with no extra commuting relations. In this package,
`NonCommutativeAlgebra` commutes operators across physical sites, so this helper
rejects objectives whose support spans more than one site.
"""
function newton_chip_basis(
    pop::PolyOpt{NonCommutativeAlgebra,T,P},
    d::Int
) where {T<:Unsigned,C<:Number,P<:Polynomial{NonCommutativeAlgebra,T,C}}
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) ||
        throw(ArgumentError("`newton_chip_basis` is only supported for unconstrained `PolyOpt` problems. Pass the resulting basis through `moment_basis` only for constraint-free problems."))
    _newton_chip_single_site_objective(pop.objective) ||
        throw(ArgumentError("`newton_chip_basis` assumes no extra commuting relations. In this package, `NonCommutativeAlgebra` commutes operators across physical sites, so the helper currently requires the objective support to use variables from at most one site."))

    return _newton_chip_basis(pop.objective, d)
end

"""
    newton_chip_basis(pop::PolyOpt, d::Int)

Construct a strict Newton cyclic chip basis for an unconstrained tracial
polynomial optimization problem over `NonCommutativeAlgebra`.

This method applies to scalar trace objectives of the form `tr(f)`, represented
in the package as `NCStatePolynomial{<:Number,MaxEntangled,NonCommutativeAlgebra}`.
The returned basis contains pure operator basis elements `<I>⊗w` whose
commutative-collapse exponent lies in half the tracial Newton polytope of the
canonical trace support carried by `pop.objective`. The intended integration
path is `SolverConfig(moment_basis=newton_chip_basis(pop, d))`. As with any
custom state/trace `moment_basis`, the injected basis must still generate every
objective/constraint moment required by the relaxation; underspecified choices
now error instead of silently dropping terms.

Warning: the Newton cyclic chip theorem is an unconstrained theorem for the
free noncommutative `*`-algebra. This helper therefore rejects constrained
problems, objectives spanning more than one site, products of traces, mixed
state/operator terms, and algebras with built-in relations beyond the free
single-site setting.
"""
function newton_chip_basis(
    pop::PolyOpt{NonCommutativeAlgebra,T,P},
    d::Int
) where {T<:Unsigned,C<:Number,P<:NCStatePolynomial{C,MaxEntangled,NonCommutativeAlgebra,T}}
    isempty(pop.eq_constraints) && isempty(pop.ineq_constraints) ||
        throw(ArgumentError("`newton_chip_basis` is only supported for unconstrained `PolyOpt` problems. Pass the resulting basis through `moment_basis` only for constraint-free problems."))
    _newton_chip_single_site_objective(pop.objective) ||
        throw(ArgumentError("`newton_chip_basis` assumes no extra commuting relations. In this package, `NonCommutativeAlgebra` commutes operators across physical sites, so the helper currently requires the objective support to use variables from at most one site."))

    return _newton_cyclic_chip_basis(pop.objective, pop.registry, d)
end

function newton_chip_basis(pop::PolyOpt{A,T,P}, _d::Int) where {A<:AlgebraType,T<:Integer,C<:Number,P<:Polynomial{A,T,C}}
    throw(ArgumentError("`newton_chip_basis` is only supported for ordinary polynomial `PolyOpt` problems over `NonCommutativeAlgebra`; got `$(nameof(A))`."))
end

function newton_chip_basis(pop::PolyOpt{NonCommutativeAlgebra,T,P}, _d::Int) where {T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,NonCommutativeAlgebra,T}}
    throw(ArgumentError("`newton_chip_basis` is only supported for tracial `PolyOpt` problems with `MaxEntangled` state type over `NonCommutativeAlgebra`; got `$(nameof(ST))`."))
end

function newton_chip_basis(pop::PolyOpt{A,T,P}, _d::Int) where {A<:AlgebraType,T<:Integer,ST<:StateType,C<:Number,P<:NCStatePolynomial{C,ST,A,T}}
    throw(ArgumentError("`newton_chip_basis` is only supported for unconstrained ordinary polynomial problems, or for unconstrained tracial `PolyOpt` problems over `NonCommutativeAlgebra`; got state type `$(nameof(ST))` over `$(nameof(A))`."))
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

        Moment equality constraints ($(length(pop.moment_eq_constraints))):
            $(isempty(pop.moment_eq_constraints) ? "(none)" : cons_str(pop.moment_eq_constraints, true))

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
# @noinline: prevents the compiler from inlining these trait methods away,
# keeping them visible to Julia's code-coverage instrumentation.
@noinline _is_complex_problem(::Type{<:Union{PauliAlgebra,FermionicAlgebra,BosonicAlgebra}}) = true
@noinline _is_complex_problem(::Type{<:Union{NonCommutativeAlgebra,ProjectorAlgebra,UnipotentAlgebra}}) = false
