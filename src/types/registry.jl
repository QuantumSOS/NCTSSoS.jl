"""
    VariableRegistry

A registry that maps between variable symbols and their indices.
Supports both contiguous (Vector-like) and sparse (Dict-based) indices.

# Fields
- `idx_to_variables::Dict{T, Symbol}`: Maps index to symbol
- `variables_to_idx::Dict{Symbol, T}`: Maps symbol to index

# Examples
```jldoctest
julia> reg, _ = create_pauli_variables(1:2);

julia> reg
VariableRegistry with 6 variables: σx₁, σy₁, σz₁, σx₂, σy₂, σz₂

julia> length(reg)
6

julia> reg[1]  # Access by index
:σx₁

julia> Int(reg[:σx₁])  # Access by symbol (convert typed index to Int for display)
1

julia> :σx₁ in reg
true
```
"""
struct VariableRegistry{A<:AlgebraType, T<:Integer}
    idx_to_variables::Dict{T, Symbol}
    variables_to_idx::Dict{Symbol,T}

    # Inner constructor to enforce invariants
    function VariableRegistry{A, T}(idx_to_vars::Dict{T, Symbol}, vars_to_idx::Dict{Symbol,T}) where {A<:AlgebraType, T<:Integer}
        # Verify consistency
        length(idx_to_vars) == length(vars_to_idx) ||
            error("Inconsistent registry: different lengths")

        for (idx, sym) in idx_to_vars
            haskey(vars_to_idx, sym) && vars_to_idx[sym] == idx ||
                error("Inconsistent registry: symbol $sym has wrong index")
        end

        new{A, T}(idx_to_vars, vars_to_idx)
    end
end

# =============================================================================
# Display Registry: module-level fallback for REPL display of monomials
# =============================================================================

"""
    _DISPLAY_REGISTRIES::Dict{DataType, VariableRegistry}

Module-level cache mapping `AlgebraType` to the most recently created
`VariableRegistry` for that algebra.  Used by `Base.show(::IO, ::NormalMonomial)`
as a fallback when no `:registry` IOContext is present, so that monomials
print with human-readable symbol names in the REPL.

Updated automatically by every `create_*_variables` call.  The IOContext
`:registry` key always takes precedence when supplied.
"""
const _DISPLAY_REGISTRIES = Dict{DataType, Any}()

"""
    _register_display_registry(reg::VariableRegistry{A}) where {A}

Store `reg` as the display-fallback registry for algebra type `A`.
"""
function _register_display_registry(reg::VariableRegistry{A}) where {A<:AlgebraType}
    _DISPLAY_REGISTRIES[A] = reg
    return nothing
end

"""
    _get_display_registry(::Type{A}) where {A<:AlgebraType} -> Union{VariableRegistry, Nothing}

Retrieve the display-fallback registry for algebra type `A`, or `nothing`.
"""
function _get_display_registry(::Type{A}) where {A<:AlgebraType}
    return get(_DISPLAY_REGISTRIES, A, nothing)
end

"""
    subregistry(reg::VariableRegistry{A,T}, subset_indices::VT) where {A,T,VT<:AbstractVector{T}}

Create a new VariableRegistry containing only the specified subset of indices.

This is useful for creating clique-local registries for basis generation in
correlative sparsity decomposition.

# Arguments
- `reg::VariableRegistry{A,T}`: The parent registry
- `subset_indices::AbstractVector{T}`: Indices to include in the sub-registry (must match registry's index type)

# Returns
A new `VariableRegistry{A,T}` containing only the variables at the specified indices.

# Examples
```jldoctest
julia> reg, (x,) = create_noncommutative_variables([("x", 1:5)]);

julia> sorted_idxs = sort(collect(keys(reg.idx_to_variables)));

julia> sub_reg = subregistry(reg, sorted_idxs[1:3]);

julia> length(sub_reg)
3
```

# Notes
- Indices not present in the parent registry are silently ignored
- The returned registry has the same algebra type as the parent
- This is a copy operation, not a view - modifications to the sub-registry
  do not affect the parent registry
"""
function subregistry(reg::VariableRegistry{A,T}, subset_indices::VT) where {A<:AlgebraType, T<:Integer, VT<:AbstractVector{T}}
    filtered_idx_to_vars = Dict{T, Symbol}()
    filtered_vars_to_idx = Dict{Symbol, T}()

    for idx in subset_indices
        if haskey(reg.idx_to_variables, idx)
            sym = reg.idx_to_variables[idx]
            filtered_idx_to_vars[idx] = sym
            filtered_vars_to_idx[sym] = idx
        end
    end

    return VariableRegistry{A,T}(filtered_idx_to_vars, filtered_vars_to_idx)
end

"""
    _basis_subregistry(reg::VariableRegistry{A,T}, subset_indices)

Create the registry used for automatic moment-basis generation on a clique.

For `PBWAlgebra`s, correlative sparsity is tracked in terms of physical modes,
so a clique over modes `{1,2}` must generate basis words over both annihilation
and creation operators `{a₁, a₂, a₁†, a₂†}`.  A plain `subregistry(reg, [1, 2])`
would drop the negative indices and silently build an incomplete basis.
"""
function _basis_subregistry(reg::VariableRegistry{A,T}, subset_indices::VT) where {A<:AlgebraType,T<:Integer,VT<:AbstractVector{T}}
    return subregistry(reg, subset_indices)
end

function _basis_subregistry(reg::VariableRegistry{A,T}, subset_indices::VT) where {A<:PBWAlgebra,T<:Signed,VT<:AbstractVector{T}}
    modes = sort!(unique(T(abs(idx)) for idx in subset_indices))
    expanded_indices = T[]
    sizehint!(expanded_indices, 2 * length(modes))

    for mode in modes
        creation_idx = -mode
        annihilation_idx = mode
        haskey(reg.idx_to_variables, creation_idx) && push!(expanded_indices, creation_idx)
        haskey(reg.idx_to_variables, annihilation_idx) && push!(expanded_indices, annihilation_idx)
    end

    return subregistry(reg, expanded_indices)
end

# Unicode subscript digits: ₀₁₂₃₄₅₆₇₈₉
const SUBSCRIPT_DIGITS = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']

"""
    _check_no_duplicate_subscripts(subscripts)

Throw `ArgumentError` if `subscripts` contains duplicate values.
"""
function _check_no_duplicate_subscripts(subscripts)
    seen = Set{eltype(subscripts)}()
    for s in subscripts
        if s in seen
            throw(ArgumentError("Duplicate subscript value: $s. Each subscript must be unique."))
        end
        push!(seen, s)
    end
end

"""
    _subscript_string(n::Integer) -> String

Convert an integer to its Unicode subscript representation.

# Example
```julia
_subscript_string(42)  # "₄₂"
_subscript_string(0)   # "₀"
```
"""
function _subscript_string(n::Integer)

    if n < 0
        error("Negative subscripts not supported: $n")
    end

    # Convert to string of digits
    digits_str = string(Int(n))

    # Map each digit to subscript
    return String([SUBSCRIPT_DIGITS[parse(Int, d)+1] for d in digits_str])
end

"""
    Base.length(reg::VariableRegistry) -> Int

Return the number of variables in the registry.
"""
Base.length(reg::VariableRegistry) = length(reg.idx_to_variables)

"""
    Base.getindex(reg::VariableRegistry{A,T}, idx::T) -> Symbol

Get the variable symbol at the given index.
"""
Base.getindex(reg::VariableRegistry{A,T}, idx::T) where {A,T} = reg.idx_to_variables[idx]

"""
    Base.getindex(reg::VariableRegistry{A,T}, idx::Integer) -> Symbol

Get the variable symbol at the given index (converts to type T).
"""
Base.getindex(reg::VariableRegistry{A,T}, idx::Integer) where {A,T} = reg.idx_to_variables[T(idx)]

"""
    Base.getindex(reg::VariableRegistry, sym::Symbol) -> Int

Get the index of the given variable symbol.
"""
Base.getindex(reg::VariableRegistry, sym::Symbol) = reg.variables_to_idx[sym]

"""
    Base.in(sym::Symbol, reg::VariableRegistry) -> Bool

Check if a variable symbol exists in the registry.
"""
Base.in(sym::Symbol, reg::VariableRegistry) = haskey(reg.variables_to_idx, sym)

"""
    symbols(reg::VariableRegistry) -> Vector{Symbol}

Return all symbols in index-sorted order.

# Examples
```jldoctest
julia> reg, _ = create_projector_variables([("P", 1:3)]);

julia> reg
VariableRegistry with 3 variables: P₁, P₂, P₃

julia> symbols(reg)
3-element Vector{Symbol}:
 :P₁
 :P₂
 :P₃
```
"""
function symbols(reg::VariableRegistry)
    sorted_pairs = sort(collect(reg.idx_to_variables), by=first)
    return [p.second for p in sorted_pairs]
end

"""
    indices(reg::VariableRegistry{A,T}) -> Vector{T}

Return all indices in sorted order.

# Examples
```jldoctest
julia> reg, _ = create_projector_variables([("P", 1:3)]);

julia> reg
VariableRegistry with 3 variables: P₁, P₂, P₃

julia> indices(reg) == UInt8[0x05, 0x09, 0x0d]
true
```
"""
indices(reg::VariableRegistry) = sort(collect(keys(reg.idx_to_variables)))

"""
    Base.show(io::IO, reg::VariableRegistry)

Display a variable registry showing its variables in order.
"""
function Base.show(io::IO, reg::VariableRegistry)
    n = length(reg)
    print(io, "VariableRegistry with $n variable")
    n != 1 && print(io, "s")

    if n > 0
        # Sort by index for consistent display
        sorted_pairs = sort(collect(reg.idx_to_variables), by=first)
        syms = [p.second for p in sorted_pairs]

        if n <= 10
            print(io, ": ", join(syms, ", "))
        else
            print(io, ": ", join(syms[1:5], ", "),
                ", ..., ", join(syms[end-2:end], ", "))
        end
    end
end

"""
    create_pauli_variables(subscripts)

Create a variable registry and monomials for Pauli spin matrices σx, σy, σz.

Pauli matrices satisfy:
- σᵢ² = 1 (each Pauli matrix squares to identity)
- {σᵢ, σⱼ} = 2δᵢⱼ (anticommutation relation)

The operators are always x, y, z components. This function creates variables
for each subscript value.

**Variable Ordering:** Variables are ordered by site first, then by Pauli type (x, y, z).
This enables efficient encoding: for index `idx`, site = `(idx - 1) ÷ 3 + 1`,
pauli_type = `(idx - 1) % 3` (0=X, 1=Y, 2=Z).

# Arguments
- `subscripts`: Subscript values for different qubits/sites

# Returns
A tuple of:
- `VariableRegistry`: Registry with Pauli variables σx, σy, σz for each subscript
- `(σx, σy, σz)`: Tuple of monomial vectors grouped by Pauli type

# Examples
```jldoctest
julia> reg, (σx, σy, σz) = create_pauli_variables(1:2);

julia> length(σx)  # One σx per site
2

julia> typeof(σx[1])  # monomial for σx₁
NormalMonomial{PauliAlgebra, UInt8}

julia> :σx₁ in reg
true

julia> Int(reg[:σx₁])  # Index 1: site 1, type X
1
```

Single qubit:
```jldoctest
julia> reg, (σx, σy, σz) = create_pauli_variables([1]);

julia> length(σx)
1
```
"""
function create_pauli_variables(subscripts)
    _check_no_duplicate_subscripts(subscripts)
    # Pauli operators ordered by site first, then by type (x, y, z)
    # This enables encoding: site = (idx-1) ÷ 3 + 1, type = (idx-1) % 3
    pauli_types = [:σx, :σy, :σz]
    all_symbols = Symbol[]
    n_sites = length(subscripts)

    for subscript in subscripts
        subscript_str = _subscript_string(subscript)
        for pauli in pauli_types
            push!(all_symbols, Symbol(string(pauli) * subscript_str))
        end
    end

    # Select smallest unsigned type that fits 3*n_sites indices
    T = _select_pauli_type(n_sites)

    # Create Dict-based bidirectional mapping (NO sorting - preserve site-first order)
    idx_to_vars = Dict{T, Symbol}()
    variables_to_idx = Dict{Symbol, T}()
    for (idx, sym) in enumerate(all_symbols)
        idx_to_vars[T(idx)] = sym
        variables_to_idx[sym] = T(idx)
    end

    reg = VariableRegistry{PauliAlgebra, T}(idx_to_vars, variables_to_idx)
    _register_display_registry(reg)

    # Build monomial vectors grouped by Pauli type (x, y, z)
    σx = Vector{NormalMonomial{PauliAlgebra,T}}(undef, n_sites)
    σy = Vector{NormalMonomial{PauliAlgebra,T}}(undef, n_sites)
    σz = Vector{NormalMonomial{PauliAlgebra,T}}(undef, n_sites)

    for (site_idx, subscript) in enumerate(subscripts)
        subscript_str = _subscript_string(subscript)
        σx[site_idx] = NormalMonomial{PauliAlgebra,T}([reg[Symbol("σx" * subscript_str)]])
        σy[site_idx] = NormalMonomial{PauliAlgebra,T}([reg[Symbol("σy" * subscript_str)]])
        σz[site_idx] = NormalMonomial{PauliAlgebra,T}([reg[Symbol("σz" * subscript_str)]])
    end

    return (reg, (σx, σy, σz))
end

"""
    _select_pauli_type(n_sites::Int) -> Type{<:Unsigned}

Select smallest unsigned type for Pauli variables.
Pauli uses contiguous indices: 3 operators per site (σx, σy, σz).
Total indices = 3 * n_sites.
"""
function _select_pauli_type(n_sites::Int)
    n_total = 3 * n_sites
    n_total <= typemax(UInt8) && return UInt8
    n_total <= typemax(UInt16) && return UInt16
    n_total <= typemax(UInt32) && return UInt32
    return UInt64
end

"""
    create_fermionic_variables(subscripts)
    create_fermionic_variables(prefix_subscripts::Vector{Tuple{String,VT}}) where {VT<:AbstractVector}

Create a variable registry and monomials for fermionic creation (a⁺) and annihilation (a) operators.

Fermionic operators satisfy anticommutation relations:
- {aᵢ, aⱼ⁺} = δᵢⱼ (creation-annihilation)
- {aᵢ, aⱼ} = 0 (annihilation-annihilation)
- {aᵢ⁺, aⱼ⁺} = 0 (creation-creation)

# Arguments
- `subscripts`: subscript values for a single operator family named `a`
- `prefix_subscripts`: vector of `(prefix, subscripts)` tuples for multiple
  fermionic species sharing one registry

# Returns
- single-family call: `(reg, (a, a⁺))`
- multi-family call: `(reg, ((a, a⁺), (b, b⁺), ...))`

!!! note "Cross-species anticommutation"
    When multiple species are created (multi-family call), **all** modes
    share one fermionic algebra and **anticommute** with each other —
    including across species: `{c_upᵢ, c_dnⱼ⁺} = 0` for i ≠ j.
    Each tuple provides a naming prefix, not an independent algebra.
    This differs from `MonoidAlgebra` constructors like
    [`create_projector_variables`](@ref), where different tuples
    correspond to different *sites* whose operators commute.

# Examples
```jldoctest
julia> reg, (a, a⁺) = create_fermionic_variables(1:2);

julia> length(a)
2

julia> :a₁ in reg
true
```

Multiple fermionic species in one registry:
```jldoctest
julia> reg, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([("c_up", 1:2), ("c_dn", 1:2)]);

julia> length(c_up), length(c_dn)
(2, 2)

julia> c_up[1] == c_dn[1]
false
```

Cross-species operators anticommute (they are **not** independent):
```jldoctest
julia> reg, ((a, a_dag), (b, b_dag)) = create_fermionic_variables([("a", 1:1), ("b", 1:1)]);

julia> a[1] * b[1] + b[1] * a[1]  # {a₁, b₁} = 0
0

julia> a[1] * b_dag[1] + b_dag[1] * a[1]  # {a₁, b₁⁺} = 0 (different species)
0
```
"""
create_fermionic_variables(subscripts) = _create_physical_variables(FermionicAlgebra, subscripts, "a")

function create_fermionic_variables(
    prefix_subscripts::Vector{Tuple{String,VT}}
) where {T<:Integer,VT<:AbstractVector{T}}
    _create_physical_variables(FermionicAlgebra, prefix_subscripts)
end

function _create_physical_variables(::Type{A}, subscripts, prefix::String) where {A<:AlgebraType}
    reg, monomial_groups = _create_physical_variables(A, [(prefix, subscripts)])
    return reg, monomial_groups[1]
end

function _create_physical_variables(
    ::Type{A},
    prefix_subscripts::Vector{Tuple{String,VT}},
) where {A<:AlgebraType,T<:Integer,VT<:AbstractVector{T}}
    for (_, subscripts) in prefix_subscripts
        _check_no_duplicate_subscripts(subscripts)
    end

    n_modes = sum(length(subscripts) for (_, subscripts) in prefix_subscripts)
    Tidx = _select_signed_index_type(n_modes)

    idx_to_vars = Dict{Tidx, Symbol}()
    variables_to_idx = Dict{Symbol, Tidx}()
    monomial_groups = Vector{Tuple{Vector{NormalMonomial{A,Tidx}},Vector{NormalMonomial{A,Tidx}}}}(
        undef,
        length(prefix_subscripts),
    )

    next_mode = 0
    for (group_idx, (prefix, subscripts)) in enumerate(prefix_subscripts)
        annihilation = Vector{NormalMonomial{A,Tidx}}(undef, length(subscripts))
        creation = Vector{NormalMonomial{A,Tidx}}(undef, length(subscripts))

        for (local_idx, subscript) in enumerate(subscripts)
            next_mode += 1
            mode = Tidx(next_mode)
            subscript_str = _subscript_string(subscript)
            ann_sym = Symbol(prefix * subscript_str)
            cre_sym = Symbol(prefix * "⁺" * subscript_str)

            if haskey(variables_to_idx, ann_sym) || haskey(variables_to_idx, cre_sym)
                throw(ArgumentError("Duplicate physical variable symbol: $ann_sym or $cre_sym"))
            end

            idx_to_vars[mode] = ann_sym
            variables_to_idx[ann_sym] = mode
            idx_to_vars[-mode] = cre_sym
            variables_to_idx[cre_sym] = -mode

            annihilation[local_idx] = NormalMonomial{A,Tidx}([mode])
            creation[local_idx] = NormalMonomial{A,Tidx}([-mode])
        end

        monomial_groups[group_idx] = (annihilation, creation)
    end

    reg = VariableRegistry{A, Tidx}(idx_to_vars, variables_to_idx)
    _register_display_registry(reg)

    return reg, Tuple(monomial_groups)
end

"""
    _select_signed_index_type(n::Int) -> Type{<:Signed}

Select smallest signed type that can hold ±n (for fermionic/bosonic).
"""
function _select_signed_index_type(n::Int)
    n <= typemax(Int8) && return Int8
    n <= typemax(Int16) && return Int16
    n <= typemax(Int32) && return Int32
    return Int64
end

"""
    create_bosonic_variables(subscripts)
    create_bosonic_variables(prefix_subscripts::Vector{Tuple{String,VT}}) where {VT<:AbstractVector}

Create a variable registry and monomials for bosonic creation (c⁺) and annihilation (c) operators.

Bosonic operators satisfy commutation relations:
- [cᵢ, cⱼ⁺] = δᵢⱼ (creation-annihilation)
- [cᵢ, cⱼ] = 0 (annihilation-annihilation)
- [cᵢ⁺, cⱼ⁺] = 0 (creation-creation)

# Arguments
- `subscripts`: subscript values for a single operator family named `c`
- `prefix_subscripts`: vector of `(prefix, subscripts)` tuples for multiple
  bosonic species sharing one registry

# Returns
- single-family call: `(reg, (c, c⁺))`
- multi-family call: `(reg, ((c, c⁺), (d, d⁺), ...))`

# Examples
```jldoctest
julia> reg, (c, c⁺) = create_bosonic_variables(1:2);

julia> length(c)
2

julia> :c₁ in reg
true
```

Multiple bosonic species in one registry:
```jldoctest
julia> reg, ((b1, b1_dag), (b2, b2_dag)) = create_bosonic_variables([("b1", 1:2), ("b2", 1:2)]);

julia> length(b1), length(b2)
(2, 2)

julia> b1[1] == b2[1]
false
```

!!! note "Cross-species commutation"
    When multiple species are created, **all** modes share one bosonic
    algebra. Cross-species operators commute: `[b1ᵢ, b2ⱼ⁺] = 0`
    for any i, j (the Kronecker delta is zero across species).
    Each tuple provides a naming prefix, not an independent algebra.
"""
create_bosonic_variables(subscripts) = _create_physical_variables(BosonicAlgebra, subscripts, "c")

function create_bosonic_variables(
    prefix_subscripts::Vector{Tuple{String,VT}}
) where {T<:Integer,VT<:AbstractVector{T}}
    _create_physical_variables(BosonicAlgebra, prefix_subscripts)
end

"""
    create_projector_variables(prefix_subscripts::Vector{Tuple{String, VT}})

Create a variable registry and monomials for projector operators.

Projector operators satisfy:
- Pi^2 = Pi (idempotency: projectors square to themselves)

Projectors are self-adjoint and commutativity is NOT enforced.

# Arguments
- `prefix_subscripts`: Vector of `(prefix, subscripts)` tuples for multi-prefix creation

# Returns
A tuple of:
- `VariableRegistry`: Registry with projector variables for each subscript
- `monomial_groups`: Tuple of monomial vectors, one per prefix group

# Examples
```jldoctest
julia> reg, (P,) = create_projector_variables([("P", 1:3)]);

julia> length(P)
3

julia> typeof(P[1])  # monomial for P₁
NormalMonomial{ProjectorAlgebra, UInt8}

julia> reg, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 3:4)]);

julia> length(P), length(Q)
(2, 2)
```
"""
function create_projector_variables(
    prefix_subscripts::Vector{Tuple{String, VT}}
) where {T<:Integer, VT<:AbstractVector{T}}
    _create_noncommutative_variables(ProjectorAlgebra, prefix_subscripts)
end

# Note: Use select_uint_type from algebra.jl for site-based index type selection

"""
    create_unipotent_variables(prefix_subscripts::Vector{Tuple{String, VT}})

Create a variable registry and monomials for unipotent operators.

Unipotent operators satisfy:
- U^2 = I (squares to identity)

This is simpler than Pauli algebra - no cyclic products or
cross-operator interactions.

# Arguments
- `prefix_subscripts`: Vector of `(prefix, subscripts)` tuples for multi-prefix creation

# Returns
A tuple of:
- `VariableRegistry`: Registry with unipotent variables for each subscript
- `monomial_groups`: Tuple of monomial vectors, one per prefix group

# Examples
```jldoctest
julia> reg, (U,) = create_unipotent_variables([("U", 1:3)]);

julia> length(U)
3

julia> typeof(U[1])  # monomial for U₁
NormalMonomial{UnipotentAlgebra, UInt8}

julia> reg, (U, V) = create_unipotent_variables([("U", 1:2), ("V", 3:4)]);

julia> length(U), length(V)
(2, 2)
```
"""
function create_unipotent_variables(
    prefix_subscripts::Vector{Tuple{String, VT}}
) where {T<:Integer, VT<:AbstractVector{T}}
    _create_noncommutative_variables(UnipotentAlgebra, prefix_subscripts)
end

"""
    create_noncommutative_variables(prefix_subscripts::Vector{Tuple{String, VT}})

Create a variable registry and monomials for generic non-commutative variables.

Non-commutative variables have no simplification rules - word order is
preserved exactly as given.

# Arguments
- `prefix_subscripts`: Vector of `(prefix, subscripts)` tuples for multi-prefix creation

# Returns
A tuple of:
- `VariableRegistry`: Registry with non-commutative variables for each subscript
- `monomial_groups`: Tuple of monomial vectors, one per prefix group

# Examples
```jldoctest
julia> reg, (x,) = create_noncommutative_variables([("x", 1:3)]);

julia> length(x)
3

julia> typeof(x[1])  # monomial for x₁
NormalMonomial{NonCommutativeAlgebra, UInt8}

julia> reg, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 3:4)]);

julia> length(x), length(y)
(2, 2)
```
"""
function create_noncommutative_variables(
    prefix_subscripts::Vector{Tuple{String, VT}}
) where {T<:Integer, VT<:AbstractVector{T}}
    _create_noncommutative_variables(NonCommutativeAlgebra, prefix_subscripts)
end

"""
    _create_noncommutative_variables(::Type{A}, prefix_subscripts::Vector{Tuple{String, VT}})

Internal function to create variables with multiple prefixes and site-encoded grouped indices.
Each tuple `(prefix, subscripts)` defines a group of variables with a specific prefix.
Variables from each group are assigned to a distinct physical site (1-indexed by group order).

# Arguments
- `A`: The algebra type for the monomials
- `prefix_subscripts`: Vector of tuples, each containing a prefix string and subscript collection

# Returns
A tuple of:
- `VariableRegistry`: Registry mapping symbols to encoded indices
- `monomial_groups`: Tuple of `Vector{NormalMonomial{A,T}}` for each prefix group

# Examples
```julia
# Create P₁, P₂, P₃ (site 1) and Q₄, Q₅ (site 2)
_create_noncommutative_variables(ProjectorAlgebra, [("P", 1:3), ("Q", 4:5)])
```
"""
function _create_noncommutative_variables(
    ::Type{A},
    prefix_subscripts::Vector{Tuple{String, VT}}
) where {A<:AlgebraType, T<:Integer, VT<:AbstractVector{T}}
    # Check for duplicate subscripts within each group
    for (prefix, subs) in prefix_subscripts
        _check_no_duplicate_subscripts(subs)
    end
    # Check for symbol collisions across groups
    all_syms = Symbol[]
    for (prefix, subs) in prefix_subscripts
        for s in subs
            sym = Symbol(prefix * _subscript_string(s))
            if sym in all_syms
                throw(ArgumentError("Duplicate variable symbol: $sym. Use distinct prefixes or subscripts."))
            end
            push!(all_syms, sym)
        end
    end
    n_operators = sum(x -> length(x[2]), prefix_subscripts; init=0)
    n_sites = length(prefix_subscripts)
    IndexT = select_uint_type(n_operators, n_sites)

    all_symbols = Vector{Symbol}(undef, n_operators)
    all_indices = Vector{IndexT}(undef, n_operators)

    global_idx = 1
    for (physical_site, (prefix, subscript_gp)) in enumerate(prefix_subscripts)
        for subscript in subscript_gp
            subscript_str = _subscript_string(subscript)
            encoded_idx = encode_index(IndexT, global_idx, physical_site)
            all_symbols[global_idx] = Symbol(prefix * subscript_str)
            all_indices[global_idx] = encoded_idx
            global_idx += 1
        end
    end

    idx_to_vars = Dict(zip(all_indices, all_symbols))
    variables_to_idx = Dict(zip(all_symbols, all_indices))

    reg = VariableRegistry{A, IndexT}(idx_to_vars, variables_to_idx)
    _register_display_registry(reg)

    # Build grouped monomial vectors
    monomial_groups = Vector{Vector{NormalMonomial{A,IndexT}}}()

    for (prefix, subscripts) in prefix_subscripts
        group = Vector{NormalMonomial{A,IndexT}}()
        for subscript in subscripts
            sym = Symbol(prefix * _subscript_string(subscript))
            idx = reg[sym]
            push!(group, NormalMonomial{A,IndexT}([idx]))
        end
        push!(monomial_groups, group)
    end

    return (reg, Tuple(monomial_groups))
end
