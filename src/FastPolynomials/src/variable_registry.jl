"""
    VariableRegistry

A registry that maps between variable symbols and their indices.
Supports both contiguous (Vector-like) and sparse (Dict-based) indices.

# Fields
- `idx_to_variables::Dict{T, Symbol}`: Maps index to symbol
- `variables_to_idx::Dict{Symbol, T}`: Maps symbol to index

# Examples
```jldoctest
julia> reg = create_pauli_variables(1:2)
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
struct VariableRegistry{T<:Integer}
    idx_to_variables::Dict{T, Symbol}
    variables_to_idx::Dict{Symbol,T}

    # Inner constructor to enforce invariants
    function VariableRegistry{T}(idx_to_vars::Dict{T, Symbol}, vars_to_idx::Dict{Symbol,T}) where {T<:Integer}
        # Verify consistency
        length(idx_to_vars) == length(vars_to_idx) ||
            error("Inconsistent registry: different lengths")

        for (idx, sym) in idx_to_vars
            haskey(vars_to_idx, sym) && vars_to_idx[sym] == idx ||
                error("Inconsistent registry: symbol $sym has wrong index")
        end

        new{T}(idx_to_vars, vars_to_idx)
    end
end

# Unicode subscript digits: ₀₁₂₃₄₅₆₇₈₉
const SUBSCRIPT_DIGITS = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']

"""
    _subscript_string(n::Int) -> String

Convert an integer to its Unicode subscript representation.

# Example
```julia
_subscript_string(42)  # "₄₂"
_subscript_string(0)   # "₀"
```
"""
function _subscript_string(n::Int)

    if n < 0
        error("Negative subscripts not supported: $n")
    end

    # Convert to string of digits
    digits_str = string(n)

    # Map each digit to subscript
    return String([SUBSCRIPT_DIGITS[parse(Int, d)+1] for d in digits_str])
end

"""
    Base.length(reg::VariableRegistry) -> Int

Return the number of variables in the registry.
"""
Base.length(reg::VariableRegistry) = length(reg.idx_to_variables)

"""
    Base.getindex(reg::VariableRegistry{T}, idx::T) -> Symbol

Get the variable symbol at the given index.
"""
Base.getindex(reg::VariableRegistry{T}, idx::T) where {T} = reg.idx_to_variables[idx]

"""
    Base.getindex(reg::VariableRegistry{T}, idx::Integer) -> Symbol

Get the variable symbol at the given index (converts to type T).
"""
Base.getindex(reg::VariableRegistry{T}, idx::Integer) where {T} = reg.idx_to_variables[T(idx)]

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
julia> reg = create_projector_variables(1:3)
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
    indices(reg::VariableRegistry{T}) -> Vector{T}

Return all indices in sorted order.

# Examples
```jldoctest
julia> reg = create_projector_variables(1:3)
VariableRegistry with 3 variables: P₁, P₂, P₃

julia> indices(reg)
3-element Vector{UInt8}:
 0x01
 0x02
 0x03
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

# AlgebraType is defined in algebra_types.jl as singleton types
# (Commutative, PauliAlgebra, FermionicAlgebra, BosonicAlgebra, ProjectorAlgebra, UnipotentAlgebra)

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

julia> σx[1]  # Monomial for σx₁
Monomial{PauliAlgebra, UInt8}(...)

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

    reg = VariableRegistry{T}(idx_to_vars, variables_to_idx)

    # Build monomial vectors grouped by Pauli type (x, y, z)
    σx = Vector{Monomial{PauliAlgebra, T}}(undef, n_sites)
    σy = Vector{Monomial{PauliAlgebra, T}}(undef, n_sites)
    σz = Vector{Monomial{PauliAlgebra, T}}(undef, n_sites)

    for (site_idx, subscript) in enumerate(subscripts)
        subscript_str = _subscript_string(subscript)
        σx[site_idx] = Monomial{PauliAlgebra}([reg[Symbol("σx" * subscript_str)]])
        σy[site_idx] = Monomial{PauliAlgebra}([reg[Symbol("σy" * subscript_str)]])
        σz[site_idx] = Monomial{PauliAlgebra}([reg[Symbol("σz" * subscript_str)]])
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

Create a variable registry and monomials for fermionic creation (a⁺) and annihilation (a) operators.

Fermionic operators satisfy anticommutation relations:
- {aᵢ, aⱼ⁺} = δᵢⱼ (creation-annihilation)
- {aᵢ, aⱼ} = 0 (annihilation-annihilation)
- {aᵢ⁺, aⱼ⁺} = 0 (creation-creation)

# Arguments
- `subscripts`: Subscript values for different modes/sites

# Returns
A tuple of:
- `VariableRegistry`: Registry with both a and a⁺ operators for each subscript
- `(a, a⁺)`: Tuple of monomial vectors (annihilation, creation)

# Examples
```jldoctest
julia> reg, (a, a⁺) = create_fermionic_variables(1:2);

julia> length(a)  # One annihilation per mode
2

julia> a[1]  # Monomial for a₁
Monomial{FermionicAlgebra, Int8}(...)

julia> :a₁ in reg  # annihilation operator
true

julia> :a⁺₁ in reg  # creation operator
true
```

Multiple modes:
```jldoctest
julia> reg, (a, a⁺) = create_fermionic_variables(1:3);

julia> length(a)
3

julia> length(a⁺)
3
```
"""
create_fermionic_variables(subscripts) = _create_physical_variables(FermionicAlgebra, subscripts, "a")

function _create_physical_variables(::Type{A}, subscripts, prefix::String) where {A<:AlgebraType}
    n_modes = length(subscripts)

    # Select signed type that can hold ±n_modes
    T = _select_signed_index_type(n_modes)

    idx_to_vars = Dict{T, Symbol}()
    variables_to_idx = Dict{Symbol, T}()

    for (i, subscript) in enumerate(subscripts)
        subscript_str = _subscript_string(subscript)

        # Annihilation operator: positive index
        ann_sym = Symbol(prefix * subscript_str)
        idx_to_vars[T(i)] = ann_sym
        variables_to_idx[ann_sym] = T(i)

        # Creation operator: negative index (same absolute value)
        cre_sym = Symbol(prefix * "⁺" * subscript_str)
        idx_to_vars[T(-i)] = cre_sym
        variables_to_idx[cre_sym] = T(-i)
    end

    reg = VariableRegistry{T}(idx_to_vars, variables_to_idx)

    # Build monomial vectors grouped by operator type (annihilation, creation)
    annihilation = Vector{Monomial{A, T}}(undef, n_modes)
    creation = Vector{Monomial{A, T}}(undef, n_modes)

    for (i, subscript) in enumerate(subscripts)
        subscript_str = _subscript_string(subscript)
        annihilation[i] = Monomial{A}([reg[Symbol(prefix * subscript_str)]])
        creation[i] = Monomial{A}([reg[Symbol(prefix * "⁺" * subscript_str)]])
    end

    return (reg, (annihilation, creation))
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

Create a variable registry and monomials for bosonic creation (c⁺) and annihilation (c) operators.

Bosonic operators satisfy commutation relations:
- [cᵢ, cⱼ⁺] = δᵢⱼ (creation-annihilation)
- [cᵢ, cⱼ] = 0 (annihilation-annihilation)
- [cᵢ⁺, cⱼ⁺] = 0 (creation-creation)

# Arguments
- `subscripts`: Subscript values for different modes/sites

# Returns
A tuple of:
- `VariableRegistry`: Registry with both c and c⁺ operators for each subscript
- `(c, c⁺)`: Tuple of monomial vectors (annihilation, creation)

# Examples
```jldoctest
julia> reg, (c, c⁺) = create_bosonic_variables(1:2);

julia> length(c)  # One annihilation per mode
2

julia> c[1]  # Monomial for c₁
Monomial{BosonicAlgebra, Int8}(...)

julia> :c₁ in reg  # annihilation operator
true

julia> :c⁺₁ in reg  # creation operator
true
```

Multiple modes:
```jldoctest
julia> reg, (c, c⁺) = create_bosonic_variables(1:3);

julia> length(c)
3

julia> length(c⁺)
3
```
"""
create_bosonic_variables(subscripts) = _create_physical_variables(BosonicAlgebra, subscripts, "c")

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

julia> P[1]  # Monomial for P₁
Monomial{ProjectorAlgebra, UInt8}(...)

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

"""
    _select_unsigned_type(n_operators::Int, n_sites::Int) -> Type{<:Unsigned}

Select the smallest unsigned integer type that can encode the given number
of operators and sites using the site-based encoding scheme.
"""
function _select_unsigned_type(n_operators::Int, n_sites::Int)
    for T in (UInt8, UInt16, UInt32, UInt64)
        if n_sites <= max_sites(T) && n_operators <= max_operators(T)
            return T
        end
    end
    error("Cannot fit $n_operators operators × $n_sites sites in any UInt type")
end

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

julia> U[1]  # Monomial for U₁
Monomial{UnipotentAlgebra, UInt8}(...)

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

julia> x[1]  # Monomial for x₁
Monomial{NonCommutativeAlgebra, UInt8}(...)

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
    _select_contiguous_type(n::Int) -> Type{<:Unsigned}

Select smallest unsigned type for contiguous 1-indexed variables.
"""
function _select_contiguous_type(n::Int)
    n <= typemax(UInt8) && return UInt8
    n <= typemax(UInt16) && return UInt16
    n <= typemax(UInt32) && return UInt32
    return UInt64
end

@inline _encode_index(::Type{T}, global_idx::Int, physical_site::Int) where {T<:Integer} = T(global_idx << (sizeof(T) * 2) | T(physical_site))

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
- `monomial_groups`: Tuple of `Vector{Monomial{A,T}}` for each prefix group

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
    n = sum(x -> length(x[2]), prefix_subscripts)
    IndexT = _select_contiguous_type(n)

    all_symbols = Vector{Symbol}(undef, n)
    all_indices = Vector{IndexT}(undef, n)

    global_idx = 1
    for (physical_site, (prefix, subscript_gp)) in enumerate(prefix_subscripts)
        for subscript in subscript_gp
            subscript_str = _subscript_string(subscript)
            encoded_idx = _encode_index(IndexT, global_idx, physical_site)
            all_symbols[global_idx] = Symbol(prefix * subscript_str)
            all_indices[global_idx] = encoded_idx
            global_idx += 1
        end
    end

    idx_to_vars = Dict(zip(all_indices, all_symbols))
    variables_to_idx = Dict(zip(all_symbols, all_indices))

    reg = VariableRegistry{IndexT}(idx_to_vars, variables_to_idx)

    # Build grouped monomial vectors
    monomial_groups = Vector{Vector{Monomial{A, IndexT}}}()

    for (prefix, subscripts) in prefix_subscripts
        group = Vector{Monomial{A, IndexT}}()
        for subscript in subscripts
            sym = Symbol(prefix * _subscript_string(subscript))
            idx = reg[sym]
            push!(group, Monomial{A}([idx]))
        end
        push!(monomial_groups, group)
    end

    return (reg, Tuple(monomial_groups))
end

