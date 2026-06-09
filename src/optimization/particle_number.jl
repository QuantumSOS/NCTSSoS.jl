"""
    particle_number_constraint(registry::VariableRegistry{A,T}, N::Integer)
    particle_number_constraint(registry::VariableRegistry{A,T}, group => N, more_groups => more_Ns...)

Construct one-sided moment-equality constraints fixing particle number for
fermionic or bosonic `polyopt` problems.

The total form returns a one-element vector containing
`N̂ - N⋅I`, where `N̂ = Σᵢ aᵢ†aᵢ` is built from every annihilation mode in
`registry`.  The grouped form returns one constraint per group, using the
annihilation-operator vectors returned by [`create_fermionic_variables`](@ref)
or [`create_bosonic_variables`](@ref).  Groups may be spin species, mode
subsets, or single-mode vectors such as `[a[1]]`.

Use the returned vector directly in [`polyopt`](@ref):

```jldoctest
julia> reg, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
           ("c_up", 1:2), ("c_dn", 1:2)
       ]);

julia> cons = particle_number_constraint(reg, c_up => 1, c_dn => 1);

julia> length(cons)
2

julia> cons[1] == sum(c_up_dag[i] * c_up[i] for i in eachindex(c_up)) - one(cons[1])
true
```

!!! note "State constraints, not algebra identities"
    These polynomials are intended for `moment_eq_constraints`, which encode
    state-sector conditions `g|ψ⟩ = 0`.  Do not pass them as `eq_constraints`.
    State/trace polynomial optimization does not yet support
    `moment_eq_constraints`.

!!! note "Bosonic truncation"
    Fixing total particle number does not impose per-mode occupation caps for
    bosons.  Add suitable `ineq_constraints` when your relaxation also needs
    finite local occupation bounds.
"""
function particle_number_constraint(
    registry::VariableRegistry{A,T},
    N::Integer,
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    modes = _particle_number_registry_modes(registry)
    isempty(modes) && throw(ArgumentError("cannot build a particle-number constraint from a registry with no annihilation modes"))
    return [_particle_number_constraint_from_modes(A, modes, N)]
end

function particle_number_constraint(
    registry::VariableRegistry{A,T},
    groups::Pair...,
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    isempty(groups) && throw(ArgumentError("pass a particle number N, or at least one `group => N` pair"))

    constraints = Polynomial{A,T,Float64}[]
    sizehint!(constraints, length(groups))

    for group_target in groups
        group, target = group_target.first, group_target.second
        target isa Integer || throw(ArgumentError("particle-number target must be an integer; got $(typeof(target))"))
        modes = _particle_number_group_modes(registry, group)
        push!(constraints, _particle_number_constraint_from_modes(A, modes, target))
    end

    return constraints
end

function particle_number_constraint(
    ::VariableRegistry{A,T},
    _args...,
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    throw(ArgumentError("expected an integer particle number or one or more `group => N` pairs"))
end

function particle_number_constraint(::VariableRegistry{A}, _args...) where {A<:AlgebraType}
    throw(ArgumentError(
        "particle_number_constraint is only supported for FermionicAlgebra and BosonicAlgebra registries; " *
        "got $(nameof(A))."
    ))
end

function _particle_number_registry_modes(registry::VariableRegistry{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    modes = sort!(T[idx for idx in keys(registry.idx_to_variables) if idx > zero(T)])
    for mode in modes
        _validate_particle_number_mode(registry, mode)
    end
    return modes
end

function _particle_number_group_modes(
    registry::VariableRegistry{A,T},
    group,
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    group isa AbstractVector || throw(ArgumentError(
        "particle-number groups must be vectors of annihilation monomials, e.g. `c_up => 2` or `[c[1]] => 1`; " *
        "got $(typeof(group))."
    ))
    isempty(group) && throw(ArgumentError("particle-number groups must not be empty"))

    modes = T[]
    sizehint!(modes, length(group))
    seen = Set{T}()

    for mono in group
        mode = _particle_number_mode(registry, mono)
        mode in seen && throw(ArgumentError("particle-number group repeats mode $mode"))
        push!(seen, mode)
        push!(modes, mode)
    end

    sort!(modes)
    return modes
end

function _particle_number_mode(
    registry::VariableRegistry{A,T},
    mono::NormalMonomial{A,T},
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    length(mono.word) == 1 || throw(ArgumentError(
        "particle-number groups must contain annihilation generators, not products; got word $(mono.word)."
    ))

    mode = only(mono.word)
    mode > zero(T) || throw(ArgumentError(
        "particle-number groups must contain annihilation operators (positive PBW indices); got creation index $mode."
    ))
    _validate_particle_number_mode(registry, mode)
    return mode
end

function _particle_number_mode(
    ::VariableRegistry{A,T},
    mono,
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    throw(ArgumentError(
        "particle-number groups must contain `NormalMonomial{$(nameof(A)),$T}` annihilation operators; got $(typeof(mono))."
    ))
end

function _validate_particle_number_mode(
    registry::VariableRegistry{A,T},
    mode::T,
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    mode > zero(T) || throw(ArgumentError("particle-number mode indices must be positive; got $mode"))
    haskey(registry.idx_to_variables, mode) || throw(ArgumentError(
        "annihilation mode $mode is not present in the registry"
    ))
    haskey(registry.idx_to_variables, -mode) || throw(ArgumentError(
        "creation mode $(-mode) paired with annihilation mode $mode is not present in the registry"
    ))
    return nothing
end

function _particle_number_constraint_from_modes(
    ::Type{A},
    modes::Vector{T},
    target::Integer,
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Signed}
    _validate_particle_number_target(A, target, length(modes))

    terms = Tuple{Float64,NormalMonomial{A,T}}[]
    sizehint!(terms, length(modes) + 1)

    for mode in modes
        push!(terms, (1.0, NormalMonomial{A,T}(T[-mode, mode])))
    end

    if !iszero(target)
        push!(terms, (-_particle_number_target_coefficient(target), one(NormalMonomial{A,T})))
    end

    return Polynomial(terms)
end

function _particle_number_target_coefficient(target::Integer)
    max_exact = BigInt(2)^precision(Float64)
    BigInt(target) <= max_exact || throw(DomainError(
        target,
        "particle number is too large to represent exactly as a Float64 coefficient",
    ))
    return Float64(target)
end

function _validate_particle_number_target(::Type{A}, target::Integer, _n_modes::Integer) where {A<:BosonicAlgebra}
    target >= 0 || throw(DomainError(target, "particle number must be non-negative"))
    _particle_number_target_coefficient(target)
    return nothing
end

function _validate_particle_number_target(::Type{A}, target::Integer, n_modes::Integer) where {A<:FermionicAlgebra}
    target >= 0 || throw(DomainError(target, "particle number must be non-negative"))
    target <= n_modes || throw(DomainError(target, "fermionic particle number cannot exceed the number of modes in the constraint group ($n_modes)"))
    _particle_number_target_coefficient(target)
    return nothing
end
