function _fermionic_irrep_table(sector::FermionicSectorSpec)
    table = sector.mode_layout.irrep_table
    isnothing(table) && throw(ArgumentError(
        "`split_irrep=true` needs an explicit `irrep_table` in the supplied `FermionicModeLayout`."
    ))
    return table
end

@inline _fermionic_irrep_identity(sector::FermionicSectorSpec) = _fermionic_irrep_table(sector).identity
@inline _fermionic_irrep_multiply(sector::FermionicSectorSpec, left, right) =
    _fermionic_irrep_table(sector).multiply(left, right)
@inline _fermionic_irrep_dual(sector::FermionicSectorSpec, irrep) =
    _fermionic_irrep_table(sector).dual(irrep)

function _fermionic_mode_orbital(mode_layout::FermionicModeLayout, mode::Int; context::AbstractString="Fermionic orbital metadata")
    haskey(mode_layout.orbital_of, mode) || throw(ArgumentError(
        "$context needs `orbital_of[$mode]` in the supplied `FermionicModeLayout`."
    ))
    return mode_layout.orbital_of[mode]
end

function _fermionic_mode_orbital(sector::FermionicSectorSpec, mode::Int)
    return _fermionic_mode_orbital(
        sector.mode_layout,
        mode;
        context="Fermionic irrep-sector splitting",
    )
end

function _fermionic_orbital_irrep(mode_layout::FermionicModeLayout, orbital::Int; context::AbstractString="Fermionic orbital metadata")
    haskey(mode_layout.irrep_of, orbital) || throw(ArgumentError(
        "$context needs `irrep_of[$orbital]` in the supplied `FermionicModeLayout`."
    ))
    return mode_layout.irrep_of[orbital]
end

function _fermionic_mode_irrep(sector::FermionicSectorSpec, mode::Int)
    orbital = _fermionic_mode_orbital(sector, mode)
    return _fermionic_orbital_irrep(
        sector.mode_layout,
        orbital;
        context="Fermionic irrep-sector splitting",
    )
end

function _validate_abelian_irrep_table(table::AbelianIrrepTable, irreps::AbstractVector)
    dual_irreps = [table.dual(irrep) for irrep in irreps]
    values = unique!(Any[table.identity; irreps; dual_irreps])

    for irrep in values
        table.multiply(table.identity, irrep) == irrep || throw(ArgumentError(
            "`AbelianIrrepTable` identity must act neutrally on the left; failed for irrep $(repr(irrep))."
        ))
        table.multiply(irrep, table.identity) == irrep || throw(ArgumentError(
            "`AbelianIrrepTable` identity must act neutrally on the right; failed for irrep $(repr(irrep))."
        ))

        dual_irrep = table.dual(irrep)
        table.dual(dual_irrep) == irrep || throw(ArgumentError(
            "`AbelianIrrepTable` dual must be involutive; failed for irrep $(repr(irrep))."
        ))
        table.multiply(irrep, dual_irrep) == table.identity || throw(ArgumentError(
            "`AbelianIrrepTable` dual must invert on the right; failed for irrep $(repr(irrep))."
        ))
        table.multiply(dual_irrep, irrep) == table.identity || throw(ArgumentError(
            "`AbelianIrrepTable` dual must invert on the left; failed for irrep $(repr(irrep))."
        ))
    end

    for left in values, right in values
        table.multiply(left, right) == table.multiply(right, left) || throw(ArgumentError(
            "`AbelianIrrepTable` must be Abelian on the supplied active irreps; failed for $(repr(left)) and $(repr(right))."
        ))
    end

    return nothing
end

function _validate_active_fermionic_sector_metadata(
    sector::FermionicSectorSpec,
    active_modes::AbstractVector{<:Integer},
)
    if sector.split_spin
        for mode in active_modes
            haskey(sector.mode_layout.spin2_of, mode) || throw(ArgumentError(
                "Fermionic spin-sector splitting needs `spin2_of[$mode]` for every active mode."
            ))
        end
    end

    sector.split_irrep || return nothing

    active_orbitals = Int[]
    sizehint!(active_orbitals, length(active_modes))
    for mode in active_modes
        push!(active_orbitals, _fermionic_mode_orbital(sector, Int(mode)))
    end
    unique!(sort!(active_orbitals))

    for orbital in active_orbitals
        haskey(sector.mode_layout.irrep_of, orbital) || throw(ArgumentError(
            "`split_irrep=true` needs `irrep_of[$orbital]` for every active orbital."
        ))
    end

    active_irreps = Any[_fermionic_orbital_irrep(sector.mode_layout, orbital) for orbital in active_orbitals]
    _validate_abelian_irrep_table(_fermionic_irrep_table(sector), active_irreps)
    return nothing
end
