"""Experimental Pauli compiler and two-phase solver, isolated from the NCTSSoS API."""
module PauliWorkflowRedesign

import NCTSSoS as NC
using NCTSSoS: JuMP
import LinearAlgebra
using LinearAlgebra: Hermitian, eigen
using SparseArrays: sparse

export PauliSymmetrySpec, PauliSparsitySpec, PauliStrengtheningSpec,
    BoundProvenance, CompiledSDPPlan, resolve_workflow_preset,
    compile_relaxation, solve_ground_state, solve_pauli_ground_state, pauli_nctssos,
    bound_observable_interval, certify_rational_sos, certify_energy_window, certify_observable_interval

Base.@kwdef struct PauliSymmetrySpec
    time_reversal_gauge::Bool = true
    klein4_sign::Bool = true
    anti_sign_realify::Bool = true
    xyz_permutation::Bool = true
    rotated_lie_detection::Bool = true
    cartan_u1_weyl::Bool = true
    su2_mz0_tree::Bool = true
    spatial_mode::Symbol = :auto
    drop_irreps::Union{Nothing,Vector{Int}} = nothing
end

Base.@kwdef struct PauliSparsitySpec
    support_mode::Symbol = :auto
    power_graph_radius::Int = 2
    local_product_depth::Int = 2
    min_cluster_size::Int = 4
    merge_threshold::Int = 8
    two_body_range_r::Int = 4
    cs_algo::Symbol = :NoElimination
    ts_algo::Symbol = :Maximal
    ts_step::Int = 0
    rewrite_complete_local_rdm::Bool = true
end

Base.@kwdef struct PauliStrengtheningSpec
    commutator_linear_order::Int = 0
    kms_psd_order::Int = 0
    kms_frustration_guard::Bool = true
    kms_null_tol::Float64 = 1e-7
    rdm_positivity_k::Int = 0
    lightcone_depth::Int = 0
    rg_coarse_grain_depth::Int = 0
    energy_upper_bound::Union{Nothing,Float64} = nothing
    target_mz_sector::Union{Nothing,Int} = nothing
    target_spin_j::Union{Nothing,Real} = nothing
    rational_sos_certify::Bool = false
end

"""A record of the relaxation and the *reported* bound, not a solver certificate."""
struct BoundProvenance
    bound_status::Symbol
    rigorous_lower_bound::Bool
    exact_wrt_dense_npa::Bool
    witnesses::Tuple
    active_relaxations::Tuple{Vararg{Symbol}}
    active_strengthenings::Tuple{Vararg{Symbol}}
    sector_targets::Tuple{Vararg{Symbol}}
    trivial_pauli_floor::Float64
end

"Resolve a named recipe into independent, immutable option sets."
function resolve_workflow_preset(preset::Symbol)
    if preset in (:ExactNPA, :ExactHierarchy)
        return (PauliSymmetrySpec(), PauliSparsitySpec(support_mode=:full), PauliStrengtheningSpec())
    elseif preset === :ScalableLattice
        return (PauliSymmetrySpec(spatial_mode=:circulant_dft),
            PauliSparsitySpec(support_mode=:contiguous_lattice),
            PauliStrengtheningSpec(commutator_linear_order=2, rdm_positivity_k=4))
    elseif preset in (:GeneralGraphFast, :FastScalable)
        return (PauliSymmetrySpec(spatial_mode=:clique_orbits),
            PauliSparsitySpec(support_mode=:correlative_cliques, cs_algo=:MF, ts_algo=:MMD, ts_step=1),
            PauliStrengtheningSpec(commutator_linear_order=2))
    elseif preset === :TightGroundState
        return (PauliSymmetrySpec(), PauliSparsitySpec(cs_algo=:MF),
            PauliStrengtheningSpec(commutator_linear_order=3, kms_psd_order=2,
                rdm_positivity_k=6, lightcone_depth=2))
    end
    throw(ArgumentError("Unsupported workflow preset: $preset"))
end

include("symmetry.jl")
include("sparsity.jl")
include("strengthening.jl")
include("wang.jl")
include("local_density.jl")
include("certificates.jl")
include("compiler.jl")
include("solve.jl")

"Compile the relaxation, returning the plan when no optimizer was supplied."
function solve_pauli_ground_state(H; optimizer=nothing, keywords...)
    plan = compile_relaxation(H; keywords...)
    return isnothing(optimizer) ? plan : solve_ground_state(plan, H; optimizer)
end

const pauli_nctssos = solve_pauli_ground_state

# Direct-include clients do not need a package environment for precompilation.
if Base.PkgId(@__MODULE__).uuid !== nothing
    include("src/precompile.jl")
end

end
