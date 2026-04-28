# H₄ periodic V2RDM active-space asset regression
#
# The core package still lacks a first-class periodic V2RDM API, but the repo
# now carries a reviewed native bridge benchmark helper under `demos/`. This
# regression therefore locks down:
#   1. the vendored Nk=2 integrals and metadata,
#   2. the HF reconstruction and additive constant shift,
#   3. the dense-graph counts showing why a naive 32-mode order-2 lift is not
#      a good CI target on the current API,
#   4. the reviewed K-only and spin-resolved native bridge block counts.

using Test, NCTSSoS, LinearAlgebra

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: H4_EXPECTATIONS_PATH,
                         H4_INTEGRALS_PATH,
                         H4_REFERENCE_PATH,
                         build_k2,
                         compound_spatial_edge_count,
                         complete_graph_edge_count,
                         fermionic_order2_basis_size,
                         fermionic_order2_nuniq,
                         hf_energy,
                         load_nk2_asset

include(joinpath(pkgdir(NCTSSoS), "demos", "H4PeriodicNativeV2RDMHelpers.jl"))
using .H4PeriodicNativeV2RDMHelpers: load_native_problem_data,
                                     native_size_summary

@testset "H4 periodic V2RDM active-space asset regression" begin
    asset = load_nk2_asset()
    (; reference, preflight, blocker, figure, nk, n_active_orb, n_active_elec,
       total_active_electrons, total_spatial_orbitals, total_spin_orbital_modes,
       h1e, eri) = asset

    k2 = build_k2(h1e, eri; nk, n_active_orb, n_total_electrons = total_active_electrons)

    @testset "Vendored Nk=2 asset stays internally consistent" begin
        @test H4_EXPECTATIONS_PATH == joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "h4_periodic_v2rdm.toml")
        @test H4_REFERENCE_PATH == joinpath(pkgdir(NCTSSoS), "test", "data", "assets", "h4_chain_reference.toml")
        @test H4_INTEGRALS_PATH == joinpath(pkgdir(NCTSSoS), "test", "data", "assets", "h4_chain_nk2_integrals.txt")

        @test reference["system"]["active_space_electrons"] == n_active_elec
        @test reference["system"]["active_space_orbitals"] == n_active_orb
        @test reference["pyscf_nk2"]["n_eri_blocks"] == preflight["eri_blocks"]
        @test reference["pyscf_nk2"]["HF_energy"] ≈ preflight["hf_total_energy"] atol = 1e-12

        @test figure["V2RDM_4_8"] == preflight["figure_v2rdm_nk2"]
        @test figure["HF"] == preflight["figure_hf_nk2"]

        @test sort(collect(keys(h1e))) == [0, 1]
        @test sort(collect(keys(eri))) == [
            (0, 0, 0, 0),
            (0, 0, 1, 1),
            (0, 1, 0, 1),
            (0, 1, 1, 0),
            (1, 0, 0, 1),
            (1, 0, 1, 0),
            (1, 1, 0, 0),
            (1, 1, 1, 1),
        ]

        for block in values(h1e)
            @test block ≈ block' atol = 1e-8
        end

        @test real.(diag(h1e[0])) ≈ preflight["h1e_diag_k0"] atol = 1e-10
        @test real.(diag(h1e[1])) ≈ preflight["h1e_diag_k1"] atol = 1e-10
        @test maximum(maximum(abs, value) for value in values(eri)) ≈ preflight["max_abs_eri"] atol = 1e-12
        @test maximum(abs, k2) ≈ preflight["max_abs_k2"] atol = 1e-12

        hf_electronic = hf_energy(h1e, eri; nk, n_active_elec)
        hf_total = hf_electronic + preflight["hf_constant_shift"]

        @test hf_electronic ≈ preflight["hf_active_space_electronic_energy"] atol = 1e-10
        @test hf_total ≈ preflight["hf_total_energy"] atol = 1e-10

        # Figure values are approximate digitizations, not exact regression targets.
        @test abs(hf_total - preflight["figure_hf_nk2"]) > preflight["figure_uncertainty_ha"]
    end

    @testset "Naive order-2 spin-orbital lift remains dense" begin
        spatial_edges = compound_spatial_edge_count(h1e, eri; n_active_orb)
        @test spatial_edges == blocker["spatial_graph_edges"]
        @test spatial_edges == complete_graph_edge_count(total_spatial_orbitals)
        @test blocker["spatial_graph_complete"] == true

        basis_size = Int(fermionic_order2_basis_size(total_spin_orbital_modes))
        nuniq = Int(fermionic_order2_nuniq(total_spin_orbital_modes))

        @test basis_size == blocker["order2_basis_size"]
        @test nuniq == blocker["order2_nuniq"]
    end

    @testset "Native bridge block counts stay reviewed" begin
        native = load_native_problem_data()
        k_only = native_size_summary(native; refinement = :k_only)
        spin_resolved = native_size_summary(native; refinement = :spin_resolved)

        @test length.(native.D_term_sparsity.block_bases) == [240, 256]
        @test length.(native.D_term_sparsity_spin.block_bases) == [56, 128, 56, 64, 128, 64]

        @test k_only.d_only_block_sizes == [240, 256]
        @test k_only.g_block_sizes == [512, 512]
        @test k_only.jump_variables == 123_136
        @test k_only.solver_psd_total_rows == 1_543_136
        @test k_only.jump_constraints_total == 11

        @test spin_resolved.d_only_block_sizes == [56, 128, 56, 64, 128, 64]
        @test spin_resolved.g_block_sizes == [128, 256, 128, 128, 256, 128]
        @test spin_resolved.jump_variables == 47_232
        @test spin_resolved.solver_psd_total_rows == 584_160
        @test spin_resolved.jump_constraints_total == 23
    end
end
