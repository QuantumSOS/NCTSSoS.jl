# H₄ periodic V2RDM active-space asset regression
#
# The current package still lacks the k-blocked / spin-adapted formulation
# needed for a faithful Nk=2 periodic V2RDM solve.  This regression therefore
# locks down the asset-backed pieces we can verify exactly today:
#   1. the vendored Nk=2 integrals and metadata,
#   2. the HF reconstruction and additive constant shift,
#   3. the dense-graph counts showing why a naive 32-mode order-2 lift is not
#      a good CI target on the current API.

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
end
