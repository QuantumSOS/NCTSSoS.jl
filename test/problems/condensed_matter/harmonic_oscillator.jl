# Quantum Harmonic Oscillator (Bosonic Sanity Check)
#
# H = ω (b† b + 1/2)
#
# Exact E₀ = ω/2 = 0.5 for ω = 1.
# This should be exactly reproduced by even a low-order SDP relaxation,
# since the objective is linear in the moment ⟨b† b⟩ and ⟨b† b⟩ ≥ 0.

using Test, NCTSSoS, JuMP

const HARMONIC_EXPECTATIONS_PATH = "expectations/harmonic_oscillator.toml"

@testset "Quantum harmonic oscillator" begin

    @testset "Single mode, ω=1, order 2" begin
        oracle = expectations_oracle(HARMONIC_EXPECTATIONS_PATH, "single_mode_order2")

        registry, (b, b_dag) = create_bosonic_variables(1:1)
        ham = b_dag[1] * b[1] + 0.5  # ω = 1

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=2))

        @test result.objective ≈ 0.5 atol = 1e-6
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
