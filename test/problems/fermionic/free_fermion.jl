# Free Fermion Hopping Benchmarks
#
# These are the simplest possible fermionic benchmarks: non-interacting tight-binding
# chains that should be solved exactly by the order-1 fermionic SDP relaxation.
#
# Two-site hopping: H = -(a†₁ a₂ + a†₂ a₁), E₀ = -1.0 (bonding orbital)
# N-site chain OBC: H = -Σ (a†_i a_{i+1} + h.c.), E₀ from diagonalizing the hopping matrix

using Test, NCTSSoS, JuMP, LinearAlgebra

const FREE_FERMION_EXPECTATIONS_PATH = "expectations/free_fermion.toml"

# ─── Exact diagonalization for tight-binding chains ───

const _FF_PAULI_Z = ComplexF64[1 0; 0 -1]
const _FF_PAULI_I = Matrix{ComplexF64}(I, 2, 2)
const _FF_JW_SM = ComplexF64[0 0; 1 0]
const _FF_JW_SP = ComplexF64[0 1; 0 0]

function _ff_jw_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = Vector{Matrix{ComplexF64}}(undef, nmodes)
    for site in 1:nmodes
        if site < mode
            mats[site] = _FF_PAULI_Z
        elseif site == mode
            mats[site] = kind === :annihilate ? _FF_JW_SM : _FF_JW_SP
        else
            mats[site] = _FF_PAULI_I
        end
    end
    return reduce(kron, mats)
end

function _ff_chain_exact_energy(N::Int; t::Real=1.0)
    a = [_ff_jw_op(i, :annihilate, N) for i in 1:N]
    ad = [_ff_jw_op(i, :create, N) for i in 1:N]
    dim = 2^N
    H = zeros(ComplexF64, dim, dim)
    for i in 1:N-1
        H .-= t .* (ad[i] * a[i+1] + ad[i+1] * a[i])
    end
    return eigmin(Hermitian(H))
end

# ─── Tests ───

@testset "Free fermion hopping benchmarks" begin

    @testset "Two-site hopping (order 1)" begin
        oracle = expectations_oracle(FREE_FERMION_EXPECTATIONS_PATH, "two_site_hopping_order1")

        registry, (a, a_dag) = create_fermionic_variables(1:2)
        ham = -(a_dag[1] * a[2] + a_dag[2] * a[1])

        # Exact check
        exact_e0 = _ff_chain_exact_energy(2)
        @test exact_e0 ≈ -1.0 atol = 1e-12

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=1))
        @test result.objective ≈ -1.0 atol = 1e-6
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Two-site hopping (order 2)" begin
        oracle = expectations_oracle(FREE_FERMION_EXPECTATIONS_PATH, "two_site_hopping_order2")

        registry, (a, a_dag) = create_fermionic_variables(1:2)
        ham = -(a_dag[1] * a[2] + a_dag[2] * a[1])

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=2))
        @test result.objective ≈ -1.0 atol = 1e-6
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Pure hopping chain N=4 OBC (order 1)" begin
        oracle = expectations_oracle(FREE_FERMION_EXPECTATIONS_PATH, "pure_hopping_chain_n4_order1")

        N = 4
        registry, (a, a_dag) = create_fermionic_variables(1:N)
        ham = -sum(a_dag[i] * a[i + 1] + a_dag[i + 1] * a[i] for i in 1:N-1)

        exact_e0 = _ff_chain_exact_energy(N)
        @test exact_e0 ≈ -sqrt(5) atol = 1e-10

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=1))
        @test result.objective ≈ exact_e0 atol = 1e-6
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
