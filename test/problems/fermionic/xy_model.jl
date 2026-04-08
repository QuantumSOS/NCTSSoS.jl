# XY Model (Jordan–Wigner → free fermions)
#
# The spin-1/2 XY Hamiltonian maps to free fermions via Jordan–Wigner.
# Spin form:  H = (j_c/4) Σ_i (σˣ_i σˣ_{i+1} + σʸ_i σʸ_{i+1})
# Fermionic:  H = (j_c/2) Σ_i (a†_i a_{i+1} + a†_{i+1} a_i)
#
# Key subtlety: spin PBC maps to fermionic APBC (even parity) or PBC (odd parity)
# via the JW string wrapping around the ring. The ground state lives in whichever
# sector gives lower energy.
#
# Single-particle dispersion: ε_k = -j_c cos(k)
#
# Reference: textbook result; Baumgratz & Plenio (2012) Table 1.

using Test, NCTSSoS, JuMP, LinearAlgebra

const XY_EXPECTATIONS_PATH = "expectations/xy_model.toml"

# ─── Exact diagonalization oracle ───

const _XY_PAULI_X = ComplexF64[0 1; 1 0]
const _XY_PAULI_Y = ComplexF64[0 -im; im 0]
const _XY_PAULI_I = Matrix{ComplexF64}(I, 2, 2)

function _xy_site_op(site::Int, op::AbstractMatrix, nsites::Int)
    mats = [_XY_PAULI_I for _ in 1:nsites]
    mats[site] = Matrix{ComplexF64}(op)
    return reduce(kron, mats)
end

function _xy_exact_energy(N::Int; jc::Real=1.0, periodic::Bool=true)
    dim = 2^N
    H = zeros(ComplexF64, dim, dim)
    bound = periodic ? N : N - 1
    for i in 1:bound
        j = mod1(i + 1, N)
        H .+= (jc / 4) .* (_xy_site_op(i, _XY_PAULI_X, N) * _xy_site_op(j, _XY_PAULI_X, N) +
                            _xy_site_op(i, _XY_PAULI_Y, N) * _xy_site_op(j, _XY_PAULI_Y, N))
    end
    return eigmin(Hermitian(H))
end

# ─── Tests ───

@testset "XY model (free fermion benchmarks)" begin
    N = 4
    jc = 1.0

    @testset "Pauli form, N=$N PBC, order 2" begin
        oracle = expectations_oracle(XY_EXPECTATIONS_PATH, "pbc_pauli_n4_order2")
        exact_e0 = _xy_exact_energy(N; jc, periodic=true)
        @test exact_e0 ≈ -sqrt(2) atol = 1e-12

        registry, (x, y, _) = create_pauli_variables(1:N)
        ham = (jc / 4) * sum(
            x[i] * x[mod1(i + 1, N)] + y[i] * y[mod1(i + 1, N)]
            for i in 1:N
        )
        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=2))

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Fermionic APBC, N=$N, order 1 (= spin PBC ground state)" begin
        oracle = expectations_oracle(XY_EXPECTATIONS_PATH, "apbc_fermionic_n4_order1")

        registry, (a, a_dag) = create_fermionic_variables(1:N)
        # APBC: sign flip on boundary hopping
        ham = (jc / 2) * (
            sum(a_dag[i] * a[i + 1] + a_dag[i + 1] * a[i] for i in 1:N-1) -
            (a_dag[N] * a[1] + a_dag[1] * a[N])
        )
        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=1))

        @test result.objective ≈ -sqrt(2) atol = 1e-6
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Fermionic PBC, N=$N, order 1" begin
        oracle = expectations_oracle(XY_EXPECTATIONS_PATH, "pbc_fermionic_n4_order1")

        registry, (a, a_dag) = create_fermionic_variables(1:N)
        # Standard PBC
        ham = (jc / 2) * sum(
            a_dag[i] * a[mod1(i + 1, N)] + a_dag[mod1(i + 1, N)] * a[i]
            for i in 1:N
        )
        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=1))

        @test result.objective ≈ -1.0 atol = 1e-6
        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
