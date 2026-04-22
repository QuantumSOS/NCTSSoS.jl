# 1D Hubbard model tests
#
# Spin is represented by doubling the fermionic mode count: the up/down species
# share one FermionicAlgebra registry but live on disjoint mode ranges. This is
# the smallest honest way to encode the spinful Hubbard Hamiltonian with the
# current polynomial interface.
#
# Grand-canonical tests (no particle-number constraint) verify the unconstrained
# lower bound. The canonical test uses `moment_eq_constraints` to fix
# (N_up, N_dn) = (2, 2) via one-sided localizing constraints g|ψ⟩ = 0.
#
# The two-site periodic ring is still useful because its unconstrained ground
# state already sits in the half-filled sector. The four-site ring exercises the
# actual N=4, t=1, U=4 model both as a grand-canonical Hamiltonian (full-Fock
# lower bound) and as a canonical half-filling problem.

using Test, NCTSSoS, JuMP, LinearAlgebra
import NCTSSoS: degree

const HUBBARD_EXPECTATIONS_PATH = "expectations/hubbard.toml"

const _PAULI_Z = ComplexF64[1 0; 0 -1]
const _PAULI_I = Matrix{ComplexF64}(I, 2, 2)
const _JW_SP = ComplexF64[0 1; 0 0]
const _JW_SM = ComplexF64[0 0; 1 0]

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0,
    )
end

function _hubbard_kron_all(mats::AbstractVector{<:AbstractMatrix})
    isempty(mats) && return Matrix{ComplexF64}(I, 1, 1)
    return reduce(kron, mats)
end

function _hubbard_jw_fermion_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = Vector{Matrix{ComplexF64}}(undef, nmodes)
    for site in 1:nmodes
        if site < mode
            mats[site] = _PAULI_Z
        elseif site == mode
            mats[site] = kind === :annihilate ? _JW_SM : _JW_SP
        else
            mats[site] = _PAULI_I
        end
    end
    return _hubbard_kron_all(mats)
end

function _hubbard_fermion_word_matrix(word::AbstractVector{<:Integer}, nmodes::Int)
    dim = 2^nmodes
    acc = Matrix{ComplexF64}(I, dim, dim)
    for op in word
        mode = abs(Int(op))
        kind = op > 0 ? :annihilate : :create
        acc = acc * _hubbard_jw_fermion_op(mode, kind, nmodes)
    end
    return acc
end

function _hubbard_fermion_poly_matrix(p::Polynomial{FermionicAlgebra,T,C}, nmodes::Int) where {T<:Integer,C<:Number}
    dim = 2^nmodes
    acc = zeros(ComplexF64, dim, dim)
    for (coef, mono) in zip(coefficients(p), monomials(p))
        acc .+= ComplexF64(coef) * _hubbard_fermion_word_matrix(mono.word, nmodes)
    end
    return acc
end

function _hubbard_problem(nsites::Int; t::Real, U::Real, periodic::Bool)
    registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
        ("c_up", 1:nsites),
        ("c_dn", 1:nsites),
    ])

    bonds = periodic ? [(i, mod1(i + 1, nsites)) for i in 1:nsites] : [(i, i + 1) for i in 1:nsites-1]

    hopping = -t * sum(
        c_up_dag[i] * c_up[j] + c_up_dag[j] * c_up[i] +
        c_dn_dag[i] * c_dn[j] + c_dn_dag[j] * c_dn[i]
        for (i, j) in bonds
    )
    interaction = U * sum(
        (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
        for i in 1:nsites
    )

    return registry, (c_up, c_up_dag), (c_dn, c_dn_dag), hopping + interaction
end

function _sector_indices(nsites::Int; nup::Int, ndn::Int)
    nmodes = 2 * nsites
    upmask = (UInt(1) << nsites) - 1
    indices = Int[]
    for state in 0:(2^nmodes - 1)
        up_occ = count_ones(UInt(state) & upmask)
        dn_occ = count_ones(UInt(state) >> nsites)
        if up_occ == nup && dn_occ == ndn
            push!(indices, state + 1)
        end
    end
    return indices
end

function _sector_ground_state_energy(H::AbstractMatrix{<:Complex}, nsites::Int; nup::Int, ndn::Int)
    keep = _sector_indices(nsites; nup, ndn)
    return eigmin(Hermitian(H[keep, keep]))
end

@testset "1D Hubbard model (spinful via doubled fermionic modes)" begin
    @testset "On-site interaction stays quartic across spin species" begin
        _, (c_up, c_up_dag), (c_dn, c_dn_dag), _ = _hubbard_problem(2; t=1.0, U=4.0, periodic=true)

        double_occupancy = Polynomial((c_up_dag[1] * c_up[1]) * (c_dn_dag[1] * c_dn[1]))

        @test degree(double_occupancy) == 4
        @test length(double_occupancy.terms) == 1
    end

    @testset "Periodic two-site ring (exact at order 2)" begin
        N = 2
        t = 1.0
        U = 4.0
        oracle = expectations_oracle(HUBBARD_EXPECTATIONS_PATH, "periodic_n2_u4_order2")

        registry, _, _, ham = _hubbard_problem(N; t, U, periodic=true)
        exact_hamiltonian = _hubbard_fermion_poly_matrix(ham, 2 * N)
        exact_e0 = eigmin(Hermitian(exact_hamiltonian))
        half_filled_e0 = _sector_ground_state_energy(exact_hamiltonian, N; nup=1, ndn=1)

        @test exact_e0 ≈ oracle.opt atol = 1e-12
        @test exact_e0 ≈ half_filled_e0 atol = 1e-12

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=2))

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Periodic four-site ring (grand-canonical lower bound)" begin
        N = 4
        t = 1.0
        U = 4.0
        oracle = expectations_oracle(HUBBARD_EXPECTATIONS_PATH, "periodic_n4_u4_order2")

        registry, _, _, ham = _hubbard_problem(N; t, U, periodic=true)
        exact_hamiltonian = _hubbard_fermion_poly_matrix(ham, 2 * N)
        exact_e0 = eigmin(Hermitian(exact_hamiltonian))
        half_filled_e0 = _sector_ground_state_energy(exact_hamiltonian, N; nup=2, ndn=2)

        @test exact_e0 ≈ -3.4185507188738447 atol = 1e-12
        @test half_filled_e0 ≈ -2.102748483462074 atol = 1e-12
        @test exact_e0 < half_filled_e0

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=2))

        @test result.objective ≈ oracle.opt atol = 5e-6
        @test result.objective ≤ exact_e0 + 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end

    @testset "Periodic four-site ring (canonical half-filling via moment_eq)" begin
        N = 4
        t = 1.0
        U = 4.0
        oracle = expectations_oracle(HUBBARD_EXPECTATIONS_PATH, "periodic_n4_u4_order2_canonical")

        registry, (c_up, c_up_dag), (c_dn, c_dn_dag), ham = _hubbard_problem(N; t, U, periodic=true)
        exact_hamiltonian = _hubbard_fermion_poly_matrix(ham, 2 * N)
        half_filled_e0 = _sector_ground_state_energy(exact_hamiltonian, N; nup=2, ndn=2)

        # Build particle-number constraints: N_up = 2, N_dn = 2
        n_up_total = 1.0 * sum(c_up_dag[i] * c_up[i] for i in 1:N)
        n_dn_total = 1.0 * sum(c_dn_dag[i] * c_dn[i] for i in 1:N)
        cons_nup = n_up_total - 2.0 * one(ham)
        cons_ndn = n_dn_total - 2.0 * one(ham)

        # Canonical: use moment_eq_constraints for state sector constraints
        pop = polyopt(ham, registry; moment_eq_constraints=[cons_nup, cons_ndn])
        result = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=2))

        @test result.objective ≈ oracle.opt atol = 1e-6
        @test result.objective ≤ half_filled_e0 + 1e-6  # valid lower bound
        @test result.objective > -3.5  # must not drop to full-Fock GS
        @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
        @test result.n_unique_moment_matrix_elements == oracle.nuniq
    end
end
