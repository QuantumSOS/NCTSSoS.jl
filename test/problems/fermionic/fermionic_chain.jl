# Kitaev Chain Sweet Spot via Majorana Fermions
#
# Requested benchmark:
#   H = i \sum_{j=1}^{N-1} γ_{2j} γ_{2j+1},   N = 4,
# with Majorana relations γ_j^2 = I and γ_j γ_k = -γ_k γ_j for j ≠ k.
#
# The current public framework does not yet solve this model directly through a
# generic Majorana/unipotent interface with complex coefficients. What it does
# support robustly is the equivalent Pauli representation obtained from the
# Jordan–Wigner image of the Majorana operators:
#
#   γ_{2j-1} = Z₁⋯Z_{j-1} X_j,
#   γ_{2j}   = -Z₁⋯Z_{j-1} Y_j,
#
# so that
#
#   i γ_{2j} γ_{2j+1} = X_j X_{j+1}.
#
# Therefore the Majorana sweet-spot chain maps exactly to the open XX chain
#
#   H = \sum_{j=1}^{N-1} X_j X_{j+1},
#
# whose ground-state energy is E₀ = -(N - 1) = -3 for N = 4.

using NCTSSoS, Test, JuMP, LinearAlgebra

const FERMIONIC_CHAIN_EXPECTATIONS_PATH = "expectations/fermionic_chain.toml"

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0,
    )
end

function _fermionic_chain_kron_all(mats::AbstractVector{<:AbstractMatrix})
    isempty(mats) && return Matrix{ComplexF64}(I, 1, 1)
    return reduce(kron, mats)
end

const _PAULI_X = ComplexF64[0 1; 1 0]
const _PAULI_Z = ComplexF64[1 0; 0 -1]
const _PAULI_I = Matrix{ComplexF64}(I, 2, 2)
const _JW_SP = ComplexF64[0 1; 0 0]
const _JW_SM = ComplexF64[0 0; 1 0]

function _jw_fermion_op(mode::Int, kind::Symbol, nmodes::Int)
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
    return _fermionic_chain_kron_all(mats)
end

function _pauli_site_op(site::Int, local_op::AbstractMatrix{<:Number}, nsites::Int)
    mats = [_PAULI_I for _ in 1:nsites]
    mats[site] = Matrix{ComplexF64}(local_op)
    return _fermionic_chain_kron_all(mats)
end

function _majorana_sweet_spot_hamiltonian_oracle(nsites::Int)
    a = [_jw_fermion_op(site, :annihilate, nsites) for site in 1:nsites]
    a_dag = [_jw_fermion_op(site, :create, nsites) for site in 1:nsites]

    gamma = Matrix{ComplexF64}[]
    for site in 1:nsites
        push!(gamma, a[site] + a_dag[site])
        push!(gamma, im * (a_dag[site] - a[site]))
    end

    dim = 2^nsites
    return sum(
        im * gamma[2 * site] * gamma[2 * site + 1]
        for site in 1:(nsites - 1);
        init=zeros(ComplexF64, dim, dim),
    )
end

function _pauli_xx_chain_hamiltonian_oracle(nsites::Int)
    dim = 2^nsites
    return sum(
        _pauli_site_op(site, _PAULI_X, nsites) * _pauli_site_op(site + 1, _PAULI_X, nsites)
        for site in 1:nsites-1;
        init=zeros(ComplexF64, dim, dim),
    )
end

@testset "Kitaev chain sweet spot (N=4, Majorana mapping)" begin
    N = 4

    H_majorana = _majorana_sweet_spot_hamiltonian_oracle(N)
    H_pauli = _pauli_xx_chain_hamiltonian_oracle(N)

    @test isapprox(H_majorana, H_pauli; atol=1e-12, rtol=0)

    exact_e0 = eigmin(Hermitian(H_majorana))
    @test exact_e0 ≈ -3.0 atol = 1e-12

    oracle = expectations_oracle(
        FERMIONIC_CHAIN_EXPECTATIONS_PATH,
        "kitaev_majorana_n4_order1",
    )

    registry, (x, _, _) = create_pauli_variables(1:N)
    ham = sum(1.0 * x[site] * x[site + 1] for site in 1:N-1)
    pop = polyopt(ham, registry)

    solver_config = SolverConfig(optimizer=SOLVER, order=1)
    result = cs_nctssos(pop, solver_config)

    @test result.objective ≈ exact_e0 atol = 1e-6
    @test result.objective ≈ oracle.opt atol = 1e-6
    @test reduce(vcat, result.moment_matrix_sizes) == oracle.sides
    @test result.n_unique_moment_matrix_elements == oracle.nuniq
end
