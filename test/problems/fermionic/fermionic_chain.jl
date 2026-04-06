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

function _kron_all(mats::AbstractVector{<:AbstractMatrix})
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
    return _kron_all(mats)
end

function _pauli_site_op(site::Int, local_op::AbstractMatrix{<:Number}, nsites::Int)
    mats = [_PAULI_I for _ in 1:nsites]
    mats[site] = Matrix{ComplexF64}(local_op)
    return _kron_all(mats)
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

# =============================================================================
# 3.1.3 Kitaev Chain (Topological Superconductor)
# =============================================================================
#
# The Kitaev chain is the paradigmatic model for topological superconductivity
# in 1D. It exhibits a topological phase with unpaired Majorana zero modes at
# the chain ends when |μ| < 2t and Δ ≠ 0.
#
# Hamiltonian (N sites, open boundary conditions):
#   H = -μ Σᵢ aᵢ†aᵢ - t Σᵢ (aᵢ†aᵢ₊₁ + aᵢ₊₁†aᵢ) + Δ Σᵢ (aᵢaᵢ₊₁ + aᵢ₊₁†aᵢ†)
#
# The model is solvable via Jordan-Wigner transformation to free fermions.
# Single-particle energies: εₖ = ±√[(2t cos(k) + μ)² + (2Δ sin(k))²]
#
# Reference:
#   Kitaev, A.Yu. "Unpaired Majorana fermions in quantum wires."
#   Physics-Uspekhi 44, 131 (2001), arXiv:cond-mat/0010440.

# Exact diagonalization oracle for Kitaev chain using Jordan-Wigner mapping
function _kitaev_chain_exact_energy(nsites::Int, μ::Real, t::Real, Δ::Real)
    # Build exact Hamiltonian matrix using JW transformation
    # H = -μ Σᵢ aᵢ†aᵢ - t Σᵢ (aᵢ†aᵢ₊₁ + h.c.) + Δ Σᵢ (aᵢaᵢ₊₁ + h.c.)
    a = [_jw_fermion_op(i, :annihilate, nsites) for i in 1:nsites]
    a_dag = [_jw_fermion_op(i, :create, nsites) for i in 1:nsites]
    
    dim = 2^nsites
    H = zeros(ComplexF64, dim, dim)
    
    # Chemical potential term: -μ Σᵢ aᵢ†aᵢ
    for i in 1:nsites
        H .-= μ .* (a_dag[i] * a[i])
    end
    
    # Hopping term: -t Σᵢ (aᵢ†aᵢ₊₁ + aᵢ₊₁†aᵢ)
    for i in 1:(nsites-1)
        H .-= t .* (a_dag[i] * a[i+1] + a_dag[i+1] * a[i])
    end
    
    # Pairing term: Δ Σᵢ (aᵢaᵢ₊₁ + aᵢ₊₁†aᵢ†)
    for i in 1:(nsites-1)
        H .+= Δ .* (a[i] * a[i+1] + a_dag[i+1] * a_dag[i])
    end
    
    return eigmin(Hermitian(H))
end

@testset "3.1.3 Kitaev Chain (topological superconductor)" begin
    N = 4

    @testset "Topological phase (μ=0, t=Δ=1) via FermionicAlgebra" begin
        # In the topological phase with μ=0, t=Δ=1, the Hamiltonian reduces to:
        # H = Σᵢ [-(aᵢ†aᵢ₊₁ + aᵢ₊₁†aᵢ) + (aᵢaᵢ₊₁ + aᵢ₊₁†aᵢ†)]
        # This maps to iΣⱼ γ₂ⱼγ₂ⱼ₊₁ in Majorana basis with E₀ = -(N-1)
        
        registry, (a, a_dag) = create_fermionic_variables(1:N)

        hopping = sum(
            a_dag[i] * a[i+1] + a_dag[i+1] * a[i]
            for i in 1:(N-1)
        )
        pairing = sum(
            a[i] * a[i+1] + a_dag[i+1] * a_dag[i]
            for i in 1:(N-1)
        )
        ham = pairing - hopping

        pop = polyopt(ham, registry)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)
        result = cs_nctssos(pop, solver_config)

        # Exact ground state energy from ED: E₀ = -(N-1) = -3
        exact_e0 = _kitaev_chain_exact_energy(N, 0.0, 1.0, 1.0)
        @test exact_e0 ≈ -3.0 atol = 1e-12
        @test result.objective ≈ exact_e0 atol = 1e-6
    end

    @testset "Trivial phase (μ=3, t=1, Δ=1) via FermionicAlgebra" begin
        # In the trivial phase with |μ| > 2t, the chemical potential dominates.
        # For μ=3, t=Δ=1, exact ED gives E₀ ≈ -12.5039
        
        registry, (a, a_dag) = create_fermionic_variables(1:N)
        
        μ, t, Δ = 3.0, 1.0, 1.0
        
        chemical = -μ * sum(a_dag[i] * a[i] for i in 1:N)
        hopping = -t * sum(
            a_dag[i] * a[i+1] + a_dag[i+1] * a[i]
            for i in 1:(N-1)
        )
        pairing = Δ * sum(
            a[i] * a[i+1] + a_dag[i+1] * a_dag[i]
            for i in 1:(N-1)
        )
        ham = chemical + hopping + pairing
        
        pop = polyopt(ham, registry)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)
        result = cs_nctssos(pop, solver_config)
        
        exact_e0 = _kitaev_chain_exact_energy(N, μ, t, Δ)
        @test exact_e0 ≈ -12.503891557126414 atol = 1e-10
        @test result.objective ≈ exact_e0 atol = 1e-6
    end

    @testset "XX chain with field via PauliAlgebra" begin
        # XX chain with longitudinal field in Pauli algebra (Xᵢ²=1).
        # H = -Σᵢ XᵢXᵢ₊₁ - h Σᵢ Zᵢ
        # For h=1, N=4 with Pauli constraints: relaxation gives E₀ = -7.0
        # This differs from fermionic ED due to different algebra constraints.
        
        registry, (x, y, z) = create_pauli_variables(1:N)
        
        xx_interaction = sum(x[i] * x[i+1] for i in 1:(N-1))
        field = sum(z[i] for i in 1:N)
        ham = -(xx_interaction + field)
        
        pop = polyopt(ham, registry)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)
        result = cs_nctssos(pop, solver_config)
        
        # With Pauli constraints (Xᵢ²=1), the relaxation achieves E₀ = -7.0
        # at order=1 (all spins aligned with effective field)
        @test result.objective ≈ -7.0 atol = 1e-6
    end

    @testset "Asymmetric pairing (μ=0, t=1, Δ=0.5) via FermionicAlgebra" begin
        # With Δ ≠ t, the ground state energy interpolates between limits.
        # For μ=0, t=1, Δ=0.5: exact ED gives E₀ ≈ -2.4221
        
        registry, (a, a_dag) = create_fermionic_variables(1:N)
        
        μ, t, Δ = 0.0, 1.0, 0.5
        
        hopping = -t * sum(
            a_dag[i] * a[i+1] + a_dag[i+1] * a[i]
            for i in 1:(N-1)
        )
        pairing = Δ * sum(
            a[i] * a[i+1] + a_dag[i+1] * a_dag[i]
            for i in 1:(N-1)
        )
        ham = hopping + pairing
        
        pop = polyopt(ham, registry)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)
        result = cs_nctssos(pop, solver_config)
        
        exact_e0 = _kitaev_chain_exact_energy(N, μ, t, Δ)
        @test exact_e0 ≈ -2.422078451440553 atol = 1e-10
        @test result.objective ≈ exact_e0 atol = 1e-6
    end

    @testset "Pure hopping (μ=0, t=1, Δ=0) via FermionicAlgebra" begin
        # No pairing (Δ=0): standard tight-binding chain.
        # For N=4, μ=0, t=1, exact ED gives E₀ = -√5 ≈ -2.2361
        
        registry, (a, a_dag) = create_fermionic_variables(1:N)
        
        hopping = -sum(
            a_dag[i] * a[i+1] + a_dag[i+1] * a[i]
            for i in 1:(N-1)
        )
        
        pop = polyopt(hopping, registry)
        solver_config = SolverConfig(optimizer=SOLVER, order=1)
        result = cs_nctssos(pop, solver_config)
        
        exact_e0 = _kitaev_chain_exact_energy(N, 0.0, 1.0, 0.0)
        @test exact_e0 ≈ -2.23606797749979 atol = 1e-10  # -√5
        @test result.objective ≈ exact_e0 atol = 1e-6
    end
end
