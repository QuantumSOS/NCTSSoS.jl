# 1D transverse-field Ising model tests
#
# The package already had a small Pauli-form script here, but it was never part
# of the curated suite. This version keeps the cheap open/periodic Pauli checks
# and adds the requested Jordan–Wigner fermionic form for the open chain.

using Test, NCTSSoS, JuMP, LinearAlgebra

const ISING_EXPECTATIONS_PATH = "expectations/ising.toml"

const _PAULI_X = ComplexF64[0 1; 1 0]
const _PAULI_Z = ComplexF64[1 0; 0 -1]
const _PAULI_I = Matrix{ComplexF64}(I, 2, 2)
const _JW_SP = ComplexF64[0 1; 0 0]
const _JW_SM = ComplexF64[0 0; 1 0]

function _ising_kron_all(mats::AbstractVector{<:AbstractMatrix})
    isempty(mats) && return Matrix{ComplexF64}(I, 1, 1)
    return reduce(kron, mats)
end

function _ising_pauli_site_op(site::Int, local_op::AbstractMatrix{<:Number}, nsites::Int)
    mats = [_PAULI_I for _ in 1:nsites]
    mats[site] = Matrix{ComplexF64}(local_op)
    return _ising_kron_all(mats)
end

function _tfim_pauli_matrix(
    nsites::Int;
    J::Real,
    h::Real,
    coupling_axis::Symbol,
    field_axis::Symbol,
    periodic::Bool,
)
    coupling_op = coupling_axis === :x ? _PAULI_X : _PAULI_Z
    field_op = field_axis === :x ? _PAULI_X : _PAULI_Z

    dim = 2^nsites
    hamiltonian = zeros(ComplexF64, dim, dim)

    for site in 1:(periodic ? nsites : nsites - 1)
        next_site = mod1(site + 1, nsites)
        hamiltonian .+= -(J / 4) * _ising_pauli_site_op(site, coupling_op, nsites) * _ising_pauli_site_op(next_site, coupling_op, nsites)
    end

    for site in 1:nsites
        hamiltonian .+= -(h / 2) * _ising_pauli_site_op(site, field_op, nsites)
    end

    return hamiltonian
end

function _ising_jw_fermion_op(mode::Int, kind::Symbol, nmodes::Int)
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
    return _ising_kron_all(mats)
end

function _ising_fermion_word_matrix(word::AbstractVector{<:Integer}, nmodes::Int)
    dim = 2^nmodes
    acc = Matrix{ComplexF64}(I, dim, dim)
    for op in word
        mode = abs(Int(op))
        kind = op > 0 ? :annihilate : :create
        acc = acc * _ising_jw_fermion_op(mode, kind, nmodes)
    end
    return acc
end

function _ising_fermion_poly_matrix(p::Polynomial{FermionicAlgebra,T,C}, nmodes::Int) where {T<:Integer,C<:Number}
    dim = 2^nmodes
    acc = zeros(ComplexF64, dim, dim)
    for (coef, mono) in zip(coefficients(p), monomials(p))
        acc .+= ComplexF64(coef) * _ising_fermion_word_matrix(mono.word, nmodes)
    end
    return acc
end

function _tfim_open_fermionic_hamiltonian(a, a_dag; J::Real, h::Real)
    coupling = sum(
        (J / 4) * (a_dag[site] - a[site]) * (a_dag[site + 1] + a[site + 1])
        for site in 1:(length(a) - 1)
    )
    field = sum(
        begin
            number_op = a_dag[site] * a[site]
            -(h / 2) * (2.0 * number_op - one(number_op))
        end
        for site in eachindex(a)
    )
    return coupling + field
end

@testset "1D transverse-field Ising model" begin
    N = 3
    J = 1.0
    h = 2.0

    pauli_open_oracle = expectations_oracle(ISING_EXPECTATIONS_PATH, "open_pauli_n3_order2")
    pauli_periodic_oracle = expectations_oracle(ISING_EXPECTATIONS_PATH, "periodic_pauli_n3_order2")
    fermionic_open_oracle = expectations_oracle(ISING_EXPECTATIONS_PATH, "open_fermionic_n3_order1")

    @testset "Pauli form (open chain)" begin
        exact_hamiltonian = _tfim_pauli_matrix(N; J, h, coupling_axis=:z, field_axis=:x, periodic=false)
        exact_e0 = eigmin(Hermitian(exact_hamiltonian))
        @test exact_e0 ≈ pauli_open_oracle.opt atol = 1e-12

        registry, (sx, _, sz) = create_pauli_variables(1:N)
        ham = sum(-(J / 4) * sz[site] * sz[site + 1] for site in 1:N-1) +
              sum(-(h / 2) * sx[site] for site in 1:N)

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=2))
        @test result.objective ≈ pauli_open_oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == pauli_open_oracle.sides
        @test result.n_unique_moment_matrix_elements == pauli_open_oracle.nuniq
    end

    @testset "Pauli form (periodic chain)" begin
        exact_hamiltonian = _tfim_pauli_matrix(N; J, h, coupling_axis=:z, field_axis=:x, periodic=true)
        exact_e0 = eigmin(Hermitian(exact_hamiltonian))
        @test exact_e0 ≈ pauli_periodic_oracle.opt atol = 1e-12

        registry, (sx, _, sz) = create_pauli_variables(1:N)
        ham = sum(-(J / 4) * sz[site] * sz[mod1(site + 1, N)] for site in 1:N) +
              sum(-(h / 2) * sx[site] for site in 1:N)

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=2))
        @test result.objective ≈ pauli_periodic_oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == pauli_periodic_oracle.sides
        @test result.n_unique_moment_matrix_elements == pauli_periodic_oracle.nuniq
    end

    @testset "Jordan–Wigner fermionic form (open chain)" begin
        # The ZZ + X convention above is globally y-rotated to XX + Z before the
        # Jordan–Wigner map. The spectrum is unchanged, but the fermionic form is
        # cleanly quadratic only in the rotated basis.
        rotated_exact_hamiltonian = _tfim_pauli_matrix(N; J, h, coupling_axis=:x, field_axis=:z, periodic=false)
        exact_e0 = eigmin(Hermitian(rotated_exact_hamiltonian))

        @test exact_e0 ≈ pauli_open_oracle.opt atol = 1e-12
        @test exact_e0 ≈ fermionic_open_oracle.opt atol = 1e-12

        registry, (a, a_dag) = create_fermionic_variables(1:N)
        ham = _tfim_open_fermionic_hamiltonian(a, a_dag; J, h)
        fermionic_matrix = _ising_fermion_poly_matrix(ham, N)

        @test isapprox(fermionic_matrix, rotated_exact_hamiltonian; atol=1e-12, rtol=0)

        result = cs_nctssos(polyopt(ham, registry), SolverConfig(optimizer=SOLVER, order=1))
        @test result.objective ≈ fermionic_open_oracle.opt atol = 1e-6
        @test reduce(vcat, result.moment_matrix_sizes) == fermionic_open_oracle.sides
        @test result.n_unique_moment_matrix_elements == fermionic_open_oracle.nuniq
    end
end
