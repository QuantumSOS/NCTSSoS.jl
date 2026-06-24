# Tests for sparse Pauli spin-chain bases and charge-compatible sign symmetry.

using Test, NCTSSoS

function _pauli_support(mono)
    return sort!(Int[NCTSSoS._pauli_site(idx) for idx in mono.word])
end

function _is_periodic_contiguous(sites::Vector{Int}, n::Int)
    isempty(sites) && return true
    target = Set(sites)
    width = length(sites)
    return any(Set(mod1(start + offset, n) for offset in 0:(width - 1)) == target for start in 1:n)
end

@testset "Pauli contiguous chain basis" begin
    reg4, (σx4, _, _) = create_pauli_variables(1:4)

    @test length(pauli_contiguous_chain_basis(reg4, 0)) == 1
    @test length(pauli_contiguous_chain_basis(reg4, 1)) == 13
    @test length(pauli_contiguous_chain_basis(reg4, 2)) == 49
    @test length(pauli_contiguous_chain_basis(reg4, 3)) == 157
    @test length(pauli_contiguous_chain_basis(reg4, 4)) == 238
    @test length(pauli_contiguous_chain_basis(reg4, 2; periodic=false)) == 40

    basis = pauli_contiguous_chain_basis(reg4, 3)
    @test one(σx4[1]) in basis
    @test all(mono -> degree(mono) <= 3, basis)
    @test all(mono -> _is_periodic_contiguous(_pauli_support(mono), 4), basis)

    reg10, _ = create_pauli_variables(1:10)
    @test length(pauli_contiguous_chain_basis(reg10, 2)) == 1 + 10 * (3 + 9)

    reg100, _ = create_pauli_variables(1:100)
    @test length(pauli_contiguous_chain_basis(reg100, 4)) == 1 + 100 * (3 + 9 + 27 + 81)

    @test_throws ArgumentError pauli_contiguous_chain_basis(reg4, -1)
end

@testset "Pauli chain basis closure under spatial and sign symmetries" begin
    n = 6
    reg, (σx, σy, σz) = create_pauli_variables(1:n)
    basis = pauli_contiguous_chain_basis(reg, 4)
    lookup = Set(basis)

    translation = pauli_site_permutation([2:n; 1])
    reflection = pauli_site_permutation(reverse(1:n))
    sign = pauli_sign_symmetry(n; integer_type=eltype(σx[1].word))

    for g in (translation, reflection, sign), mono in basis
        _, image = NCTSSoS._act_monomial(g, mono)
        @test image in lookup
    end

    @test NCTSSoS._act_monomial(sign, σx[1]) == (-1, σx[1])
    @test NCTSSoS._act_monomial(sign, σy[1]) == (-1, σy[1])
    @test NCTSSoS._act_monomial(sign, σz[1]) == (1, σz[1])
end

@testset "Sparse Pauli charge words follow the supplied chain basis" begin
    n = 6
    reg, _ = create_pauli_variables(1:n)
    basis = pauli_contiguous_chain_basis(reg, 4)

    charge_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=n, max_degree=4),
        nothing,
    )
    blocks = collect(Iterators.flatten(charge_groups))

    @test sum(size(block.row_basis, 1) for block in blocks) == length(basis)
    @test Set(block.label.charge for block in blocks) == Set(-4:4)
    @test all(block -> block.provenance == :charge_sector, blocks)

    sign_group = CliffordSymmetryGroup(
        pauli_sign_symmetry(n; integer_type=eltype(basis[1].word));
        nqubits=n,
        integer_type=eltype(basis[1].word),
    )
    signed_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=n, max_degree=4),
        sign_group,
    )
    signed_blocks = collect(Iterators.flatten(signed_groups))

    @test sum(size(block.row_basis, 1) for block in signed_blocks) == length(basis)
    @test all(block -> block.label.group_order == 2, signed_blocks)
    @test all(block -> block.provenance == :charge_sector, signed_blocks)
end


# Tests for translation-invariant Pauli chain relaxations.

using Test, NCTSSoS, JuMP

if !@isdefined(SOLVER)
    using COSMO
    const SOLVER = optimizer_with_attributes(
        COSMO.Optimizer,
        "verbose" => false,
        "eps_abs" => 1e-8,
        "eps_rel" => 1e-8,
        "max_iter" => 50_000,
    )
end

if !@isdefined(flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

@testset "Translation-invariant Pauli chain relaxation" begin
    quiet(f) = redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            f()
        end
    end

    @testset "contiguous basis and Heisenberg helpers" begin
        n = 8
        registry, ops = create_pauli_variables(1:n)

        basis = pauli_contiguous_chain_basis(ops, 2)
        hamiltonian = heisenberg_chain_hamiltonian(ops)

        @test length(basis) == 1 + n * (3 + 9)
        @test length(terms(hamiltonian)) == 3n
        @test polyopt(hamiltonian, registry).objective == hamiltonian
    end

    @testset "normalized momentum blocks and N=100-size block structure" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        mp, report = pauli_translation_invariant_moment_relaxation(pop, ops, 1; sign_symmetry=false)
        identity_cross = mp.constraints[1][2][1, 2]
        @test only(coefficients(identity_cross)) ≈ sqrt(n) + 0im atol = 1e-12
        @test report.psd_block_sizes == [8, 6, 6]
        @test report.real_moment_matrix

        n_large = 20
        registry_large, ops_large = create_pauli_variables(1:n_large)
        pop_large = polyopt(heisenberg_chain_hamiltonian(ops_large), registry_large)
        _, report_large = pauli_translation_invariant_moment_relaxation(pop_large, ops_large, 4)

        @test report_large.basis_size == 1 + n_large * sum(3^ℓ for ℓ in 1:4)
        @test report_large.orbit_basis_size == 1 + sum(3^ℓ for ℓ in 1:4)
        @test maximum(report_large.psd_block_sizes) == 62
        @test length(report_large.psd_block_sizes) == 4 * (fld(n_large, 2) + 1)
    end

    @testset "guardrails reject invalid reductions" begin
        registry, ops = create_pauli_variables(1:4)
        heisenberg_pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(heisenberg_pop, ops, 0)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(heisenberg_pop, ops, 1; momenta=[1, 2], real_moment_matrix=false)

        field_pop = polyopt(sum(ops[1]), registry)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(field_pop, ops, 1)
        @test pauli_translation_invariant_moment_relaxation(field_pop, ops, 1; sign_symmetry=false)[2].psd_block_sizes == [8, 6, 6]

        σx, σy, σz = ops
        @test_throws ArgumentError pauli_contiguous_chain_basis((σy, σx, σz), 1)

        registry8, ops8 = create_pauli_variables(1:8)
        mismatched_pop = polyopt(ops8[1][6] * ops8[1][7], registry8)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(mismatched_pop, ops, 2; sign_symmetry=false)

        registry100, ops100 = create_pauli_variables(1:100)
        type_mismatched_pop = polyopt(ops100[1][1] * ops100[1][2], registry100)
        @test_throws ArgumentError pauli_translation_invariant_moment_relaxation(type_mismatched_pop, ops, 2; sign_symmetry=false)
    end

    @testset "small chain agrees with dense order-1 relaxation" begin
        n = 4
        registry, ops = create_pauli_variables(1:n)
        pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)

        dense = quiet() do
            cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=1); dualize=false)
        end
        reduced = quiet() do
            pauli_translation_invariant_nctssos(pop, ops, 1, SOLVER; dualize=false)
        end

        @test termination_status(dense.model) == JuMP.MOI.OPTIMAL
        @test termination_status(reduced.model) == JuMP.MOI.OPTIMAL
        @test reduced.objective ≈ dense.objective atol = 1e-6
        @test maximum(reduced.report.psd_block_sizes) < only(flatten_sizes(dense.moment_matrix_sizes))
    end
end
