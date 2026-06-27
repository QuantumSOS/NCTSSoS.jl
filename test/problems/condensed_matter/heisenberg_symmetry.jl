# Public API regression for Clifford symmetry reduction on a 2-site Heisenberg model.

using Test, NCTSSoS, JuMP, SparseArrays

if !@isdefined(SOLVER)
    using COSMO
    const SOLVER = optimizer_with_attributes(
        COSMO.Optimizer,
        "verbose" => false,
        "eps_abs" => 1e-8,
        "eps_rel" => 1e-8,
        "eps_prim_inf" => 1e-6,
        "eps_dual_inf" => 1e-6,
        "max_iter" => 50_000,
        "rho" => 1.0,
        "adaptive_rho" => true,
        "alpha" => 1.0,
        "scaling" => 10,
    )
end

if !@isdefined(flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes)
end

@testset "2-site Heisenberg model with SWAP Clifford symmetry" begin
    registry, (σx, σy, σz) = create_pauli_variables(1:2)
    heisenberg = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
    pop = polyopt(heisenberg, registry)
    basis = [one(σx[1]); σx; σy; σz]

    plain_config = SolverConfig(
        optimizer=SOLVER,
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
    )
    symmetric_config = SolverConfig(
        optimizer=SOLVER,
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2)),
    )

    plain = cs_nctssos(pop, plain_config)
    symmetric = cs_nctssos(pop, symmetric_config)

    @test plain.objective ≈ -0.75 atol = 1e-6
    @test symmetric.objective ≈ plain.objective atol = 1e-6
    @test !isnothing(symmetric.symmetry)
    @test symmetric.symmetry.group_order == 2
    @test symmetric.symmetry.psd_block_sizes == [4, 3]
    @test maximum(symmetric.symmetry.psd_block_sizes) < only(flatten_sizes(plain.moment_matrix_sizes))
end

@testset "4-site Heisenberg ring with charge/spatial/singlet symmetry" begin
    n_ring = 4
    registry, (σx, σy, σz) = create_pauli_variables(1:n_ring)
    heisenberg = sum(
        ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n_ring)]
        for op in (σx, σy, σz) for i in 1:n_ring
    )
    pop = polyopt(heisenberg, registry)

    translation = pauli_site_permutation([2, 3, 4, 1])
    reflection = pauli_site_permutation([4, 3, 2, 1])
    symmetry = SymmetrySpec(
        [translation, reflection];
        pauli_charge=PauliChargeSectorSpec(nqubits=n_ring),
        pauli_singlet=PauliSingletConstraintSpec(nqubits=n_ring),
    )

    plain = cs_nctssos(
        pop,
        SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        ),
    )
    reduced = cs_nctssos(
        pop,
        SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=symmetry,
        ),
    )

    @test reduced.objective ≈ plain.objective atol = 1e-6
    @test !isnothing(reduced.symmetry)
    @test reduced.symmetry.group_order == 8
    @test maximum(reduced.symmetry.psd_block_sizes) < only(flatten_sizes(plain.moment_matrix_sizes))
    @test any(label -> label isa PauliChargeBlockLabel && label.charge == 0, reduced.symmetry.block_labels)
end

@testset "16-site order-2 Heisenberg symmetry size evidence" begin
    n_ring = 16
    registry, (σx, σy, σz) = create_pauli_variables(1:n_ring)
    heisenberg = sum(
        ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n_ring)]
        for op in (σx, σy, σz) for i in 1:n_ring
    )
    pop = polyopt(heisenberg, registry)

    translation = pauli_site_permutation([2:16; 1])
    reflection = pauli_site_permutation(reverse(1:16))
    symmetry = SymmetrySpec(
        [translation, reflection];
        pauli_charge=PauliChargeSectorSpec(nqubits=n_ring),
        pauli_singlet=PauliSingletConstraintSpec(nqubits=n_ring),
        offblock_check=:off,
    )
    cfg = SolverConfig(
        optimizer=nothing,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )
    sparsity = compute_sparsity(pop, cfg)
    basis = only(sparsity.corr_sparsity.clq_mom_mtx_bases)

    support_domain = NCTSSoS._symmetry_domain(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
    )
    sw_group = CliffordSymmetryGroup(
        symmetry.clifford_generators;
        nqubits=n_ring,
        integer_type=eltype(basis).parameters[2],
        domain=support_domain,
    )
    charge_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        symmetry.pauli_charge,
        sw_group,
    )
    block_sizes = [size(block.row_basis, 1) for group in charge_groups for block in group]
    block_variable_count = sum(size * (size + 1) ÷ 2 for size in block_sizes)

    @test all(issparse(block.row_basis) for group in charge_groups for block in group)
    @test length(basis) == 1129
    @test 1129 * 1130 ÷ 2 == 637885
    @test maximum(block_sizes) <= 24
    @test block_variable_count == 5108
    @test block_variable_count < 637885 ÷ 100
    @test Set(block.label.charge for group in charge_groups for block in group) == Set(-2:2)
end

@testset "16-site sparse degree-4 Heisenberg charge/sign size evidence" begin
    n_ring = 16
    registry, _ = create_pauli_variables(1:n_ring)
    basis = pauli_contiguous_chain_basis(registry, 4)

    translation = pauli_site_permutation([2:n_ring; 1])
    reflection = pauli_site_permutation(reverse(1:n_ring))
    sign = pauli_sign_symmetry(n_ring; integer_type=eltype(basis[1].word))
    sw_group = CliffordSymmetryGroup(
        [translation, reflection, sign];
        nqubits=n_ring,
        integer_type=eltype(basis[1].word),
    )
    charge_groups = NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=n_ring, max_degree=4),
        sw_group,
    )
    block_sizes = [size(block.row_basis, 1) for group in charge_groups for block in group]
    block_variable_count = sum(size * (size + 1) ÷ 2 for size in block_sizes)

    @test length(basis) == 1 + n_ring * (3 + 9 + 27 + 81)
    @test length(sw_group) == 64
    @test maximum(block_sizes) == 19
    @test block_variable_count < length(basis) * (length(basis) + 1) ÷ 200
    @test Set(block.label.charge for group in charge_groups for block in group) == Set(-4:4)
end

@testset "32-site Heisenberg large charge-sector fallback" begin
    n_ring = 32
    registry, (σx, σy, σz) = create_pauli_variables(1:n_ring)
    heisenberg = sum(
        ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n_ring)]
        for op in (σx, σy, σz) for i in 1:n_ring
    )
    pop = polyopt(heisenberg, registry)

    translation = pauli_site_permutation([2:n_ring; 1])
    reflection = pauli_site_permutation(reverse(1:n_ring))
    symmetry = SymmetrySpec(
        [translation, reflection];
        pauli_charge=PauliChargeSectorSpec(nqubits=n_ring),
        pauli_singlet=PauliSingletConstraintSpec(nqubits=n_ring),
        offblock_check=:off,
    )
    cfg = SolverConfig(
        optimizer=nothing,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )
    sparsity = compute_sparsity(pop, cfg)
    mp, report = NCTSSoS.moment_relax_symmetric(
        pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
        symmetry,
    )

    @test report.basis_half_size == 4561
    @test report.group_order == 64
    @test report.psd_block_sizes == [16, 17, 34, 17, 16]
    @test all(==(:charge_orbit_representative), report.block_provenance)
    @test report.invariant_moment_count == 769
    @test length(mp.constraints) == 88
    @test length(mp.linear.zero_constraints) == 83
end
