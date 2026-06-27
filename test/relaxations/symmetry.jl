# test/relaxations/symmetry.jl
# Tests: fail-fast behavior for dense symmetry MVP

using Test, NCTSSoS, LinearAlgebra

if !@isdefined(SOLVER)
    using JuMP, COSMO
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

if !@isdefined(expectations_oracle)
    include("../Expectations.jl")
    using .TestExpectations: expectations_oracle
end

const SW = NCTSSoS.SymbolicWedderburn

function _capture_exception(f)
    try
        f()
        return nothing
    catch err
        return err
    end
end

function _create_chsh_problem()
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    return polyopt(-f, reg), reg, x, y
end

function _create_chsh_symmetry(x, y)
    alice_swap = SignedPermutation(
        x[1].word[1] => x[2].word[1],
        x[2].word[1] => x[1].word[1],
        y[2].word[1] => (-1, y[2].word[1]),
    )
    bob_swap = SignedPermutation(
        x[2].word[1] => (-1, x[2].word[1]),
        y[1].word[1] => y[2].word[1],
        y[2].word[1] => y[1].word[1],
    )
    party_swap = SignedPermutation(
        x[1].word[1] => y[1].word[1],
        x[2].word[1] => y[2].word[1],
        y[1].word[1] => x[1].word[1],
        y[2].word[1] => x[2].word[1],
    )
    return SymmetrySpec(alice_swap, bob_swap, party_swap)
end

function _create_spinless_fermion_swap_problem()
    reg, (a, a_dag) = create_fermionic_variables(1:2)
    ham = -(a_dag[1] * a[2] + a_dag[2] * a[1])
    basis = [one(a[1]), a[1], a[2], a_dag[1], a_dag[2]]
    symmetry = SymmetrySpec(
        fermionic_generators=[FermionicModePermutation(1 => 2, 2 => 1)],
        sector=FermionicSectorSpec(split_parity=true, split_number=true),
    )
    return polyopt(ham, reg), reg, a, a_dag, basis, symmetry
end

@testset "Symmetry MVP Guards" begin
    @testset "non-invariant objective errors" begin
        _, reg, x, y = _create_chsh_problem()
        pop = polyopt(-(x[1] * y[1]), reg)
        T = typeof(x[1].word[1])
        spec = NCTSSoS._convert_symmetry_spec(T, _create_chsh_symmetry(x, y))
        domain = sort!(T[x[1].word[1], x[2].word[1], y[1].word[1], y[2].word[1]])
        group = NCTSSoS._enumerate_symmetry_group(spec, domain)

        err = _capture_exception(() -> NCTSSoS._check_symmetry_invariance(pop, group))

        @test err isa ArgumentError
        @test occursin("objective", sprint(showerror, err))
    end

    @testset "basis closure errors" begin
        _, _, x, y = _create_chsh_problem()
        T = typeof(x[1].word[1])
        spec = NCTSSoS._convert_symmetry_spec(T, _create_chsh_symmetry(x, y))
        domain = sort!(T[x[1].word[1], x[2].word[1], y[1].word[1], y[2].word[1]])
        group = NCTSSoS._enumerate_symmetry_group(spec, domain)
        bad_basis = [one(x[1]), x[1], y[1]]

        err = _capture_exception(() -> NCTSSoS._check_basis_closure("test basis", bad_basis, group))

        @test err isa ArgumentError
        @test occursin("closure", sprint(showerror, err))
    end

    @testset "SignedPermutation validates malformed pair inputs" begin
        err = _capture_exception(() -> SignedPermutation(1.0 => 2))
        @test err isa ArgumentError
        @test occursin("source index must be an integer", sprint(showerror, err))

        err = _capture_exception(() -> SignedPermutation(1 => 2.0))
        @test err isa ArgumentError
        @test occursin("target index must be an integer", sprint(showerror, err))

        err = _capture_exception(() -> SignedPermutation(1 => (1, 2, 3)))
        @test err isa ArgumentError
        @test occursin("(sign, target_index)", sprint(showerror, err))

        err = _capture_exception(() -> SignedPermutation(1 => (1.0, 2)))
        @test err isa ArgumentError
        @test occursin("signs must be integer ±1", sprint(showerror, err))

        err = _capture_exception(() -> SignedPermutation(1 => (0, 2)))
        @test err isa ArgumentError
        @test occursin("±1", sprint(showerror, err))
    end

    @testset "CliffordSymmetry named gates act on Pauli words" begin
        _, (σx, σy, σz) = create_pauli_variables(1:2)
        act(g, mono) = NCTSSoS._act_monomial(g, mono)
        mono(poly) = only(monomials(poly))

        h = CliffordSymmetry(:H, 1)
        @test h.nqubits == 1
        @test act(h, σx[1]) == (1, σz[1])
        @test act(h, σy[1]) == (-1, σy[1])
        @test act(h, σz[1]) == (1, σx[1])

        s = CliffordSymmetry(:S, 1)
        @test act(s, σx[1]) == (1, σy[1])
        @test act(s, σy[1]) == (-1, σx[1])
        @test act(s, σz[1]) == (1, σz[1])

        cnot = CliffordSymmetry(:CNOT, 1, 2)
        @test cnot.nqubits == 2
        @test act(cnot, σx[1]) == (1, mono(σx[1] * σx[2]))
        @test act(cnot, σy[1]) == (1, mono(σy[1] * σx[2]))
        @test act(cnot, σz[2]) == (1, mono(σz[1] * σz[2]))
        @test act(cnot, σy[2]) == (1, mono(σz[1] * σy[2]))
        @test act(cnot, mono(σx[1] * σx[2])) == (1, σx[1])
        @test act(cnot, mono(σx[1] * σz[2])) == (-1, mono(σy[1] * σy[2]))

        swap = CliffordSymmetry(:SWAP, 1, 2)
        @test act(swap, σx[1]) == (1, σx[2])
        @test act(swap, σy[2]) == (1, σy[1])
        @test act(swap, mono(σx[1] * σz[2])) == (1, mono(σz[1] * σx[2]))

        heisenberg = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
        @test NCTSSoS._act_polynomial(swap, heisenberg) == heisenberg

        err = _capture_exception(() -> CliffordSymmetry(:CNOT, 1, 1))
        @test err isa ArgumentError
        @test occursin("distinct sites", sprint(showerror, err))

        err = _capture_exception(() -> CliffordSymmetry(Dict(mono(σx[1] * σx[2]) => (1, σx[1]))))
        @test err isa ArgumentError
        @test occursin("single-Pauli-letter", sprint(showerror, err))

        err = _capture_exception(() -> CliffordSymmetry(Dict(σx[1] => (1, one(σx[1]))); nqubits=1))
        @test err isa ArgumentError
        @test occursin("non-identity", sprint(showerror, err))

        raw_images = Dict{NormalMonomial{PauliAlgebra,UInt8},Tuple{Int,NormalMonomial{PauliAlgebra,UInt8}}}(
            NormalMonomial{PauliAlgebra,UInt8}(UInt8[1]) =>
                (2, NormalMonomial{PauliAlgebra,UInt8}(UInt8[3])),
        )
        err = _capture_exception(() -> CliffordSymmetry{UInt8}(1, raw_images))
        @test err isa ArgumentError
        @test occursin("signs must be", sprint(showerror, err))

        bad_pauli_index = NormalMonomial{PauliAlgebra,Int}(Int[-1])
        good_pauli_index = NormalMonomial{PauliAlgebra,Int}(Int[1])
        err = _capture_exception(() -> CliffordSymmetry(Dict(bad_pauli_index => (1, good_pauli_index)); nqubits=1))
        @test err isa ArgumentError
        @test occursin("non-positive Pauli index", sprint(showerror, err))

        err = _capture_exception(() -> CliffordSymmetry(Dict(
            σx[1] => (1, mono(σx[1] * σx[2])),
            σz[1] => (1, mono(σz[1] * σz[2])),
        ); nqubits=2))
        @test err isa ArgumentError
        @test occursin("commutation/anticommutation", sprint(showerror, err))

        _, (τx, τy, τz) = create_pauli_variables(1:4)
        err = _capture_exception(() -> CliffordSymmetry(Dict(
            τx[1] => (1, mono(τx[2] * τx[3])),
            τy[1] => (1, mono(τy[2] * τx[3] * τx[4])),
            τz[1] => (1, mono(τz[2] * τx[4])),
        ); nqubits=4))
        @test err isa ArgumentError
        @test occursin("commutation/anticommutation", sprint(showerror, err))

        err = _capture_exception(() -> CliffordSymmetry(Dict(
            σy[1] => (1, σz[1]),
            σz[1] => (1, σy[1]),
        ); nqubits=1))
        @test err isa ArgumentError
        @test occursin("Pauli multiplication", sprint(showerror, err))
    end

    @testset "CliffordSymmetryGroup closure matches SymbolicWedderburn action convention" begin
        _, (σx, _, σz) = create_pauli_variables(1:2)
        mono(poly) = only(monomials(poly))

        h_group = CliffordSymmetryGroup(CliffordSymmetry(:H, 1); integer_type=UInt8)
        s_group = CliffordSymmetryGroup(CliffordSymmetry(:S, 1); integer_type=UInt8)
        hs_group = CliffordSymmetryGroup(CliffordSymmetry(:H, 1), CliffordSymmetry(:S, 1); integer_type=UInt8)
        cnot_group = CliffordSymmetryGroup(CliffordSymmetry(:CNOT, 1, 2); integer_type=UInt8)
        swap_group = CliffordSymmetryGroup(CliffordSymmetry(:SWAP, 1, 2); integer_type=UInt8)

        @test length(h_group) == 2
        @test length(s_group) == 4
        @test length(hs_group) == 24
        @test length(cnot_group) == 2
        @test length(swap_group) == 2
        @test NCTSSoS.GroupsCore.order(Int, only(NCTSSoS.GroupsCore.gens(h_group))) == 2

        action = NCTSSoS.NCPauliCliffordAction{UInt8}()
        probe = mono(σx[1] * σz[2])
        for group in (cnot_group, hs_group), g in collect(group), h in collect(group)
            lhs = SW.action(action, g * h, probe)
            mid_mono, mid_sign = SW.action(action, g, probe)
            rhs_mono, rhs_sign = SW.action(action, h, mid_mono)
            @test lhs == (rhs_mono, mid_sign * rhs_sign)
        end

        high_swap = CliffordSymmetry(:SWAP, 999, 1000; nqubits=1000, integer_type=UInt16)
        high_domain = UInt16[NCTSSoS._pauli_index(999, 0)]
        high_group = CliffordSymmetryGroup(high_swap; nqubits=1000, integer_type=UInt16, domain=high_domain)
        @test length(high_group) == 2
        @test length(high_group.inner.domain) == 6

        spec = SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2))
        @test length(spec.clifford_generators) == 1
    end

    @testset "Pauli Clifford symmetry wires into moment relaxation" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:2)
        heisenberg = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
        pop = polyopt(heisenberg, reg)
        basis = [one(σx[1]); σx; σy; σz]
        cfg = SolverConfig(
            optimizer=nothing,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2)),
        )
        sparsity = compute_sparsity(pop, cfg)

        @test NCTSSoS._check_symmetry_mvp_support(pop, cfg, sparsity) === nothing
        _, report = NCTSSoS.moment_relax_symmetric(
            pop,
            sparsity.corr_sparsity,
            sparsity.cliques_term_sparsities,
            cfg.symmetry,
        )

        @test report.group_order == 2
        @test report.psd_block_sizes == [4, 3]
        @test maximum(report.psd_block_sizes) < length(basis)
        @test all(==(:wedderburn), report.block_provenance)

        single_reg, (τx, _, τz) = create_pauli_variables(1:1)
        z_pop = polyopt(1.0 * τz[1], single_reg)
        z_basis = [one(τx[1]); τz]
        z_cfg = SolverConfig(
            optimizer=nothing,
            moment_basis=z_basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=SymmetrySpec(CliffordSymmetry(:S, 1)),
        )
        z_sparsity = compute_sparsity(z_pop, z_cfg)
        _, z_report = NCTSSoS.moment_relax_symmetric(
            z_pop,
            z_sparsity.corr_sparsity,
            z_sparsity.cliques_term_sparsities,
            z_cfg.symmetry,
        )
        @test z_report.group_order == 4
        @test z_report.psd_block_sizes == [2]
    end

    @testset "Pauli charge sectors and singlet constraints" begin
        n_ring = 4
        reg, (σx, σy, σz) = create_pauli_variables(1:n_ring)
        heisenberg = sum(
            ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n_ring)]
            for op in (σx, σy, σz) for i in 1:n_ring
        )
        pop = polyopt(heisenberg, reg)

        charge_spec = SymmetrySpec(pauli_charge=PauliChargeSectorSpec(nqubits=n_ring))
        charge_cfg = SolverConfig(
            optimizer=nothing,
            order=2,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=charge_spec,
        )
        charge_sparsity = compute_sparsity(pop, charge_cfg)
        _, charge_report = NCTSSoS.moment_relax_symmetric(
            pop,
            charge_sparsity.corr_sparsity,
            charge_sparsity.cliques_term_sparsities,
            charge_spec,
        )

        @test charge_report.psd_block_sizes == [6, 16, 23, 16, 6]
        @test [label.charge for label in charge_report.block_labels] == [-2, -1, 0, 1, 2]
        @test all(==(:charge_sector), charge_report.block_provenance)
        @test charge_report.basis_half_size == 67

        translation = pauli_site_permutation([2, 3, 4, 1])
        reflection = pauli_site_permutation([4, 3, 2, 1])
        spatial_singlet_spec = SymmetrySpec(
            [translation, reflection];
            pauli_charge=PauliChargeSectorSpec(nqubits=n_ring),
            pauli_singlet=PauliSingletConstraintSpec(nqubits=n_ring),
        )
        spatial_mp, spatial_report = NCTSSoS.moment_relax_symmetric(
            pop,
            charge_sparsity.corr_sparsity,
            charge_sparsity.cliques_term_sparsities,
            spatial_singlet_spec,
        )

        @test spatial_report.group_order == 8
        @test maximum(spatial_report.psd_block_sizes) < maximum(charge_report.psd_block_sizes)
        @test Set(label.charge for label in spatial_report.block_labels) == Set(-2:2)
        @test all(label.group_order == 8 for label in spatial_report.block_labels)
        @test all(==(:charge_wedderburn), spatial_report.block_provenance)
        @test count(con -> con[1] == :Zero, spatial_mp.constraints) > 0

        spatial_singlet_off_spec = SymmetrySpec(
            [translation, reflection];
            pauli_charge=PauliChargeSectorSpec(nqubits=n_ring),
            pauli_singlet=PauliSingletConstraintSpec(nqubits=n_ring),
            offblock_check=:off,
        )
        spatial_off_mp, spatial_off_report = NCTSSoS.moment_relax_symmetric(
            pop,
            charge_sparsity.corr_sparsity,
            charge_sparsity.cliques_term_sparsities,
            spatial_singlet_off_spec,
        )
        @test spatial_off_report.psd_block_sizes == spatial_report.psd_block_sizes
        @test spatial_off_report.block_labels == spatial_report.block_labels
        @test spatial_off_report.invariant_moment_count == spatial_report.invariant_moment_count
        @test length(spatial_off_mp.linear.moments) == length(spatial_mp.linear.moments)

        single_reg, (τx, τy, τz) = create_pauli_variables(1:1)
        bad_spec = SymmetrySpec(
            CliffordSymmetry(:H, 1);
            pauli_charge=PauliChargeSectorSpec(nqubits=1),
        )
        bad_pop = polyopt(1.0 * τz[1], single_reg)
        bad_cfg = SolverConfig(
            optimizer=nothing,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=bad_spec,
        )
        bad_sparsity = compute_sparsity(bad_pop, bad_cfg)
        err = _capture_exception(() -> NCTSSoS.moment_relax_symmetric(
            bad_pop,
            bad_sparsity.corr_sparsity,
            bad_sparsity.cliques_term_sparsities,
            bad_spec,
        ))
        @test err isa ArgumentError
        @test occursin("charge-compatible Clifford", sprint(showerror, err))

        @test_throws ArgumentError pauli_site_permutation(Int[])
        @test_throws ArgumentError pauli_site_permutation([1, 1])
        @test sprint(show, PauliChargeBlockLabel(0, nothing, 1, 3)) ==
            "PauliChargeBlockLabel(charge=0, finite_block=nothing, group_order=1, sector_dimension=3)"
        degree3_charge_spec = SymmetrySpec(
            pauli_charge=PauliChargeSectorSpec(nqubits=1, max_degree=3),
        )
        @test degree3_charge_spec.pauli_charge.max_degree == 3
        @test_throws ArgumentError SymmetrySpec(
            pauli_singlet=PauliSingletConstraintSpec(nqubits=1, max_degree=3),
        )
        @test_throws ArgumentError NCTSSoS._pauli_charge_letter_expansion(Int, 1, Int8(42))

        charge_half_basis = [one(τx[1]), τz[1], τx[1], τy[1]]
        identity_group = NCTSSoS.CliffordSymmetryGroup(pauli_site_permutation([1]))
        inferred_identity_groups = NCTSSoS._pauli_charge_transform_groups(
            charge_half_basis,
            PauliChargeSectorSpec(),
            identity_group,
        )
        inferred_identity_blocks = collect(Iterators.flatten(inferred_identity_groups))
        @test all(block -> block.label.group_order == 1, inferred_identity_blocks)
        @test all(block -> block.provenance == :charge_sector, inferred_identity_blocks)

        err = _capture_exception(() -> NCTSSoS._pauli_charge_transform_groups(
            [one(τz[1])],
            PauliChargeSectorSpec(),
            nothing,
        ))
        @test err isa ArgumentError
        @test occursin("needs `nqubits`", sprint(showerror, err))

        charge_non_neutral_spec = SymmetrySpec(
            pauli_charge=PauliChargeSectorSpec(nqubits=1),
        )
        charge_non_neutral_cfg = SolverConfig(
            optimizer=nothing,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=charge_non_neutral_spec,
        )
        charge_non_neutral_pop = polyopt(1.0 * τx[1], single_reg)
        charge_non_neutral_sparsity = compute_sparsity(charge_non_neutral_pop, charge_non_neutral_cfg)
        err = _capture_exception(() -> NCTSSoS.moment_relax_symmetric(
            charge_non_neutral_pop,
            charge_non_neutral_sparsity.corr_sparsity,
            charge_non_neutral_sparsity.cliques_term_sparsities,
            charge_non_neutral_spec,
        ))
        @test err isa ArgumentError
        @test occursin("charge-neutral", sprint(showerror, err))

        constrained_pop = polyopt(
            1.0 * τz[1],
            single_reg;
            eq_constraints=[1.0 * τz[1]],
            ineq_constraints=[1.0 * τz[1]],
            moment_eq_constraints=[1.0 * τz[1]],
        )
        @test NCTSSoS._check_pauli_charge_neutral(constrained_pop) === nothing

        bad_singlet_spec = SymmetrySpec(
            pauli_singlet=PauliSingletConstraintSpec(nqubits=0),
        )
        bad_singlet_cfg = SolverConfig(
            optimizer=nothing,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=bad_singlet_spec,
        )
        bad_singlet_sparsity = compute_sparsity(bad_pop, bad_singlet_cfg)
        err = _capture_exception(() -> NCTSSoS.moment_relax_symmetric(
            bad_pop,
            bad_singlet_sparsity.corr_sparsity,
            bad_singlet_sparsity.cliques_term_sparsities,
            bad_singlet_spec,
        ))
        @test err isa ArgumentError
        @test occursin("nqubits", sprint(showerror, err))

        nonpauli_reg, (u,) = create_unipotent_variables([("u", 1:1)])
        nonpauli_pop = polyopt(1.0 * u[1], nonpauli_reg)
        nonpauli_charge_cfg = SolverConfig(
            optimizer=nothing,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=charge_non_neutral_spec,
        )
        nonpauli_charge_sparsity = compute_sparsity(nonpauli_pop, nonpauli_charge_cfg)
        err = _capture_exception(() -> NCTSSoS.moment_relax_symmetric(
            nonpauli_pop,
            nonpauli_charge_sparsity.corr_sparsity,
            nonpauli_charge_sparsity.cliques_term_sparsities,
            charge_non_neutral_spec,
        ))
        @test err isa ArgumentError
        @test occursin("Pauli charge-sector splitting is only supported", sprint(showerror, err))

        nonpauli_singlet_spec = SymmetrySpec(
            pauli_singlet=PauliSingletConstraintSpec(nqubits=1),
        )
        nonpauli_singlet_cfg = SolverConfig(
            optimizer=nothing,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=nonpauli_singlet_spec,
        )
        nonpauli_singlet_sparsity = compute_sparsity(nonpauli_pop, nonpauli_singlet_cfg)
        err = _capture_exception(() -> NCTSSoS.moment_relax_symmetric(
            nonpauli_pop,
            nonpauli_singlet_sparsity.corr_sparsity,
            nonpauli_singlet_sparsity.cliques_term_sparsities,
            nonpauli_singlet_spec,
        ))
        @test err isa ArgumentError
        @test occursin("Pauli SU(2) singlet constraints are only supported", sprint(showerror, err))

        empty_clifford_charge = SymmetrySpec(
            CliffordSymmetry[];
            pauli_charge=PauliChargeSectorSpec(nqubits=1),
        )
        @test isempty(empty_clifford_charge.clifford_generators)
        @test empty_clifford_charge.pauli_charge isa PauliChargeSectorSpec

        overspecified_singlet = SymmetrySpec(
            pauli_singlet=PauliSingletConstraintSpec(nqubits=2),
        )
        overspecified_cfg = SolverConfig(
            optimizer=nothing,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=overspecified_singlet,
        )
        overspecified_sparsity = compute_sparsity(bad_pop, overspecified_cfg)
        err = _capture_exception(() -> NCTSSoS.moment_relax_symmetric(
            bad_pop,
            overspecified_sparsity.corr_sparsity,
            overspecified_sparsity.cliques_term_sparsities,
            overspecified_singlet,
        ))
        @test err isa ArgumentError
        @test occursin("supported by the relaxation basis", sprint(showerror, err))
    end

    @testset "SympleQ recognizes Pauli Hamiltonian Clifford symmetry" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:1)
        hamiltonian = 1.0 * (σx[1] + σz[1])
        spec = sympleq_symmetry_spec(hamiltonian)

        swap_idx = findfirst(spec.clifford_generators) do g
            NCTSSoS._act_polynomial(g, hamiltonian) == hamiltonian &&
                NCTSSoS._act_monomial(g, σx[1]) == (1, σz[1]) &&
                NCTSSoS._act_monomial(g, σz[1]) == (1, σx[1]) &&
                NCTSSoS._act_monomial(g, σy[1]) == (-1, σy[1])
        end
        @test !isnothing(swap_idx)
        generator = spec.clifford_generators[swap_idx]

        fixed_term = 1.0 * σx[1]
        fixed_spec = sympleq_symmetry_spec(fixed_term)
        @test any(fixed_spec.clifford_generators) do g
            NCTSSoS._act_polynomial(g, fixed_term) == fixed_term &&
                NCTSSoS._act_monomial(g, σx[1]) == (1, σx[1]) &&
                NCTSSoS._act_monomial(g, σz[1]) == (-1, σz[1])
        end

        signed_hamiltonian = 1.0 * (σx[1] - σz[1])
        signed_generators = sympleq_symmetry_spec(signed_hamiltonian).clifford_generators
        @test any(signed_generators) do g
            NCTSSoS._act_polynomial(g, signed_hamiltonian) == signed_hamiltonian &&
                NCTSSoS._act_monomial(g, σx[1]) == (-1, σz[1]) &&
                NCTSSoS._act_monomial(g, σz[1]) == (-1, σx[1])
        end

        pop = polyopt(hamiltonian, reg)
        basis = [one(σx[1]), σx[1], σz[1]]
        cfg = SolverConfig(
            optimizer=nothing,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=spec,
        )
        sparsity = compute_sparsity(pop, cfg)
        @test NCTSSoS._check_symmetry_mvp_support(pop, cfg, sparsity) === nothing
        _, report = NCTSSoS.moment_relax_symmetric(
            pop,
            sparsity.corr_sparsity,
            sparsity.cliques_term_sparsities,
            cfg.symmetry,
        )

        @test report.group_order == 2
        @test sort(report.psd_block_sizes) == [1, 2]

        _, (τx, τy, _) = create_pauli_variables(1:2)
        xy_pair = 1.0 * (τx[1] * τx[2] + τy[1] * τy[2])
        phase_spec = sympleq_symmetry_spec(xy_pair)
        @test any(phase_spec.clifford_generators) do g
            NCTSSoS._act_polynomial(g, xy_pair) == xy_pair &&
                NCTSSoS._act_monomial(g, τx[1]) == (1, τy[1]) &&
                NCTSSoS._act_monomial(g, τy[1]) == (-1, τx[1])
        end
    end

    @testset "SympleQ deterministic internals" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:3)
        hamiltonian = 2.0 * σx[1] + 2.0 * σz[1] + 3.0 * σx[2] * σz[3]
        tab = SymplecticTableau(hamiltonian)

        @test length(tab) == 3
        @test size(tab) == (3, 6)
        @test sprint(show, tab) == "SymplecticTableau(3 terms, 3 qubits)"
        @test tab.coeffs == [2.0 + 0.0im, 2.0 + 0.0im, 3.0 + 0.0im]
        @test tab.nqubits == 3

        y_row, y_eta = NCTSSoS._pauli_word_to_symplectic(σy[2], 3)
        @test y_eta == 1
        @test y_row == UInt8[0, 1, 0, 0, 1, 0]
        y_word, eta_back = NCTSSoS._symplectic_row_to_pauli_word(y_row, UInt8)
        @test eta_back == 1
        @test y_word == σy[2].word
        @test NCTSSoS._single_pauli_letter_from_symplectic(y_row, UInt8) == (only(σy[2].word), Int8(1))
        @test isnothing(NCTSSoS._single_pauli_letter_from_symplectic(UInt8[1, 1, 0, 0], UInt8))

        products = NCTSSoS.symplectic_product_matrix(tab)
        @test products == UInt8[0 1 0; 1 0 0; 0 0 0]
        graph = NCTSSoS.anticommutation_graph(tab)
        @test sprint(show, graph) == "SympleQGraph(3 vertices, 1 edges, 3 Pauli vertices, 0 cycle vertices)"
        @test graph.colors == [(:pauli, 2.0), (:pauli, 2.0), (:pauli, 3.0)]

        cycle_tab = SymplecticTableau(1.0 * (σx[1] + σx[2] + σx[1] * σx[2]))
        @test NCTSSoS.pauli_cycle_basis(cycle_tab) == [[1, 2, 3]]
        augmented = NCTSSoS.cycle_augmented_graph(cycle_tab)
        @test augmented.cycles == [[1, 2, 3]]
        @test augmented.auxiliary_vertices == [4]
        @test_throws ArgumentError NCTSSoS.pauli_cycle_basis(cycle_tab; cycle_strategy=:minimal_circuits)

        @test_throws ArgumentError NCTSSoS.TermPermutation([1, 1])
        swap_perm = NCTSSoS.TermPermutation([2, 1, 3])
        @test length(swap_perm) == 3
        @test swap_perm[1] == 2
        @test !isone(swap_perm)
        @test sprint(show, swap_perm) == "TermPermutation([2, 1, 3])"

        automorphisms = NCTSSoS.automorphism_generators(graph; backend=:backtracking)
        @test any(perm -> perm.images == [2, 1, 3], automorphisms)
        @test NCTSSoS.automorphism_generators(graph; backend=:bliss) == automorphisms
        @test_throws ArgumentError NCTSSoS.automorphism_generators(graph; backend=:nauty)

        S = NCTSSoS.symplectic_matrix_from_permutation(tab, swap_perm)
        @test sprint(show, S) == "SymplecticMatrix(6×6)"
        @test size(S) == (6, 6)
        @test Matrix(S) == S.data
        @test NCTSSoS.is_symplectic_matrix(S)
        @test NCTSSoS._gf2_matmul(tab.paulis, S.data) == tab.paulis[swap_perm.images, :]
        @test_throws ArgumentError NCTSSoS.symplectic_matrix_from_permutation(tab, NCTSSoS.TermPermutation([1, 2]))
        @test_throws ArgumentError NCTSSoS.SymplecticMatrix(ones(Int, 2, 3))
        @test_throws ArgumentError NCTSSoS.SymplecticMatrix(ones(Int, 3, 3))

        phase = PhaseVector(zeros(Int8, 6), true)
        @test length(phase) == 6
        @test phase[1] == 0
        generator = NCTSSoS.SympleQGenerator(swap_perm, S, phase, true)
        clifford = sympleq_clifford_symmetry(generator; integer_type=UInt8)
        @test NCTSSoS._act_polynomial(clifford, hamiltonian) == hamiltonian
        @test NCTSSoS._act_monomial(clifford, σx[1]) == (1, σz[1])
        @test NCTSSoS._act_monomial(clifford, σz[1]) == (1, σx[1])

        recovered_phase = NCTSSoS.recover_phase_vector(tab, swap_perm, S)
        @test recovered_phase.verified
        @test NCTSSoS.sympleq_generators(hamiltonian) isa Vector{NCTSSoS.SympleQGenerator}
        @test_throws ArgumentError sympleq_symmetry_spec(1.0 * one(σx[1]))
    end

    @testset "SympleQ Witt extension scales past 10 qubits" begin
        # GF(2) affine solver behind the Witt extension step.
        A = UInt8[1 0 1; 0 1 1]
        b = UInt8[1, 0]
        particular, kernel = NCTSSoS._gf2_affine_solution_space(A, b)
        @test mod.(Int.(A) * Int.(particular), 2) == Int.(b)
        @test length(kernel) == 1
        @test mod.(Int.(A) * Int.(only(kernel)), 2) == [0, 0]
        # Inconsistent system: x₁ = 0 and x₁ = 1.
        @test NCTSSoS._gf2_affine_solution_space(UInt8[1 0; 1 0], UInt8[0, 1]) === nothing

        # Rank-deficient target extension satisfies the pairing constraints and
        # stays independent of the target rows (dim 6, rank-2 partial basis).
        D = UInt8[1 0 0 0 0 0; 0 1 0 0 0 0]
        R = UInt8[0 1 0 0 0 0; 1 0 0 0 0 0]
        d = NCTSSoS._choose_domain_extension(D)
        r = NCTSSoS._choose_target_extension(d, D, R)
        @test !NCTSSoS._gf2_row_in_span(r, R)
        for i in axes(D, 1)
            @test NCTSSoS._gf2_pairing(d, view(D, i, :)) == NCTSSoS._gf2_pairing(r, view(R, i, :))
        end

        # 12-site Heisenberg ring: the binary symplectic dimension is 24 and the
        # tableau has GF(2) rank deficiency 2, so every detected automorphism
        # exercises the Witt extension. The pre-linear-solve implementation
        # enumerated all 2^dim GF(2) vectors and refused dim > 20 outright,
        # which silently collapsed the detected group on ≥ 11 qubits.
        n_ring = 12
        _, (ρx, ρy, ρz) = create_pauli_variables(1:n_ring)
        ring_h = sum(
            ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, n_ring)]
            for op in (ρx, ρy, ρz) for i in 1:n_ring
        )
        ring_spec = sympleq_symmetry_spec(ring_h)
        @test length(ring_spec.clifford_generators) >= 3
        @test all(
            g -> NCTSSoS._act_polynomial(g, ring_h) == ring_h,
            ring_spec.clifford_generators,
        )
    end

    @testset "symmetry helpers cover invariant constraints and scalar reductions" begin
        reg, (x,) = create_unipotent_variables([("x", 1:2)])
        invariant_poly = 1.0 * (x[1] + x[2])
        pop = polyopt(
            invariant_poly,
            reg;
            eq_constraints=[invariant_poly],
            ineq_constraints=[1.0 - invariant_poly],
            moment_eq_constraints=[invariant_poly],
        )
        cfg = SolverConfig(
            optimizer=SOLVER,
            order=1,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        sparsity = compute_sparsity(pop, cfg)

        T = typeof(x[1].word[1])
        spec = NCTSSoS._convert_symmetry_spec(
            T,
            SymmetrySpec(SignedPermutation(
                x[1].word[1] => x[2].word[1],
                x[2].word[1] => x[1].word[1],
            )),
        )
        domain = NCTSSoS._symmetry_domain(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        group = NCTSSoS._enumerate_symmetry_group(spec, domain)
        sw_group = NCTSSoS._sw_signed_permutation_group(spec, group, domain)

        @test domain == sort!(T[x[1].word[1], x[2].word[1]])
        @test NCTSSoS._check_symmetry_invariance(pop, group) === nothing

        @test NCTSSoS._assert_poly_equal(invariant_poly, invariant_poly, "same") === nothing
        err = _capture_exception(() -> NCTSSoS._assert_poly_equal(invariant_poly, invariant_poly + 1.0, "different"))
        @test err isa ArgumentError
        @test occursin("different", sprint(showerror, err))

        reducer = NCTSSoS._build_orbit_reducer([one(x[1])], group)
        scalar_blocks = NCTSSoS._reduce_constraint_matrix_symmetric(
            reshape([1.0 * one(x[1])], 1, 1),
            [one(x[1])],
            sw_group,
            reducer,
        )
        @test scalar_blocks == [reshape([1.0 * one(x[1])], 1, 1)]

        total_basis, moment_eq_row_bases, moment_eq_row_basis_degrees =
            NCTSSoS._polynomial_total_basis(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        reducer_full = NCTSSoS._build_orbit_reducer(total_basis, group)
        constraints = Tuple{Symbol, Matrix{Polynomial{UnipotentAlgebra,T,Float64}}}[]
        NCTSSoS._add_moment_eq_constraints_symmetric!(
            constraints,
            pop,
            moment_eq_row_bases,
            moment_eq_row_basis_degrees,
            reducer_full,
            Polynomial{UnipotentAlgebra,T,Float64},
        )

        @test length(constraints) == 3
        @test all(cone == :Zero for (cone, _) in constraints)
        @test sort([string(only(mat)) for (_, mat) in constraints]) == ["1 + x₁x₂", "1 + x₁x₂", "2.0 * x₁"]
    end

    @testset "lazy orbit reduction helper branches" begin
        reg, (x,) = create_unipotent_variables([("x", 1:2)])
        T = typeof(x[1].word[1])
        signed_swap = SignedPermutation(x[1].word[1] => (-1, x[2].word[1]))
        sign, image = NCTSSoS._orbit_reducer_action(signed_swap, x[1])
        @test sign == -1
        @test image == x[2]

        _, (σx, _, σz) = create_pauli_variables(1:2)
        pauli_T = typeof(σx[1].word[1])
        spatial_swap = pauli_site_permutation([2, 1])
        spatial_sign, spatial_image = NCTSSoS._orbit_reducer_action(spatial_swap, σx[1])
        @test spatial_sign == 1
        @test spatial_image == σx[2]

        hadamard = CliffordSymmetry(:H, 1; nqubits=2)
        h_sign, h_image = NCTSSoS._orbit_reducer_action(hadamard, σx[1])
        @test h_sign == 1
        @test h_image == σz[1]

        lazy_h = NCTSSoS._lazy_orbit_reducer(PauliAlgebra, pauli_T, [hadamard])
        @test lazy_h.spatial_permutations === nothing

        @test NCTSSoS._wedderburn_row_basis([1.0 0.0; 0.0 1.0]) == [1.0 0.0; 0.0 1.0]

        nonpauli_pop = polyopt(1.0 * x[1], reg)
        lazy_charge_spec = SymmetrySpec(
            pauli_charge=PauliChargeSectorSpec(nqubits=1),
            offblock_check=:off,
        )
        @test !NCTSSoS._can_use_lazy_pauli_charge_reduction(nonpauli_pop, lazy_charge_spec)

        P = Polynomial{UnipotentAlgebra,T,Float64}
        M = typeof(one(x[1]))
        multi_clique_cs = NCTSSoS.CorrelativeSparsity{UnipotentAlgebra,T,P,M,Nothing}(
            [T[x[1].word[1]], T[x[2].word[1]]],
            reg,
            P[],
            [Int[], Int[]],
            Int[],
            [[one(x[1]), x[2]], [x[1], x[2]]],
            [Vector{M}[], Vector{M}[]],
        )
        half_basis = NCTSSoS._half_basis_vector(multi_clique_cs)
        @test issorted(half_basis)
        @test Set(half_basis) == Set([one(x[1]), x[1], x[2]])
    end

    @testset "signed-permutation group multiplication matches SymbolicWedderburn's action convention" begin
        _, _, x, y = _create_chsh_problem()
        T = typeof(x[1].word[1])
        spec = NCTSSoS._convert_symmetry_spec(T, _create_chsh_symmetry(x, y))
        domain = sort!(T[x[1].word[1], x[2].word[1], y[1].word[1], y[2].word[1]])
        group = NCTSSoS._enumerate_symmetry_group(spec, domain)
        sw_group = NCTSSoS._sw_signed_permutation_group(spec, group, domain)
        action = NCTSSoS.NCWordSignedPermutationAction{UnipotentAlgebra,T}(UnipotentAlgebra)
        probe = only(monomials(x[1] * y[2]))

        for g in collect(sw_group), h in collect(sw_group)
            lhs = SW.action(action, g * h, probe)
            mid_mono, mid_sign = SW.action(action, g, probe)
            rhs_mono, rhs_sign = SW.action(action, h, mid_mono)
            @test lhs == (rhs_mono, mid_sign * rhs_sign)
        end
    end

    @testset "fermionic signed support keeps creation and annihilation distinct" begin
        pop, _, a, a_dag, basis, _ = _create_spinless_fermion_swap_problem()
        cfg = SolverConfig(
            optimizer=SOLVER,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        sparsity = compute_sparsity(pop, cfg)

        domain = NCTSSoS._symmetry_domain(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        mode_domain = NCTSSoS._mode_symmetry_domain(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        T = typeof(a[1].word[1])

        @test domain == sort!(T[-2, -1, 1, 2])
        @test mode_domain == [1, 2]
        @test Set(domain) != Set(T[1, 2])
        @test only(monomials(NCTSSoS._act_polynomial(FermionicModePermutation(1 => 2, 2 => 1), a_dag[1] * a[2]))) == only(monomials(a_dag[2] * a[1]))
    end

    @testset "signed-permutation keys are collision-free on signed domains" begin
        domain = [-1, 1]
        g1 = SignedPermutation(Dict(-1 => (1, -1), 1 => (1, 1)))
        g2 = SignedPermutation(Dict(-1 => (-1, 1), 1 => (-1, -1)))

        @test g1 != g2
        @test NCTSSoS._symmetry_key(g1, domain) != NCTSSoS._symmetry_key(g2, domain)
    end

    @testset "complex transformed blocks use unitary conjugation" begin
        reg, (σx, _, _) = create_pauli_variables(1:1)
        one_poly = 1.0 * one(σx[1])
        zero_poly = zero(one_poly)
        mat = [one_poly zero_poly; zero_poly 2.0 * one_poly]
        U = ComplexF64[1 im; 1 -im] / sqrt(2)

        transformed = NCTSSoS._transform_polynomial_block(mat, U, U)
        expected = U * ComplexF64[1 0; 0 2] * U'

        for j in 1:2, i in 1:2
            @test only(monomials(transformed[i, j])) == one(σx[1])
            @test only(coefficients(transformed[i, j])) ≈ expected[i, j] atol = 1e-12
        end
    end

    @testset "fermionic irrep labels compose exactly" begin
        reg, (a, a_dag) = create_fermionic_variables(1:3)
        c3_table = AbelianIrrepTable(0; multiply=(x, y) -> mod(x + y, 3), dual=x -> mod(-x, 3))
        layout = FermionicModeLayout(
            Dict(1 => 1, 2 => 2, 3 => 3);
            irrep_of=Dict(1 => 1, 2 => 2, 3 => 1),
            irrep_table=c3_table,
        )
        sector = FermionicSectorSpec(mode_layout=layout, split_irrep=true, split_parity=true)

        label_create_annihilate = NCTSSoS._fermionic_sector_label(
            only(monomials(a_dag[1] * a[2])),
            sector,
        )
        label_double_create = NCTSSoS._fermionic_sector_label(
            only(monomials(a_dag[3] * a_dag[1])),
            sector,
        )
        neutral_label = NCTSSoS._fermionic_sector_label(one(a[1]), sector)

        @test label_create_annihilate == FermionicSectorLabel(0, nothing, nothing, 2)
        @test label_double_create == FermionicSectorLabel(0, nothing, nothing, 2)
        @test neutral_label == FermionicSectorLabel(0, nothing, nothing, 0)

        missing_mode_layout = FermionicModeLayout(
            Dict(1 => 1);
            irrep_of=Dict(1 => 1),
            irrep_table=c3_table,
        )
        err = _capture_exception(() -> NCTSSoS._validate_active_fermionic_sector_metadata(
            FermionicSectorSpec(mode_layout=missing_mode_layout, split_irrep=true),
            [1, 2],
        ))
        @test err isa ArgumentError
        @test occursin("orbital_of[2]", sprint(showerror, err))
    end

    @testset "spin Casimir transform splits a spin-zero density sector deterministically" begin
        reg, ((c_up, c_up_dag), (c_dn, c_dn_dag)) =
            create_fermionic_variables([("c_up", 1:1), ("c_dn", 1:1)])
        up_mode = Int(c_up[1].word[1])
        dn_mode = Int(c_dn[1].word[1])
        layout = FermionicModeLayout(
            Dict(up_mode => 1, dn_mode => 1);
            spin2_of=Dict(up_mode => 1, dn_mode => -1),
        )
        sector = FermionicSectorSpec(
            mode_layout=layout,
            split_parity=true,
            split_number=true,
            split_spin=true,
        )
        spin = FermionicSpinAdaptationSpec(mode_layout=layout)
        density_basis = [
            only(monomials(c_up_dag[1] * c_up[1])),
            only(monomials(c_dn_dag[1] * c_dn[1])),
        ]
        label = NCTSSoS._fermionic_sector_label(first(density_basis), sector)

        blocks = NCTSSoS._fermionic_sector_transform_blocks(
            density_basis,
            label,
            nothing,
            spin,
            [up_mode, dn_mode],
        )

        @test length(blocks) == 2
        @test [size(block.row_basis, 1) for block in blocks] == [1, 1]
        @test [block.label.total_spin2 for block in blocks] == [0, 2]
        @test all(block.provenance == :sector_spin for block in blocks)
        @test blocks[1].row_basis ≈ ComplexF64[1 1] / sqrt(2) atol = 1e-12
        @test blocks[2].row_basis ≈ ComplexF64[1 -1] / sqrt(2) atol = 1e-12
    end

    @testset "two-orbital pair sector resolves into singlet and triplet blocks" begin
        oracle = expectations_oracle(
            "expectations/fermionic_symmetry.toml",
            "two_orbital_pair_sector_spin_blocks",
        )
        reg, ((c_up, c_up_dag), (c_dn, c_dn_dag)) =
            create_fermionic_variables([("c_up", 1:2), ("c_dn", 1:2)])
        up_modes = [Int(op.word[1]) for op in c_up]
        dn_modes = [Int(op.word[1]) for op in c_dn]
        layout = FermionicModeLayout(
            Dict(
                up_modes[1] => 1,
                up_modes[2] => 2,
                dn_modes[1] => 1,
                dn_modes[2] => 2,
            );
            spin2_of=Dict(
                up_modes[1] => 1,
                up_modes[2] => 1,
                dn_modes[1] => -1,
                dn_modes[2] => -1,
            ),
        )
        sector = FermionicSectorSpec(
            mode_layout=layout,
            split_parity=true,
            split_number=true,
            split_spin=true,
        )
        spin = FermionicSpinAdaptationSpec(mode_layout=layout)
        pair_basis = [
            only(monomials(c_up_dag[1] * c_dn_dag[1])),
            only(monomials(c_up_dag[2] * c_dn_dag[2])),
            only(monomials(c_up_dag[1] * c_dn_dag[2])),
            only(monomials(c_up_dag[2] * c_dn_dag[1])),
        ]
        sector_label = NCTSSoS._fermionic_sector_label(first(pair_basis), sector)

        blocks = NCTSSoS._fermionic_sector_transform_blocks(
            pair_basis,
            sector_label,
            nothing,
            spin,
            sort!(vcat(up_modes, dn_modes)),
        )
        sort!(blocks; by=block -> block.label.total_spin2)

        @test [size(block.row_basis, 1) for block in blocks] == oracle.sides
        @test [block.label.total_spin2 for block in blocks] == [0, 2]
        @test all(block.provenance == :sector_spin for block in blocks)

        singlet_block = only(filter(block -> block.label.total_spin2 == 0, blocks))
        triplet_block = only(filter(block -> block.label.total_spin2 == 2, blocks))
        singlet_projector = singlet_block.row_basis' * singlet_block.row_basis
        triplet_projector = triplet_block.row_basis' * triplet_block.row_basis

        triplet_vector = ComplexF64[0, 0, 1, -1] / sqrt(2)
        expected_triplet_projector = triplet_vector * triplet_vector'
        singlet_mixed = ComplexF64[0, 0, 1, 1] / sqrt(2)
        expected_singlet_projector = Diagonal(ComplexF64[1, 1, 0, 0]) + singlet_mixed * singlet_mixed'

        @test triplet_projector ≈ expected_triplet_projector atol = 1e-12
        @test singlet_projector ≈ expected_singlet_projector atol = 1e-12
    end

    @testset "fermionic symmetric moment relaxation keeps parity constraints and supports sector blocks" begin
        pop, _, a, _, basis, symmetry = _create_spinless_fermion_swap_problem()
        cfg = SolverConfig(
            optimizer=SOLVER,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        sparsity = compute_sparsity(pop, cfg)

        group_only = SymmetrySpec(FermionicModePermutation(1 => 2, 2 => 1))
        mp_group, report_group = NCTSSoS.moment_relax_symmetric(
            pop,
            sparsity.corr_sparsity,
            sparsity.cliques_term_sparsities,
            group_only,
        )
        mp_sector, report_sector = NCTSSoS.moment_relax_symmetric(
            pop,
            sparsity.corr_sparsity,
            sparsity.cliques_term_sparsities,
            symmetry,
        )

        @test count(c -> c[1] == :Zero, mp_group.constraints) > 0
        @test report_group.psd_block_sizes == [3, 2]
        @test report_sector.psd_block_sizes == [1, 1, 1, 1, 1]
        @test count(==(:sector_wedderburn), report_sector.block_provenance) == 4
        @test count(==(:sector_split), report_sector.block_provenance) == 1
        @test Set(report_sector.block_labels) == Set([
            FermionicSectorLabel(1, -1, nothing, nothing),
            FermionicSectorLabel(1, 1, nothing, nothing),
            FermionicSectorLabel(0, 0, nothing, nothing),
        ])
    end

    @testset "spin-adapted symmetry path adds zero constraints for non-scalar off-blocks" begin
        reg, ((c_up, c_up_dag), (c_dn, c_dn_dag)) =
            create_fermionic_variables([("c_up", 1:1), ("c_dn", 1:1)])
        up_mode = Int(c_up[1].word[1])
        dn_mode = Int(c_dn[1].word[1])
        layout = FermionicModeLayout(
            Dict(up_mode => 1, dn_mode => 1);
            spin2_of=Dict(up_mode => 1, dn_mode => -1),
        )
        sector = FermionicSectorSpec(
            mode_layout=layout,
            split_parity=true,
            split_number=true,
            split_spin=true,
        )
        spin = FermionicSpinAdaptationSpec(mode_layout=layout)
        basis = [
            one(c_up[1]),
            only(monomials(c_up_dag[1] * c_up[1])),
            only(monomials(c_dn_dag[1] * c_dn[1])),
        ]
        ham = -(c_up_dag[1] * c_up[1] + c_dn_dag[1] * c_dn[1])
        pop = polyopt(ham, reg)
        cfg = SolverConfig(
            optimizer=SOLVER,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
        )
        plain = cs_nctssos(pop, cfg)
        sector_only = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                moment_basis=basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
                symmetry=SymmetrySpec(sector=sector),
            ),
        )
        spin_res = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                moment_basis=basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
                symmetry=SymmetrySpec(sector=sector, spin_adaptation=spin),
            ),
        )

        @test spin_res.objective ≈ plain.objective atol = 1e-6
        @test sector_only.objective ≈ plain.objective atol = 1e-6
        @test !isnothing(spin_res.symmetry)
        @test spin_res.symmetry.psd_block_sizes == [2, 1]
        @test Set(label.total_spin2 for label in spin_res.symmetry.block_labels if label isa FermionicSpinBlockLabel) == Set([0, 2])
        @test maximum(spin_res.symmetry.psd_block_sizes) < only(flatten_sizes(plain.moment_matrix_sizes))
    end
end

@testset "Off-diagonal block certificate" begin
    # Shared fixture: 2-site Heisenberg, order-1 basis, SWAP Clifford symmetry.
    function _offblock_fixture()
        reg, (σx, σy, σz) = create_pauli_variables(1:2)
        heisenberg = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
        pop = polyopt(heisenberg, reg)
        basis = [one(σx[1]); σx; σy; σz]
        spec = SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2))
        cfg = SolverConfig(
            optimizer=nothing,
            moment_basis=basis,
            cs_algo=NoElimination(),
            ts_algo=NoElimination(),
            symmetry=spec,
        )
        sparsity = compute_sparsity(pop, cfg)
        return pop, spec, sparsity
    end

    @testset "all modes produce the same reduction on a valid spec" begin
        pop, _, sparsity = _offblock_fixture()
        reports = map((:full, :randomized, :off)) do mode
            spec = SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2); offblock_check=mode)
            mp, report = NCTSSoS.moment_relax_symmetric(
                pop,
                sparsity.corr_sparsity,
                sparsity.cliques_term_sparsities,
                spec,
            )
            (mp, report)
        end
        for (mp, report) in reports
            @test report.group_order == 2
            @test report.psd_block_sizes == [4, 3]
            @test mp.n_unique_moment_matrix_elements == reports[1][1].n_unique_moment_matrix_elements
        end
    end

    @testset "solved objectives agree across modes" begin
        pop, _, _ = _offblock_fixture()
        reg, (σx, σy, σz) = create_pauli_variables(1:2)
        heisenberg = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
        pop = polyopt(heisenberg, reg)
        basis = [one(σx[1]); σx; σy; σz]
        objectives = map((:full, :randomized, :off)) do mode
            cfg = SolverConfig(
                optimizer=SOLVER,
                moment_basis=basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
                symmetry=SymmetrySpec(CliffordSymmetry(:SWAP, 1, 2); offblock_check=mode),
            )
            cs_nctssos(pop, cfg).objective
        end
        for obj in objectives
            @test obj ≈ -0.75 atol = 1e-4
        end
    end

    @testset "corrupted symmetry-adapted basis is caught by :full and :randomized but not :off" begin
        pop, spec, sparsity = _offblock_fixture()
        corr = sparsity.corr_sparsity
        cts = sparsity.cliques_term_sparsities
        basis = only(corr.clq_mom_mtx_bases)
        T = eltype(basis).parameters[2]  # NormalMonomial{PauliAlgebra,T} index type

        support_domain = NCTSSoS._symmetry_domain(pop, corr, cts)
        nqubits = max(
            NCTSSoS._clifford_nqubits_from_domain(support_domain),
            NCTSSoS._clifford_max_generator_nqubits(spec.clifford_generators),
        )
        sw_group = CliffordSymmetryGroup(
            spec.clifford_generators;
            nqubits, integer_type=T, domain=support_domain,
        )
        group = NCTSSoS.CliffordSymmetry{T}[NCTSSoS._clifford_group_value(el) for el in sw_group]
        total_basis, _, _ = NCTSSoS._polynomial_total_basis(pop, corr, cts)
        reducer = NCTSSoS._build_orbit_reducer(total_basis, group)

        MP_C = NCTSSoS._moment_problem_coeff_type(PauliAlgebra, ComplexF64)
        MP_P = NCTSSoS.Polynomial{PauliAlgebra,T,MP_C}
        _, raw_mat = NCTSSoS._build_constraint_matrix(one(pop.objective), basis, :HPSD)
        mat = NCTSSoS._convert_polynomial_matrix(MP_P, raw_mat)

        blocks = NCTSSoS._sw_decompose_half_basis(basis, sw_group)
        row_bases = [Matrix(b) for b in blocks]
        @test length(row_bases) == 2

        # Sanity: the genuine decomposition passes every mode.
        for mode in (:full, :randomized, :off)
            result = NCTSSoS._reduce_transformed_blocks(
                mat, row_bases, reducer; offblock_check=mode,
            )
            @test length(result) == 2
        end

        # Corrupt the basis: swap a row between the two isotypic blocks. The
        # diagonal blocks then no longer carry the full PSD constraint and the
        # off-diagonal blocks are nonzero.
        corrupted = deepcopy(row_bases)
        tmp = copy(corrupted[1][end, :])
        corrupted[1][end, :] = corrupted[2][1, :]
        corrupted[2][1, :] = tmp

        @test_throws ArgumentError NCTSSoS._reduce_transformed_blocks(
            mat, corrupted, reducer; offblock_check=:full,
        )
        @test_throws ArgumentError NCTSSoS._reduce_transformed_blocks(
            mat, corrupted, reducer; offblock_check=:randomized,
        )
        # :off skips the certificate — this documents the risk of disabling it.
        @test length(NCTSSoS._reduce_transformed_blocks(
            mat, corrupted, reducer; offblock_check=:off,
        )) == 2
    end

    @testset "invalid offblock_check mode is rejected" begin
        @test_throws ArgumentError SymmetrySpec(
            CliffordSymmetry(:SWAP, 1, 2); offblock_check=:bogus,
        )
    end
end
