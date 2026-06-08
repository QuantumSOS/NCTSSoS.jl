# test/relaxations/symmetry.jl
# Tests: fail-fast behavior for dense symmetry MVP

using Test, NCTSSoS, LinearAlgebra

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
