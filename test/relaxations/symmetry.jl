# test/relaxations/symmetry.jl
# Tests: fail-fast behavior for dense symmetry MVP

using Test, NCTSSoS

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

    @testset "multiplicity-bearing decompositions error" begin
        reg, (x,) = create_unipotent_variables([("x", 1:2)])
        T = typeof(x[1].word[1])
        spec = NCTSSoS._convert_symmetry_spec(
            T,
            SymmetrySpec(
                SignedPermutation(
                    x[1].word[1] => x[2].word[1],
                    x[2].word[1] => x[1].word[1],
                ),
            ),
        )
        domain = sort!(T[x[1].word[1], x[2].word[1]])
        group = NCTSSoS._enumerate_symmetry_group(spec, domain)
        sw_group = NCTSSoS._sw_signed_permutation_group(spec, group, domain)
        basis = [one(x[1]), x[1], x[2]]

        err = _capture_exception(() -> NCTSSoS._sw_decompose_half_basis(basis, sw_group))

        @test err isa ArgumentError
        @test occursin("multiplicity-free", sprint(showerror, err))
    end
end
