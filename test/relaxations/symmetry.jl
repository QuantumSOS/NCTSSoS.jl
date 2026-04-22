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
