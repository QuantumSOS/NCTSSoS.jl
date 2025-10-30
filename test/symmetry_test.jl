using Test
using NCTSSoS
using NCTSSoS: NCVariablePermutation, normalform
using PermutationGroups

@testset "Symmetry Support" begin
    @testset "Group Actions on Monomials" begin
        @ncpolyvar x[1:3]

        # Test identity permutation
        action = NCVariablePermutation(vec(x))
        mon = x[1] * x[2]
        g = one(PermGroup([perm"(1,2)"]))
        @test SymbolicWedderburn.action(action, g, mon) == mon

        # Test simple transposition (1,2)
        g = perm"(1,2)"
        result = SymbolicWedderburn.action(action, g, mon)
        expected = x[2] * x[1]
        @test result == expected

        # Test on monomial with exponents
        mon2 = x[1]^2 * x[2]
        result2 = SymbolicWedderburn.action(action, g, mon2)
        expected2 = x[2]^2 * x[1]
        @test result2 == expected2

        # Test cyclic permutation (1,2,3)
        g_cycle = perm"(1,2,3)"
        mon3 = x[1] * x[2] * x[3]
        result3 = SymbolicWedderburn.action(action, g_cycle, mon3)
        expected3 = x[2] * x[3] * x[1]
        @test result3 == expected3

        # Test on constant monomial
        mon_const = one(Monomial)
        @test SymbolicWedderburn.action(action, g, mon_const) == mon_const
    end

    @testset "Normal Form Computation" begin
        @ncpolyvar x[1:3]

        # Test with S_3 group
        G = PermGroup([perm"(1,2)", perm"(1,2,3)"])
        action = NCVariablePermutation(vec(x))
        sa = SimplifyAlgorithm(comm_gps=[vec(x)], is_unipotent=false, is_projective=false)

        # All permutations of x[1]*x[2] should have the same normal form in commutative case
        mon1 = x[1] * x[2]
        mon2 = x[2] * x[1]
        nf1 = normalform(mon1, G, action, sa)
        nf2 = normalform(mon2, G, action, sa)
        @test nf1 == nf2  # Both should map to x[1]*x[2] (lexicographically smaller)

        # Test that normal form is idempotent
        nf_again = normalform(nf1, G, action, sa)
        @test nf1 == nf_again
    end
end
