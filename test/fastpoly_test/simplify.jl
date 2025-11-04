using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: simplify, get_state_basis, NCStateWord, get_basis
using NCTSSoS.FastPolynomials:  symmetric_canonicalize, Arbitrary, is_symmetric

@testset "Simplification Interface" begin
    @ncpolyvar x[1:3] y[1:3]

    sa1 = SimplifyAlgorithm(; comm_gps=[x, y], is_unipotent=true, is_projective=false)
    sa2 = SimplifyAlgorithm(; comm_gps=[x, y], is_unipotent=false, is_projective=true)

    @testset "Simplify Monomial" begin
        @test simplify(x[1] * y[1] * x[1] * y[2], sa1) == y[1] * y[2]
        @test simplify(x[1] * y[1] * x[1] * y[2], sa2) == x[1] * y[1] * y[2]
    end

	@testset "Simplify StateWord" begin
        sw = ς(x[1]^2) * ς(y[1] * y[2]) * ς(x[1] * y[1] * x[1])
        @test simplify(sw, sa1) == ς(y[1] * y[2]) * ς(y[1])
        @test simplify(sw, sa2) == ς(x[1])* ς(y[1] * y[2]) * ς(x[1] *y[1])
	end

	@testset "Simplify NCStateWord" begin
        ncsw =
            ς(x[1]^2) * ς(y[1] * y[2]) * ς(x[1] * y[1] * x[1]) * (x[1]^2 * y[1] * x[2] * y[1])

		@test simplify(ncsw, sa1) == 
            ς(y[1] * y[2]) * ς(y[1]) * monomial(x[2])

		@test simplify(ncsw, sa2) == 
            ς(x[1]) * ς(y[1] * y[2]) * ς(x[1] * y[1]) * (x[1] * x[2] * y[1])
	end

    @testset "Get Basis" begin
		@ncpolyvar a b

        sa1 = SimplifyAlgorithm(;
            comm_gps=[[a, b]], is_unipotent=true, is_projective=false
        )

        @test get_basis(Polynomial{ComplexF64}, [a, b], 2, sa1) == [
            one(Monomial),
            a,
            b,
            a * b,
            b * a,
        ]

		sa2 = SimplifyAlgorithm(;
			comm_gps=[[a, b]], is_unipotent=false, is_projective=true
		)

        @test get_basis(Polynomial{ComplexF64}, [a, b], 2, sa2) == [
            one(Monomial),
            a,
            b,
            a * b, 
            b * a
        ]


        sa3 = SimplifyAlgorithm(;
            comm_gps=[[a], [b]], is_unipotent=false, is_projective=false
        )

        @test get_basis(Polynomial{ComplexF64}, [a, b], 2, sa3) == [
            one(Monomial),
            a,
            b,
            a*b,
            a^2,
            b^2
        ]

        sa4 = SimplifyAlgorithm(;
            comm_gps=[[a, b]], is_unipotent=false, is_projective=false
        )

        @test get_basis(Polynomial{ComplexF64}, [a, b], 2, sa4) == [
            one(Monomial),
            a,
            b,
            a * b, 
            a^2,
            b * a,
            b^2
        ]
    end

	@testset "Get State Basis" begin
		@ncpolyvar x y
        sa1 = SimplifyAlgorithm(;
            comm_gps=[[x], [y]], is_unipotent=true, is_projective=false
        )
		sa2 = SimplifyAlgorithm(;
			comm_gps=[[x], [y]], is_unipotent=false, is_projective=true
		)

		sa3 = SimplifyAlgorithm(;
			comm_gps=[[x], [y]], is_unipotent=false, is_projective=false
		)

        target_sbasis_1 = [
            one(NCStateWord{Arbitrary}),
            ς(x) * one(Monomial),
            ς(y) * one(Monomial),
            ς(x * y) * one(Monomial),
            ς(x) * ς(x) * one(Monomial),
            ς(x) * ς(y) * one(Monomial),
            ς(y) * ς(y) * one(Monomial),
            ς(one(Monomial)) * monomial(x),
            ς(x) * monomial(x),
            ς(y) * monomial(x),
            ς(one(Monomial)) * monomial(y),
            ς(x) * monomial(y),
            ς(y) * monomial(y),
            ς(one(Monomial)) * (x * y),
        ]

        @test sort(get_state_basis(Arbitrary,[x, y], 2, sa1)) == sort(target_sbasis_1)

		target_sbasis_2 = [
            one(NCStateWord{Arbitrary}),
            ς(x) * one(Monomial),
            ς(y) * one(Monomial),
            ς(x * y) * one(Monomial),
            ς(one(Monomial)) * monomial(x),
            ς(one(Monomial)) * monomial(y),
            ς(x) * ς(x) * one(Monomial),
            ς(x) * ς(y) * one(Monomial),
            ς(y) * ς(y) * one(Monomial),
            ς(x) * monomial(x),
            ς(y) * monomial(x),
            ς(x) * monomial(y),
            ς(y) * monomial(y),
            ς(one(Monomial)) * (x * y),
        ]

        @test sort(get_state_basis(Arbitrary,[x, y], 2, sa2)) == sort(target_sbasis_2)

		target_sbasis_3 = [
            one(NCStateWord{Arbitrary}),
            ς(x) * one(Monomial),
            ς(y) * one(Monomial),
            ς(one(Monomial)) * monomial(x),
            ς(one(Monomial)) * monomial(y),
            ς(x) * ς(x) * one(Monomial),
            ς(x) * ς(y) * one(Monomial),
            ς(y) * ς(y) * one(Monomial),
            ς(x * y) * one(Monomial),
            ς(x^2) * one(Monomial),
            ς(y^2) * one(Monomial),
            ς(x) * monomial(x),
            ς(y) * monomial(x),
            ς(x) * monomial(y),
            ς(y) * monomial(y),
            ς(one(Monomial)) * (x * y),
            ς(one(Monomial)) * (x^2 ),
            ς(one(Monomial)) * (y^2),
        ]

        @test sort(get_state_basis(Arbitrary,[x, y], 2, sa3)) == sort(target_sbasis_3)
	end
end

@testset "Symmetric Canonicalie" begin
    @ncpolyvar x[1:2] y[1:2] 

    sa1 = SimplifyAlgorithm(; comm_gps=[x, y], is_unipotent=false, is_projective=false)

    sa2 = SimplifyAlgorithm(; comm_gps=[x, y], is_unipotent=true, is_projective=false)
    
    sa3 = SimplifyAlgorithm(; comm_gps=[x, y], is_unipotent=false, is_projective=true)

    @testset "Monomial" begin
        @test symmetric_canonicalize(x[2] * y[2]^2 * x[1] * y[1], sa1) ==
            x[1] * x[2] * y[1] * y[2]^2

        @test symmetric_canonicalize(x[2] * y[2]^2 * x[1] * y[1], sa2) ==
            x[1] * x[2] * y[1]

        @test symmetric_canonicalize(x[2] * y[2]^2 * x[1] * y[1], sa3) ==
            x[1] * x[2] * y[1] * y[2]
    end

    @testset "StateWord" begin
        @test canonicalize(
            ς(x[2] * y[1] * x[1]) * ς(y[2] * x[2] * x[1] * y[2]), sa1
        ) == ς(x[1] * x[2] * y[1]) * ς(x[1] * x[2] * y[2]^2)

        @test canonicalize(
            ς(x[2] * y[1] * x[1]) * ς(y[2] * x[2] * x[1] * y[2]), sa2
        ) == ς(x[1] * x[2] * y[1]) * ς(x[1] * x[2])

        @test canonicalize(
            ς(x[2] * y[1] * x[1]) * ς(y[2] * x[2] * x[1] * y[2]), sa3
        ) == ς(x[1] * x[2] * y[1]) * ς(x[1] * x[2]*y[2])
    end

    @testset "NCStateWord" begin
        @test canonicalize(
            ς(x[2] * y[1] * x[1]) *
            ς(y[2] * x[2] * x[1] * y[2]) *
            (x[2] * y[1] * x[1] * y[1]),
            sa1,
        ) == ς(x[1] * x[2] * y[1]) * ς(x[1] * x[2] * y[2]^2) * (x[1]*x[2]*y[1]^2)

        @test canonicalize(
            ς(x[2] * y[1] * x[1]) *
            ς(y[2] * x[2] * x[1] * y[2]) *
            (x[2] * y[1] * x[1] * y[1]),
            sa2,
        ) == ς(x[1] * x[2] * y[1]) * ς(x[1] * x[2]) * (x[1]*x[2])

        @test canonicalize(
            ς(x[2] * y[1] * x[1]) *
            ς(y[2] * x[2] * x[1] * y[2]) *
            (x[2] * y[1] * x[1] * y[1]),
            sa3,
        ) == ς(x[1] * x[2] * y[1]) * ς(x[1] * x[2] * y[2]) * (x[1]*x[2]*y[1])
    end
end

@testset "Test if polynomial is symmetric" begin
    @ncpolyvar x[1:2] y[1:2] z[1:2]
    sa = SimplifyAlgorithm(comm_gps=[[x[i],y[i],z[i]] for i in 1:2])

    sym_poly = sum(one(ComplexF64)*op[1]*op[2] for op in [x,y,z]) 

    @test is_symmetric(sym_poly, sa)
    unsym_poly = one(ComplexF64)*x[1]*y[1]- im* z[1]
    @test !is_symmetric(unsym_poly, sa)
end