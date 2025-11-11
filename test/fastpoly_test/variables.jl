using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: Variable
using NCTSSoS.FastPolynomials:
    polyarrayvar, buildpolyvar, get_basis, monomials

@testset "Variable" begin
    @testset "Creation by Macros" begin
        my_vars = @ncpolyvar x y z
        @test my_vars isa Tuple
        @test x isa Variable
        @test y isa Variable
        @test x.name == :x
        @test y.name == :y

        @ncpolyvar x[1:10]
        @test x isa Vector{Variable}
        @test length(x) == 10
        @test x[1].name == Symbol("x₁")

        @ncpolyvar x[1:10, 1:5]
        @test x isa Matrix{Variable}
        @test size(x) == (10, 5)
        @test x[1,3].name == Symbol("x₁,₃")
    end

    @testset "String" begin
        @ncpolyvar xyd[1:10]
        @test xyd[2].name == Symbol("xyd₂")

        @ncpolyvar xyz
        @test xyz.name == :xyz
    end

    @testset "Complex Conversion" begin
        @ncpolyvar x_real
        x_complex = Variable(:x_real; iscomplex=true)
        @test x_complex isa Variable
        @test x_complex.name == :x_real
        @test x_complex.iscomplex == true
    end

    @testset "Utils: polyarrayvar" begin
        x = polyarrayvar(:x, 1:3, 1:2; iscomplex=false)
        @test x isa Matrix{Variable}
        @test size(x) == (3, 2)
    end

    @testset "Utils: buildpolyvar" begin
        @test_throws ErrorException buildpolyvar(Expr(:call, :x, +), false)
    end

    @testset "Hash" begin
        @ncpolyvar x y z
        @test hash(x) == hash(x)
        @test hash(x) != hash(y)
        @test hash(x) != hash(z)

        @ncpolyvar x[1:10]
        @test hash(x[1]) == hash(x[1])
        @test hash(x[1]) != hash(x[2])
    end

    @testset "^" begin
        @ncpolyvar x

        @test x^2 == monomial([x], [2])
        @test x^0 == monomial([], [])

        @test_throws AssertionError Base.:(^)(x, -1)
    end

    @testset "Get basis" begin
        @ncpolyvar x y z

        monomials_deg2 = monomials([x, y, z], Val(2))
        @test sort(monomials_deg2) == sort([
            monomial([x], [2]),
            monomial([y], [2]),
            monomial([z], [2]),
            monomial([x, y], [1, 1]),
            monomial([x, z], [1, 1]),
            monomial([y, z], [1, 1]),
            monomial([z, y], [1, 1]),
            monomial([z, x], [1, 1]),
            monomial([y, x], [1, 1]),
        ])

        nc_basis_deg2 = get_basis([x, y, z], 2)

        @test sort(nc_basis_deg2) == sort([
            one(x),
            monomial([x], [1]),
            monomial([y], [1]),
            monomial([z], [1]),
            x^2,
            y^2,
            z^2,
            x * y,
            x * z,
            y * z,
            z * x,
            z * y,
            y * x,
        ])
    end
end
