# Note: FastPolynomials is loaded by setup.jl
using NCTSSoS.FastPolynomials:
    StateWord,
    NCStateWord,
    StatePolynomial,
    NCStatePolynomial,
    Arbitrary,
    MaxEntangled,
    monomials,
    coefficients,
    ς,
    tr

@testset "StatePolynomial{MaxEntangled}" begin
    @testset "Creation" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([3])

        sw1 = tr(m1)
        sw2 = tr(m2)

        sp = StatePolynomial([1.0, 2.0], [sw1, sw2])

        @test sp isa StatePolynomial{Float64,MaxEntangled}
        @test length(sp.state_words) == 2
    end

    @testset "Automatic Deduplication" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        sw = tr(m)

        # Same state word appears twice with different coefficients
        sp = StatePolynomial([1.0, 2.0], [sw, sw])

        @test length(sp.state_words) == 1
        @test sp.coeffs == [3.0]  # 1.0 + 2.0
    end

    @testset "Zero Coefficient Removal" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        sw1 = tr(m1)
        sw2 = tr(m2)

        # Zero coefficient should be removed
        sp = StatePolynomial([0.0, 1.0], [sw1, sw2])
        @test length(sp.state_words) == 1

        # Cancellation
        sp2 = StatePolynomial([2.0, -2.0], [sw1, sw1])
        @test iszero(sp2)
    end

    @testset "Degree" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])  # degree 2
        m2 = Monomial{NonCommutativeAlgebra}([3, 4, 5])  # degree 3

        sw1 = tr(m1)
        sw2 = tr(m2)

        sp = StatePolynomial([1.0, 1.0], [sw1, sw2])
        @test degree(sp) == 3  # max degree
    end

    @testset "Variables" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([3, 4])

        sw1 = tr(m1)
        sw2 = tr(m2)

        sp = StatePolynomial([1.0, 1.0], [sw1, sw2])
        vars = variables(sp)

        @test all(i in vars for i in 1:4)
    end

    @testset "Zero and One" begin
        T = StatePolynomial{Float64,MaxEntangled,NonCommutativeAlgebra,Int64}

        sp_zero = zero(T)
        @test iszero(sp_zero)

        sp_one = one(T)
        @test isone(sp_one)
    end

    @testset "Addition" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        sw1 = tr(m1)
        sw2 = tr(m2)

        sp1 = StatePolynomial([1.0], [sw1])
        sp2 = StatePolynomial([2.0], [sw2])

        sp_sum = sp1 + sp2
        @test length(sp_sum.state_words) == 2
    end

    @testset "Subtraction" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        sw = tr(m)

        sp1 = StatePolynomial([3.0], [sw])
        sp2 = StatePolynomial([1.0], [sw])

        sp_diff = sp1 - sp2
        @test sp_diff.coeffs == [2.0]
    end

    @testset "Scalar Multiplication" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        sw = tr(m)

        sp = StatePolynomial([2.0], [sw])

        sp_scaled = 3.0 * sp
        @test sp_scaled.coeffs == [6.0]

        sp_scaled_right = sp * 3.0
        @test sp_scaled_right.coeffs == [6.0]
    end

    @testset "Polynomial Multiplication" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        sw1 = tr(m1)
        sw2 = tr(m2)

        sp1 = StatePolynomial([1.0], [sw1])
        sp2 = StatePolynomial([2.0], [sw2])

        sp_prod = sp1 * sp2
        @test sp_prod.coeffs == [2.0]
        @test length(sp_prod.state_words[1].state_monos) == 2
    end

    @testset "Equality" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        sw = tr(m)

        sp1 = StatePolynomial([1.0], [sw])
        sp2 = StatePolynomial([1.0], [sw])
        sp3 = StatePolynomial([2.0], [sw])

        @test sp1 == sp2
        @test sp1 != sp3
    end
end

@testset "StatePolynomial{Arbitrary}" begin
    @testset "Creation with ς" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([3])

        sw1 = ς(m1)
        sw2 = ς(m2)

        sp = StatePolynomial([1.0, 2.0], [sw1, sw2])

        @test sp isa StatePolynomial{Float64,Arbitrary}
    end

    @testset "Arithmetic" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        # Create polynomials
        sp1 = StatePolynomial([1.0], [ς(m1)])
        sp2 = StatePolynomial([2.0], [ς(m2)])

        # Addition
        sp_sum = sp1 + sp2
        @test length(sp_sum.state_words) == 2

        # Multiplication
        sp_prod = sp1 * sp2
        @test sp_prod.coeffs == [2.0]
    end
end

@testset "NCStatePolynomial{MaxEntangled}" begin
    @testset "Creation" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}([3])

        sw = tr(m_state)
        ncsw = NCStateWord(sw, m_nc)

        ncsp = NCStatePolynomial([1.0], [ncsw])

        @test ncsp isa NCStatePolynomial{Float64,MaxEntangled}
    end

    @testset "Automatic Deduplication" begin
        m_state = Monomial{NonCommutativeAlgebra}([1])
        m_nc = Monomial{NonCommutativeAlgebra}([2])

        sw = tr(m_state)
        ncsw = NCStateWord(sw, m_nc)

        ncsp = NCStatePolynomial([1.0, 2.0], [ncsw, ncsw])
        @test length(ncsp.nc_state_words) == 1
        @test ncsp.coeffs == [3.0]
    end

    @testset "Degree" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])  # degree 2
        m_nc = Monomial{NonCommutativeAlgebra}([3])        # degree 1

        sw = tr(m_state)
        ncsw = NCStateWord(sw, m_nc)

        ncsp = NCStatePolynomial([1.0], [ncsw])
        @test degree(ncsp) == 3
    end

    @testset "Monomials Accessor" begin
        m_state = Monomial{NonCommutativeAlgebra}([1])
        m_nc = Monomial{NonCommutativeAlgebra}([2])

        sw = tr(m_state)
        ncsw = NCStateWord(sw, m_nc)

        ncsp = NCStatePolynomial([1.0], [ncsw])
        monos = monomials(ncsp)

        @test length(monos) == 1
        @test monos[1] == ncsw
    end

    @testset "Zero and One" begin
        T = NCStatePolynomial{Float64,MaxEntangled,NonCommutativeAlgebra,Int64}

        ncsp_zero = zero(T)
        @test iszero(ncsp_zero)

        ncsp_one = one(T)
        @test isone(ncsp_one)
    end

    @testset "Addition" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        ncsw1 = NCStateWord(tr(m1), m1)
        ncsw2 = NCStateWord(tr(m2), m2)

        ncsp1 = NCStatePolynomial([1.0], [ncsw1])
        ncsp2 = NCStatePolynomial([2.0], [ncsw2])

        ncsp_sum = ncsp1 + ncsp2
        @test length(ncsp_sum.nc_state_words) == 2
    end

    @testset "Scalar Multiplication" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        ncsw = NCStateWord(tr(m), m)

        ncsp = NCStatePolynomial([2.0], [ncsw])

        ncsp_scaled = 3.0 * ncsp
        @test ncsp_scaled.coeffs == [6.0]
    end

    @testset "Equality" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        ncsw = NCStateWord(tr(m), m)

        ncsp1 = NCStatePolynomial([1.0], [ncsw])
        ncsp2 = NCStatePolynomial([1.0], [ncsw])
        ncsp3 = NCStatePolynomial([2.0], [ncsw])

        @test ncsp1 == ncsp2
        @test ncsp1 != ncsp3
    end
end

@testset "NCStatePolynomial{Arbitrary}" begin
    @testset "Creation" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}([3])

        sw = ς(m_state)
        ncsw = NCStateWord(sw, m_nc)

        ncsp = NCStatePolynomial([1.0], [ncsw])

        @test ncsp isa NCStatePolynomial{Float64,Arbitrary}
    end

    @testset "Arithmetic" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        ncsw1 = NCStateWord(ς(m1), m1)
        ncsw2 = NCStateWord(ς(m2), m2)

        ncsp1 = NCStatePolynomial([1.0], [ncsw1])
        ncsp2 = NCStatePolynomial([2.0], [ncsw2])

        # Addition
        ncsp_sum = ncsp1 + ncsp2
        @test length(ncsp_sum.nc_state_words) == 2

        # Subtraction
        ncsp_diff = ncsp1 - ncsp2
        @test length(ncsp_diff.nc_state_words) == 2
    end

    @testset "Variables" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}([3, 4])

        ncsw = NCStateWord(ς(m_state), m_nc)
        ncsp = NCStatePolynomial([1.0], [ncsw])

        vars = variables(ncsp)
        # NCStatePolynomial.variables() returns Set{T} of indices (new API)
        # Check that we have 4 variables with indices 1-4
        @test length(vars) == 4
        @test all(v in 1:4 for v in vars)
    end
end
