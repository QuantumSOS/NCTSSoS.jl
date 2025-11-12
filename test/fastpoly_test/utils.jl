using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials:
    star, symmetric_canonicalize, cyclic_canonicalize, _comm!, simplify!

@testset "Utilities" begin
    @ncpolyvar x y z

    @testset "Symmetric Canonical Form" begin
        sa = SimplifyAlgorithm(
            comm_gps=[[x, y, z]],
            is_projective=false,
            is_unipotent=false,
        )
        mono1 = monomial([z, y, x], [1, 1, 2])

        mono1_sym = symmetric_canonicalize(mono1, sa)
        @test mono1_sym.vars == [x, y, z]
        @test mono1_sym.z == [2, 1, 1]

        mono2 = monomial([x, y, z], [2, 1, 1])
        mono2_sym = symmetric_canonicalize(mono2, sa)
        @test mono2_sym.vars == [x, y, z]
        @test mono2_sym.z == [2, 1, 1]

        poly1 = Polynomial(
            [0.1, 0.2, 0.3],
            [
                monomial([x, y, z], [2, 1, 1]),
                monomial([z, y, x], [1, 1, 2]),
                monomial([], []),
            ],
        )
        poly1_sym = canonicalize(poly1, sa)

        @test poly1_sym.coeffs ≈ [0.3, 0.3]
        @test poly1_sym.monos == [monomial([], []), monomial([x, y, z], [2, 1, 1])]

        n = 3
        @ncpolyvar a[1:n]

        sa = SimplifyAlgorithm(
            comm_gps=[a],
            is_projective=false,
            is_unipotent=false,
        )

        poly3 = Polynomial(
            [1, -1, -1, 3, -2, 2, -1, -1, 6, 9, 9, -54, 142],
            [
                monomial([a[1]], [2]),
                monomial([a[1], a[2]], [1, 1]),
                monomial([a[2], a[1]], [1, 1]),
                monomial([a[2]], [2]),
                monomial([a[1], a[2], a[1]], [1, 1, 1]),
                monomial([a[1], a[2], a[1]], [1, 2, 1]),
                monomial([a[2], a[3]], [1, 1]),
                monomial([a[3], a[2]], [1, 1]),
                monomial([a[3]], [2]),
                monomial([a[2], a[3]], [2, 1]),
                monomial([a[3], a[2]], [1, 2]),
                monomial([a[3], a[2], a[3]], [1, 1, 1]),
                monomial([a[3], a[2], a[3]], [1, 2, 1]),
            ],
        )

        supp = sort(
            map([
                [1, 2, 2, 1],
                [3, 2, 2, 3],
                [1, 2, 1],
                [2, 2, 3],
                [3, 2, 3],
                [1, 1],
                [1, 2],
                [2, 2],
                [2, 3],
                [3, 3],
            ]) do v
                monomial(a[v], ones(length(v)))
            end,
        )

        coe = [1, -2, 3, -2, 6, -2, 18, -54, 2, 142]

        poly3_sym = canonicalize(poly3, sa)

        @test poly3_sym.coeffs == coe
        @test poly3_sym.monos == supp
    end

    @testset "Cyclic Canonical Form" begin
        sa = SimplifyAlgorithm(
            comm_gps=[[x, y, z]],
            is_projective=false,
            is_unipotent=false,
        )
        @test cyclic_canonicalize(x^2 * y^2, sa) == x^2 * y^2
        @test cyclic_canonicalize(z * y^2, sa) == y^2 * z
        @test cyclic_canonicalize(z * z, sa) == z^2
        @test isone(cyclic_canonicalize(one(x), sa))

        n = 3
        @ncpolyvar a[1:n]
        sa = SimplifyAlgorithm(
            comm_gps=[a],
            is_projective=false,
            is_unipotent=false,
        )
        f = mapreduce(
            +,
            zip(
                [1, -1, -1, 3, -2, 2, -1, -1, 6, 9, 9, -54, 142, 5, 5, 5],
                [
                    monomial([a[1]], [2]),
                    monomial([a[1], a[2]], [1, 1]),
                    monomial([a[2], a[1]], [1, 1]),
                    monomial([a[2]], [2]),
                    monomial([a[1], a[2], a[1]], [1, 1, 1]),
                    monomial([a[1], a[2], a[1]], [1, 2, 1]),
                    monomial([a[2], a[3]], [1, 1]),
                    monomial([a[3], a[2]], [1, 1]),
                    monomial([a[3]], [2]),
                    monomial([a[2], a[3]], [2, 1]),
                    monomial([a[3], a[2]], [1, 2]),
                    monomial([a[3], a[2], a[3]], [1, 1, 1]),
                    monomial([a[3], a[2], a[3]], [1, 2, 1]),
                    monomial([a[1], a[2]], [2, 2]),
                    monomial([a[1], a[2], a[3]], [1, 1, 1]),
                    monomial([a[3], a[2], a[1]], [1, 1, 1]),
                ],
            ),
        ) do (coef, mono)
            coef * tr(mono)
        end

        f2 = mapreduce(+,
            zip([7, 142, -2, 18, -54, 1, -2, 3, -2, 6, 10],
                [
                    monomial([a[1], a[2], a[1]], [1, 2, 1]),
                    monomial([a[3], a[2], a[3]], [1, 2, 1]),
                    monomial([a[1], a[2], a[1]], [1, 1, 1]),
                    monomial([a[2], a[3]], [2, 1]),
                    monomial([a[3], a[2], a[3]], [1, 1, 1]),
                    monomial([a[1]], [2]),
                    monomial([a[1], a[2]], [1, 1]),
                    monomial([a[2]], [2]),
                    monomial([a[2], a[3]], [1, 1]),
                    monomial([a[3]], [2]),
                    monomial([a[1], a[2], a[3]], [1, 1, 1]),
                ])
        ) do (coef, mono)
            coef * tr(mono)
        end
        @test canonicalize(f, sa) == canonicalize(f2, sa)
    end

    @testset "_comm" begin
        @ncpolyvar a[1:3]
        @ncpolyvar b[1:3]

        mono = a[1]^2 * b[2]^2 * a[2] * b[1]^3
        comm_gp_dict = Dict(zip([a; b], [fill(1, 3); fill(2, 3)]))
        _comm!(mono, comm_gp_dict)
        @test mono == a[1]^2 * a[2] * b[2]^2 * b[1]^3
        mono = a[1]^3 * a[3]
        _comm!(mono, comm_gp_dict)
        @test mono == a[1]^3 * a[3]

        @ncpolyvar x[1:2] y[1:2]
        comm_gp_dict = Dict(zip([x; y], [fill(1, 2); fill(2, 2)]))
        mono = y[1] * x[1] * y[2]^2 * x[2]
        _comm!(mono, comm_gp_dict)
        @test mono == x[1] * x[2] * y[1] * y[2]^2
    end

    @testset "_projective" begin
        @ncpolyvar x y z
        mono = y * x^3 * y * z^3
        sa = SimplifyAlgorithm(
            comm_gps=[[x, y, z]],
            is_projective=true,
            is_unipotent=false,
        )
        @test simplify!(mono, sa) == y * x * y * z
    end

    @testset "_unipotent" begin
        @ncpolyvar x y z
        mono = z * x * y * z^2 * y * x * z
        sa = SimplifyAlgorithm(
            comm_gps=[[x, y, z]],
            is_projective=false,
            is_unipotent=true,
        )
        @test simplify!(mono, sa) == one(mono)
    end
end
