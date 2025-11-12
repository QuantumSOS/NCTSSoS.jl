using Test, NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials:
    StateWord,
    NCStateWord,
    degree,
    NCStatePolynomial,
    StatePolynomial,
    get_state_basis,
    neat_dot,
    expval,
    monomials,
    Arbitrary,
    MaxEntangled

@testset "NCTraceStatePolynomial Components" begin
    @ncpolyvar x[1:2] y[1:2]

    @testset "NCTraceStatePolynomial" begin
        sws =
            NCStateWord.(
                Ref(MaxEntangled),
                [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3]],
                Ref(one(x[1])),
            )
        sp = ncstatepoly([1.0, 2.0, 5.0], sws)
        @test string(sp) ==
              "2.0 * tr(x₁¹x₂¹) * 1 + 5.0 * tr(x₂³) * 1 + 1.0 * tr(x₁¹x₂¹) * tr(x₂²) * 1"
        sws_rep =
            NCStateWord.(
                Ref(MaxEntangled),
                [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^2, x[1] * x[2]], [x[2]^3]],
                Ref(one(x[1])),
            )
        sp_rep = ncstatepoly([0.5, 2.0, 0.5, 5.0], sws_rep)
        @test sp == sp_rep

        sws_diff =
            NCStateWord.(
                Ref(MaxEntangled),
                [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3, x[1]]],
                Ref(one(x[1])),
            )
        sp_diff = ncstatepoly([1.0, 2.0, 5.0], sws_diff)
        @test sp_diff != sp

        @test degree(sp) == 4

        @test one(sp) == ncstatepoly([1.0], [one(sws[1])])

        @test sort(monomials(sp)) == sort(
            NCStateWord.(
                Ref(MaxEntangled),
                [[x[2]^2, x[1] * x[2]], [x[1] * x[2]], [x[2]^3]],
                Ref(one(x[1])),
            ),
        )

        sp_ez =
            -1.0 * tr(x[1] * y[1]) * one(x[1]) - 1.0 * tr(x[1] * y[2]) * one(x[1]) -
            (1.0 * tr(x[2] * y[1]) * one(x[1])) + 1.0 * tr(x[2] * y[2]) * one(x[1])

        sp_hard = mapreduce(
            a -> a[1] * NCStateWord(MaxEntangled, [a[2]], a[3]),
            +,
            zip(
                [-1.0, -1.0, -1.0, 1.0],
                [x[1] * y[1], x[1] * y[2], x[2] * y[1], x[2] * y[2]],
                fill(one(x[1]), 4),
            ),
        )
        @test sp_ez == sp_hard
    end
end

@testset "Trace State Polynomial" begin
    @testset "free trace state polynomial" begin
        @ncpolyvar x[1:2]
        ncsws = [
            prod(tr.(wd)) for
            wd in [[x[1] * x[2], x[2]^2], [x[2] * x[1], x[1]^2], [x[1] * x[2]], [x[2]^3]]
        ]
        spop = StatePolynomial([1.0, 0.5, 2.0, 3.0], ncsws)
        @test string(spop) == "2.0 * tr(x₁¹x₂¹) + 3.0 * tr(x₂³) + 0.5 * tr(x₁²) * tr(x₂¹x₁¹) + 1.0 * tr(x₁¹x₂¹) * tr(x₂²)"

        p1 = 1.0 + tr(x[1] * x[2])
        @test string(p1) == "1.0 * tr(1) + 1.0 * tr(x₁¹x₂¹)"
    end

    @testset "nc trace state polynomial" begin
        @ncpolyvar x[1:2]
        ncsws = [
            NCStateWord(MaxEntangled, wd, ncw) for (wd, ncw) in
            zip([[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3]], [x[1], x[2]^2, one(x[1])])
        ]

        spop = ncstatepoly([1.0, 2.0, 3.0], ncsws)

        @test string(spop) ==
              "3.0 * tr(x₂³) * 1 + 2.0 * tr(x₁¹x₂¹) * x₂² + 1.0 * tr(x₁¹x₂¹) * tr(x₂²) * x₁¹"

        @test degree(spop) == 5
        @test sort(monomials(spop)) == sort(
            map(
                a -> NCStateWord(MaxEntangled, a[1], monomial(a[2])),
                [
                    ([x[2]^3], one(x[1])),
                    ([x[2]^2, x[1] * x[2]], x[1]),
                    ([x[1] * x[2]], x[2]^2),
                ],
            ),
        )
    end

    @testset "Trace State Polynomial Arithmetic" begin
        @ncpolyvar x[1:2]

        sps = [
            1.0 * tr(x[1] * x[2]) + 2.0 * tr(x[2] * x[1]),
            2.0 * tr(x[2]^2) + 3.0 * tr(x[1]^2),
        ]

        sp_prod =
            2.0 * tr(x[1] * x[2]) * tr(x[2]^2) +
            4.0 * tr(x[2] * x[1]) * tr(x[2]^2) +
            3.0 * tr(x[1] * x[2]) * tr(x[1]^2) +
            6.0 * tr(x[2] * x[1]) * tr(x[1]^2)
        @test sps[1] * sps[2] == sp_prod
        sp_sum =
            1.0 * tr(x[1] * x[2]) +
            2.0 * tr(x[2] * x[1]) +
            2.0 * tr(x[2]^2) +
            3.0 * tr(x[1]^2)
        @test sps[1] + sps[2] == sp_sum
    end

    @testset "variables of NC Trace State Polynomial" begin
        @ncpolyvar x[1:10]

        sp =
            1.0 * (tr(x[4] * x[5]) * tr(x[6] * x[7])) +
            2.0 * tr(x[8]) +
            3.0 * tr(x[9] * x[10])
        @test sort(variables(sp)) == sort(x[4:10])

        words = monomial.([x[1], x[2], x[3]])
        spop =
            1.0 * NCStateWord(MaxEntangled, [x[4] * x[5], x[6] * x[7]], x[1]) +
            2.0 * NCStateWord(MaxEntangled, [x[8]], x[2]) +
            3.0 * NCStateWord(MaxEntangled, [x[9], x[10]], x[3])
        @test variables(spop) == sort(x)
    end

    @testset "Trace State Polynomial Op comparison" begin
        @ncpolyvar a b c
        words_order1 = [a * b, b * c, c^2]
        words_order2 = [b * c, c^2, a * b]
        state_word_order1 =
            NCStateWord.(Ref(MaxEntangled), [[a^2], [b^2], [c^2]], Ref(one(a)))
        state_word_order2 =
            NCStateWord.(Ref(MaxEntangled), [[b^2], [c^2], [a^2]], Ref(one(a)))
        state_poly_order1 = ncstatepoly([1.0, 2.0, 3.0], state_word_order1)
        state_poly_order2 = ncstatepoly([2.0, 3.0, 1.0], state_word_order2)
        @test state_poly_order1 == state_poly_order2
    end
end

@testset "get state basis" begin
    @ncpolyvar x y

    sa = SimplifyAlgorithm(; comm_gps=[[x, y]], is_projective=false, is_unipotent=false)
    c_words = [[one(x)], [y], [x], [one(x)], [one(x)]]
    nc_words = [one(x), one(x), one(x), y, x]
    @test sort(get_state_basis(MaxEntangled, [x, y], 1, sa)) ==
          sort(map(a -> NCStateWord(MaxEntangled, a[1], a[2]), zip(c_words, nc_words)))

    c_words = [
        [one(x)],
        [y],
        [x],
        [x * y],
        [y^2],
        [x^2],
        [y, y],
        [y, x],
        [x, x],
        [one(x)],
        [y],
        [x],
        [one(x)],
        [y],
        [x],
        fill([one(x)], 4)...,
    ]
    nc_words = [fill(one(x), 9); fill(y, 3); fill(x, 3); [y * x, y^2, x * y, x^2]]

    @test sort(get_state_basis(MaxEntangled, [x, y], 2, sa)) ==
          sort(map(a -> NCStateWord(MaxEntangled, a[1], a[2]), zip(c_words, nc_words)))

    nc_words = [fill(one(x), 7); fill(x, 4); fill(x^2, 2); [x^3]]
    c_words = [
        [one(x)],
        [x],
        [x^2],
        [x^3],
        [x, x],
        [x, x^2],
        [x, x, x],
        [one(x)],
        [x],
        [x^2],
        [x, x],
        [one(x)],
        [x],
        [one(x)],
    ]
    @test sort(get_state_basis(MaxEntangled, [x], 3, sa)) ==
          sort(map(a -> NCStateWord(MaxEntangled, a[1], a[2]), zip(c_words, nc_words)))
end

@testset "NCStatePolynomial Components" begin
    @ncpolyvar x[1:2] y[1:2]

    @testset "NCStatePolynomial" begin
        sws =
            NCStateWord.(
                Ref(Arbitrary),
                [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3]],
                Ref(one(x[1])),
            )
        sp = ncstatepoly([1.0, 2.0, 5.0], sws)
        @test string(sp) ==
              "2.0 * <x₁¹x₂¹> * 1 + 5.0 * <x₂³> * 1 + 1.0 * <x₁¹x₂¹> * <x₂²> * 1"
        sws_rep =
            NCStateWord.(
                Ref(Arbitrary),
                [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^2, x[1] * x[2]], [x[2]^3]],
                Ref(one(x[1])),
            )
        sp_rep = ncstatepoly([0.5, 2.0, 0.5, 5.0], sws_rep)
        @test sp == sp_rep

        sws_diff =
            NCStateWord.(
                Ref(Arbitrary),
                [[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3, x[1]]],
                Ref(one(x[1])),
            )
        sp_diff = ncstatepoly([1.0, 2.0, 5.0], sws_diff)
        @test sp_diff != sp

        @test degree(sp) == 4

        @test one(sp) == ncstatepoly([1.0], [one(sws[1])])

        @test one(NCStatePolynomial{Float64,Arbitrary}) == ncstatepoly([1.0], [one(sws[1])])

        @test zero(sp) == ncstatepoly([0.0], [one(sws[1])])
        @test zero(NCStatePolynomial{Float64,Arbitrary}) == ncstatepoly([0.0], [one(sws[1])])

        @test sort(monomials(sp)) == sort(
            NCStateWord.(
                Ref(Arbitrary),
                [[x[2]^2, x[1] * x[2]], [x[1] * x[2]], [x[2]^3]],
                Ref(one(x[1])),
            ),
        )

        sp_ez =
            -1.0 * ς(x[1] * y[1]) * one(x[1]) - 1.0 * ς(x[1] * y[2]) * one(x[1]) -
            (1.0 * ς(x[2] * y[1]) * one(x[1])) + 1.0 * ς(x[2] * y[2]) * one(x[1])
        sp_hard = mapreduce(
            a -> a[1] * NCStateWord(Arbitrary, [a[2]], a[3]),
            +,
            zip(
                [-1.0, -1.0, -1.0, 1.0],
                [x[1] * y[1], x[1] * y[2], x[2] * y[1], x[2] * y[2]],
                fill(one(x[1]), 4),
            ),
        )
        @test sp_ez == sp_hard
    end
end

@testset "State Polynomial" begin
    @testset "free state polynomial" begin
        @ncpolyvar x[1:2]
        ncsws = [
            prod(ς.(wd)) for
            wd in [[x[1] * x[2], x[2]^2], [x[2] * x[1], x[1]^2], [x[1] * x[2]], [x[2]^3]]
        ]
        spop = StatePolynomial([1.0, 0.5, 2.0, 3.0], ncsws)
        @test string(spop) == "2.0 * <x₁¹x₂¹> + 3.0 * <x₂³> + 0.5 * <x₁²> * <x₂¹x₁¹> + 1.0 * <x₁¹x₂¹> * <x₂²>"
    end

    @testset "nc state polynomial" begin
        @ncpolyvar x[1:2]
        ncsws = [
            NCStateWord(Arbitrary, wd, ncw) for (wd, ncw) in
            zip([[x[1] * x[2], x[2]^2], [x[1] * x[2]], [x[2]^3]], [x[1], x[2]^2, one(x[1])])
        ]

        spop = ncstatepoly([1.0, 2.0, 3.0], ncsws)

        @test string(spop) ==
              "3.0 * <x₂³> * 1 + 2.0 * <x₁¹x₂¹> * x₂² + 1.0 * <x₁¹x₂¹> * <x₂²> * x₁¹"

        @test degree(spop) == 5
        @test sort(monomials(spop)) == sort(
            map(
                a -> NCStateWord(Arbitrary, a[1], monomial(a[2])),
                [
                    ([x[2]^3], one(x[1])),
                    ([x[2]^2, x[1] * x[2]], x[1]),
                    ([x[1] * x[2]], x[2]^2),
                ],
            ),
        )
    end

    @testset "State Polynomial Arithmetic" begin
        @ncpolyvar x[1:2]

        sps = [
            1.0 * ς(x[1] * x[2]) + 2.0 * ς(x[2] * x[1]), 2.0 * ς(x[2]^2) + 3.0 * ς(x[1]^2)
        ]

        sp_prod =
            2.0 * ς(x[1] * x[2]) * ς(x[2]^2) +
            4.0 * ς(x[2] * x[1]) * ς(x[2]^2) +
            3.0 * ς(x[1] * x[2]) * ς(x[1]^2) +
            6.0 * ς(x[2] * x[1]) * ς(x[1]^2)
        @test sps[1] * sps[2] == sp_prod
        sp_sum =
            1.0 * ς(x[1] * x[2]) + 2.0 * ς(x[2] * x[1]) + 2.0 * ς(x[2]^2) + 3.0 * ς(x[1]^2)
        @test sps[1] + sps[2] == sp_sum
    end

    @testset "variables of NC State Polynomial" begin
        @ncpolyvar x[1:10]

        sp = 1.0 * (ς(x[4] * x[5]) * ς(x[6] * x[7])) + 2.0 * ς(x[8]) + 3.0 * ς(x[9] * x[10])
        @test sort(variables(sp)) == sort(x[4:10])

        words = monomial.([x[1], x[2], x[3]])
        spop =
            1.0 * NCStateWord(Arbitrary, [x[4] * x[5], x[6] * x[7]], x[1]) +
            2.0 * NCStateWord(Arbitrary, [x[8]], x[2]) +
            3.0 * NCStateWord(Arbitrary, [x[9], x[10]], x[3])
        @test variables(spop) == sort(x)
    end

    @testset "State Polynomial Op comparison" begin
        @ncpolyvar a b c
        words_order1 = [a * b, b * c, c^2]
        words_order2 = [b * c, c^2, a * b]
        state_word_order1 = NCStateWord.(Ref(Arbitrary), [[a^2], [b^2], [c^2]], Ref(one(a)))
        state_word_order2 = NCStateWord.(Ref(Arbitrary), [[b^2], [c^2], [a^2]], Ref(one(a)))
        state_poly_order1 = ncstatepoly([1.0, 2.0, 3.0], state_word_order1)
        state_poly_order2 = ncstatepoly([2.0, 3.0, 1.0], state_word_order2)
        @test state_poly_order1 == state_poly_order2
    end
end

@testset "get state basis" begin
    @ncpolyvar x y

    sa = SimplifyAlgorithm(; comm_gps=[[x, y]], is_projective=false, is_unipotent=false)
    c_words = [[one(x)], [y], [x], [one(x)], [one(x)]]
    nc_words = [one(x), one(x), one(x), y, x]
    @test sort(get_state_basis(Arbitrary, [x, y], 1, sa)) ==
          sort(map(x -> NCStateWord(Arbitrary, x[1], x[2]), zip(c_words, nc_words)))

    c_words = [
        [one(x)],
        [y],
        [x],
        [x * y],
        [y^2],
        [x^2],
        [y, y],
        [y, x],
        [x, x],
        [one(x)],
        [y],
        [x],
        [one(x)],
        [y],
        [x],
        fill([one(x)], 4)...,
    ]
    nc_words = [fill(one(x), 9); fill(y, 3); fill(x, 3); [y * x, y^2, x * y, x^2]]

    @test sort(get_state_basis(Arbitrary, [x, y], 2, sa)) ==
          sort(map(x -> NCStateWord(Arbitrary, x[1], x[2]), zip(c_words, nc_words)))

    nc_words = [fill(one(x), 7); fill(x, 4); fill(x^2, 2); [x^3]]
    c_words = [
        [one(x)],
        [x],
        [x^2],
        [x^3],
        [x, x],
        [x, x^2],
        [x, x, x],
        [one(x)],
        [x],
        [x^2],
        [x, x],
        [one(x)],
        [x],
        [one(x)],
    ]
    @test sort(get_state_basis(Arbitrary, [x], 3, sa)) ==
          sort(map(x -> NCStateWord(Arbitrary, x[1], x[2]), zip(c_words, nc_words)))
end

@testset "Arithmetic" begin
    @ncpolyvar x[1:3]
    tp = tr(x[1] * x[2]) - tr(x[1]) * tr(x[2])
    @test string(tp) == "-1.0 * tr(x₁¹) * tr(x₂¹) + 1.0 * tr(x₁¹x₂¹)"

    sp = ς(x[1] * x[2]) - ς(x[1]) * ς(x[2])
    @test string(sp) == "-1.0 * <x₁¹> * <x₂¹> + 1.0 * <x₁¹x₂¹>"
end
