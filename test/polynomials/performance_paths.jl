using NCTSSoS, Test

function _reference_product_terms(
    scale,
    m1::NormalMonomial{A,T},
    m2::NormalMonomial{A,T},
    ::Type{C},
) where {A<:AlgebraType,T<:Integer,C<:Number}
    out = Tuple{C,NormalMonomial{A,T}}[]
    if isempty(m1.word)
        iszero(scale) || push!(out, (C(scale), m2))
        return out
    elseif isempty(m2.word)
        iszero(scale) || push!(out, (C(scale), m1))
        return out
    end

    simplified = simplify(A, vcat(m1.word, m2.word))
    for (coef, mono) in NCTSSoS._simplified_to_terms(A, simplified, T)
        total = C(scale) * C(coef)
        iszero(total) || push!(out, (total, mono))
    end
    return out
end

function _reference_product(
    m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    C = coeff_type(A)
    return Polynomial(_reference_product_terms(one(C), m1, m2, C))
end

function _reference_left_product(
    m::NormalMonomial{A,T}, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    NC = promote_type(C, coeff_type(A))
    terms = Tuple{NC,NormalMonomial{A,T}}[]
    for (coef, mono) in p.terms
        append!(terms, _reference_product_terms(coef, m, mono, NC))
    end
    return Polynomial(terms)
end

function _reference_right_product(
    p::Polynomial{A,T,C}, m::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    NC = promote_type(C, coeff_type(A))
    terms = Tuple{NC,NormalMonomial{A,T}}[]
    for (coef, mono) in p.terms
        append!(terms, _reference_product_terms(coef, mono, m, NC))
    end
    return Polynomial(terms)
end

function _reference_polynomial_product(
    p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}
) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    NC = promote_type(C1, C2, coeff_type(A))
    terms = Tuple{NC,NormalMonomial{A,T}}[]
    for (coef1, mono1) in p1.terms
        for (coef2, mono2) in p2.terms
            append!(terms, _reference_product_terms(coef1 * coef2, mono1, mono2, NC))
        end
    end
    return Polynomial(terms)
end

function _canonical_monomials(::Type{A}, raw_words::Vector{Vector{T}}) where {A<:AlgebraType,T<:Integer}
    monos = NormalMonomial{A,T}[one(NormalMonomial{A,T})]
    for raw in raw_words
        for (coef, mono) in NCTSSoS._simplified_to_terms(A, simplify(A, raw), T)
            iszero(coef) || push!(monos, mono)
        end
    end
    return unique(monos)
end

function _exercise_multiplication_family(::Type{A}, monos::Vector{NormalMonomial{A,T}}) where {A<:AlgebraType,T<:Integer}
    @test length(monos) >= 4

    for m1 in monos, m2 in monos
        @test m1 * m2 == _reference_product(m1, m2)
    end

    C = coeff_type(A)
    input_terms = Tuple{C,NormalMonomial{A,T}}[
        (C(2), monos[2]),
        (C(-1), monos[2]),
        (zero(C), monos[3]),
        (C(3), monos[3]),
        (C(-3), monos[3]),
        (C(4), monos[4]),
    ]
    input_snapshot = copy(input_terms)
    p = Polynomial(input_terms)
    @test input_terms == input_snapshot  # public constructor must not sort! caller storage

    parametric_input = copy(input_terms)
    p_parametric = Polynomial{A,T,C}(parametric_input)
    @test parametric_input == input_snapshot

    expected_p = Polynomial(Tuple{C,NormalMonomial{A,T}}[(C(1), monos[2]), (C(4), monos[4])])
    @test p == expected_p
    @test p_parametric == expected_p

    q = Polynomial(Tuple{C,NormalMonomial{A,T}}[
        (C(-2), monos[4]),
        (C(5), monos[2]),
        (C(1), monos[1]),
    ])

    for m in monos[1:4]
        @test m * p == _reference_left_product(m, p)
        @test p * m == _reference_right_product(p, m)
    end
    @test p * q == _reference_polynomial_product(p, q)
end

function _hash_allocated(m)
    hash(m, UInt(0x12345678))
    return @allocated hash(m, UInt(0x12345678))
end

@testset "Performance-sensitive monomial paths" begin
    @testset "Differential multiplication across algebras" begin
        nc_words = Vector{UInt16}[
            UInt16[nc_idx(UInt16, 2, 1)],
            UInt16[nc_idx(UInt16, 2, 2), nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 1, 1), nc_idx(UInt16, 3, 1)],
            UInt16[nc_idx(UInt16, 1, 3), nc_idx(UInt16, 2, 1)],
        ]
        _exercise_multiplication_family(NonCommutativeAlgebra, _canonical_monomials(NonCommutativeAlgebra, nc_words))

        proj_words = Vector{UInt16}[
            UInt16[nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 1, 1), nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 1, 2), nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 2, 1), nc_idx(UInt16, 1, 1)],
        ]
        _exercise_multiplication_family(ProjectorAlgebra, _canonical_monomials(ProjectorAlgebra, proj_words))

        unip_words = Vector{UInt16}[
            UInt16[nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 1, 1), nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 1, 2), nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 2, 1), nc_idx(UInt16, 1, 1)],
            UInt16[nc_idx(UInt16, 1, 3), nc_idx(UInt16, 2, 1)],
        ]
        _exercise_multiplication_family(UnipotentAlgebra, _canonical_monomials(UnipotentAlgebra, unip_words))

        pauli_words = Vector{UInt8}[
            UInt8[1],
            UInt8[2, 1],
            UInt8[1, 4],
            UInt8[3, 5, 1],
            UInt8[2, 5, 6],
        ]
        _exercise_multiplication_family(PauliAlgebra, _canonical_monomials(PauliAlgebra, pauli_words))

        fermion_words = Vector{Int8}[
            Int8[-1],
            Int8[1],
            Int8[-1, 1],
            Int8[1, -1],
            Int8[-2, -1, 1],
            Int8[2, -2, 1],
        ]
        _exercise_multiplication_family(FermionicAlgebra, _canonical_monomials(FermionicAlgebra, fermion_words))

        boson_words = Vector{Int8}[
            Int8[-1],
            Int8[1],
            Int8[-1, 1],
            Int8[1, -1],
            Int8[-2, -1, 1],
            Int8[2, -2, 1],
        ]
        _exercise_multiplication_family(BosonicAlgebra, _canonical_monomials(BosonicAlgebra, boson_words))
    end

    @testset "NormalMonomial hash contract" begin
        seed1 = UInt(0x12345678)
        seed2 = UInt(0x87654321)
        m8 = NormalMonomial{PauliAlgebra,UInt8}(UInt8[1, 4, 7])
        m16 = NormalMonomial{PauliAlgebra,UInt16}(UInt16[1, 4, 7])
        @test m8 == m16
        @test isequal(m8, m16)
        @test hash(m8, seed1) == hash(m16, seed1)
        @test hash(m8, seed1) != hash(m8, seed2)

        # Warm the method before asserting runtime allocation behavior.
        _hash_allocated(m8)
        @test _hash_allocated(m8) == 0

        same_unsigned_word = UInt16[1]
        unsigned_hashes = [
            hash(NormalMonomial{NonCommutativeAlgebra,UInt16}(same_unsigned_word), seed1),
            hash(NormalMonomial{ProjectorAlgebra,UInt16}(same_unsigned_word), seed1),
            hash(NormalMonomial{UnipotentAlgebra,UInt16}(same_unsigned_word), seed1),
            hash(NormalMonomial{PauliAlgebra,UInt16}(same_unsigned_word), seed1),
        ]
        @test length(unique(unsigned_hashes)) == length(unsigned_hashes)

        f = NormalMonomial{FermionicAlgebra,Int32}(Int32[-1, 1])
        b = NormalMonomial{BosonicAlgebra,Int32}(Int32[-1, 1])
        @test f != b
        @test hash(f, seed1) != hash(b, seed1)
    end
end
