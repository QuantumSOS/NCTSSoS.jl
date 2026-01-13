using Test
using NCTSSoS
using NCTSSoS: coeff_type

# Test data for algebra types
const ALGEBRA_TYPES = [
    NonCommutativeAlgebra,
    PauliAlgebra,
    FermionicAlgebra,
    BosonicAlgebra,
    ProjectorAlgebra,
    UnipotentAlgebra,
]

# Expected coeff_type for each algebra (PauliAlgebra is special)
const ALGEBRA_COEFF_TYPES = Dict(
    PauliAlgebra => ComplexF64,
    NonCommutativeAlgebra => Float64,
    FermionicAlgebra => Float64,
    BosonicAlgebra => Float64,
    ProjectorAlgebra => Float64,
    UnipotentAlgebra => Float64,
)

# Index types for NormalMonomial (signed for PBW algebras)
const ALGEBRA_INDEX_TYPES = Dict(
    PauliAlgebra => UInt16,
    NonCommutativeAlgebra => UInt16,
    FermionicAlgebra => Int32,
    BosonicAlgebra => Int32,
    ProjectorAlgebra => UInt16,
    UnipotentAlgebra => UInt16,
)

# UInt types and their encoding parameters
const UINT_TYPES = [UInt8, UInt16, UInt32, UInt64]
const UINT_PARAMS = Dict(
    UInt8 => (site_bits=2, max_sites=3, max_ops=63),
    UInt16 => (site_bits=4, max_sites=15, max_ops=4095),
    UInt32 => (site_bits=8, max_sites=255, max_ops=16777215),
    UInt64 => (site_bits=16, max_sites=65535, max_ops=281474976710655),
)

# Test cases for encode_index
const ENCODE_TEST_CASES = Dict(
    UInt8 => [
        (op=1, site=1, expected=0x05),
        (op=5, site=2, expected=0x16),
        (op=63, site=3, expected=0xFF),
    ],
    UInt16 => [
        (op=1, site=1, expected=0x0011),
        (op=100, site=10, expected=0x064A),
        (op=4095, site=15, expected=0xFFFF),
    ],
    UInt32 => [
        (op=1, site=1, expected=0x00000101),
        (op=1000, site=100, expected=0x0003E864),
        (op=16777215, site=255, expected=0xFFFFFFFF),
    ],
    UInt64 => [
        (op=1, site=1, expected=0x0000000000010001),
        (op=10000, site=5000, expected=0x0000000027101388),
        (op=281474976710655, site=65535, expected=0xFFFFFFFFFFFFFFFF),
    ],
)

# Round-trip test values
const ROUNDTRIP_TEST_VALUES = Dict(
    UInt8 => (ops=[1, 10, 32, 63], sites=[1, 2, 3]),
    UInt16 => (ops=[1, 100, 1000, 4095], sites=[1, 7, 15]),
    UInt32 => (ops=[1, 1000, 100000, 16777215], sites=[1, 100, 255]),
    UInt64 => (ops=[1, 10000, 1000000, 100000000, 281474976710655], sites=[1, 1000, 10000, 65535]),
)

@testset "coeff_type trait for algebras" begin
    for (A, expected) in ALGEBRA_COEFF_TYPES
        @test coeff_type(A) == expected
    end
end

@testset "coeff_type trait for types" begin
    # NormalMonomial: returns coeff_type(A)
    for (A, expected) in ALGEBRA_COEFF_TYPES
        T = ALGEBRA_INDEX_TYPES[A]
        @test coeff_type(NormalMonomial{A,T}) == expected
    end

    # Polynomial: returns C (the coefficient type)
    @test coeff_type(Polynomial{PauliAlgebra,UInt16,ComplexF64}) == ComplexF64
    @test coeff_type(Polynomial{NonCommutativeAlgebra,UInt16,Float64}) == Float64
    @test coeff_type(Polynomial{BosonicAlgebra,Int32,Float64}) == Float64

    T = UInt16
    # Instance methods
    m = NormalMonomial{PauliAlgebra,T}(T[1, 4])
    @test coeff_type(m) == ComplexF64

    p = Polynomial([(1.0 + 0im, NormalMonomial{PauliAlgebra,T}(T[1]))])
    @test coeff_type(p) == ComplexF64
end

@testset "Index Encoding: site_bits" begin
    for T in UINT_TYPES
        @test NCTSSoS.site_bits(T) == UINT_PARAMS[T].site_bits
    end
end

@testset "Index Encoding: max_sites" begin
    for T in UINT_TYPES
        @test NCTSSoS.max_sites(T) == UINT_PARAMS[T].max_sites
    end
end

@testset "Index Encoding: max_operators" begin
    for T in UINT_TYPES
        @test NCTSSoS.max_operators(T) == UINT_PARAMS[T].max_ops
    end
end

@testset "Index Encoding: encode_index" begin
    for T in UINT_TYPES
        @testset "$T encoding" begin
            for tc in ENCODE_TEST_CASES[T]
                idx = NCTSSoS.encode_index(T, tc.op, tc.site)
                @test idx == tc.expected
            end
        end
    end
end

@testset "Index Encoding: encode_index boundary tests" begin
    for T in UINT_TYPES
        params = UINT_PARAMS[T]
        @testset "$T boundaries" begin
            # Valid boundaries
            @test_nowarn NCTSSoS.encode_index(T, 1, 1)
            @test_nowarn NCTSSoS.encode_index(T, params.max_ops, params.max_sites)

            # Invalid: site out of range
            @test_throws ArgumentError NCTSSoS.encode_index(T, 1, 0)
            @test_throws ArgumentError NCTSSoS.encode_index(T, 1, params.max_sites + 1)

            # Invalid: operator_id out of range
            @test_throws ArgumentError NCTSSoS.encode_index(T, 0, 1)
            @test_throws ArgumentError NCTSSoS.encode_index(T, params.max_ops + 1, 1)
        end
    end
end

@testset "Index Encoding: decode_site" begin
    for T in UINT_TYPES
        @testset "$T decode_site" begin
            for tc in ENCODE_TEST_CASES[T]
                idx = NCTSSoS.encode_index(T, tc.op, tc.site)
                @test NCTSSoS.decode_site(idx) == tc.site
            end
        end
    end
end

@testset "Index Encoding: decode_operator_id" begin
    for T in UINT_TYPES
        @testset "$T decode_operator_id" begin
            for tc in ENCODE_TEST_CASES[T]
                idx = NCTSSoS.encode_index(T, tc.op, tc.site)
                @test NCTSSoS.decode_operator_id(idx) == tc.op
            end
        end
    end
end

@testset "Index Encoding: round-trip tests" begin
    for T in UINT_TYPES
        @testset "$T round-trip" begin
            vals = ROUNDTRIP_TEST_VALUES[T]
            for op_id in vals.ops
                for site in vals.sites
                    idx = NCTSSoS.encode_index(T, op_id, site)
                    @test NCTSSoS.decode_operator_id(idx) == op_id
                    @test NCTSSoS.decode_site(idx) == site
                end
            end
        end
    end
end

@testset "Index Encoding: select_uint_type" begin
    @testset "UInt8 selection" begin
        @test NCTSSoS.select_uint_type(10, 3) == UInt8
        @test NCTSSoS.select_uint_type(63, 3) == UInt8
        @test NCTSSoS.select_uint_type(1, 1) == UInt8
        @test NCTSSoS.select_uint_type(50, 2) == UInt8
    end

    @testset "UInt16 selection" begin
        @test NCTSSoS.select_uint_type(100, 4) == UInt16
        @test NCTSSoS.select_uint_type(64, 3) == UInt16   # exceeds UInt8 operators
        @test NCTSSoS.select_uint_type(10, 4) == UInt16   # exceeds UInt8 sites
        @test NCTSSoS.select_uint_type(4095, 15) == UInt16
        @test NCTSSoS.select_uint_type(1000, 10) == UInt16
    end

    @testset "UInt32 selection" begin
        @test NCTSSoS.select_uint_type(5000, 16) == UInt32    # exceeds UInt16 sites
        @test NCTSSoS.select_uint_type(4096, 10) == UInt32    # exceeds UInt16 operators
        @test NCTSSoS.select_uint_type(100000, 100) == UInt32
        @test NCTSSoS.select_uint_type(16777215, 255) == UInt32
    end

    @testset "UInt64 selection" begin
        @test NCTSSoS.select_uint_type(100, 256) == UInt64        # exceeds UInt32 sites
        @test NCTSSoS.select_uint_type(16777216, 100) == UInt64   # exceeds UInt32 operators
        @test NCTSSoS.select_uint_type(1000000, 1000) == UInt64
        @test NCTSSoS.select_uint_type(100000000, 10000) == UInt64
    end

    @testset "select_uint_type error handling" begin
        @test_throws ArgumentError NCTSSoS.select_uint_type(100, 65536)
        @test_throws ArgumentError NCTSSoS.select_uint_type(281474976710656, 100)
        @test_throws ArgumentError NCTSSoS.select_uint_type(281474976710656, 65536)
    end

    @testset "select_uint_type error message validation" begin
        try
            NCTSSoS.select_uint_type(100, 65536)
            @test false
        catch e
            @test e isa ArgumentError
            @test occursin("Cannot fit", e.msg)
            @test occursin("65536", e.msg)
        end

        try
            NCTSSoS.select_uint_type(281474976710656, 100)
            @test false
        catch e
            @test e isa ArgumentError
            @test occursin("Cannot fit", e.msg)
            @test occursin("281474976710656", e.msg)
        end
    end
end

@testset "Index Encoding: edge case combinations" begin
    @testset "Maximum values" begin
        for T in UINT_TYPES
            params = UINT_PARAMS[T]
            idx = NCTSSoS.encode_index(T, params.max_ops, params.max_sites)
            @test NCTSSoS.decode_operator_id(idx) == params.max_ops
            @test NCTSSoS.decode_site(idx) == params.max_sites
        end
    end

    @testset "Mixed boundary values" begin
        idx = NCTSSoS.encode_index(UInt16, 4095, 1)
        @test NCTSSoS.decode_operator_id(idx) == 4095
        @test NCTSSoS.decode_site(idx) == 1

        idx = NCTSSoS.encode_index(UInt16, 1, 15)
        @test NCTSSoS.decode_operator_id(idx) == 1
        @test NCTSSoS.decode_site(idx) == 15
    end
end
