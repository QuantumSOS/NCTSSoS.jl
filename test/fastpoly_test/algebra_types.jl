using Test
using NCTSSoS.FastPolynomials

@testset "AlgebraType Singleton Types" begin
    @testset "NonCommutativeAlgebra" begin
        alg = NonCommutativeAlgebra()
        @test alg isa AlgebraType
        @test alg isa NonCommutativeAlgebra
        # Singleton type - all instances equal
        @test NonCommutativeAlgebra() === NonCommutativeAlgebra()
    end

    @testset "PauliAlgebra" begin
        alg = PauliAlgebra()
        @test alg isa AlgebraType
        @test alg isa PauliAlgebra
        @test PauliAlgebra() === PauliAlgebra()
    end

    @testset "FermionicAlgebra" begin
        alg = FermionicAlgebra()
        @test alg isa AlgebraType
        @test alg isa FermionicAlgebra
        @test FermionicAlgebra() === FermionicAlgebra()
    end

    @testset "BosonicAlgebra" begin
        alg = BosonicAlgebra()
        @test alg isa AlgebraType
        @test alg isa BosonicAlgebra
        @test BosonicAlgebra() === BosonicAlgebra()
    end

    @testset "ProjectorAlgebra" begin
        alg = ProjectorAlgebra()
        @test alg isa AlgebraType
        @test alg isa ProjectorAlgebra
        @test ProjectorAlgebra() === ProjectorAlgebra()
    end

    @testset "UnipotentAlgebra" begin
        alg = UnipotentAlgebra()
        @test alg isa AlgebraType
        @test alg isa UnipotentAlgebra
        @test UnipotentAlgebra() === UnipotentAlgebra()
    end
end

@testset "AlgebraType Show Methods" begin
    @testset "NonCommutativeAlgebra show" begin
        alg = NonCommutativeAlgebra()
        @test sprint(show, alg) == "NonCommutativeAlgebra()"
    end

    @testset "PauliAlgebra show" begin
        alg = PauliAlgebra()
        @test sprint(show, alg) == "PauliAlgebra()"
    end

    @testset "FermionicAlgebra show" begin
        alg = FermionicAlgebra()
        @test sprint(show, alg) == "FermionicAlgebra()"
    end

    @testset "BosonicAlgebra show" begin
        alg = BosonicAlgebra()
        @test sprint(show, alg) == "BosonicAlgebra()"
    end

    @testset "ProjectorAlgebra show" begin
        alg = ProjectorAlgebra()
        @test sprint(show, alg) == "ProjectorAlgebra()"
    end

    @testset "UnipotentAlgebra show" begin
        alg = UnipotentAlgebra()
        @test sprint(show, alg) == "UnipotentAlgebra()"
    end
end

@testset "Index Encoding: site_bits" begin
    @test FastPolynomials.site_bits(UInt8) == 2
    @test FastPolynomials.site_bits(UInt16) == 4
    @test FastPolynomials.site_bits(UInt32) == 8
    @test FastPolynomials.site_bits(UInt64) == 16
end

@testset "Index Encoding: max_sites" begin
    @testset "UInt8 max_sites" begin
        # 2 bits for site: 2^2 - 1 = 3
        @test FastPolynomials.max_sites(UInt8) == 3
    end

    @testset "UInt16 max_sites" begin
        # 4 bits for site: 2^4 - 1 = 15
        @test FastPolynomials.max_sites(UInt16) == 15
    end

    @testset "UInt32 max_sites" begin
        # 8 bits for site: 2^8 - 1 = 255
        @test FastPolynomials.max_sites(UInt32) == 255
    end

    @testset "UInt64 max_sites" begin
        # 16 bits for site: 2^16 - 1 = 65535
        @test FastPolynomials.max_sites(UInt64) == 65535
    end
end

@testset "Index Encoding: max_operators" begin
    @testset "UInt8 max_operators" begin
        # 6 bits for operator: 2^6 - 1 = 63
        @test FastPolynomials.max_operators(UInt8) == 63
    end

    @testset "UInt16 max_operators" begin
        # 12 bits for operator: 2^12 - 1 = 4095
        @test FastPolynomials.max_operators(UInt16) == 4095
    end

    @testset "UInt32 max_operators" begin
        # 24 bits for operator: 2^24 - 1 = 16777215
        @test FastPolynomials.max_operators(UInt32) == 16777215
    end

    @testset "UInt64 max_operators" begin
        # 48 bits for operator: 2^48 - 1 = 281474976710655
        @test FastPolynomials.max_operators(UInt64) == 281474976710655
    end
end

@testset "Index Encoding: encode_index" begin
    @testset "UInt8 encoding" begin
        # UInt8: 2 bits for site, 6 bits for operator_id
        idx = FastPolynomials.encode_index(UInt8, 1, 1)
        @test idx isa UInt8
        @test idx == 0x05  # (1 << 2) | 1 = 5

        idx = FastPolynomials.encode_index(UInt8, 5, 2)
        @test idx == 0x16  # (5 << 2) | 2 = 22

        idx = FastPolynomials.encode_index(UInt8, 63, 3)
        @test idx == 0xFF  # (63 << 2) | 3 = 255
    end

    @testset "UInt16 encoding" begin
        # UInt16: 4 bits for site, 12 bits for operator_id
        idx = FastPolynomials.encode_index(UInt16, 1, 1)
        @test idx isa UInt16
        @test idx == 0x0011  # (1 << 4) | 1

        idx = FastPolynomials.encode_index(UInt16, 100, 10)
        @test idx == 0x064A  # (100 << 4) | 10

        idx = FastPolynomials.encode_index(UInt16, 4095, 15)
        @test idx == 0xFFFF  # (4095 << 4) | 15
    end

    @testset "UInt32 encoding" begin
        # UInt32: 8 bits for site, 24 bits for operator_id
        idx = FastPolynomials.encode_index(UInt32, 1, 1)
        @test idx isa UInt32
        @test idx == 0x00000101  # (1 << 8) | 1

        idx = FastPolynomials.encode_index(UInt32, 1000, 100)
        @test idx == 0x0003E864  # (1000 << 8) | 100

        idx = FastPolynomials.encode_index(UInt32, 16777215, 255)
        @test idx == 0xFFFFFFFF  # (16777215 << 8) | 255
    end

    @testset "UInt64 encoding" begin
        # UInt64: 16 bits for site, 48 bits for operator_id
        idx = FastPolynomials.encode_index(UInt64, 1, 1)
        @test idx isa UInt64
        @test idx == 0x0000000000010001  # (1 << 16) | 1

        idx = FastPolynomials.encode_index(UInt64, 10000, 5000)
        @test idx == 0x0000000027101388  # (10000 << 16) | 5000

        idx = FastPolynomials.encode_index(UInt64, 281474976710655, 65535)
        @test idx == 0xFFFFFFFFFFFFFFFF
    end
end

@testset "Index Encoding: encode_index boundary tests" begin
    @testset "UInt8 boundaries" begin
        # Valid: site 1 to 3, operator_id 1 to 63
        @test_nowarn FastPolynomials.encode_index(UInt8, 1, 1)
        @test_nowarn FastPolynomials.encode_index(UInt8, 63, 3)

        # Invalid: site out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt8, 1, 0)
        @test_throws AssertionError FastPolynomials.encode_index(UInt8, 1, 4)

        # Invalid: operator_id out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt8, 0, 1)
        @test_throws AssertionError FastPolynomials.encode_index(UInt8, 64, 1)
    end

    @testset "UInt16 boundaries" begin
        # Valid: site 1 to 15, operator_id 1 to 4095
        @test_nowarn FastPolynomials.encode_index(UInt16, 1, 1)
        @test_nowarn FastPolynomials.encode_index(UInt16, 4095, 15)

        # Invalid: site out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt16, 1, 0)
        @test_throws AssertionError FastPolynomials.encode_index(UInt16, 1, 16)

        # Invalid: operator_id out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt16, 0, 1)
        @test_throws AssertionError FastPolynomials.encode_index(UInt16, 4096, 1)
    end

    @testset "UInt32 boundaries" begin
        # Valid: site 1 to 255, operator_id 1 to 16777215
        @test_nowarn FastPolynomials.encode_index(UInt32, 1, 1)
        @test_nowarn FastPolynomials.encode_index(UInt32, 16777215, 255)

        # Invalid: site out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt32, 1, 0)
        @test_throws AssertionError FastPolynomials.encode_index(UInt32, 1, 256)

        # Invalid: operator_id out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt32, 0, 1)
        @test_throws AssertionError FastPolynomials.encode_index(UInt32, 16777216, 1)
    end

    @testset "UInt64 boundaries" begin
        # Valid: site 1 to 65535, operator_id 1 to 281474976710655
        @test_nowarn FastPolynomials.encode_index(UInt64, 1, 1)
        @test_nowarn FastPolynomials.encode_index(UInt64, 281474976710655, 65535)

        # Invalid: site out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt64, 1, 0)
        @test_throws AssertionError FastPolynomials.encode_index(UInt64, 1, 65536)

        # Invalid: operator_id out of range
        @test_throws AssertionError FastPolynomials.encode_index(UInt64, 0, 1)
        @test_throws AssertionError FastPolynomials.encode_index(UInt64, 281474976710656, 1)
    end
end

@testset "Index Encoding: decode_site" begin
    @testset "UInt8 decode_site" begin
        idx = FastPolynomials.encode_index(UInt8, 10, 2)
        @test FastPolynomials.decode_site(idx) == 2

        idx = FastPolynomials.encode_index(UInt8, 63, 3)
        @test FastPolynomials.decode_site(idx) == 3

        idx = FastPolynomials.encode_index(UInt8, 1, 1)
        @test FastPolynomials.decode_site(idx) == 1
    end

    @testset "UInt16 decode_site" begin
        idx = FastPolynomials.encode_index(UInt16, 100, 10)
        @test FastPolynomials.decode_site(idx) == 10

        idx = FastPolynomials.encode_index(UInt16, 4095, 15)
        @test FastPolynomials.decode_site(idx) == 15

        idx = FastPolynomials.encode_index(UInt16, 1, 1)
        @test FastPolynomials.decode_site(idx) == 1
    end

    @testset "UInt32 decode_site" begin
        idx = FastPolynomials.encode_index(UInt32, 1000, 100)
        @test FastPolynomials.decode_site(idx) == 100

        idx = FastPolynomials.encode_index(UInt32, 16777215, 255)
        @test FastPolynomials.decode_site(idx) == 255

        idx = FastPolynomials.encode_index(UInt32, 1, 1)
        @test FastPolynomials.decode_site(idx) == 1
    end

    @testset "UInt64 decode_site" begin
        idx = FastPolynomials.encode_index(UInt64, 10000, 5000)
        @test FastPolynomials.decode_site(idx) == 5000

        idx = FastPolynomials.encode_index(UInt64, 281474976710655, 65535)
        @test FastPolynomials.decode_site(idx) == 65535

        idx = FastPolynomials.encode_index(UInt64, 1, 1)
        @test FastPolynomials.decode_site(idx) == 1
    end
end

@testset "Index Encoding: decode_operator_id" begin
    @testset "UInt8 decode_operator_id" begin
        idx = FastPolynomials.encode_index(UInt8, 10, 2)
        @test FastPolynomials.decode_operator_id(idx) == 10

        idx = FastPolynomials.encode_index(UInt8, 63, 3)
        @test FastPolynomials.decode_operator_id(idx) == 63

        idx = FastPolynomials.encode_index(UInt8, 1, 1)
        @test FastPolynomials.decode_operator_id(idx) == 1
    end

    @testset "UInt16 decode_operator_id" begin
        idx = FastPolynomials.encode_index(UInt16, 100, 10)
        @test FastPolynomials.decode_operator_id(idx) == 100

        idx = FastPolynomials.encode_index(UInt16, 4095, 15)
        @test FastPolynomials.decode_operator_id(idx) == 4095

        idx = FastPolynomials.encode_index(UInt16, 1, 1)
        @test FastPolynomials.decode_operator_id(idx) == 1
    end

    @testset "UInt32 decode_operator_id" begin
        idx = FastPolynomials.encode_index(UInt32, 1000, 100)
        @test FastPolynomials.decode_operator_id(idx) == 1000

        idx = FastPolynomials.encode_index(UInt32, 16777215, 255)
        @test FastPolynomials.decode_operator_id(idx) == 16777215

        idx = FastPolynomials.encode_index(UInt32, 1, 1)
        @test FastPolynomials.decode_operator_id(idx) == 1
    end

    @testset "UInt64 decode_operator_id" begin
        idx = FastPolynomials.encode_index(UInt64, 10000, 5000)
        @test FastPolynomials.decode_operator_id(idx) == 10000

        idx = FastPolynomials.encode_index(UInt64, 281474976710655, 65535)
        @test FastPolynomials.decode_operator_id(idx) == 281474976710655

        idx = FastPolynomials.encode_index(UInt64, 1, 1)
        @test FastPolynomials.decode_operator_id(idx) == 1
    end
end

@testset "Index Encoding: round-trip tests" begin
    @testset "UInt8 round-trip" begin
        for op_id in [1, 10, 32, 63]
            for site in [1, 2, 3]
                idx = FastPolynomials.encode_index(UInt8, op_id, site)
                @test FastPolynomials.decode_operator_id(idx) == op_id
                @test FastPolynomials.decode_site(idx) == site
            end
        end
    end

    @testset "UInt16 round-trip" begin
        for op_id in [1, 100, 1000, 4095]
            for site in [1, 7, 15]
                idx = FastPolynomials.encode_index(UInt16, op_id, site)
                @test FastPolynomials.decode_operator_id(idx) == op_id
                @test FastPolynomials.decode_site(idx) == site
            end
        end
    end

    @testset "UInt32 round-trip" begin
        for op_id in [1, 1000, 100000, 16777215]
            for site in [1, 100, 255]
                idx = FastPolynomials.encode_index(UInt32, op_id, site)
                @test FastPolynomials.decode_operator_id(idx) == op_id
                @test FastPolynomials.decode_site(idx) == site
            end
        end
    end

    @testset "UInt64 round-trip" begin
        for op_id in [1, 10000, 1000000, 100000000, 281474976710655]
            for site in [1, 1000, 10000, 65535]
                idx = FastPolynomials.encode_index(UInt64, op_id, site)
                @test FastPolynomials.decode_operator_id(idx) == op_id
                @test FastPolynomials.decode_site(idx) == site
            end
        end
    end
end

@testset "Index Encoding: select_uint_type" begin
    @testset "UInt8 selection" begin
        # UInt8: max 63 operators, max 3 sites
        @test FastPolynomials.select_uint_type(10, 3) == UInt8
        @test FastPolynomials.select_uint_type(63, 3) == UInt8
        @test FastPolynomials.select_uint_type(1, 1) == UInt8
        @test FastPolynomials.select_uint_type(50, 2) == UInt8
    end

    @testset "UInt16 selection" begin
        # UInt16: max 4095 operators, max 15 sites
        @test FastPolynomials.select_uint_type(100, 4) == UInt16
        @test FastPolynomials.select_uint_type(64, 3) == UInt16  # exceeds UInt8 operators
        @test FastPolynomials.select_uint_type(10, 4) == UInt16  # exceeds UInt8 sites
        @test FastPolynomials.select_uint_type(4095, 15) == UInt16
        @test FastPolynomials.select_uint_type(1000, 10) == UInt16
    end

    @testset "UInt32 selection" begin
        # UInt32: max 16777215 operators, max 255 sites
        @test FastPolynomials.select_uint_type(5000, 16) == UInt32  # exceeds UInt16 sites
        @test FastPolynomials.select_uint_type(4096, 10) == UInt32  # exceeds UInt16 operators
        @test FastPolynomials.select_uint_type(100000, 100) == UInt32
        @test FastPolynomials.select_uint_type(16777215, 255) == UInt32
    end

    @testset "UInt64 selection" begin
        # UInt64: max 281474976710655 operators, max 65535 sites
        @test FastPolynomials.select_uint_type(100, 256) == UInt64  # exceeds UInt32 sites
        @test FastPolynomials.select_uint_type(16777216, 100) == UInt64  # exceeds UInt32 operators
        @test FastPolynomials.select_uint_type(1000000, 1000) == UInt64
        @test FastPolynomials.select_uint_type(100000000, 10000) == UInt64
    end

    @testset "select_uint_type error handling" begin
        # Too many sites for UInt64
        @test_throws ErrorException FastPolynomials.select_uint_type(100, 65536)

        # Too many operators for UInt64 (using large value that doesn't overflow max_operators calculation)
        @test_throws ErrorException FastPolynomials.select_uint_type(281474976710656, 100)

        # Both exceed UInt64 capacity
        @test_throws ErrorException FastPolynomials.select_uint_type(281474976710656, 65536)
    end

    @testset "select_uint_type error message validation" begin
        # Verify error messages are informative
        try
            FastPolynomials.select_uint_type(100, 65536)
            @test false  # Should not reach here
        catch e
            @test e isa ErrorException
            @test occursin("Cannot fit", e.msg)
            @test occursin("65536", e.msg)  # sites value in message
        end

        try
            FastPolynomials.select_uint_type(281474976710656, 100)
            @test false  # Should not reach here
        catch e
            @test e isa ErrorException
            @test occursin("Cannot fit", e.msg)
            @test occursin("281474976710656", e.msg)  # operators value in message
        end
    end
end

@testset "Index Encoding: edge case combinations" begin
    @testset "Minimum values" begin
        # Test minimum valid values (1, 1) for all types
        for T in [UInt8, UInt16, UInt32, UInt64]
            idx = FastPolynomials.encode_index(T, 1, 1)
            @test FastPolynomials.decode_operator_id(idx) == 1
            @test FastPolynomials.decode_site(idx) == 1
        end
    end

    @testset "Maximum values" begin
        # Test maximum valid values for each type
        idx = FastPolynomials.encode_index(UInt8, 63, 3)
        @test FastPolynomials.decode_operator_id(idx) == 63
        @test FastPolynomials.decode_site(idx) == 3

        idx = FastPolynomials.encode_index(UInt16, 4095, 15)
        @test FastPolynomials.decode_operator_id(idx) == 4095
        @test FastPolynomials.decode_site(idx) == 15

        idx = FastPolynomials.encode_index(UInt32, 16777215, 255)
        @test FastPolynomials.decode_operator_id(idx) == 16777215
        @test FastPolynomials.decode_site(idx) == 255

        idx = FastPolynomials.encode_index(UInt64, 281474976710655, 65535)
        @test FastPolynomials.decode_operator_id(idx) == 281474976710655
        @test FastPolynomials.decode_site(idx) == 65535
    end

    @testset "Mixed boundary values" begin
        # Test max operator with min site and vice versa
        idx = FastPolynomials.encode_index(UInt16, 4095, 1)
        @test FastPolynomials.decode_operator_id(idx) == 4095
        @test FastPolynomials.decode_site(idx) == 1

        idx = FastPolynomials.encode_index(UInt16, 1, 15)
        @test FastPolynomials.decode_operator_id(idx) == 1
        @test FastPolynomials.decode_site(idx) == 15
    end
end

@testset "Index Encoding: integration with simplification" begin
    @testset "ProjectorAlgebra encoding preserved through simplification" begin
        # Create encoded indices for projector operators on different sites
        p1_site1 = FastPolynomials.encode_index(UInt16, 1, 1)  # P₁ on site 1
        p1_site2 = FastPolynomials.encode_index(UInt16, 1, 2)  # P₁ on site 2

        # Create monomial P₁¹ P₁¹ P₁² (two P₁ on site 1, one on site 2)
        m = Monomial{ProjectorAlgebra}(UInt16[p1_site1, p1_site1, p1_site2])

        # Simplify: P² = P, so P₁¹ P₁¹ → P₁¹
        result = simplify(m)

        # Verify result preserves encoding (simplify returns Monomial for ProjectorAlgebra)
        @test result.word[1] == p1_site1  # First should be site 1
        @test result.word[2] == p1_site2  # Second should be site 2
        @test FastPolynomials.decode_site(result.word[1]) == 1
        @test FastPolynomials.decode_site(result.word[2]) == 2
    end

    @testset "UnipotentAlgebra encoding preserved through simplification" begin
        # Create encoded indices for unipotent operators on different sites
        u1_site1 = FastPolynomials.encode_index(UInt16, 1, 1)  # U₁ on site 1
        u1_site2 = FastPolynomials.encode_index(UInt16, 1, 2)  # U₁ on site 2

        # Create monomial U₁¹ U₁¹ U₁² (two U₁ on site 1, one on site 2)
        m = Monomial{UnipotentAlgebra}(UInt16[u1_site1, u1_site1, u1_site2])

        # Simplify: U² = I, so U₁¹ U₁¹ → identity (removed)
        result = simplify(m)

        # After U₁¹ U₁¹ cancels, only U₁² remains (simplify returns Monomial for UnipotentAlgebra)
        @test length(result.word) == 1
        @test result.word[1] == u1_site2
        @test FastPolynomials.decode_site(result.word[1]) == 2
    end

    @testset "NonCommutativeAlgebra site-based commutation" begin
        # Create encoded indices for operators on different sites
        x1_site1 = FastPolynomials.encode_index(UInt16, 1, 1)  # x₁ on site 1
        x2_site2 = FastPolynomials.encode_index(UInt16, 2, 2)  # x₂ on site 2
        x3_site1 = FastPolynomials.encode_index(UInt16, 3, 1)  # x₃ on site 1

        # Create monomial x₂² x₁¹ x₃¹ (out of site order)
        m = Monomial{NonCommutativeAlgebra}(UInt16[x2_site2, x1_site1, x3_site1])

        # Simplify should sort by site (operators on different sites commute)
        result = simplify(m)

        # Site 1 operators should come before site 2 (simplify returns Monomial for NonCommutativeAlgebra)
        @test FastPolynomials.decode_site(result.word[1]) == 1
        @test FastPolynomials.decode_site(result.word[2]) == 1
        @test FastPolynomials.decode_site(result.word[3]) == 2
    end
end
