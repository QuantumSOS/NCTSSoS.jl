using Test, NCTSSoS.FastPolynomials
# Import internal functions for testing
using NCTSSoS.FastPolynomials:
    _subscript_string,
    _select_pauli_type,
    _select_signed_index_type,
    _select_contiguous_type

@testset "Variable Registry" begin

    @testset "VariableRegistry Invariant Validation" begin
        # Should error: Different lengths
        @test_throws ErrorException VariableRegistry{NonCommutativeAlgebra, UInt8}(
            Dict{UInt8,Symbol}(0x01 => :x),
            Dict{Symbol,UInt8}(:x => 0x01, :y => 0x02)
        )

        # Should error: Mismatched indices
        @test_throws ErrorException VariableRegistry{NonCommutativeAlgebra, UInt8}(
            Dict{UInt8,Symbol}(0x01 => :x),
            Dict{Symbol,UInt8}(:x => 0x02)  # Wrong index
        )

        # Should succeed: Consistent registry
        reg = VariableRegistry{NonCommutativeAlgebra, UInt8}(
            Dict{UInt8,Symbol}(0x01 => :x, 0x02 => :y),
            Dict{Symbol,UInt8}(:x => 0x01, :y => 0x02)
        )
        @test length(reg) == 2

        # Empty registry should be valid
        empty_reg = VariableRegistry{NonCommutativeAlgebra, UInt8}(
            Dict{UInt8,Symbol}(),
            Dict{Symbol,UInt8}()
        )
        @test length(empty_reg) == 0
    end

    @testset "_subscript_string Edge Cases" begin
        # Zero
        @test _subscript_string(0) == "₀"

        # Single digits
        @test _subscript_string(1) == "₁"
        @test _subscript_string(9) == "₉"

        # Multi-digit numbers
        @test _subscript_string(10) == "₁₀"
        @test _subscript_string(42) == "₄₂"
        @test _subscript_string(123) == "₁₂₃"

        # Large numbers
        @test _subscript_string(9999) == "₉₉₉₉"

        # Negative numbers should error
        @test_throws ErrorException _subscript_string(-1)
        @test_throws ErrorException _subscript_string(-42)
    end

    @testset "Base.getindex Error Handling" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Valid access by symbol then by index
        idx = reg[:x₁]
        @test reg[idx] == :x₁

        # Integer type conversion works for valid indices
        @test reg[Int(idx)] isa Symbol

        # Missing symbol should throw KeyError
        @test_throws KeyError reg[:nonexistent]

        # Invalid index should throw KeyError
        @test_throws KeyError reg[UInt8(0)]
    end

    @testset "Base.show Display" begin
        # Single variable (singular)
        reg1, _ = create_noncommutative_variables([("x", 1:1)])
        str1 = string(reg1)
        @test contains(str1, "1 variable")
        @test !contains(str1, "variables")  # Not plural

        # Multiple variables (plural)
        reg3, _ = create_noncommutative_variables([("x", 1:3)])
        str3 = string(reg3)
        @test contains(str3, "3 variables")

        # More than 10 variables (ellipsis)
        reg15, _ = create_noncommutative_variables([("x", 1:15)])
        str15 = string(reg15)
        @test contains(str15, "15 variables")
        @test contains(str15, "...")
        @test contains(str15, "x₁")   # First shown
        @test contains(str15, "x₁₅")  # Last shown

        # Empty registry
        empty_reg = VariableRegistry{NonCommutativeAlgebra, UInt8}(Dict{UInt8,Symbol}(), Dict{Symbol,UInt8}())
        str_empty = string(empty_reg)
        @test contains(str_empty, "0 variables")
    end

    @testset "Type Selection Functions" begin
        # _select_pauli_type: 3 operators per site
        @test _select_pauli_type(85) == UInt8    # 3*85 = 255 ≤ typemax(UInt8)
        @test _select_pauli_type(86) == UInt16   # 3*86 = 258 > typemax(UInt8)

        # _select_signed_index_type
        @test _select_signed_index_type(127) == Int8
        @test _select_signed_index_type(128) == Int16

        # _select_contiguous_type
        @test _select_contiguous_type(255) == UInt8
        @test _select_contiguous_type(256) == UInt16
    end

    @testset "symbols() and indices() Ordering" begin
        # Create registry and verify ordering
        reg, _ = create_pauli_variables(1:3)

        # symbols() should return in index-sorted order
        syms = symbols(reg)
        @test length(syms) == 9  # 3 sites × 3 Pauli types
        @test syms[1] == :σx₁
        @test syms[2] == :σy₁
        @test syms[3] == :σz₁

        # indices() should return sorted
        idxs = indices(reg)
        @test issorted(idxs)
        @test length(idxs) == 9
    end

    @testset "Empty and Edge Case Collections" begin
        # Single element
        reg1, (σx, σy, σz) = create_pauli_variables([1])
        @test length(reg1) == 3
        @test length(σx) == 1

        # Non-contiguous subscripts
        reg_nc, (P,) = create_projector_variables([("P", [1, 5, 100])])
        @test :P₁ in reg_nc
        @test :P₅ in reg_nc
        @test :P₁₀₀ in reg_nc
        @test length(P) == 3

        # Large subscripts
        reg_large, (U,) = create_unipotent_variables([("U", [9999])])
        @test :U₉₉₉₉ in reg_large
    end
    @testset "Non-Commutative Variable Creation" begin
        # Create variables with single prefix
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        @test reg isa VariableRegistry
        @test length(x) == 3

        # Check that monomials are created correctly
        @test x[1] isa Monomial{NonCommutativeAlgebra}
        @test x[2] isa Monomial{NonCommutativeAlgebra}
        @test x[3] isa Monomial{NonCommutativeAlgebra}

        # Check registry contains symbols
        @test :x₁ in reg
        @test :x₂ in reg
        @test :x₃ in reg
    end

    @testset "Multi-Prefix Variable Creation" begin
        reg, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 3:4)])

        @test length(x) == 2
        @test length(y) == 2

        @test :x₁ in reg
        @test :x₂ in reg
        @test :y₃ in reg
        @test :y₄ in reg
    end

    @testset "Pauli Variable Creation" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:2)

        @test length(σx) == 2
        @test length(σy) == 2
        @test length(σz) == 2

        @test σx[1] isa Monomial{PauliAlgebra}
        @test σy[1] isa Monomial{PauliAlgebra}
        @test σz[1] isa Monomial{PauliAlgebra}

        # Check registry symbols
        @test :σx₁ in reg
        @test :σy₁ in reg
        @test :σz₁ in reg
        @test :σx₂ in reg
        @test :σy₂ in reg
        @test :σz₂ in reg
    end

    @testset "Projector Variable Creation" begin
        reg, (P,) = create_projector_variables([("P", 1:3)])

        @test length(P) == 3
        @test P[1] isa Monomial{ProjectorAlgebra}

        @test :P₁ in reg
        @test :P₂ in reg
        @test :P₃ in reg
    end

    @testset "Unipotent Variable Creation" begin
        reg, (U,) = create_unipotent_variables([("U", 1:3)])

        @test length(U) == 3
        @test U[1] isa Monomial{UnipotentAlgebra}

        @test :U₁ in reg
        @test :U₂ in reg
        @test :U₃ in reg
    end

    @testset "Fermionic Variable Creation" begin
        reg, (a, a_dag) = create_fermionic_variables(1:2)

        @test length(a) == 2
        @test length(a_dag) == 2
        @test a[1] isa Monomial{FermionicAlgebra}
        @test a_dag[1] isa Monomial{FermionicAlgebra}

        # Check both annihilation and creation operators
        @test :a₁ in reg
        @test Symbol("a⁺₁") in reg
        @test :a₂ in reg
        @test Symbol("a⁺₂") in reg
    end

    @testset "Bosonic Variable Creation" begin
        reg, (c, c_dag) = create_bosonic_variables(1:2)

        @test length(c) == 2
        @test length(c_dag) == 2
        @test c[1] isa Monomial{BosonicAlgebra}
        @test c_dag[1] isa Monomial{BosonicAlgebra}

        # Check both annihilation and creation operators
        @test :c₁ in reg
        @test Symbol("c⁺₁") in reg
    end

    @testset "Registry Access" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:5)])

        # Test symbols function
        syms = symbols(reg)
        @test length(syms) == 5
        @test :x₁ in syms
        @test :x₅ in syms

        # Test indices function
        idxs = indices(reg)
        @test length(idxs) == 5
    end

    @testset "Monomial Power" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # x^2 should create a monomial with word [idx, idx]
        x2 = x[1] * x[1]
        @test x2 isa Monomial  # Multiplication now returns Monomial
        @test degree(x2) == 2

        # x^0 should be identity
        x0 = one(x[1])
        @test isone(x0)
    end
end
