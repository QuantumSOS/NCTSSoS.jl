using Test, NCTSSoS
# Import internal functions for testing
using NCTSSoS:
    _subscript_string,
    _select_pauli_type,
    _select_signed_index_type,
    max_operators,
    max_sites

# Test helper: verify all expected symbols exist in registry
function check_symbols_in_registry(reg, expected_symbols)
    for sym in expected_symbols
        @test sym in reg
    end
end

# Test helper: verify variable arrays have correct type and length
function check_variable_arrays(vars, expected_type, expected_lengths)
    for (var_array, expected_len) in zip(vars, expected_lengths)
        @test length(var_array) == expected_len
        @test all(v -> v isa NormalMonomial{expected_type}, var_array)
    end
end

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
        # Multi-digit numbers
        @test _subscript_string(42) == "₄₂"
        @test _subscript_string(123) == "₁₂₃"
        @test _subscript_string(9999) == "₉₉₉₉"

        # Negative numbers should error
        @test_throws ErrorException _subscript_string(-1)
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
        @test _select_pauli_type(85) == UInt8    # 3*85 = 255 <= typemax(UInt8)
        @test _select_pauli_type(86) == UInt16   # 3*86 = 258 > typemax(UInt8)

        # _select_signed_index_type
        @test _select_signed_index_type(127) == Int8
        @test _select_signed_index_type(128) == Int16

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
        # Non-contiguous subscripts
        reg_nc, (P,) = create_projector_variables([("P", [1, 5, 100])])
        check_symbols_in_registry(reg_nc, [:P₁, :P₅, :P₁₀₀])

        # Large subscripts
        reg_large, (U,) = create_unipotent_variables([("U", [9999])])
        @test :U₉₉₉₉ in reg_large
    end
    @testset "Non-Commutative Variable Creation" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        check_variable_arrays([x], NonCommutativeAlgebra, [3])
        check_symbols_in_registry(reg, [:x₁, :x₂, :x₃])
    end

    @testset "Multi-Prefix Variable Creation" begin
        reg, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 3:4)])
        check_variable_arrays([x, y], NonCommutativeAlgebra, [2, 2])
        check_symbols_in_registry(reg, [:x₁, :x₂, :y₃, :y₄])
    end

    @testset "Pauli Variable Creation" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:2)
        check_variable_arrays([σx, σy, σz], PauliAlgebra, [2, 2, 2])
        check_symbols_in_registry(reg, [:σx₁, :σy₁, :σz₁, :σx₂, :σy₂, :σz₂])
    end

    @testset "Projector Variable Creation" begin
        reg, (P,) = create_projector_variables([("P", 1:3)])
        check_variable_arrays([P], ProjectorAlgebra, [3])
        check_symbols_in_registry(reg, [:P₁, :P₂, :P₃])
    end

    @testset "Unipotent Variable Creation" begin
        reg, (U,) = create_unipotent_variables([("U", 1:3)])
        check_variable_arrays([U], UnipotentAlgebra, [3])
        check_symbols_in_registry(reg, [:U₁, :U₂, :U₃])
    end

    @testset "Fermionic Variable Creation" begin
        reg, (a, a_dag) = create_fermionic_variables(1:2)
        check_variable_arrays([a, a_dag], FermionicAlgebra, [2, 2])
        check_symbols_in_registry(reg, [:a₁, Symbol("a⁺₁"), :a₂, Symbol("a⁺₂")])
    end

    @testset "Bosonic Variable Creation" begin
        reg, (c, c_dag) = create_bosonic_variables(1:2)
        check_variable_arrays([c, c_dag], BosonicAlgebra, [2, 2])
        check_symbols_in_registry(reg, [:c₁, Symbol("c⁺₁"), :c₂, Symbol("c⁺₂")])
    end

    @testset "Monomial Power" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # x^2 should create a polynomial with degree 2
        x2 = x[1] * x[1]
        @test x2 isa Polynomial
        @test degree(x2) == 2

        # x^0 should be identity
        x0 = one(x[1])
        @test isone(x0)
    end

    @testset "subregistry" begin
        # Basic functionality: create subset of registry
        reg, (x,) = create_noncommutative_variables([("x", 1:5)])
        sorted_idxs = sort(collect(keys(reg.idx_to_variables)))

        sub_reg = subregistry(reg, sorted_idxs[1:3])
        @test length(sub_reg) == 3
        @test sub_reg isa VariableRegistry{NonCommutativeAlgebra}

        # Check symbols are preserved
        for idx in sorted_idxs[1:3]
            @test reg[idx] == sub_reg[idx]
        end

        # Check symbols not in subset are missing
        @test_throws KeyError sub_reg[sorted_idxs[4]]
        @test_throws KeyError sub_reg[sorted_idxs[5]]

        # Subregistry with Pauli algebra preserves type
        reg_pauli, _ = create_pauli_variables(1:3)
        pauli_idxs = sort(collect(keys(reg_pauli.idx_to_variables)))
        sub_pauli = subregistry(reg_pauli, pauli_idxs[1:4])
        @test sub_pauli isa VariableRegistry{PauliAlgebra}
        @test length(sub_pauli) == 4

        # Empty subset yields empty registry
        empty_sub = subregistry(reg, eltype(sorted_idxs)[])
        @test length(empty_sub) == 0

        # Indices not in parent are silently ignored
        fake_idxs = [sorted_idxs[1], eltype(sorted_idxs)(255)]  # 255 unlikely to exist
        sub_with_fake = subregistry(reg, fake_idxs)
        @test length(sub_with_fake) == 1
        @test sub_with_fake[sorted_idxs[1]] == reg[sorted_idxs[1]]

        # Full subset equals original (same content)
        full_sub = subregistry(reg, sorted_idxs)
        @test length(full_sub) == length(reg)
        for idx in sorted_idxs
            @test full_sub[idx] == reg[idx]
        end
    end
end
