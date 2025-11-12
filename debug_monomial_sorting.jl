#!/usr/bin/env julia

"""
Debug script for monomial sorting and comparison logic.

This script investigates:
1. The order in which get_basis generates monomials
2. Whether the comparison function implements correct lexicographical ordering
3. Discrepancies between generated order and expected order
"""

using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: Variable, get_basis, monomial

function print_monomial_details(m::NCTSSoS.FastPolynomials.Monomial)
    println("  Monomial: $(m)")
    println("    vars: $(m.vars)")
    println("    exps: $(m.z)")
    println("    degree: $(sum(m.z))")
    println()
end

function main()
    println("="^80)
    println("MONOMIAL SORTING DEBUG SCRIPT")
    println("="^80)
    println()

    # Create variables
    @ncpolyvar x y z

    println("Variables created: x, y, z")
    println("Variable comparison results:")
    println("  x < y: $(cmp(x, y) < 0)")
    println("  y < z: $(cmp(y, z) < 0)")
    println("  x < z: $(cmp(x, z) < 0)")
    println()

    # Get basis up to degree 2
    println("-"^80)
    println("PART 1: Generated Basis from get_basis([x, y, z], 2)")
    println("-"^80)
    println()

    basis = get_basis([x, y, z], 2)

    println("Generated order ($(length(basis)) monomials):")
    for (i, m) in enumerate(basis)
        print("  [$i] ")
        print_monomial_details(m)
    end

    # Expected lexicographical order
    println("-"^80)
    println("PART 2: Expected Lexicographical Order")
    println("-"^80)
    println()

    expected = [
        one(x),                    # 1 (degree 0)
        monomial([x], [1]),        # x (degree 1)
        monomial([y], [1]),        # y
        monomial([z], [1]),        # z
        x^2,                       # x² (degree 2)
        x * y,                     # xy
        x * z,                     # xz
        y * x,                     # yx
        y^2,                       # y²
        y * z,                     # yz
        z * x,                     # zx
        z * y,                     # zy
        z^2,                       # z²
    ]

    println("Expected order ($(length(expected)) monomials):")
    for (i, m) in enumerate(expected)
        print("  [$i] ")
        print_monomial_details(m)
    end

    # Check if basis is sorted
    println("-"^80)
    println("PART 3: Sorting Analysis")
    println("-"^80)
    println()

    println("Is generated basis sorted? $(issorted(basis))")
    println("Does generated basis match expected? $(basis == expected)")
    println()

    # Compare each position
    println("Position-by-position comparison:")
    for i in 1:min(length(basis), length(expected))
        match = basis[i] == expected[i]
        status = match ? "✓" : "✗"
        println("  [$i] $status  Generated: $(basis[i])  |  Expected: $(expected[i])")
        if !match
            println("      cmp(generated, expected) = $(cmp(basis[i], expected[i]))")
        end
    end
    println()

    # Test comparison function with specific pairs
    println("-"^80)
    println("PART 4: Testing Comparison Function")
    println("-"^80)
    println()

    test_pairs = [
        (monomial([x], [1]), monomial([y], [1]), "x vs y"),
        (monomial([x], [2]), monomial([x, y], [1, 1]), "x² vs xy"),
        (monomial([x, y], [1, 1]), monomial([x, z], [1, 1]), "xy vs xz"),
        (monomial([x, y], [1, 1]), monomial([y, x], [1, 1]), "xy vs yx"),
        (monomial([y], [2]), monomial([y, z], [1, 1]), "y² vs yz"),
        (monomial([z, x], [1, 1]), monomial([z, y], [1, 1]), "zx vs zy"),
    ]

    println("Testing specific monomial pairs:")
    for (a, b, desc) in test_pairs
        cmp_result = cmp(a, b)
        cmp_str = cmp_result < 0 ? "a < b" : (cmp_result > 0 ? "a > b" : "a == b")
        println("  $(desc): cmp($(a), $(b)) = $cmp_result  ($cmp_str)")
    end
    println()

    # Sort the generated basis and compare
    println("-"^80)
    println("PART 5: Sorted Generated Basis")
    println("-"^80)
    println()

    sorted_basis = sort(basis)
    println("Sorted generated basis:")
    for (i, m) in enumerate(sorted_basis)
        print("  [$i] ")
        print_monomial_details(m)
    end

    println("Does sorted basis match expected? $(sorted_basis == expected)")
    println()

    # Analyze comparison logic issues
    println("-"^80)
    println("PART 6: Deep Dive into Comparison Logic")
    println("-"^80)
    println()

    println("Testing degree-2 monomials with same degree:")
    deg2_monos = filter(m -> sum(m.z) == 2, expected)
    println("\nDegree-2 monomials in expected order:")
    for (i, m) in enumerate(deg2_monos)
        println("  [$i] $(m)")
    end

    println("\nPairwise comparisons of consecutive degree-2 monomials:")
    for i in 1:(length(deg2_monos) - 1)
        a, b = deg2_monos[i], deg2_monos[i + 1]
        cmp_result = cmp(a, b)
        is_correct = cmp_result < 0
        status = is_correct ? "✓" : "✗ WRONG!"
        println("  $status  $(a) < $(b): cmp = $cmp_result")

        if !is_correct
            println("      Analysis:")
            println("        $(a).vars = $(a.vars), $(a).z = $(a.z)")
            println("        $(b).vars = $(b.vars), $(b).z = $(b.z)")
        end
    end
    println()

    println("="^80)
    println("END OF DEBUG SCRIPT")
    println("="^80)
end

# Run the main function
main()
