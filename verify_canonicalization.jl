println("=" ^ 80)
println("CANONICALIZATION VERIFICATION")
println("=" ^ 80)

using NCTSSoS

issues_found = String[]

# ============================================================================
# 1. SYMMETRIC CANONICALIZATION
# ============================================================================
println("\n" * "=" ^ 40)
println("1. SYMMETRIC CANONICALIZATION")
println("=" ^ 40)

println("\nTest: symmetric_canon picks lexicographically smaller of word vs reverse")

# Test case 1: [1,2,3] < [3,2,1], should return [1,2,3]
word1 = [1, 2, 3]
result1 = symmetric_canon(word1)
if result1 == [1, 2, 3]
    println("  ✓ symmetric_canon([1,2,3]) = [1,2,3]")
else
    msg = "FAIL: symmetric_canon([1,2,3]) = $result1 (expected [1,2,3])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test case 2: [3,2,1] should return [1,2,3] (reverse is smaller)
word2 = [3, 2, 1]
result2 = symmetric_canon(word2)
if result2 == [1, 2, 3]
    println("  ✓ symmetric_canon([3,2,1]) = [1,2,3]")
else
    msg = "FAIL: symmetric_canon([3,2,1]) = $result2 (expected [1,2,3])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test case 3: [1,2,1] palindrome - should return as-is
word3 = [1, 2, 1]
result3 = symmetric_canon(word3)
if result3 == [1, 2, 1]
    println("  ✓ symmetric_canon([1,2,1]) = [1,2,1] (palindrome)")
else
    msg = "FAIL: symmetric_canon([1,2,1]) = $result3 (expected [1,2,1])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test case 4: [2,1,3] vs [3,1,2] - [2,1,3] < [3,1,2]
word4 = [2, 1, 3]
result4 = symmetric_canon(word4)
if result4 == [2, 1, 3]  # 2 < 3
    println("  ✓ symmetric_canon([2,1,3]) = [2,1,3] (word < reverse)")
else
    msg = "FAIL: symmetric_canon([2,1,3]) = $result4 (expected [2,1,3])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Physical interpretation: For Hermitian observables, ⟨ψ|ABC|ψ⟩ = ⟨ψ|C†B†A†|ψ⟩*
# For Pauli matrices (self-adjoint), this means σxσyσz ≡ σzσyσx
println("\nPhysical test: Pauli observable canonicalization")
reg, (σx, σy, σz) = create_pauli_variables(1:3)
m_xyz = Monomial{PauliAlgebra}([UInt8(1), UInt8(2), UInt8(3)])  # σx₁σy₁σz₁ encoding
m_zyx = Monomial{PauliAlgebra}([UInt8(3), UInt8(2), UInt8(1)])  # σz₁σy₁σx₁ encoding

canon_xyz = symmetric_canon(m_xyz)
canon_zyx = symmetric_canon(m_zyx)

if canon_xyz == canon_zyx
    println("  ✓ symmetric_canon(σx₁σy₁σz₁) = symmetric_canon(σz₁σy₁σx₁)")
else
    msg = "FAIL: Symmetric monomials should have same canonical form"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# 2. CYCLIC CANONICALIZATION
# ============================================================================
println("\n" * "=" ^ 40)
println("2. CYCLIC CANONICALIZATION")
println("=" ^ 40)

println("\nTest: cyclic_canon finds minimum cyclic rotation")

# Test case 1: [2,3,1] rotations are [2,3,1], [3,1,2], [1,2,3] - minimum is [1,2,3]
word_c1 = [2, 3, 1]
result_c1 = cyclic_canon(word_c1)
if result_c1 == [1, 2, 3]
    println("  ✓ cyclic_canon([2,3,1]) = [1,2,3]")
else
    msg = "FAIL: cyclic_canon([2,3,1]) = $result_c1 (expected [1,2,3])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test case 2: [3,1,2] rotations are [3,1,2], [1,2,3], [2,3,1] - minimum is [1,2,3]
word_c2 = [3, 1, 2]
result_c2 = cyclic_canon(word_c2)
if result_c2 == [1, 2, 3]
    println("  ✓ cyclic_canon([3,1,2]) = [1,2,3]")
else
    msg = "FAIL: cyclic_canon([3,1,2]) = $result_c2 (expected [1,2,3])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test case 3: [1,2,3] already minimal
word_c3 = [1, 2, 3]
result_c3 = cyclic_canon(word_c3)
if result_c3 == [1, 2, 3]
    println("  ✓ cyclic_canon([1,2,3]) = [1,2,3] (already minimal)")
else
    msg = "FAIL: cyclic_canon([1,2,3]) = $result_c3 (expected [1,2,3])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test case 4: [1,1,2] rotations are [1,1,2], [1,2,1], [2,1,1]
word_c4 = [2, 1, 1]
result_c4 = cyclic_canon(word_c4)
if result_c4 == [1, 1, 2]
    println("  ✓ cyclic_canon([2,1,1]) = [1,1,2]")
else
    msg = "FAIL: cyclic_canon([2,1,1]) = $result_c4 (expected [1,1,2])"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Physical interpretation: For trace optimization, Tr(ABC) = Tr(BCA) = Tr(CAB)
println("\nPhysical test: Trace equivalence")
println("  Tr(σx₁σy₂σz₃) = Tr(σy₂σz₃σx₁) = Tr(σz₃σx₁σy₂)")
m_trace1 = Monomial{PauliAlgebra}([UInt8(1), UInt8(5), UInt8(9)])  # σx₁ σy₂ σz₃
m_trace2 = Monomial{PauliAlgebra}([UInt8(5), UInt8(9), UInt8(1)])  # σy₂ σz₃ σx₁
m_trace3 = Monomial{PauliAlgebra}([UInt8(9), UInt8(1), UInt8(5)])  # σz₃ σx₁ σy₂

canon_t1 = cyclic_canon(m_trace1)
canon_t2 = cyclic_canon(m_trace2)
canon_t3 = cyclic_canon(m_trace3)

if canon_t1 == canon_t2 == canon_t3
    println("  ✓ All cyclic rotations have the same canonical form")
else
    msg = "FAIL: Cyclic rotations should have same canonical form"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# 3. CYCLIC-SYMMETRIC CANONICALIZATION
# ============================================================================
println("\n" * "=" ^ 40)
println("3. CYCLIC-SYMMETRIC CANONICALIZATION")
println("=" ^ 40)

println("\nTest: cyclic_symmetric_canon combines both transformations")

# Test case: [3,2,1] and [1,2,3] should have same canonical form
# Also [2,1,3], [1,3,2], [3,2,1] etc.
using NCTSSoS: cyclic_symmetric_canon

word_cs1 = [3, 2, 1]
word_cs2 = [1, 2, 3]
result_cs1 = cyclic_symmetric_canon(word_cs1)
result_cs2 = cyclic_symmetric_canon(word_cs2)

if result_cs1 == result_cs2
    println("  ✓ cyclic_symmetric_canon([3,2,1]) = cyclic_symmetric_canon([1,2,3])")
else
    msg = "FAIL: Should be equal: $result_cs1 vs $result_cs2"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test [3,1,4,2] as in docstring
word_cs3 = [3, 1, 4, 2]
result_cs3 = cyclic_symmetric_canon(word_cs3)
expected_cs3 = [1, 3, 2, 4]  # from docstring
if result_cs3 == expected_cs3
    println("  ✓ cyclic_symmetric_canon([3,1,4,2]) = [1,3,2,4]")
else
    msg = "FAIL: cyclic_symmetric_canon([3,1,4,2]) = $result_cs3 (expected $expected_cs3)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# 4. POLYNOMIAL CANONICALIZATION
# ============================================================================
println("\n" * "=" ^ 40)
println("4. POLYNOMIAL CANONICALIZATION")
println("=" ^ 40)

println("\nTest: Canonicalizing polynomial combines equivalent terms")

# Create two monomials that are symmetric equivalents
m_a = Monomial{PauliAlgebra}([UInt8(1), UInt8(2), UInt8(3)])
m_b = Monomial{PauliAlgebra}([UInt8(3), UInt8(2), UInt8(1)])  # reverse

# Create polynomial with both terms
p = Polynomial([Term(1.0+0im, m_a), Term(2.0+0im, m_b)])
println("  Original: $p (2 terms)")

p_canon = canonicalize(p)
println("  Canonicalized: $p_canon")

if length(terms(p_canon)) == 1
    coeff = coefficients(p_canon)[1]
    if coeff == 3.0+0im
        println("  ✓ Terms combined: coefficient = 3.0")
    else
        msg = "FAIL: Expected coefficient 3.0, got $coeff"
        println("  ✗ $msg")
        push!(issues_found, msg)
    end
else
    msg = "FAIL: Expected 1 term after canonicalization, got $(length(terms(p_canon)))"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test cyclic canonicalization on polynomial
println("\nTest: Cyclic canonicalization on polynomial")
m_c1 = Monomial{PauliAlgebra}([UInt8(1), UInt8(2), UInt8(3)])  # ABC
m_c2 = Monomial{PauliAlgebra}([UInt8(2), UInt8(3), UInt8(1)])  # BCA (cyclic rotation)

p_cyc = Polynomial([Term(1.0+0im, m_c1), Term(2.0+0im, m_c2)])
println("  Original: $p_cyc (2 terms)")

p_cyc_canon = canonicalize(p_cyc; cyclic=true)
println("  Cyclic canonicalized: $p_cyc_canon")

if length(terms(p_cyc_canon)) == 1
    println("  ✓ Cyclic equivalent terms combined")
else
    msg = "FAIL: Cyclic equivalent terms should combine"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# SUMMARY
# ============================================================================
println("\n" * "=" ^ 80)
println("CANONICALIZATION VERIFICATION SUMMARY")
println("=" ^ 80)

if isempty(issues_found)
    println("✅ ALL CANONICALIZATION TESTS PASSED")
else
    println("❌ ISSUES FOUND:")
    for issue in issues_found
        println("  - $issue")
    end
end

exit(isempty(issues_found) ? 0 : 1)
