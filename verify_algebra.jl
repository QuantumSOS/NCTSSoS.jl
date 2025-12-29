println("=" ^ 80)
println("COMPREHENSIVE VERIFICATION: NCTSSoS.jl")
println("=" ^ 80)

using NCTSSoS
using LinearAlgebra

issues_found = String[]

# ============================================================================
# 1. PAULI ALGEBRA VERIFICATION
# ============================================================================
println("\n" * "=" ^ 40)
println("1. PAULI ALGEBRA")
println("=" ^ 40)

reg, (σx, σy, σz) = create_pauli_variables(1:3)
println("✓ Created 3-site Pauli variables")

# Test: σ² = I (identity)
println("\nTest: σᵢ² = I")
for (name, σ) in [("σx", σx), ("σy", σy), ("σz", σz)]
    poly = Polynomial(σ[1]) * Polynomial(σ[1])
    if poly == one(poly)
        println("  ✓ $(name)₁² = 1")
    else
        msg = "FAIL: $(name)₁² = $poly (expected 1)"
        println("  ✗ $msg")
        push!(issues_found, msg)
    end
end

# Test: σx σy = i σz (at same site)
println("\nTest: σx₁ σy₁ = i σz₁")
prod_xy = Polynomial(σx[1]) * Polynomial(σy[1])
expected_xy = im * Polynomial(σz[1])
if prod_xy == expected_xy
    println("  ✓ σx₁ σy₁ = i σz₁")
else
    msg = "FAIL: σx₁ σy₁ = $prod_xy (expected $expected_xy)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: σy σx = -i σz (anticommutation)
println("\nTest: σy₁ σx₁ = -i σz₁")
prod_yx = Polynomial(σy[1]) * Polynomial(σx[1])
expected_yx = -im * Polynomial(σz[1])
if prod_yx == expected_yx
    println("  ✓ σy₁ σx₁ = -i σz₁")
else
    msg = "FAIL: σy₁ σx₁ = $prod_yx (expected $expected_yx)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: Different sites commute
println("\nTest: [σx₁, σy₂] = 0 (different sites commute)")
comm = Polynomial(σx[1]) * Polynomial(σy[2]) - Polynomial(σy[2]) * Polynomial(σx[1])
if iszero(comm)
    println("  ✓ [σx₁, σy₂] = 0")
else
    msg = "FAIL: [σx₁, σy₂] = $comm (expected 0)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: σx σy σz = i I (cycle)
println("\nTest: σx₁ σy₁ σz₁ = i")
cycle = Polynomial(σx[1]) * Polynomial(σy[1]) * Polynomial(σz[1])
if cycle == im * one(cycle)
    println("  ✓ σx₁ σy₁ σz₁ = i")
else
    msg = "FAIL: σx₁ σy₁ σz₁ = $cycle (expected i)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# 2. FERMIONIC ALGEBRA VERIFICATION
# ============================================================================
println("\n" * "=" ^ 40)
println("2. FERMIONIC ALGEBRA")
println("=" ^ 40)

reg_f, (a, a⁺) = create_fermionic_variables(1:3)
println("✓ Created 3-mode fermionic variables")

# Test: a² = 0 (nilpotent)
println("\nTest: aᵢ² = 0 (nilpotent)")
a_sq = Polynomial(a[1]) * Polynomial(a[1])
if iszero(a_sq)
    println("  ✓ a₁² = 0")
else
    msg = "FAIL: a₁² = $a_sq (expected 0)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: (a⁺)² = 0
println("\nTest: (a⁺ᵢ)² = 0")
adag_sq = Polynomial(a⁺[1]) * Polynomial(a⁺[1])
if iszero(adag_sq)
    println("  ✓ (a⁺₁)² = 0")
else
    msg = "FAIL: (a⁺₁)² = $adag_sq (expected 0)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: {aᵢ, a⁺ᵢ} = 1 (anticommutation at same site)
println("\nTest: {a₁, a⁺₁} = 1")
anticomm = Polynomial(a[1]) * Polynomial(a⁺[1]) + Polynomial(a⁺[1]) * Polynomial(a[1])
if anticomm == one(anticomm)
    println("  ✓ {a₁, a⁺₁} = 1")
else
    msg = "FAIL: {a₁, a⁺₁} = $anticomm (expected 1)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: {aᵢ, a⁺ⱼ} = 0 for i ≠ j
println("\nTest: {a₁, a⁺₂} = 0 (different modes)")
anticomm_diff = Polynomial(a[1]) * Polynomial(a⁺[2]) + Polynomial(a⁺[2]) * Polynomial(a[1])
if iszero(anticomm_diff)
    println("  ✓ {a₁, a⁺₂} = 0")
else
    msg = "FAIL: {a₁, a⁺₂} = $anticomm_diff (expected 0)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: Normal ordering - a⁺a (number operator)
println("\nTest: Normal ordering preserved for a⁺₁ a₁")
num_op = Polynomial(a⁺[1]) * Polynomial(a[1])
println("  a⁺₁ a₁ = $num_op")

# ============================================================================
# 3. BOSONIC ALGEBRA VERIFICATION
# ============================================================================
println("\n" * "=" ^ 40)
println("3. BOSONIC ALGEBRA")
println("=" ^ 40)

reg_b, (c, c⁺) = create_bosonic_variables(1:3)
println("✓ Created 3-mode bosonic variables")

# Test: [cᵢ, c⁺ᵢ] = 1 (commutation at same site)
println("\nTest: [c₁, c⁺₁] = 1")
comm_b = Polynomial(c[1]) * Polynomial(c⁺[1]) - Polynomial(c⁺[1]) * Polynomial(c[1])
if comm_b == one(comm_b)
    println("  ✓ [c₁, c⁺₁] = 1")
else
    msg = "FAIL: [c₁, c⁺₁] = $comm_b (expected 1)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: [cᵢ, c⁺ⱼ] = 0 for i ≠ j
println("\nTest: [c₁, c⁺₂] = 0 (different modes)")
comm_b_diff = Polynomial(c[1]) * Polynomial(c⁺[2]) - Polynomial(c⁺[2]) * Polynomial(c[1])
if iszero(comm_b_diff)
    println("  ✓ [c₁, c⁺₂] = 0")
else
    msg = "FAIL: [c₁, c⁺₂] = $comm_b_diff (expected 0)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: [cᵢ, cⱼ] = 0
println("\nTest: [c₁, c₂] = 0")
comm_cc = Polynomial(c[1]) * Polynomial(c[2]) - Polynomial(c[2]) * Polynomial(c[1])
if iszero(comm_cc)
    println("  ✓ [c₁, c₂] = 0")
else
    msg = "FAIL: [c₁, c₂] = $comm_cc (expected 0)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# 4. PROJECTOR ALGEBRA VERIFICATION
# ============================================================================
println("\n" * "=" ^ 40)
println("4. PROJECTOR ALGEBRA")
println("=" ^ 40)

reg_p, (P,) = create_projector_variables([("P", 1:3)])
println("✓ Created 3 projector variables")

# Test: P² = P (idempotent)
println("\nTest: P₁² = P₁ (idempotent)")
P_sq = Polynomial(P[1]) * Polynomial(P[1])
if P_sq == Polynomial(P[1])
    println("  ✓ P₁² = P₁")
else
    msg = "FAIL: P₁² = $P_sq (expected P₁)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: Different projectors don't simplify (non-commutative)
println("\nTest: P₁ P₂ ≠ P₂ P₁ (non-commutative)")
P12 = Polynomial(P[1]) * Polynomial(P[2])
P21 = Polynomial(P[2]) * Polynomial(P[1])
if P12 != P21
    println("  ✓ P₁ P₂ ≠ P₂ P₁")
else
    println("  ⚠ P₁ P₂ = P₂ P₁ (unexpected but not necessarily wrong)")
end

# ============================================================================
# 5. UNIPOTENT ALGEBRA VERIFICATION
# ============================================================================
println("\n" * "=" ^ 40)
println("5. UNIPOTENT ALGEBRA")
println("=" ^ 40)

reg_u, (U,) = create_unipotent_variables([("U", 1:3)])
println("✓ Created 3 unipotent variables")

# Test: U² = I (involutory)
println("\nTest: U₁² = I (involutory)")
U_sq = Polynomial(U[1]) * Polynomial(U[1])
if U_sq == one(U_sq)
    println("  ✓ U₁² = 1")
else
    msg = "FAIL: U₁² = $U_sq (expected 1)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# 6. NONCOMMUTATIVE ALGEBRA VERIFICATION
# ============================================================================
println("\n" * "=" ^ 40)
println("6. NONCOMMUTATIVE ALGEBRA (generic)")
println("=" ^ 40)

reg_nc, (x,) = create_noncommutative_variables([("x", 1:3)])
println("✓ Created 3 generic NC variables")

# Test: No simplification (x₁ x₁ stays as-is at polynomial level)
println("\nTest: x₁ x₁ has no automatic simplification")
x_sq = Polynomial(x[1]) * Polynomial(x[1])
println("  x₁² = $x_sq")
# In generic NC, x₁² is NOT simplified
deg = degree(x_sq)
if deg == 2
    println("  ✓ Degree is 2 (no simplification)")
else
    msg = "FAIL: degree(x₁²) = $deg (expected 2)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# Test: Order preserved
println("\nTest: x₁ x₂ ≠ x₂ x₁ (order preserved)")
x12 = Polynomial(x[1]) * Polynomial(x[2])
x21 = Polynomial(x[2]) * Polynomial(x[1])
if x12 != x21
    println("  ✓ x₁ x₂ ≠ x₂ x₁")
else
    msg = "FAIL: x₁ x₂ = x₂ x₁ (should be different)"
    println("  ✗ $msg")
    push!(issues_found, msg)
end

# ============================================================================
# SUMMARY
# ============================================================================
println("\n" * "=" ^ 80)
println("VERIFICATION SUMMARY")
println("=" ^ 80)

if isempty(issues_found)
    println("✅ ALL ALGEBRA TESTS PASSED")
else
    println("❌ ISSUES FOUND:")
    for issue in issues_found
        println("  - $issue")
    end
end

# Return exit code
exit(isempty(issues_found) ? 0 : 1)
