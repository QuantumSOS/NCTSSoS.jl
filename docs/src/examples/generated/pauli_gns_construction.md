<!-- nctssos-literate-source: pauli_gns_construction.jl sha256: b76158188e144f8ef743079f82a0ea66d63fe0ab1ecd60c95e1d7602e5d04423 -->

```@meta
EditURL = "../literate/pauli_gns_construction.jl"
```

# GNS Construction for Pauli Operator Reconstruction

The Gelfand-Naimark-Segal (GNS) construction is a fundamental mathematical tool
in quantum mechanics that allows us to represent abstract quantum states and operators
as concrete matrices acting on a Hilbert space. For physicists, this provides a systematic
way to reconstruct operator representations from expectation values (moments).

## Background: From Expectation Values to Matrix Representations

In quantum mechanics, we often know the expectation values of operators in a given state:
```math
\langle A  \rangle = \text{Tr}(\rho A)
```
where ``\rho`` is the density matrix and ``A`` is an operator.

The GNS construction answers the question: **Can we reconstruct the actual matrices
representing operators from just these expectation values?**

The key insight is that the collection of all expectation values defines a **moment matrix**
(also called a Hankel matrix in the context of polynomial optimization):
```math
H_{ij} = \langle b_i^\dagger b_j \rangle
```
where ``\{b_i\}`` is a basis of operators (monomials in our variables).

## Pauli Operators: A Simple Quantum System

Let's start with the simplest non-trivial quantum system: a single qubit with Pauli operators.
The Pauli matrices are:
- ``\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0  \end{pmatrix}``
- ``\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}``
- ``\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}```

These satisfy the fundamental Pauli algebra:
1. ``\sigma_x^2 = \sigma_y^2 = \sigma_z^2 = I`` (squares equal identity)
2. ``\{\sigma_i, \sigma_j\} = 0`` for ``i \neq j`` (anti-commutation)
3. ``[\sigma_i, \sigma_j] = 2i\epsilon_{ijk}\sigma_k`` (commutation relations)

## Step 1: Define Non-commuting Variables

````@example pauli_gns_construction
using NCTSSoS
using LinearAlgebra
using LinearAlgebra: tr
````

In `NCTSSoS.jl`, we represent Pauli operators as non-commuting polynomial variables:
Create non-commuting variables using the typed algebra system

````@example pauli_gns_construction
registry, (x, y, z) = create_noncommutative_variables([("x", 1:1), ("y", 1:1), ("z", 1:1)])
````

Extract single variables from arrays

````@example pauli_gns_construction
x, y, z = x[1], y[1], z[1]
````

These variables x, y, z will represent σₓ, σᵧ, σ_z respectively

````@example pauli_gns_construction
vars = [x, y, z];
nothing #hide
````

## Step 2: Choose a Quantum State

````@example pauli_gns_construction
zero_state = ComplexF64[1; 0];    # |0⟩
one_state = ComplexF64[0; 1];
nothing #hide
````

Define quantum states for testing
For clear reconstruction, use a pure state

````@example pauli_gns_construction
ρ =  zero_state * zero_state'
````

## Step 3: Compute Expectation Values

We need a function to compute the expectation value of any monomial in our variables:

````@example pauli_gns_construction
function expval_pauli(mono::Monomial, ρ::Matrix)
    if isone(mono)
        return 1.0 + 0.0im  # ⟨I⟩ = 1 for normalized states
    end

    # Start with identity matrix
    mat = Matrix{ComplexF64}(I, 2, 2)

    # Multiply the appropriate Pauli matrices
    for (var, exp) in zip(mono.vars, mono.z)
        pauli_mat = if var == x
            ComplexF64[0 1; 1 0]      # σₓ
        elseif var == y
            ComplexF64[0 -im; im 0]   # σᵧ
        elseif var == z
            ComplexF64[1 0; 0 -1]     # σ_z
        else
            error("Unknown variable: $var")
        end

        # Raise to the appropriate power
        for _ in 1:exp
            mat = mat * pauli_mat
        end
    end

    # Compute Tr(ρ * mat)
    return tr(ρ * mat)
end
````

## Step 4: Build the Moment (Hankel) Matrix

The moment matrix encodes all expectation values of products of our basis operators:

Choose the degree of our polynomial basis

````@example pauli_gns_construction
degree = 4;
nothing #hide
````

Generate the basis of monomials up to the specified degree

````@example pauli_gns_construction
basis = NCTSSoS.get_basis(vars, degree)

println("Basis operators (monomials):")
for (i, b) in enumerate(basis)
    println("$i: $b")
end
````

Build the moment matrix H where H[i,j] = ⟨b_i† * b_j⟩

````@example pauli_gns_construction
n = length(basis)
H = zeros(ComplexF64, n, n)

for i in 1:n
    for j in 1:n
        mono_i = basis[i]
        mono_j = basis[j]

        # Compute the product monomial
        product = NCTSSoS.neat_dot(mono_i, mono_j)

        # Get the expectation value
        H[i, j] = expval_pauli(product, ρ)
    end
end
````

Verify that H is Hermitian (as required by quantum mechanics)

````@example pauli_gns_construction
println("Moment matrix H is Hermitian: ", H ≈ H')
````

## Step 5: GNS Reconstruction

Now we use the `reconstruct` function to perform the GNS construction and obtain
concrete matrix representations of our abstract operators:

````@example pauli_gns_construction
X_recon, Y_recon, Z_recon = reconstruct(H, vars, degree; atol=0.001)
````

````@example pauli_gns_construction
println("Reconstructed Pauli operators:")
println("σₓ (reconstructed):")
round.(X_recon, digits=6)
````

````@example pauli_gns_construction
println("σᵧ (reconstructed):")
round.(Y_recon, digits=6)
````

````@example pauli_gns_construction
println("σ_z (reconstructed):")
round.(Z_recon, digits=6)
````

## Step 6: Verify Pauli Algebra

The true test of our reconstruction is whether the recovered operators satisfy
the Pauli algebra relations:

### Test 1: Squares should equal identity (X² = Y² = Z² = I)

````@example pauli_gns_construction
X2 = X_recon * X_recon
Y2 = Y_recon * Y_recon
Z2 = Z_recon * Z_recon

println("   ||X² - I|| = $(norm(X2 - I))")
println("   ||Y² - I|| = $(norm(Y2 - I))")
println("   ||Z² - I|| = $(norm(Z2 - I))")
````

### Test 2: Anti-commutation relations ({σᵢ, σⱼ} = 0 for i ≠ j)

````@example pauli_gns_construction
anticomm_XY = X_recon * Y_recon + Y_recon * X_recon
anticomm_YZ = Y_recon * Z_recon + Z_recon * Y_recon
anticomm_ZX = Z_recon * X_recon + X_recon * Z_recon

println("||{X,Y}|| = $(norm(anticomm_XY))")
println("||{Y,Z}|| = $(norm(anticomm_YZ))")
println("||{Z,X}|| = $(norm(anticomm_ZX))")
````

### Test 3: Commutation relations ([σᵢ, σⱼ] = 2iε_ijkσₖ)

````@example pauli_gns_construction
comm_XY = X_recon * Y_recon - Y_recon * X_recon
comm_YZ = Y_recon * Z_recon - Z_recon * Y_recon
comm_ZX = Z_recon * X_recon - X_recon * Z_recon

println("||[X,Y] - 2iZ|| = $(norm(comm_XY - 2im * Z_recon))")
println("||[Y,Z] - 2iX|| = $(norm(comm_YZ - 2im * X_recon))")
println("||[Z,X] - 2iY|| = $(norm(comm_ZX - 2im * Y_recon))")
````

## Understanding Unitary Freedom in GNS Reconstruction

**Important observation**: Looking at the reconstructed matrices above, you may notice that
the reconstructed X looks like the usual Pauli Y matrix (with an overall negative sign),
and the reconstructed Y looks like the usual Pauli X matrix. This might seem surprising,
but it's completely fine! Here's why:

1. **The Pauli algebra is still satisfied**: Even though `X_recon` ≈ -σ_y (negative Pauli-Y),
   it still satisfies all the required algebraic relations:
   - X² = I (since (-σ_y)² = σ_y² = I)
   - {X, Y} = 0 (anticommutation still holds)
   - [X, Y] = 2iZ (commutation relations preserved)

2. **Unitary equivalence**: The reconstructed matrices differ from standard Pauli matrices
   by a unitary transformation. This is fundamental to GNS construction - it only guarantees
   reconstruction **up to unitary equivalence**, not exact reconstruction of a specific
   matrix representation.

## Example 2: Random Pure State

Let's verify that GNS reconstruction works for an arbitrary pure state:

Create a random normalized pure state (not aligned with computational basis)

````@example pauli_gns_construction
using Random  # For reproducibility in random state example
Random.seed!(42)  # For reproducibility
ψ_random = normalize(randn(ComplexF64, 2))
ρ_random = ψ_random * ψ_random'

@show ψ_random;
nothing #hide
````

Build moment matrix and reconstruct

````@example pauli_gns_construction
H_random = zeros(ComplexF64, n, n)
for i in 1:n, j in 1:n
    product = NCTSSoS.neat_dot(basis[i], basis[j])
    H_random[i, j] = expval_pauli(product, ρ_random)
end

X_rand, Y_rand, Z_rand = reconstruct(H_random, vars, degree; atol=0.001)

@show rank(H_random, atol=1e-6);
@show size(X_rand);
nothing #hide
````

The reconstructed operators look very different from standard Pauli matrices:

````@example pauli_gns_construction
X_rand
````

````@example pauli_gns_construction
Y_rand
````

````@example pauli_gns_construction
Z_rand
````

The reconstructed operators are not in the familar form but still satisfy the Pauli algebra:

````@example pauli_gns_construction
@show norm(X_rand * X_rand - I);
@show norm(X_rand * Y_rand + Y_rand * X_rand);
@show norm((X_rand * Y_rand - Y_rand * X_rand) - 2im * Z_rand);
nothing #hide
````

Now, let's diagonalize Z_rand and use its eigenvectors to transform all three operators:

### Transforming to Eigenbasis of Z

The key insight: If we diagonalize Z_rand, the resulting unitary transformation
will bring X, Y, and Z into the standard Pauli matrix form (up to ordering and phases).

````@example pauli_gns_construction
F_rand = eigen(Z_rand)  # Diagonalize Z, not Y!

@show F_rand.values;
nothing #hide
````

Important: eigen() may return eigenvectors in any order. The standard Pauli Z matrix
is σ_z = diag(1, -1), so we want eigenvalue +1 first and -1 second.

````@example pauli_gns_construction
if real(F_rand.values[1]) < real(F_rand.values[2])
    # Need to swap eigenvectors to get correct ordering (eigenvalue +1 first)
    U_rand = F_rand.vectors[:, [2, 1]]
    @info "Reordered eigenvectors to match σ_z = diag(1, -1) convention"
else
    U_rand = F_rand.vectors
end
````

Apply the unitary transformation U_rand† · (operator) · U_rand to all three operators:

````@example pauli_gns_construction
Z_transformed = U_rand' * Z_rand * U_rand;
X_transformed = U_rand' * X_rand * U_rand;
Y_transformed = U_rand' * Y_rand * U_rand;
nothing #hide
````

Display the transformed operators:

````@example pauli_gns_construction
println("\nZ after transformation (should be diagonal):")
round.(Z_transformed, digits=4)
````

````@example pauli_gns_construction
println("\nX after transformation:")
round.(X_transformed, digits=4)
````

````@example pauli_gns_construction
println("\nY after transformation:")
round.(Y_transformed, digits=4)
````

## Second Unitary Transformation to Fix X and Y

The first transformation (diagonalizing Z) gave us a diagonal Z, but X and Y might still
differ from standard Pauli matrices by permutations or sign flips. We need a second
unitary transformation that:
1. Preserves the diagonal form of Z (must be diagonal in the Z eigenbasis)
2. Brings X and Y to their standard forms

The most general unitary that preserves a diagonal matrix is another diagonal unitary
(phase matrix) or a permutation followed by phases. Let's construct this explicitly:

## Systematic Approach: Determine U₂ from X_transformed

**Goal**: Find U₂ such that:
1. U₂† · X_transformed · U₂ = σ_x = `[0 1; 1 0]`
2. U₂† · Z_transformed · U₂ = Z_transformed (preserve Z's diagonal form)
3. U₂† · Y_transformed · U₂ = σ_y = `[0 -im; im 0]` (automatically satisfied if 1 & 2 hold)

**Key insight**: Since Z is diagonal, any unitary that preserves Z must be diagonal:
U₂ = diag(e^(iα), e^(iβ))

If `X_transformed` = `[0  a; b  0]`, we want:
U₂† · `X_transformed` · U₂ = `[0  e^(-iα)·a·e^(iβ); e^(-iβ)·b·e^(iα)  0]` = `[0  1; 1  0]`

This requires: e^(i(β-α))·a = 1, so β - α = -arg(a)
And for Hermiticity (b = a*): e^(i(α-β))·a* = 1, which is automatically satisfied

We can choose α = 0, then β = -arg(a)

Let's compute this systematically:

````@example pauli_gns_construction
σ_x = ComplexF64[0 1; 1 0]
σ_y = ComplexF64[0 -im; im 0]
σ_z = ComplexF64[1 0; 0 -1]
if abs(Z_transformed[1,1] - 1) < 0.1 || abs(Z_transformed[1,1] + 1) < 0.1  # Z is diagonal
    # Extract the off-diagonal elements of X_transformed
    # X should have form [0  a; a*  0] for Hermitian matrix
    a = X_transformed[1, 2]

    # We want U₂† [0 a; a* 0] U₂ = [0 1; 1 0]
    # With U₂ = diag(e^(iα), e^(iβ)), this gives:
    # [0  e^(-iα)ae^(iβ); e^(-iβ)a*e^(iα)  0] = [0 1; 1 0]
    # So we need: e^(i(β-α)) = 1/a, or β - α = -arg(a)

    # Choose α = 0, then β = -arg(a)
    α = 0.0
    β = -angle(a)

    # Construct U₂
    U2_systematic = diagm([exp(im*α), exp(im*β)])

    println("\nPhases determined: α = $α, β = $β")
    println("This should make X_transformed[1,2] real and positive: $(exp(im*(β-α)) * a)")

    # Apply the systematic transformation
    Z_corrected = U2_systematic' * Z_transformed * U2_systematic
    X_corrected = U2_systematic' * X_transformed * U2_systematic
    Y_corrected = U2_systematic' * Y_transformed * U2_systematic

    println("\n=== After Systematic Phase Correction ===")

    println("\nZ (should remain unchanged since U₂ is diagonal):")
    println(round.(Z_corrected, digits=6))
    println("Verification: ||Z_corrected - Z_transformed|| = $(norm(Z_corrected - Z_transformed))")

    println("\nX (should now match σ_x):")
    println(round.(X_corrected, digits=6))

    println("\nY (should now match σ_y):")
    println(round.(Y_corrected, digits=6))

    println("\n=== Final Verification ===")
    println("Distance from standard Pauli matrices:")
    println("||Z_corrected - σ_z|| = $(norm(Z_corrected - σ_z))")
    println("||X_corrected - σ_x|| = $(norm(X_corrected - σ_x))")
    println("||Y_corrected - σ_y|| = $(norm(Y_corrected - σ_y))")

    # If still not perfect, check for permutations
    if norm(X_corrected - σ_x) > 0.01
        println("\nChecking if X and Y are swapped:")
        println("||X_corrected - σ_y|| = $(norm(X_corrected - σ_y))")
        println("||Y_corrected - σ_x|| = $(norm(Y_corrected - σ_x))")
    end

    # Verify Pauli algebra is still satisfied
    println("\n=== Verify Pauli Algebra Still Holds ===")
    println("||X² - I|| = $(norm(X_corrected * X_corrected - I))")
    println("||Y² - I|| = $(norm(Y_corrected * Y_corrected - I))")
    println("||Z² - I|| = $(norm(Z_corrected * Z_corrected - I))")
    println("||[X,Y] - 2iZ|| = $(norm((X_corrected * Y_corrected - Y_corrected * X_corrected) - 2im * Z_corrected))")
end
````

## Example 3: Mixed State with Higher Rank

Now let's see what happens with a mixed state that is not a pure state.
Mixed state: ½|0⟩⟨0| + ½|+⟩⟨+| where |+⟩ = (|0⟩ + |1⟩)/√2

````@example pauli_gns_construction
plus_state = normalize(ComplexF64[1; 1])
ρ_mixed = 0.5 * (zero_state * zero_state') + 0.5 * (plus_state * plus_state')
````

Build moment matrix and reconstruct

````@example pauli_gns_construction
H_mixed = zeros(ComplexF64, n, n)
for i in 1:n, j in 1:n
    product = NCTSSoS.neat_dot(basis[i], basis[j])
    H_mixed[i, j] = expval_pauli(product, ρ_mixed)
end

X_mixed, Y_mixed, Z_mixed = reconstruct(H_mixed, vars, degree; atol=0.001)

@show rank(H_mixed, atol=1e-6);
@show size(X_mixed);
nothing #hide
````

For a mixed state, the rank can be higher! The reconstructed operators now act on a
**direct sum** of smaller Hilbert spaces. But remarkably, they still satisfy the
Pauli algebra:

### Complete Pauli Algebra Verification

````@example pauli_gns_construction
dim_mixed = size(X_mixed, 1)
````

**Test 1: Squares should equal identity (σᵢ² = I)**

````@example pauli_gns_construction
println("||X² - I|| = $(norm(X_mixed * X_mixed - I(dim_mixed)))")
println("||Y² - I|| = $(norm(Y_mixed * Y_mixed - I(dim_mixed)))")
println("||Z² - I|| = $(norm(Z_mixed * Z_mixed - I(dim_mixed)))")
````

**Test 2: Anti-commutation relations ({σᵢ, σⱼ} = 0 for i ≠ j)**

````@example pauli_gns_construction
anticomm_XY_mixed = X_mixed * Y_mixed + Y_mixed * X_mixed
anticomm_YZ_mixed = Y_mixed * Z_mixed + Z_mixed * Y_mixed
anticomm_ZX_mixed = Z_mixed * X_mixed + X_mixed * Z_mixed

println("||{X,Y}|| = $(norm(anticomm_XY_mixed))")
println("||{Y,Z}|| = $(norm(anticomm_YZ_mixed))")
println("||{Z,X}|| = $(norm(anticomm_ZX_mixed))")
````

**Test 3: Commutation relations ([σᵢ, σⱼ] = 2iε_ijkσₖ)**

````@example pauli_gns_construction
comm_XY_mixed = X_mixed * Y_mixed - Y_mixed * X_mixed
comm_YZ_mixed = Y_mixed * Z_mixed - Z_mixed * Y_mixed
comm_ZX_mixed = Z_mixed * X_mixed - X_mixed * Z_mixed

println("||[X,Y] - 2iZ|| = $(norm(comm_XY_mixed - 2im * Z_mixed))")
println("||[Y,Z] - 2iX|| = $(norm(comm_YZ_mixed - 2im * X_mixed))")
println("||[Z,X] - 2iY|| = $(norm(comm_ZX_mixed - 2im * Y_mixed))")
````

This demonstrates that GNS reconstruction can handle both pure and mixed states,
automatically determining the appropriate Hilbert space dimension!

## Understanding the Flat Extension Property

A crucial condition for successful GNS reconstruction is the **flat extension property**.
This requires that the rank of the full moment matrix equals the rank of its principal
submatrix (the Hankel block).

The `reconstruct` function automatically checks this condition and will issue a warning
if it's not satisfied. When the flatness condition holds, we can be confident that our
reconstruction is valid and complete.

## Summary

✓ GNS reconstruction successfully extracts operator representations from moments

✓ Reconstructed operators satisfy the correct algebraic relations

✓ Works for both pure states (low rank) and mixed states (higher rank)

✓ Operators are unique up to unitary transformations

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

