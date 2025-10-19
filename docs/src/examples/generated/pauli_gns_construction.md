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

````julia
using NCTSSoS
using NCTSSoS.FastPolynomials
using LinearAlgebra
using LinearAlgebra: tr
````

In `NCTSSoS.jl`, we represent Pauli operators as non-commuting polynomial variables:
Declare non-commuting variables for Pauli operators

````julia
@ncpolyvar x y z;
````

These variables x, y, z will represent σₓ, σᵧ, σ_z respectively

````julia
vars = [x, y, z];
````

## Step 2: Choose a Quantum State

````julia
zero_state = ComplexF64[1; 0];    # |0⟩
one_state = ComplexF64[0; 1];
````

Define quantum states for testing
For clear reconstruction, use a pure state

````julia
ρ =  zero_state * zero_state'
````

````
2×2 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im
````

## Step 3: Compute Expectation Values

We need a function to compute the expectation value of any monomial in our variables:

````julia
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

````
expval_pauli (generic function with 1 method)
````

## Step 4: Build the Moment (Hankel) Matrix

The moment matrix encodes all expectation values of products of our basis operators:

Choose the degree of our polynomial basis

````julia
degree = 4;
````

Generate the basis of monomials up to the specified degree

````julia
basis = NCTSSoS.get_basis(vars, degree)

println("Basis operators (monomials):")
for (i, b) in enumerate(basis)
    println("$i: $b")
end
````

````
Basis operators (monomials):
1: 1
2: x¹
3: y¹
4: z¹
5: x¹y¹
6: x¹z¹
7: x²
8: y¹x¹
9: y¹z¹
10: y²
11: z¹x¹
12: z¹y¹
13: z²
14: x¹y¹x¹
15: x¹y¹z¹
16: x¹y²
17: x¹z¹x¹
18: x¹z¹y¹
19: x¹z²
20: x²y¹
21: x²z¹
22: x³
23: y¹x¹y¹
24: y¹x¹z¹
25: y¹x²
26: y¹z¹x¹
27: y¹z¹y¹
28: y¹z²
29: y²x¹
30: y²z¹
31: y³
32: z¹x¹y¹
33: z¹x¹z¹
34: z¹x²
35: z¹y¹x¹
36: z¹y¹z¹
37: z¹y²
38: z²x¹
39: z²y¹
40: z³
41: x¹y¹x¹y¹
42: x¹y¹x¹z¹
43: x¹y¹x²
44: x¹y¹z¹x¹
45: x¹y¹z¹y¹
46: x¹y¹z²
47: x¹y²x¹
48: x¹y²z¹
49: x¹y³
50: x¹z¹x¹y¹
51: x¹z¹x¹z¹
52: x¹z¹x²
53: x¹z¹y¹x¹
54: x¹z¹y¹z¹
55: x¹z¹y²
56: x¹z²x¹
57: x¹z²y¹
58: x¹z³
59: x²y¹x¹
60: x²y¹z¹
61: x²y²
62: x²z¹x¹
63: x²z¹y¹
64: x²z²
65: x³y¹
66: x³z¹
67: x⁴
68: y¹x¹y¹x¹
69: y¹x¹y¹z¹
70: y¹x¹y²
71: y¹x¹z¹x¹
72: y¹x¹z¹y¹
73: y¹x¹z²
74: y¹x²y¹
75: y¹x²z¹
76: y¹x³
77: y¹z¹x¹y¹
78: y¹z¹x¹z¹
79: y¹z¹x²
80: y¹z¹y¹x¹
81: y¹z¹y¹z¹
82: y¹z¹y²
83: y¹z²x¹
84: y¹z²y¹
85: y¹z³
86: y²x¹y¹
87: y²x¹z¹
88: y²x²
89: y²z¹x¹
90: y²z¹y¹
91: y²z²
92: y³x¹
93: y³z¹
94: y⁴
95: z¹x¹y¹x¹
96: z¹x¹y¹z¹
97: z¹x¹y²
98: z¹x¹z¹x¹
99: z¹x¹z¹y¹
100: z¹x¹z²
101: z¹x²y¹
102: z¹x²z¹
103: z¹x³
104: z¹y¹x¹y¹
105: z¹y¹x¹z¹
106: z¹y¹x²
107: z¹y¹z¹x¹
108: z¹y¹z¹y¹
109: z¹y¹z²
110: z¹y²x¹
111: z¹y²z¹
112: z¹y³
113: z²x¹y¹
114: z²x¹z¹
115: z²x²
116: z²y¹x¹
117: z²y¹z¹
118: z²y²
119: z³x¹
120: z³y¹
121: z⁴

````

Build the moment matrix H where H[i,j] = ⟨b_i† * b_j⟩

````julia
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

````julia
println("Moment matrix H is Hermitian: ", H ≈ H')
````

````
Moment matrix H is Hermitian: true

````

## Step 5: GNS Reconstruction

Now we use the `reconstruct` function to perform the GNS construction and obtain
concrete matrix representations of our abstract operators:

````julia
X_recon, Y_recon, Z_recon = reconstruct(H, vars, degree; atol=0.001)
````

````
3-element Vector{Matrix{ComplexF64}}:
 [0.0 + 0.0im 3.574334513209998e-48 + 1.0000000000000004im; 2.753963633492892e-48 - 1.0000000000000004im -8.310098917615953e-33 + 6.366149832128084e-51im]
 [0.0 + 0.0im 1.0000000000000004 - 3.574334513209998e-48im; 1.0000000000000004 + 2.753963633492892e-48im -2.347184397671815e-50 - 6.1530573799022466e-33im]
 [1.0000000000000002 + 0.0im -2.353834293849223e-50 - 5.7353166945488984e-33im; -2.51012310344404e-50 + 4.651031097596431e-33im -1.0000000000000009 + 2.3381282238640562e-48im]
````

````julia
println("Reconstructed Pauli operators:")
println("σₓ (reconstructed):")
round.(X_recon, digits=6)
````

````
2×2 Matrix{ComplexF64}:
 0.0+0.0im   0.0+1.0im
 0.0-1.0im  -0.0+0.0im
````

````julia
println("σᵧ (reconstructed):")
round.(Y_recon, digits=6)
````

````
2×2 Matrix{ComplexF64}:
 0.0+0.0im   1.0-0.0im
 1.0+0.0im  -0.0-0.0im
````

````julia
println("σ_z (reconstructed):")
round.(Z_recon, digits=6)
````

````
2×2 Matrix{ComplexF64}:
  1.0+0.0im  -0.0-0.0im
 -0.0+0.0im  -1.0+0.0im
````

## Step 6: Verify Pauli Algebra

The true test of our reconstruction is whether the recovered operators satisfy
the Pauli algebra relations:

### Test 1: Squares should equal identity (X² = Y² = Z² = I)

````julia
X2 = X_recon * X_recon
Y2 = Y_recon * Y_recon
Z2 = Z_recon * Z_recon

println("   ||X² - I|| = $(norm(X2 - I))")
println("   ||Y² - I|| = $(norm(Y2 - I))")
println("   ||Z² - I|| = $(norm(Z2 - I))")
````

````
   ||X² - I|| = 1.2560739669470201e-15
   ||Y² - I|| = 1.2560739669470201e-15
   ||Z² - I|| = 1.831026719408895e-15

````

### Test 2: Anti-commutation relations ({σᵢ, σⱼ} = 0 for i ≠ j)

````julia
anticomm_XY = X_recon * Y_recon + Y_recon * X_recon
anticomm_YZ = Y_recon * Z_recon + Z_recon * Y_recon
anticomm_ZX = Z_recon * X_recon + X_recon * Z_recon

println("||{X,Y}|| = $(norm(anticomm_XY))")
println("||{Y,Z}|| = $(norm(anticomm_YZ))")
println("||{Z,X}|| = $(norm(anticomm_ZX))")
````

````
||{X,Y}|| = 1.4623122726759252e-32
||{Y,Z}|| = 9.42055475210265e-16
||{Z,X}|| = 9.42055475210265e-16

````

### Test 3: Commutation relations ([σᵢ, σⱼ] = 2iε_ijkσₖ)

````julia
comm_XY = X_recon * Y_recon - Y_recon * X_recon
comm_YZ = Y_recon * Z_recon - Z_recon * Y_recon
comm_ZX = Z_recon * X_recon - X_recon * Z_recon

println("||[X,Y] - 2iZ|| = $(norm(comm_XY - 2im * Z_recon))")
println("||[Y,Z] - 2iX|| = $(norm(comm_YZ - 2im * X_recon))")
println("||[Z,X] - 2iY|| = $(norm(comm_ZX - 2im * Y_recon))")
````

````
||[X,Y] - 2iZ|| = 1.3322676295501878e-15
||[Y,Z] - 2iX|| = 1.2560739669470201e-15
||[Z,X] - 2iY|| = 1.2560739669470201e-15

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

````julia
using Random  # For reproducibility in random state example
Random.seed!(42)  # For reproducibility
ψ_random = normalize(randn(ComplexF64, 2))
ρ_random = ψ_random * ψ_random'

@show ψ_random;
````

````
ψ_random = ComplexF64[0.48005951599229074 - 0.5357791469433075im, -0.532085869997182 - 0.44650665589136673im]

````

Build moment matrix and reconstruct

````julia
H_random = zeros(ComplexF64, n, n)
for i in 1:n, j in 1:n
    product = NCTSSoS.neat_dot(basis[i], basis[j])
    H_random[i, j] = expval_pauli(product, ρ_random)
end

X_rand, Y_rand, Z_rand = reconstruct(H_random, vars, degree; atol=0.001)

@show rank(H_random, atol=1e-6);
@show size(X_rand);
````

````
Rank of full Hankel matrix H: 2 (using atol=0.001)
Rank of hankel_block (degree 3): 2 (using atol=0.001)
✓ Flatness condition satisfied: rank(H) = rank(hankel_block) = 2
GNS reconstruction: keeping 2 singular values > 0.001, reconstructed matrices will be 2×2
Variable x: constructed (2, 2) matrix representation
Variable y: constructed (2, 2) matrix representation
Variable z: constructed (2, 2) matrix representation
rank(H_random, atol = 1.0e-6) = 2
size(X_rand) = (2, 2)

````

The reconstructed operators look very different from standard Pauli matrices:

````julia
X_rand
````

````
2×2 Matrix{ComplexF64}:
 -0.0324079+1.23262e-19im   0.994615-0.0984462im
   0.994615+0.0984462im    0.0324079-3.68775e-18im
````

````julia
Y_rand
````

````
2×2 Matrix{ComplexF64}:
  -0.998861-1.45507e-19im  -0.0356829-0.0316907im
 -0.0356829+0.0316907im      0.998861+1.01064e-17im
````

````julia
Z_rand
````

````
2×2 Matrix{ComplexF64}:
 0.0350329-3.53414e-18im   -0.097307-0.994638im
 -0.097307+0.994638im     -0.0350329-8.01641e-18im
````

The reconstructed operators are not in the familar form but still satisfy the Pauli algebra:

````julia
@show norm(X_rand * X_rand - I);
@show norm(X_rand * Y_rand + Y_rand * X_rand);
@show norm((X_rand * Y_rand - Y_rand * X_rand) - 2im * Z_rand);
````

````
norm(X_rand * X_rand - I) = 4.9544522472129334e-17
norm(X_rand * Y_rand + Y_rand * X_rand) = 9.504511772915181e-16
norm((X_rand * Y_rand - Y_rand * X_rand) - (2im) * Z_rand) = 5.162078360867593e-16

````

Now, let's diagonalize Z_rand and use its eigenvectors to transform all three operators:

### Transforming to Eigenbasis of Z

The key insight: If we diagonalize Z_rand, the resulting unitary transformation
will bring X, Y, and Z into the standard Pauli matrix form (up to ordering and phases).

````julia
F_rand = eigen(Z_rand)  # Diagonalize Z, not Y!

@show F_rand.values;
````

````
F_rand.values = ComplexF64[-1.0 - 1.672081166228218e-17im, 1.0 + 6.247779567868549e-18im]

````

Important: eigen() may return eigenvectors in any order. The standard Pauli Z matrix
is σ_z = diag(1, -1), so we want eigenvalue +1 first and -1 second.

````julia
if real(F_rand.values[1]) < real(F_rand.values[2])
    # Need to swap eigenvectors to get correct ordering (eigenvalue +1 first)
    U_rand = F_rand.vectors[:, [2, 1]]
    @info "Reordered eigenvectors to match σ_z = diag(1, -1) convention"
else
    U_rand = F_rand.vectors
end
````

````
[ Info: Reordered eigenvectors to match σ_z = diag(1, -1) convention

````

Apply the unitary transformation U_rand† · (operator) · U_rand to all three operators:

````julia
Z_transformed = U_rand' * Z_rand * U_rand;
X_transformed = U_rand' * X_rand * U_rand;
Y_transformed = U_rand' * Y_rand * U_rand;
````

Display the transformed operators:

````julia
println("\nZ after transformation (should be diagonal):")
round.(Z_transformed, digits=4)
````

````
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0-0.0im
````

````julia
println("\nX after transformation:")
round.(X_transformed, digits=4)
````

````
2×2 Matrix{ComplexF64}:
    0.0+0.0im     0.9916-0.1296im
 0.9916+0.1296im     0.0-0.0im
````

````julia
println("\nY after transformation:")
round.(Y_transformed, digits=4)
````

````
2×2 Matrix{ComplexF64}:
     0.0+0.0im     -0.1296-0.9916im
 -0.1296+0.9916im      0.0+0.0im
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

````julia
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

````

Phases determined: α = 0.0, β = 0.12995470935070005
This should make X_transformed[1,2] real and positive: 0.9999999999999997 + 2.7755575615628914e-17im

=== After Systematic Phase Correction ===

Z (should remain unchanged since U₂ is diagonal):
ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 - 0.0im -1.0 - 0.0im]
Verification: ||Z_corrected - Z_transformed|| = 2.4360324397437158e-17

X (should now match σ_x):
ComplexF64[0.0 + 0.0im 1.0 + 0.0im; 1.0 - 0.0im 0.0 - 0.0im]

Y (should now match σ_y):
ComplexF64[0.0 + 0.0im 0.0 - 1.0im; 0.0 + 1.0im 0.0 + 0.0im]

=== Final Verification ===
Distance from standard Pauli matrices:
||Z_corrected - σ_z|| = 2.8052198857519395e-16
||X_corrected - σ_x|| = 4.922826142722048e-16
||Y_corrected - σ_y|| = 6.521822047694349e-16

=== Verify Pauli Algebra Still Holds ===
||X² - I|| = 9.499064443774055e-16
||Y² - I|| = 1.2948801525946436e-15
||Z² - I|| = 4.478422697720281e-16
||[X,Y] - 2iZ|| = 1.1309194191732256e-15

````

## Example 3: Mixed State with Higher Rank

Now let's see what happens with a mixed state that is not a pure state.
Mixed state: ½|0⟩⟨0| + ½|+⟩⟨+| where |+⟩ = (|0⟩ + |1⟩)/√2

````julia
plus_state = normalize(ComplexF64[1; 1])
ρ_mixed = 0.5 * (zero_state * zero_state') + 0.5 * (plus_state * plus_state')
````

````
2×2 Matrix{ComplexF64}:
 0.75+0.0im  0.25+0.0im
 0.25+0.0im  0.25+0.0im
````

Build moment matrix and reconstruct

````julia
H_mixed = zeros(ComplexF64, n, n)
for i in 1:n, j in 1:n
    product = NCTSSoS.neat_dot(basis[i], basis[j])
    H_mixed[i, j] = expval_pauli(product, ρ_mixed)
end

X_mixed, Y_mixed, Z_mixed = reconstruct(H_mixed, vars, degree; atol=0.001)

@show rank(H_mixed, atol=1e-6);
@show size(X_mixed);
````

````
Rank of full Hankel matrix H: 4 (using atol=0.001)
Rank of hankel_block (degree 3): 4 (using atol=0.001)
✓ Flatness condition satisfied: rank(H) = rank(hankel_block) = 4
GNS reconstruction: keeping 4 singular values > 0.001, reconstructed matrices will be 4×4
Variable x: constructed (4, 4) matrix representation
Variable y: constructed (4, 4) matrix representation
Variable z: constructed (4, 4) matrix representation
rank(H_mixed, atol = 1.0e-6) = 4
size(X_mixed) = (4, 4)

````

For a mixed state, the rank can be higher! The reconstructed operators now act on a
**direct sum** of smaller Hilbert spaces. But remarkably, they still satisfy the
Pauli algebra:

### Complete Pauli Algebra Verification

````julia
dim_mixed = size(X_mixed, 1)
````

````
4
````

**Test 1: Squares should equal identity (σᵢ² = I)**

````julia
println("||X² - I|| = $(norm(X_mixed * X_mixed - I(dim_mixed)))")
println("||Y² - I|| = $(norm(Y_mixed * Y_mixed - I(dim_mixed)))")
println("||Z² - I|| = $(norm(Z_mixed * Z_mixed - I(dim_mixed)))")
````

````
||X² - I|| = 2.4124410454036057e-15
||Y² - I|| = 2.7614386208727306e-15
||Z² - I|| = 2.6211266463523374e-15

````

**Test 2: Anti-commutation relations ({σᵢ, σⱼ} = 0 for i ≠ j)**

````julia
anticomm_XY_mixed = X_mixed * Y_mixed + Y_mixed * X_mixed
anticomm_YZ_mixed = Y_mixed * Z_mixed + Z_mixed * Y_mixed
anticomm_ZX_mixed = Z_mixed * X_mixed + X_mixed * Z_mixed

println("||{X,Y}|| = $(norm(anticomm_XY_mixed))")
println("||{Y,Z}|| = $(norm(anticomm_YZ_mixed))")
println("||{Z,X}|| = $(norm(anticomm_ZX_mixed))")
````

````
||{X,Y}|| = 1.81421523891841e-15
||{Y,Z}|| = 2.6164352461025136e-15
||{Z,X}|| = 2.5553316515572013e-15

````

**Test 3: Commutation relations ([σᵢ, σⱼ] = 2iε_ijkσₖ)**

````julia
comm_XY_mixed = X_mixed * Y_mixed - Y_mixed * X_mixed
comm_YZ_mixed = Y_mixed * Z_mixed - Z_mixed * Y_mixed
comm_ZX_mixed = Z_mixed * X_mixed - X_mixed * Z_mixed

println("||[X,Y] - 2iZ|| = $(norm(comm_XY_mixed - 2im * Z_mixed))")
println("||[Y,Z] - 2iX|| = $(norm(comm_YZ_mixed - 2im * X_mixed))")
println("||[Z,X] - 2iY|| = $(norm(comm_ZX_mixed - 2im * Y_mixed))")
````

````
||[X,Y] - 2iZ|| = 3.560277133665956e-15
||[Y,Z] - 2iX|| = 3.5287866252028285e-15
||[Z,X] - 2iY|| = 2.4546810564734144e-15

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

