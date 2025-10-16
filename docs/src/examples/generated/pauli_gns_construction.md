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

In `NCTSSoS.jl`, we represent Pauli operators as non-commuting polynomial variables:

````julia
using NCTSSoS
using NCTSSoS.FastPolynomials
using LinearAlgebra
using LinearAlgebra: tr
````

Declare non-commuting variables for Pauli operators

````julia
@ncpolyvar x y z;
````

These variables x, y, z will represent σₓ, σᵧ, σ_z respectively

````julia
vars = [x, y, z];
````

## Step 2: Choose a Quantum State

Define quantum states for testing

````julia
zero_state = ComplexF64[1; 0];    # |0⟩
````

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
        # For Hermitian operators, b_i† = b_i
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
 [0.0 + 0.0im 3.363354988591067e-48 + 1.0000000000000004im; -6.119919185539761e-49 - 1.0000000000000004im -8.310098917615953e-33 + 5.098721737023209e-51im]
 [0.0 + 0.0im 1.0000000000000004 - 3.363354988591067e-48im; 1.0000000000000004 - 6.119919185539761e-49im -7.67288159674941e-50 - 6.1530573799022466e-33im]
 [1.0000000000000002 + 0.0im -1.9243339383232715e-50 - 5.7353166945488984e-33im; -5.289187967971371e-50 + 4.651031097596431e-33im -1.0000000000000009 + 2.4324615273979e-48im]
````

````julia
println("Reconstructed Pauli operators:")
println("σₓ (reconstructed):")
round.(X_recon, digits=6)
````

````
2×2 Matrix{ComplexF64}:
  0.0+0.0im   0.0+1.0im
 -0.0-1.0im  -0.0+0.0im
````

````julia
println("σᵧ (reconstructed):")
round.(Y_recon, digits=6)
````

````
2×2 Matrix{ComplexF64}:
 0.0+0.0im   1.0-0.0im
 1.0-0.0im  -0.0-0.0im
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

````julia
println("=== PAULI ALGEBRA VERIFICATION ===")
````

````
=== PAULI ALGEBRA VERIFICATION ===

````

Test 1: Squares should equal identity

````julia
println("1. Testing X² = Y² = Z² = I:")
X2 = X_recon * X_recon
Y2 = Y_recon * Y_recon
Z2 = Z_recon * Z_recon

println("   ||X² - I|| = $(norm(X2 - I(2)))")
println("   ||Y² - I|| = $(norm(Y2 - I(2)))")
println("   ||Z² - I|| = $(norm(Z2 - I(2)))")
````

````
1. Testing X² = Y² = Z² = I:
   ||X² - I|| = 1.2560739669470201e-15
   ||Y² - I|| = 1.2560739669470201e-15
   ||Z² - I|| = 1.831026719408895e-15

````

Test 2: Anti-commutators should be zero

````julia
println("2. Testing anti-commutation {σ_i, σ_j} = 0 for i ≠ j:")
anticomm_XY = X_recon * Y_recon + Y_recon * X_recon
anticomm_YZ = Y_recon * Z_recon + Z_recon * Y_recon
anticomm_ZX = Z_recon * X_recon + X_recon * Z_recon

println("   ||{X,Y}|| = $(norm(anticomm_XY))")
println("   ||{Y,Z}|| = $(norm(anticomm_YZ))")
println("   ||{Z,X}|| = $(norm(anticomm_ZX))")
````

````
2. Testing anti-commutation {σ_i, σ_j} = 0 for i ≠ j:
   ||{X,Y}|| = 1.4623122726759252e-32
   ||{Y,Z}|| = 9.42055475210265e-16
   ||{Z,X}|| = 9.42055475210265e-16

````

Test 3: Commutation relations

````julia
println("3. Testing commutation [σ_i, σ_j] = 2iε_ijkσ_k:")
comm_XY = X_recon * Y_recon - Y_recon * X_recon
comm_YZ = Y_recon * Z_recon - Z_recon * Y_recon
comm_ZX = Z_recon * X_recon - X_recon * Z_recon

println("   ||[X,Y] - 2iZ|| = $(norm(comm_XY - 2im * Z_recon))")
println("   ||[Y,Z] - 2iX|| = $(norm(comm_YZ - 2im * X_recon))")
println("   ||[Z,X] - 2iY|| = $(norm(comm_ZX - 2im * Y_recon))")
````

````
3. Testing commutation [σ_i, σ_j] = 2iε_ijkσ_k:
   ||[X,Y] - 2iZ|| = 1.3322676295501878e-15
   ||[Y,Z] - 2iX|| = 1.2560739669470201e-15
   ||[Z,X] - 2iY|| = 1.2560739669470201e-15

````

## Understanding the Flat Extension Property

A crucial condition for successful GNS reconstruction is the **flat extension property**.
This requires that the rank of the full moment matrix equals the rank of its principal
submatrix (the Hankel block).

The `reconstruct` function automatically checks this condition:

````julia
println("=== FLAT EXTENSION CHECK ===")
println("(This check is performed automatically in the reconstruct function)")
````

````
=== FLAT EXTENSION CHECK ===
(This check is performed automatically in the reconstruct function)

````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

