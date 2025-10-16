# # GNS Construction for Pauli Operator Reconstruction

# The Gelfand-Naimark-Segal (GNS) construction is a fundamental mathematical tool
# in quantum mechanics that allows us to represent abstract quantum states and operators
# as concrete matrices acting on a Hilbert space. For physicists, this provides a systematic
# way to reconstruct operator representations from expectation values (moments).

# ## Background: From Expectation Values to Matrix Representations

# In quantum mechanics, we often know the expectation values of operators in a given state:
# ```math
# \langle A  \rangle = \text{Tr}(\rho A)
# ```
# where ``\rho`` is the density matrix and ``A`` is an operator.

# The GNS construction answers the question: **Can we reconstruct the actual matrices
# representing operators from just these expectation values?**

# The key insight is that the collection of all expectation values defines a **moment matrix**
# (also called a Hankel matrix in the context of polynomial optimization):
# ```math
# H_{ij} = \langle b_i^\dagger b_j \rangle
# ```
# where ``\{b_i\}`` is a basis of operators (monomials in our variables).

# ## Pauli Operators: A Simple Quantum System

# Let's start with the simplest non-trivial quantum system: a single qubit with Pauli operators.
# The Pauli matrices are:
# - ``\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0  \end{pmatrix}``
# - ``\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}``
# - ``\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}```

# These satisfy the fundamental Pauli algebra:
# 1. ``\sigma_x^2 = \sigma_y^2 = \sigma_z^2 = I`` (squares equal identity)
# 2. ``\{\sigma_i, \sigma_j\} = 0`` for ``i \neq j`` (anti-commutation)
# 3. ``[\sigma_i, \sigma_j] = 2i\epsilon_{ijk}\sigma_k`` (commutation relations)

# ## Step 1: Define Non-commuting Variables

# In `NCTSSoS.jl`, we represent Pauli operators as non-commuting polynomial variables:

using NCTSSoS
using NCTSSoS.FastPolynomials
using LinearAlgebra
using LinearAlgebra: tr

# Declare non-commuting variables for Pauli operators
@ncpolyvar x y z;

# These variables x, y, z will represent σₓ, σᵧ, σ_z respectively
vars = [x, y, z];

# ## Step 2: Choose a Quantum State

# Define quantum states for testing
zero_state = ComplexF64[1; 0];    # |0⟩

# For clear reconstruction, use a pure state
ρ =  zero_state * zero_state'

# ## Step 3: Compute Expectation Values

# We need a function to compute the expectation value of any monomial in our variables:

function expval_pauli(mono::Monomial, ρ::Matrix)
    if isone(mono)
        return 1.0 + 0.0im  # ⟨I⟩ = 1 for normalized states
    end

    ## Start with identity matrix
    mat = Matrix{ComplexF64}(I, 2, 2)

    ## Multiply the appropriate Pauli matrices
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

        ## Raise to the appropriate power
        for _ in 1:exp
            mat = mat * pauli_mat
        end
    end

    ## Compute Tr(ρ * mat)
    return tr(ρ * mat)
end

# ## Step 4: Build the Moment (Hankel) Matrix

# The moment matrix encodes all expectation values of products of our basis operators:

# Choose the degree of our polynomial basis
degree = 4;

# Generate the basis of monomials up to the specified degree
basis = NCTSSoS.get_basis(vars, degree)

println("Basis operators (monomials):")
for (i, b) in enumerate(basis)
    println("$i: $b")
end

# Build the moment matrix H where H[i,j] = ⟨b_i† * b_j⟩
n = length(basis)
H = zeros(ComplexF64, n, n)

for i in 1:n
    for j in 1:n
        ## For Hermitian operators, b_i† = b_i
        mono_i = basis[i]
        mono_j = basis[j]

        ## Compute the product monomial
        product = NCTSSoS.neat_dot(mono_i, mono_j)

        ## Get the expectation value
        H[i, j] = expval_pauli(product, ρ)
    end
end

# Verify that H is Hermitian (as required by quantum mechanics)
println("Moment matrix H is Hermitian: ", H ≈ H')

# ## Step 5: GNS Reconstruction

# Now we use the `reconstruct` function to perform the GNS construction and obtain
# concrete matrix representations of our abstract operators:

X_recon, Y_recon, Z_recon = reconstruct(H, vars, degree; atol=0.001)

#
println("Reconstructed Pauli operators:")
println("σₓ (reconstructed):")
round.(X_recon, digits=6)

#
println("σᵧ (reconstructed):")
round.(Y_recon, digits=6)

#
println("σ_z (reconstructed):")
round.(Z_recon, digits=6)

# ## Step 6: Verify Pauli Algebra

# The true test of our reconstruction is whether the recovered operators satisfy
# the Pauli algebra relations:

println("=== PAULI ALGEBRA VERIFICATION ===")

# Test 1: Squares should equal identity
println("1. Testing X² = Y² = Z² = I:")
X2 = X_recon * X_recon
Y2 = Y_recon * Y_recon
Z2 = Z_recon * Z_recon

println("   ||X² - I|| = $(norm(X2 - I(2)))")
println("   ||Y² - I|| = $(norm(Y2 - I(2)))")
println("   ||Z² - I|| = $(norm(Z2 - I(2)))")

# Test 2: Anti-commutators should be zero
println("2. Testing anti-commutation {σ_i, σ_j} = 0 for i ≠ j:")
anticomm_XY = X_recon * Y_recon + Y_recon * X_recon
anticomm_YZ = Y_recon * Z_recon + Z_recon * Y_recon
anticomm_ZX = Z_recon * X_recon + X_recon * Z_recon

println("   ||{X,Y}|| = $(norm(anticomm_XY))")
println("   ||{Y,Z}|| = $(norm(anticomm_YZ))")
println("   ||{Z,X}|| = $(norm(anticomm_ZX))")

# Test 3: Commutation relations
println("3. Testing commutation [σ_i, σ_j] = 2iε_ijkσ_k:")
comm_XY = X_recon * Y_recon - Y_recon * X_recon
comm_YZ = Y_recon * Z_recon - Z_recon * Y_recon
comm_ZX = Z_recon * X_recon - X_recon * Z_recon

println("   ||[X,Y] - 2iZ|| = $(norm(comm_XY - 2im * Z_recon))")
println("   ||[Y,Z] - 2iX|| = $(norm(comm_YZ - 2im * X_recon))")
println("   ||[Z,X] - 2iY|| = $(norm(comm_ZX - 2im * Y_recon))")

# ## Understanding the Flat Extension Property

# A crucial condition for successful GNS reconstruction is the **flat extension property**.
# This requires that the rank of the full moment matrix equals the rank of its principal
# submatrix (the Hankel block).

# The `reconstruct` function automatically checks this condition:
println("=== FLAT EXTENSION CHECK ===")
println("(This check is performed automatically in the reconstruct function)")