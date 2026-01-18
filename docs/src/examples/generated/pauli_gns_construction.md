# GNS Construction for Operator Reconstruction

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

## The `reconstruct` Function

NCTSSoS provides a `reconstruct` function that performs GNS construction. Given a moment
matrix `H` and a registry of variables, it returns matrix representations of the variables.

The function signature is:
```julia
reconstruct(H::Matrix, registry::VariableRegistry, degree::Int; atol=1e-3)
```

## Example: Simple Non-Commutative Variables

Let's demonstrate with a simple example using two non-commuting variables.

````julia
using NCTSSoS
using LinearAlgebra
````

Create two non-commutative variables

````julia
registry, (x,) = create_noncommutative_variables([("x", 1:2)])
````

````
(VariableRegistry with 2 variables: x₁, x₂, (NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}[UInt8[0x05], UInt8[0x09]],))
````

Generate a basis of monomials up to degree 2

````julia
basis = get_ncbasis(registry, 2)
println("Basis monomials (degree ≤ 2): ", length(basis))
````

````
Basis monomials (degree ≤ 2): 7

````

For GNS reconstruction, we need a moment matrix that encodes expectation values
of all products of basis elements. The matrix entry H[i,j] = ⟨basis[i]† · basis[j]⟩.

Here we'll use a simple positive definite moment matrix as an example:

````julia
n = length(basis)
H = zeros(Float64, n, n)
````

````
7×7 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
````

Create a simple valid moment matrix (identity-like with some correlations)

````julia
for i in 1:n
    H[i,i] = 1.0
end
````

Add small off-diagonal terms for correlation structure

````julia
H[1,2] = H[2,1] = 0.3
H[1,3] = H[3,1] = 0.3
H[2,3] = H[3,2] = 0.2
````

````
0.2
````

The moment matrix should be positive semidefinite

````julia
@assert isposdef(Hermitian(H)) "Moment matrix must be positive definite"
````

Perform GNS reconstruction

````julia
matrices = reconstruct(H, registry, 2; atol=0.001)
````

````
Dict{UInt8, Matrix{Float64}} with 2 entries:
  0x05 => [0.0 0.0 0.0; -0.0 0.0 0.0; -0.0 0.0 0.0]
  0x09 => [0.0 0.0 0.0; -0.0 0.0 0.0; -0.0 0.0 0.0]
````

Access the reconstructed matrix for each variable

````julia
var_indices = collect(NCTSSoS.indices(registry))
X1 = matrices[var_indices[1]]
X2 = matrices[var_indices[2]]

println("Reconstructed X₁ matrix:")
display(round.(X1, digits=4))

println("\nReconstructed X₂ matrix:")
display(round.(X2, digits=4))
````

````
Reconstructed X₁ matrix:

Reconstructed X₂ matrix:

````

## Verifying the Reconstruction

The reconstructed matrices should be consistent with the moment data.
Specifically, the GNS construction guarantees that these matrices satisfy
the algebraic relations encoded by the algebra type.

````julia
@show size(X1)
@show size(X2)
````

````
(3, 3)
````

## Key Properties of GNS Reconstruction

1. **Unitary Freedom**: Reconstructed matrices are unique only up to unitary transformations.
   Different moment matrices from the same state can yield unitarily equivalent matrices.

2. **Rank Determines Dimension**: The rank of the Hankel matrix determines the dimension
   of the reconstructed Hilbert space. Pure states give lower-rank matrices.

3. **Flatness Condition**: For valid reconstruction, the moment matrix should satisfy
   the "flat extension" property (rank should stabilize at lower degrees).

## Using GNS with Polynomial Optimization

The primary use of GNS reconstruction in NCTSSoS is to extract optimal solutions
from semidefinite programming relaxations. When solving polynomial optimization:

1. The SDP relaxation produces a moment matrix as part of its dual solution
2. `reconstruct` extracts matrix representations of variables
3. These matrices can be used to verify or extract the optimal point

See the optimization examples for practical applications of GNS reconstruction
in the context of polynomial optimization over non-commutative variables.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

