# [Overview of Optimizer](@id overview-of-optimizers)
A brief overview of semidefinite program optimizers is provided below.
Choose the problem and relaxation first, then pass a solver to `solve`:

```julia
using NCTSSoS.Next
using JuMP

_, (_, _, z) = pauli_variables(1)
problem = Problem(z[1])
relaxation = Relaxation(order=1)
```

Each example below solves this problem. Its minimum is −1.
`MosekBackend` writes the solver form directly. For other solvers,
`JuMPBackend` accepts a JuMP optimizer or `optimizer_with_attributes(...)`.

## [Mosek](@id mosek)
Mosek [mosek](@cite) is a high-performance commercial solver for convex
optimization problems, including linear, conic, and semidefinite programming. It
is known for its speed and reliability, making it a popular choice for academic
and industrial applications. While it is a closed-source product, it offers free
licenses for academic use.

```julia
using NCTSSoS.Next
using Mosek
result = solve(problem, relaxation; solver=MosekBackend())
@assert isapprox(result.bound, -1; atol=1e-6)
```

## [Clarabel](@id clarabel)
Clarabel [Clarabel_2024](@cite) is a modern, open-source interior-point solver
for convex optimization. It is written in Rust and Julia and is designed to be
fast, reliable, and easy to use. It has a growing community of users and
developers and is a promising alternative to commercial solvers like Mosek. For
performance, you may need to install closed-source [linear system
solvers](https://clarabel.org/stable/user_guide_linsolvers/).

The example disables chordal decomposition because Clarabel 0.11.1 with
OrderedCollections 2.0.3 fails in dual PSD completion for sparse PSD constraints.

```julia
using NCTSSoS.Next
using Clarabel, JuMP
optimizer = optimizer_with_attributes(Clarabel.Optimizer,
    "chordal_decomposition_enable" => false)
result = solve(problem, relaxation; solver=JuMPBackend(optimizer))
@assert isapprox(result.bound, -1; atol=1e-6)
```

## [COSMO](@id cosmo)
COSMO [Garstka_2021](@cite) (Conic Operator Splitting Method) is a first-order
solver for convex conic optimization problems. It is based on the alternating
direction method of multipliers (ADMM) and is designed to be fast and scalable,
especially for large-scale problems. However, it may not be as accurate as
interior-point methods like Mosek and Clarabel.

```julia
using NCTSSoS.Next
using COSMO
result = solve(problem, relaxation; solver=JuMPBackend(COSMO.Optimizer))
@assert isapprox(result.bound, -1; atol=1e-3)
```

Check `result.status` before you use a result. If the solver does not return
a usable solution, `result.bound` is `nothing`. A floating-point SDP result
is subject to solver tolerances; it is not a proof in exact arithmetic.
