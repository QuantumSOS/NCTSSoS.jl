# [Overview of Optimizer](@id overview-of-optimizers)
A brief overview of different Semideinite Program optimizers are provided below.

## [Mosek](@id mosek)
Mosek [mosek](@cite) is a high-performance commercial solver for convex
optimization problems, including linear, conic, and semidefinite programming. It
is known for its speed and reliability, making it a popular choice for academic
and industrial applications. While it is a closed-source product, it offers free
licenses for academic use.

```julia
using MosekTools
SolverConfig(optimizer=Mosek.Optimizer) # using Mosek as the optimizer
```

## [Clarabel](@id clarabel)
Clarabel [Clarabel_2024](@cite) is a modern, open-source interior-point solver
for convex optimization. It is written in Rust and Julia and is designed to be
fast, reliable, and easy to use. It has a growing community of users and
developers and is a promising alternative to commercial solvers like Mosek. For
performance, you may need to install closed-source [linear system
solvers](https://clarabel.org/stable/user_guide_linsolvers/).

```julia
using Clarabel
SolverConfig(optimizer=Clarabel.Optimizer) # using Clarabel as the optimizer
```

## [COSMO](@id cosmo)
COSMO [Garstka_2021](@cite) (Conic Operator Splitting Method) is a first-order
solver for convex conic optimization problems. It is based on the alternating
direction method of multipliers (ADMM) and is designed to be fast and scalable,
especially for large-scale problems. However, it may not be as accurate as
interior-point methods like Mosek and Clarabel.

```julia
using COSMO
SolverConfig(optimizer=COSMO.Optimizer) # using COSMO as the optimizer
```