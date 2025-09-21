```@meta
EditURL = "../literate/trace_poly.jl"
```

# Tracial Polynomial Optimization

## Toy Example
Let's learn how to do [tracial polynomial optimization](@ref
tracial-polynomial) from a toy example.

We use [`NCTSSoS.FastPolynomials.tr`](@ref) to declare a part of a term in
tracial polynomial.

```julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial
@ncpolyvar x[1:3]

p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * one(Monomial)
```

Polynomial Optimization declaration and solving interface is the same as regular
polynomial optimization.

```julia
spop = polyopt(p; is_projective=true, comm_gps=[x])

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)

result = cs_nctssos(spop, solver_config)

@assert isapprox(result.objective , -0.046717378455438933, atol = 1e-6)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=3)

result = cs_nctssos(spop, solver_config)

@assert isapprox(result.objective, -0.03124998978001017, atol = 1e-6)
```

The results matches within $10^{-6}$ absolute tolerance comparing to answer in
[klep2022Optimization](@cite)!

## Polynomial Bell Inequalities

Polynomial Bell inequalities provide a powerful framework for detecting quantum
entanglement and non-locality in bipartite quantum systems. These inequalities
impose constraints on the correlations that can be achieved by local hidden
variable models, and their violation serves as a signature of quantum mechanical
behavior. For maximally entangled bipartite states, such as Bell states, the
quantum correlations can exceed the classical bounds imposed by these polynomial
inequalities, demonstrating the non-local nature of quantum entanglement. The
following examples illustrate how tracial polynomial optimization can be used to
compute the maximum violation of specific Bell inequalities, revealing the
extent to which quantum mechanics transcends classical limitations.

```julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial

@ncpolyvar x[1:2] y[1:2]

p = -1.0 * tr(x[1] * y[1]) - 1.0 * tr(x[1] * y[2]) - 1.0 * tr(x[2] * y[1]) + 1.0 * tr(x[2] * y[2])

tpop = polyopt(p * one(Monomial); is_unipotent=true)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=1, ts_algo=MaximalElimination())

result = cs_nctssos(tpop, solver_config)

@assert isapprox(result.objective, -2.8284271157283083, atol = 1e-5)
```

Our computation matches with the theoretical prediction for maximally entangled
bipartite state with $10^{-6}$ absolute tolerance [klep2022Optimization](@cite)!

## Covariance of quantum correlation

As introduced in [Bell Inequalities example](@ref bell-inequalities), we may
also compute the covariance of quantum correlations while limiting the state to
maximally entangled bipartite state.

```julia
using NCTSSoS, MosekTools
using NCTSSoS.FastPolynomials:tr, Monomial

@ncpolyvar x[1:3] y[1:3]

cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))
tpop = polyopt(p * one(Monomial); is_unipotent=true)

solver_config = SolverConfig(; optimizer=Mosek.Optimizer, order=2)

result = cs_nctssos(tpop, solver_config)

@assert isapprox(result.objective,-5.0, atol = 1e-5)
```

Again, the result matches the theoretical prediction for maximally entangled
bipartite state with $10^{-6}$ absolute tolerance [klep2022Optimization](@cite)!

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

