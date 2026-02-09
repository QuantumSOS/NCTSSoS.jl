# Tracial Polynomial Optimization

## Toy Example
Let's learn how to do [tracial polynomial optimization](@ref
tracial-polynomial) from a toy example.

We use `NCTSSoS.tr` to declare a part of a term in
tracial polynomial.

````julia
using NCTSSoS, MosekTools
using NCTSSoS: tr
````

We use MathOptInterface (MOI) to configure the SDP solver backend (here: Mosek)
and silence its output.

````julia
const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)
````

Create projector variables using the typed algebra system.

- `registry` stores the variable mapping and the algebra constraints.
- `x` is a length-3 vector of noncommutative projector variables satisfying
  `x[i]^2 = x[i]` (idempotency is enforced by the `ProjectorAlgebra` type).

````julia
registry, (x,) = create_projector_variables([("x", 1:3)])
````

`tr(·)` constructs a *tracial state polynomial* (a symbolic moment in a tracial
state).

Internally, NCTSSoS represents noncommutative words using `NormalMonomial{A,T}`
(where `A` is the algebra type and `T` is the index integer type). To pass a
tracial state polynomial to `polyopt`, we multiply by the identity monomial to
embed it as a noncommutative polynomial objective.

Identity monomial for converting StatePolynomial to NCStatePolynomial.

````julia
const ID_PROJECTOR = one(NormalMonomial{ProjectorAlgebra,UInt8})

p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * ID_PROJECTOR
````

Polynomial Optimization declaration and solving interface is the same as regular
polynomial optimization. No need for is_projective or comm_gps - the registry
encodes all algebra constraints!

````julia
spop = polyopt(p, registry)
````

The SDP relaxation is configured via `SolverConfig`.

`order` is the moment/SOS relaxation order: higher order typically yields a
tighter bound at the cost of a larger SDP.

We solve the relaxation with `cs_nctssos`; the returned `result.objective` is a
certified bound from the SDP relaxation.

````julia
solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2)

result = cs_nctssos(spop, solver_config)

@show result.objective

@assert isapprox(result.objective , -0.046717378455438933, atol = 1e-6)

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=3)

result = cs_nctssos(spop, solver_config)

@show result.objective

@assert isapprox(result.objective, -0.03124998978001017, atol = 1e-6)
````

The literature values are $-0.046717378455438933$ (order 2) and
$-0.03124998978001017$ (order 3), and the outputs above match within
$10^{-6}$ absolute tolerance [klep2022Optimization](@cite).

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

````julia
using NCTSSoS, MosekTools
using NCTSSoS: tr
````

Create unipotent variables (operators that square to identity).

Here we model Alice's observables as `x[1:2]` and Bob's observables as
`y[1:2]`. Variables from different labels ("x" vs "y") commute by construction,
which encodes the bipartite/spatial separation assumption.

````julia
registry, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
````

Identity monomial for converting StatePolynomial to NCStatePolynomial

````julia
const ID_UNIPOTENT = one(NormalMonomial{UnipotentAlgebra,UInt8})
````

`cs_nctssos` minimizes the objective by default. To compute the *maximum*
quantum value of the CHSH expression, we minimize its negative.

````julia
p = -1.0 * tr(x[1] * y[1]) - 1.0 * tr(x[1] * y[2]) - 1.0 * tr(x[2] * y[1]) + 1.0 * tr(x[2] * y[2])

tpop = polyopt(p * ID_UNIPOTENT, registry)
````

`ts_algo` selects a term-sparsity heuristic to reduce SDP size (optional).

````julia
solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=1, ts_algo=MaximalElimination())

result = cs_nctssos(tpop, solver_config)

@show result.objective

@assert isapprox(result.objective, -2.8284271157283083, atol = 1e-5)
````

Our computation matches the theoretical prediction $-2.8284271247461903$ to
within $10^{-6}$ absolute tolerance [klep2022Optimization](@cite).

## Covariance of quantum correlation

This section follows the "covariance Bell inequality" setup: the covariance of
two observables is
```math
\operatorname{Cov}(A, B) = \langle AB \rangle - \langle A \rangle\langle B \rangle.
```
We construct a linear combination of such covariances and compute its optimum
via tracial polynomial optimization.

As in the Bell example, we use unipotent variables to represent $\pm 1$ valued
observables.

Here we introduce two commuting groups of variables, Alice `x[1:3]` and Bob
`y[1:3]`.

````julia
using NCTSSoS, MosekTools
using NCTSSoS: tr

registry, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])

cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
````

As above, we minimize the negative of the covariance Bell expression to obtain
the maximum quantum value.

````julia
p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))
tpop = polyopt(p * ID_UNIPOTENT, registry)

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2)

result = cs_nctssos(tpop, solver_config)

@show result.objective
abs_error = abs(result.objective + 5.0)
@show abs_error
@assert abs_error < 1e-3
````

Again, the literature value is $-5$, and the printed absolute error above
shows the match [klep2022Optimization](@cite).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

