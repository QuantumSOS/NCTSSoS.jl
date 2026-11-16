# [Polynomial Optimization](@id polynomial-optimization)

Polynomial optimization is a [mathematical optimization
problem](https://en.wikipedia.org/wiki/Mathematical_optimization) where the
objective function is a polynomial of variable unknowns and the constraints are
polynomial inequalities and equalities. In general, finding the exact solution
to a polynomial optimization problem is
NP-hard[wangExploitingTermSparsity2021](@cite). However, [moment](@ref
moment-problem) and [sum-of-squares](@ref sohs-problem) relaxations can be
applied to obtain lower/upper bounds to polynomial optimization problems
efficiently.

## Solver results

`Problem` defines the objective and operator constraints. `Relaxation`
selects the order, sparsity and strengthening. `solve` takes both and a solver:

```julia
using NCTSSoS
using Mosek

_, (x,) = site_variables(Free(), 1, 1)
problem = Problem(x[1]; constraints=(Positive(1 - x[1]^2),))
relaxation = Relaxation(order=1)
solver = MosekBackend()
result = solve(problem, relaxation; solver)
@assert isapprox(result.bound, -1; atol=1e-6)
@assert isapprox(expectation(result, x[1]), result.bound; atol=1e-6)
```

The constraint is the operator inequality ``1-x_1^2\succeq0``. Its minimum
is −1, attained by the scalar operator ``x_1=-1``.

`result.bound` is a numerical relaxation optimum, not an interval-certified
bound. `result.status` gives the solver termination status. If the solver
does not return an accepted solution, `result.bound` is `nothing`, not `NaN`.
Check it before you use the bound or call `expectation(result, p)`.
`result.sizes` describes the SDP, and `result.valid_for` states which states
the relaxation covers. These fields apply to dense, sparse and symmetric
relaxations. Solver failures can also raise errors; a returned result is not
a guarantee that a bound is available.

## [Noncommutative Polynomial Optimization](@id noncommutative-polynomial-optimization)

Noncommutative polynomial optimization concerns when variables are noncommuting.
These variables can be thought of as matrices or operators acting on an
(infinite dimension) separable Hilbert space $\mathcal{H}$. The numerical value
of the objective function may take two meanings.

- The eigenvalue of the polynomial of operators, giving a general form of minimization problem as
```math
\mathrm{inf}_{\mathbf{x}\in\mathcal{B}(\mathcal{H})^n}\ \lambda_{\min}(f(\mathbf{x}))\ \ \text{ s.t. }\ g_1(\mathbf{x})\ge0,\ldots,g_m(\mathbf{x})\ge0,h_1(\mathbf{x})=0,\ldots,h_{\ell}(\mathbf{x})=0
```
- The trace of the polynomial of operators, giving a general form of minimization problem as

```math
\mathrm{inf}_{\mathbf{x}\in\mathcal{B}(\mathcal{H})^n}\ \mathrm{tr}(f(\mathbf{x})) \ \text{ s.t. }\ g_1(\mathbf{x})\ge0,\ldots,g_m(\mathbf{x})\ge0,h_1(\mathbf{x})=0,\ldots,h_{\ell}(\mathbf{x})=0
```

where
```math
f,g_1,\ldots,g_m,h_1,\ldots,h_{\ell}\in\mathbb{R}\langle\mathbf{x}\rangle
```
are noncommutative polynomials in noncommuting variables $\mathbf{x}$.
