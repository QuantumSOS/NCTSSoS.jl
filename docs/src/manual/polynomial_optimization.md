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
