# [Quick Start](@id quick-start)

The non-commutative Broyden banded function is a generalization of the classical
Broyden banded function to non-commuting variables. It is often used in
optimization and numerical analysis to test the performance of algorithms. We
will use it as an example.

The function is defined as:
```math
f(x_1, \dots, x_n) = \sum_{i=1}^n f_i(x_1, \dots, x_n)^\dagger f_i(x_1, \dots, x_n)
```
where
```math
f_i(x_1, \dots, x_n) = 2x_i + 5x_i^3 + 1 - \sum_{j \in J_i} (x_j + x_j^2)
```
with
```math
J_i = \{j | j \neq i, \max(1, i-5) \leq j \leq \min(n, i+1)\}
```

The variables $x_i$ are non-commuting. You may think of them as matrices or
operators that is assigned with a representation. It's possible to put constraints
on the variables in the form of polynomial inequalities. For example, we may require:

```math
1 - x_i^2 \geq 0 \quad \text{and} \quad x_i - \frac{1}{3} \geq 0 \quad \forall i \in [1,n]
```

## Setting Up the Problem

First, create non-commutative variables using `site_variables(Free(), 1, n)`.
All variables belong to one site, so no commutation between them is imposed.
Then use `Problem` with `Positive` to define the operator inequalities.

```julia
using NCTSSoS.Next, Mosek
import NCTSSoS.Next as Next

function broyden_banded(n::Int)
    # Create non-commutative variables
    registry, (x,) = site_variables(Free(), 1, n)

    # Build the objective function using sum
    f = sum(1:n) do i
        jset = setdiff(max(1, i-5):min(n, i+1), i)
        g = isempty(jset) ? 0.0 * x[1] : sum(x[j] + x[j]^2 for j in jset)
        (2.0*x[i] + 5.0*x[i]^3 + 1 - g)^2
    end

    # Define inequality constraints
    ineq_cons = [[1.0 - x[i]^2 for i in 1:n]; [x[i] - 1//3 for i in 1:n]]

    return registry, Problem(f; constraints=Positive.(ineq_cons))
end

registry, problem = broyden_banded(6)
```

The `site_variables` function returns a registry and
variable arrays. All algebraic constraints are automatically encoded through
the type system.

[Polynomial Optimization](@ref polynomial-optimization) is NP-hard, therefore it
is considered impossible to solve them efficiently in general. However, it is
possible to relax the problem into [Semidefinite Programming](@ref
semidefinite-programming). The solution of the Semidefinite program will be the
lower/upper bound of the original minimization/maximization polynomial
optimziation problem.

The relaxation can be done in two different forms, [moment relaxation](@ref
moment-problem) and [sum of hermitian square relaxation](@ref sohs-problem). The
relaxation is tight in the limit of a parameter, moment order, reaching
infinity. However, you may get lucky and be able to sovle the problem at finite
moment order.

Information about the relaxation is encoded in `Relaxation`.
Choose the moment order here and pass a [SDP Solver](@ref
overview-of-optimizers) to `solve`. Different optimizers may have different performance
characteristics and may be more or less suitable for a given problem.

```julia
relaxation = Relaxation(order=3)
solver = MosekBackend()
```

If you want full control over the relaxation basis, you can provide
`rows` directly instead of `order`:

```julia
rows = collect(Next.basis(registry, 3))
custom_relaxation = Relaxation(rows=rows)
```

## Solving the Problem

Finally, we are ready to solve this problem using `solve`.

```julia
result = solve(problem, relaxation; solver)
# Objective: ~2.9796579831938192
@assert isapprox(result.bound, 2.97965798; atol=1e-6)
```

Although we have reached a tight bound, time to solve the problem can still be
significant. As a remedy, [Sparsities](@ref sparsities) can be utilized to
reduce the problem size. Choose a correlative grouping and a term split,
each with an
[`EliminationAlgorithm`](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#Elimination-Algorithms)
on `Relaxation`. Name `MMD()` explicitly: `TermSplit` defaults to
`MaximalElimination()`.

```julia
relaxation_sparse = Relaxation(order=3,
    grouping=Correlative(algorithm=MF()),
    split=TermSplit(algorithm=MMD(), level=1))

result_sparse = solve(problem, relaxation_sparse; solver)
# Objective: ~2.9796579733401436
@assert isapprox(result_sparse.bound, result.bound; atol=1e-6)
```

Sparse blocks can reduce the cost of solving the problem.
However, the result may no longer be a tight lower bound since sparsity is itself
a kind of relaxation. Increasing `TermSplit(level=...)` can sometimes tighten
the relaxation in the [Term Sparsity](@ref term-sparsity) sense by extending
the completed moment and localizing supports at the same order.

This is not guaranteed to improve the bound: if the raw term-sparsity support is
already at a fixed point, the next sparse step will return the same relaxation.
In that case, use a higher `order` or a less aggressive sparsity configuration.

```julia
relaxation_higher = Relaxation(order=3,
    grouping=Correlative(algorithm=MF()),
    split=TermSplit(algorithm=MMD(), level=2))
result_higher = solve(problem, relaxation_higher; solver)
# Objective: ~2.979657987721403
@assert isapprox(result_higher.bound, result.bound; atol=1e-6)
```

## Workflow

To summarize, the workflow for solving polynomial optimization can be summarized as

![`Workflow for solving Polynomial Optimization problem`](assets/overall_workflow.typ.svg)

If you would like to understand more please refer to [examples section](@ref "bell-inequalities") in the document.
