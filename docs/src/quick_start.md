# [Quick Start](@id quick-start)

The non-commutative Broyden banded function is a generalization of the classical
Broyden banded function to non-commuting variables. It is often used in
optimization and numerical analysis to test the performance of algorithms. We
will use it as an example.

The function is defined as:
```math
f(x_1, \dots, x_n) = \sum_{i=1}^n f_i(x_1, \dots, x_n)^2
```
where
```math
f_i(x_1, \dots, x_n) = 2x_i + 5x_i^3 + 1 - \sum_{j \in J_i} (x_j +x_j)^2
```
with
```math
J_i = \{j | j \neq i, \max(1, i-5) \leq j \leq \min(n, i+1)\}
```

The variables $x_i$ are non-commuting. You may think of them as matrices or
operators that is assigned with a representation. It's possible puts constraints
on the variables in the form of polynomials of equalities and inequalities. For
example, we may require

```math
1 - x_i^2 \geq 0 \quad \text{and} \quad x_i - \frac{1}{3} \geq 0 \quad \forall i \in [1,n]
```

Firstly, use [`NCTSSoS.PolyOpt`](@ref) object to represent this problem. Since
the constraints are inequalities, we pass it to `ineq_constraints` argument to
constructor [`NCTSSoS.polyopt`](@ref).

```julia quick-start
using NCTSSoS, MosekTools

function broyden_banded(n::Int)
	@ncpolyvar x[1:n]
	f = 0.
	for i = 1:n
	    jset = max(1, i-5) : min(n, i+1)
	    jset = setdiff(jset, i)
	    g = sum(x[j] + x[j]^2 for j in jset)
	    f += (2*x[i] + 5*x[i]^3 + 1 - g)^2
	end
	ineq_cons = [[1 - x[i]^2 for i in 1:n];[x[i] - 1/3 for i in 1:n]]

	return polyopt(f; ineq_constraints=ineq_cons)
end

pop = broyden_banded(6)
```

```julia
obj:

    6.0 * 1 + -6.0 * x₁¹ + -6.0 * x₂¹ + -4.0 * x₃¹ + -2.0 * x₄¹ + 2.0 * x₆¹ + 2.0 * x₁¹x₃¹ + 1.0 * x₁¹x₄¹ + -1.0 * x₁¹x₆¹ + -1.0 * x₁² + -1.0 * x₂¹x₃¹ + 1.0 * x₂¹x₄¹ + -1.0 * x₂¹x₆¹ + -1.0 * x₂² + 2.0 * x₃¹x₁¹ + -1.0 * x₃¹x₂¹ + -2.0 * x₃¹x₄¹ + -1.0 * x₃¹x₆¹ + 1.0 * x₄¹x₁¹ + 1.0 * x₄¹x₂¹ + -2.0 * x₄¹x₃¹ + -3.0 * x₄¹x₅¹ + -1.0 * x₄¹x₆¹ + 1.0 * x₄² + -3.0 * x₅¹x₄¹ + -4.0 * x₅¹x₆¹ + 2.0 * x₅² + -1.0 * x₆¹x₁¹ + -1.0 * x₆¹x₂¹ + -1.0 * x₆¹x₃¹ + -1.0 * x₆¹x₄¹ + -4.0 * x₆¹x₅¹ + 3.0 * x₆² + 2.0 * x₁¹x₂² + 4.0 * x₁¹x₃² + 3.0 * x₁¹x₄² + 2.0 * x₁¹x₅² + 1.0 * x₁¹x₆² + 2.0 * x₁²x₂¹ + 2.0 * x₁²x₃¹ + 1.0 * x₁²x₄¹ + -1.0 * x₁²x₆¹ + 20.0 * x₁³ + 2.0 * x₂¹x₁² + 1.0 * x₂¹x₃² + 3.0 * x₂¹x₄² + 2.0 * x₂¹x₅² + 1.0 * x₂¹x₆² + 2.0 * x₂²x₁¹ + 1.0 * x₂²x₃¹ + 1.0 * x₂²x₄¹ + -1.0 * x₂²x₆¹ + 20.0 * x₂³ + 2.0 * x₃¹x₁² + 1.0 * x₃¹x₂² + 2.0 * x₃¹x₅² + 1.0 * x₃¹x₆² + 4.0 * x₃²x₁¹ + 1.0 * x₃²x₂¹ + -1.0 * x₃²x₆¹ + 18.0 * x₃³ + 1.0 * x₄¹x₁² + 1.0 * x₄¹x₂² + -1.0 * x₄¹x₅² + 1.0 * x₄¹x₆² + 3.0 * x₄²x₁¹ + 3.0 * x₄²x₂¹ + -1.0 * x₄²x₅¹ + -1.0 * x₄²x₆¹ + 16.0 * x₄³ + -1.0 * x₅¹x₄² + -2.0 * x₅¹x₆² + 2.0 * x₅²x₁¹ + 2.0 * x₅²x₂¹ + 2.0 * x₅²x₃¹ + -1.0 * x₅²x₄¹ + -2.0 * x₅²x₆¹ + 14.0 * x₅³ + -1.0 * x₆¹x₁² + -1.0 * x₆¹x₂² + -1.0 * x₆¹x₃² + -1.0 * x₆¹x₄² + -2.0 * x₆¹x₅² + 1.0 * x₆²x₁¹ + 1.0 * x₆²x₂¹ + 1.0 * x₆²x₃¹ + 1.0 * x₆²x₄¹ + -2.0 * x₆²x₅¹ + 12.0 * x₆³ + -5.0 * x₁¹x₂³ + -5.0 * x₁¹x₃³ + -5.0 * x₁¹x₄³ + -5.0 * x₁¹x₅³ + -5.0 * x₁¹x₆³ + 4.0 * x₁²x₂² + 4.0 * x₁²x₃² + 3.0 * x₁²x₄² + 2.0 * x₁²x₅² + 1.0 * x₁²x₆² + -5.0 * x₁³x₂¹ + 25.0 * x₁⁴ + -5.0 * x₂¹x₁³ + -5.0 * x₂¹x₃³ + -5.0 * x₂¹x₄³ + -5.0 * x₂¹x₅³ + -5.0 * x₂¹x₆³ + 4.0 * x₂²x₁² + 3.0 * x₂²x₃² + 3.0 * x₂²x₄² + 2.0 * x₂²x₅² + 1.0 * x₂²x₆² + -5.0 * x₂³x₁¹ + -5.0 * x₂³x₃¹ + 25.0 * x₂⁴ + -5.0 * x₃¹x₂³ + -5.0 * x₃¹x₄³ + -5.0 * x₃¹x₅³ + -5.0 * x₃¹x₆³ + 4.0 * x₃²x₁² + 3.0 * x₃²x₂² + 2.0 * x₃²x₄² + 2.0 * x₃²x₅² + 1.0 * x₃²x₆² + -5.0 * x₃³x₁¹ + -5.0 * x₃³x₂¹ + -5.0 * x₃³x₄¹ + 24.0 * x₃⁴ + -5.0 * x₄¹x₃³ + -5.0 * x₄¹x₅³ + -5.0 * x₄¹x₆³ + 3.0 * x₄²x₁² + 3.0 * x₄²x₂² + 2.0 * x₄²x₃² + 1.0 * x₄²x₅² + 1.0 * x₄²x₆² + -5.0 * x₄³x₁¹ + -5.0 * x₄³x₂¹ + -5.0 * x₄³x₃¹ + -5.0 * x₄³x₅¹ + 23.0 * x₄⁴ + -5.0 * x₅¹x₄³ + -5.0 * x₅¹x₆³ + 2.0 * x₅²x₁² + 2.0 * x₅²x₂² + 2.0 * x₅²x₃² + 1.0 * x₅²x₄² + -5.0 * x₅³x₁¹ + -5.0 * x₅³x₂¹ + -5.0 * x₅³x₃¹ + -5.0 * x₅³x₄¹ + -5.0 * x₅³x₆¹ + 22.0 * x₅⁴ + -5.0 * x₆¹x₅³ + 1.0 * x₆²x₁² + 1.0 * x₆²x₂² + 1.0 * x₆²x₃² + 1.0 * x₆²x₄² + -5.0 * x₆³x₁¹ + -5.0 * x₆³x₂¹ + -5.0 * x₆³x₃¹ + -5.0 * x₆³x₄¹ + -5.0 * x₆³x₅¹ + 21.0 * x₆⁴ + -5.0 * x₁²x₂³ + -5.0 * x₁²x₃³ + -5.0 * x₁²x₄³ + -5.0 * x₁²x₅³ + -5.0 * x₁²x₆³ + -5.0 * x₁³x₂² + -5.0 * x₂²x₁³ + -5.0 * x₂²x₃³ + -5.0 * x₂²x₄³ + -5.0 * x₂²x₅³ + -5.0 * x₂²x₆³ + -5.0 * x₂³x₁² + -5.0 * x₂³x₃² + -5.0 * x₃²x₂³ + -5.0 * x₃²x₄³ + -5.0 * x₃²x₅³ + -5.0 * x₃²x₆³ + -5.0 * x₃³x₁² + -5.0 * x₃³x₂² + -5.0 * x₃³x₄² + -5.0 * x₄²x₃³ + -5.0 * x₄²x₅³ + -5.0 * x₄²x₆³ + -5.0 * x₄³x₁² + -5.0 * x₄³x₂² + -5.0 * x₄³x₃² + -5.0 * x₄³x₅² + -5.0 * x₅²x₄³ + -5.0 * x₅²x₆³ + -5.0 * x₅³x₁² + -5.0 * x₅³x₂² + -5.0 * x₅³x₃² + -5.0 * x₅³x₄² + -5.0 * x₅³x₆² + -5.0 * x₆²x₅³ + -5.0 * x₆³x₁² + -5.0 * x₆³x₂² + -5.0 * x₆³x₃² + -5.0 * x₆³x₄² + -5.0 * x₆³x₅² + 25.0 * x₁⁶ + 25.0 * x₂⁶ + 25.0 * x₃⁶ + 25.0 * x₄⁶ + 25.0 * x₅⁶ + 25.0 * x₆⁶

constraints:


    1.0 * 1 + -1.0 * x₁² >= 0
    1.0 * 1 + -1.0 * x₂² >= 0
    1.0 * 1 + -1.0 * x₃² >= 0
    1.0 * 1 + -1.0 * x₄² >= 0
    1.0 * 1 + -1.0 * x₅² >= 0
    1.0 * 1 + -1.0 * x₆² >= 0
    -0.3333333333333333 * 1 + 1.0 * x₁¹ >= 0
    -0.3333333333333333 * 1 + 1.0 * x₂¹ >= 0
    -0.3333333333333333 * 1 + 1.0 * x₃¹ >= 0
    -0.3333333333333333 * 1 + 1.0 * x₄¹ >= 0
    -0.3333333333333333 * 1 + 1.0 * x₅¹ >= 0
    -0.3333333333333333 * 1 + 1.0 * x₆¹ >= 0

variables:
    x₁ x₂ x₃ x₄ x₅ x₆

is_unipotent:
    false

is_projective:
    false
```

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

Information about the relaxation is encoded in [`NCTSSoS.SolverConfig`](@ref).
Besides moment order, it also needs to be provided with a [SDP Solver](@ref
overview-of-optimizers). Different optimizers may have different performance
characteristics and may be more or less suitable for a given problem.

```julia quick-start
solver_config = SolverConfig(optimizer=Mosek.Optimizer;order=3)
```

Finally, we are ready to solve this problem. This is accomplishded with [`NCTSSoS.cs_nctssos`](@ref).

```julia quick-start
@time result = cs_nctssos(pop, solver_config)

@assert isapprox(result.objective, 2.979657980133734; atol=1e-5) # find correct value
```

```julia
131.787849 seconds (80.18 M allocations: 4.424 GiB, 0.67% gc time, 4.86% compilation time)
Objective: 2.979657977586228
```

Although we have reached a tight bound, time to solve the problem can still be
significant. As a remedy, [Sparsities](@ref sparsities) can be utilized to
reduce the problem size. This is achieved by supplying an
[`EliminationAlgorithm`](https://algebraicjulia.github.io/CliqueTrees.jl/stable/api/#Elimination-Algorithms)
to [`NCTSSoS.SolverConfig`](@ref).

```julia
solver_config = SolverConfig(optimizer=Mosek.Optimizer; order=3, cs_algo=MF(), ts_algo=MMD())

@time result_ts_cs = cs_nctssos(pop, solver_config)

@assert result_ts_cs.objective <= result.objective
```

```julia
1.761473 seconds (46.53 M allocations: 1.740 GiB, 13.13% gc time, 9.54% compilation time)
Objective: 2.979657981888441
```

As expected, the time taken to solve the problem has been significantly reduced.
However, the result may no longer be a tight lower bound. Sparsity is itself a kind
of relaxation. Luckily, we may tighten the relaxation in the [Term
Sparsity](@ref term-sparsity) sense, this is done with
[`NCTSSoS.cs_nctssos_higher`](@ref).

```julia quick-start
@time result_higher = cs_nctssos_higher(pop, result_ts_cs, solver_config)

@assert result_higher.objective ≈ result.objective atol=1e-5
```

```julia
14.902031 seconds (444.36 M allocations: 15.470 GiB, 8.69% gc time, 0.46% compilation time)
Objective: 2.9796579998271047
```

## Workflow

To summarize, the workflow for solving polynomial optimization can be summarized as

![`Workflow for solving Polynomial Optimization problem`](assets/overall_workflow.typ.svg)

If you would like to understand more please refer to [examples section](@ref "bell-inequalities") in the document.
