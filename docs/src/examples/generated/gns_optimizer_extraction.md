```@meta
EditURL = "../literate/gns_optimizer_extraction.jl"
```

# [GNS Optimizer Extraction](@id gns-optimizer-extraction)

In this page, we reproduce Example 4.23 from
[burgdorfOptimizationPolynomialsNonCommuting2016](@cite) — the GNS extraction
workflow applied to

```math
f = X Y X
```

on the noncommutative polydisc

```math
1-X^2\succeq0,\qquad 1-Y^2\succeq0.
```

The optimal eigenvalue is ``\lambda_{\min}(f)=-1``. The order-2 relaxation
already gives this bound. For GNS we solve at order 3, with moments up to
degree 6, and check the resulting `5 × 5` model.

**Prerequisites**: familiarity with the
[polynomial optimization API](@ref polynomial-optimization) and the
[SDP relaxation workflow](@ref moment-sohs-hierarchy).

## Setup and problem

The two free variables share one site. `Relaxation(order=3)` requests
the dense moment relaxation. The direct Mosek backend is silent by default.

````julia
using NCTSSoS.Next, Mosek, LinearAlgebra

registry, (vars,) = site_variables(Free(), 1, 2);
X, Y = vars;
f = 1.0 * X * Y * X;
g1 = 1.0 - X^2;
g2 = 1.0 - Y^2;
problem = Problem(f; constraints=(Positive(g1), Positive(g2)));
````

## Solve the order-3 moment relaxation

`solve` returns both the bound and the moments needed by `gns`.

````julia
result = solve(problem, Relaxation(order=3); solver=MosekBackend());
[
    :objective => result.bound,
    :objective_moment => expectation(result, f),
    :sizes => result.sizes,
]
````

````
3-element Vector{Pair{Symbol, Any}}:
        :objective => -1.0000000037277186
 :objective_moment => -0.9999999934607839
            :sizes => 15 rows, 78 moments, 77 real variables, 3 blocks (largest 15); localizing rows: [7, 7]
````

The bound agrees with -1 at tolerance `1e-5`.

````julia
@assert isapprox(result.bound, -1.0; atol=1e-5)
````

## Run the checked GNS reconstruction

`gns` selects the principal and full degree pair and checks the model.
Inspect the verdict before reading operators. These are numerical checks,
not a rigorous certificate.

````julia
g = gns(result; atol=1e-6, check_atol=5e-5);
g.verdict
````

````
:optimizer_certified
````

This solve returns an accepted optimizer. Moment checks follow below.

````julia
@assert g.verdict != :withheld

A1 = g[X];
A2 = g[Y];
extraction_summary = (
    degrees = (g.principal_degree, g.full_degree),
    rank = g.rank,
    full_rank = g.full_rank,
    state = round.(g.state, digits=6),
    A1 = round.(A1, digits=3),
    A2 = round.(A2, digits=3),
)
extraction_summary
````

````
(degrees = (2, 3), rank = 5, full_rank = 9, state = ComplexF64[-0.995358 - 0.0im, 0.0 - 0.0im, -0.0 - 0.0im, 0.0 - 0.0im, 0.096241 - 0.0im], A1 = ComplexF64[-0.0 + 0.0im -0.995 + 0.0im 0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im; -0.995 + 0.0im 0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im 0.096 + 0.0im; 0.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im -0.739 + 0.0im 0.0 + 0.0im; -0.0 + 0.0im 0.0 + 0.0im -0.739 + 0.0im 0.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.096 + 0.0im 0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im], A2 = ComplexF64[0.0 + 0.0im -0.0 + 0.0im -0.785 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; -0.0 + 0.0im -1.0 + 0.0im 0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im; -0.785 + 0.0im 0.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im -0.448 + 0.0im; 0.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im -0.448 + 0.0im -0.0 + 0.0im 0.0 + 0.0im])
````

The principal degree is 2 and the full degree is 3. The principal
rank gives `5 × 5` operators. Matrix entries depend on the SVD basis;
the operator dimension, objective and constraint values do not.

````julia
@assert (g.principal_degree, g.full_degree) == (2, 3)
@assert g.rank == 5
@assert size(A1) == size(A2) == (5, 5)
@assert norm(A1 - A1') <= 1e-8
@assert norm(A2 - A2') <= 1e-8
@assert isapprox(norm(g.state), 1.0; atol=1e-8)
````

## Verify the extracted optimizer

`evaluate` substitutes the extracted operators into a polynomial;
`expectation` also uses the reconstructed state.

````julia
F = Hermitian(evaluate(g, f));
eigenvalues = eigvals(F);
round.(eigenvalues, digits=6)
````

````
5-element Vector{Float64}:
 -1.0
 -0.545614
 -0.0
  0.0
  0.545614
````

The nonminimal levels are not fixed by the problem: the optimal moments
are not unique. This run with the default SOS form gives levels near
±0.545614. A separate run with `MosekBackend(form=:moment)` gives levels
near ±0.536868. These are two examples of optimal moment choices, not
properties of the forms. The models are not merely different bases of
the same matrices. Both attain -1 within `1e-5`.
Check that value and both constraints:

````julia
@assert isapprox(minimum(eigenvalues), -1.0; atol=1e-5)
@assert isapprox(expectation(g, f), -1.0; atol=1e-5)
@assert eigmin(Hermitian(evaluate(g, g1))) >= -5e-5
@assert eigmin(Hermitian(evaluate(g, g2))) >= -5e-5
````

The model is feasible and attains the bound -1 within the stated tolerance.
This checks the optimizer, not just an eigenvalue of an unconstrained matrix.

## Check moments separately

The checked result reports errors against the solved moments. The full
program check includes moments beyond the row words of the selected Hankel.

````julia
(
    verdict = g.verdict,
    hankel_moment_error = g.checks.hankel_moment_error,
    program_moment_error = g.checks.program_moment_error,
    objective_error = g.checks.objective_error,
)
````

````
(verdict = :optimizer_certified, hankel_moment_error = 0.14884636138021173, program_moment_error = 0.14884636138021173, objective_error = 7.135449120454496e-9)
````

As a direct comparison, check every word through degree 3 with the public
polynomial interface. Start with the identity and append `X` or `Y`.

````julia
words = [one(X)];
layer = [one(X)];
for degree in 1:3
    global layer = [w * x for w in layer for x in (X, Y)]
    append!(words, layer)
end
errors = [abs(expectation(g, w) - expectation(result, w)) for w in words];
(words_checked=length(words), max_error=maximum(errors))
````

````
(words_checked = 15, max_error = 3.400235204509272e-9)
````

The 15 words pass the original degree-3 tolerance.

````julia
@assert maximum(errors) < 5e-5
````

Here the verdict is `:optimizer_certified`: the ranks are 5 → 9 and the
full program moment error is about 0.148846. The smaller degree-3 set agrees,
but this does not establish reproduction of all program moments. See the
[GNS Construction Guide](@ref gns-construction-guide) for the distinction.

**See also:**
- [Newton Chip Method](@ref newton-chip-method): basis reduction.
- [Polynomial Optimization](@ref polynomial-optimization): modeling entry point.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

