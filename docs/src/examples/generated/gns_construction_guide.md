```@meta
EditURL = "../literate/gns_construction_guide.jl"
```

# [GNS Construction Interface Guide](@id gns-construction-guide)

## What is the GNS construction?

A moment relaxation gives numerical expectation values ``L(w)`` for words
``w``. The GNS (Gelfand-Naimark-Segal) construction seeks matrices ``A_i``
and a unit vector ``\xi`` such that
``\langle\xi,w(A)\xi\rangle = L(w)``.

Solved moments need not admit this finite-dimensional representation.
`gns(result)` checks the reconstruction before it returns operators.
It can also return a feasible optimizer that attains the bound but does
**not** reproduce all the moments. Read its verdict before using the model.
These are numerical checks at stated tolerances, not rigorous certificates.

## Key concepts

**Moment matrix (Hankel matrix).** For an ordered basis of words ``b_i``,
its entries are

```math
H_{ij} = L(b_i^\dagger b_j).
```

**Relaxation order and Hankel degree.** A dense order-``d`` relaxation
provides moments up to degree ``2d``. A Hankel matrix with row words up to
degree ``k`` needs moments up to degree ``2k``. Here the relaxation has
order 3, and the selected full Hankel has row degree 3.

**Principal block and full Hankel.** `gns` selects an adjacent degree pair.
`basis` contains the principal row words; `full_basis` contains the full
row words. The principal numerical rank sets the operator dimension.

**Flatness.** Equal principal and full ranks, with small extension
residuals, mean that the higher-degree rows add no new directions. The
returned diagnostics use the actual solved moments; the API does not
replace their extension block to force equal ranks.

## API summary

| Call or field | Purpose |
|:---|:---|
| `gns(result)` | Reconstruct and check a solved result |
| `g.verdict`, `g.status`, `g.advice` | Read what was accepted or withheld |
| `g[x]` | Read an accepted variable matrix |
| `g.state` | Read the cyclic vector |
| `evaluate(g, p)` | Evaluate a polynomial as a matrix |
| `expectation(g, p)` | Evaluate a polynomial in the reconstructed state |
| `g.flatness`, `g.checks` | Inspect numerical checks |

**Prerequisites**: the [polynomial optimization API](@ref polynomial-optimization)
and the [SDP relaxation workflow](@ref moment-sohs-hierarchy).

````julia
using NCTSSoS.Next, Mosek, LinearAlgebra
````

## Step 1 — Solve the noncommutative unit-ball problem

We minimize

```math
f = 2-X^2+XY^2X-Y^2
```

subject to ``1-X^2-Y^2\succeq0``. The optimal value is 1.

````julia
registry, (vars,) = site_variables(Free(), 1, 2);
X, Y = vars;
f = 2.0 - X^2 + X * Y^2 * X - Y^2;
ball = 1.0 - X^2 - Y^2;
problem = Problem(f; constraints=(Positive(ball),));
result = solve(problem, Relaxation(order=3); solver=MosekBackend());
result.bound
````

````
0.9999999943399142
````

## Step 2 — Read the checked reconstruction

The rank cutoff is explicit here. `check_atol` sets the absolute tolerance
for matrix, moment and objective checks.

````julia
g = gns(result; atol=1e-6, check_atol=5e-5);
g.verdict
````

````
:optimizer_certified
````

This example returns `:optimizer_certified`, not `:moments_reproduced`.
The operator model attains the bound and satisfies the ball constraint,
but it does not reproduce every solved moment.

The three possible verdicts have different meanings:

| Verdict | Meaning |
|:---|:---|
| `:moments_reproduced` | A positive, flat pair passes the moment and model checks |
| `:optimizer_certified` | The model satisfies the problem and attains the bound; moment reproduction is not certified |
| `:withheld` | Operators and state are `nothing`; inspect `status`, `checks` and `advice` |

Do not treat an optimizer-only verdict as a moment reconstruction.
`expectation(result, p)` reads solved moments, while `expectation(g, p)`
reads the returned model. They can differ.

## Step 3 — Inspect degrees, bases and ranks

````julia
(
    solved_order = g.relaxation_order,
    degrees = (g.principal_degree, g.full_degree),
    basis_sizes = (length(g.basis), length(g.full_basis)),
    ranks = (g.rank, g.full_rank),
)
````

````
(solved_order = 3, degrees = (2, 3), basis_sizes = (7, 15), ranks = (5, 9))
````

The selected degrees are 2 and 3, with 7 and 15 row words. The actual
ranks are **5 → 9**, not 5 → 5. The full matrix is not flat.
Earlier versions of this guide replaced the extension block before
extraction. The checked API keeps the original moments and reports this
difference instead.

Singular values show which directions the cutoff retains:

````julia
(principal=g.singular_values, full=g.full_singular_values)
````

````
(principal = [1.5012343702087998, 0.523614468909496, 0.4763855260778388, 0.1447750348094701, 0.11914702764208825, 1.4518008115105776e-9, 4.397520930571975e-10], full = [1.501234370208799, 0.8422166125249535, 0.7310212130853533, 0.14477503480947027, 0.1191470276420883, 0.04009663824207139, 0.03308724840421216, 0.03117658729766389, 0.024772450148146257, 2.5033927668499987e-9, 1.5543787634358816e-9, 1.4518008089914379e-9, 4.80869466570765e-10, 4.397520881391201e-10, 3.622578578784968e-10])
````

The flatness record includes positivity, cutoffs, singular-value gaps,
the minimum eigenvalue, and range and extension residuals:

````julia
collect(pairs(g.flatness))
````

````
11-element Vector{Pair{Symbol, Any}}:
            :is_flat => false
        :is_positive => true
             :cutoff => 1.0e-6
        :full_cutoff => 1.0e-6
                :gap => (retained = 0.11914702764208825, discarded = 1.4518008115105776e-9, resolved = true, automatic = false)
           :full_gap => (retained = 0.024772450148146257, discarded = 2.5033927668499987e-9, resolved = true, automatic = false)
        :range_error => 1.0779543221227156e-16
    :extension_error => 0.03673409612098843
 :minimum_eigenvalue => -1.5543787634358816e-9
      :psd_tolerance => 1.501234370208799e-8
         :check_atol => 5.0e-5
````

## Step 4 — Read operators and evaluate the objective

The principal rank still gives a `5 × 5` model. Matrix entries depend on
the SVD basis, but the objective and constraint checks do not.

````julia
A1 = g[X];
A2 = g[Y];
(A1=round.(A1, digits=3), A2=round.(A2, digits=3), state=round.(g.state, digits=6))


F = evaluate(g, f);
(
    objective = expectation(g, f),
    objective_error = abs(expectation(g, f) - result.bound),
    constraint_min_eigenvalue = eigmin(Hermitian(evaluate(g, ball))),
)
````

````
(objective = 1.0000000041331605, objective_error = 9.793246324996119e-9, constraint_min_eigenvalue = 4.23202513267292e-9)
````

## Step 5 — Separate objective checks from moment checks

````julia
collect(pairs(g.checks))
````

````
14-element Vector{Pair{Symbol, Any}}:
        :hankel_moment_error => 0.03308724840421217
       :program_moment_error => 0.03308724840421217
         :program_components => 78
            :objective_error => 9.793246324996119e-9
 :constraint_min_eigenvalues => [4.232025094404645e-9]
            :equality_errors => [0.0]
              :sector_errors => Float64[]
              :algebra_error => 2.2887833992611187e-16
        :normalization_error => 1.4658474434270374e-10
           :distance_to_flat => 0.04127091054790795
                  :sigma_min => 0.11914702764208825
                  :sigma_max => 1.5012343702087998
           :condition_number => 12.599847431514895
    :distance_over_sigma_min => 0.34638640480301125
````

`program_moment_error` checks all stored program moment components;
`hankel_moment_error` checks the selected full Hankel. The program error
is about 0.033 here, far above `5e-5`. The small objective error and the
nonnegative constraint eigenvalue still support the optimizer verdict.

Conditioning is also reported: `sigma_min`, `sigma_max`, and
`condition_number` describe the retained principal subspace.
`distance_to_flat` and `distance_over_sigma_min` describe the departure
from flatness. A small objective error alone says nothing about these.

## Tolerances and withheld models

- `atol=nothing` selects a resolved singular-value gap above the numerical
  tail. An explicit `atol` sets the absolute rank cutoff.
- `rtol * σmax` is a lower bound on that cutoff (default `rtol=0`).
- `check_atol` controls the absolute matrix, moment and objective errors.
- `psd_rtol * max(1, σmax)` controls allowed full-Hankel negativity
  (default `psd_rtol=1e-8`).

If a result is withheld, inspect its singular values, residuals and advice.
Improve solver accuracy or increase the relaxation order when needed.
Do not raise tolerances solely to change the verdict.

**See also:**
- [GNS Optimizer Extraction](@ref gns-optimizer-extraction) — a textbook
  example with the noncommutative polydisc.
- [Recovering the Textbook CHSH Optimizer](@ref chsh-gns-reconstruction) —
  a flat model with explicit gauge fixing.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

