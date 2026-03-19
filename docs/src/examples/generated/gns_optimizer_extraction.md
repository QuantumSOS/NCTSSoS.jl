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
1 - X^2 \ge 0, \quad 1 - Y^2 \ge 0.
```

The optimal eigenvalue is ``\lambda_{\min}(f) = -1``. The order-2 relaxation
already certifies this bound, but the GNS step needs moments up to degree `6`,
so we solve the order-3 moment relaxation, extract a degree-2 quotient space,
and verify the resulting `5 × 5` model.

**Prerequisites**: familiarity with the
[polynomial optimization API](@ref polynomial-optimization) and the
[SDP relaxation workflow](@ref moment-sohs-hierarchy).

## Setup

We use Mosek in silent mode so the page stays compact.

````julia
using NCTSSoS, MosekTools, LinearAlgebra, Logging

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(
    Mosek.Optimizer,
    MOI.Silent() => true,
);

function evaluate_monomial(matrices::Dict, mono::NormalMonomial)
    sample = first(values(matrices))
    T = eltype(sample)
    dim = size(sample, 1)
    value = Matrix{T}(I, dim, dim)
    for idx in mono.word
        value *= matrices[idx]
    end
    return value
end

moment_from_model(matrices::Dict, xi::AbstractVector, mono::NormalMonomial) = real(dot(xi, evaluate_monomial(matrices, mono) * xi))
moment_from_sdp(monomap::AbstractDict, mono::NormalMonomial) = get(monomap, symmetric_canon(NCTSSoS.expval(mono)), 0.0)
````

`evaluate_monomial` maps a word to its matrix product. The two `moment_*`
helpers compare the extracted model with the solved SDP moments.

## Define the polynomial optimization problem

We keep `X` and `Y` in one noncommutative label group, so there is a single
dense moment matrix. That matches the textbook presentation of the example.

````julia
registry, (vars,) = create_noncommutative_variables([("X", 1:2)]);
X, Y = vars;
````

`registry` is a `VariableRegistry`; `X` and `Y` are the two operator symbols.

````julia
f = 1.0 * X * Y * X;
g1 = 1.0 - X^2;
g2 = 1.0 - Y^2;
````

`f` is the objective polynomial. `g1`, `g2` define the nc polydisc.

````julia
pop = polyopt(f, registry; ineq_constraints=[g1, g2]);
````

`pop` is an ordinary noncommutative `PolyOpt`.

````julia
pop
````

````
Optimization Problem (NonCommutativeAlgebra)
────────────────────────────────────
Objective:
    X₁X₂X₁

Equality constraints (0):
    (none)

Inequality constraints (2):
    1 + -X₁² >= 0
            1 + -X₂² >= 0

Variables (2):
    X₁, X₂

````

## Solve the order-3 moment relaxation directly

The high-level [`cs_nctssos`](@ref) interface returns the bound, but GNS needs
the solved moments themselves. For that reason we use the low-level symbolic
moment workflow:

1. `compute_sparsity`
2. `NCTSSoS.moment_relax`
3. `NCTSSoS.solve_moment_problem`

````julia
solver_config = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=3,
    cs_algo=NoElimination(),
    ts_algo=NoElimination(),
);
````

`solver_config` requests the dense order-3 relaxation.

````julia
sparsity = compute_sparsity(pop, solver_config);
moment_problem = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities);
moment_result = NCTSSoS.solve_moment_problem(moment_problem, SILENT_MOSEK);
````

`moment_result` stores the objective value and the solved numeric moment table
in `moment_result.monomap`.

````julia
moment_summary = (
    objective = moment_result.objective,
    unique_moments = moment_result.n_unique_elements,
    solved_moments = length(moment_result.monomap),
)
moment_summary
````

````
(objective = -0.999999938049795, unique_moments = 78, solved_moments = 78)
````

`moment_summary` shows the certified bound together with the size of the
solved moment table.

````julia
@assert isapprox(moment_result.objective, -1.0; atol=1e-5)
````

The relaxation already certifies the optimal value `-1`.

## Run the GNS reconstruction

We reconstruct from the solved moments with:

- full Hankel degree `3`,
- principal quotient degree `2`.

The returned `GNSResult` contains the multiplication operators, the
distinguished vector `xi`, the basis, and rank information.

````julia
gns = with_logger(Logging.SimpleLogger(devnull, Logging.Error)) do
    gns_reconstruct(moment_result.monomap, registry, 3; hankel_deg=2, atol=1e-6)
end;
````

`gns` is a `GNSResult`; `gns.matrices` stores the extracted operators and
`gns.xi` is the class of the identity word.

````julia
A1 = gns.matrices[registry[:X₁]];
A2 = gns.matrices[registry[:X₂]];
````

`A1` and `A2` are the extracted `5 × 5` symmetric operators.

````julia
extraction_summary = (
    rank = gns.rank,
    full_rank = gns.full_rank,
    xi = round.(gns.xi, digits=6),
    A1 = round.(A1, digits=3),
    A2 = round.(A2, digits=3),
)
extraction_summary
````

````
(rank = 5, full_rank = 9, xi = [-0.995501, -0.0, 0.0, -0.0, -0.094753], A1 = [-0.0 0.996 0.0 -0.0 -0.0; 0.996 0.0 0.0 0.0 0.095; 0.0 0.0 -0.0 0.745 0.0; -0.0 0.0 0.745 0.0 0.0; -0.0 0.095 0.0 0.0 0.0], A2 = [0.0 0.0 -0.769 -0.0 0.0; 0.0 -1.0 -0.0 -0.0 0.0; -0.769 -0.0 -0.0 0.0 0.473; -0.0 -0.0 0.0 -0.0 0.0; 0.0 0.0 0.473 0.0 -0.0])
````

`extraction_summary` exposes the `5 × 5` operators in a compact form.

````julia
@assert gns.rank == 5
@assert size(A1) == (5, 5)
@assert size(A2) == (5, 5)
@assert norm(A1 - A1') ≤ 1e-8
@assert norm(A2 - A2') ≤ 1e-8
@assert isapprox(norm(gns.xi), 1.0; atol=1e-8)
````

The exact matrix entries depend on the orthonormal basis chosen by the SVD, so
they are unique only up to an orthogonal change of basis. The invariant facts
are the dimension, symmetry, and reproduced moments.

## Verify the extracted optimizer

First evaluate the objective polynomial on the extracted operators.

````julia
F = Matrix(A1 * A2 * A1);
F = Symmetric((F + F') / 2);
````

`F` is the symmetric matrix `f(A1, A2) = A1 * A2 * A1`.

````julia
eigenvalues = eigvals(F)
round.(eigenvalues, digits=6)
````

````
5-element Vector{Float64}:
 -1.0
 -0.53687
 -0.0
  0.0
  0.53687
````

`eigenvalues` shows that the extracted model attains the lower bound `-1`.

````julia
@assert isapprox(minimum(eigenvalues), -1.0; atol=1e-5)
@assert isapprox(real(dot(gns.xi, F * gns.xi)), -1.0; atol=1e-5)
````

The smallest eigenvalue is `-1`, so the extracted model attains the SDP bound.

## Verify the moments

The GNS construction should reproduce the solved moments on all words up to the
full Hankel degree.

````julia
basis_deg3 = get_ncbasis(registry, 3);
errors = [
    abs(moment_from_model(gns.matrices, gns.xi, mono) - moment_from_sdp(moment_result.monomap, mono))
    for mono in basis_deg3
];
````

`errors` stores the reconstruction error for every basis word up to degree `3`.

````julia
moment_check = (
    basis_size = length(basis_deg3),
    max_error = maximum(errors),
)
moment_check
````

````
(basis_size = 15, max_error = 3.045198293527562e-8)
````

`moment_check` confirms that the reconstructed model reproduces the SDP
moments up to the requested degree.

````julia
@assert maximum(errors) < 5e-5
````

## Next steps

- `gns_reconstruct`: full GNS API, including `xi`, basis data, and
  rank information.
- [Newton Chip Method](@ref newton-chip-method): basis reduction for ordinary
  and tracial relaxations.
- [Polynomial Optimization](@ref polynomial-optimization): the higher-level
  modeling entry point.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

