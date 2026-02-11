# Tracial Polynomial Optimization

Tracial polynomial optimization minimizes polynomial expressions involving
traces of noncommutative operators. This arises naturally in quantum
information when optimizing over quantum states via the moment-SOS
hierarchy [klep2022Optimization](@cite).

This example covers three problems of increasing complexity:

1. A **toy problem** with projector variables — minimal setup to learn the API.
2. The **CHSH Bell inequality** — finding the Tsirelson bound $2\sqrt{2}$.
3. A **covariance Bell inequality** — nonlinear tracial objective involving
   products of expectation values.

**Prerequisites**: familiarity with
tracial polynomial concepts and
the basic polynomial optimization API.

## Setup

We use MosekTools as the SDP solver backend. To silence solver output,
we wrap it with `MOI.OptimizerWithAttributes`.

````julia
using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true);
````

---
## Toy Example: Projector Trace Polynomial

We minimize a tracial polynomial over three projector variables $P_1, P_2, P_3$
satisfying $P_i^2 = P_i$.

The objective is:
```math
f = \operatorname{tr}(P_1 P_2 P_3) + \operatorname{tr}(P_1 P_2)\,\operatorname{tr}(P_3)
```

#### Step 1: Create projector variables

````julia
registry, (x,) = create_projector_variables([("x", 1:3)]);
````

registry: stores variable-to-index mapping and the ProjectorAlgebra constraint
x: length-3 vector of NormalMonomial{ProjectorAlgebra, UInt8}


#### Step 2: Build the tracial objective

`tr` converts a `Polynomial` into a `StatePolynomial` — a symbolic
expression of trace moments.  To pass it to `polyopt`, we embed it as
an `NCStatePolynomial` by multiplying with the identity monomial.

````julia
ID_PROJ = one(NormalMonomial{ProjectorAlgebra, UInt8})
````

````
𝟙
````

ID_PROJ: identity element 𝟙 — the multiplication converts StatePolynomial → NCStatePolynomial

````julia
p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * ID_PROJ;
````

#### Step 3: Create the optimization problem

````julia
spop = polyopt(p, registry);
````

spop: PolyOpt wrapping objective + registry (algebra constraints encoded automatically)


#### Step 4: Solve at relaxation order 2

`SolverConfig` specifies the SDP backend and the moment/SOS hierarchy
`order`. Higher order yields tighter bounds at the cost of larger SDPs.

````julia
solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2);

result = cs_nctssos(spop, solver_config);
````

result.objective: certified SDP lower bound on f

````julia
@show result.objective
@assert isapprox(result.objective, -0.046717378455438933, atol=1e-6)
````

````
result.objective = -0.046717378455481205

````

#### Step 5: Tighten the bound at order 3

````julia
solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=3);

result = cs_nctssos(spop, solver_config);

@show result.objective
@assert isapprox(result.objective, -0.03124998978001017, atol=1e-6)
````

````
result.objective = -0.03124998978003755

````

The literature values are $-0.0467$ (order 2) and $-0.0312$ (order 3); our
results match within $10^{-6}$ absolute tolerance [klep2022Optimization](@cite).

---
## CHSH Bell Inequality (Tracial Form)

The CHSH inequality bounds correlations between two parties with $\pm 1$-valued
observables $A_1, A_2$ (Alice) and $B_1, B_2$ (Bob):

```math
\mathcal{B}_{\text{CHSH}}
  = \operatorname{tr}(A_1 B_1) + \operatorname{tr}(A_1 B_2)
  + \operatorname{tr}(A_2 B_1) - \operatorname{tr}(A_2 B_2).
```

Classical bound: $\mathcal{B} \leq 2$. Quantum bound (Tsirelson):
$\mathcal{B} \leq 2\sqrt{2} \approx 2.828$.

Here we use `UnipotentAlgebra` ($U^2 = I$) to model the observables. Since
`cs_nctssos` *minimizes*, we minimize $-\mathcal{B}$ and expect
$\approx -2\sqrt{2}$.

#### Step 1: Create unipotent variables

Variables from different labels (`"x"` vs `"y"`) commute, encoding the
bipartite/spatial-separation assumption.

````julia
registry, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)]);
````

x: Alice's observables [A₁, A₂]
y: Bob's observables [B₁, B₂]

````julia

ID_UNIP = one(NormalMonomial{UnipotentAlgebra, UInt8});
````

ID_UNIP: identity element for StatePolynomial → NCStatePolynomial conversion


#### Step 2: Define the (negated) CHSH tracial expression

````julia
p = -1.0 * tr(x[1] * y[1]) +  ## −tr(A₁B₁)
    -1.0 * tr(x[1] * y[2]) +  ## −tr(A₁B₂)
    -1.0 * tr(x[2] * y[1]) +  ## −tr(A₂B₁)
     1.0 * tr(x[2] * y[2]);   ## +tr(A₂B₂)
````

#### Step 3: Solve with term-sparsity exploitation

`ts_algo = MaximalElimination()` decomposes the SDP into independent blocks
via connected components of the variable interaction graph.

````julia
tpop = polyopt(p * ID_UNIP, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=1, ts_algo=MaximalElimination());

result = cs_nctssos(tpop, solver_config);

@show result.objective
@assert isapprox(result.objective, -2 * sqrt(2), atol=1e-5)
````

````
result.objective = -2.8284271321623193

````

The result matches the Tsirelson bound $-2\sqrt{2} \approx -2.828$
[klep2022Optimization](@cite).

---
## Covariance Bell Inequality

Nonlinear Bell inequalities involve products of expectation values, not just
single traces. The covariance of observables $A_i, B_j$ is:

```math
\operatorname{Cov}(A_i, B_j)
  = \langle A_i B_j \rangle
  - \langle A_i \rangle\,\langle B_j \rangle.
```

We compute the maximum quantum value of the covariance Bell expression
[pozsgay2017Covariance](@cite):

```math
f = \sum_{(i,j)\in S^+} \operatorname{Cov}(A_i, B_j)
  - \sum_{(i,j)\in S^-} \operatorname{Cov}(A_i, B_j)
```

where $S^+ = \{(1,1),(1,2),(1,3),(2,1),(2,2),(3,1)\}$ and
$S^- = \{(2,3),(3,2)\}$.
Classical bound: $f \leq 4.5$.  Quantum bound: $f = 5$.

> **tr vs ς**
>
> `tr` builds a `StatePolynomial{MaxEntangled}` — a moment in the
> maximally entangled (tracial) state.  `ς` (type `\varsigma` + Tab)
> builds a `StatePolynomial{Arbitrary}` — a moment optimized over all states.
> The covariance Bell inequality requires optimizing over arbitrary states,
> so we use `ς` here.

#### Step 1: Create unipotent variables (3 per party)

````julia
registry, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)]);
````

x: Alice's observables [A₁, A₂, A₃]
y: Bob's observables [B₁, B₂, B₃]

````julia

ID_UNIP = one(NormalMonomial{UnipotentAlgebra, UInt8});
````

#### Step 2: Define covariance using arbitrary-state expectation values

Each factor `ς(·) * ID_UNIP` embeds one expectation value as an
`NCStatePolynomial`.  Their product captures the nonlinear structure
$\langle\cdot\rangle - \langle\cdot\rangle\langle\cdot\rangle$ symbolically.

````julia
cov(a, b) = 1.0 * ς(x[a] * y[b]) * ID_UNIP -
            1.0 * ς(x[a]) * ς(y[b]) * ID_UNIP;
````

#### Step 3: Build and solve the (negated) objective

````julia
p = -1.0 * (
    cov(1, 1) + cov(1, 2) + cov(1, 3) +   ## Cov(A₁,B₁) + Cov(A₁,B₂) + Cov(A₁,B₃)
    cov(2, 1) + cov(2, 2) - cov(2, 3) +   ## Cov(A₂,B₁) + Cov(A₂,B₂) − Cov(A₂,B₃)
    cov(3, 1) - cov(3, 2)                  ## Cov(A₃,B₁) − Cov(A₃,B₂)
);

tpop = polyopt(p, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2);

result = cs_nctssos(tpop, solver_config);

@show result.objective
abs_error = abs(result.objective + 5.0)
@show abs_error
@assert abs_error < 1e-3
````

````
result.objective = -4.99999999982401
abs_error = 1.7598988932832071e-10

````

The quantum value $-5$ is recovered within $10^{-3}$ absolute tolerance
[pozsgay2017Covariance](@cite).

---
## Summary

| Problem | Algebra | Key function | Classical bound | Quantum bound |
|:--------|:--------|:-------------|:----------------|:--------------|
| Toy trace polynomial | `ProjectorAlgebra` ($P^2 = P$) | `tr` | — | −0.0312 (order 3) |
| CHSH (tracial) | `UnipotentAlgebra` ($U^2 = I$) | `tr` | 2 | $2\sqrt{2} \approx 2.828$ |
| Covariance Bell | `UnipotentAlgebra` ($U^2 = I$) | `ς(·) − ς(·)·ς(·)` | 4.5 | 5 |

Key API:
- `create_projector_variables`: create $P^2 = P$ operators
- `create_unipotent_variables`: create $U^2 = I$ operators
- `tr`: build tracial (maximally entangled) state polynomials
- `ς` / `varsigma`: build arbitrary-state polynomials
- `polyopt`: formulate optimization problem
- `SolverConfig`: configure solver and sparsity options
- `cs_nctssos`: solve via moment-SOS hierarchy

For linear (non-tracial) Bell inequalities, see the
Bell inequalities example.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

