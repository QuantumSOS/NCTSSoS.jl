# # [Tracial Polynomial Optimization](@id tracial-polynomial-optimization)
#
# Tracial polynomial optimization minimizes polynomial expressions involving
# traces of noncommutative operators — a natural formulation for optimizing
# over quantum states via the moment-SOS hierarchy
# [klep2022Optimization](@cite).
#
# This example covers three problems of increasing complexity:
#
# 1. A **toy problem** with projector variables — minimal setup to learn the API.
# 2. The **CHSH Bell inequality** — recovering the Tsirelson bound $2\sqrt{2}$.
# 3. A **covariance Bell inequality** — nonlinear objective involving products
#    of trace moments.
#
# **Prerequisites**: familiarity with
# [tracial polynomial concepts](@ref tracial-polynomial) and the
# [polynomial optimization API](@ref polynomial-optimization).

# ## Setup
#
# We use Mosek as the SDP solver. The `MOI.Silent()` attribute suppresses
# solver output.

using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true);
nothing #hide

# ---
# ## Toy Example: Projector Trace Polynomial
#
# Minimize a tracial polynomial over three projector variables $P_1, P_2, P_3$
# satisfying $P_i^2 = P_i$:
#
# ```math
# f = \operatorname{tr}(P_1 P_2 P_3)
#   + \operatorname{tr}(P_1 P_2)\,\operatorname{tr}(P_3)
# ```

# #### Step 1 — Create projector variables

registry, (x,) = create_projector_variables([("x", 1:3)]);
nothing #hide

# #### Step 2 — Build the tracial objective
#
# [`tr`](@ref) wraps a `Polynomial` into a `StatePolynomial` — a symbolic
# expression over trace moments. To feed it to [`polyopt`](@ref), multiply
# by the identity monomial `𝟙` to obtain an `NCStatePolynomial`.

𝟙 = one(NormalMonomial{ProjectorAlgebra, UInt8})
p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * 𝟙;
nothing #hide

# #### Step 3 — Formulate and solve
#
# [`SolverConfig`](@ref) sets the SDP backend and hierarchy `order`.
# Higher orders yield tighter bounds at the cost of larger SDPs.

spop = polyopt(p, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2);
result = cs_nctssos(spop, solver_config);

@show result.objective
@assert isapprox(result.objective, -0.046717378455438933, atol=1e-6)

# #### Step 4 — Tighten the bound at order 3

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=3);
result = cs_nctssos(spop, solver_config);

@show result.objective
@assert isapprox(result.objective, -0.03124998978001017, atol=1e-6)

# The literature values are $-0.0467$ (order 2) and $-0.0312$ (order 3); our
# results match within $10^{-6}$ [klep2022Optimization](@cite).

# ---
# ## CHSH Bell Inequality (Tracial Form)
#
# The CHSH inequality bounds correlations between two parties whose $\pm 1$
# observables are $A_1, A_2$ (Alice) and $B_1, B_2$ (Bob):
#
# ```math
# \mathcal{B}_{\text{CHSH}}
#   = \operatorname{tr}(A_1 B_1) + \operatorname{tr}(A_1 B_2)
#   + \operatorname{tr}(A_2 B_1) - \operatorname{tr}(A_2 B_2).
# ```
#
# Classical bound: $\mathcal{B} \leq 2$. Quantum bound (Tsirelson):
# $\mathcal{B} \leq 2\sqrt{2} \approx 2.828$.
#
# We model the observables with `UnipotentAlgebra` ($U^2 = I$). Since
# [`cs_nctssos`](@ref) *minimizes*, we negate the Bell expression and expect
# $\approx -2\sqrt{2}$.

# #### Step 1 — Create unipotent variables
#
# Separate label groups (`"x"` vs `"y"`) commute, encoding the bipartite
# locality assumption.

registry, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)]);
nothing #hide

# #### Step 2 — Define the negated CHSH expression

𝟙 = one(NormalMonomial{UnipotentAlgebra, UInt8})

p = -1.0 * tr(x[1] * y[1]) +  ## −tr(A₁B₁)
    -1.0 * tr(x[1] * y[2]) +  ## −tr(A₁B₂)
    -1.0 * tr(x[2] * y[1]) +  ## −tr(A₂B₁)
     1.0 * tr(x[2] * y[2]);   ## +tr(A₂B₂)
nothing #hide

# #### Step 3 — Solve with term-sparsity exploitation
#
# Setting `ts_algo = MaximalElimination()` decomposes the SDP into smaller
# blocks via connected components of the variable interaction graph.

tpop = polyopt(p * 𝟙, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=1, ts_algo=MaximalElimination());
result = cs_nctssos(tpop, solver_config);

@show result.objective
@assert isapprox(result.objective, -2 * sqrt(2), atol=1e-5)

# The result recovers the Tsirelson bound $-2\sqrt{2} \approx -2.828$
# [klep2022Optimization](@cite).

# ---
# ## Covariance Bell Inequality
#
# Nonlinear Bell inequalities involve products of trace moments. The
# covariance of observables $A_i, B_j$ is:
#
# ```math
# \operatorname{Cov}(A_i, B_j)
#   = \operatorname{tr}(A_i B_j)
#   - \operatorname{tr}(A_i)\,\operatorname{tr}(B_j).
# ```
#
# We maximize the covariance Bell expression from
# [pozsgay2017Covariance](@cite):
#
# ```math
# f = \sum_{(i,j)\in S^+} \operatorname{Cov}(A_i, B_j)
#   - \sum_{(i,j)\in S^-} \operatorname{Cov}(A_i, B_j)
# ```
#
# where $S^+ = \{(1,1),(1,2),(1,3),(2,1),(2,2),(3,1)\}$ and
# $S^- = \{(2,3),(3,2)\}$.
# Classical bound: $f \leq 4.5$.  Quantum bound: $f = 5$.

# #### Step 1 — Create unipotent variables (single group)
#
# Variables within the same label group do **not** commute. To keep all six
# observables non-commuting, we place them in one group and split afterward.

registry, (vars,) = create_unipotent_variables([("v", 1:6)]);
x = vars[1:3];  # Alice: A₁, A₂, A₃
y = vars[4:6];  # Bob:   B₁, B₂, B₃
nothing #hide

𝟙 = one(typeof(x[1]));
nothing #hide

# #### Step 2 — Define the covariance helper

cov(a, b) = 1.0 * tr(x[a] * y[b]) - 1.0 * tr(x[a]) * tr(y[b]);
nothing #hide

# #### Step 3 — Build and solve the negated objective

p = -1.0 * (
    cov(1, 1) + cov(1, 2) + cov(1, 3) +   ## S⁺ terms
    cov(2, 1) + cov(2, 2) - cov(2, 3) +   ## mixed signs
    cov(3, 1) - cov(3, 2)                  ## S⁻ term: −Cov(A₃,B₂)
);
nothing #hide

tpop = polyopt(p * 𝟙, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2);
result = cs_nctssos(tpop, solver_config);

@show result.objective
abs_error = abs(result.objective + 5.0)
@show abs_error
@assert abs_error < 1e-3

# The quantum value $-5$ is recovered within $10^{-3}$
# [pozsgay2017Covariance](@cite).

# ---
# ## Summary
#
# | Problem | Algebra | Classical | Quantum |
# |:--------|:--------|:----------|:--------|
# | Toy trace polynomial | `ProjectorAlgebra` ($P^2 = P$) | — | −0.0312 (order 3) |
# | CHSH (tracial) | `UnipotentAlgebra` ($U^2 = I$) | 2 | $2\sqrt{2} \approx 2.828$ |
# | Covariance Bell | `UnipotentAlgebra` ($U^2 = I$) | 4.5 | 5 |
#
# **Key API**:
# - [`create_projector_variables`](@ref) / [`create_unipotent_variables`](@ref):
#   create operators with $P^2 = P$ or $U^2 = I$ constraints
# - [`tr`](@ref): build trace moments (tracial state)
# - [`polyopt`](@ref): formulate the optimization problem
# - [`SolverConfig`](@ref): set solver backend, hierarchy order, and sparsity options
# - [`cs_nctssos`](@ref): solve via the moment-SOS hierarchy
#
# For linear (non-tracial) Bell inequalities, see the
# [Bell inequalities example](@ref bell-inequalities).
