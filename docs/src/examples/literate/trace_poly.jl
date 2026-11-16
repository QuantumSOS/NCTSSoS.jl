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
# 2. The **CHSH Bell inequality** — constructing its tracial objective.
# 3. A **covariance Bell inequality** — nonlinear objective involving products
#    of trace moments.
#
# **Prerequisites**: familiarity with
# [tracial polynomial concepts](@ref tracial-polynomial) and the
# [polynomial optimization API](@ref polynomial-optimization).

# ## Setup
#
# We use the native Mosek backend as the SDP solver.

using NCTSSoS.Next
import NCTSSoS, Mosek
using JSON3

const TRACE_POLY_REFS = JSON3.read(read(joinpath(pkgdir(NCTSSoS), "docs", "src", "examples", "literate", "data", "trace_poly_refs.json"), String))

const SOLVER = MosekBackend();
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
#
# One site contains three noncommuting projectors. The registry stores their
# algebra and names. Each polynomial carries that registry into `Problem`.
# `Next.RealTrace()` selects normalized real traces.

registry, (x,) = site_variables(Projector(), 1, 3; functional=Next.RealTrace());
nothing #hide

# #### Step 2 — Build the tracial objective
#
# `trace_symbol` builds a symbolic trace. Products of these symbols describe
# products of trace moments. No identity factor is needed.

p = trace_symbol(x[1] * x[2] * x[3]) + trace_symbol(x[1] * x[2]) * trace_symbol(x[3]);
nothing #hide

# #### Step 3 — Formulate and solve
#
# `Relaxation` sets the hierarchy `order`; `solve` selects the backend.
# Higher orders yield tighter bounds at the cost of larger SDPs.

spop = Problem(p);

result = solve(spop, Relaxation(order=2); solver=SOLVER);

@show result.bound
@assert isapprox(result.bound, TRACE_POLY_REFS["toy_projector_order2_objective"], atol=1e-6)

# #### Step 4 — Tighten the bound at order 3

result = solve(spop, Relaxation(order=3); solver=SOLVER);

@show result.bound
@assert isapprox(result.bound, TRACE_POLY_REFS["toy_projector_order3_objective"], atol=1e-6)

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
# We model the observables with `Unipotent()` ($U^2 = I$). Since
# `solve` *minimizes*, we negate the Bell expression and expect
# $\approx -2\sqrt{2}$.

# #### Step 1 — Create unipotent variables (single site)
#
# Variables at the same site do **not** commute. We place all four
# observables at one site, then split into Alice/Bob symbols.
#
# !!! note "Why non-commuting? The transpose trick"
#     In the tracial formulation, bipartite expectations over a maximally
#     entangled state are rewritten via the identity
#     ``\langle\phi^+|A\otimes B|\phi^+\rangle = \tfrac{1}{k}\operatorname{Tr}(A\,B^{\mathsf T})``
#     [klep2022Optimization](@cite).
#     When we write `trace_symbol(xᵢ * yⱼ)`, the variable `yⱼ` represents the
#     **transposed** Bob operator ``B_j^{\mathsf T}``. Placing Alice and Bob
#     at separate sites would impose ``[A_i, B_j^{\mathsf T}] = 0``,
#     which is **stronger** than the physical tensor-product commutation
#     ``[A_i \otimes I,\, I \otimes B_j] = 0``. Using a single site leaves
#     the variables non-commuting, matching the correct constraint set.
#     See the [Bell inequalities example](@ref bell-inequalities) for the
#     state-polynomial formulation, where separate sites *are* appropriate.

registry, (vars,) = site_variables(Unipotent(), 1, 4; functional=Next.RealTrace());
x = vars[1:2];  # Alice: A₁, A₂
y = vars[3:4];  # Bob:   B₁, B₂
nothing #hide

# #### Step 2 — Define the negated CHSH expression

p = -1.0 * trace_symbol(x[1] * y[1]) +  ## −tr(A₁B₁)
    -1.0 * trace_symbol(x[1] * y[2]) +  ## −tr(A₁B₂)
    -1.0 * trace_symbol(x[2] * y[1]) +  ## −tr(A₂B₁)
     1.0 * trace_symbol(x[2] * y[2]);   ## +tr(A₂B₂)
nothing #hide

# Term sparsity for state and trace polynomials is not available in this release.

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

# #### Step 1 — Create unipotent variables (single site)
#
# As with CHSH above, the tracial formulation requires non-commuting
# variables (see the note on the transpose trick). We place all six
# observables at one site and split afterward.

registry, (vars,) = site_variables(Unipotent(), 1, 6; functional=Next.RealTrace());
x = vars[1:3];  # Alice: A₁, A₂, A₃
y = vars[4:6];  # Bob:   B₁, B₂, B₃
nothing #hide

# #### Step 2 — Define the covariance helper

cov(a, b) = trace_symbol(x[a] * y[b]) - trace_symbol(x[a]) * trace_symbol(y[b]);
nothing #hide

# #### Step 3 — Build and solve the negated objective

p = -1.0 * (
    cov(1, 1) + cov(1, 2) + cov(1, 3) +   ## S⁺ terms
    cov(2, 1) + cov(2, 2) - cov(2, 3) +   ## mixed signs
    cov(3, 1) - cov(3, 2)                  ## S⁻ term: −Cov(A₃,B₂)
);
nothing #hide

tpop = Problem(p);

result = solve(tpop, Relaxation(order=2); solver=SOLVER);

@show result.bound
abs_error = abs(result.bound + 5.0)
@show abs_error
@assert abs_error < 1e-3

# The quantum value $-5$ is recovered within $10^{-3}$
# [pozsgay2017Covariance](@cite).

# ---
# ## Summary
#
# | Problem | Algebra | Classical | Quantum |
# |:--------|:--------|:----------|:--------|
# | Toy trace polynomial | `Projector()` ($P^2 = P$) | — | −0.0312 (order 3) |
# | CHSH (tracial) | `Unipotent()` ($U^2 = I$) | 2 | $2\sqrt{2} \approx 2.828$ |
# | Covariance Bell | `Unipotent()` ($U^2 = I$) | 4.5 | 5 |
#
# **Key API**:
# - `site_variables` with `Projector()` or `Unipotent()`:
#   create operators with $P^2 = P$ or $U^2 = I$ constraints
# - `trace_symbol`: build trace moments (tracial state)
# - `Problem`: formulate the optimization problem
# - `Relaxation`: set the hierarchy order
# - `solve`: solve via the moment-SOS hierarchy
#
# For linear (non-tracial) Bell inequalities, see the
# [Bell inequalities example](@ref bell-inequalities).
