# # [Tracial Polynomial Optimization](@id tracial-polynomial-optimization)
#
# Can quantum mechanics really outperform classical physics in a
# *provable*, *quantitative* way?  Bell inequalities answer yes:
# they define sharp thresholds that separate classical correlations from
# quantum ones.  Computing these quantum bounds is, at heart, a
# **tracial polynomial optimization** problem — minimize a polynomial in
# noncommutative operators under a trace (expectation) functional
# [klep2022Optimization](@cite).
#
# In this tutorial we build up to that punchline in three stages:
#
# 1. **Warm-up** — a toy projector problem to learn the API in isolation.
# 2. **CHSH inequality** — the most celebrated Bell test, where quantum
#    mechanics beats the classical limit by a factor of $\sqrt{2}$.
# 3. **Covariance Bell inequality** — a *nonlinear* generalization that
#    shows the quantum advantage persists even when we strip away marginal
#    biases.
#
# Each section first tells the *physics story* (what we measure, what
# classical and quantum theory predict), then shows how NCTSSoS turns
# that story into code.
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
# ## Warm-up: Projector Trace Polynomial
#
# Before we tackle Bell inequalities, let us get comfortable with the API
# on a small, self-contained problem.
#
# **Physical picture.**  Projectors ($P^2 = P$) appear everywhere in
# quantum mechanics — they represent yes/no measurement outcomes.
# A trace polynomial in projectors, such as
#
# ```math
# f = \operatorname{tr}(P_1 P_2 P_3)
#   + \operatorname{tr}(P_1 P_2)\,\operatorname{tr}(P_3),
# ```
#
# measures how three projective measurements interact.  The first term
# captures a *joint* three-body correlation; the second is a product of
# two lower-order moments.  We want the minimum of $f$ over *all*
# possible quantum realizations of $P_1, P_2, P_3$.
#
# **Goal.**  Find the tightest lower bound on $f$ using the moment-SOS
# hierarchy, and see how the bound improves as we increase the
# relaxation order.

# We declare three projector variables.  Each tuple `("x", 1:3)` creates
# a **label group**: the string is a name prefix and the range gives the
# indices.  The `registry` stores symbol ↔ index mappings and algebra
# constraints, which [`polyopt`](@ref) needs later to build the SDP.

registry, (x,) = create_projector_variables([("x", 1:3)]);
nothing #hide

# Translating the math into code is direct: [`tr`](@ref) wraps a
# polynomial into a trace moment, and multiplying by the identity
# monomial `𝟙` produces the `NCStatePolynomial` that [`polyopt`](@ref) expects.

𝟙 = one(NormalMonomial{ProjectorAlgebra, UInt8})
p = (tr(x[1] * x[2] * x[3]) + tr(x[1] * x[2]) * tr(x[3])) * 𝟙;
nothing #hide

# At **order 2** the SDP is small and fast:

spop = polyopt(p, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2);
result = cs_nctssos(spop, solver_config);

@show result.objective
@assert isapprox(result.objective, -0.046717378455438933, atol=1e-6)

# Raising the order to **3** adds more moment constraints and tightens
# the bound:

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=3);
result = cs_nctssos(spop, solver_config);

@show result.objective
@assert isapprox(result.objective, -0.03124998978001017, atol=1e-6)

# The literature values are $-0.0467$ (order 2) and $-0.0312$ (order 3);
# our results match within $10^{-6}$ [klep2022Optimization](@cite).
# With the API under our belt, we are ready for physics.

# ---
# ## CHSH Bell Inequality (Tracial Form)
#
# ### The physics: a tale of two parties
#
# The CHSH inequality [clauser1969Proposed](@cite) is the most famous
# Bell test.  The setup is deceptively simple:
#
# > Alice and Bob share a bipartite state.  Each independently chooses
# > one of two measurements with outcomes $\pm 1$.
#
# Denote Alice's observables $A_1, A_2$ and Bob's $B_1, B_2$.
# The CHSH score is the combination
#
# ```math
# \mathcal{B}_{\text{CHSH}}
#   = \operatorname{tr}(A_1 B_1) + \operatorname{tr}(A_1 B_2)
#   + \operatorname{tr}(A_2 B_1) - \operatorname{tr}(A_2 B_2).
# ```
#
# What makes this remarkable:
# - **Classical prediction.**  If outcomes are governed by shared
#   randomness (local hidden variables), then
#   $\mathcal{B}_{\text{CHSH}} \leq 2$ — always.
# - **Quantum prediction.**  Tsirelson showed that entangled quantum
#   states can reach $\mathcal{B}_{\text{CHSH}} = 2\sqrt{2} \approx 2.828$,
#   violating the classical limit by a factor of $\sqrt{2}$
#   [tsirelson2007Extremal](@cite).
#
# This gap is not a rounding error — it is a *proof* that quantum
# correlations are fundamentally stronger than classical ones.
# Let us verify the Tsirelson bound with NCTSSoS.
#
# ### Encoding observables as code
#
# Observables with outcomes $\pm 1$ satisfy $A_i^2 = I$, which is
# exactly the `UnipotentAlgebra` constraint.  We create four variables
# — $A_1, A_2, B_1, B_2$ — in a **single label group**.
#
# !!! note "Why one label group? The transpose trick"
#     In the standard bipartite picture, correlations are expectations like
#     ``\langle A_i \otimes B_j\rangle``. In the tracial formulation, for a
#     maximally entangled state these become
#     ``\langle\phi^+|A\otimes B|\phi^+\rangle = \tfrac{1}{k}\operatorname{Tr}(A\,B^{\mathsf T})``
#     [klep2022Optimization](@cite).
#     When we write `tr(Aᵢ * Bⱼ)`, the symbol `Bⱼ` represents the
#     **transposed** Bob operator ``B_j^{\mathsf T}``. Placing Alice and Bob
#     in separate label groups would impose ``[A_i, B_j^{\mathsf T}] = 0``,
#     which is **stronger** than the physical tensor-product commutation
#     ``[A_i \otimes I,\, I \otimes B_j] = 0``. Using a single group leaves
#     the variables non-commuting, matching the intended relaxation.
#     See the [Bell inequalities example](@ref bell-inequalities) for the
#     state-polynomial formulation, where separate groups *are* appropriate.

registry, (vars,) = create_unipotent_variables([("v", 1:4)]);
A1, A2, B1, B2 = vars
nothing #hide

# The code variable `A1` is Alice's $A_1$, `B2` is Bob's $B_2$, and so
# on — the mapping to the math is one-to-one.  Now we write the
# **negated** CHSH expression (since [`cs_nctssos`](@ref) minimizes):

𝟙 = one(NormalMonomial{UnipotentAlgebra, UInt8})

p = -1.0 * (
    tr(A1 * B1) + tr(A1 * B2) +
    tr(A2 * B1) - tr(A2 * B2)
);
nothing #hide

# Finally, we solve.  Setting `ts_algo = MaximalElimination()` exploits
# term sparsity — the SDP decomposes into smaller blocks, which is
# faster while giving the same bound.

tpop = polyopt(p * 𝟙, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=1, ts_algo=MaximalElimination());
result = cs_nctssos(tpop, solver_config);

@show result.objective
@assert isapprox(result.objective, -2 * sqrt(2), atol=1e-5)

# The optimum is $-2\sqrt{2} \approx -2.828$ — the Tsirelson bound,
# recovered at the *first* level of the hierarchy.  Quantum mechanics
# really does beat the classical limit of $2$.

# ---
# ## Covariance Bell Inequality
#
# ### Beyond linear correlations
#
# CHSH captures *linear* correlations: each term is a single
# $\operatorname{tr}(A_i B_j)$.  But what if Alice's and Bob's
# individual marginals are biased — what if
# $\operatorname{tr}(A_i) \neq 0$?  A biased coin can fake correlation.
#
# **Covariance** strips away that marginal bias:
#
# ```math
# \operatorname{Cov}(A_i, B_j)
#   = \operatorname{tr}(A_i B_j)
#   - \operatorname{tr}(A_i)\,\operatorname{tr}(B_j).
# ```
#
# This is a *nonlinear* function of trace moments (it involves a
# product of two traces), making the optimization harder.
# Pozsgay *et al.* [pozsgay2017Covariance](@cite) constructed a
# covariance Bell inequality for three measurements per party:
#
# ```math
# f = \sum_{(i,j)\in S^+} \operatorname{Cov}(A_i, B_j)
#   - \sum_{(i,j)\in S^-} \operatorname{Cov}(A_i, B_j)
# ```
#
# where $S^+ = \{(1,1),(1,2),(1,3),(2,1),(2,2),(3,1)\}$ and
# $S^- = \{(2,3),(3,2)\}$.
#
# The result is striking:
# - **Classical bound:** $f \leq 4.5$.
# - **Quantum bound:** $f = 5$.
#
# Even after removing marginal bias, quantum correlations still exceed
# the classical limit.  Let us verify this with NCTSSoS.
#
# ### From math to code
#
# As with CHSH, we use `UnipotentAlgebra` for $\pm 1$ observables.
# The six variables split into Alice ($A_1, A_2, A_3$) and Bob
# ($B_1, B_2, B_3$), all in a single label group (same transpose-trick
# reasoning as above).

registry, (vars,) = create_unipotent_variables([("v", 1:6)]);
A = vars[1:3];  # Alice: A₁, A₂, A₃
B = vars[4:6];  # Bob:   B₁, B₂, B₃
nothing #hide

𝟙 = one(typeof(A[1]));
nothing #hide

# The covariance helper translates the formula directly — note the
# product of two trace moments in the second term:

cov(i, j) = 1.0 * tr(A[i] * B[j]) - 1.0 * tr(A[i]) * tr(B[j]);
nothing #hide

# We negate the objective (minimization) and solve at order 2:

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
# [pozsgay2017Covariance](@cite).  NCTSSoS handles the nonlinear
# trace-moment products without any special treatment — the same
# `tr` / `polyopt` / `cs_nctssos` pipeline works for both linear
# and nonlinear Bell inequalities.

# ---
# ## Summary
#
# We told three stories of increasing complexity, all solved by the
# same pipeline:
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
