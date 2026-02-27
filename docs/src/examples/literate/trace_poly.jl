# # [Tracial Polynomial Optimization](@id tracial-polynomial-optimization)
#
# Can quantum mechanics outperform classical physics in a *provable*,
# *quantitative* way?  Bell inequalities answer yes: they set sharp
# thresholds separating classical from quantum correlations.  Computing
# quantum bounds is a **tracial polynomial optimization** problem —
# minimize a polynomial in noncommutative operators under a
# [tracial state](@ref tracial-polynomial) (normalized trace).
# For a maximally entangled bipartite state
# $\psi_k = \tfrac{1}{\sqrt{k}}\sum_{i=1}^k |ii\rangle$,
# ``\langle\psi_k|X\otimes Y|\psi_k\rangle = \tfrac{1}{k}\operatorname{Tr}(X\,Y^{\mathsf T})``,
# so Bell-inequality optimization over such states reduces to tracial form
# [klep2022Optimization](@cite).
#
# This tutorial covers three problems of increasing complexity:
#
# 1. **Warm-up** — a toy projector ($P^2 = P$) problem to learn the API.
# 2. **CHSH inequality** — a two-party Bell test where quantum beats
#    classical by $\sqrt{2}$.
# 3. **Covariance Bell inequality** — a *nonlinear* generalization using
#    bias-corrected correlations.
#
# !!! note "Prerequisites"
#     [Trace polynomials](@ref tracial-polynomial),
#     [polynomial optimization](@ref polynomial-optimization), and the
#     [moment-SOHS hierarchy](@ref moment-sohs-hierarchy).
#     For the same inequalities in *state-polynomial* form (bipartite,
#     separate sites), see [Bell inequalities](@ref bell-inequalities).

# ## Setup
#
# We use [Mosek](@ref mosek) as the SDP solver, silenced via MathOptInterface.

using NCTSSoS, MosekTools
using JSON3

const TRACE_POLY_REFS = JSON3.read(read(joinpath(@__DIR__, "data", "trace_poly_refs.json"), String))

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true);
nothing #hide

# ---
# ## Warm-up: Three Projectors ($P^2 = P$)
#
# A [projector](@ref monoid-algebras-showcase) ($P^2 = P$) models a yes/no
# measurement in quantum mechanics.
# Given three projectors $P_1, P_2, P_3$, define
#
# ```math
# f = \operatorname{tr}(P_1 P_2 P_3)
#   + \operatorname{tr}(P_1 P_2)\,\operatorname{tr}(P_3).
# ```
#
# The first term is a three-body trace moment; the second mixes a pairwise
# moment with a single-body moment.  If the projectors commuted (a classical
# model), both terms would be nonnegative joint probabilities and $f \ge 0$.
# A negative value witnesses genuinely noncommuting structure.
#
# How negative can $f$ get over all quantum realizations?
# We approximate this minimum with the
# [moment-SOHS hierarchy](@ref moment-sohs-hierarchy): a sequence of
# [SDP](@ref semidefinite-programming) relaxations whose bounds tighten
# with increasing **order**.

# Each tuple `("P", 1:3)` creates one **label group** (one physical site):
# operators across groups commute; within a group they need not.
# The `registry` stores symbol ↔ index mappings and the constraint
# $P_i^2 = P_i$.

registry, (P,) = create_projector_variables([("P", 1:3)]);
nothing #hide

# [`tr`](@ref) wraps an operator word in the trace functional, producing a
# scalar trace moment.  Since the objective is purely tracial (no remaining
# operator word), we multiply by the identity `𝟙` to form a state
# polynomial that [`polyopt`](@ref) can optimize.

𝟙 = one(typeof(P[1]))
p = (tr(P[1] * P[2] * P[3]) + tr(P[1] * P[2]) * tr(P[3])) * 𝟙;
nothing #hide

# At **order 2** the SDP is small and fast:

spop = polyopt(p, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2);
result = cs_nctssos(spop, solver_config);

@show result.objective
@assert isapprox(result.objective, TRACE_POLY_REFS["toy_projector_order2_objective"], atol=1e-6)

# Raising to **order 3** tightens the bound:

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=3);
result = cs_nctssos(spop, solver_config);

@show result.objective
@assert isapprox(result.objective, TRACE_POLY_REFS["toy_projector_order3_objective"], atol=1e-6)

# Literature values: $-0.0467$ (order 2) and $-0.0312$ (order 3); our
# results match within $10^{-6}$ [klep2022Optimization](@cite).

# ---
# ## CHSH Bell Inequality (Tracial Form)
#
# The [CHSH inequality](@ref bell-inequalities) [clauser1969Proposed](@cite)
# is the most famous Bell test: Alice and Bob each choose one of two $\pm 1$
# measurements.  The CHSH score is
#
# ```math
# \mathcal{B}_{\text{CHSH}}
#   = \operatorname{tr}(A_1 B_1) + \operatorname{tr}(A_1 B_2)
#   + \operatorname{tr}(A_2 B_1) - \operatorname{tr}(A_2 B_2).
# ```
#
# Classically $\mathcal{B}_{\text{CHSH}} \leq 2$; quantum mechanics reaches
# the Tsirelson bound $2\sqrt{2} \approx 2.828$
# [cirelson1980Quantum](@cite).
# Any correlation above $2$ rules out a local hidden-variable (shared-randomness)
# model. We now recover the quantum maximum with NCTSSoS.
#
# ### The transpose trick
#
# In the tracial formulation, a maximally entangled state of local dimension
# $k$ converts bipartite expectations to single-system traces:
# ``\langle\psi_k|A\otimes B|\psi_k\rangle = \tfrac{1}{k}\operatorname{Tr}(A\,B^{\mathsf T})``.
# Each `Bⱼ` below therefore represents the *transposed* operator
# $B_j^{\mathsf T}$, and all four variables share **one label group** —
# they act on the same space and need not commute
# [klep2022Optimization](@cite).
#
# !!! tip "Tracial vs. state-polynomial formulation"
#     Separate label groups would enforce
#     ``[A_i, B_j^{\mathsf T}] = 0``, changing the feasible set.
#     For that bipartite formulation, see
#     [Bell inequalities](@ref bell-inequalities).

# Observables with $\pm 1$ outcomes satisfy $U^2 = I$
# ([UnipotentAlgebra](@ref monoid-algebras-showcase)):

registry, (vars,) = create_unipotent_variables([("v", 1:4)]);
A1, A2, B1, B2 = vars
# In code, `A1`/`A2` represent $A_1$/$A_2$ and `B1`/`B2` represent
# $B_1^{\mathsf T}$/$B_2^{\mathsf T}$ (transpose trick).
nothing #hide

# Negate for minimization ([`cs_nctssos`](@ref) minimizes):

𝟙 = one(typeof(A1))

p = -1.0 * (
    tr(A1 * B1) + tr(A1 * B2) +
    tr(A2 * B1) - tr(A2 * B2)
);
nothing #hide

# `MaximalElimination()` exploits [term sparsity](@ref term-sparsity)
# to decompose the SDP into smaller independent blocks:

tpop = polyopt(p * 𝟙, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=1, ts_algo=MaximalElimination());
result = cs_nctssos(tpop, solver_config);

@show result.objective
@assert isapprox(result.objective, -2 * sqrt(2), atol=1e-5)

# The optimum $-2\sqrt{2} \approx -2.828$ recovers the Tsirelson bound
# at the first hierarchy level.

# ---
# ## Covariance Bell Inequality (Bias-Corrected Correlations)
#
# CHSH uses *linear* correlations $\operatorname{tr}(A_i B_j)$.  When
# marginals are biased ($\operatorname{tr}(A_i) \neq 0$), the
# **covariance** removes that bias:
#
# ```math
# \operatorname{Cov}(A_i, B_j)
#   = \operatorname{tr}(A_i B_j)
#   - \operatorname{tr}(A_i)\,\operatorname{tr}(B_j).
# ```
#
# This involves a *product of two traces* — a nonlinear trace-moment
# expression.  Pozsgay *et al.* [pozsgay2017Covariance](@cite)
# constructed a covariance Bell inequality for three measurements per
# party with sign sets
# $S^+ = \{(1,1),(1,2),(1,3),(2,1),(2,2),(3,1)\}$ and
# $S^- = \{(2,3),(3,2)\}$:
#
# - **Classical bound:** $f \leq 4.5$
# - **Quantum bound:** $f = 5$
#
# See [Bell inequalities](@ref bell-inequalities) for an alternative
# implementation using [state polynomials](@ref state-polynomial).

# Six unipotent variables in a single label group (same transpose-trick
# reasoning as CHSH):

registry, (vars,) = create_unipotent_variables([("v", 1:6)]);
A = vars[1:3];  # Alice: A₁, A₂, A₃
B = vars[4:6];  # Bob:   B₁, B₂, B₃

𝟙 = one(typeof(A[1]));
nothing #hide

# Covariance helper — note the product of two trace moments:

cov(i, j) = 1.0 * tr(A[i] * B[j]) - 1.0 * tr(A[i]) * tr(B[j]);
nothing #hide

# Negate and solve at order 2:

p = -1.0 * (
    cov(1, 1) + cov(1, 2) + cov(1, 3) +   ## (1,·) ∈ S⁺
    cov(2, 1) + cov(2, 2) - cov(2, 3) +   ## (2,1),(2,2) ∈ S⁺; (2,3) ∈ S⁻
    cov(3, 1) - cov(3, 2)                  ## (3,1) ∈ S⁺; (3,2) ∈ S⁻
);
nothing #hide

tpop = polyopt(p * 𝟙, registry);

solver_config = SolverConfig(; optimizer=SILENT_MOSEK, order=2);
result = cs_nctssos(tpop, solver_config);

@show result.objective
abs_error = abs(result.objective + 5.0)
@show abs_error
@assert abs_error < 1e-3

# Quantum value $-5$ recovered within $10^{-3}$
# [pozsgay2017Covariance](@cite).  The same
# [`tr`](@ref) / [`polyopt`](@ref) / [`cs_nctssos`](@ref) pipeline handles
# both linear and nonlinear trace-moment objectives.

# ---
# ## Summary
#
# Three problems of increasing complexity, all solved by the same pipeline:
#
# | Problem | Algebra | Classical | Quantum |
# |:--------|:--------|:----------|:--------|
# | Toy trace polynomial | [ProjectorAlgebra](@ref monoid-algebras-showcase) ($P^2 = P$) | — | −0.0312 (order 3) |
# | CHSH (tracial) | [UnipotentAlgebra](@ref monoid-algebras-showcase) ($U^2 = I$) | 2 | $2\sqrt{2} \approx 2.828$ |
# | Covariance Bell | [UnipotentAlgebra](@ref monoid-algebras-showcase) ($U^2 = I$) | 4.5 | 5 |
#
# **API pipeline**:
# [`create_projector_variables`](@ref) /
# [`create_unipotent_variables`](@ref) →
# [`tr`](@ref) →
# [`polyopt`](@ref) →
# [`SolverConfig`](@ref) →
# [`cs_nctssos`](@ref)
#
# For the *state-polynomial* formulation (bipartite, separate sites), see
# [Bell inequalities](@ref bell-inequalities).
