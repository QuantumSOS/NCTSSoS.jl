# # [GNS Construction Interface Guide](@id gns-construction-guide)
#
# ## What is the GNS construction?
#
# A moment relaxation gives numerical expectation values ``L(w)`` for words
# ``w``. The GNS (Gelfand-Naimark-Segal) construction seeks matrices ``A_i``
# and a unit vector ``\xi`` such that
# ``\langle\xi,w(A)\xi\rangle = L(w)``.
#
# Solved moments need not admit this finite-dimensional representation.
# `gns(result)` checks the reconstruction before it returns operators.
# It can also return a feasible optimizer that attains the bound but does
# **not** reproduce all the moments. Read its verdict before using the model.
# These are numerical checks at stated tolerances, not rigorous certificates.
#
# ## Key concepts
#
# **Moment matrix (Hankel matrix).** For an ordered basis of words ``b_i``,
# its entries are
#
# ```math
# H_{ij} = L(b_i^\dagger b_j).
# ```
#
# **Relaxation order and Hankel degree.** A dense order-``d`` relaxation
# provides moments up to degree ``2d``. A Hankel matrix with row words up to
# degree ``k`` needs moments up to degree ``2k``. Here the relaxation has
# order 3, and the selected full Hankel has row degree 3.
#
# **Principal block and full Hankel.** `gns` selects an adjacent degree pair.
# `basis` contains the principal row words; `full_basis` contains the full
# row words. The principal numerical rank sets the operator dimension.
#
# **Flatness.** Equal principal and full ranks, with small extension
# residuals, mean that the higher-degree rows add no new directions. The
# returned diagnostics use the actual solved moments; the API does not
# replace their extension block to force equal ranks.
#
# ## API summary
#
# | Call or field | Purpose |
# |:---|:---|
# | `gns(result)` | Reconstruct and check a solved result |
# | `g.verdict`, `g.status`, `g.advice` | Read what was accepted or withheld |
# | `g[x]` | Read an accepted variable matrix |
# | `g.state` | Read the cyclic vector |
# | `evaluate(g, p)` | Evaluate a polynomial as a matrix |
# | `expectation(g, p)` | Evaluate a polynomial in the reconstructed state |
# | `g.flatness`, `g.checks` | Inspect numerical checks |
#
# **Prerequisites**: the [polynomial optimization API](@ref polynomial-optimization)
# and the [SDP relaxation workflow](@ref moment-sohs-hierarchy).

using NCTSSoS.Next, Mosek, LinearAlgebra

# ## Step 1 — Solve the noncommutative unit-ball problem
#
# We minimize
#
# ```math
# f = 2-X^2+XY^2X-Y^2
# ```
#
# subject to ``1-X^2-Y^2\succeq0``. The optimal value is 1.

registry, (vars,) = site_variables(Free(), 1, 2);
X, Y = vars;
f = 2.0 - X^2 + X * Y^2 * X - Y^2;
ball = 1.0 - X^2 - Y^2;
problem = Problem(f; constraints=(Positive(ball),));
result = solve(problem, Relaxation(order=3); solver=MosekBackend());
result.bound
@assert isapprox(result.bound, 1.0; atol=1e-5) #src

# ## Step 2 — Read the checked reconstruction
#
# The rank cutoff is explicit here. `check_atol` sets the absolute tolerance
# for matrix, moment and objective checks.

g = gns(result; atol=1e-6, check_atol=5e-5);
g.verdict

# This example returns `:optimizer_certified`, not `:moments_reproduced`.
# The operator model attains the bound and satisfies the ball constraint,
# but it does not reproduce every solved moment.

@assert g.verdict == :optimizer_certified #src

# The three possible verdicts have different meanings:
#
# | Verdict | Meaning |
# |:---|:---|
# | `:moments_reproduced` | A positive, flat pair passes the moment and model checks |
# | `:optimizer_certified` | The model satisfies the problem and attains the bound; moment reproduction is not certified |
# | `:withheld` | Operators and state are `nothing`; inspect `status`, `checks` and `advice` |
#
# Do not treat an optimizer-only verdict as a moment reconstruction.
# `expectation(result, p)` reads solved moments, while `expectation(g, p)`
# reads the returned model. They can differ.

# ## Step 3 — Inspect degrees, bases and ranks

(
    solved_order = g.relaxation_order,
    degrees = (g.principal_degree, g.full_degree),
    basis_sizes = (length(g.basis), length(g.full_basis)),
    ranks = (g.rank, g.full_rank),
)

# The selected degrees are 2 and 3, with 7 and 15 row words. The actual
# ranks are **5 → 9**, not 5 → 5. The full matrix is not flat.
# Earlier versions of this guide replaced the extension block before
# extraction. The checked API keeps the original moments and reports this
# difference instead.

@assert (g.principal_degree, g.full_degree) == (2, 3) #src
@assert (length(g.basis), length(g.full_basis)) == (7, 15) #src
@assert (g.rank, g.full_rank) == (5, 9) #src
@assert !g.flatness.is_flat #src

# Singular values show which directions the cutoff retains:

(principal=g.singular_values, full=g.full_singular_values)

# The flatness record includes positivity, cutoffs, singular-value gaps,
# the minimum eigenvalue, and range and extension residuals:

collect(pairs(g.flatness))

# ## Step 4 — Read operators and evaluate the objective
#
# The principal rank still gives a `5 × 5` model. Matrix entries depend on
# the SVD basis, but the objective and constraint checks do not.

A1 = g[X];
A2 = g[Y];
(A1=round.(A1, digits=3), A2=round.(A2, digits=3), state=round.(g.state, digits=6))

@assert size(A1) == size(A2) == (5, 5) #src
@assert norm(A1 - A1') <= 1e-8 #src
@assert isapprox(norm(g.state), 1.0; atol=1e-8) #src

F = evaluate(g, f);
(
    objective = expectation(g, f),
    objective_error = abs(expectation(g, f) - result.bound),
    constraint_min_eigenvalue = eigmin(Hermitian(evaluate(g, ball))),
)

@assert isapprox(expectation(g, f), 1.0; atol=1e-5) #src
@assert eigmin(Hermitian(evaluate(g, ball))) >= -5e-5 #src

# ## Step 5 — Separate objective checks from moment checks

collect(pairs(g.checks))

# `program_moment_error` checks all stored program moment components;
# `hankel_moment_error` checks the selected full Hankel. The program error
# is about 0.033 here, far above `5e-5`. The small objective error and the
# nonnegative constraint eigenvalue still support the optimizer verdict.

@assert g.checks.program_moment_error > 5e-5 #src
@assert g.checks.objective_error <= 5e-5 #src
@assert all(v -> v >= -5e-5, g.checks.constraint_min_eigenvalues) #src

# Conditioning is also reported: `sigma_min`, `sigma_max`, and
# `condition_number` describe the retained principal subspace.
# `distance_to_flat` and `distance_over_sigma_min` describe the departure
# from flatness. A small objective error alone says nothing about these.
#
# ## Tolerances and withheld models
#
# - `atol=nothing` selects a resolved singular-value gap above the numerical
#   tail. An explicit `atol` sets the absolute rank cutoff.
# - `rtol * σmax` is a lower bound on that cutoff (default `rtol=0`).
# - `check_atol` controls the absolute matrix, moment and objective errors.
# - `psd_rtol * max(1, σmax)` controls allowed full-Hankel negativity
#   (default `psd_rtol=1e-8`).
#
# If a result is withheld, inspect its singular values, residuals and advice.
# Improve solver accuracy or increase the relaxation order when needed.
# Do not raise tolerances solely to change the verdict.
#
# **See also:**
# - [GNS Optimizer Extraction](@ref gns-optimizer-extraction) — a textbook
#   example with the noncommutative polydisc.
# - [Recovering the Textbook CHSH Optimizer](@ref chsh-gns-reconstruction) —
#   a flat model with explicit gauge fixing.
