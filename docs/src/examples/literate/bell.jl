# # [Bell inequalities](@id bell-inequalities)
#
# Bell inequalities test whether quantum mechanics can be explained by local hidden variable
# theories. They are linear combinations of expectation values with bounds that differ between
# classical and quantum theories.
#
# The general form of a Bell inequality is:
#
# ```math
# \sum_{i,j} c_{ij} \langle A_i B_j \rangle \leq C
# ```
#
# where $A_i$ and $B_j$ are observables measured by Alice and Bob, $c_{ij}$ are coefficients,
# and $C$ is the classical bound. Quantum mechanics can exceed this bound.

# ## Setup
#
# We use `NCTSSoS.jl` for polynomial optimization and `Mosek` as the SDP solver backend.
# Use `JuMPBackend(COSMO.Optimizer)` for an open-source alternative.

using NCTSSoS.Next
import Mosek

# ## Key Concepts: Unipotent and Projector Variables
#
# Bell inequalities use two types of measurement operators:
#
# 1. **Unipotent operators** ($U^2 = I$): Model ±1-valued observables like Pauli matrices
# 2. **Projector operators** ($P^2 = P$): Model projection measurements like $|0\rangle\langle 0|$
#
# Let's demonstrate both:

# ### Unipotent Variables (U² = I)
#
# Create operators that square to identity:

reg_unip, (A, B) = site_variables(Unipotent(), 2, 2)
# reg_unip: registry storing variable names and algebra type
# A: Alice's measurement operators [A₁, A₂] on site 1
# B: Bob's measurement operators [B₁, B₂] on site 2

# Verify the unipotent property (U² = I):
A[1] * A[1]  # should simplify to identity

# Check that operators on different sites commute:
A[1] * B[1] == B[1] * A[1]  # true: different sites commute

# ### Projector Variables (P² = P)
#
# Create operators that are idempotent:

reg_proj, (P, Q) = site_variables(Projector(), 2, 2)
# reg_proj: registry for projector algebra
# P: Alice's projectors [P₁, P₂] on site 1
# Q: Bob's projectors [Q₁, Q₂] on site 2

# Verify the idempotent property (P² = P):
P[1] * P[1] == P[1]  # true

# ---
# ## Linear Bell Inequalities
# ---

# ### CHSH Inequality
#
# The CHSH inequality involves two parties with two ±1-valued observables each.
# The objective function is:
#
# ```math
# f(A_1, A_2, B_1, B_2) = \langle A_1 B_1 \rangle + \langle A_1 B_2 \rangle + \langle A_2 B_1 \rangle - \langle A_2 B_2 \rangle
# ```
#
# Classical bound: $f \leq 2$. Quantum bound (Tsirelson): $f \leq 2\sqrt{2} \approx 2.828$.

# #### Step 1: Create unipotent variables for CHSH

registry, (x, y) = site_variables(Unipotent(), 2, 2)
# registry: variable registry encoding U² = I constraint
# x: Alice's observables [x₁, x₂] = [A₁, A₂]
# y: Bob's observables [y₁, y₂] = [B₁, B₂]

# #### Step 2: Define the CHSH objective function

f = 1.0 * x[1] * y[1] +  # ⟨A₁B₁⟩ term
    1.0 * x[1] * y[2] +  # ⟨A₁B₂⟩ term
    1.0 * x[2] * y[1] -  # ⟨A₂B₁⟩ term
    1.0 * x[2] * y[2]    # -⟨A₂B₂⟩ term
# f: polynomial representing the CHSH Bell operator

# Inspect the polynomial structure:
[(c, sprint(show, w; context=:registry => registry)) for (c, w) in Next.terms(f)]
# Each tuple contains a coefficient and its word, printed with operator names.

# #### Step 3: Create the optimization problem
#
# `Problem` creates a **minimization** problem. To find the maximum quantum
# violation we minimize `-f` and negate the result (same pattern as I₃₃₂₂ below).

pop = Problem(-f)
# pop: minimize -f ≡ maximize f, subject to algebraic constraints (U² = I)

# #### Step 4: Configure and run the SDP solver

relaxation = Relaxation(order=1)
# relaxation: specifies the hierarchy level

result = solve(pop, relaxation; solver=MosekBackend())
# result: optimization result containing objective value and solver info

# #### Step 5: Extract the upper bound

chsh_bound = -result.bound
# chsh_bound: upper bound on maximal quantum violation (negate since we minimized -f)

# Compare with Tsirelson's bound:
tsirelson_bound = 2 * sqrt(2)
# tsirelson_bound: theoretical maximum = 2√2 ≈ 2.828

abs(chsh_bound - tsirelson_bound)  # difference (should be ~1e-7)
@assert isapprox(chsh_bound, tsirelson_bound; atol=1e-7) #hide

# !!! tip "Going further: shrink this SDP with symmetry"
#     The CHSH operator is invariant under a 16-element symmetry group. The
#     [CHSH with Symmetry Reduction](@ref chsh-symmetry) example shows how to
#     use that group to replace the dense `5\times 5` PSD block with three
#     independent `1\times 1` PSD blocks while still recovering `2\sqrt{2}`.
#     The architectural picture is on the [Symmetry-Adapted Basis](@ref symmetry-adapted-basis)
#     manual page.

# ---
# ### $I_{3322}$ Inequality
#
# The $I_{3322}$ inequality uses **projector** observables (P² = P) with three measurements
# per party [pal2010maximal](@cite).
#
# ```math
# f = \langle A_1(B_1+B_2+B_3) \rangle + \langle A_2(B_1+B_2-B_3) \rangle + \langle A_3(B_1-B_2) \rangle - \langle A_1 \rangle - 2\langle B_1 \rangle - \langle B_2 \rangle
# ```
#
# Classical bound: $f \leq 0$. Quantum bound: $f \leq 0.25$.

# #### Step 1: Create projector variables

registry, (x, y) = site_variables(Projector(), 2, 3)
# registry: variable registry encoding P² = P constraint
# x: Alice's projectors [x₁, x₂, x₃] = [A₁, A₂, A₃]
# y: Bob's projectors [y₁, y₂, y₃] = [B₁, B₂, B₃]

# #### Step 2: Define the I₃₃₂₂ objective function

f = 1.0 * x[1] * (y[1] + y[2] + y[3]) +  # A₁(B₁+B₂+B₃)
    1.0 * x[2] * (y[1] + y[2] - y[3]) +  # A₂(B₁+B₂-B₃)
    1.0 * x[3] * (y[1] - y[2]) -         # A₃(B₁-B₂)
    1.0 * x[1] -                          # -A₁
    2.0 * y[1] -                          # -2B₁
    1.0 * y[2]                            # -B₂
# f: I₃₃₂₂ Bell polynomial

# Check the number of terms:
length(Next.terms(f))  # number of monomials

# #### Step 3: Solve (minimizing -f to find maximum of f)

pop = Problem(-f)
# pop: minimize -f (equivalent to maximize f)

relaxation = Relaxation(order=2)
# order=2: second level of the moment hierarchy

result = solve(pop, relaxation; solver=MosekBackend())
i3322_bound = -result.bound
# i3322_bound: upper bound on I₃₃₂₂ violation (negate since we minimized -f)

i3322_bound  # should be close to 0.25

# ---
# ### Exploiting Sparsity for Larger Problems
#
# Higher relaxation orders improve bounds but increase SDP size.
# **Sparsity exploitation** reduces computational cost:
#
# 1. **Correlative Sparsity (CS)**: Decomposes problem by variable interactions
# 2. **Term Sparsity (TS)**: Removes unnecessary monomials from moment matrices
#
# Let's solve I₃₃₂₂ at order=6 using correlative sparsity:

# #### Without sparsity (for comparison, order=3)

registry, (x, y) = site_variables(Projector(), 2, 3)
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]
pop = Problem(-f)

relaxation_dense = Relaxation(order=3)
# relaxation_dense: no sparsity exploitation

@time result_dense = solve(pop, relaxation_dense; solver=MosekBackend())
bound_dense = -result_dense.bound
# bound_dense: bound without sparsity

bound_dense

# #### With correlative sparsity (order=6)

relaxation_sparse = Relaxation(
    order = 6,             # higher order for better bound
    grouping = Correlative(MF())  # minimum-fill chordal completion
)
# Correlative(MF()): enables correlative sparsity via chordal graph decomposition

@time result_sparse = solve(pop, relaxation_sparse; solver=MosekBackend())
bound_sparse = -result_sparse.bound
# bound_sparse: improved bound using sparsity

bound_sparse  # closer to theoretical 0.25

# Improvement in bound:
bound_dense - bound_sparse  # positive = improvement

# ---
# ## Nonlinear Bell Inequalities
#
# Nonlinear Bell inequalities involve polynomial functions of expectation values,
# not just linear combinations. They can detect non-locality where linear inequalities fail.
#
# ### Covariance Bell Inequality
#
# The covariance between observables A and B is:
#
# ```math
# \text{Cov}(A, B) = \langle AB \rangle - \langle A \rangle \langle B \rangle
# ```
#
# This is **nonlinear** because it involves products of expectation values.
#
# The covariance Bell inequality [pozsgay2017Covariance](@cite):
#
# ```math
# f = \sum_{i,j} s_{ij} \text{Cov}(A_i, B_j)
# ```
#
# with signs $s_{ij} \in \{+1, -1\}$. Classical bound: $f \leq 4.5$. Quantum bound: $f = 5$.

# #### Step 1: Create unipotent variables

registry, (x, y) = site_variables(Unipotent(), 2, 3; functional=Next.RealState())
# x: Alice's observables [A₁, A₂, A₃]
# y: Bob's observables [B₁, B₂, B₃]

# #### Step 2: Inspect the identity
#
# The polynomial carries its algebra and functional type. Its identity belongs
# to the same registry. No identity factor is needed to build state polynomials.

ID = one(x[1])
# ID: identity monomial (𝟙) inferred from the variable type

ID  # display the identity

# #### Step 3: Define the covariance function using state polynomials
#
# `state_symbol` denotes a symbolic expectation value ⟨·⟩. It returns a
# polynomial in the same registry, with functional words for its terms.

cov(a, b) = state_symbol(x[a] * y[b]) -  # ⟨AᵢBⱼ⟩
            state_symbol(x[a]) * state_symbol(y[b])  # -⟨Aᵢ⟩⟨Bⱼ⟩
# cov(a,b): covariance Cov(Aₐ, Bᵦ) as a state polynomial

# Example: Cov(A₁, B₁)
cov(1, 1)

# #### Step 4: Build the objective function

sp = cov(1,1) + cov(1,2) + cov(1,3) +  # Cov(A₁, B₁) + Cov(A₁, B₂) + Cov(A₁, B₃)
     cov(2,1) + cov(2,2) - cov(2,3) +  # Cov(A₂, B₁) + Cov(A₂, B₂) - Cov(A₂, B₃)
     cov(3,1) - cov(3,2)               # Cov(A₃, B₁) - Cov(A₃, B₂)
# sp: state polynomial for covariance Bell inequality

# #### Step 5: Create optimization problem and solve

# Flipping all of Alice's observables changes `sp` to `-sp` and preserves the
# constraints. Thus minus the minimum of `sp` also bounds its maximum.
spop = Problem(sp)
# spop: state polynomial optimization problem

relaxation = Relaxation(order=2)

result = solve(spop, relaxation; solver=MosekBackend())
cov_bound = -result.bound
# cov_bound: upper bound on covariance Bell violation

cov_bound  # should be close to 5.0

# Compare with known quantum value:
abs(cov_bound - 5.0)  # difference from theoretical value

# Term sparsity for state and trace polynomials is not available in this release.

# ---
# ## Summary
#
# | Inequality | Operator Type | Classical Bound | Quantum Bound | API |
# |:-----------|:--------------|:----------------|:--------------|:----|
# | CHSH | Unipotent (U²=I) | 2 | 2√2 ≈ 2.828 | `site_variables(Unipotent(), 2, 2)` |
# | I₃₃₂₂ | Projector (P²=P) | 0 | 0.25 | `site_variables(Projector(), 2, 3)` |
# | Covariance | Unipotent + State | 4.5 | 5 | `state_symbol` |
#
# Key functions:
# - `site_variables`: Create operators with the selected algebra
# - `Problem`: Create optimization problem
# - `Relaxation`: Configure hierarchy order and sparsity options
# - `solve`: Solve using moment-SOS hierarchy
# - `state_symbol`: Symbolic expectation value for state polynomials
