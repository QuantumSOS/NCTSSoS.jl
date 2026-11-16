# # [Recovering the Textbook CHSH Optimizer from GNS](@id chsh-gns-reconstruction)
#
# The CHSH inequality is the standard Bell-inequality warm-up:
#
# ```math
# S = A_1 B_1 + A_1 B_2 + A_2 B_1 - A_2 B_2,
# ```
#
# with `A_i^2 = B_j^2 = I` and `[A_i, B_j] = 0`.
# The quantum maximum is Tsirelson's bound `2\sqrt{2}`.
#
# The [Bell inequalities](@ref bell-inequalities) page shows how to certify that
# bound with `NCTSSoS.jl`. Here we take the next step: **recover an actual
# optimizer** from the solved moments.
#
# Concretely, this page shows how to:
#
# 1. solve a dense CHSH moment relaxation,
# 2. run a flat GNS reconstruction,
# 3. fix the unitary gauge explicitly,
# 4. recover the textbook two-qubit model
#    ```math
#    A_1 = \sigma_z \otimes I,
#    \quad
#    A_2 = \sigma_x \otimes I,
#    \quad
#    B_1 = I \otimes \frac{\sigma_z + \sigma_x}{\sqrt{2}},
#    \quad
#    B_2 = I \otimes \frac{\sigma_z - \sigma_x}{\sqrt{2}},
#    ```
#    together with the Bell state
#    ```math
#    |\Phi^+\rangle = \frac{|00\rangle + |11\rangle}{\sqrt{2}}.
#    ```
#
# **Prerequisites**: familiarity with the [Bell inequalities](@ref bell-inequalities)
# example and the [GNS construction interface](@ref gns-construction-guide).
#
# One subtlety matters up front: **order 1 already certifies Tsirelson's bound,
# but it does not give a flat enough Hankel matrix for a closed finite-dimensional
# GNS model**. In this example, the first clean dense reconstruction comes from an
# **order-4** relaxation with a degree-3 Hankel and a degree-2 principal block.

using NCTSSoS, Mosek, LinearAlgebra

function blockdiag(A::AbstractMatrix, B::AbstractMatrix)
    T = promote_type(eltype(A), eltype(B))
    m, n = size(A)
    p, q = size(B)
    M = zeros(T, m + p, n + q)
    M[1:m, 1:n] = A
    M[m+1:end, n+1:end] = B
    return M
end
nothing #hide

function involution_eigenspaces(H::AbstractMatrix; atol::Real=1e-8)
    F = eigen(Hermitian((H + H') / 2))
    idx_plus = findall(λ -> λ > 1 - atol, F.values)
    idx_minus = findall(λ -> λ < -1 + atol, F.values)
    return F.vectors[:, idx_plus], F.vectors[:, idx_minus]
end
nothing #hide

# ## Step 1 — Solve a dense order-4 CHSH relaxation
#
# We use unipotent variables, so the relations `A_i^2 = I`, `B_j^2 = I`, and the
# bipartite commutation rules are built into the algebra.
#
# We optimize `-S`, because `solve` minimizes.

registry, (A, B) = site_variables(Unipotent(), 2, 2);

chsh = 1.0 * A[1] * B[1] +
       1.0 * A[1] * B[2] +
       1.0 * A[2] * B[1] -
       1.0 * A[2] * B[2];

problem = Problem(-chsh);
result = solve(problem, Relaxation(order=4); solver=MosekBackend());

solve_summary = (
    objective = result.bound,
    tsirelson = -result.bound,
    expected = 2 * sqrt(2),
    sizes = result.sizes,
)
solve_summary

# The objective should match `-2√2` up to solver tolerance.

@assert isapprox(result.bound, -2 * sqrt(2); atol=1e-8)

# ## Step 2 — Check the flat Hankel condition
#
# `gns` selects the degree pair and checks the reconstruction. Here it selects
# a degree-3 full Hankel with a degree-2 principal block. Flatness at this pair
# supports a closed `4 × 4` model, within the stated numerical tolerances.

g = gns(result; atol=1e-8)

# The verdict reports reproduced moments. The full and principal ranks agree.
# This is a numerical check, not a rigorous certificate.

@assert g.verdict == :moments_reproduced
@assert g.flatness.is_flat
@assert (g.principal_degree, g.full_degree) == (2, 3)
@assert g.rank == g.full_rank == 4

# ## Step 3 — Read the GNS operators
#
# Index the checked result with a variable to obtain its operator matrix.
# `g.state` is the cyclic vector.

A1_raw = g[A[1]];
A2_raw = g[A[2]];
B1_raw = g[B[1]];
B2_raw = g[B[2]];

raw_summary = (
    rank = g.rank,
    full_rank = g.full_rank,
    xi = round.(g.state, digits=6),
    A₁ = round.(A1_raw, digits=3),
    B₁ = round.(B1_raw, digits=3),
)
raw_summary

# The raw matrices are already a valid CHSH model — just not in a friendly basis.

I4 = Matrix{Float64}(I, 4, 4)
raw_errors = (
    A₁² = opnorm(A1_raw * A1_raw - I4),
    A₂² = opnorm(A2_raw * A2_raw - I4),
    B₁² = opnorm(B1_raw * B1_raw - I4),
    B₂² = opnorm(B2_raw * B2_raw - I4),
    comm_A₁B₁ = opnorm(A1_raw * B1_raw - B1_raw * A1_raw),
    comm_A₁B₂ = opnorm(A1_raw * B2_raw - B2_raw * A1_raw),
    comm_A₂B₁ = opnorm(A2_raw * B1_raw - B1_raw * A2_raw),
    comm_A₂B₂ = opnorm(A2_raw * B2_raw - B2_raw * A2_raw),
    anticomm_A = opnorm(A1_raw * A2_raw + A2_raw * A1_raw),
    anticomm_B = opnorm(B1_raw * B2_raw + B2_raw * B1_raw),
)
raw_errors

# `raw_errors` should stay at the solver-noise level.

@assert all(err -> err < 5e-8, values(raw_errors))

g.checks

# The built-in checks agree with the hand-checked algebra facts above.

@assert g.checks.algebra_error < 5e-8
@assert g.checks.program_moment_error < 1e-8
@assert g.checks.objective_error < 1e-8

# The CHSH operator should have spectrum `{-2√2, 0, 0, 2√2}`.

chsh_raw = A1_raw * B1_raw + A1_raw * B2_raw + A2_raw * B1_raw - A2_raw * B2_raw;
raw_spectrum = round.(eigvals(Hermitian((chsh_raw + chsh_raw') / 2)), digits=6)
raw_spectrum

# ## Step 4 — Fix Alice's gauge
#
# GNS is only defined up to unitary change of basis. We remove that ambiguity in
# two local steps.
#
# First diagonalize `A₁`. Because `A₁² = I`, its `+1` and `-1` eigenspaces are
# both two-dimensional. In that basis, `A₂` must be off-diagonal because it
# anticommutes with `A₁`. A block-diagonal correction then turns `A₂` into
# `σx ⊗ I`.

Vp, Vm = involution_eigenspaces(A1_raw; atol=1e-8)
U_a1 = hcat(Vp, Vm)

A1_split = U_a1' * A1_raw * U_a1
A2_split = U_a1' * A2_raw * U_a1
B1_split = U_a1' * B1_raw * U_a1
B2_split = U_a1' * B2_raw * U_a1
xi_split = U_a1' * g.state

alice_block = A2_split[1:2, 3:4]
U_a2 = blockdiag(Matrix{ComplexF64}(I, 2, 2), alice_block')

U_alice = U_a1 * U_a2

A1_alice = U_alice' * A1_raw * U_alice
A2_alice = U_alice' * A2_raw * U_alice
B1_alice = U_alice' * B1_raw * U_alice
B2_alice = U_alice' * B2_raw * U_alice
xi_alice = U_alice' * g.state

alice_summary = (
    A₁ = round.(A1_alice, digits=6),
    A₂ = round.(A2_alice, digits=6),
    B₁_offdiag = opnorm(B1_alice[1:2, 3:4]),
    B₂_offdiag = opnorm(B2_alice[1:2, 3:4]),
    xi = round.(xi_alice, digits=6),
)
alice_summary

# After Alice's gauge-fixing, the Bob operators should be block-diagonal with
# identical `2 × 2` blocks, i.e. of the form `I ⊗ b_j`.

@assert opnorm(A1_alice - kron(ComplexF64[1 0; 0 -1], Matrix{ComplexF64}(I, 2, 2))) < 5e-8
@assert opnorm(A2_alice - kron(ComplexF64[0 1; 1 0], Matrix{ComplexF64}(I, 2, 2))) < 5e-8
@assert opnorm(B1_alice[1:2, 3:4]) < 5e-8
@assert opnorm(B2_alice[1:2, 3:4]) < 5e-8

# ## Step 5 — Fix Bob's gauge
#
# The two Bob observables commute with Alice, so only their common `2 × 2` block
# remains to be aligned.
#
# Instead of matching `B₁` and `B₂` directly, form the Pauli pair
#
# ```math
# Z_B = \frac{B_1 + B_2}{\sqrt{2}},
# \qquad
# X_B = \frac{B_1 - B_2}{\sqrt{2}}.
# ```
#
# At the optimum, these are just `σz` and `σx` on Bob's qubit.

b1 = (B1_alice[1:2, 1:2] + B1_alice[3:4, 3:4]) / 2
b2 = (B2_alice[1:2, 1:2] + B2_alice[3:4, 3:4]) / 2

Zb = (b1 + b2) / sqrt(2)
Xb = (b1 - b2) / sqrt(2)

Fz = eigen(Hermitian((Zb + Zb') / 2))
perm = sortperm(Fz.values; rev=true)
Vb = Fz.vectors[:, perm]

Xb_diag = Vb' * Xb * Vb
phase = Xb_diag[1, 2] / abs(Xb_diag[1, 2])
Vb = Vb * Diagonal(ComplexF64[1, conj(phase)])

U_bob = kron(Matrix{ComplexF64}(I, 2, 2), Vb)
U_total = U_alice * U_bob

A1_aligned = U_total' * A1_raw * U_total
A2_aligned = U_total' * A2_raw * U_total
B1_aligned = U_total' * B1_raw * U_total
B2_aligned = U_total' * B2_raw * U_total
xi_aligned = U_total' * g.state

σx = ComplexF64[0 1; 1 0]
σz = ComplexF64[1 0; 0 -1]
I2 = Matrix{ComplexF64}(I, 2, 2)

A1_ref = kron(σz, I2)
A2_ref = kron(σx, I2)
B1_ref = kron(I2, (σz + σx) / sqrt(2))
B2_ref = kron(I2, (σz - σx) / sqrt(2))

alignment_errors = (
    A₁ = opnorm(A1_aligned - A1_ref),
    A₂ = opnorm(A2_aligned - A2_ref),
    B₁ = opnorm(B1_aligned - B1_ref),
    B₂ = opnorm(B2_aligned - B2_ref),
)
alignment_errors

# These errors should collapse to roundoff once the gauge is fixed.

@assert all(err -> err < 5e-8, values(alignment_errors))

# The operator gauge is now fixed exactly to the textbook CHSH observables.

aligned_summary = (
    A₁ = round.(A1_aligned, digits=3),
    A₂ = round.(A2_aligned, digits=3),
    B₁ = round.(B1_aligned, digits=3),
    B₂ = round.(B2_aligned, digits=3),
)
aligned_summary

# ## Step 6 — Recover the Bell state
#
# The cyclic vector should now be the maximally entangled state `|Φ⁺⟩`, up to a
# harmless global phase.

phi_plus = ComplexF64[1 / sqrt(2), 0, 0, 1 / sqrt(2)]
phase = dot(phi_plus, xi_aligned)
xi_aligned *= conj(phase) / abs(phase)

state_summary = (
    xi = round.(xi_aligned, digits=6),
    phi_plus = round.(phi_plus, digits=6),
    overlap = abs(dot(phi_plus, xi_aligned)),
)
state_summary

# Global phase is gone; the cyclic vector is now the textbook Bell state.

@assert norm(xi_aligned - phi_plus) < 5e-8

# The aligned model should still attain Tsirelson's bound.

chsh_aligned = A1_aligned * B1_aligned + A1_aligned * B2_aligned + A2_aligned * B1_aligned - A2_aligned * B2_aligned
final_check = (
    expectation = real(dot(xi_aligned, chsh_aligned * xi_aligned)),
    spectrum = round.(eigvals(Hermitian((chsh_aligned + chsh_aligned') / 2)), digits=6),
)
final_check

# The aligned state still saturates Tsirelson's bound.

@assert isapprox(final_check.expectation, 2 * sqrt(2); atol=1e-8)

# ## Summary
#
# For CHSH, the dense order-4 relaxation gives a flat degree-3 Hankel matrix,
# and GNS recovers a closed `4 × 4` optimizer.
#
# After an explicit gauge-fixing,
#
# - `A₁ = σz ⊗ I`,
# - `A₂ = σx ⊗ I`,
# - `B₁ = I ⊗ (σz + σx)/√2`,
# - `B₂ = I ⊗ (σz - σx)/√2`,
# - `ξ = |Φ⁺⟩`.
#
# So the SDP moments do not just certify the Tsirelson bound — they recover the
# standard two-qubit realization of the optimal CHSH strategy.
#
# **See also**:
# - [Bell inequalities](@ref bell-inequalities) — basic CHSH and `I_{3322}` workflows.
# - [GNS Construction Guide](@ref gns-construction-guide) — dense GNS API and diagnostics.
# - [GNS Construction for Pauli Operators](@ref pauli-gns-construction) — another explicit gauge-fixing example.
