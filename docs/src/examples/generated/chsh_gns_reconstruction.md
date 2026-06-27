```@meta
EditURL = "../literate/chsh_gns_reconstruction.jl"
```

# [Recovering the Textbook CHSH Optimizer from GNS](@id chsh-gns-reconstruction)

The CHSH inequality is the standard Bell-inequality warm-up:

```math
S = A_1 B_1 + A_1 B_2 + A_2 B_1 - A_2 B_2,
```

with `A_i^2 = B_j^2 = I` and `[A_i, B_j] = 0`.
The quantum maximum is Tsirelson's bound `2\sqrt{2}`.

The [Bell inequalities](@ref bell-inequalities) page shows how to certify that
bound with `NCTSSoS.jl`. Here we take the next step: **recover an actual
optimizer** from the solved moments.

Concretely, this page shows how to:

1. solve a dense CHSH moment relaxation,
2. run a flat GNS reconstruction,
3. fix the unitary gauge explicitly,
4. recover the textbook two-qubit model
   ```math
   A_1 = \sigma_z \otimes I,
   \quad
   A_2 = \sigma_x \otimes I,
   \quad
   B_1 = I \otimes \frac{\sigma_z + \sigma_x}{\sqrt{2}},
   \quad
   B_2 = I \otimes \frac{\sigma_z - \sigma_x}{\sqrt{2}},
   ```
   together with the Bell state
   ```math
   |\Phi^+\rangle = \frac{|00\rangle + |11\rangle}{\sqrt{2}}.
   ```

**Prerequisites**: familiarity with the [Bell inequalities](@ref bell-inequalities)
example and the [GNS construction interface](@ref gns-construction-guide).

One subtlety matters up front: **order 1 already certifies Tsirelson's bound,
but it does not give a flat enough Hankel matrix for a closed finite-dimensional
GNS model**. In this example, the first clean dense reconstruction comes from an
**order-4** relaxation with a degree-3 Hankel and a degree-2 principal block.

````julia
using NCTSSoS, MosekTools, LinearAlgebra, Logging

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(
    Mosek.Optimizer,
    MOI.Silent() => true,
)

function blockdiag(A::AbstractMatrix, B::AbstractMatrix)
    T = promote_type(eltype(A), eltype(B))
    m, n = size(A)
    p, q = size(B)
    M = zeros(T, m + p, n + q)
    M[1:m, 1:n] = A
    M[m+1:end, n+1:end] = B
    return M
end

function involution_eigenspaces(H::AbstractMatrix; atol::Real=1e-8)
    F = eigen(Hermitian((H + H') / 2))
    idx_plus = findall(λ -> λ > 1 - atol, F.values)
    idx_minus = findall(λ -> λ < -1 + atol, F.values)
    return F.vectors[:, idx_plus], F.vectors[:, idx_minus]
end
````

## Step 1 — Solve a dense order-4 CHSH relaxation

We use unipotent variables, so the relations `A_i^2 = I`, `B_j^2 = I`, and the
bipartite commutation rules are built into the algebra.

We still optimize `-S`, because [`cs_nctssos`](@ref) minimizes.

````julia
registry, (A, B) = create_unipotent_variables([("A", 1:2), ("B", 1:2)]);

chsh = 1.0 * A[1] * B[1] +
       1.0 * A[1] * B[2] +
       1.0 * A[2] * B[1] -
       1.0 * A[2] * B[2];

pop = polyopt(-chsh, registry);

solver_config = SolverConfig(
    optimizer = SILENT_MOSEK,
    order = 4,
    cs_algo = NoElimination(),
    ts_algo = NoElimination(),
);

sparsity = compute_sparsity(pop, solver_config);
moment_problem = NCTSSoS.moment_relax(
    pop,
    sparsity.corr_sparsity,
    sparsity.cliques_term_sparsities,
);
moment_result = NCTSSoS.solve_moment_problem(moment_problem, SILENT_MOSEK);

solve_summary = (
    objective = moment_result.objective,
    tsirelson = -moment_result.objective,
    expected = 2 * sqrt(2),
    n_unique_moments = moment_result.n_unique_elements,
    solved_moments = length(moment_result.monomap),
)
solve_summary
````

````
(objective = -2.828427124682018, tsirelson = 2.828427124682018, expected = 2.8284271247461903, n_unique_moments = 101, solved_moments = 101)
````

The objective should match `-2√2` up to solver tolerance.

````julia
@assert isapprox(moment_result.objective, -2 * sqrt(2); atol=1e-8)
````

## Step 2 — Check the flat Hankel condition

We reconstruct on a degree-3 Hankel matrix with a degree-2 principal block.
Flatness at this pair is what makes the final `4 × 4` model close exactly.

````julia
using NCTSSoS: get_ncbasis

full_basis = get_ncbasis(registry, 3);
basis = get_ncbasis(registry, 2);
hankel = NCTSSoS.hankel_matrix(moment_result.monomap, full_basis);
flatness = test_flatness(hankel, full_basis, basis; atol=1e-8)
````

````
NCTSSoS.FlatnessResult(true, 4, 4, 2.2140277378315492e-11)
````

`flatness` should report equal numerical ranks for the full Hankel and the
principal block.

````julia
@assert flatness.is_flat
@assert flatness.rank_principal == 4
@assert flatness.rank_full == 4
````

## Step 3 — Run the GNS reconstruction

The flatness check tells us the degree-2 quotient already closes, so GNS should
recover a genuine `4 × 4` representation.

````julia
gns = with_logger(Logging.SimpleLogger(devnull, Logging.Error)) do
    gns_reconstruct(moment_result.monomap, registry, 3; hankel_deg=2, atol=1e-8)
end;

A1_raw = gns.matrices[registry[:A₁]];
A2_raw = gns.matrices[registry[:A₂]];
B1_raw = gns.matrices[registry[:B₁]];
B2_raw = gns.matrices[registry[:B₂]];

raw_summary = (
    rank = gns.rank,
    full_rank = gns.full_rank,
    xi = round.(gns.xi, digits=6),
    A₁ = round.(A1_raw, digits=3),
    B₁ = round.(B1_raw, digits=3),
)
raw_summary
````

````
(rank = 4, full_rank = 4, xi = [-0.0, 1.0, -0.0, 0.0], A₁ = [-0.0 -0.0 1.0 -0.0; -0.0 -0.0 -0.0 -1.0; 1.0 -0.0 -0.0 -0.0; -0.0 -1.0 -0.0 0.0], B₁ = [0.0 -0.0 -0.707 -0.707; -0.0 -0.0 0.707 -0.707; -0.707 0.707 0.0 0.0; -0.707 -0.707 0.0 0.0])
````

The raw matrices are already a valid CHSH model — just not in a friendly basis.

````julia
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
````

````
(A₁² = 2.1183753326918557e-11, A₂² = 2.116762248081691e-11, B₁² = 2.1184447245239306e-11, B₂² = 2.118462772283887e-11, comm_A₁B₁ = 1.0737410962186246e-11, comm_A₁B₂ = 1.0725531578241554e-11, comm_A₂B₁ = 1.0725420558495212e-11, comm_A₂B₂ = 1.0739631415612551e-11, anticomm_A = 1.5177535868931984e-11, anticomm_B = 1.520554749857501e-11)
````

`raw_errors` should stay at the solver-noise level.

````julia
@assert all(err -> err < 5e-8, values(raw_errors))

verification = verify_gns(
    gns,
    moment_result.monomap,
    registry;
    poly = -chsh,
    f_star = moment_result.objective,
    atol = 1e-5,
)
verification
````

````
NCTSSoS.VerificationReport(true, 1.0682343898338331e-11, 1.6981083206246694e-11, Float64[], true)
````

The generic verification suite should agree with the hand-checked algebra facts above.

````julia
@assert verification.is_symmetric
@assert verification.moment_max_error < 1e-8
@assert verification.objective_error < 1e-8
````

The CHSH operator should have spectrum `{-2√2, 0, 0, 2√2}`.

````julia
chsh_raw = A1_raw * B1_raw + A1_raw * B2_raw + A2_raw * B1_raw - A2_raw * B2_raw;
raw_spectrum = round.(eigvals(Hermitian((chsh_raw + chsh_raw') / 2)), digits=6)
raw_spectrum
````

````
4-element Vector{Float64}:
 -2.828427
  0.0
  0.0
  2.828427
````

## Step 4 — Fix Alice's gauge

GNS is only defined up to unitary change of basis. We remove that ambiguity in
two local steps.

First diagonalize `A₁`. Because `A₁² = I`, its `+1` and `-1` eigenspaces are
both two-dimensional. In that basis, `A₂` must be off-diagonal because it
anticommutes with `A₁`. A block-diagonal correction then turns `A₂` into
`σx ⊗ I`.

````julia
Vp, Vm = involution_eigenspaces(A1_raw; atol=1e-8)
U_a1 = hcat(Vp, Vm)

A1_split = U_a1' * A1_raw * U_a1
A2_split = U_a1' * A2_raw * U_a1
B1_split = U_a1' * B1_raw * U_a1
B2_split = U_a1' * B2_raw * U_a1
xi_split = U_a1' * gns.xi

alice_block = A2_split[1:2, 3:4]
U_a2 = blockdiag(Matrix{ComplexF64}(I, 2, 2), alice_block')

U_alice = U_a1 * U_a2

A1_alice = U_alice' * A1_raw * U_alice
A2_alice = U_alice' * A2_raw * U_alice
B1_alice = U_alice' * B1_raw * U_alice
B2_alice = U_alice' * B2_raw * U_alice
xi_alice = U_alice' * gns.xi

alice_summary = (
    A₁ = round.(A1_alice, digits=6),
    A₂ = round.(A2_alice, digits=6),
    B₁_offdiag = opnorm(B1_alice[1:2, 3:4]),
    B₂_offdiag = opnorm(B2_alice[1:2, 3:4]),
    xi = round.(xi_alice, digits=6),
)
alice_summary
````

````
(A₁ = ComplexF64[1.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im; -0.0 + 0.0im 1.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im; -0.0 + 0.0im 0.0 + 0.0im -1.0 + 0.0im 0.0 + 0.0im; -0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im -1.0 + 0.0im], A₂ = ComplexF64[-0.0 + 0.0im -0.0 + 0.0im 1.0 + 0.0im -0.0 + 0.0im; -0.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im -0.0 + 0.0im; -0.0 + 0.0im 1.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im], B₁_offdiag = 2.6823790630177104e-12, B₂_offdiag = 2.6801523132166027e-12, xi = ComplexF64[0.000236 + 0.0im, -0.707107 + 0.0im, 0.707107 + 0.0im, 0.000236 + 0.0im])
````

After Alice's gauge-fixing, the Bob operators should be block-diagonal with
identical `2 × 2` blocks, i.e. of the form `I ⊗ b_j`.

````julia
@assert opnorm(A1_alice - kron(ComplexF64[1 0; 0 -1], Matrix{ComplexF64}(I, 2, 2))) < 5e-8
@assert opnorm(A2_alice - kron(ComplexF64[0 1; 1 0], Matrix{ComplexF64}(I, 2, 2))) < 5e-8
@assert opnorm(B1_alice[1:2, 3:4]) < 5e-8
@assert opnorm(B2_alice[1:2, 3:4]) < 5e-8
````

## Step 5 — Fix Bob's gauge

The two Bob observables commute with Alice, so only their common `2 × 2` block
remains to be aligned.

Instead of matching `B₁` and `B₂` directly, form the Pauli pair

```math
Z_B = \frac{B_1 + B_2}{\sqrt{2}},
\qquad
X_B = \frac{B_1 - B_2}{\sqrt{2}}.
```

At the optimum, these are just `σz` and `σx` on Bob's qubit.

````julia
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
xi_aligned = U_total' * gns.xi

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
````

````
(A₁ = 2.4177994380540508e-11, A₂ = 1.7377504887508332e-11, B₁ = 2.3440019234320677e-11, B₂ = 2.34431598498403e-11)
````

These errors should collapse to roundoff once the gauge is fixed.

````julia
@assert all(err -> err < 5e-8, values(alignment_errors))
````

The operator gauge is now fixed exactly to the textbook CHSH observables.

````julia
aligned_summary = (
    A₁ = round.(A1_aligned, digits=3),
    A₂ = round.(A2_aligned, digits=3),
    B₁ = round.(B1_aligned, digits=3),
    B₂ = round.(B2_aligned, digits=3),
)
aligned_summary
````

````
(A₁ = ComplexF64[1.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im -0.0 + 0.0im; -0.0 + 0.0im 1.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im -0.0 + 0.0im -1.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im -1.0 + 0.0im], A₂ = ComplexF64[-0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im], B₁ = ComplexF64[0.707 + 0.0im 0.707 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.707 + 0.0im -0.707 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.707 + 0.0im 0.707 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.707 + 0.0im -0.707 + 0.0im], B₂ = ComplexF64[0.707 + 0.0im -0.707 + 0.0im -0.0 + 0.0im -0.0 + 0.0im; -0.707 + 0.0im -0.707 + 0.0im -0.0 + 0.0im -0.0 + 0.0im; -0.0 + 0.0im -0.0 + 0.0im 0.707 + 0.0im -0.707 + 0.0im; -0.0 + 0.0im -0.0 + 0.0im -0.707 + 0.0im -0.707 + 0.0im])
````

## Step 6 — Recover the Bell state

The cyclic vector should now be the maximally entangled state `|Φ⁺⟩`, up to a
harmless global phase.

````julia
phi_plus = ComplexF64[1 / sqrt(2), 0, 0, 1 / sqrt(2)]
phase = dot(phi_plus, xi_aligned)
xi_aligned *= conj(phase) / abs(phase)

state_summary = (
    xi = round.(xi_aligned, digits=6),
    phi_plus = round.(phi_plus, digits=6),
    overlap = abs(dot(phi_plus, xi_aligned)),
)
state_summary
````

````
(xi = ComplexF64[0.707107 + 0.0im, -0.0 - 0.0im, 0.0 + 0.0im, 0.707107 + 0.0im], phi_plus = ComplexF64[0.707107 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.707107 + 0.0im], overlap = 0.9999999999912621)
````

Global phase is gone; the cyclic vector is now the textbook Bell state.

````julia
@assert norm(xi_aligned - phi_plus) < 5e-8
````

The aligned model should still attain Tsirelson's bound.

````julia
chsh_aligned = A1_aligned * B1_aligned + A1_aligned * B2_aligned + A2_aligned * B1_aligned - A2_aligned * B2_aligned
final_check = (
    expectation = real(dot(xi_aligned, chsh_aligned * xi_aligned)),
    spectrum = round.(eigvals(Hermitian((chsh_aligned + chsh_aligned') / 2)), digits=6),
)
final_check
````

````
(expectation = 2.828427124641356, spectrum = [-2.828427, 0.0, 0.0, 2.828427])
````

The aligned state still saturates Tsirelson's bound.

````julia
@assert isapprox(final_check.expectation, 2 * sqrt(2); atol=1e-8)
````

## Summary

For CHSH, the dense order-4 relaxation gives a flat degree-3 Hankel matrix,
and GNS recovers a closed `4 × 4` optimizer.

After an explicit gauge-fixing,

- `A₁ = σz ⊗ I`,
- `A₂ = σx ⊗ I`,
- `B₁ = I ⊗ (σz + σx)/√2`,
- `B₂ = I ⊗ (σz - σx)/√2`,
- `ξ = |Φ⁺⟩`.

So the SDP moments do not just certify the Tsirelson bound — they recover the
standard two-qubit realization of the optimal CHSH strategy.

**See also**:
- [Bell inequalities](@ref bell-inequalities) — basic CHSH and `I_{3322}` workflows.
- [GNS Construction Guide](@ref gns-construction-guide) — dense GNS API and diagnostics.
- [GNS Construction for Pauli Operators](@ref pauli-gns-construction) — another explicit gauge-fixing example.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

