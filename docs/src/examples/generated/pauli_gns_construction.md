```@meta
EditURL = "../literate/pauli_gns_construction.jl"
```

# [Recovering Textbook Pauli Matrices and the XXX Singlet from GNS](@id pauli-gns-construction)

The GNS construction reconstructs operators only **up to unitary change of
basis**. That is the right mathematics, but it hides the familiar physics.
For the two-qubit XXX model we can do better: solve the moment relaxation,
run GNS, and then **fix the basis explicitly** so the reconstructed operators
become the standard matrices

```math
\sigma_x \otimes I, \quad \sigma_y \otimes I, \quad \sigma_z \otimes I,
\quad
I \otimes \sigma_x, \quad I \otimes \sigma_y, \quad I \otimes \sigma_z,
```

while the GNS cyclic vector becomes the singlet ground state

```math
|\psi_-\rangle = \frac{|01\rangle - |10\rangle}{\sqrt{2}}.
```

We use the two-qubit Heisenberg XXX Hamiltonian

```math
H = \frac{1}{4}\left(\sigma^x_1\sigma^x_2 + \sigma^y_1\sigma^y_2 + \sigma^z_1\sigma^z_2\right),
```

whose exact spectrum is one singlet level at ``-3/4`` and a triply-degenerate
triplet level at ``+1/4``.

**Prerequisites**: familiarity with the [polynomial optimization API](@ref polynomial-optimization)
and the [GNS construction interface](@ref gns-construction-guide).

Concretely, this page shows:

1. that the dense moment and SOS forms give matching bounds and moments,
2. the checked GNS reconstruction from each form,
3. an explicit gauge-fixing that recovers the textbook computational basis.

````julia
using NCTSSoS, Mosek, LinearAlgebra
import NCTSSoS
````

## Step 1 — Solve the order-2 relaxation in moment and SOS forms

Build one dense relaxation and solve it in both forms. The backend recovers
moments in either case, including the equality duals of the SOS form.
The user does not need to read JuMP equality multipliers.

````julia
registry, (σx, σy, σz) = pauli_variables(2);

ham = ComplexF64(1 / 4) * (
    σx[1] * σx[2] +
    σy[1] * σy[2] +
    σz[1] * σz[2]
);
built = NCTSSoS.build(Problem(ham), Relaxation(order=2));
moment_result = solve(built; solver=MosekBackend(form=:moment));
sos_result = solve(built; solver=MosekBackend(form=:sos));
````

Compare all 16 two-site Pauli basis operators, including the identity.
Every word in this algebra reduces to one of these, times a phase.

````julia
pauli_words = [a * b
    for a in (one(σx[1]), σx[1], σy[1], σz[1])
    for b in (one(σx[2]), σx[2], σy[2], σz[2])];
max_moment_recovery_error = maximum(
    abs(expectation(moment_result, w) - expectation(sos_result, w))
    for w in pauli_words
);

solve_summary = (
    moment_bound = moment_result.bound,
    sos_bound = sos_result.bound,
    words_checked = length(pauli_words),
    max_recovered_moment_error = max_moment_recovery_error,
)
solve_summary
````

````
(moment_bound = -0.7499999966623365, sos_bound = -0.749999999999977, words_checked = 16, max_recovered_moment_error = 4.44889936002113e-9)
````

Check bounds at `5e-6` and individual moments separately at `1e-8`.
An accurate bound alone does not imply accurate moments.

````julia
@assert isapprox(moment_result.bound, -0.75; atol=5e-6)
@assert isapprox(sos_result.bound, moment_result.bound; atol=5e-6)
@assert max_moment_recovery_error < 1e-8
````

## Step 2 — Checked GNS reconstruction from both forms

Each solve result goes directly to `gns`. Both must reproduce the moments
within `1e-6`. These checks are numerical, not rigorous certificates.

````julia
g_moment = gns(moment_result; atol=1e-6, check_atol=1e-6);
g = gns(sos_result; atol=1e-6, check_atol=1e-6);
(moment=g_moment.verdict, sos=g.verdict)
````

````
(moment = :moments_reproduced, sos = :moments_reproduced)
````

Compare verdicts, ranks and reconstructed expectations, not raw matrix
entries: the two GNS bases can differ by a unitary transformation.

````julia
@assert g_moment.verdict == g.verdict == :moments_reproduced
@assert g_moment.rank == g_moment.full_rank == g.rank == g.full_rank == 4
@assert maximum(abs(expectation(g_moment, w) - expectation(g, w)) for w in pauli_words) < 1e-6
````

Use the SOS result for the gauge alignment below. The raw operators form a
`4 × 4` Pauli representation, but not yet in the computational basis.

````julia
raw_X1 = g[σx[1]];
raw_Y1 = g[σy[1]];
raw_Z1 = g[σz[1]];
raw_X2 = g[σx[2]];
raw_Y2 = g[σy[2]];
raw_Z2 = g[σz[2]];

raw_summary = (
    rank = g.rank,
    full_rank = g.full_rank,
    xi = round.(g.state, digits=6),
    σx₁ = round.(raw_X1, digits=3),
    σz₁ = round.(raw_Z1, digits=3),
)
raw_summary
````

````
(rank = 4, full_rank = 4, xi = ComplexF64[-0.0 - 0.0im, -0.0 + 0.0im, 0.0 - 0.0im, -1.0 - 0.0im], σx₁ = ComplexF64[-0.0 + 0.0im -0.027 + 0.715im -0.009 - 0.633im -0.296 + 0.0im; -0.027 - 0.715im 0.036 - 0.0im 0.02 + 0.293im -0.631 - 0.054im; -0.009 + 0.633im 0.02 - 0.293im -0.036 - 0.0im -0.708 - 0.098im; -0.296 - 0.0im -0.631 + 0.054im -0.708 + 0.098im 0.0 + 0.0im], σz₁ = ComplexF64[0.0 - 0.0im -0.032 - 0.092im -0.048 - 0.503im 0.858 + 0.0im; -0.032 + 0.092im -0.105 - 0.0im -0.057 - 0.849im -0.505 - 0.002im; -0.048 + 0.503im -0.057 + 0.849im 0.105 - 0.0im 0.095 - 0.023im; 0.858 - 0.0im -0.505 + 0.002im 0.095 + 0.023im -0.0 - 0.0im])
````

The matrices look unfamiliar, but they already satisfy the Pauli algebra.

````julia
I4 = Matrix{ComplexF64}(I, 4, 4)
raw_representation_errors = (
    σx₁² = opnorm(raw_X1 * raw_X1 - I4),
    σy₁² = opnorm(raw_Y1 * raw_Y1 - I4),
    σz₁² = opnorm(raw_Z1 * raw_Z1 - I4),
    anticommutator_σx₁σy₁ = opnorm(raw_X1 * raw_Y1 + raw_Y1 * raw_X1),
    commutator_σz₁σz₂ = opnorm(raw_Z1 * raw_Z2 - raw_Z2 * raw_Z1),
)
raw_representation_errors
````

````
(σx₁² = 6.590489063477322e-13, σy₁² = 6.590478726880507e-13, σz₁² = 6.59366876213847e-13, anticommutator_σx₁σy₁ = 7.774285125744733e-16, commutator_σz₁σz₂ = 5.932219962321572e-16)
````

This really is a four-dimensional representation with a normalized cyclic
vector; it is just written in an inconvenient basis.

````julia
@assert g.rank == 4
@assert g.full_rank == 4
@assert isapprox(norm(g.state), 1.0; atol=1e-8)
@assert all(err -> err < 1e-6, values(raw_representation_errors))
````

The matrices above look nothing like the textbook Pauli matrices, but that is
just basis freedom. GNS only reconstructs the representation up to a unitary.

## Step 3 — Fix the computational basis

We now remove that gauge freedom in two steps.

1. `σz₁` and `σz₂` commute, so their joint eigenspaces define the four
   computational basis states `|00⟩`, `|01⟩`, `|10⟩`, `|11⟩`.
2. Each basis vector is still free up to a phase. We choose those phases so
   `σx₁` and `σx₂` become the standard bit-flip matrices.

After that, the `σy` matrices and the GNS cyclic vector fall into place.

````julia
function leading_projector_vector(P; atol=1e-6)
    F = eigen(Hermitian((P + P') / 2))
    idx = argmax(F.values)
    λ = F.values[idx]
    λ ≥ 1 - atol || error("Expected a rank-1 projector, got leading eigenvalue $λ.")
    v = F.vectors[:, idx]
    return v / norm(v)
end

function pauli_alignment_unitary(X1, X2, Z1, Z2, ξ; atol=1e-6)
    T = promote_type(eltype(X1), eltype(X2), eltype(Z1), eltype(Z2), eltype(ξ))
    I4 = Matrix{T}(I, 4, 4)
    labels = ((1, 1), (1, -1), (-1, 1), (-1, -1))

    U_joint = hcat([
        leading_projector_vector(((I4 + s1 * Z1) / 2) * ((I4 + s2 * Z2) / 2); atol=atol)
        for (s1, s2) in labels
    ]...)

    X1_joint = U_joint' * X1 * U_joint
    X2_joint = U_joint' * X2 * U_joint

    phases = Vector{ComplexF64}(undef, 4)
    phases[1] = 1.0 + 0im
    phases[2] = conj(X2_joint[1, 2]) / abs(X2_joint[1, 2])
    phases[3] = conj(X1_joint[1, 3]) / abs(X1_joint[1, 3])

    phase_from_X1 = conj(X1_joint[2, 4]) / abs(X1_joint[2, 4]) * phases[2]
    phase_from_X2 = conj(X2_joint[3, 4]) / abs(X2_joint[3, 4]) * phases[3]
    @assert isapprox(phase_from_X1, phase_from_X2; atol=100atol)
    phases[4] = phase_from_X1 / abs(phase_from_X1)

    U = U_joint * Diagonal(phases)

    ψ_singlet_ref = ComplexF64[0, 1 / sqrt(2), -1 / sqrt(2), 0]
    ξ_aligned = U' * ξ
    overlap = dot(ψ_singlet_ref, ξ_aligned)
    global_phase = abs(overlap) ≤ atol ? (1.0 + 0im) : overlap / abs(overlap)

    return U * global_phase
end

U = pauli_alignment_unitary(raw_X1, raw_X2, raw_Z1, raw_Z2, g.state)

X1 = U' * raw_X1 * U
Y1 = U' * raw_Y1 * U
Z1 = U' * raw_Z1 * U
X2 = U' * raw_X2 * U
Y2 = U' * raw_Y2 * U
Z2 = U' * raw_Z2 * U
ψ_gns = U' * g.state

joint_eigenbasis_summary = (
    σz₁ = round.(Z1, digits=6),
    σz₂ = round.(Z2, digits=6),
    ψ_gns = round.(ψ_gns, digits=6),
)
joint_eigenbasis_summary
````

````
(σz₁ = ComplexF64[1.0 + 0.0im 0.0 - 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 - 0.0im -0.0 - 0.0im; 0.0 - 0.0im 0.0 + 0.0im -1.0 + 0.0im 0.0 - 0.0im; 0.0 - 0.0im -0.0 - 0.0im 0.0 + 0.0im -1.0 + 0.0im], σz₂ = ComplexF64[1.0 - 0.0im -0.0 - 0.0im 0.0 - 0.0im 0.0 + 0.0im; -0.0 - 0.0im -1.0 - 0.0im 0.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -0.0 - 0.0im 1.0 - 0.0im -0.0 + 0.0im; 0.0 - 0.0im 0.0 + 0.0im 0.0 - 0.0im -1.0 - 0.0im], ψ_gns = ComplexF64[0.0 + 0.0im, 0.707107 + 0.0im, -0.707107 - 0.0im, -0.0 + 0.0im])
````

## Step 4 — Compare with the textbook Pauli matrices

Now we build the exact computational-basis matrices and check that the aligned
GNS operators agree with them numerically.

````julia
σx_ref = ComplexF64[0 1; 1 0]
σy_ref = ComplexF64[0 -im; im 0]
σz_ref = ComplexF64[1 0; 0 -1]
I2 = Matrix{ComplexF64}(I, 2, 2)

X1_ref = kron(σx_ref, I2)
Y1_ref = kron(σy_ref, I2)
Z1_ref = kron(σz_ref, I2)
X2_ref = kron(I2, σx_ref)
Y2_ref = kron(I2, σy_ref)
Z2_ref = kron(I2, σz_ref)

site_1_matrices = (
    σx₁ = round.(X1, digits=3),
    σy₁ = round.(Y1, digits=3),
    σz₁ = round.(Z1, digits=3),
)
site_1_matrices
````

````
(σx₁ = ComplexF64[0.0 - 0.0im -0.0 + 0.0im 1.0 + 0.0im -0.0 - 0.0im; -0.0 - 0.0im -0.0 - 0.0im -0.0 + 0.0im 1.0 + 0.0im; 1.0 - 0.0im -0.0 - 0.0im 0.0 + 0.0im 0.0 + 0.0im; -0.0 - 0.0im 1.0 - 0.0im 0.0 - 0.0im -0.0 + 0.0im], σy₁ = ComplexF64[0.0 - 0.0im -0.0 - 0.0im 0.0 - 1.0im -0.0 + 0.0im; -0.0 + 0.0im -0.0 + 0.0im 0.0 + 0.0im -0.0 - 1.0im; 0.0 + 1.0im 0.0 - 0.0im 0.0 - 0.0im -0.0 - 0.0im; 0.0 - 0.0im -0.0 + 1.0im -0.0 + 0.0im -0.0 - 0.0im], σz₁ = ComplexF64[1.0 + 0.0im 0.0 - 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im 0.0 - 0.0im -0.0 - 0.0im; 0.0 - 0.0im 0.0 + 0.0im -1.0 + 0.0im 0.0 - 0.0im; 0.0 - 0.0im -0.0 - 0.0im 0.0 + 0.0im -1.0 + 0.0im])
````

Site 2 should match the second-qubit Pauli operators just as cleanly.

````julia
site_2_matrices = (
    σx₂ = round.(X2, digits=3),
    σy₂ = round.(Y2, digits=3),
    σz₂ = round.(Z2, digits=3),
)
site_2_matrices
````

````
(σx₂ = ComplexF64[0.0 - 0.0im 1.0 + 0.0im -0.0 + 0.0im -0.0 - 0.0im; 1.0 - 0.0im -0.0 + 0.0im 0.0 + 0.0im -0.0 + 0.0im; -0.0 - 0.0im 0.0 - 0.0im -0.0 + 0.0im 1.0 + 0.0im; -0.0 - 0.0im 0.0 - 0.0im 1.0 - 0.0im -0.0 - 0.0im], σy₂ = ComplexF64[0.0 - 0.0im 0.0 - 1.0im -0.0 - 0.0im -0.0 + 0.0im; 0.0 + 1.0im 0.0 + 0.0im -0.0 + 0.0im -0.0 + 0.0im; -0.0 + 0.0im -0.0 - 0.0im -0.0 - 0.0im 0.0 - 1.0im; 0.0 - 0.0im -0.0 + 0.0im 0.0 + 1.0im -0.0 - 0.0im], σz₂ = ComplexF64[1.0 - 0.0im -0.0 - 0.0im 0.0 - 0.0im 0.0 + 0.0im; -0.0 - 0.0im -1.0 - 0.0im 0.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -0.0 - 0.0im 1.0 - 0.0im -0.0 + 0.0im; 0.0 - 0.0im 0.0 + 0.0im 0.0 - 0.0im -1.0 - 0.0im])
````

The actual numerical check is basis-independent: every aligned operator should
agree with its textbook reference up to tiny solver noise.

````julia
alignment_errors = (
    σx₁ = opnorm(X1 - X1_ref),
    σy₁ = opnorm(Y1 - Y1_ref),
    σz₁ = opnorm(Z1 - Z1_ref),
    σx₂ = opnorm(X2 - X2_ref),
    σy₂ = opnorm(Y2 - Y2_ref),
    σz₂ = opnorm(Z2 - Z2_ref),
)
alignment_errors
````

````
(σx₁ = 3.2935064941233763e-13, σy₁ = 3.2971354623823036e-13, σz₁ = 3.295462068981634e-13, σx₂ = 3.2936172282349783e-13, σy₂ = 3.2945129820172017e-13, σz₂ = 3.2962417400920433e-13)
````

Those errors should stay at the solver-noise level.

````julia
@assert all(err -> err < 1e-6, values(alignment_errors))
````

This is the point of the page: the raw GNS output was merely unitary-equivalent
to the usual Pauli representation, while the aligned output is the actual
textbook matrix model.

## Step 5 — Reconstruct the XXX Hamiltonian and its spectrum

With the aligned Pauli matrices in hand, we rebuild the Hamiltonian as an
ordinary `4 × 4` matrix.

````julia
H_gns = (X1 * X2 + Y1 * Y2 + Z1 * Z2) / 4
H_ref = (X1_ref * X2_ref + Y1_ref * Y2_ref + Z1_ref * Z2_ref) / 4

hamiltonian_summary = (
    H = round.(H_gns, digits=3),
    deviation_from_reference = opnorm(H_gns - H_ref),
)
hamiltonian_summary
````

````
(H = ComplexF64[0.25 + 0.0im -0.0 + 0.0im -0.0 - 0.0im -0.0 + 0.0im; -0.0 - 0.0im -0.25 + 0.0im 0.5 - 0.0im 0.0 - 0.0im; 0.0 + 0.0im 0.5 + 0.0im -0.25 - 0.0im -0.0 - 0.0im; -0.0 - 0.0im -0.0 - 0.0im 0.0 + 0.0im 0.25 + 0.0im], deviation_from_reference = 4.935774282274562e-13)
````

Diagonalizing the aligned Hamiltonian now exposes the familiar spectrum.

````julia
F = eigen(Hermitian((H_gns + H_gns') / 2))
spectrum = round.(F.values, digits=6)
spectrum
````

````
4-element Vector{Float64}:
 -0.75
  0.25
  0.25
  0.25
````

We also pin the exact spectrum numerically.

````julia
@assert opnorm(H_gns - H_ref) < 1e-6
@assert isapprox(F.values[1], -0.75; atol=1e-6)
@assert all(λ -> isapprox(λ, 0.25; atol=1e-6), F.values[2:4])
````

The GNS-reconstructed Hamiltonian has exactly the singlet/triplet structure we
expect: one ground-state energy `-3/4` and three excited states at `+1/4`.

## Step 6 — Recover the singlet ground state

Because the ground state is non-degenerate, the aligned GNS cyclic vector must
be that ground state.

````julia
ψ_singlet_ref = ComplexF64[0, 1 / sqrt(2), -1 / sqrt(2), 0]
ground_state_check = (
    ψ_gns = round.(ψ_gns, digits=6),
    ψ_singlet = round.(ψ_singlet_ref, digits=6),
    overlap = abs(dot(ψ_singlet_ref, ψ_gns)),
    eigen_residual = norm(H_gns * ψ_gns + (3 / 4) * ψ_gns),
)
ground_state_check
````

````
(ψ_gns = ComplexF64[0.0 + 0.0im, 0.707107 + 0.0im, -0.707107 - 0.0im, -0.0 + 0.0im], ψ_singlet = ComplexF64[0.0 + 0.0im, 0.707107 + 0.0im, -0.707107 + 0.0im, 0.0 + 0.0im], overlap = 1.0000000000000002, eigen_residual = 4.934800691006673e-13)
````

The aligned cyclic vector should agree with the singlet directly, not just up
to energy expectation value.

````julia
@assert norm(ψ_gns - ψ_singlet_ref) < 1e-6
@assert abs(real(dot(ψ_gns, H_gns * ψ_gns)) + 0.75) < 1e-6
@assert norm(H_gns * ψ_gns + (3 / 4) * ψ_gns) < 1e-6
````

The GNS vector is not just some optimizer in an abstract quotient space. After
fixing the basis, it is the familiar singlet

```math
\frac{|01\rangle - |10\rangle}{\sqrt{2}}.
```

## Summary

For the two-qubit XXX model, the dense order-2 relaxation gives:

1. matching moments from the moment and SOS forms, recovered by the backend,
2. a `4 × 4` operator representation of the Pauli algebra,
3. the exact textbook matrices after basis alignment,
4. the full XXX Hamiltonian matrix with spectrum `{-3/4, 1/4, 1/4, 1/4}`,
5. the singlet ground state as the aligned cyclic vector.

The important lesson is simple: **raw GNS output is only defined up to unitary
gauge, but for a small full-rank example you can fix that gauge explicitly and
recover the usual computational-basis physics — with either solver form.**

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

