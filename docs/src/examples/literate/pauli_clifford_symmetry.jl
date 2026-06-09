# # [Pauli Symmetry Reduction: Manual Clifford Gates and Automatic SympleQ Detection](@id pauli-clifford-symmetry)
#
# The [CHSH symmetry](@ref chsh-symmetry) example shows how signed permutations
# on unipotent (measurement) variables can shrink a Bell-inequality SDP. This
# page covers the **Pauli** counterpart: Clifford gate conjugation on spin
# operators, the natural symmetry language for qubit Hamiltonians.
#
# The two features demonstrated here are:
#
# 1. **Manual [`CliffordSymmetry`](@ref)**: you specify which Clifford gate
#    leaves the Hamiltonian invariant (e.g. a SWAP gate on two qubits).
# 2. **Automatic SympleQ detection**: [`sympleq_symmetry_spec`](@ref) builds a
#    binary symplectic tableau from a Pauli Hamiltonian, discovers its
#    colour-preserving graph automorphisms, and hands back a ready-to-use
#    [`SymmetrySpec`](@ref) — no manual bookkeeping required.
#
# **Prerequisites:** familiarity with the
# [CHSH symmetry](@ref chsh-symmetry) example (for the symmetry workflow) and
# the [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) manual page (for
# the architectural picture: what NCTSSoS does locally, what
# [`SymbolicWedderburn`](https://github.com/kalmarek/SymbolicWedderburn.jl)
# does, and what is not yet supported).

using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(
    Mosek.Optimizer,
    MOI.Silent() => true,
)
nothing #hide

# ## Step 1 — Build the 2-site Heisenberg model
#
# The isotropic Heisenberg Hamiltonian on two qubits is
#
# ```math
# H = \frac{1}{4}\bigl(\sigma_1^x \sigma_2^x
#   + \sigma_1^y \sigma_2^y
#   + \sigma_1^z \sigma_2^z\bigr).
# ```
#
# Its ground-state energy is ``-3/4``. We will certify that lower bound with a
# moment relaxation and then show that Clifford symmetry makes it cheaper.

registry, (σx, σy, σz) = create_pauli_variables(1:2);

# Pauli Hamiltonians use `ComplexF64` coefficients (the algebra's phase
# structure needs complex arithmetic internally, even for Hermitian operators):

H = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))

# Package the Hamiltonian for optimization — the Pauli commutation and
# anticommutation relations are encoded in the algebra type, so no explicit
# constraints are needed:

pop = polyopt(H, registry);

# Fix an explicit order-1 half-basis so every run below uses exactly the same
# monomials. The seven elements are: the identity, plus all single-site Pauli
# operators.

basis = [one(σx[1]), σx[1], σx[2], σy[1], σy[2], σz[1], σz[2]];
length(basis)

# ## Step 2 — Dense baseline (no symmetry)
#
# Without symmetry or sparsity reduction, the moment matrix is a single dense
# block indexed by all seven basis monomials.

dense_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
)

dense_result = cs_nctssos(pop, dense_config);

# The objective recovers the exact ground-state energy:

dense_result.objective

# Sanity-check against the analytical value:

abs(dense_result.objective - (-0.75))

# The dense layout produces a single `7×7` PSD block:

dense_result.moment_matrix_sizes

# with this many distinct moment variables:

dense_result.n_unique_moment_matrix_elements

# No symmetry was applied, so the report field is empty:

dense_result.symmetry === nothing

# ## Step 3 — Manual Clifford symmetry: the SWAP gate
#
# The 2-site Heisenberg Hamiltonian is invariant under qubit exchange. In Pauli
# language, the SWAP gate conjugates every single-site Pauli operator on qubit 1
# into the corresponding operator on qubit 2, and vice versa:
#
# ```math
# \text{SWAP}\;\sigma_i^a\;\text{SWAP}^\dagger = \sigma_j^a, \quad
# \{i,j\} = \{1,2\},\; a \in \{x,y,z\}.
# ```
#
# Since each term ``\sigma_1^a \sigma_2^a`` is symmetric under this exchange,
# ``H`` is invariant.
#
# [`CliffordSymmetry`](@ref) provides named constructors for the standard
# Clifford gates — `:H` (Hadamard), `:S` (phase), `:CNOT`, and `:SWAP`:

swap = CliffordSymmetry(:SWAP, 1, 2)

# Wrap the generator in a [`SymmetrySpec`](@ref). The package will enumerate
# the full group (here: order 2) and verify that the objective is truly
# invariant before solving:

manual_spec = SymmetrySpec(swap)

# Build a config identical to the dense one, except for the `symmetry` keyword:

manual_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = manual_spec,
)

manual_result = cs_nctssos(pop, manual_config);

# Same objective — the symmetry reduction is exact, not an additional relaxation:

manual_result.objective

# Difference from the dense baseline:

abs(manual_result.objective - dense_result.objective)

# ### Reading the [`SymmetryReport`](@ref)
#
# The `symmetry` field now contains a diagnostic summary of what the reduction
# did:

report = manual_result.symmetry

# The SWAP generator generates a group of order 2 (identity + swap):

report.group_order

# The original half-basis had 7 monomials:

(report.basis_half_size, report.basis_full_size)

# The dense `7×7` PSD block has been split into two independent blocks — one
# of size 4 (symmetric subspace) and one of size 3 (antisymmetric subspace):

report.psd_block_sizes

# That is the headline win: the solver now handles two smaller PSD constraints
# instead of one large one — same answer, less work.

# ## Step 4 — Automatic SympleQ detection
#
# For the 2-site Heisenberg model, writing down SWAP by hand is trivial. For
# larger or less symmetric Hamiltonians, manual enumeration becomes tedious and
# error-prone.
#
# [`sympleq_symmetry_spec`](@ref) automates the process. It builds a binary
# symplectic tableau from the Pauli polynomial, encodes commutation and
# coefficient structure as a coloured graph, enumerates colour-preserving
# graph automorphisms, lifts them to symplectic Clifford actions (including a
# GF(2) phase solve to ensure the signs work out), and returns only those
# generators that genuinely leave `H` invariant.

auto_spec = sympleq_symmetry_spec(H)

# How many independent Clifford generators did it find?

length(auto_spec.clifford_generators)

# Solve with the auto-detected symmetry:

auto_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = auto_spec,
)

auto_result = cs_nctssos(pop, auto_config);

# Same objective again:

auto_result.objective

# Same block structure as the manual SWAP case:

auto_report = auto_result.symmetry
auto_report.group_order

#

auto_report.psd_block_sizes

# The SympleQ path found the same symmetry automatically — no Clifford-gate
# bookkeeping required.

# ## Step 5 — Beyond two qubits: 4-site Heisenberg chain
#
# The real payoff of automatic detection shows up when the system gets larger.
# Consider a 4-site periodic Heisenberg chain:
#
# ```math
# H_4 = \frac{1}{4}\sum_{i=1}^{4}\sum_{a \in \{x,y,z\}}
#        \sigma_i^a\,\sigma_{i+1\,(\mathrm{mod}\,4)}^a.
# ```

N = 4
registry4, (σx4, σy4, σz4) = create_pauli_variables(1:N);

H4 = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (σx4, σy4, σz4) for i in 1:N
)

pop4 = polyopt(H4, registry4);

# SympleQ detects the full point-group symmetry of the ring automatically:

auto_spec4 = sympleq_symmetry_spec(H4)

# Number of independent Clifford generators for the 4-site ring:

length(auto_spec4.clifford_generators)

# Build an order-1 half-basis (identity + all single-site Paulis):

basis4 = [one(σx4[1]); σx4; σy4; σz4];
length(basis4)

# Solve with symmetry:

config4 = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis4,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = auto_spec4,
)

result4 = cs_nctssos(pop4, config4);

result4.objective

# The symmetry report shows the group order and the reduced PSD block sizes:

report4 = result4.symmetry
report4.group_order

#

report4.psd_block_sizes

# Compare: without symmetry, we would have a single ``13 \times 13`` PSD block.
# With symmetry, it splits into several smaller blocks — and the solver does
# correspondingly less work. On larger systems, this reduction is the
# difference between feasible and intractable.

# ## Summary
#
# | run | Hamiltonian | PSD block sizes | group order | objective |
# |:----|:------------|:----------------|:------------|:----------|
# | dense baseline | 2-site | `[7]` | — | ``-0.75`` |
# | manual `CliffordSymmetry(:SWAP, 1, 2)` | 2-site | `[4, 3]` | 2 | ``-0.75`` |
# | `sympleq_symmetry_spec(H)` | 2-site | `[4, 3]` | 2 | ``-0.75`` |
# | `sympleq_symmetry_spec(H₄)` | 4-site ring | *(auto-detected)* | *(auto-detected)* | lower bound |
#
# The pattern: same answer, smaller SDP. Manual Clifford gates give you
# control; SympleQ gives you automation. For production use on large Pauli
# Hamiltonians, start with `sympleq_symmetry_spec` and fall back to manual
# gates only if you need to enforce a specific subgroup.
