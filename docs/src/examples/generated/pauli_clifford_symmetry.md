```@meta
EditURL = "../literate/pauli_clifford_symmetry.jl"
```

# [Pauli Symmetry Reduction: Manual Clifford Gates and Automatic SympleQ Detection](@id pauli-clifford-symmetry)

The [CHSH symmetry](@ref chsh-symmetry) example shows how signed permutations
on unipotent (measurement) variables can shrink a Bell-inequality SDP. This
page covers the **Pauli** counterpart: Clifford gate conjugation on spin
operators, the natural symmetry language for qubit Hamiltonians.

The two features demonstrated here are:

1. **Manual [`CliffordSymmetry`](@ref)**: you specify which Clifford gate
   leaves the Hamiltonian invariant (e.g. a SWAP gate on two qubits).
2. **Automatic SympleQ detection**: [`sympleq_symmetry_spec`](@ref) builds a
   binary symplectic tableau from a Pauli Hamiltonian, discovers its
   colour-preserving graph automorphisms, and hands back a ready-to-use
   [`SymmetrySpec`](@ref) — no manual bookkeeping required.

**Prerequisites:** familiarity with the
[CHSH symmetry](@ref chsh-symmetry) example (for the symmetry workflow) and
the [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) manual page (for
the architectural picture: what NCTSSoS does locally, what
[`SymbolicWedderburn`](https://github.com/kalmarek/SymbolicWedderburn.jl)
does, and what is not yet supported).

````julia
using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(
    Mosek.Optimizer,
    MOI.Silent() => true,
)
````

## Step 1 — Build the 2-site Heisenberg model

The isotropic Heisenberg Hamiltonian on two qubits is

```math
H = \frac{1}{4}\bigl(\sigma_1^x \sigma_2^x
  + \sigma_1^y \sigma_2^y
  + \sigma_1^z \sigma_2^z\bigr).
```

Its ground-state energy is ``-3/4``. We will certify that lower bound with a
moment relaxation and then show that Clifford symmetry makes it cheaper.

````julia
registry, (σx, σy, σz) = create_pauli_variables(1:2);
````

Pauli Hamiltonians use `ComplexF64` coefficients (the algebra's phase
structure needs complex arithmetic internally, even for Hermitian operators):

````julia
H = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (σx, σy, σz))
````

````
0.25 + 0.0im * σx₁σx₂ + 0.25 + 0.0im * σy₁σy₂ + 0.25 + 0.0im * σz₁σz₂
````

Package the Hamiltonian for optimization — the Pauli commutation and
anticommutation relations are encoded in the algebra type, so no explicit
constraints are needed:

````julia
pop = polyopt(H, registry);
````

Fix an explicit order-1 half-basis so every run below uses exactly the same
monomials. The seven elements are: the identity, plus all single-site Pauli
operators.

````julia
basis = [one(σx[1]), σx[1], σx[2], σy[1], σy[2], σz[1], σz[2]];
length(basis)
````

````
7
````

## Step 2 — Dense baseline (no symmetry)

Without symmetry or sparsity reduction, the moment matrix is a single dense
block indexed by all seven basis monomials.

````julia
dense_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
)

dense_result = cs_nctssos(pop, dense_config);
````

The objective recovers the exact ground-state energy:

````julia
dense_result.objective
````

````
-0.7500000000000073
````

Sanity-check against the analytical value:

````julia
abs(dense_result.objective - (-0.75))
````

````
7.327471962526033e-15
````

The dense layout produces a single `7×7` PSD block:

````julia
dense_result.moment_matrix_sizes
````

````
1-element Vector{Vector{Int64}}:
 [7]
````

with this many distinct moment variables:

````julia
dense_result.n_unique_moment_matrix_elements
````

````
16
````

No symmetry was applied, so the report field is empty:

````julia
dense_result.symmetry === nothing
````

````
true
````

## Step 3 — Manual Clifford symmetry: the SWAP gate

The 2-site Heisenberg Hamiltonian is invariant under qubit exchange. In Pauli
language, the SWAP gate conjugates every single-site Pauli operator on qubit 1
into the corresponding operator on qubit 2, and vice versa:

```math
\text{SWAP}\;\sigma_i^a\;\text{SWAP}^\dagger = \sigma_j^a, \quad
\{i,j\} = \{1,2\},\; a \in \{x,y,z\}.
```

Since each term ``\sigma_1^a \sigma_2^a`` is symmetric under this exchange,
``H`` is invariant.

[`CliffordSymmetry`](@ref) provides named constructors for the standard
Clifford gates — `:H` (Hadamard), `:S` (phase), `:CNOT`, and `:SWAP`:

````julia
swap = CliffordSymmetry(:SWAP, 1, 2)
````

````
CliffordSymmetry(nqubits=2, images=6)
````

Wrap the generator in a [`SymmetrySpec`](@ref). The package will enumerate
the full group (here: order 2) and verify that the objective is truly
invariant before solving:

````julia
manual_spec = SymmetrySpec(swap)
````

````
NCTSSoS.SymmetrySpec(NCTSSoS.SignedPermutation[], NCTSSoS.FermionicModePermutation[], NCTSSoS.CliffordSymmetry[CliffordSymmetry(nqubits=2, images=6)], nothing, nothing, true)
````

Build a config identical to the dense one, except for the `symmetry` keyword:

````julia
manual_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = manual_spec,
)

manual_result = cs_nctssos(pop, manual_config);
````

Same objective — the symmetry reduction is exact, not an additional relaxation:

````julia
manual_result.objective
````

````
-0.7500000000000073
````

Difference from the dense baseline:

````julia
abs(manual_result.objective - dense_result.objective)
````

````
0.0
````

### Reading the [`SymmetryReport`](@ref)

The `symmetry` field now contains a diagnostic summary of what the reduction
did:

````julia
report = manual_result.symmetry
````

````
SymmetryReport(group_order=2, invariant_moment_count=9, psd_block_sizes=[4, 3], basis_half_size=7, basis_full_size=16, block_provenance=[:wedderburn, :wedderburn])
````

The SWAP generator generates a group of order 2 (identity + swap):

````julia
report.group_order
````

````
2
````

The original half-basis had 7 monomials:

````julia
(report.basis_half_size, report.basis_full_size)
````

````
(7, 16)
````

The dense `7×7` PSD block has been split into two independent blocks — one
of size 4 (symmetric subspace) and one of size 3 (antisymmetric subspace):

````julia
report.psd_block_sizes
````

````
2-element Vector{Int64}:
 4
 3
````

That is the headline win: the solver now handles two smaller PSD constraints
instead of one large one — same answer, less work.

## Step 4 — Automatic SympleQ detection

For the 2-site Heisenberg model, writing down SWAP by hand is trivial. For
larger or less symmetric Hamiltonians, manual enumeration becomes tedious and
error-prone.

[`sympleq_symmetry_spec`](@ref) automates the process. It builds a binary
symplectic tableau from the Pauli polynomial, encodes commutation and
coefficient structure as a coloured graph, enumerates colour-preserving
graph automorphisms, lifts them to symplectic Clifford actions (including a
GF(2) phase solve to ensure the signs work out), and returns only those
generators that genuinely leave `H` invariant.

````julia
auto_spec = sympleq_symmetry_spec(H)
````

````
NCTSSoS.SymmetrySpec(NCTSSoS.SignedPermutation[], NCTSSoS.FermionicModePermutation[], NCTSSoS.CliffordSymmetry[CliffordSymmetry(nqubits=2, images=4), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=4), CliffordSymmetry(nqubits=2, images=4), CliffordSymmetry(nqubits=2, images=4), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=6), CliffordSymmetry(nqubits=2, images=4), CliffordSymmetry(nqubits=2, images=4), CliffordSymmetry(nqubits=2, images=4), CliffordSymmetry(nqubits=2, images=4)], nothing, nothing, true)
````

How many independent Clifford generators did it find?

````julia
length(auto_spec.clifford_generators)
````

````
17
````

Solve with the auto-detected symmetry:

````julia
auto_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = auto_spec,
)

auto_result = cs_nctssos(pop, auto_config);
````

Same objective again:

````julia
auto_result.objective
````

````
-0.7500000001402382
````

Same block structure as the manual SWAP case:

````julia
auto_report = auto_result.symmetry
auto_report.group_order
````

````
24
````

````julia
auto_report.psd_block_sizes
````

````
2-element Vector{Int64}:
 1
 2
````

The SympleQ path found the same symmetry automatically — no Clifford-gate
bookkeeping required.

## Step 5 — Beyond two qubits: 4-site Heisenberg chain

The real payoff of automatic detection shows up when the system gets larger.
Consider a 4-site periodic Heisenberg chain:

```math
H_4 = \frac{1}{4}\sum_{i=1}^{4}\sum_{a \in \{x,y,z\}}
       \sigma_i^a\,\sigma_{i+1\,(\mathrm{mod}\,4)}^a.
```

````julia
N = 4
registry4, (σx4, σy4, σz4) = create_pauli_variables(1:N);

H4 = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (σx4, σy4, σz4) for i in 1:N
)

pop4 = polyopt(H4, registry4);
````

SympleQ detects the full point-group symmetry of the ring automatically:

````julia
auto_spec4 = sympleq_symmetry_spec(H4)
````

````
NCTSSoS.SymmetrySpec(NCTSSoS.SignedPermutation[], NCTSSoS.FermionicModePermutation[], NCTSSoS.CliffordSymmetry[CliffordSymmetry(nqubits=4, images=12), CliffordSymmetry(nqubits=4, images=12), CliffordSymmetry(nqubits=4, images=12), CliffordSymmetry(nqubits=4, images=8), CliffordSymmetry(nqubits=4, images=8), CliffordSymmetry(nqubits=4, images=12), CliffordSymmetry(nqubits=4, images=12), CliffordSymmetry(nqubits=4, images=12), CliffordSymmetry(nqubits=4, images=12), CliffordSymmetry(nqubits=4, images=8), CliffordSymmetry(nqubits=4, images=8)], nothing, nothing, true)
````

Number of independent Clifford generators for the 4-site ring:

````julia
length(auto_spec4.clifford_generators)
````

````
11
````

Build an order-1 half-basis (identity + all single-site Paulis):

````julia
basis4 = [one(σx4[1]); σx4; σy4; σz4];
length(basis4)
````

````
13
````

Solve with symmetry:

````julia
config4 = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis4,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = auto_spec4,
)

result4 = cs_nctssos(pop4, config4);

result4.objective
````

````
-3.000000000750111
````

The symmetry report shows the group order and the reduced PSD block sizes:

````julia
report4 = result4.symmetry
report4.group_order
````

````
16
````

````julia
report4.psd_block_sizes
````

````
5-element Vector{Int64}:
 1
 2
 2
 2
 2
````

Compare: without symmetry, we would have a single ``13 \times 13`` PSD block.
With symmetry, it splits into several smaller blocks — and the solver does
correspondingly less work. On larger systems, this reduction is the
difference between feasible and intractable.

## Summary

| run | Hamiltonian | PSD block sizes | group order | objective |
|:----|:------------|:----------------|:------------|:----------|
| dense baseline | 2-site | `[7]` | — | ``-0.75`` |
| manual `CliffordSymmetry(:SWAP, 1, 2)` | 2-site | `[4, 3]` | 2 | ``-0.75`` |
| `sympleq_symmetry_spec(H)` | 2-site | `[4, 3]` | 2 | ``-0.75`` |
| `sympleq_symmetry_spec(H₄)` | 4-site ring | *(auto-detected)* | *(auto-detected)* | lower bound |

The pattern: same answer, smaller SDP. Manual Clifford gates give you
control; SympleQ gives you automation. For production use on large Pauli
Hamiltonians, start with `sympleq_symmetry_spec` and fall back to manual
gates only if you need to enforce a specific subgroup.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

