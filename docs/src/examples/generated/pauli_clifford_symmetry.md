```@meta
EditURL = "../literate/pauli_clifford_symmetry.jl"
```

# [Pauli Symmetry Reduction: Manual Clifford Gates and Automatic SympleQ Detection](@id pauli-clifford-symmetry)

The cost of a moment relaxation is driven by the size of its PSD blocks: an
order-2 half-basis on 10 qubits already produces a ``436 \times 436`` moment
matrix with thousands of distinct moment variables, and interior-point
solvers pay roughly cubic cost in the block side. Symmetry reduction is one
of the few reductions that shrinks those blocks **exactly** — the optimum is
unchanged, unlike sparsity-based relaxations which may weaken the bound.

Using symmetry takes two separate steps, and this page motivates both:

1. **Exploit.** Given a group that provably leaves the problem invariant,
   block-diagonalize the moment matrix and merge moment variables into
   orbits. NCTSSoS does this whenever you pass a [`SymmetrySpec`](@ref); the
   machinery is described in
   [Symmetry-Adapted Basis](@ref symmetry-adapted-basis).
2. **Find.** Produce that group in the first place — either write
   [`CliffordSymmetry`](@ref) generators by hand, or let
   [`sympleq_symmetry_spec`](@ref) detect them from the Hamiltonian
   automatically (the algorithm is described in
   [Clifford Symmetry Detection](@ref clifford-symmetry-detection)).

Why bother with automatic detection when a model's symmetry looks obvious?
Because **the symmetry you can see is usually not all there is** — and the
two-site Heisenberg model, the smallest example with a known answer, already
makes the point sharply:

- the symmetry everyone writes down (qubit swap) is a group of order **2**
  and splits the ``7\times 7`` moment block into sizes ``[4, 3]``;
- the model's *internal* symmetry — global rotations of the Pauli axes,
  which fix no individual term of ``H`` — is a group of order **24** that
  SympleQ finds automatically, collapsing the SDP to blocks ``[1, 2]`` with
  a *single* nontrivial moment variable;
- the union of both is a group of order **48** that reduces the whole
  ``7 \times 7`` PSD constraint to three scalars.

Detection and inspection find *different* parts of the symmetry group, so
the best spec is their union. We demonstrate this on 2 and 4 sites where
every number can be checked by hand, then time the payoff on a 10-site ring
at relaxation order 2.

**Prerequisites:** the [CHSH symmetry](@ref chsh-symmetry) example (the
same workflow with signed permutations on measurement variables) and the
[Symmetry-Adapted Basis](@ref symmetry-adapted-basis) manual page (the
architectural picture).

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
moment relaxation and then shrink the SDP with progressively more symmetry.

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
-0.7500000000000072
````

Sanity-check against the analytical value:

````julia
abs(dense_result.objective - (-0.75))
````

````
7.216449660063518e-15
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

## Step 3 — The symmetry you can see: qubit SWAP

The 2-site Heisenberg Hamiltonian is invariant under qubit exchange. In Pauli
language, the SWAP gate conjugates every single-site Pauli operator on qubit 1
into the corresponding operator on qubit 2, and vice versa:

```math
\text{SWAP}\;\sigma_i^a\;\text{SWAP}^\dagger = \sigma_j^a, \quad
\{i,j\} = \{1,2\},\; a \in \{x,y,z\}.
```

Each term ``\sigma_1^a \sigma_2^a`` is symmetric under this exchange, so
``H`` is invariant — this is the symmetry anyone spots in five seconds.

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
manual_spec = SymmetrySpec(swap);
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
1.1102230246251565e-16
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

and the 16 distinct moment variables of the dense run have merged into 9
nontrivial orbit representatives (plus the identity):

````julia
report.invariant_moment_count
````

````
9
````

Two smaller PSD constraints, fewer variables, same answer. Useful — but
this is only the beginning, because SWAP is not the symmetry that does the
heavy lifting here.

## Step 4 — The symmetry you'd probably miss: global Pauli-axis rotations

The Heisenberg coupling is **isotropic**: ``\vec{\sigma}_1 \cdot
\vec{\sigma}_2`` treats the ``x``, ``y``, ``z`` axes identically. Any global
spin rotation — the *same* single-qubit unitary applied to every site —
leaves it invariant. The Clifford subgroup of those rotations is the
24-element octahedral group: the rotations that permute the coordinate axes
``\pm x, \pm y, \pm z`` among themselves.

Unlike SWAP, these symmetries move no qubit. They act in the *internal*
Pauli-axis space, permuting the terms of ``H`` among each other. Two
generators suffice. The **global Hadamard** exchanges the ``x`` and ``z``
axes on both sites:

```math
\sigma_i^x \mapsto \sigma_i^z, \quad
\sigma_i^z \mapsto \sigma_i^x, \quad
\sigma_i^y \mapsto -\sigma_i^y \qquad (i = 1, 2),
```

which maps ``\sigma_1^x\sigma_2^x \leftrightarrow \sigma_1^z\sigma_2^z`` and
fixes ``\sigma_1^y\sigma_2^y`` (the two ``-`` signs cancel). Note what just
happened: this Clifford **permutes the terms of ``H`` without fixing them**
— no relabeling of qubits can do that. Hamiltonian-level invariance holds
even though no individual term is invariant.

Single-gate constructors like `CliffordSymmetry(:H, 1)` won't express this —
a Hadamard on qubit 1 *alone* is not a symmetry (it would map
``\sigma_1^x\sigma_2^x`` to ``\sigma_1^z\sigma_2^x``, which is not in ``H``).
For gates acting on several sites at once, pass an image dictionary
`source => (sign, target)`:

````julia
M = typeof(σx[1])
global_hadamard = CliffordSymmetry(Dict{M,Tuple{Int,M}}(
    σx[1] => (1, σz[1]), σz[1] => (1, σx[1]), σy[1] => (-1, σy[1]),
    σx[2] => (1, σz[2]), σz[2] => (1, σx[2]), σy[2] => (-1, σy[2]),
))
````

````
CliffordSymmetry(nqubits=2, images=6)
````

The **global S gate** is a quarter turn about the ``z`` axis,
``\sigma^x \mapsto \sigma^y \mapsto -\sigma^x`` on both sites; it swaps the
``\sigma^x\sigma^x`` and ``\sigma^y\sigma^y`` terms and fixes
``\sigma^z\sigma^z``:

````julia
global_s = CliffordSymmetry(Dict{M,Tuple{Int,M}}(
    σx[1] => (1, σy[1]), σy[1] => (-1, σx[1]),
    σx[2] => (1, σy[2]), σy[2] => (-1, σx[2]),
))
````

````
CliffordSymmetry(nqubits=2, images=4)
````

Together with SWAP, these two generators produce a much larger group —
NCTSSoS enumerates it and re-verifies invariance of every element:

````julia
full_manual_spec = SymmetrySpec(swap, global_hadamard, global_s);

full_manual_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = full_manual_spec,
)

full_manual_result = cs_nctssos(pop, full_manual_config);
full_manual_result.objective
````

````
-0.7500000003352729
````

The three generators close into a group of order **48** —
``2 \;(\text{SWAP}) \times 24 \;(\text{octahedral axis rotations})``:

````julia
full_report = full_manual_result.symmetry
full_report.group_order
````

````
48
````

And the SDP collapses: the ``7 \times 7`` PSD block becomes three scalars,

````julia
full_report.psd_block_sizes
````

````
3-element Vector{Int64}:
 1
 1
 1
````

with a *single* nontrivial invariant moment left — the common bond
correlator ``\langle \sigma_1^a \sigma_2^a \rangle``, the same for all three
axes by symmetry. Every other moment (single-site averages, mixed-axis
correlators) is forced to zero by some group element:

````julia
full_report.invariant_moment_count
````

````
1
````

Same exact objective, from an SDP that is barely an SDP anymore. The catch:
would you have written `global_hadamard` and `global_s` down unprompted —
and would you trust yourself to have found *all* such generators on a model
that isn't a two-site textbook exercise? This is the problem automatic
detection solves.

## Step 5 — Automatic SympleQ detection — and its blind spot

[`sympleq_symmetry_spec`](@ref) inspects the Hamiltonian itself. It encodes
the Pauli terms as a binary symplectic tableau, builds a coloured graph
whose automorphisms preserve coefficients, commutation relations, and
``\mathbb{F}_2``-linear structure, lifts each automorphism to a Clifford
(including a GF(2) phase solve so the signs work out), and returns only
generators that verifiably leave ``H`` invariant. Details:
[Clifford Symmetry Detection](@ref clifford-symmetry-detection).

````julia
auto_spec = sympleq_symmetry_spec(H);
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
auto_result.objective
````

````
-0.7500000001402383
````

The detected generators close into a group of order **24**:

````julia
auto_report = auto_result.symmetry
auto_report.group_order
````

````
24
````

That is precisely the global octahedral rotation group from Step 4 — the
part of the symmetry a human is likely to *miss*, found automatically with
no hints. It is already enough to collapse the moment matrix almost
completely:

````julia
auto_report.psd_block_sizes
````

````
2-element Vector{Int64}:
 1
 2
````

````julia
auto_report.invariant_moment_count
````

````
1
````

But notice what is **not** there: the group order is 24, not 48. SympleQ did
not return SWAP. The reason is structural, not a bug. The detection graph
sees a Clifford only through the *permutation it induces on the terms of*
``H`` — and SWAP fixes all three terms ``\sigma_1^a\sigma_2^a``
individually, so at the term level it is indistinguishable from the
identity. Cliffords that move operators while fixing every term can
therefore escape detection (see
[Scaling and limitations](@ref sympleq-limits)).

The fix is one line: **union the generators you know with the generators
SympleQ finds.** A [`SymmetrySpec`](@ref) accepts any mix:

````julia
combined_spec = SymmetrySpec(auto_spec.clifford_generators..., swap);

combined_config = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = combined_spec,
)

combined_result = cs_nctssos(pop, combined_config);
combined_result.objective
````

````
-0.7500000003352729
````

The union recovers the full order-48 group and the three-scalar SDP of
Step 4 — without anyone having to invent the axis rotations by hand:

````julia
combined_report = combined_result.symmetry
(combined_report.group_order, combined_report.psd_block_sizes)
````

````
(48, [1, 1, 1])
````

## Step 6 — Bigger ring, same lesson: the 4-site Heisenberg chain

On a 4-site periodic chain

```math
H_4 = \frac{1}{4}\sum_{i=1}^{4}\sum_{a \in \{x,y,z\}}
       \sigma_i^a\,\sigma_{i+1\,(\mathrm{mod}\,4)}^a
```

the "visible" symmetry is the spatial point group of the ring: translations
and a reflection (the dihedral group ``D_4``, order 8). The "invisible" one
is again the global axis-rotation tower. Let's see who finds what.

````julia
N = 4
registry4, (σx4, σy4, σz4) = create_pauli_variables(1:N);

H4 = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (σx4, σy4, σz4) for i in 1:N
)

pop4 = polyopt(H4, registry4);
````

Order-1 half-basis: identity plus all single-site Paulis (13 monomials).

````julia
basis4 = [one(σx4[1]); σx4; σy4; σz4];
length(basis4)
````

````
13
````

Write the spatial generators by hand, as image dictionaries: a lattice
translation ``i \mapsto i+1 \pmod 4`` and a reflection ``i \mapsto 5-i``,
both acting identically on the three Pauli types:

````julia
M4 = typeof(σx4[1])
translation = CliffordSymmetry(Dict{M4,Tuple{Int,M4}}(
    op[i] => (1, op[mod1(i + 1, N)]) for op in (σx4, σy4, σz4), i in 1:N
))
reflection = CliffordSymmetry(Dict{M4,Tuple{Int,M4}}(
    op[i] => (1, op[N + 1 - i]) for op in (σx4, σy4, σz4), i in 1:N
))
````

Solve with the spatial point group only:

````julia
spatial_config4 = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis4,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = SymmetrySpec(translation, reflection),
)

spatial_result4 = cs_nctssos(pop4, spatial_config4);
spatial_result4.objective
````

````
-2.999999999999893
````

The full ``D_4`` (order 8) splits the ``13 \times 13`` block into three:

````julia
spatial_report4 = spatial_result4.symmetry
(spatial_report4.group_order, spatial_report4.psd_block_sizes,
 spatial_report4.invariant_moment_count)
````

````
(8, [4, 3, 3], 15)
````

Now automatic detection — the Hamiltonian has 12 terms, and SympleQ returns
its verified generators in a fraction of a second:

````julia
auto_spec4 = sympleq_symmetry_spec(H4);
length(auto_spec4.clifford_generators)

auto_config4 = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis4,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = auto_spec4,
)

auto_result4 = cs_nctssos(pop4, auto_config4);
auto_result4.objective
````

````
-3.00000000075011
````

A different group of order 16 — again containing axis-mixing Cliffords the
spatial spec lacks, again missing part of the point group the spatial spec
has:

````julia
auto_report4 = auto_result4.symmetry
(auto_report4.group_order, auto_report4.psd_block_sizes,
 auto_report4.invariant_moment_count)
````

````
(16, [1, 2, 2, 2, 2], 8)
````

Neither alone is the whole story; their union is strictly bigger than both:

````julia
combined_config4 = SolverConfig(
    optimizer    = SILENT_MOSEK,
    moment_basis = basis4,
    cs_algo      = NoElimination(),
    ts_algo      = NoElimination(),
    symmetry     = SymmetrySpec(
        auto_spec4.clifford_generators..., translation, reflection,
    ),
)

combined_result4 = cs_nctssos(pop4, combined_config4);
combined_result4.objective
````

````
-2.999999999999906
````

Order 64 — and the ``13 \times 13`` moment matrix dissolves into **seven
scalars** with just four invariant moments:

````julia
combined_report4 = combined_result4.symmetry
(combined_report4.group_order, combined_report4.psd_block_sizes,
 combined_report4.invariant_moment_count)
````

````
(64, [1, 1, 1, 1, 1, 1, 1], 4)
````

Recap of the 4-site runs (all certify the same bound):

| spec | group order | PSD blocks | invariant moments |
|:-----|:------------|:-----------|:------------------|
| dense | — | `[13]` | 67 |
| manual ``D_4`` (translation + reflection) | 8 | `[4, 3, 3]` | 15 |
| `sympleq_symmetry_spec` | 16 | `[1, 2, 2, 2, 2]` | 8 |
| union of both | 64 | `[1, 1, 1, 1, 1, 1, 1]` | 4 |

## Step 7 — The payoff at scale: 10-site Heisenberg ring at order 2

Steps 1–6 used **order-1** relaxations, where the moment matrices are small
enough that symmetry overhead can exceed solver time. The real payoff of
symmetry reduction shows at **higher relaxation orders**, where PSD block
sizes grow quadratically and the solver’s ``O(k^3)`` cost per block
dominates.

Consider a 10-site periodic Heisenberg chain, now at **order 2**:

```math
H_{10} = \frac{1}{4}\sum_{i=1}^{10}\sum_{a \in \{x,y,z\}}
       \sigma_i^a\,\sigma_{i+1\,(\mathrm{mod}\,10)}^a.
```

At order 2, the half-basis includes all products of up to two single-site
Pauli operators — 436 monomials instead of 31. To see *where* the time goes,
we will also report the raw interior-point time inside each total via
`JuMP.solve_time`.

````julia
import JuMP

N10 = 10
registry10, (σx10, σy10, σz10) = create_pauli_variables(1:N10);

H10 = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N10)]
    for op in (σx10, σy10, σz10) for i in 1:N10
)

pop10 = polyopt(H10, registry10);
````

Use `order=2` to let NCTSSoS auto-generate the full degree-2 half-basis.

````julia
dense_config10 = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
);
````

### Warm up before timing

The stopwatch must measure problem solving, not Julia compiling methods for
the first time. `make examples` already precompiles the docs environment, but
first use of this exact order-2 dense and symmetry-reduced path can still pay
JIT and one-time solver-library setup costs. We therefore run each code path
once and discard those results before recording any timings.

The measured numbers below still include SDP assembly, SympleQ detection,
model construction, and the optimizer solve. What they exclude is cold-start
compilation noise.

````julia
warm_spec10 = sympleq_symmetry_spec(H10);
warm_config10 = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = warm_spec10,
);
cs_nctssos(pop10, dense_config10);
cs_nctssos(pop10, warm_config10);
GC.gc()
````

### Dense baseline (timed)

A single ``436 \times 436`` PSD block with over twenty thousand distinct
moment variables — Mosek needs a couple of minutes:

````julia
GC.gc()
t_dense10 = @elapsed dense_result10 = cs_nctssos(pop10, dense_config10);

dense_result10.objective
````

````
-4.51815932435283
````

````julia
dense_result10.moment_matrix_sizes
````

````
1-element Vector{Vector{Int64}}:
 [436]
````

````julia
t_dense10
````

````
88.598676167
````

Nearly all of that wall-clock time is the interior-point solver itself:

````julia
t_dense_mosek = JuMP.solve_time(dense_result10.model)
````

````
88.2443618774414
````

### SympleQ detection + symmetry-reduced solve (timed)

SympleQ detection depends on the *Hamiltonian* (30 Pauli terms), not on the
relaxation order, so it stays cheap relative to the solve.

````julia
GC.gc()
t_detect10 = @elapsed auto_spec10 = sympleq_symmetry_spec(H10);

length(auto_spec10.clifford_generators)
````

````
11
````

````julia
t_detect10
````

````
0.001293625
````

Solve with the auto-detected symmetry at order 2:

````julia
config10 = SolverConfig(
    optimizer = SILENT_MOSEK,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = auto_spec10,
)

GC.gc()
t_sym10 = @elapsed result10 = cs_nctssos(pop10, config10);

result10.objective
````

````
-4.518159330977996
````

````julia
report10 = result10.symmetry
report10.group_order
````

````
16
````

````julia
report10.psd_block_sizes
````

````
10-element Vector{Int64}:
 51
 40
 25
 20
 25
 20
 25
 30
 50
 50
````

````julia
t_sym10
````

````
0.828725625
````

And the solver time hiding inside that total:

````julia
t_sym_mosek = JuMP.solve_time(result10.model)
````

````
0.592249870300293
````

### Comparison

The ``436 \times 436`` dense PSD block splits into **10 smaller blocks**
(largest ~ 51), and the distinct moment variables drop from ~20,000 to
~1,400. The objectives match — the reduction is exact.

````julia
println("Steady-state timing; warm-up solves discarded")
println("run       total        of which Mosek   moment vars   PSD blocks")
println("dense     $(round(t_dense10; digits=1))s      $(round(t_dense_mosek; digits=1))s           $(dense_result10.n_unique_moment_matrix_elements)         $(dense_result10.moment_matrix_sizes)")
println("symmetry  $(round(t_detect10 + t_sym10; digits=1))s       $(round(t_sym_mosek; digits=1))s             $(result10.n_unique_moment_matrix_elements)          $(report10.psd_block_sizes)")
println("  symmetry time = detect $(round(t_detect10; sigdigits=2))s + solve $(round(t_sym10; digits=1))s")
println("End-to-end speedup: $(round(t_dense10 / (t_detect10 + t_sym10); digits=1))x")
println("Solver-only speedup: $(round(t_dense_mosek / t_sym_mosek; digits=0))x")
println("Moment variable reduction: $(dense_result10.n_unique_moment_matrix_elements) → $(result10.n_unique_moment_matrix_elements)")
println("Objective difference: $(round(abs(dense_result10.objective - result10.objective); sigdigits=3))")
println("Group order: $(report10.group_order), generators: $(length(auto_spec10.clifford_generators))")
````

````
Steady-state timing; warm-up solves discarded
run       total        of which Mosek   moment vars   PSD blocks
dense     88.6s      88.2s           20686         [[436]]
symmetry  0.8s       0.6s             1421          [51, 40, 25, 20, 25, 20, 25, 30, 50, 50]
  symmetry time = detect 0.0013s + solve 0.8s
End-to-end speedup: 106.7x
Solver-only speedup: 149.0x
Moment variable reduction: 20686 → 1421
Objective difference: 6.63e-9
Group order: 16, generators: 11

````

### Where does the remaining time go?

Look at the two speedup numbers separately. The **solver-only** speedup is
roughly two orders of magnitude: Mosek's share collapses from ~130 s to
under a second, because ten blocks of side ≤ 51 are vastly cheaper than one
``436 \times 436`` block. That is the genuine, exact payoff of symmetry
reduction — and it grows with system size and relaxation order, since the
dense solver cost explodes while the reduced blocks stay manageable.

The **end-to-end** speedup is smaller because the symmetry run pays two
overheads the dense run does not: SympleQ detection (well under a second
here, and independent of relaxation order) and constructing the
symmetry-adapted relaxation — decomposing the 436-dimensional representation into isotypic
blocks and congruence-transforming the moment matrix into them. That
construction cost grows with the order of the group and the number of
blocks, not with what the solver does afterwards.

Part of that construction is a **soundness certificate**: replacing one
big PSD block by small diagonal blocks is only exact if the off-diagonal
blocks of the transformed moment matrix vanish, and NCTSSoS verifies this
at runtime rather than trusting the supplied generators blindly. The
`offblock_check` option of [`SymmetrySpec`](@ref) selects the verification
mode: the default `:randomized` certificate evaluates every off-diagonal
block at deterministic pseudo-random values of the moment variables (each
entry is linear in them, so vanishing at generic points certifies the
symbolic property); `:full` expands every entry symbolically — an
order-of-magnitude more construction time, useful for diagnosing a failing
certificate entry-by-entry; `:off` skips the check.

What about the “union the specs” advice of Step 6? On this ring you can
union `auto_spec10` with the 10-site translation and reflection (the
dihedral group ``D_{10}``, exactly as in Step 6). The combined group has
order **160** and shreds the SDP to 37 blocks of side ≤ 11 with only 179
moment variables — Mosek then solves it in ~0.1 s. But on the same
hardware as the run above, building that finer decomposition takes
roughly twice as long as building the order-16 one, so the total is
*slower* than `auto_spec10` alone despite the smaller SDP. A bigger group
always means a smaller SDP, not always a faster total.

Rules of thumb:

- When the solver dominates (large basis, high order, or many repeated
  solves of the same relaxation), symmetry reduction wins decisively.
- When a one-shot solve is cheap, or the group is enormous relative to the
  moment matrix, detection and construction can eat the gains — measure
  both parts (`JuMP.solve_time` vs. total) before scaling up.

### Scaling beyond 10 sites

The gap between dense and symmetry-reduced widens rapidly with system
size. The numbers below were measured once with the same protocol as the
live run above (same Heisenberg ring at order 2, same machine — 28 cores,
111 GB RAM — Mosek, warm-up runs discarded); they are quoted statically
because the dense baselines take from twenty minutes to two hours:

| sites | dense total (Mosek) | symmetry total (Mosek) | end-to-end | solver-only | objective difference |
|:------|:--------------------|:-----------------------|:-----------|:------------|:---------------------|
| 12 | 1098 s (1093 s) | 31 s (3.1 s) | 35× | 354× | ``8.7 \times 10^{-9}`` |
| 14 | 6969 s (6957 s) | 102 s (11.8 s) | 68× | 588× | ``4.1 \times 10^{-8}`` |
| 16 | intractable | 330 s (31.9 s) | — | — | — |

The structure behind those totals: at 12 sites the dense relaxation is a
single ``631 \times 631`` block with 46,666 moment variables, against 10
symmetry-adapted blocks of side ≤ 73 with 3,178; at 14 sites,
``862 \times 862`` with 91,771 against blocks ≤ 99 with 6,224; at 16
sites the symmetry-reduced SDP has 10 blocks of side ≤ 129 with 11,077
moment variables and solves in about five minutes, while the dense
``1129 \times 1129`` baseline was not attempted — the 14-site dense solve
already consumed ~99 of the machine's 111 GB, and both its time and
memory grow steeply from there. The detected group is the same order-16
group at every size, found in ≈ 0.5 s.

With no dense run to compare against at 16 sites, is the bound still
trustworthy? Physics provides the check. Per site, the order-2 bounds are

```math
E_{10}/10 = -0.451816,\quad E_{12}/12 = -0.449407,\quad
E_{14}/14 = -0.447994,\quad E_{16}/16 = -0.447093,
```

a monotone family approaching the exact thermodynamic-limit ground-state
energy per site of the Heisenberg chain, ``1/4 - \ln 2 \approx
-0.443147``, from below — just as the true finite-ring energies do, with
every bound below its ring's true energy, as a lower bound must be. The symmetry-reduced numbers agree with the dense
objectives to ``\sim 10^{-8}`` wherever both exist, and the certificate
described above guards the block structure at every size.

At these sizes the remaining cost is almost entirely the symmetry-adapted
*construction* (about 90% of the symmetry total at 16 sites), not the
solve — so the gap between the two columns is also a roadmap for future
optimization.

## Summary

The 2-site numbers, side by side:

| spec | group order | PSD blocks | invariant moments | objective |
|:-----|:------------|:-----------|:------------------|:----------|
| dense baseline | — | `[7]` | 16 | ``-0.75`` |
| manual SWAP | 2 | `[4, 3]` | 9 | ``-0.75`` |
| manual SWAP + axis rotations | 48 | `[1, 1, 1]` | 1 | ``-0.75`` |
| `sympleq_symmetry_spec` | 24 | `[1, 2]` | 1 | ``-0.75`` |
| auto ∪ manual SWAP | 48 | `[1, 1, 1]` | 1 | ``-0.75`` |

Three takeaways:

1. **Symmetry reduction is exact** — every row certifies the same bound,
   only the SDP shrinks. At low orders the overhead is not worth it; by
   order 2 on 10 qubits the pure solver time drops by two orders of
   magnitude (~130 s → under a second), and the advantage widens with
   size — 68× end-to-end at 14 sites, and at 16 sites the symmetry-reduced
   relaxation solves in minutes where the dense one is out of reach. The
   remaining overhead — detection plus symmetry-adapted construction —
   grows with the group order, so the largest group is not automatically
   the fastest run.
2. **Inspection and detection see different things.** Humans spot spatial
   symmetries (SWAP, translations); SympleQ finds internal, axis-mixing
   Clifford symmetries that fix no individual term — but can miss Cliffords
   that fix *every* term, because those induce the identity permutation on
   the detection graph.
3. **So use both.** For production runs, start from
   `sympleq_symmetry_spec(H)` and splice in the structural generators you
   know: `SymmetrySpec(auto.clifford_generators..., your_gates...)`. The
   union is often strictly stronger than either part — as a *reduction*;
   check Step 7's caveat before assuming it is also faster end-to-end.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

