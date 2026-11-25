```@meta
EditURL = "../literate/pauli_clifford_symmetry.jl"
```

# [Pauli Symmetry Reduction: Manual Clifford Gates and Automatic SympleQ Detection](@id pauli-clifford-symmetry)

An order-2 relaxation on 10 qubits has 436 rows. Symmetry can replace its
large PSD block with smaller blocks without weakening the invariant
relaxation. There are two separate tasks:

| Task | Input | Call | Output |
|:-----|:------|:-----|:-------|
| Find symmetry | Hamiltonian | `detect_symmetry(H)` | Verified generators |
| Use symmetry | Problem and rows | `NCTSSoS.build(problem, relaxation)` | Reduced conic program |
| Solve | Built program | `solve(built; solver)` | Bound and moments |

A symmetry that is easy to see need not give the smallest blocks. On two
sites, qubit exchange gives sides `[4, 3]`, while global Pauli-axis rotations
can reduce the problem further. We first compare those choices, then study
four- and ten-site rings.

Prerequisites: [CHSH symmetry](@ref chsh-symmetry),
[Symmetry-Adapted Basis](@ref symmetry-adapted-basis), and
[Clifford Symmetry Detection](@ref clifford-symmetry-detection).

````julia
using NCTSSoS, Mosek
import NCTSSoS

solver = MosekBackend();
````

## Step 1 — The two-site Heisenberg model

The isotropic Hamiltonian is

```math
H = \frac14(X_1X_2 + Y_1Y_2 + Z_1Z_2).
```

Its ground-state energy is `-0.75`. The Pauli algebra supplies the
commutation and anticommutation rules; no operator constraints are needed.

````julia
registry, (X, Y, Z) = pauli_variables(2);
H = sum(op[1] * op[2] / 4 for op in (X, Y, Z))
problem = Problem(H);
````

Use the same seven rows in every two-site run: identity and the six
single-site operators.

````julia
rows = [one(X[1]); X; Y; Z];
relaxation = Relaxation(rows=rows);
length(rows)
````

````
7
````

## Step 2 — Dense baseline

````julia
dense = NCTSSoS.build(problem, relaxation);
dense_result = solve(dense; solver);
dense_result.bound
````

````
-0.7500000000000073
````

`NCTSSoS.block_sides` reports the final PSD/HPSD sides before realification.
`sizes.moments` includes the identity; `real_variables` counts free real
coordinates and does not count its fixed value.

````julia
(NCTSSoS.block_sides(dense), dense_result.sizes.moments, dense_result.sizes.real_variables)
````

````
([7], 16, 15)
````

## Step 3 — The visible symmetry: qubit SWAP

Qubit exchange sends every Pauli on site 1 to the same Pauli on site 2.
Each term of `H` is unchanged. `SitePermutation` expresses this without
separate maps for the three axes.

````julia
swap = SitePermutation([2, 1]);
swapped = NCTSSoS.build(Problem(H; symmetry=(swap,)), relaxation);
swap_result = solve(swapped; solver);
(swap_result.bound, NCTSSoS.block_sides(swapped), swap_result.sizes.real_variables)
````

````
(-0.7500000000000072, [4, 3], 9)
````

The order-two exchange group separates symmetric and antisymmetric rows.
It keeps the bound but leaves nine free coordinates. We can do more.

## Step 4 — Global Pauli-axis rotations

The coupling treats all three axes equally. The same rotation on both
sites preserves `H`, even when it permutes the three terms. A Hadamard
on only one site does not preserve `H`.

A global Hadamard sends `X → Z`, `Z → X`, and `Y → -Y`. The two minus
signs in `Y₁Y₂` cancel. A global phase gate sends `X → Y`, `Y → -X`, and
`Z → Z`. Together they generate the 24-element octahedral rotation group.

`CliffordMap` takes images of `X₁,…,Xₙ,Z₁,…,Zₙ`, as Pauli words in
that order. The images of Y and their signs follow from Pauli multiplication.
`NCTSSoS.terms` gives coefficient-word pairs; each operator below has one term.

````julia
global_hadamard = CliffordMap([last(only(NCTSSoS.terms(p))) for p in [Z; X]]);
global_s = CliffordMap([last(only(NCTSSoS.terms(p))) for p in [Y; Z]]);
manual_actions = (swap, global_hadamard, global_s);
manual = NCTSSoS.build(Problem(H; symmetry=manual_actions), relaxation);
manual_result = solve(manual; solver);
(manual_result.bound, NCTSSoS.block_sides(manual), manual_result.sizes.real_variables)
````

````
(-0.7500000003352736, [1, 1, 1], 1)
````

SWAP and axis rotations generate a group of order 48. Its reduced SDP has
three scalar blocks and one free coordinate: the common bond correlator.
Single-site and mixed-axis moments vanish by symmetry.

````julia
expectation(manual_result, X[1] * X[2])
````

````
-0.9999999990976801
````

## Step 5 — Automatic detection and its blind spot

`detect_symmetry` encodes the Hamiltonian terms in a symplectic graph. It
finds graph automorphisms, lifts them to Clifford maps, and checks their
phases and invariance. See [Clifford Symmetry Detection](@ref clifford-symmetry-detection).

````julia
detected = detect_symmetry(H);
(length(detected.generators), length(detected.moment_only))
````

````
(3, 0)
````

Build and solve with the detected symmetry:

````julia
automatic = NCTSSoS.build(Problem(H; symmetry=detected), relaxation);
auto_result = solve(automatic; solver);
(auto_result.bound, NCTSSoS.block_sides(automatic), auto_result.sizes.real_variables)
````

````
(-0.7500000001402383, [1, 2], 1)
````

Detection sees maps through their action on the terms of `H`. A map such
as SWAP that fixes every term can be invisible at that level. Add known
spatial generators to the detected block generators, and retain any
moment-only or continuous generators.

````julia
combined_symmetry = NCTSSoS.Symmetry(
    detected.generators..., swap;
    moment_only=detected.moment_only, continuous=detected.continuous,
);
combined = NCTSSoS.build(Problem(H; symmetry=combined_symmetry), relaxation);
combined_result = solve(combined; solver);
(combined_result.bound, NCTSSoS.block_sides(combined), combined_result.sizes.real_variables)
````

````
(-0.750000000335274, [1, 1, 1], 1)
````

Group orders above describe the specified mathematical groups. The public
size report does not compute their orders; a generator count is not a group
order. The order-48 subgroup is not the full Clifford stabilizer of `H`:
for example, `exp(-iπ X₁X₂/4)` preserves `H` but changes some word degrees.

## Step 6 — Four-site ring: Hamiltonian symmetry and row closure

For a periodic ring,

```math
H_N = \frac14\sum_{i=1}^N\sum_{a=x,y,z}\sigma_i^a\sigma_{i+1}^a.
```

Translation and reflection generate the spatial dihedral group, of order
`2N`. Global axis rotations give a degree-preserving subgroup of order
`48N`. We construct it explicitly rather than discard individual detected
generators: products of degree-changing maps can preserve degree.

````julia
function ring_actions(ops)
    x, y, z = ops
    n = length(x)
    translation = SitePermutation([mod1(i + 1, n) for i in 1:n])
    reflection = SitePermutation(collect(n:-1:1))
    hadamard = CliffordMap([last(only(NCTSSoS.terms(p))) for p in [z; x]])
    phase = CliffordMap([last(only(NCTSSoS.terms(p))) for p in [y; z]])
    return (translation, reflection, hadamard, phase)
end

registry4, ops4 = pauli_variables(4);
H4 = sum(op[i] * op[mod1(i + 1, 4)] / 4 for op in ops4 for i in 1:4);
actions4 = ring_actions(ops4);
relaxation4 = Relaxation(order=1);
spatial4 = NCTSSoS.build(Problem(H4; symmetry=actions4[1:2]), relaxation4);
spatial_result4 = solve(spatial4; solver);
(spatial_result4.bound, NCTSSoS.block_sides(spatial4), spatial_result4.sizes.real_variables)
````

````
(-2.9999999999998943, [4, 3, 3], 15)
````

The spatial group splits the 13-row moment matrix. Add axis rotations:

````julia
full4 = NCTSSoS.build(Problem(H4; symmetry=actions4), relaxation4);
full_result4 = solve(full4; solver);
(full_result4.bound, NCTSSoS.block_sides(full4), full_result4.sizes.real_variables)
````

````
(-2.999999999998245, [1, 1, 1, 1], 2)
````

Detection separates degree-preserving block generators from moment-only
generators. A degree-changing symmetry is not a licence to rotate rows
outside the selected family. The latter generators instead impose relations
on represented moments; see the manual's supported scope.

````julia
detected4 = detect_symmetry(H4);
(length(detected4.generators), length(detected4.moment_only))
````

````
(6, 1)
````

This family can be passed directly to the problem:

````julia
automatic4 = NCTSSoS.build(Problem(H4; symmetry=detected4), relaxation4);
auto_result4 = solve(automatic4; solver);
(auto_result4.bound, NCTSSoS.block_sides(automatic4), auto_result4.sizes.real_variables)
````

````
(-2.999999999998246, [1, 1, 1, 1], 2)
````

Adding spatial generators already in the explicit subgroup changes nothing:

````julia
redundant4 = NCTSSoS.build(Problem(H4; symmetry=(actions4..., actions4[1:2]...)), relaxation4);
redundant_result4 = solve(redundant4; solver);
(redundant_result4.bound, NCTSSoS.block_sides(redundant4), redundant_result4.sizes.real_variables)
````

````
(-2.999999999998245, [1, 1, 1, 1], 2)
````

## Step 7 — Ten-site ring at order 2

The dense relaxation has 436 rows. Small blocks can reduce optimization
time, but detection and block construction cost time too. Measure the
stages separately before choosing a symmetry group.

The first two timers measure wall time for build and backend preparation.
The third value is Mosek's optimizer time for that call, including presolve;
it excludes package build, preparation and result extraction. A single run
is an illustration, not a speedup benchmark. Package wall times can include
first-use compilation, and all measurements depend on the machine.

````julia
registry10, ops10 = pauli_variables(10);
H10 = sum(op[i] * op[mod1(i + 1, 10)] / 4 for op in ops10 for i in 1:10);
relaxation10 = Relaxation(order=2);

function staged_run(problem, relaxation; solver)
    build_seconds = @elapsed built = NCTSSoS.build(problem, relaxation)
    prepare_seconds = @elapsed NCTSSoS.prepare(built; solver)
    result = solve(built; solver)
    solve_seconds = result.solution.diagnostics.solve_seconds
    return (; built, result, build_seconds, prepare_seconds, solve_seconds)
end

dense10 = staged_run(Problem(H10), relaxation10; solver);
(dense10.result.bound, NCTSSoS.block_sides(dense10.built), dense10.result.sizes.real_variables)
````

````
(-4.518159324352234, [436], 20685)
````

Seconds, in order: build wall time, prepare wall time, Mosek optimizer time.

````julia
(dense10.build_seconds, dense10.prepare_seconds, dense10.solve_seconds)
````

````
(0.007168266, 0.119272876, 313.87293162999777)
````

Detection depends on the 30 Hamiltonian terms, not on the relaxation order.
Inspect its two generator classes, then use the explicit order-480
lattice-and-axis subgroup for the timed reduced solve.

````julia
detect_seconds = @elapsed detected10 = detect_symmetry(H10);
(detect_seconds, length(detected10.generators), length(detected10.moment_only))
````

````
(0.280939249, 6, 1)
````

Inspect the detected family's final block sides without solving it:

````julia
detected_build10 = NCTSSoS.build(Problem(H10; symmetry=detected10), relaxation10);
NCTSSoS.block_sides(detected_build10)
````

````
29-element Vector{Int64}:
 6
 5
 1
 4
 5
 4
 4
 6
 4
 5
 5
 5
 5
 4
 4
 6
 4
 5
 5
 5
 5
 2
 2
 2
 2
 4
 2
 2
 2
````

The timed run below uses the explicit subgroup, not this detected family.

````julia
actions10 = ring_actions(ops10);
reduced10 = staged_run(Problem(H10; symmetry=actions10), relaxation10; solver);
(reduced10.result.bound, NCTSSoS.block_sides(reduced10.built), reduced10.result.sizes.real_variables)
````

````
(-4.518159328030607, [6, 5, 5, 1, 4, 4, 4, 4, 6, 5, 5, 5, 5, 4, 4, 4, 6, 5, 5, 5, 5, 2, 2, 2, 2, 2, 4, 2, 2], 71)
````

Seconds, in the same order:

````julia
(reduced10.build_seconds, reduced10.prepare_seconds, reduced10.solve_seconds)
````

````
(0.01355881, 0.015115917, 0.04896847900090506)
````

Compare bounds with tolerance `1e-6`. This checks numerical agreement,
not an exact arithmetic certificate of the computed bound.

````julia
abs(dense10.result.bound - reduced10.result.bound)
````

````
3.6783731616196746e-9
````

`NCTSSoS.build` checks invariance of the objective and constraints exactly and
has no option to skip it. The reduction uses group projectors and needs no
numerical off-block test; see [Symmetry-Adapted Basis](@ref symmetry-adapted-basis).

### Historical scaling beyond ten sites

These are historical measurements from the former implementation, not
measurements of this API. They used an incomplete order-16 detected
subgroup, Mosek, 28 cores and 111 GB RAM, with warm-up runs discarded.
They illustrate why both construction and solver time matter; they do not
predict the cost of the order-480 subgroup above.

| Sites | Dense total (solver) | Reduced total (solver) | Bound difference |
|:------|:---------------------|:-----------------------|:-----------------|
| 12 | 1098 s (1093 s) | 31 s (3.1 s) | `8.7e-9` |
| 14 | 6969 s (6957 s) | 102 s (11.8 s) | `4.1e-8` |
| 16 | Not attempted | 330 s (31.9 s) | — |

The old dense sides were 631 and 862 at 12 and 14 sites; at 16 sites the
side would be 1129. The 14-site run used about 99 GB. The historical
per-site bounds were `-0.451816`, `-0.449407`, `-0.447994`, and `-0.447093`
at 10, 12, 14, and 16 sites. They approach the infinite-chain value
`1/4 - log(2) ≈ -0.443147` from below. No new large-ring run is made here.

## Summary

1. Symmetry can keep the bound while shrinking the PSD blocks. Read current
   sizes from `NCTSSoS.block_sides` and `result.sizes`, not historical timings.
2. Spatial symmetries and detected axis symmetries complement each other.
   Detection can miss maps that fix every Hamiltonian term.
3. Row closure and Hamiltonian invariance are different requirements. Use
   block generators on rows and moment-only generators for moment relations.
4. Compare build and prepare wall times with Mosek optimizer time separately. A larger group
   need not give a faster end-to-end run.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

