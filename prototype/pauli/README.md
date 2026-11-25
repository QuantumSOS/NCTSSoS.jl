# Pauli relaxation prototype

`PauliWorkflowRedesign` is an experimental front end to NCTSSoS. It does not change the default solver or export a new API from `NCTSSoS`.

## Compile the relaxation before choosing the solver form

1. **Select the Gram rows.** The sparsity policy specifies which Pauli words can appear in a test operator. The graph-distance route builds them directly from the interaction graph. A sparse row pool changes the relaxation; it is not an exact reduction of the dense hierarchy.
2. **Verify symmetry on what will be solved.** Check candidate actions against the Hamiltonian, observable and selected rows before using them to identify moments or split cones. The bounded detector need not find every symmetry. In a benchmark that requires a periodic layout, a missing translation is a failed audit, not a successful smaller computation.
3. **Add strengthening in the same moment coordinates.** Record the actual stationarity equations, positive-semidefinite blocks and local density constraints. Requested flags are not a specification of the resulting model: adding Gram rows can also make previously unsupported stationarity equations available.
4. **Choose how to submit the plan.** `compile_relaxation` constructs the plan; `solve_ground_state` submits it. Explicit moment and sum-of-squares forms can have different numerical conditioning even for the same relaxation. Keep the model fingerprint, native status, primal/dual objectives and residuals in the original units. An accepted numerical estimate is not a rational certificate.

For a one-factor row-pool comparison, match the objective, symmetry actions, density blocks and stationarity coefficients by physical Pauli word, not just by count. On the 10×10 square, changing the pair stencil from a diamond to a box doubled the site pairs (2,000 to 4,000). Unrestricted assembly also added 115 stationarity equations. The isolated experiment retained the original 49 equations and five density blocks; it was **audited, not solved**.

## Options that do not require a square lattice

| Need | Interface | Limit |
| --- | --- | --- |
| Reuse known Pauli symmetries | `NCTSSoS.detect_symmetry(H; ga_backend=:seeded, known_generators=...)` | Unconstrained Pauli objective; checks every supplied action. The default is still backtracking. |
| Split a Gram matrix by diagonal Pauli signs | Supply the verified sign actions in `NCTSSoS.Symmetry(...)`. | The fast path applies to sign-only unitary stabilizers. Other actions use the existing projector route. |
| Choose rows from the interaction graph | `PauliSparsitySpec(support_mode=:graph_motifs, power_graph_radius=r)` | Pair distance uses the original graph; orders 3 and 4 also include induced two-edge paths and chordless four-cycles. This changes the relaxation. |
| Avoid a full triplet copy on transfer | `NCTSSoS.prepare_standard(form, NCTSSoS.MosekBackend(form=:moment))` | Writes CSC coefficients in batches. The cone constraints and objective do not change. |
| Shorten a supplied SOS proof | `prune_sos_factors(H, level; squares, loss_budget)` | One absolute loss budget for all removed terms. Replay with `certify_rational_sos(H, level; squares=...)` and the same sector multipliers, if any. No solver needed. |

`solve_ground_state(plan, H; optimizer, prune_budget=...)` can extract factors for
an intact, unprojected moment cone. When the lift is unavailable, `certified_bound`
is `nothing`; `bound` remains a numerical solver result.

## Compile the periodic square XXX witness relaxation

For the historical 16×16 form, start Julia from the repository root with
`OPENBLAS_NUM_THREADS=16 julia --project --threads=16,0` (Julia 1.12.7).
The compiler does not load an optimizer or need a license, saved point, or coordinate map.
Different BLAS reduction orders can make a near-real Fourier sector fail its
reality check. The compiler then keeps the full Hermitian cone, which uses more rows.

```julia
import NCTSSoS
include("prototype/pauli/PauliWorkflowRedesign.jl")
using .PauliWorkflowRedesign

L = 16
_, axes = NCTSSoS.pauli_variables(L^2)
east(i) = i % L == 0 ? i-L+1 : i+1
north(i) = mod1(i+L,L^2)
H = sum(a[i]*(a[east(i)]+a[north(i)]) for i in 1:L^2 for a in axes)/(4L^2)
plan = compile_square_xxx_witness(H, L)
form = plan.standard
```

This is the Tencent `cut-recovery` recipe: corner and straight triplet motifs,
pair shell 2, the 2×3 RDM, both four-site lines and the five-site cross, and a
single scalar PSD row from the [exact 3×3 witness](square_witness3x3.tsv). The
TSV gives products of disjoint Pauli dots (XX + YY + ZZ), in row-major tile
order, and the exact local operator floor −7607443585/68719476736. Positive
definiteness above that floor was checked with integer arithmetic in all ten
magnetization sectors. To replay that proof without Julia or MOSEK, run
`uv run --with python-flint python prototype/pauli/benchmark/verify_square_witness.py prototype/pauli/square_witness3x3.tsv`.
Projection into the compiled moments adds one row and
no variables. Use `plan.standard` for submission; `plan.program` is the
pre-cut program and has no Gram lifts for the witness. Other square recipes
can opt in with `compile_minimal_square(H,L; physical_witness3x3=true)` if
their compiled basis covers every witness term. Missing coverage is an error.

To solve with the benchmark driver in an environment with the prototype,
JuMP, Mosek and MosekTools, and a valid MOSEK license:

```sh
OPENBLAS_NUM_THREADS=16 julia --project=<env> --threads=16,0 prototype/pauli/benchmark/minimal_square.jl 16 0 <out> --xxx-witness
```

`--xxx-witness` uses native tolerances of 1e−9 and a strict acceptance gate: `OPTIMAL`,
both `FEASIBLE_POINT`, and all seven original-form metrics (including the
relative objective gap) at most 1e−8. The round TOML reports status acceptance
and the original-form gate separately. Other recipes keep status-only
acceptance unless `--strict` is specified. `--xxx-witness` does not combine with
pricing rounds or `--certify`. The historical strict dual at L = 16
was −0.6790106635450809 per spin, a numerical bound from the SDP, **not** an
exact energy certificate.

Run `julia --project --threads=16,0 prototype/pauli/test_square_xxx16.jl --full`
to check the compiled form without a solver. Allow several minutes and about
25 GiB of RAM for the full build. Omit `--full` for the small checks only.

See [Wang comparison and assumptions](WANG_MATCH.md), [benchmark protocol](benchmark/README.md), and [recorded results](benchmark/RESULTS.md).
