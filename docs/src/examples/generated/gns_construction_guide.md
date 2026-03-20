```@meta
EditURL = "../literate/gns_construction_guide.jl"
```

# [GNS Construction Interface Guide](@id gns-construction-guide)

## What is the GNS construction?

In noncommutative polynomial optimization we solve a semidefinite program (SDP)
that produces **moments** — numerical expectation values ``L(w)`` for every
monomial word ``w`` up to some degree. These moments are abstract numbers; the
**GNS (Gelfand-Naimark-Segal) construction** turns them into concrete matrices
acting on a finite-dimensional Hilbert space that *reproduce* those moments.

Concretely, given the solved moment data the GNS construction returns:

- a matrix ``A_i`` for each variable ``X_i``, and
- a distinguished unit vector ``\xi`` (called the **cyclic vector**),

such that ``\langle \xi,\, w(A_1,\dots,A_n)\, \xi \rangle = L(w)`` for every
word ``w`` in the moment table.

## Key concepts used throughout this page

**Moment matrix (Hankel matrix).**
Choose an ordered list of monomial words ``\{b_1, b_2, \dots\}`` up to some
degree. The **moment matrix** (also called **Hankel matrix**) is the square
matrix whose ``(i,j)`` entry is the moment of the product
``b_i^\dagger \cdot b_j``:
```math
H_{ij} \;=\; L(b_i^\dagger\, b_j),
```
where ``L`` is the linear functional on words defined by the SDP solution.

**Relaxation order, moment degree, and Hankel degree.**
An SDP relaxation of **order** ``d`` produces moments for all words up to
degree ``2d``. To build a Hankel matrix whose rows and columns are indexed by
words of degree ``\le k`` we need moments up to degree ``2k``, because
``\deg(b_i^\dagger b_j) \le 2k``. In this guide we solve at order 3
(moments up to degree 6) and build a Hankel matrix indexed by words up to
degree 3 — we call this degree `H_deg = 3`.

**Principal block and full Hankel.**
The GNS construction partitions the Hankel matrix into two nested levels:

- `full_basis` — all words up to degree `H_deg` (rows/columns of the full
  Hankel matrix),
- `basis` — all words up to degree `hankel_deg` (rows/columns of the
  **principal block**, a sub-matrix of the full Hankel).

`hankel_deg` is always strictly less than `H_deg`; the standard choice is
`hankel_deg = H_deg - 1`. The principal block determines the dimension of
the reconstructed Hilbert space (its numerical rank = the dimension of the
output matrices).

**Flat extension.**
When the numerical rank of the full Hankel matrix equals the rank of its
principal block, we say the Hankel matrix is **flat**. Flatness guarantees
that the GNS quotient space — the vector space spanned by the columns of the
principal block — already captures all the information in the larger matrix.
In other words, the extension rows/columns are linearly dependent on the
principal block, so no new directions appear at higher degree.

## API summary

| Function | Purpose |
|----------|---------|
| `gns_reconstruct` | Full reconstruction → `GNSResult` |
| `reconstruct` | Convenience wrapper → matrices only |
| `test_flatness` | Check the flat-extension rank condition |
| `flat_extend` | Project the extension block to its canonical flat form |
| `robustness_report` | Conditioning and distance-to-flatness diagnostics |
| `verify_gns` | Moment reproduction, symmetry, and constraint checks |

**Prerequisites**: the [polynomial optimization API](@ref polynomial-optimization)
and the [SDP relaxation workflow](@ref moment-sohs-hierarchy).

## Setup

````julia
using NCTSSoS, MosekTools, LinearAlgebra, Logging

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(
    Mosek.Optimizer,
    MOI.Silent() => true,
);
````

## Problem: noncommutative unit ball

We minimize

```math
f \;=\; 2 - X^2 + X Y^2 X - Y^2
```

subject to the ball constraint ``1 - X^2 - Y^2 \ge 0``.
The optimal value is ``f^* = 1``.

````julia
registry, (vars,) = create_noncommutative_variables([("X", 1:2)]);
X, Y = vars;

f = 2.0 - X^2 + X * Y^2 * X - Y^2;
g = 1.0 - X^2 - Y^2;
pop = polyopt(f, registry; ineq_constraints=[g]);
pop
````

````
Optimization Problem (NonCommutativeAlgebra)
────────────────────────────────────
Objective:
    2.0 + -X₁² + -X₂² + X₁X₂²X₁

Equality constraints (0):
    (none)

Inequality constraints (1):
    1 + -X₁² + -X₂² >= 0

Variables (2):
    X₁, X₂

````

## Step 1 — Solve the order-3 moment relaxation

We want a Hankel matrix with `H_deg = 3` and a principal block with
`hankel_deg = 2`. Building the full Hankel requires moments up to degree
``2 \times 3 = 6``, so we solve the SDP at **order 3** (which produces
moments up to degree ``2 \times 3 = 6``).

````julia
solver_config = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=3,
    cs_algo=NoElimination(),
    ts_algo=NoElimination(),
);

sparsity = compute_sparsity(pop, solver_config);
moment_problem = NCTSSoS.moment_relax(
    pop,
    sparsity.corr_sparsity,
    sparsity.cliques_term_sparsities,
);
moment_result = NCTSSoS.solve_moment_problem(moment_problem, SILENT_MOSEK);

@assert isapprox(moment_result.objective, 1.0; atol=1e-5)

(objective=moment_result.objective, n_moments=length(moment_result.monomap))
````

````
(objective = 1.000000000697771, n_moments = 78)
````

## Step 2 — Build the Hankel matrix

We now build the two monomial bases introduced in "Key concepts":

- `full_basis`: all words up to degree `H_deg = 3` — indexes rows and columns
  of the full Hankel matrix.
- `basis`: all words up to degree `hankel_deg = 2` — indexes the principal
  block (a sub-matrix sitting inside the full Hankel).

`NCTSSoS.hankel_matrix` fills the matrix entry-by-entry from the solved
moment table `moment_result.monomap`. It is not exported, so we qualify the
call.

````julia
using NCTSSoS: get_ncbasis

full_basis = get_ncbasis(registry, 3);
basis = get_ncbasis(registry, 2);
hankel = NCTSSoS.hankel_matrix(moment_result.monomap, full_basis);

(hankel_size=size(hankel), basis_sizes=(full=length(full_basis), principal=length(basis)))
````

````
(hankel_size = (15, 15), basis_sizes = (full = 15, principal = 7))
````

## Step 3 — Test flatness

Recall from "Key concepts": the Hankel matrix is **flat** when its full rank
equals the rank of the principal block. `test_flatness` checks this
numerically.

````julia
flatness = test_flatness(hankel, full_basis, basis; atol=1e-6)
````

````
NCTSSoS.FlatnessResult(false, 5, 9, 0.03707312397043376)
````

`FlatnessResult` fields:

- `is_flat`: whether the rank condition holds and the residual is small.
- `rank_principal` / `rank_full`: numerical ranks of the principal block and
  the full Hankel matrix, respectively.
- `err_flat`: measures how well the extension block is predicted by the
  principal block. Internally, the full Hankel is partitioned as
  ```math
  H = \begin{pmatrix} \tilde H & B \\ B^\top & C \end{pmatrix},
  ```
  where ``\tilde H`` is the principal block. The residual is
  ``\|C - Z^\top \tilde H Z\| / (1 + \|C\| + \|Z^\top \tilde H Z\|)``
  with ``\tilde H Z = B``.

SDP solutions are numerical, so the Hankel matrix is only *approximately*
flat. We can canonically project it to the nearest flat extension.

## Step 4 — Flat extension

`flat_extend` replaces the lower-right block ``C`` (see the partition above)
by its least-squares prediction ``Z^\top \tilde H Z`` (where
``\tilde H Z = B``). The result is a matrix whose rank equals
`rank_principal`, i.e. it is exactly flat.

````julia
hankel_flat = flat_extend(hankel, full_basis, basis; atol=1e-6);
````

Verify that the flattened matrix is now exactly flat:

````julia
flatness_after = test_flatness(hankel_flat, full_basis, basis; atol=1e-6)
````

````
NCTSSoS.FlatnessResult(true, 5, 5, 0.0)
````

## Step 5 — GNS reconstruction (SVD method)

`gns_reconstruct` is the main entry point. It accepts three input forms:

1. **Moment map** + registry — most common after `solve_moment_problem`;
   builds the Hankel matrix internally.
2. **Dense Hankel matrix** + registry — use when you want to pre-process the
   Hankel (e.g. flat-extend it first).
3. **Moment map** + `SparsityResult` — for problems with pairwise-disjoint
   variable cliques (sparse GNS).

Here we use form 2 with the flat-extended Hankel. The default `method=:svd`
takes the SVD of the principal block, keeps the columns whose singular values
exceed `atol`, and uses those columns as an orthonormal basis for the
reconstructed Hilbert space. The operator matrix for each variable ``X_i`` is
then computed by projecting the **localizing matrix** (the moment matrix of
``b_k^\dagger X_i b_j``) into this basis.

````julia
gns_svd = with_logger(Logging.SimpleLogger(devnull, Logging.Error)) do
    gns_reconstruct(hankel_flat, registry, 3; method=:svd, hankel_deg=2, atol=1e-6)
end;
````

The returned `GNSResult` stores everything about the reconstruction:

````julia
(
    type = typeof(gns_svd),
    rank = gns_svd.rank,
    full_rank = gns_svd.full_rank,
    n_basis = length(gns_svd.basis),
    n_full_basis = length(gns_svd.full_basis),
    n_singular_values = length(gns_svd.singular_values),
    xi_norm = norm(gns_svd.xi),
)
````

````
(type = NCTSSoS.GNSResult{Float64, Float64, NCTSSoS.NonCommutativeAlgebra, UInt8, NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}}, rank = 5, full_rank = 5, n_basis = 7, n_full_basis = 15, n_singular_values = 7, xi_norm = 1.0000000000181841)
````

### `GNSResult` fields

| Field | Description |
|-------|-------------|
| `matrices` | `Dict{VarIdx, Matrix}` — one operator matrix per variable |
| `xi` | The **cyclic vector**: the unit vector corresponding to the identity word ``\mathbb{1}`` in the orthonormal basis. Moments are recovered as ``L(w) = \langle \xi, w(A)\, \xi \rangle``. |
| `basis` | Monomials of degree ≤ `hankel_deg` (index the principal block) |
| `full_basis` | Monomials of degree ≤ `H_deg` (index the full Hankel) |
| `singular_values` | Singular values of the principal block |
| `rank` | Numerical rank of the principal block (= dimension of reconstructed space) |
| `full_rank` | Numerical rank of the full Hankel |

````julia
A1_svd = gns_svd.matrices[registry[:X₁]];
A2_svd = gns_svd.matrices[registry[:X₂]];

@assert gns_svd.rank == 5
@assert size(A1_svd) == (5, 5)
@assert norm(A1_svd - A1_svd') ≤ 1e-8   # NCTSSoS symmetrizes operators for self-adjoint algebras (NC, Pauli, …)
````

## Step 6 — GNS reconstruction (Cholesky method)

The `:cholesky` method is an alternative to SVD. Instead of using all
singular vectors, it picks a subset of basis words whose columns in the
principal block are numerically linearly independent (via QR with column
pivoting). It then Cholesky-factors the corresponding sub-block. This can be
faster for large problems, but is less tolerant of near-singular blocks.

````julia
gns_chol = with_logger(Logging.SimpleLogger(devnull, Logging.Error)) do
    gns_reconstruct(hankel_flat, registry, 3; method=:cholesky, hankel_deg=2, atol=1e-6)
end;

A1_chol = gns_chol.matrices[registry[:X₁]];
A2_chol = gns_chol.matrices[registry[:X₂]];
````

Both methods reproduce the same moments ``L(w)``. The operator matrices may
differ — they are related by an orthogonal change of basis of the
reconstructed Hilbert space — but all expectation values
``\langle \xi, w(A)\, \xi \rangle`` agree:

````julia
function eval_moment(gns, mono)
    sample = first(values(gns.matrices))
    dim = size(sample, 1)
    T = eltype(sample)
    op = Matrix{T}(I, dim, dim)
    for idx in mono.word
        op = op * gns.matrices[idx]
    end
    return real(dot(gns.xi, op * gns.xi))
end

check_basis = get_ncbasis(registry, 3);
max_diff = maximum(
    abs(eval_moment(gns_svd, m) - eval_moment(gns_chol, m)) for m in check_basis
)
@assert max_diff < 1e-6

(cholesky_rank=gns_chol.rank, max_moment_diff=max_diff)
````

````
(cholesky_rank = 5, max_moment_diff = 7.273553981335112e-11)
````

## Step 7 — Monomap shortcut

When you don't need the flat-extended Hankel explicitly, pass the moment map
directly. `gns_reconstruct` builds the Hankel internally.

````julia
gns_direct = with_logger(Logging.SimpleLogger(devnull, Logging.Error)) do
    gns_reconstruct(moment_result.monomap, registry, 3; hankel_deg=2, atol=1e-6)
end;

@assert gns_direct.rank == 5
````

If you only need the operator matrices (no `xi`, no diagnostics), use the
convenience wrapper `reconstruct` with a Hankel matrix:

````julia
matrices_only = with_logger(Logging.SimpleLogger(devnull, Logging.Error)) do
    reconstruct(hankel_flat, registry, 3; hankel_deg=2, atol=1e-6)
end;

@assert haskey(matrices_only, registry[:X₁])
@assert haskey(matrices_only, registry[:X₂])

(direct_rank=gns_direct.rank, matrices_only_keys=collect(keys(matrices_only)))
````

````
(direct_rank = 5, matrices_only_keys = UInt8[0x05, 0x09])
````

## Step 8 — Robustness diagnostics

`robustness_report` quantifies how well-conditioned the reconstruction is.

````julia
report = robustness_report(gns_svd, hankel, full_basis, basis);
````

Or, using the shorter two-argument form (bases are taken from the `GNSResult`):

````julia
report2 = robustness_report(gns_svd, hankel);
````

`RobustnessReport` fields:

| Field | Meaning |
|-------|---------|
| `sigma_min` | Smallest singular value above the `atol` cutoff (i.e. the weakest direction retained in the basis) |
| `sigma_max` | Largest singular value of the principal block |
| `condition_number` | ``\kappa = \sigma_{\max} / \sigma_{\min}`` — large values signal ill-conditioning |
| `dist_to_flat` | Operator-norm distance between the original Hankel and its canonical flat extension |
| `operator_error_bound` | ``\text{dist\_to\_flat} / \sigma_{\min}`` — first-order bound on how much the reconstructed operators can change due to non-flatness |

````julia
(
    σ_min = report.sigma_min,
    σ_max = report.sigma_max,
    κ = report.condition_number,
    dist_to_flat = report.dist_to_flat,
    op_error_bound = report.operator_error_bound,
)
````

````
(σ_min = 0.12315958149571547, σ_max = 1.501343449377641, κ = 12.190228572917572, dist_to_flat = 0.040074940795970314, op_error_bound = 0.3253903619132098)
````

## Step 9 — Full verification suite

`verify_gns` checks moment reproduction, symmetry, and (optionally)
objective value, constraint feasibility, and ball containment in one call.

````julia
gns = gns_svd   # pick either method
verification = verify_gns(
    gns,
    moment_result.monomap,
    registry;
    poly=f,
    f_star=moment_result.objective,
    constraints=[g],
    ball=true,
    atol=5e-5,
);
````

`VerificationReport` fields:

| Field | Meaning |
|-------|---------|
| `is_symmetric` | All operator matrices satisfy ``\|A_i - A_i^\top\| \le \text{atol}`` |
| `moment_max_error` | Largest ``|L_{\text{model}}(w) - L_{\text{SDP}}(w)|`` over all words in `full_basis` |
| `objective_error` | ``|\langle \xi,\, f(A)\, \xi \rangle - f^*|`` — how closely the model attains the SDP bound |
| `constraint_min_eigenvalues` | Minimum eigenvalue of each constraint polynomial evaluated on the model; non-negative means feasible |
| `ball_contained` | Whether ``I - \sum_i A_i^2 \succeq 0``, i.e. the operators lie inside the unit ball (checked only when `ball=true`) |

````julia
@assert verification.is_symmetric
@assert verification.moment_max_error ≤ 5e-5
@assert verification.objective_error ≤ 5e-5
@assert all(v -> v ≥ -5e-5, verification.constraint_min_eigenvalues)
@assert verification.ball_contained

verification
````

````
NCTSSoS.VerificationReport(true, 7.273570634680482e-11, 1.0130301042465817e-9, [1.6736739682933458e-9], true)
````

## Summary

The typical GNS workflow is:

```
solve SDP  →  (optional) hankel_matrix + test_flatness + flat_extend
           →  gns_reconstruct  →  robustness_report / verify_gns
```

**Key choices:**
- `method=:svd` (default) is more numerically robust; `:cholesky` can be
  faster for large problems.
- `H_deg`: maximum degree of words indexing the full Hankel. `hankel_deg`:
  maximum degree of words indexing the principal block (determines the
  dimension of the output). The standard choice is `hankel_deg = H_deg - 1`.
- `atol`: singular values below this threshold are treated as zero when
  computing the numerical rank.

**See also:**
- [GNS Optimizer Extraction](@ref gns-optimizer-extraction) — a textbook
  example with the noncommutative polydisc.
- [Newton Chip Method](@ref newton-chip-method) — basis reduction for
  ordinary and tracial relaxations.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

