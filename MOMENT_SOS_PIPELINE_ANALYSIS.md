# NCTSSoS polynomial optimization → moment → SOS → JuMP lowering analysis

Scope: ordinary polynomial optimization, especially complex Hermitian problems over `PauliAlgebra`, `FermionicAlgebra`, and `BosonicAlgebra`. State/trace moment problems are noted only where they differ.

## Executive take

The current library has a clean symbolic `MomentProblem`, but the JuMP lowering is hard-wired in two places:

1. `MomentProblem -> JuMP` in `src/optimization/moment.jl` is **moment-variable-first**.
   - real problems: one real free variable per canonical moment;
   - complex Hermitian problems: two real free variables `y_re/y_im` per canonical moment;
   - `:HPSD` constraints are immediately realified as `[Re(H) -Im(H); Im(H) Re(H)] in PSDCone()`.

2. `MomentProblem -> SOSProblem` in `src/optimization/sos.jl` is **not symbolic dualization** in the reusable sense. `SOSProblem` already stores a concrete `JuMP.GenericModel`. For complex Hermitian problems it also realifies the Hermitian cones into lifted real PSD variables immediately.

So if we want solver-native lowering choices, the right cut is:

- keep `PolyOpt`, sparsity, and `MomentProblem` symbolic;
- factor JuMP construction out of `solve_moment_problem` / `sos_dualize`;
- add a new symbolic SOS-dual data type, because current `SOSProblem` cannot be re-lowered.

For BPSDP specifically, the production target should be PSD/HPSD-block-first with native `HermitianPSDCone()` when possible. The existing variable-first H2/Nk=2 path creates 18,786 free real moment variables and BPSDP sees 406,706 scalar `1x1` blocks after bridges. That is the thing to kill.

## Source map

Primary files read:

- `src/types/algebra.jl`
- `src/types/registry.jl`
- `src/types/monomial.jl`
- `src/types/polynomial.jl`
- `src/types/arithmetic.jl`
- `src/algorithms/canonicalization.jl`
- `src/algorithms/basis.jl`
- `src/states/word.jl`
- `src/states/polynomial.jl`
- `src/util/helpers.jl`
- `src/optimization/problem.jl`
- `src/optimization/sparsity.jl`
- `src/optimization/moment.jl`
- `src/optimization/sos.jl`
- `src/optimization/interface.jl`
- `demos/SDPALibsdpExport.jl`
- `demos/h2_periodic_nk2_moment_sos.jl`

## Core representation

### Algebra and variable layer

`AlgebraType` is a type-level dispatch tag:

```text
AlgebraType
├─ MonoidAlgebra
│  ├─ NonCommutativeAlgebra
│  ├─ ProjectorAlgebra
│  └─ UnipotentAlgebra
├─ TwistedGroupAlgebra
│  └─ PauliAlgebra
└─ PBWAlgebra
   ├─ FermionicAlgebra
   └─ BosonicAlgebra
```

Storage:

- `VariableRegistry{A,T}`
  - `idx_to_variables::Dict{T,Symbol}`
  - `variables_to_idx::Dict{Symbol,T}`
- `NormalMonomial{A,T}`
  - `word::Vector{T}`
  - identity is empty word;
  - multiplication delegates to algebra-specific `simplify`.
- `Polynomial{A,T,C}`
  - `terms::Vector{Tuple{C,NormalMonomial{A,T}}}`
  - terms sorted, duplicate monomials combined, zero coefficients removed.

Coefficient defaults:

- `PauliAlgebra` defaults to `ComplexF64`.
- Other algebras default to `Float64`, but `moment_relax` promotes Fermionic/Bosonic moment-problem coefficients to complex because `_is_complex_problem(FermionicAlgebra/BosonicAlgebra) == true`.

Adjoints:

- `Polynomial'` conjugates coefficients and adjoints each monomial.
- PBW polynomial adjoint reverses and sign-flips signed indices, then simplifies.
- `NormalMonomial` PBW adjoint alone is deliberately not implemented; polynomial adjoint handles PBW correctly.

### Canonical moment keys

JuMP moment variables are keyed by canonical expectation words, not by `NormalMonomial` objects.

Important functions:

- `expval(m::NormalMonomial)` wraps a monomial in `StateSymbol{Arbitrary}`.
- `symmetric_canon(...)` returns a canonical `Vector{T}`.
- `_sorted_symmetric_basis(xs)` returns sorted unique canonical vectors.

For complex algebras (`Pauli`, `Fermionic`, `Bosonic`), `symmetric_canon(::NormalMonomial)` currently returns `copy(m.word)`. It does **not** identify `w` with `w†`. Instead complex moments are represented by independent real/imag coordinates and Hermitian constraints enforce structure.

## Stage 1: `PolyOpt`

Defined in `src/optimization/problem.jl`:

```julia
struct PolyOpt{A<:AlgebraType,T<:Integer,P<:AbstractPolynomial} <: OptimizationProblem{A,P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    moment_eq_constraints::Vector{P}
    registry::VariableRegistry{A,T}
end
```

Meaning:

- `objective`: polynomial to minimize.
- `eq_constraints`: operator polynomial equalities `p = 0`; later become `:Zero` localizing matrices.
- `ineq_constraints`: operator polynomial inequalities `p >= 0`; later become PSD/HPSD localizing matrices.
- `moment_eq_constraints`: one-sided state-sector constraints `g|ψ⟩ = 0`, modeled as rows `⟨b†g⟩ = 0`; important for fermionic particle-number/trace constraints.
- `registry`: symbol/index map.

Complex Hermitian validation happens here:

- If `_is_complex_problem(A)` is true, objective must be Hermitian.
- Inequality constraints must be Hermitian.
- Equality constraints may be non-Hermitian; they are later split into Hermitian real/imag components.
- `moment_eq_constraints` are not Hermitian-validated; they become scalar zero constraints and are split if needed.

## Stage 2: sparsity

User-facing result in `src/optimization/interface.jl`:

```julia
struct SparsityResult{A,TI,P,M,ST}
    corr_sparsity::CorrelativeSparsity{A,TI,P,M,ST}
    initial_activated_supps::Vector{Vector{M}}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
end
```

### `CorrelativeSparsity`

Defined in `src/optimization/sparsity.jl`:

```julia
struct CorrelativeSparsity{A,T,P,M,ST}
    cliques::Vector{Vector{T}}
    registry::VariableRegistry{A,T}
    cons::Vector{P}
    clq_cons::Vector{Vector{Int}}
    global_cons::Vector{Int}
    clq_mom_mtx_bases::Vector{Vector{M}}
    clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}
end
```

Meaning:

- `cliques`: variable-index cliques.
- `cons`: `vcat(pop.eq_constraints, pop.ineq_constraints)`. `moment_eq_constraints` are used for graph connectivity but are not stored in `cons`.
- `clq_cons`: per-clique indices into `cons`.
- `global_cons`: constraints not contained in any clique.
- `clq_mom_mtx_bases`: per-clique monomial basis for the moment matrix.
- `clq_localizing_mtx_bases`: per-clique, per-constraint basis for localizing matrices.

Basis generation:

- If `SolverConfig(moment_basis=nothing)`, `compute_relaxation_order` chooses order, then `get_ncbasis(subregistry, order)` enumerates monomials.
- If a custom `moment_basis` is supplied, it is normalized, validated, and filtered by clique.
- For PBW algebras, `_basis_subregistry` expands physical modes to include creation and annihilation indices.

### `TermSparsity`

```julia
struct TermSparsity{M}
    term_sparse_graph_supp::Vector{M}
    block_bases::Vector{Vector{M}}
end
```

Meaning:

- `term_sparse_graph_supp`: support activated by term sparsity.
- `block_bases`: basis blocks after term-sparsity graph clique decomposition.

For each clique, `cliques_term_sparsities[i]` is ordered as:

1. moment matrix term sparsity for multiplier `1`;
2. localizing term sparsities for each constraint assigned to that clique.

## Stage 3: `MomentProblem`

Defined in `src/optimization/moment.jl`:

```julia
struct MomentProblem{A<:AlgebraType,T<:Integer,M<:NormalMonomial{A,T},P<:Polynomial{A,T}}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}
    total_basis::Vector{M}
    n_unique_moment_matrix_elements::Int
end
```

This is symbolic. It has no JuMP model.

Fields:

- `objective`: coefficient-promoted objective polynomial.
- `constraints`: polynomial-valued matrix constraints.
  - `:PSD`: real positive semidefinite matrix.
  - `:HPSD`: complex Hermitian positive semidefinite matrix.
  - `:Zero`: polynomial-valued matrix equality.
- `total_basis`: all monomials needed by moment, localizing, global, and one-sided moment equality constraints.
- `n_unique_moment_matrix_elements`: count of unique canonical moments appearing in moment matrices only, not necessarily all of `total_basis`.

### Polynomial optimization → moment problem

Conversion point: `moment_relax(pop, corr_sparsity, cliques_term_sparsities)`.

Key steps:

1. Build `moment_matrix_basis` from the first term-sparsity entry of every clique.
2. Compute `n_unique_moment_matrix_elements = length(_sorted_symmetric_basis(moment_matrix_basis))`.
3. Build `total_basis` via `_polynomial_total_basis(...)` from every actual symbolic constraint matrix entry; include one-sided `moment_eq_constraints` support.
4. Validate objective/eq/ineq support against `total_basis`.
5. Choose cone type:
   - real algebra: `psd_cone = :PSD`;
   - complex algebra: `psd_cone = :HPSD`.
6. Promote coefficient type:
   - complex problem: `MP_C = Complex{typeof(real(zero(C)))}`;
   - real problem: `MP_C = C`.
7. For each clique block, build symbolic matrix with `_build_constraint_matrix(poly, block_basis, cone)`.
8. Append via `_append_constraint!`.
   - For `:Zero` in complex problems, `_zero_constraint_components` splits a non-Hermitian zero matrix into Hermitian part and skew-Hermitian-as-Hermitian part:
     - `(mat + mat') / 2`
     - `(mat - mat') / (2im)`
9. Add fermionic parity superselection zero constraints via `_add_parity_constraints!`.
10. Add one-sided moment equality rows via `_add_moment_eq_constraints!`.

### Constraint matrix entry semantics

`_build_constraint_matrix(poly, local_basis, cone)` constructs:

```text
M[i,j] = Σ conj(c_row) * c_col * coef * simplify(row_word† * mono * col_word)
```

where `row_mono` and `col_mono` may expand over multiple terms for PBW-style normalized monomials.

For the bare moment matrix, `poly = 1`. For localizing constraints, `poly` is the constraint polynomial.

## Stage 4A: current `MomentProblem -> JuMP` direct lowering

Conversion point: `solve_moment_problem(mp, optimizer)`.

Dispatch:

- `_solve_real_moment_problem` if `_is_complex_problem(A) == false`.
- `_solve_complex_moment_problem` otherwise.

### Real path

Storage in JuMP:

- `model = GenericModel{C}()`.
- `basis = sorted unique symmetric canonical words from mp.total_basis`.
- `@variable(model, y[1:length(basis)])`.
- normalization: identity moment equals 1.
- `monovars = Dict(canonical_word => y[i])`.
- each polynomial matrix entry is substituted by `_substitute_poly`.
- constraints:
  - `:Zero`: `jump_mat in Zeros()`;
  - `:PSD`: `jump_mat in PSDCone()`.
- objective is `_substitute_poly(mp.objective, monovars)`.

### Complex Hermitian path

Storage in JuMP:

- `model = GenericModel{real(C)}()`.
- `basis = sorted unique canonical words from mp.total_basis`.
- `@variable(model, y_re[1:n_basis])`.
- `@variable(model, y_im[1:n_basis])`.
- identity normalization:
  - `y_re[one] == 1`
  - `y_im[one] == 0`
- `basis_to_idx = Dict(canonical_word => i)`.
- each polynomial entry is converted by `_substitute_complex_poly` into `(real_affine, imag_affine)` using:

```text
coef = a + ib, y = x + iz
coef*y = (a*x - b*z) + i(b*x + a*z)
```

- `:Zero` constraints become elementwise real and imaginary equalities.
- `:HPSD` constraints are immediately realified:

```julia
embedded = [ Re(H)  -Im(H)
             Im(H)   Re(H) ]
@constraint(model, embedded in PSDCone())
```

- objective uses only the real part of the Hermitian objective.
- solution extraction returns `Dict(canonical_word => Complex(value(y_re[i]), value(y_im[i])))`.

This is the current BPSDP pain point: actual H2/Nk=2 creates 18,786 free real moment variables, and BPSDP sees 406,706 scalar `1x1` blocks after bridging.

## Stage 4B: current `MomentProblem -> SOSProblem`

Conversion point: `sos_dualize(mp)` in `src/optimization/sos.jl`.

Current type:

```julia
struct SOSProblem{T}
    model::GenericModel{T}
    n_unique_elements::Int
end
```

This is already a JuMP model. It is not a symbolic SOS problem.

### Real SOS path

`_sos_dualize_real(mp)`:

- Creates `dual_model = GenericModel{C}()`.
- For each `mp.constraints` entry:
  - `:Zero`: free symmetric matrix variable via `SymmetricMatrixSpace()`;
  - `:PSD`: PSD matrix variable via `PSDCone()`.
- Creates scalar `b`, objective `Max b`.
- Builds `symmetric_basis = _sorted_symmetric_basis(mp.total_basis)`.
- Initializes one coefficient-matching affine expression per canonical moment.
- Adds objective coefficients and subtracts `b` from constant term.
- For each constraint matrix, `get_Cαj(sorted_basis, mat)` maps `(basis_idx,row,col) => coefficient`.
- Adds coefficient-matching equations:

```text
objective_coeff[α] - δ_{α,1} b - Σ_j C_{α,j,row,col} * G_j[row,col] = 0
```

- Adds all constraints with `@constraint(dual_model, fα_constraints .== 0)`.

### Complex Hermitian SOS path

`_sos_dualize_hermitian(mp)`:

- Creates real model `GenericModel{RC}`.
- For each constraint:
  - `:Zero`: requires Hermitian matrix; creates lifted real symmetric free matrix of size `2n x 2n`.
  - `:HPSD`: creates lifted real PSD matrix of size `2n x 2n`.
- Creates scalar `b`, objective `Max b`.
- Builds real and imaginary coefficient equations for every canonical moment.
- For lifted dual matrix `Z`, uses the adjoint of the realification map:

```text
X1 = Z11 + Z22
X2 = Z21 - Z12
G  = X1 + i X2
```

Then for coefficient `c = c_re + i c_im`:

```text
Real(c*G) = c_re*X1 - c_im*X2
Imag(c*G) = c_im*X1 + c_re*X2
```

The code explicitly notes the factor-of-2 trap: summing both diagonal lifted blocks is intentional.

## Stage 5: top-level solving

`cs_nctssos(pop, solver_config; dualize=true)`:

1. `compute_sparsity(pop, solver_config)`
2. `moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)`
3. `solve_sdp(moment_problem, solver_config.optimizer; dualize)`
4. wrap into `PolyOptResult`.

`solve_sdp` either:

- `dualize=true`: calls `sos_dualize`, attaches optimizer, optimizes;
- `dualize=false`: calls `solve_moment_problem`, which builds and solves the moment JuMP model.

## Complex/Hermitian conversion points

Precise points where complex/Hermitian behavior is decided:

1. `polyopt`: validates complex objective and inequality constraints are Hermitian.
2. `_is_complex_problem(A)`: returns true for Pauli, Fermionic, Bosonic.
3. `moment_relax`:
   - chooses `:HPSD` instead of `:PSD`;
   - promotes coefficient type to complex;
   - splits complex non-Hermitian `:Zero` matrices into Hermitian components.
4. `_solve_complex_moment_problem`:
   - introduces `y_re/y_im`;
   - turns every `:HPSD` into real `2n x 2n` `PSDCone()`.
5. `_sos_dualize_hermitian`:
   - turns every Hermitian primal cone into lifted real dual variables;
   - emits real and imaginary coefficient matching equations.
6. `demos/h2_periodic_nk2_moment_sos.jl` duplicates the direct complex moment lowering for labeled diagnostics.

The symbolic `MomentProblem` itself has already marked cones as `:HPSD`, but it has **not** realified them. Realification currently happens only in JuMP model construction.

## PSD cone shape: square vs triangular JuMP constraints

JuMP distinguishes two real PSD encodings:

- `matrix in PSDCone()` with a plain `Matrix` becomes `MOI.PositiveSemidefiniteConeSquare`.
- `Symmetric(matrix) in PSDCone()` becomes `MOI.PositiveSemidefiniteConeTriangle`.

The square cone is not just a storage choice. MOI's `SquareBridge` converts it to:

1. triangular PSD cone on one triangle;
2. explicit equality constraints forcing off-diagonal symmetry.

So wrapping in `Symmetric(...)` is only safe when the matrix is already syntactically symmetric, or when we add those off-diagonal equality constraints ourselves.

Patch implication:

- Real moment path and state moment path can assert symmetry, then use `Symmetric(jump_mat)`.
- Current complex variable-first real-lift should **not** be blindly wrapped. Because complex algebras store `y[w]` and `y[w†]` as separate real/imag coordinates, the realified matrix `[Re(H) -Im(H); Im(H) Re(H)]` is not always syntactically symmetric before the PSD cone imposes conjugacy equalities. Adding those equalities manually and then wrapping in `Symmetric(...)` merely reproduces MOI's `SquareBridge` row burden. That is not a meaningful improvement.

The real fix for complex Hermitian problems is not this micro-patch. It is either native Hermitian modeling with correct star-orbit moment coordinates, or the PSD/HPSD-block-first formulation described below.

## H2/Nk=2 observed pathology

From `output/phase2/h2_nk2_bpsdp_bridge_inspect.json`:

```text
symbolic complex moments        9393
current direct real variables   18786
JuMP variables before bridge    18786
BPSDP blocks after bridge       406724
BPSDP scalar 1x1 blocks         406706
BPSDP non-scalar PSD blocks     18
BPSDP A shape                   333405 x 479158
BPSDP A nnz                     1279400
```

The toy smoke test showed JuMP affine PSD syntax is not automatically the villain. The actual NCTSSoS complex moment lowering/bridge stack is the villain. In particular, the free moment coordinates are toxic for BPSDP at this scale.

## Interface design

### Principle

Separate **symbolic SDP data** from **JuMP lowering**.

Do not hide this in `moment_relax`. `moment_relax` should stay symbolic. The lowering function should be explicit and inspectable.

### Proposed public surface

Use a result wrapper, not a bare model:

```julia
struct ModelBuildResult{T}
    model::JuMP.GenericModel{T}
    problem_kind::Symbol          # :moment or :sos
    formulation::Symbol           # :moment_variables or :psd_blocks, etc.
    representation::Symbol        # :real or :complex (see below). Real algebras: always :real.
    basis
    variable_refs
    block_refs
    eq_refs
    row_labels::Vector{String}
    pivot_map
    extraction
end
```

Possible API:

```julia
build_jump_model(mp::MomentProblem;
    formulation = :moment_variables,
    representation = :real,
    orphan_policy = :error,
    strict = true,
) -> ModelBuildResult
```

For solving:

```julia
solve_moment_problem(mp, optimizer;
    silent = true,
    formulation = :moment_variables,
    representation = :real,
    kwargs...,
)
```

Keep existing defaults unchanged.

Avoid vague `:auto` initially. If added later, it must resolve deterministically and record the decision in `ModelBuildResult`.

#### Why one `representation` knob, not two

A naive design would expose two orthogonal axes:

- `cone ∈ {:real_lift, :native_hermitian}` — shape of the PSD constraint;
- `moment_coordinate ∈ {:split_re_im, :complex_plane}` — shape of each scalar moment variable.

These are logically independent at model construction. In practice only two of the four combinations carry their weight:

| cone | coord | use |
|---|---|---|
| real_lift | split_re_im | current behavior, all-real model |
| native_hermitian | complex_plane | clean Hermitian model when backend supports it |
| native_hermitian | split_re_im | only useful as a workaround for buggy `ComplexPlane` bridges |
| real_lift | complex_plane | pointless: complex scalars get split anyway |

Furthermore, in `formulation = :psd_blocks` mode the moment-coordinate question dissolves entirely — moments are block entries, not free variables — so a second knob would only do work in the formulation we are trying to move away from.

We therefore collapse to a single knob:

- `representation = :real` — real-lift cone + (where applicable) split `y_re/y_im` scalars.
- `representation = :complex` — `HermitianPSDCone()` + (where applicable) JuMP-native complex scalars.

Real algebras ignore this knob. If the buggy-bridge workaround is ever needed, add a third value (e.g. `:complex_cone_real_scalars`) at that time, not preemptively.

### Moment lowering modes

#### 1. `formulation = :moment_variables`

This is current behavior, factored out.

Options (only meaningful for complex algebras; real algebras always behave as `:real`):

- `representation = :real`: realified `[Re(H) -Im(H); Im(H) Re(H)] in PSDCone()` for `:HPSD` blocks, plus split `y_re[w]`/`y_im[w]` scalars per moment. Current complex behavior.
- `representation = :complex`: `LinearAlgebra.Hermitian(affine_matrix) in HermitianPSDCone()` for `:HPSD` blocks, plus JuMP-native complex scalars (`ComplexPlane()`) per moment. Cleaner model when the backend supports it; still leaves free scalar coordinates underneath, so does not by itself fix BPSDP clutter.

This mode is easy to trust because it preserves existing semantics. It does not fix BPSDP scalar clutter unless the bridge stack can eliminate free coordinates.

#### 2. `formulation = :psd_blocks`

This is the important new mode.

Goal:

```text
Declare PSD/HPSD block variables first.
Tie block entries together with linear equalities.
Do not declare a standalone free moment vector.
```

Algorithm for native Hermitian blocks:

1. Scan `mp.constraints`.
2. For every `:HPSD` matrix of size `n`, create:

```julia
@variable(model, X[k][1:n, 1:n] in HermitianPSDCone())
```

For `:PSD`, create real PSD block variables.

3. Discover pivots for canonical moments. Reuse and harden the logic from `demos/SDPALibsdpExport.jl`:

```text
if mat[i,j] = phase * y_α and |phase| = 1,
then define y_α := X_block[i,j] / phase.
```

Pivot metadata must store:

```julia
(canonical_key, block_index, row, col, phase, cone_kind)
```

4. Default orphan policy should be `:error`, not silent fallback. Optional `:aux_psd_free` can allocate auxiliary Hermitian blocks for moments that appear in objective/zero constraints but have no pivot, as the libsdp exporter does.

5. Add normalization:

```text
y_identity = 1
```

through the pivot expression.

6. For each HPSD block entry, add binding constraints:

```text
X_b[i,j] == evaluate_polynomial(mat[i,j], pivot_map)
```

Skip only the exact pivot-defining row when it is provably tautological.

For native Hermitian variables, emit upper-triangle bindings and split scalar complex equality into real and imaginary equalities. Lower-triangle rows are redundant if the symbolic matrix is Hermitian.

7. For every `:Zero` matrix, emit zero equations in terms of pivot expressions. For Hermitian zero matrices, upper triangle is enough; for manually-built non-Hermitian zero matrices, either split first or error.

8. Objective:

```text
minimize real(evaluate_polynomial(mp.objective, pivot_map))
```

9. Extraction callback reconstructs a `monomap` by reading pivot block values and dividing by phase.

For `representation = :real` in `:psd_blocks` mode, each complex Hermitian block is represented by a real PSD variable of size `2n`. Then all four real-lift pieces must be bound:

```text
R[i,j]       =  Re(H[i,j])
R[i,n+j]     = -Im(H[i,j])
R[n+i,j]     =  Im(H[i,j])
R[n+i,n+j]   =  Re(H[i,j])
```

Because the PSD variable is otherwise a generic symmetric real matrix, these equalities are not optional; they are what impose the real-lift structure.

### SOS relowering requires a new symbolic type

Current `SOSProblem` cannot support lowering choices because it is already a JuMP model.

Add a symbolic dual representation, for example:

```julia
struct SOSDualForm{A,T,M,P}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}
    basis::Vector
    coefficient_maps
    n_unique_elements::Int
end
```

Then:

```julia
sos_dual_form(mp::MomentProblem) -> SOSDualForm
build_jump_model(form::SOSDualForm; representation = :real) -> ModelBuildResult
sos_dualize(mp) -> SOSProblem   # compatibility wrapper using current default lowering
```

For complex Hermitian SOS lowering:

- `representation = :real`: reproduce current `_sos_dualize_hermitian` exactly (lifted real symmetric / real PSD `2n × 2n` dual variables).
- `representation = :complex`: use `HermitianPSDCone()` and free Hermitian dual variables directly, but do not implement by hand-waving. The off-diagonal scaling and Hermitian inner product must be specified and tested against the existing factor-of-2 regression before trusting it.

## Backend policy

Recommended explicit choices:

| backend | complex moment default worth testing | reason |
|---|---|---|
| BPSDP | `formulation=:psd_blocks, representation=:complex` | BPSDP MOI supports `HermitianPositiveSemidefiniteConeTriangle`; avoids real doubling and free scalar moment clutter. |
| COSMO / Clarabel | current `formulation=:moment_variables, representation=:real` | conservative compatibility path. |
| SCS | `representation=:complex` worth testing; verify bridge result | SCS has complex cone bridges, but check actual received cone form. |

Do not rely on advertised support alone. The build result and smoke tests should inspect actual backend cone/block metadata where possible.

## Tests needed before implementation is trusted

1. **Structural equivalence**
   - Build current variable-first and new PSD-first models for tiny real and complex `MomentProblem`s.
   - Compare objective/equality maps by evaluating random consistent moments or by coefficient extraction.

2. **Optimization equivalence**
   - Small Pauli/Bosonic/Fermionic Hermitian cases.
   - Existing factor-of-2 test in `test/relaxations/interface.jl` must pass for any native-Hermitian SOS path.

3. **Pivot correctness**
   - Duplicate pivots with phases `1`, `-1`, `im`, `-im`.
   - Moment appearing in multiple blocks.
   - Missing/orphan moments must error by default and report exact keys.

4. **Zero constraint handling**
   - Non-Hermitian equality split into Hermitian components.
   - Hermitian zero matrix emits only independent rows in PSD-first mode.

5. **BPSDP bridge regression**
   - H2/Nk=2 build with `max_iter=0`.
   - Assert scalar `1x1` block count does not scale with moment count.
   - Assert Hermitian block count/dimensions match expected symbolic blocks when `representation=:complex`.

6. **Extraction**
   - From a solved PSD-first model, reconstruct `monomap` and verify all symbolic HPSD matrices evaluate to PSD and all zero constraints have small residual.

## Concrete implementation order

1. Extract current direct moment lowering into `build_jump_model(mp; formulation=:moment_variables, representation=:real)` without changing behavior.
2. Add `ModelBuildResult` metadata and make `solve_moment_problem` call the builder.
3. Factor pivot discovery from `demos/SDPALibsdpExport.jl` into library internals with strict orphan reporting.
4. Implement `formulation=:psd_blocks, representation=:complex` for complex `MomentProblem` first, because that is the BPSDP target.
5. Add H2/Nk=2 bridge-count smoke test on HAI.
6. Only then add `representation=:real` PSD-first if needed for non-Hermitian-capable solvers.
7. Split SOS dualization into symbolic `SOSDualForm` plus lowering. Keep old `SOSProblem` as compatibility wrapper.

## Main pitfalls

- Treating `SOSProblem` as symbolic. It is not.
- Letting `:auto` hide solver-dependent model changes.
- Pivot selection without a formal phase/canonicalization contract.
- Real-lift PSD-first blocks without binding both duplicated real blocks and skew imaginary blocks.
- Assuming native `HermitianPSDCone()` survives MOI bridges for every solver.
- Keeping `_substitute_complex_poly`'s silent skip behavior in new lowering for manually constructed bad `MomentProblem`s. New builders should default to strict missing-moment errors.

## Bottom line

The library already has the right symbolic object for the primal moment side. The missing abstraction is a proper lowering layer. For BPSDP, the correct direction is not more clever bridging of free moments; it is to make PSD/HPSD blocks the primary JuMP variables and express moments as entries of those blocks.

For SOS, there is no comparable symbolic object yet. Build that first, or any “SOS lowering options” will just be lipstick on an already-lowered JuMP model.
