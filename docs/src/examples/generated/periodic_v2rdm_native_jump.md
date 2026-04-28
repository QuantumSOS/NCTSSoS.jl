```@meta
EditURL = "../literate/periodic_v2rdm_native_jump.jl"
```

# [Periodic V2RDM: D-Only Bridge Model, Spin-Resolved Native Blocks, and Solver Probes](@id periodic-v2rdm-native-jump)

This page keeps the **honest** bridge formulation for the periodic H₄
`Nk = 2` benchmark:

1. use NCTSSoS to build the fermionic pair basis for `²D`,
2. stop before the full PQG moment-vector route,
3. treat `²D` as the only free PSD variable,
4. build `¹D`, `²Q(D)`, and `²G(D)` manually as affine images in JuMP,
5. keep only the 5 scalar number / trace equalities.

The important refinement on this page is that we do **not** stop at the
paper's `K`-only blocks. The Hamiltonian also conserves spin projection, so we
can split the native `²D`, `²Q`, and `²G` cones further **structurally**:

- `²D`, `²Q` split by `(K, 2S_z(pair))`,
- `²G` splits by `(Δk, Δ(2S_z))`.

That is a real block decomposition, not a pile of extra zero constraints.

See also:

- [Periodic V2RDM Benchmark (H₄ Chain)](@ref periodic-v2rdm-h4)
- [Periodic V2RDM: Building the Paper's PQG Relaxation](@ref periodic-v2rdm-pqg-construction)
- [Periodic V2RDM: PQG Basis, Trace Constraints, and COSMO](@ref periodic-v2rdm-pqg-cosmo)

## Load the native-model helper

The heavy lifting lives under `demos/`. The docs page just reads the block
counts, builds the spin-resolved JuMP model once, and records the reviewed
solver behavior.

````julia
using NCTSSoS
using JuMP

include(joinpath(pkgdir(NCTSSoS), "demos", "H4PeriodicNativeV2RDMHelpers.jl"))
using .H4PeriodicNativeV2RDMHelpers
````

## Step 1 — Build the `²D` basis once, then inspect the two blockings

The helper still starts from the same `²D` pair basis as the earlier K-only
bridge page. The difference is that it now exposes **two** native partitions:

- `:k_only` — the paper's original pair-momentum blocking,
- `:spin_resolved` — pair momentum plus pair-spin / spin-transfer labels.

````julia
data = load_native_problem_data()
k_only = native_size_summary(data; refinement = :k_only)
spin_resolved = native_size_summary(data; refinement = :spin_resolved)

block_partition_summary = (
    k_only_pair_sectors = k_only.pair_sector_rows,
    spin_resolved_pair_sectors = spin_resolved.pair_sector_rows,
    spin_resolved_particle_hole_sectors = spin_resolved.ph_sector_rows,
    d_term_sparsity_k_only = length.(data.D_term_sparsity.block_bases),
    d_term_sparsity_spin_resolved = length.(data.D_term_sparsity_spin.block_bases),
)

@assert k_only.d_only_block_sizes == [240, 256]
@assert spin_resolved.d_only_block_sizes == [56, 128, 56, 64, 128, 64]
@assert spin_resolved.g_block_sizes == [128, 256, 128, 128, 256, 128]

block_partition_summary
````

````
(k_only_pair_sectors = [(K = 0, dim = 240), (K = 1, dim = 256)], spin_resolved_pair_sectors = [(K = 0, two_Sz = -2, dim = 56), (K = 0, two_Sz = 0, dim = 128), (K = 0, two_Sz = 2, dim = 56), (K = 1, two_Sz = -2, dim = 64), (K = 1, two_Sz = 0, dim = 128), (K = 1, two_Sz = 2, dim = 64)], spin_resolved_particle_hole_sectors = [(delta_k = 0, two_Sz_transfer = -2, dim = 128), (delta_k = 0, two_Sz_transfer = 0, dim = 256), (delta_k = 0, two_Sz_transfer = 2, dim = 128), (delta_k = 1, two_Sz_transfer = -2, dim = 128), (delta_k = 1, two_Sz_transfer = 0, dim = 256), (delta_k = 1, two_Sz_transfer = 2, dim = 128)], d_term_sparsity_k_only = [240, 256], d_term_sparsity_spin_resolved = [56, 128, 56, 64, 128, 64])
````

Read that tuple as follows:

- the old native bridge used `²D` / `²Q` blocks `[240, 256]` and `²G` blocks
  `[512, 512]`,
- the spin-resolved bridge splits the pair spaces into
  - `K = 0`: `2S_z = -2, 0, +2` with sizes `[56, 128, 56]`,
  - `K = 1`: `2S_z = -2, 0, +2` with sizes `[64, 128, 64]`,
- and the particle-hole spaces into
  - `Δk = 0`: `Δ(2S_z) = -2, 0, +2` with sizes `[128, 256, 128]`,
  - `Δk = 1`: `Δ(2S_z) = -2, 0, +2` with sizes `[128, 256, 128]`.

The structural caveat is in `²G`.
We **do** split `G^{αβ}` and `G^{βα}` because they carry different spin
transfer, but we **do not** split the zero-transfer sector into separate
`αα` and `ββ` blocks. Those still couple through the `D` term in

```math
G_{pq,st} = δ_{qt} \, {}^{1}D_{ps} - {}^{2}D_{pt,sq}.
```

So the right structural block for `²G` is spin transfer, not naive spin label.

## Step 2 — Compare the three formulations again, now with the spin split

There are still three distinct coordinate systems worth comparing:

1. **current NCTSSoS paper-route** — full PQG moment vector plus lifted
   `moment_eq_constraints`,
2. **D-only moment-vector bridge** — only `²D` basis kept, but entries stored
   as generic complex coordinates,
3. **native Hermitian `²D` model** — only `²D` kept, stored directly as
   Hermitian PSD blocks.

The spin refinement shrinks both the D-only bridge and the native model,
because both inherit their free coordinates from the `²D` block structure.

````julia
current_moment_route = (
    real_variables = 340_546,
    psd_rows = 1_545_187,
    equality_rows = 16_012,
    total_rows = 1_561_199,
)

coordinate_comparison = (
    current_nctssos_real_variables = current_moment_route.real_variables,
    k_only_d_moment_real_variables = k_only.d_only_moment_vector_real_variables,
    spin_resolved_d_moment_real_variables = spin_resolved.d_only_moment_vector_real_variables,
    k_only_native_real_variables = k_only.jump_variables,
    spin_resolved_native_real_variables = spin_resolved.jump_variables,
)

coordinate_comparison
````

````
(current_nctssos_real_variables = 340546, k_only_d_moment_real_variables = 246272, spin_resolved_d_moment_real_variables = 94464, k_only_native_real_variables = 123136, spin_resolved_native_real_variables = 47232)
````

Numerically that means:

- current NCTSSoS paper-route: `340,546` real variables,
- K-only D-only bridge: `246,272` real variables,
- spin-resolved D-only bridge: `94,464` real variables,
- K-only native Hermitian `²D`: `123,136` real variables,
- spin-resolved native Hermitian `²D`: `47,232` real variables.

The spin split therefore buys another clean **`≈ 2.61×`** reduction on the
native free variables relative to the K-only native bridge, and **`≈ 7.21×`**
relative to the current paper-route's direct conic variables.

### Which “V2RDM path” do we mean?

This is where people get sloppy.
There are two different direct-V2RDM comparisons, and only one of them is an
apples-to-apples match to the eliminated native bridge.

1. **Eliminated direct V2RDM**: only `²D` is free; `²Q(D)` and `²G(D)` are
   affine images.
2. **Primal DQG V2RDM**: `²D`, `²Q`, and `²G` are all independent primal PSD
   variables tied together by linear consistency constraints.

The first comparison is exact; the second is only a same-order-of-magnitude
comparison.

````julia
direct_primal_k_only = sum(n^2 for n in k_only.d_only_block_sizes) +
                       sum(n^2 for n in k_only.q_block_sizes) +
                       sum(n^2 for n in k_only.g_block_sizes)
direct_primal_spin_resolved = sum(n^2 for n in spin_resolved.d_only_block_sizes) +
                              sum(n^2 for n in spin_resolved.q_block_sizes) +
                              sum(n^2 for n in spin_resolved.g_block_sizes)

v2rdm_variable_count_comparison = (
    current_nctssos_paper_route = current_moment_route.real_variables,
    eliminated_direct_v2rdm_k_only = k_only.jump_variables,
    eliminated_direct_v2rdm_spin_resolved = spin_resolved.jump_variables,
    primal_dqg_v2rdm_k_only = direct_primal_k_only,
    primal_dqg_v2rdm_spin_resolved = direct_primal_spin_resolved,
)

@assert direct_primal_k_only == 770_560
@assert direct_primal_spin_resolved == 291_072

v2rdm_variable_count_comparison
````

````
(current_nctssos_paper_route = 340546, eliminated_direct_v2rdm_k_only = 123136, eliminated_direct_v2rdm_spin_resolved = 47232, primal_dqg_v2rdm_k_only = 770560, primal_dqg_v2rdm_spin_resolved = 291072)
````

The two honest conclusions are:

- **eliminated vs eliminated:** the spin-resolved native NCTSSoS path and a
  direct eliminated V2RDM path have the **same** free-variable count,
  `47,232`, because they are the same `²D`-only parameterization;
- **full current route vs primal DQG:** the current NCTSSoS paper-route at
  `340,546` variables and a spin-resolved primal DQG V2RDM model at `291,072`
  variables are **on the same order** and differ by only about **17%**.

So if the claim is “the best eliminated NCTSSoS path is on par with an
eliminated V2RDM path,” then the statement is actually stronger than that:
the counts are **identical**.

If the claim is “the current paper-route is already on par with a classical
primal DQG V2RDM model,” that is also defensible — but it is a different
comparison.

## Step 3 — Why the PSD burden drops this time

The K-only bridge killed the equality spam but left the PSD cones almost as
large as the current route.

The spin-resolved bridge is different. It actually shrinks the solver-facing
cones because the sectors themselves are smaller.

````julia
cone_comparison = (
    k_only_hermitian_blocks = k_only.hermitian_block_sizes,
    spin_resolved_hermitian_blocks = spin_resolved.hermitian_block_sizes,
    k_only_real_lift_blocks = k_only.solver_real_lift_block_sizes,
    spin_resolved_real_lift_blocks = spin_resolved.solver_real_lift_block_sizes,
    k_only_psd_rows = k_only.solver_psd_total_rows,
    spin_resolved_psd_rows = spin_resolved.solver_psd_total_rows,
    k_only_total_rows = k_only.solver_total_rows,
    spin_resolved_total_rows = spin_resolved.solver_total_rows,
)

cone_comparison
````

````
(k_only_hermitian_blocks = [240, 256, 240, 256, 512, 512], spin_resolved_hermitian_blocks = [56, 128, 56, 64, 128, 64, 56, 128, 56, 64, 128, 64, 128, 256, 128, 128, 256, 128], k_only_real_lift_blocks = [480, 512, 480, 512, 1024, 1024], spin_resolved_real_lift_blocks = [112, 256, 112, 128, 256, 128, 112, 256, 112, 128, 256, 128, 256, 512, 256, 256, 512, 256], k_only_psd_rows = 1543136, spin_resolved_psd_rows = 584160, k_only_total_rows = 1543141, spin_resolved_total_rows = 584165)
````

So the spin-resolved native model has:

- Hermitian PQG blocks
  `[56, 128, 56, 64, 128, 64, 56, 128, 56, 64, 128, 64, 128, 256, 128, 128, 256, 128]`,
- real-lift PSD blocks
  `[112, 256, 112, 128, 256, 128, 112, 256, 112, 128, 256, 128, 256, 512, 256, 256, 512, 256]`,
- **PSD rows:** `584,160`,
- **total rows including 5 scalar equalities:** `584,165`.

Compared with the K-only native bridge, that is a real **`≈ 2.64×`** drop in
PSD rows. Compared with the current NCTSSoS paper-route, total rows drop by
**`≈ 2.67×`**.

This is the part the K-only bridge could not deliver.
Here the symmetry split hits the actual cone geometry.

## Step 4 — Build the spin-resolved eliminated JuMP model

The helper can now build either refinement. For the docs page we build the
**spin-resolved** model, because that is the interesting artifact.

````julia
built_spin = build_native_jump_model(refinement = :spin_resolved)
spin_model_summary = jump_model_summary(built_spin.model)

spin_native_jump_summary = (
    build_times = built_spin.build_times,
    n_variables = spin_model_summary.n_variables,
    n_constraints = spin_model_summary.n_constraints,
    constraint_counts = spin_model_summary.constraint_counts,
)

@assert spin_model_summary.n_variables == spin_resolved.jump_variables == 47_232
@assert spin_model_summary.n_constraints == spin_resolved.jump_constraints_total == 23

spin_native_jump_summary
````

````
(build_times = (one_rdm = 0.063053042, q_blocks = 0.084957291, g_blocks = 0.297064417, scalar_constraints = 0.005899583, objective = 3.397501542, total = 3.899168584), n_variables = 47232, n_constraints = 23, constraint_counts = Dict("Vector{JuMP.AffExpr} in MathOptInterface.HermitianPositiveSemidefiniteConeTriangle" => 12, "JuMP.AffExpr in MathOptInterface.EqualTo{Float64}" => 5, "Vector{JuMP.VariableRef} in MathOptInterface.HermitianPositiveSemidefiniteConeTriangle" => 6))
````

The native spin-resolved JuMP model therefore contains:

- `6` Hermitian PSD variable blocks for `²D`,
- `12` affine Hermitian PSD constraints for `²Q(D)` and `²G(D)`,
- `5` scalar equalities,
- `47,232` scalar real variables.

That is still a serious SDP, but it is no longer absurd.

## Step 5 — Solver behavior: COSMO vs Mosek

The benchmark scripts are:

- `demos/h4_periodic_native_v2rdm_jump_benchmark.jl`
- `demos/h4_periodic_native_v2rdm_mosek_benchmark.jl`

Reviewed probe outcomes:

````julia
solver_probe_summary = (
    cosmo_k_only = (
        reached_banner = false,
        reached_iterations = false,
        note = "K-only native bridge still stalled before any COSMO banner or iteration output.",
    ),
    cosmo_spin_resolved = (
        log = "demos/results/h4_periodic_native_v2rdm_jump_spin_resolved.log",
        build_wall_s = 56.40,
        reached_banner = false,
        reached_iterations = false,
        note = "Spin-resolved COSMO still printed nothing after the explicit 'begin COSMO output' line before the outer wall budget killed the run.",
    ),
    mosek_spin_resolved = (
        log = "demos/results/h4_periodic_native_v2rdm_mosek_benchmark.log",
        build_wall_s = 57.22,
        reached_banner = true,
        reached_presolve = true,
        reached_iterations = false,
        note = "Mosek parsed the model, printed the conic summary, and completed presolve, but still did not reach a printed interior-point iteration line before the outer wall timeout.",
    ),
)

solver_probe_summary
````

````
(cosmo_k_only = (reached_banner = false, reached_iterations = false, note = "K-only native bridge still stalled before any COSMO banner or iteration output."), cosmo_spin_resolved = (log = "demos/results/h4_periodic_native_v2rdm_jump_spin_resolved.log", build_wall_s = 56.4, reached_banner = false, reached_iterations = false, note = "Spin-resolved COSMO still printed nothing after the explicit 'begin COSMO output' line before the outer wall budget killed the run."), mosek_spin_resolved = (log = "demos/results/h4_periodic_native_v2rdm_mosek_benchmark.log", build_wall_s = 57.22, reached_banner = true, reached_presolve = true, reached_iterations = false, note = "Mosek parsed the model, printed the conic summary, and completed presolve, but still did not reach a printed interior-point iteration line before the outer wall timeout."))
````

The blunt summary:

- **COSMO** still stalls before visible solver setup progress, even after the
  spin refinement cuts the cones to `≈ 5.84e5` real PSD rows.
- **Mosek** behaves materially better: it recognizes the conic model and gets
  through presolve.
- But on the tested wall budget it still did **not** reach the printed
  iteration table.

So Mosek is the first solver here that behaves like a solver rather than a
black hole, but the problem is still not cheap.

Free the large JuMP model before finishing the page.

````julia
built_spin = nothing
GC.gc()
````

## Bottom line

Three conclusions survived contact with the spin-resolved native model.

1. **The honest bridge is still the eliminated native formulation.**
   We keep NCTSSoS for algebra and basis bookkeeping, then switch to direct
   `²D` / `Q(D)` / `G(D)` cones in JuMP.
2. **Spin resolution is finally a structural win, not just aesthetic tidying.**
   Native variables drop from `123,136` to `47,232`, and solver-facing rows
   drop from `1,543,141` to `584,165`.
3. **Mosek is better than COSMO on this shape, but not magically so.**
   COSMO still stalls before visible progress; Mosek at least gets through the
   model summary and presolve.

That is useful engineering progress.
We are finally reducing the right thing.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

