```@meta
EditURL = "../literate/newton_chip_method.jl"
```

# [Newton Chip Method](@id newton-chip-method)

This page shows how the Newton chip idea appears in the current `NCTSSoS`
API for both ordinary and tracial noncommutative optimization. The package
exposes [`newton_chip_basis`](@ref) for unconstrained ordinary `PolyOpt`
problems, and the same helper dispatches to the Newton cyclic chip reduction
for unconstrained tracial `PolyOpt` problems over single-site
`NonCommutativeAlgebra`. The background for these two reductions is developed
in [burgdorfOptimizationPolynomialsNonCommuting2016](@cite) and
[klep2022Optimization](@cite).

**Prerequisites**: familiarity with the
[polynomial optimization API](@ref polynomial-optimization) and the
[trace polynomial example](@ref tracial-polynomial-optimization).

## Setup

We use `COSMO` with `MOI.Silent()` so the generated page stays compact.

````julia
using NCTSSoS, COSMO

const MOI = NCTSSoS.MOI
const SILENT_COSMO = MOI.OptimizerWithAttributes(COSMO.Optimizer, MOI.Silent() => true);

function pretty_basis(basis, registry)
    io = IOBuffer()
    print(io, "[")
    for (i, elem) in pairs(basis)
        i > 1 && print(io, ", ")
        show(IOContext(io, :registry => registry), elem)
    end
    print(io, "]")
    return String(take!(io))
end
````

````
Precompiling packages...
   2329.3 ms  ✓ COSMO
  1 dependency successfully precompiled in 2 seconds. 72 already precompiled.

````

`pretty_basis` returns a registry-aware string for basis inspection.

---
## Eigenvalue Problems: `newton_chip_basis`

In the unconstrained eigenvalue setting, the certified polynomial is
`f - lambda`. The Newton chip method keeps only the right chips of words
`w` whose Hermitian square `w' * w` appears in the support of that
certificate [burgdorfOptimizationPolynomialsNonCommuting2016](@cite).

The package exposes that reduction through [`newton_chip_basis`](@ref).

````julia
registry, (x,) = create_noncommutative_variables([("x", 1:2)]);
````

`registry` is a `VariableRegistry`; `x` stores two noncommuting operators.

````julia
objective = 1.0 * x[1]^2 + 2.0 * x[2] * x[1] * x[1] * x[2] + 1.0;
pop = polyopt(objective, registry);
````

`pop` is an unconstrained ordinary `PolyOpt`; this is the scope accepted by
`newton_chip_basis`.

````julia
dense_basis = get_ncbasis(registry, 2);
chip_basis = newton_chip_basis(pop, 2);
````

`dense_basis` and `chip_basis` are vectors of monomials. The Newton-chip
basis keeps only the chips needed by the support of `objective`.

````julia
basis_summary = (
    dense = pretty_basis(dense_basis, registry),
    chip = pretty_basis(chip_basis, registry),
    sizes = (length(dense_basis), length(chip_basis)),
)
basis_summary
````

````
(dense = "[𝟙, x₁, x₂, x₁², x₁x₂, x₂x₁, x₂²]", chip = "[𝟙, x₁, x₂, x₁x₂]", sizes = (7, 4))
````

#### Configure dense and custom relaxations

````julia
dense_cfg = SolverConfig(; optimizer=SILENT_COSMO, order=2);
chip_cfg = SolverConfig(; optimizer=SILENT_COSMO, moment_basis=chip_basis);
````

`dense_cfg` uses the full order-2 basis. `chip_cfg` injects the custom
Newton-chip basis through `moment_basis`.

````julia
dense_sparsity = compute_sparsity(pop, dense_cfg);
chip_sparsity = compute_sparsity(pop, chip_cfg);
````

Each `SparsityResult` stores the moment-matrix basis used by the solver.

````julia
@assert chip_sparsity.corr_sparsity.clq_mom_mtx_bases[1] == chip_basis

ordinary_size_summary = (
    length(dense_sparsity.corr_sparsity.clq_mom_mtx_bases[1]),
    length(chip_sparsity.corr_sparsity.clq_mom_mtx_bases[1]),
)
ordinary_size_summary
````

````
(7, 4)
````

#### Solve both relaxations

````julia
dense_result = cs_nctssos(pop, dense_cfg);
chip_result = cs_nctssos(pop, chip_cfg);
````

Each result is a `PolyOptResult`; compare the bound and the SDP size.

````julia
@assert isapprox(chip_result.objective, dense_result.objective; atol=1e-6)

(
    dense_result.objective,
    chip_result.objective,
    dense_result.n_unique_moment_matrix_elements,
    chip_result.n_unique_moment_matrix_elements,
)
````

````
(1.0000000001164153, 1.0000000001164153, 22, 9)
````

The ordinary workflow is therefore:
1. build an unconstrained `PolyOpt`,
2. call [`newton_chip_basis`](@ref),
3. pass the result through `SolverConfig(moment_basis=...)`.

---
## Tracial Problems: `newton_chip_basis`

In the tracial setting, the certified polynomial is handled modulo cyclic
equivalence because `tr(AB) = tr(BA)`. For unconstrained scalar trace
objectives `tr(f)` over single-site `NonCommutativeAlgebra`,
[`newton_chip_basis`](@ref) dispatches to the Newton cyclic chip reduction and
returns a reduced `Vector{NCStateWord{MaxEntangled}}` that can be injected
through `moment_basis` [klep2022Optimization](@cite).

The helper is intentionally narrower than the full trace-polynomial API: it
covers the free, unconstrained trace setting from the theorem, not products of
traces, mixed state/operator objectives, built-in algebraic relations like
`U^2 = I`, or constrained trace problems.

````julia
trace_registry, (y,) = create_noncommutative_variables([("y", 1:2)]);
````

`trace_registry` stores a single-site free-algebra trace problem.

````julia
trace_objective = tr(1.0 * y[1] * y[1] + 1.0 * y[2] * y[1] * y[1] * y[2] + 1.0) * one(typeof(y[1]));
trace_pop = polyopt(trace_objective, trace_registry);
````

`trace_pop` is an unconstrained tracial `PolyOpt`, so it is within the scope
accepted by [`newton_chip_basis`](@ref).

````julia
dense_trace_basis = get_state_basis(trace_registry, 2; state_type=MaxEntangled);
chip_trace_basis = newton_chip_basis(trace_pop, 2);
````

`dense_trace_basis` is the full order-2 trace basis. `chip_trace_basis`
keeps only the pure operator words needed by the tracial Newton polytope.

````julia
trace_chip_words = [elem.nc_word for elem in chip_trace_basis];
trace_basis_summary = (
    dense_size = length(dense_trace_basis),
    chip_words = pretty_basis(trace_chip_words, trace_registry),
    chip_size = length(chip_trace_basis),
)
trace_basis_summary
````

````
(dense_size = 19, chip_words = "[𝟙, y₁, y₁y₂, y₂y₁]", chip_size = 4)
````

#### Configure dense and custom relaxations

````julia
dense_trace_cfg = SolverConfig(; optimizer=SILENT_COSMO, order=2);
chip_trace_cfg = SolverConfig(; optimizer=SILENT_COSMO, moment_basis=chip_trace_basis);
````

`dense_trace_cfg` builds the full order-2 trace basis. `chip_trace_cfg`
injects the Newton cyclic chip basis through `moment_basis`.

For state/trace problems the package now validates custom bases before solving:
the injected basis must generate every objective/constraint moment that the
relaxation needs. If you choose `d` too small, `compute_sparsity`/`cs_nctssos`
will throw instead of silently dropping trace words from the SDP.

````julia
dense_trace_sparsity = compute_sparsity(trace_pop, dense_trace_cfg);
chip_trace_sparsity = compute_sparsity(trace_pop, chip_trace_cfg);
````

The custom trace basis is stored in the same `SparsityResult` structure used
by the ordinary polynomial solver.

````julia
@assert chip_trace_sparsity.corr_sparsity.clq_mom_mtx_bases[1] == chip_trace_basis

trace_size_summary = (
    length(dense_trace_sparsity.corr_sparsity.clq_mom_mtx_bases[1]),
    length(chip_trace_sparsity.corr_sparsity.clq_mom_mtx_bases[1]),
)
trace_size_summary
````

````
(19, 4)
````

#### Solve both relaxations

````julia
dense_trace_result = cs_nctssos(trace_pop, dense_trace_cfg);
chip_trace_result = cs_nctssos(trace_pop, chip_trace_cfg);
````

The Newton cyclic chip basis is a basis reduction only; the relaxation value
should match the dense order-2 tracial solve.

````julia
@assert isapprox(chip_trace_result.objective, dense_trace_result.objective; atol=1e-6)

(
    dense_trace_result.objective,
    chip_trace_result.objective,
    dense_trace_result.n_unique_moment_matrix_elements,
    chip_trace_result.n_unique_moment_matrix_elements,
)
````

````
(0.9999999998835847, 1.0000000001164153, 57, 7)
````

The tracial workflow is therefore:
1. build an unconstrained scalar trace objective with [`tr`](@ref),
2. call [`newton_chip_basis`](@ref),
3. pass the result through `SolverConfig(moment_basis=...)`.

For general trace polynomials or trace problems with relations/constraints,
fall back to the dense trace basis via [`get_state_basis`](@ref) and solve
with `order`.

## Next steps

- [Trace Polynomial](@ref tracial-polynomial-optimization): larger tracial
  examples, including Bell inequalities and nonlinear trace objectives.
- [`newton_chip_basis`](@ref): API reference for both the ordinary Newton-chip
  and tracial Newton cyclic chip helpers.
- [Quick Start](@ref quick-start): dense and sparse polynomial workflows.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

