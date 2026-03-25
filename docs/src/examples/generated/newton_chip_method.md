```@meta
EditURL = "../literate/newton_chip_method.jl"
```

# [Newton Chip Method](@id newton-chip-method)

This page shows how the Newton chip idea appears in the current `NCTSSoS`
API for both ordinary and tracial noncommutative optimization. For
unconstrained ordinary `PolyOpt` problems, the package exposes
[`newton_chip_basis`](@ref) directly. For tracial problems, the package
exposes the dense trace basis through [`get_state_basis`](@ref), while the
solver still uses `order` to build that basis internally. The background for
these two reductions is developed in
[burgdorfOptimizationPolynomialsNonCommuting2016](@cite) and
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
chip_cfg = SolverConfig(; optimizer=SILENT_COSMO, monomial_basis=chip_basis);
````

`dense_cfg` uses the full order-2 basis. `chip_cfg` injects the custom
Newton-chip basis through `monomial_basis`.

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
3. pass the result through `SolverConfig(monomial_basis=...)`.

---
## Tracial Problems: inspect with `get_state_basis`, solve with `order`

In the tracial setting, the certified polynomial is handled modulo cyclic
equivalence because `tr(AB) = tr(BA)`. The paper's Newton cyclic chip method
prunes the trace basis using the tracial Newton polytope
[klep2022Optimization](@cite). The current package does **not** expose a
public `monomial_basis` hook for trace problems, so the honest workflow is:
inspect the dense trace basis with [`get_state_basis`](@ref), then solve with
`order`.

````julia
trace_registry, (vars,) = create_unipotent_variables([("v", 1:4)]);
alice = vars[1:2];
bob = vars[3:4];
````

`alice` and `bob` are unipotent observables (`U^2 = I`) used in the CHSH
trace problem.

````julia
trace_basis = get_state_basis(trace_registry, 1; state_type=MaxEntangled);
````

`trace_basis` is a `Vector{NCStateWord{MaxEntangled}}`; it is the public
basis constructor for tracial relaxations.

````julia
trace_basis_summary = (
    kind = "NCStateWord{MaxEntangled}",
    size = length(trace_basis),
)
trace_basis_summary
````

````
(kind = "NCStateWord{MaxEntangled}", size = 9)
````

#### Build the CHSH trace problem

````julia
trace_objective = (
    -1.0 * tr(alice[1] * bob[1]) -
    tr(alice[1] * bob[2]) -
    tr(alice[2] * bob[1]) +
    tr(alice[2] * bob[2])
) * one(typeof(alice[1]));
trace_pop = polyopt(trace_objective, trace_registry);
````

`trace_pop` is a tracial `PolyOpt`; use `order`, not `monomial_basis`.

````julia
trace_cfg = SolverConfig(; optimizer=SILENT_COSMO, order=1);
trace_sparsity = compute_sparsity(trace_pop, trace_cfg);
````

`trace_sparsity` shows the internal order-1 basis used by the trace solver.

````julia
@assert trace_sparsity.corr_sparsity.clq_mom_mtx_bases[1] == trace_basis

trace_size_summary = (
    length(trace_basis),
    length(trace_sparsity.corr_sparsity.clq_mom_mtx_bases[1]),
)
trace_size_summary
````

````
(9, 9)
````

#### Solve at order 1

````julia
trace_result = cs_nctssos(trace_pop, trace_cfg);
````

`trace_result` is a `PolyOptResult`; order 1 already recovers the Tsirelson
bound for CHSH in this formulation.

````julia
@assert isapprox(trace_result.objective, -2 * sqrt(2); atol=1e-5)

trace_result.objective
````

````
-2.8284271246604606
````

The tracial workflow is therefore:
1. build the objective with [`tr`](@ref),
2. inspect the dense trace basis with [`get_state_basis`](@ref) if needed,
3. solve with `SolverConfig(order=...)`.

The missing piece, compared with the paper, is a public Newton cyclic chip
helper that would replace the dense trace basis with a tracially pruned one.

## Next steps

- [Trace Polynomial](@ref tracial-polynomial-optimization): larger tracial
  examples with Bell inequalities.
- [`newton_chip_basis`](@ref): API reference for the ordinary Newton-chip
  helper.
- [Quick Start](@ref quick-start): dense and sparse polynomial workflows.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

