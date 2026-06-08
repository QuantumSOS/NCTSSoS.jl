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

````
Maximize ScalarAffineFunction{Float64}:
 1.0 - 1.0 v[28]

Subject to:

VectorAffineFunction{Float64}-in-Zeros
 ┌               ┐
 │0.0 - 2.0 v[27]│
 └               ┘ ∈ Zeros(1)
 ┌                           ┐
 │1.0 - 1.0 v[21] - 2.0 v[25]│
 └                           ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[19]│
 └               ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 1.0 v[10]│
 └               ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[9]│
 └              ┘ ∈ Zeros(1)
 ┌                           ┐
 │0.0 - 2.0 v[14] - 2.0 v[18]│
 └                           ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[8]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[7]│
 └              ┘ ∈ Zeros(1)
 ┌                                       ┐
 │0.0 - 2.0 v[20] - 2.0 v[23] - 2.0 v[24]│
 └                                       ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[17]│
 └               ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[5]│
 └              ┘ ∈ Zeros(1)
 ┌                           ┐
 │0.0 - 2.0 v[12] - 2.0 v[16]│
 └                           ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 1.0 v[3]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[2]│
 └              ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[26]│
 └               ┘ ∈ Zeros(1)
 ┌              ┐
 │2.0 - 1.0 v[6]│
 └              ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[13]│
 └               ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[4]│
 └              ┘ ∈ Zeros(1)
 ┌                           ┐
 │0.0 - 1.0 v[15] - 2.0 v[22]│
 └                           ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[11]│
 └               ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 1.0 v[1]│
 └              ┘ ∈ Zeros(1)

VectorAffineFunction{Float64}-in-Scaled{PositiveSemidefiniteConeTriangle}
 ┌                              ┐
 │0.0 + 1.0 v[1]                │
 │0.0 + 1.4142135623730951 v[2] │
 │0.0 + 1.0 v[3]                │
 │0.0 + 1.4142135623730951 v[4] │
 │0.0 + 1.4142135623730951 v[5] │
 │0.0 + 1.0 v[6]                │
 │0.0 + 1.4142135623730951 v[7] │
 │0.0 + 1.4142135623730951 v[8] │
 │0.0 + 1.4142135623730951 v[9] │
 │0.0 + 1.0 v[10]               │
 │0.0 + 1.4142135623730951 v[11]│
 │0.0 + 1.4142135623730951 v[12]│
 │0.0 + 1.4142135623730951 v[13]│
 │0.0 + 1.4142135623730951 v[14]│
 │0.0 + 1.0 v[15]               │
 │0.0 + 1.4142135623730951 v[16]│
 │0.0 + 1.4142135623730951 v[17]│
 │0.0 + 1.4142135623730951 v[18]│
 │0.0 + 1.4142135623730951 v[19]│
 │0.0 + 1.4142135623730951 v[20]│
 │0.0 + 1.0 v[21]               │
 │0.0 + 1.4142135623730951 v[22]│
 │0.0 + 1.4142135623730951 v[23]│
 │0.0 + 1.4142135623730951 v[24]│
 │0.0 + 1.4142135623730951 v[25]│
 │0.0 + 1.4142135623730951 v[26]│
 │0.0 + 1.4142135623730951 v[27]│
 │0.0 + 1.0 v[28]               │
 └                              ┘ ∈ Scaled{PositiveSemidefiniteConeTriangle}(PositiveSemidefiniteConeTriangle(7))

Maximize ScalarAffineFunction{Float64}:
 1.0 - 1.0 v[10]

Subject to:

VectorAffineFunction{Float64}-in-Zeros
 ┌              ┐
 │0.0 - 2.0 v[9]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │1.0 - 1.0 v[6]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[4]│
 └              ┘ ∈ Zeros(1)
 ┌                         ┐
 │0.0 - 2.0 v[5] - 2.0 v[7]│
 └                         ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[8]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │2.0 - 1.0 v[1]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[2]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 1.0 v[3]│
 └              ┘ ∈ Zeros(1)

VectorAffineFunction{Float64}-in-Scaled{PositiveSemidefiniteConeTriangle}
 ┌                             ┐
 │0.0 + 1.0 v[1]               │
 │0.0 + 1.4142135623730951 v[2]│
 │0.0 + 1.0 v[3]               │
 │0.0 + 1.4142135623730951 v[4]│
 │0.0 + 1.4142135623730951 v[5]│
 │0.0 + 1.0 v[6]               │
 │0.0 + 1.4142135623730951 v[7]│
 │0.0 + 1.4142135623730951 v[8]│
 │0.0 + 1.4142135623730951 v[9]│
 │0.0 + 1.0 v[10]              │
 └                             ┘ ∈ Scaled{PositiveSemidefiniteConeTriangle}(PositiveSemidefiniteConeTriangle(4))


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
(1.0, 0.9999999999606235, 22, 9)
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

````
Maximize ScalarAffineFunction{Float64}:
 1.0 - 1.0 v[190]

Subject to:

VectorAffineFunction{Float64}-in-Zeros
 ┌                             ┐
 │0.0 - 2.0 v[187] - 2.0 v[189]│
 └                             ┘ ∈ Zeros(1)
 ┌                             ┐
 │0.0 - 2.0 v[186] - 2.0 v[188]│
 └                             ┘ ∈ Zeros(1)
 ┌                                                       ┐
 │0.0 - 2.0 v[169] - 1.0 v[171] - 2.0 v[179] - 2.0 v[185]│
 └                                                       ┘ ∈ Zeros(1)
 ┌                                                                                 ┐
 │0.0 - 2.0 v[152] - 2.0 v[168] - 2.0 v[170] - 2.0 v[177] - 2.0 v[178] - 2.0 v[184]│
 └                                                                                 ┘ ∈ Zeros(1)
 ┌                                                       ┐
 │0.0 - 2.0 v[151] - 1.0 v[153] - 2.0 v[176] - 2.0 v[183]│
 └                                                       ┘ ∈ Zeros(1)
 ┌                                          ┐
 │1.0 - 1.0 v[136] - 2.0 v[175] - 2.0 v[182]│
 └                                          ┘ ∈ Zeros(1)
 ┌                                                       ┐
 │0.0 - 2.0 v[135] - 2.0 v[173] - 2.0 v[174] - 2.0 v[181]│
 └                                                       ┘ ∈ Zeros(1)
 ┌                                          ┐
 │0.0 - 1.0 v[120] - 2.0 v[172] - 2.0 v[180]│
 └                                          ┘ ∈ Zeros(1)
 ┌                                          ┐
 │0.0 - 2.0 v[134] - 2.0 v[161] - 2.0 v[167]│
 └                                          ┘ ∈ Zeros(1)
 ┌                                                                                              ┐
 │0.0 - 2.0 v[119] - 2.0 v[133] - 2.0 v[144] - 2.0 v[150] - 2.0 v[159] - 2.0 v[160] - 2.0 v[166]│
 └                                                                                              ┘ ∈ Zeros(1)
 ┌                                                                                              ┐
 │0.0 - 2.0 v[118] - 2.0 v[132] - 2.0 v[142] - 2.0 v[143] - 2.0 v[149] - 2.0 v[158] - 2.0 v[165]│
 └                                                                                              ┘ ∈ Zeros(1)
 ┌                                                       ┐
 │0.0 - 2.0 v[128] - 2.0 v[131] - 2.0 v[157] - 2.0 v[164]│
 └                                                       ┘ ∈ Zeros(1)
 ┌                                                                                 ┐
 │0.0 - 2.0 v[113] - 2.0 v[126] - 2.0 v[130] - 2.0 v[155] - 2.0 v[156] - 2.0 v[163]│
 └                                                                                 ┘ ∈ Zeros(1)
 ┌                                                       ┐
 │0.0 - 2.0 v[111] - 2.0 v[129] - 2.0 v[154] - 2.0 v[162]│
 └                                                       ┘ ∈ Zeros(1)
 ┌                                          ┐
 │0.0 - 2.0 v[117] - 2.0 v[141] - 2.0 v[148]│
 └                                          ┘ ∈ Zeros(1)
 ┌                                                       ┐
 │0.0 - 2.0 v[116] - 2.0 v[127] - 2.0 v[140] - 2.0 v[147]│
 └                                                       ┘ ∈ Zeros(1)
 ┌                                                                                 ┐
 │0.0 - 2.0 v[112] - 2.0 v[115] - 2.0 v[125] - 2.0 v[138] - 2.0 v[139] - 2.0 v[146]│
 └                                                                                 ┘ ∈ Zeros(1)
 ┌                                                       ┐
 │0.0 - 2.0 v[110] - 2.0 v[114] - 2.0 v[137] - 2.0 v[145]│
 └                                                       ┘ ∈ Zeros(1)
 ┌                ┐
 │0.0 - 2.0 v[124]│
 └                ┘ ∈ Zeros(1)
 ┌                                          ┐
 │0.0 - 2.0 v[109] - 2.0 v[122] - 2.0 v[123]│
 └                                          ┘ ∈ Zeros(1)
 ┌                                          ┐
 │0.0 - 2.0 v[107] - 2.0 v[108] - 2.0 v[121]│
 └                                          ┘ ∈ Zeros(1)
 ┌                ┐
 │0.0 - 2.0 v[106]│
 └                ┘ ∈ Zeros(1)
 ┌                            ┐
 │0.0 - 2.0 v[99] - 1.0 v[105]│
 └                            ┘ ∈ Zeros(1)
 ┌                                                    ┐
 │0.0 - 2.0 v[86] - 2.0 v[97] - 2.0 v[98] - 2.0 v[104]│
 └                                                    ┘ ∈ Zeros(1)
 ┌                                                                            ┐
 │0.0 - 2.0 v[74] - 2.0 v[84] - 2.0 v[85] - 1.0 v[91] - 2.0 v[96] - 2.0 v[103]│
 └                                                                            ┘ ∈ Zeros(1)
 ┌                                                    ┐
 │0.0 - 1.0 v[36] - 2.0 v[63] - 2.0 v[95] - 2.0 v[102]│
 └                                                    ┘ ∈ Zeros(1)
 ┌                                                                ┐
 │0.0 - 2.0 v[34] - 2.0 v[53] - 2.0 v[93] - 2.0 v[94] - 2.0 v[101]│
 └                                                                ┘ ∈ Zeros(1)
 ┌                                                    ┐
 │0.0 - 1.0 v[21] - 2.0 v[44] - 2.0 v[92] - 2.0 v[100]│
 └                                                    ┘ ∈ Zeros(1)
 ┌                                                   ┐
 │0.0 - 2.0 v[72] - 2.0 v[73] - 2.0 v[83] - 2.0 v[90]│
 └                                                   ┘ ∈ Zeros(1)
 ┌                                                               ┐
 │0.0 - 2.0 v[35] - 2.0 v[61] - 2.0 v[62] - 2.0 v[82] - 2.0 v[89]│
 └                                                               ┘ ∈ Zeros(1)
 ┌                                                                                       ┐
 │0.0 - 2.0 v[27] - 2.0 v[33] - 2.0 v[51] - 2.0 v[52] - 2.0 v[80] - 2.0 v[81] - 2.0 v[88]│
 └                                                                                       ┘ ∈ Zeros(1)
 ┌                                                               ┐
 │0.0 - 2.0 v[20] - 2.0 v[42] - 2.0 v[43] - 2.0 v[79] - 2.0 v[87]│
 └                                                               ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[32]│
 └               ┘ ∈ Zeros(1)
 ┌                                       ┐
 │0.0 - 2.0 v[19] - 2.0 v[30] - 2.0 v[31]│
 └                                       ┘ ∈ Zeros(1)
 ┌                                       ┐
 │0.0 - 2.0 v[17] - 2.0 v[18] - 2.0 v[29]│
 └                                       ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[16]│
 └               ┘ ∈ Zeros(1)
 ┌                           ┐
 │0.0 - 2.0 v[71] - 1.0 v[78]│
 └                           ┘ ∈ Zeros(1)
 ┌                                                   ┐
 │0.0 - 1.0 v[28] - 2.0 v[60] - 2.0 v[70] - 2.0 v[77]│
 └                                                   ┘ ∈ Zeros(1)
 ┌                                                               ┐
 │0.0 - 2.0 v[26] - 2.0 v[50] - 2.0 v[68] - 2.0 v[69] - 2.0 v[76]│
 └                                                               ┘ ∈ Zeros(1)
 ┌                                                   ┐
 │0.0 - 1.0 v[15] - 2.0 v[41] - 2.0 v[67] - 2.0 v[75]│
 └                                                   ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[25]│
 └               ┘ ∈ Zeros(1)
 ┌                                       ┐
 │0.0 - 2.0 v[14] - 2.0 v[23] - 2.0 v[24]│
 └                                       ┘ ∈ Zeros(1)
 ┌                                       ┐
 │0.0 - 2.0 v[12] - 2.0 v[13] - 2.0 v[22]│
 └                                       ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 2.0 v[11]│
 └               ┘ ∈ Zeros(1)
 ┌                           ┐
 │0.0 - 2.0 v[59] - 1.0 v[66]│
 └                           ┘ ∈ Zeros(1)
 ┌                                                   ┐
 │0.0 - 2.0 v[49] - 2.0 v[57] - 2.0 v[58] - 2.0 v[65]│
 └                                                   ┘ ∈ Zeros(1)
 ┌                                       ┐
 │0.0 - 2.0 v[40] - 2.0 v[56] - 2.0 v[64]│
 └                                       ┘ ∈ Zeros(1)
 ┌                                       ┐
 │0.0 - 2.0 v[47] - 2.0 v[48] - 1.0 v[55]│
 └                                       ┘ ∈ Zeros(1)
 ┌                                                   ┐
 │0.0 - 2.0 v[38] - 2.0 v[39] - 2.0 v[46] - 2.0 v[54]│
 └                                                   ┘ ∈ Zeros(1)
 ┌                           ┐
 │0.0 - 2.0 v[37] - 1.0 v[45]│
 └                           ┘ ∈ Zeros(1)
 ┌               ┐
 │0.0 - 1.0 v[10]│
 └               ┘ ∈ Zeros(1)
 ┌                         ┐
 │0.0 - 2.0 v[8] - 2.0 v[9]│
 └                         ┘ ∈ Zeros(1)
 ┌                                    ┐
 │1.0 - 1.0 v[3] - 1.0 v[6] - 2.0 v[7]│
 └                                    ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[5]│
 └              ┘ ∈ Zeros(1)
 ┌                         ┐
 │0.0 - 2.0 v[2] - 2.0 v[4]│
 └                         ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 1.0 v[1]│
 └              ┘ ∈ Zeros(1)

VectorAffineFunction{Float64}-in-Scaled{PositiveSemidefiniteConeTriangle}
 ┌                               ┐
 │0.0 + 1.0 v[1]                 │
 │0.0 + 1.4142135623730951 v[2]  │
 │0.0 + 1.0 v[3]                 │
 │0.0 + 1.4142135623730951 v[4]  │
 │0.0 + 1.4142135623730951 v[5]  │
 │0.0 + 1.0 v[6]                 │
 │0.0 + 1.4142135623730951 v[7]  │
 │0.0 + 1.4142135623730951 v[8]  │
 │0.0 + 1.4142135623730951 v[9]  │
 │0.0 + 1.0 v[10]                │
 │0.0 + 1.4142135623730951 v[11] │
 │0.0 + 1.4142135623730951 v[12] │
 │0.0 + 1.4142135623730951 v[13] │
 │0.0 + 1.4142135623730951 v[14] │
 │0.0 + 1.0 v[15]                │
 │0.0 + 1.4142135623730951 v[16] │
 │0.0 + 1.4142135623730951 v[17] │
 │0.0 + 1.4142135623730951 v[18] │
 │0.0 + 1.4142135623730951 v[19] │
 │0.0 + 1.4142135623730951 v[20] │
 │0.0 + 1.0 v[21]                │
 │0.0 + 1.4142135623730951 v[22] │
 │0.0 + 1.4142135623730951 v[23] │
 │0.0 + 1.4142135623730951 v[24] │
 │0.0 + 1.4142135623730951 v[25] │
 │0.0 + 1.4142135623730951 v[26] │
 │0.0 + 1.4142135623730951 v[27] │
 │0.0 + 1.0 v[28]                │
 │0.0 + 1.4142135623730951 v[29] │
 │0.0 + 1.4142135623730951 v[30] │
 │0.0 + 1.4142135623730951 v[31] │
 │0.0 + 1.4142135623730951 v[32] │
 │0.0 + 1.4142135623730951 v[33] │
 │0.0 + 1.4142135623730951 v[34] │
 │0.0 + 1.4142135623730951 v[35] │
 │0.0 + 1.0 v[36]                │
 │0.0 + 1.4142135623730951 v[37] │
 │0.0 + 1.4142135623730951 v[38] │
 │0.0 + 1.4142135623730951 v[39] │
 │0.0 + 1.4142135623730951 v[40] │
 │0.0 + 1.4142135623730951 v[41] │
 │0.0 + 1.4142135623730951 v[42] │
 │0.0 + 1.4142135623730951 v[43] │
 │0.0 + 1.4142135623730951 v[44] │
 │0.0 + 1.0 v[45]                │
 │0.0 + 1.4142135623730951 v[46] │
 │0.0 + 1.4142135623730951 v[47] │
 │0.0 + 1.4142135623730951 v[48] │
 │0.0 + 1.4142135623730951 v[49] │
 │0.0 + 1.4142135623730951 v[50] │
 │0.0 + 1.4142135623730951 v[51] │
 │0.0 + 1.4142135623730951 v[52] │
 │0.0 + 1.4142135623730951 v[53] │
 │0.0 + 1.4142135623730951 v[54] │
 │0.0 + 1.0 v[55]                │
 │0.0 + 1.4142135623730951 v[56] │
 │0.0 + 1.4142135623730951 v[57] │
 │0.0 + 1.4142135623730951 v[58] │
 │0.0 + 1.4142135623730951 v[59] │
 │0.0 + 1.4142135623730951 v[60] │
 │0.0 + 1.4142135623730951 v[61] │
 │0.0 + 1.4142135623730951 v[62] │
 │0.0 + 1.4142135623730951 v[63] │
 │0.0 + 1.4142135623730951 v[64] │
 │0.0 + 1.4142135623730951 v[65] │
 │0.0 + 1.0 v[66]                │
 │0.0 + 1.4142135623730951 v[67] │
 │0.0 + 1.4142135623730951 v[68] │
 │0.0 + 1.4142135623730951 v[69] │
 │0.0 + 1.4142135623730951 v[70] │
 │0.0 + 1.4142135623730951 v[71] │
 │0.0 + 1.4142135623730951 v[72] │
 │0.0 + 1.4142135623730951 v[73] │
 │0.0 + 1.4142135623730951 v[74] │
 │0.0 + 1.4142135623730951 v[75] │
 │0.0 + 1.4142135623730951 v[76] │
 │0.0 + 1.4142135623730951 v[77] │
 │0.0 + 1.0 v[78]                │
 │0.0 + 1.4142135623730951 v[79] │
 │0.0 + 1.4142135623730951 v[80] │
 │0.0 + 1.4142135623730951 v[81] │
 │0.0 + 1.4142135623730951 v[82] │
 │0.0 + 1.4142135623730951 v[83] │
 │0.0 + 1.4142135623730951 v[84] │
 │0.0 + 1.4142135623730951 v[85] │
 │0.0 + 1.4142135623730951 v[86] │
 │0.0 + 1.4142135623730951 v[87] │
 │0.0 + 1.4142135623730951 v[88] │
 │0.0 + 1.4142135623730951 v[89] │
 │0.0 + 1.4142135623730951 v[90] │
 │0.0 + 1.0 v[91]                │
 │0.0 + 1.4142135623730951 v[92] │
 │0.0 + 1.4142135623730951 v[93] │
 │0.0 + 1.4142135623730951 v[94] │
 │0.0 + 1.4142135623730951 v[95] │
 │0.0 + 1.4142135623730951 v[96] │
 │0.0 + 1.4142135623730951 v[97] │
 │0.0 + 1.4142135623730951 v[98] │
 │0.0 + 1.4142135623730951 v[99] │
 │0.0 + 1.4142135623730951 v[100]│
 │0.0 + 1.4142135623730951 v[101]│
 │0.0 + 1.4142135623730951 v[102]│
 │0.0 + 1.4142135623730951 v[103]│
 │0.0 + 1.4142135623730951 v[104]│
 │0.0 + 1.0 v[105]               │
 │0.0 + 1.4142135623730951 v[106]│
 │0.0 + 1.4142135623730951 v[107]│
 │0.0 + 1.4142135623730951 v[108]│
 │0.0 + 1.4142135623730951 v[109]│
 │0.0 + 1.4142135623730951 v[110]│
 │0.0 + 1.4142135623730951 v[111]│
 │0.0 + 1.4142135623730951 v[112]│
 │0.0 + 1.4142135623730951 v[113]│
 │0.0 + 1.4142135623730951 v[114]│
 │0.0 + 1.4142135623730951 v[115]│
 │0.0 + 1.4142135623730951 v[116]│
 │0.0 + 1.4142135623730951 v[117]│
 │0.0 + 1.4142135623730951 v[118]│
 │0.0 + 1.4142135623730951 v[119]│
 │0.0 + 1.0 v[120]               │
 │0.0 + 1.4142135623730951 v[121]│
 │0.0 + 1.4142135623730951 v[122]│
 │0.0 + 1.4142135623730951 v[123]│
 │0.0 + 1.4142135623730951 v[124]│
 │0.0 + 1.4142135623730951 v[125]│
 │0.0 + 1.4142135623730951 v[126]│
 │0.0 + 1.4142135623730951 v[127]│
 │0.0 + 1.4142135623730951 v[128]│
 │0.0 + 1.4142135623730951 v[129]│
 │0.0 + 1.4142135623730951 v[130]│
 │0.0 + 1.4142135623730951 v[131]│
 │0.0 + 1.4142135623730951 v[132]│
 │0.0 + 1.4142135623730951 v[133]│
 │0.0 + 1.4142135623730951 v[134]│
 │0.0 + 1.4142135623730951 v[135]│
 │0.0 + 1.0 v[136]               │
 │0.0 + 1.4142135623730951 v[137]│
 │0.0 + 1.4142135623730951 v[138]│
 │0.0 + 1.4142135623730951 v[139]│
 │0.0 + 1.4142135623730951 v[140]│
 │0.0 + 1.4142135623730951 v[141]│
 │0.0 + 1.4142135623730951 v[142]│
 │0.0 + 1.4142135623730951 v[143]│
 │0.0 + 1.4142135623730951 v[144]│
 │0.0 + 1.4142135623730951 v[145]│
 │0.0 + 1.4142135623730951 v[146]│
 │0.0 + 1.4142135623730951 v[147]│
 │0.0 + 1.4142135623730951 v[148]│
 │0.0 + 1.4142135623730951 v[149]│
 │0.0 + 1.4142135623730951 v[150]│
 │0.0 + 1.4142135623730951 v[151]│
 │0.0 + 1.4142135623730951 v[152]│
 │0.0 + 1.0 v[153]               │
 │0.0 + 1.4142135623730951 v[154]│
 │0.0 + 1.4142135623730951 v[155]│
 │0.0 + 1.4142135623730951 v[156]│
 │0.0 + 1.4142135623730951 v[157]│
 │0.0 + 1.4142135623730951 v[158]│
 │0.0 + 1.4142135623730951 v[159]│
 │0.0 + 1.4142135623730951 v[160]│
 │0.0 + 1.4142135623730951 v[161]│
 │0.0 + 1.4142135623730951 v[162]│
 │0.0 + 1.4142135623730951 v[163]│
 │0.0 + 1.4142135623730951 v[164]│
 │0.0 + 1.4142135623730951 v[165]│
 │0.0 + 1.4142135623730951 v[166]│
 │0.0 + 1.4142135623730951 v[167]│
 │0.0 + 1.4142135623730951 v[168]│
 │0.0 + 1.4142135623730951 v[169]│
 │0.0 + 1.4142135623730951 v[170]│
 │0.0 + 1.0 v[171]               │
 │0.0 + 1.4142135623730951 v[172]│
 │0.0 + 1.4142135623730951 v[173]│
 │0.0 + 1.4142135623730951 v[174]│
 │0.0 + 1.4142135623730951 v[175]│
 │0.0 + 1.4142135623730951 v[176]│
 │0.0 + 1.4142135623730951 v[177]│
 │0.0 + 1.4142135623730951 v[178]│
 │0.0 + 1.4142135623730951 v[179]│
 │0.0 + 1.4142135623730951 v[180]│
 │0.0 + 1.4142135623730951 v[181]│
 │0.0 + 1.4142135623730951 v[182]│
 │0.0 + 1.4142135623730951 v[183]│
 │0.0 + 1.4142135623730951 v[184]│
 │0.0 + 1.4142135623730951 v[185]│
 │0.0 + 1.4142135623730951 v[186]│
 │0.0 + 1.4142135623730951 v[187]│
 │0.0 + 1.4142135623730951 v[188]│
 │0.0 + 1.4142135623730951 v[189]│
 │0.0 + 1.0 v[190]               │
 └                               ┘ ∈ Scaled{PositiveSemidefiniteConeTriangle}(PositiveSemidefiniteConeTriangle(19))

Maximize ScalarAffineFunction{Float64}:
 1.0 - 1.0 v[10]

Subject to:

VectorAffineFunction{Float64}-in-Zeros
 ┌              ┐
 │0.0 - 2.0 v[9]│
 └              ┘ ∈ Zeros(1)
 ┌              ┐
 │1.0 - 1.0 v[6]│
 └              ┘ ∈ Zeros(1)
 ┌                         ┐
 │0.0 - 2.0 v[7] - 2.0 v[8]│
 └                         ┘ ∈ Zeros(1)
 ┌                         ┐
 │0.0 - 2.0 v[4] - 2.0 v[5]│
 └                         ┘ ∈ Zeros(1)
 ┌                         ┐
 │1.0 - 1.0 v[1] - 1.0 v[3]│
 └                         ┘ ∈ Zeros(1)
 ┌              ┐
 │0.0 - 2.0 v[2]│
 └              ┘ ∈ Zeros(1)

VectorAffineFunction{Float64}-in-Scaled{PositiveSemidefiniteConeTriangle}
 ┌                             ┐
 │0.0 + 1.0 v[1]               │
 │0.0 + 1.4142135623730951 v[2]│
 │0.0 + 1.0 v[3]               │
 │0.0 + 1.4142135623730951 v[4]│
 │0.0 + 1.4142135623730951 v[5]│
 │0.0 + 1.0 v[6]               │
 │0.0 + 1.4142135623730951 v[7]│
 │0.0 + 1.4142135623730951 v[8]│
 │0.0 + 1.4142135623730951 v[9]│
 │0.0 + 1.0 v[10]              │
 └                             ┘ ∈ Scaled{PositiveSemidefiniteConeTriangle}(PositiveSemidefiniteConeTriangle(4))


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
(1.0000000001164153, 0.9999999997876812, 57, 7)
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

