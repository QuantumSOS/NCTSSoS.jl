# Moment Matrix Extraction Feature - Context

## Task Overview
Implement a feature to extract the moment matrix from a `PolyOptResult` model. The extraction method depends on whether the model came from a Moment problem or a Sum-of-Squares (SOS) problem.

## Current Code Structure

### PolyOptResult (src/interface.jl:1-6)
```julia
struct PolyOptResult{T,P,M}
    objective::T
    corr_sparsity::CorrelativeSparsity{P,M}
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
    model::GenericModel{T}
end
```

### MomentProblem (src/moment_solver.jl:3-8)
```julia
struct MomentProblem{T,M,CR<:ConstraintRef,JS<:AbstractJuMPScalar}
    model::GenericModel{T}
    constraints::Vector{CR}
    monomap::Dict{M,JS}
    sa::SimplifyAlgorithm
end
```

### SOSProblem (src/sos_solver.jl:1-3)
```julia
struct SOSProblem{T}
    model::GenericModel{T}
end
```

## Key Findings

1. **Problem Creation Flow** (src/interface.jl:89-96):
   - `moment_relax()` creates a `MomentProblem`
   - If `dualize=true`: `sos_dualize(moment_problem)` creates an `SOSProblem`
   - If `dualize=false`: uses `MomentProblem` directly
   - Only the `model` field is stored in `PolyOptResult`, not the full problem structure

2. **Moment Problem Structure**:
   - Has constraints created via `@constraint(model, moment_mtx in cone)` (moment_solver.jl:104)
   - The `moment_mtx` is a matrix of JuMP expressions (moment_solver.jl:100-103)
   - For moment problems: use `dual()` on constraints to get moment matrix values

3. **SOS Problem Structure**:
   - Has `dual_variables` that are matrix variables representing moment matrices (sos_solver.jl:54-57)
   - Has equality constraints `fα_constraints .== 0` (sos_solver.jl:87)
   - For SOS problems: use `dual()` on equality constraints to get entries

4. **GNS Reconstruction** (src/gns.jl:51-133):
   - The `reconstruct()` function takes a Hankel matrix (moment matrix) as input
   - The moment matrix should be indexed by a monomial basis
   - Current usage requires manually constructing the moment matrix

## User Requirements

1. If `PolyOptResult` has model from Moment problem:
   - Directly use `dual()` function to get the moment matrix

2. If `PolyOptResult` has model from SOS problem:
   - Use `dual()` on equality constraints to get each entry in the moment matrix

## Problem to Solve

Since `PolyOptResult` only stores `model` (not the full `MomentProblem` or `SOSProblem`), we need to:
1. Either distinguish between Moment and SOS problems by inspecting the model
2. Or add metadata to `PolyOptResult` to track which problem type was solved

## Planning Complete - Key Decisions

### 1. Problem Type Tracking
**Decision**: Add `problem_type::ProblemType` field to `PolyOptResult` with values `MomentPrimal` or `SOSDual`.
**Rationale**: More reliable and maintainable than model inspection.

### 2. API Design
```julia
moment_matrix(result::PolyOptResult; clique::Int=1) -> MomentMatrix
moment_matrices(result::PolyOptResult) -> Vector{MomentMatrix}
```
**Return Type**: `MomentMatrix{T,M}` struct containing:
- `matrix::Matrix{T}`: The moment matrix values
- `basis::Vector{M}`: Monomial basis indexing the matrix
- `clique_vars::Vector{Variable}`: Variables in this clique
- `clique_index::Int`: Index of the clique

### 3. Sparsity Handling
- **Correlative Sparsity**: Extract one moment matrix per clique (user selects via `clique` parameter)
- **Term Sparsity**: Reconstruct full block-diagonal matrix from `block_bases`
- **Localizing Matrices**: Future extension (not in initial implementation)

### 4. Extraction Methods
- **Primal (Moment)**: Use `dual()` on PSD constraints to get moment matrix values
- **Dual (SOS)**: Use `dual()` on equality constraints, map back via monomap

### 5. Implementation Approach
- Add `ProblemType` enum and field to `PolyOptResult`
- Create `MomentMatrix` struct in `src/gns.jl`
- Implement extraction functions for primal and dual cases
- Handle term sparsity by reconstructing full basis and assembling blocks
- Comprehensive tests including GNS integration

### 6. Scope Limitations (Initial Version)
- Real problems only (Complex support in future PR)
- Moment matrices only (not localizing matrices)
- Assumes block-diagonal structure for term sparsity

See `plan.md` for complete implementation details.

## Related Code Locations

- `src/interface.jl`: PolyOptResult definition and cs_nctssos function
- `src/moment_solver.jl`: MomentProblem and moment_relax function
- `src/sos_solver.jl`: SOSProblem and sos_dualize function
- `src/gns.jl`: GNS reconstruction that uses moment matrices
