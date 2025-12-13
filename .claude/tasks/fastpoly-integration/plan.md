# Implementation Plan: FastPolynomials Integration (Aggressive Refactor)

## Approach

**Strategy:** Full API redesign exposing `AlgebraType` throughout the optimization layer. Eliminate legacy `Variable` type, use `VariableRegistry` exclusively. Break backward compatibility for cleaner, more Julian design.

## Design Principles

1. **Type-parametric on algebra**: `PolyOpt{A<:AlgebraType, P}` instead of `PolyOpt{P}`
2. **Registry-first**: All variable creation through `VariableRegistry`, no standalone `Variable`
3. **Polynomial-based basis**: Basis is `Vector{Polynomial}` (not Monomial) for algebras with multi-term simplification
4. **DRY simplification**: All simplification delegated to FastPolynomials, remove duplicate logic
5. **Type-stable paths**: Concrete types at call sites, no runtime dispatch

## Key Design Decisions (from Q&A)

### Decision 1: Merge PolyOpt and ComplexPolyOpt
- Single `PolyOpt{A, P}` type handles both real and complex
- Dispatch real vs complex moment problem on `(AlgebraType, CoefficientType)` tuple
- Example: `is_complex_problem(PauliAlgebra, Float64) → true` (Pauli produces i phases)

### Decision 2: Keep Monomial type parameter M in CorrelativeSparsity
- While derivable from P, keeping M explicit improves clarity
- `monomap::Dict{M, JuMP.Variable}` benefits from explicit M

### Decision 3: Sub-registry for cliques
- Implement `subregistry(reg::VariableRegistry{A,T}, indices::Vector{T})`
- Returns new VariableRegistry with filtered variables
- No lazy view needed - simple filter operation

### Decision 4: No ComposedRegistry
- NCTSSoS doesn't use ComposedMonomial currently
- If needed, Cartesian product of separate registries suffices

### Decision 5: Signed indices are separate variables
- For Fermionic/Bosonic: `+i` (annihilation) and `-i` (creation) are DISTINCT
- Graph nodes use position-based mapping, not abs(index)

### Decision 6: Add to FastPolynomials
- `get_ncbasis_monomials(registry, order) → Vector{Monomial}` for convenience
- `variable_indices(poly) → Vector{T}` returns native index type

### Decision 7: Basis contains Polynomials for multi-term algebras
- NonCommutative, Pauli, Projector, Unipotent: 1 word → 1 term (can use Monomial)
- Fermionic, Bosonic: 1 word → potentially multiple terms (must use Polynomial)

## New Type Hierarchy

```julia
# FastPolynomials types (unchanged)
Monomial{A<:AlgebraType, T<:Integer}
Term{M<:AbstractMonomial, C<:Number}
Polynomial{A<:AlgebraType, T<:Integer, C<:Number}
VariableRegistry{A<:AlgebraType, T<:Integer}

# NCTSSoS types (refactored)
OptimizationProblem{A<:AlgebraType, P<:Polynomial{A}}
PolyOpt{A, P} <: OptimizationProblem{A, P}
ComplexPolyOpt{A, P} <: OptimizationProblem{A, P}
CorrelativeSparsity{A, P, M<:Monomial{A}}
TermSparsity{A, M<:Monomial{A}}
MomentProblem{A, T, M<:Monomial{A}}
```

## Implementation Steps

### Phase 1: Core Type Refactoring

#### Step 1.1: Refactor pop.jl [COMPLETE]
- [x] Add algebra type parameter to `OptimizationProblem{A, P}`
- [x] Change `variables::Vector{Variable}` -> `registry::VariableRegistry{A}`
- [x] Remove `comm_gps` (use algebra type for commutation rules)
- [x] Remove `is_unipotent`, `is_projective` (inferred from algebra type)
- [x] Update `polyopt()` signature to accept registry
- [x] Add `_is_complex_problem(::Type{A})` trait for dispatch
- [x] Add backward compatibility accessors (variables, is_unipotent, is_projective, comm_gps)
- [x] Update dependent files with minimal type signature changes (sparse.jl, interface.jl, moment_solver.jl, complex_moment_solver.jl)

**New API:**
```julia
struct PolyOpt{A<:AlgebraType, P<:Polynomial{A}} <: OptimizationProblem{A, P}
    objective::P
    eq_constraints::Vector{P}
    ineq_constraints::Vector{P}
    registry::VariableRegistry{A}
end

function polyopt(objective::Polynomial{A}, registry::VariableRegistry{A};
                 eq_constraints=Polynomial{A}[],
                 ineq_constraints=Polynomial{A}[]) where A
```

#### Step 1.2: Refactor sparse.jl [COMPLETE]
- [x] Add algebra type to `CorrelativeSparsity{A, T, P, M}` (with index type T)
- [x] Replace `cliques::Vector{Vector{Variable}}` -> `cliques::Vector{Vector{T}}` (indices)
- [x] Implement `subregistry(reg, indices)` for clique-local basis generation
- [x] Replace `get_basis(vars, order)` -> `get_ncbasis(sub_registry, order)`
- [x] Add `variable_indices(poly)` and `variable_indices(mono)` -> returns `Set{T}`
- [x] Update graph construction with position-based index mapping
- [x] Handle signed indices using `abs(idx)` for consistency
- [x] Add `extract_monomials_from_basis()` helper
- [x] Add `clique_variables()` backward compatibility function
- [x] Update `PolyOptResult` type signature
- [x] Update `cs_nctssos` to use `variable_indices()` for clique objective
- [x] Update `moment_relax` and `complex_moment_relax` signatures

**Key changes:**
```julia
# Cliques store indices, not Variables
cliques::Vector{Vector{T}}  # T is registry index type

# Sub-registry for clique basis generation
sub_reg = subregistry(registry, clique_indices)
basis = get_ncbasis(sub_reg, order)  # Returns Vector{Polynomial}

# For simple algebras, can extract monomials:
basis_monomials = get_ncbasis_monomials(sub_reg, order)  # New FastPolynomials function

# Graph uses position mapping
all_indices = sorted_unique(variable_indices(objective)...)
idx_to_node = Dict(idx => pos for (pos, idx) in enumerate(all_indices))
```

#### Step 1.3: Refactor & Merge moment_solver.jl + complex_moment_solver.jl [COMPLETE]
- [x] Create single symbolic `MomentProblem{A, T, M, P}` (no JuMP model)
- [x] Store constraints as `Vector{Tuple{Symbol, Matrix{P}}}` (:Zero, :PSD, :HPSD)
- [x] Basis is `Vector{M}` (monomials as basis elements)
- [x] Use existing `_is_complex_problem(A)` trait for dispatch
- [x] Add `solve_moment_problem(mp, optimizer)` interface for direct solving
- [x] Delete complex_moment_solver.jl after merge

**New structure:**
```julia
struct MomentProblem{A<:AlgebraType, T, M, P<:Polynomial{A}}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}  # :Zero, :PSD, :HPSD
    total_basis::Vector{M}
end

# Trait dispatch
is_complex_problem(::Type{A}, ::Type{C}) where {A,C} = C <: Complex
is_complex_problem(::Type{PauliAlgebra}, ::Type{<:Number}) = true

# User interfaces
moment_problem = moment_relax(pop, corr_sparsity, term_sparsities)  # Symbolic
sos_problem = sos_dualize(moment_problem)                           # Common path
result = solve(sos_problem, optimizer)
# OR
result = solve_moment_problem(moment_problem, optimizer)            # Direct (rare)
```

### Phase 2: Algebra Constructor Refactoring [COMPLETE]

#### Step 2.1: Refactor algebra_constructors.jl [COMPLETE]
- [x] Replace `@ncpolyvar` with `create_*_variables` for each algebra
- [x] Remove manual equality constraints (algebra simplification handles Pauli rules, etc.)
- [x] Remove `is_unipotent`, `is_projective`, `comm_gps` from return type
- [x] Return simplified NamedTuple: `(registry=..., var1=..., var2=...)`
- [x] Users add problem-specific constraints via `polyopt(...; eq_constraints=[...])`

**New API:**
```julia
function pauli_algebra(N::Int)
    registry, (σx, σy, σz) = create_pauli_variables(1:N)
    return (registry = registry, σx = σx, σy = σy, σz = σz)
end

function fermionic_algebra(N::Int)
    registry, (a, a†) = create_fermionic_variables(1:N)
    return (registry = registry, a = a, a_dag = a†)
end

function bosonic_algebra(N::Int)
    registry, (c, c†) = create_bosonic_variables(1:N)
    return (registry = registry, c = c, c_dag = c†)
end

function projector_algebra(symbols::Vector{Symbol}, N::Int)
    registry, projs = create_projector_variables(symbols, 1:N)
    return (registry = registry, projectors = projs)
end

function unipotent_algebra(symbols::Vector{Symbol}, N::Int)
    registry, vars = create_unipotent_variables(symbols, 1:N)
    return (registry = registry, variables = vars)
end

function noncommutative_algebra(symbols::Vector{Symbol}, N::Int)
    registry, vars = create_noncommutative_variables(symbols, 1:N)
    return (registry = registry, variables = vars)
end
```

### Phase 3: Interface Refactoring

#### Step 3.1: Refactor interface.jl
- [ ] Add algebra type `A` to `PolyOptResult{A, T, P, M}`
- [ ] Update `cs_nctssos` signature: `cs_nctssos(pop::PolyOpt{A, P}, ...)`
- [ ] Remove `ComplexPolyOpt` special case handling
- [ ] Use `variable_indices()` instead of `variables()`

#### Step 3.2: Refactor sos_solver.jl
- [ ] Add `ProblemKind` type hierarchy for dispatch:
  ```julia
  abstract type ProblemKind end
  struct RealProblem <: ProblemKind end
  struct HermitianProblem <: ProblemKind end
  ```
- [ ] Add `K<:ProblemKind` parameter to `MomentProblem{A, T, M, P, K}`
- [ ] Dispatch `sos_dualize` on problem kind:
  ```julia
  sos_dualize(mp::MomentProblem{A,T,M,P,RealProblem}) = _sos_dualize_real(mp)
  sos_dualize(mp::MomentProblem{A,T,M,P,HermitianProblem}) = _sos_dualize_hermitian(mp)
  ```
- [ ] Remove `SimplifyAlgorithm`, use `symmetric_canon` from FastPolynomials

#### Step 3.3: Refactor gns.jl [COMPLETE]
- [x] Change signature: `reconstruct(H, registry::VariableRegistry{A,TI}, H_deg)`
- [x] Use `get_ncbasis(registry, deg)` + `extract_monomials_from_basis()` instead of `get_basis`
- [x] Return `Dict{TI, Matrix{T}}` mapping variable indices to matrices
- [x] Replace `monomial(var)` with `Monomial{A}([TM(var_idx)])`

### Phase 4: Remove Legacy Code [COMPLETE]

**Decision:** Complete removal, no deprecation warnings or backward compatibility shims.

#### Step 4.1: Clean up utils.jl in FastPolynomials
- [x] Remove `Variable` struct
- [x] Remove `@ncpolyvar` macro
- [x] Remove `get_basis` wrapper (keep `get_ncbasis` only)
- [x] Remove `variables(poly) -> Vector{Variable}`
- [x] Remove `_VAR_INDEX_TO_NAME` and `_VAR_NAME_TO_INDEX` globals
- [x] Keep `AbstractPolynomial` as union type (for flexibility/extensions)

#### Step 4.2: Update NCTSSoS.jl exports
- [x] Remove legacy type exports (`Variable`, `@ncpolyvar`)
- [x] Export new algebra constructors (`pauli_algebra`, `fermionic_algebra`, etc.)
- [x] Update `using .FastPolynomials` imports
- [x] Clean up re-exports

### Phase 5: Test Migration

**Strategy:** Comment out ALL tests in `runtests.jl`. Re-enable incrementally as each phase completes. User will manually verify results.

#### Step 5.0: Prepare test infrastructure
- [ ] Comment out all `include()` statements in `test/runtests.jl`
- [ ] Add comment block explaining incremental re-enablement plan

#### Step 5.1: FastPolynomials tests (should work immediately)
- [ ] `test/fastpoly_test/` - Already uses new API, should pass
- [ ] Suggest: Uncomment `include("fastpoly_test/runtests.jl")`

#### Step 5.2: Core optimization tests (after Phase 1)
- [ ] Update `test/pop.jl` to registry-based API
- [ ] Suggest: Uncomment after Phase 1 Step 1.1 complete

#### Step 5.3: Sparsity tests (after Phase 1)
- [ ] Update `test/sparse.jl` to new clique format
- [ ] Suggest: Uncomment after Phase 1 Step 1.2 complete

#### Step 5.4: Moment solver tests (after Phase 1)
- [ ] Update `test/moment_solver.jl` to unified MomentProblem
- [ ] Suggest: Uncomment after Phase 1 Step 1.3 complete

#### Step 5.5: Interface tests (after Phase 3)
- [ ] Update `test/interface.jl` - end-to-end cs_nctssos
- [ ] Suggest: Uncomment after Phase 3 complete

#### Step 5.6: Advanced features (after Phase 3)
- [ ] `test/sos_solver.jl` - SOS dualization
- [ ] `test/state_poly_opt.jl` - State polynomial optimization
- [ ] `test/trace_poly_opt.jl` - Trace polynomial optimization
- [ ] `test/test_gns.jl` - GNS construction
- [ ] `test/algebra_constructors.jl` - Algebra helpers

#### Test Re-enablement Order
```
Phase 1.1 complete → pop.jl
Phase 1.2 complete → sparse.jl
Phase 1.3 complete → moment_solver.jl
Phase 2 complete   → algebra_constructors.jl
Phase 3 complete   → interface.jl, sos_solver.jl
Phase 4 complete   → All remaining tests
```

## File Migration Order

**Critical Path (must be sequential):**
1. `pop.jl` - Core type definitions
2. `sparse.jl` - Basis generation (depends on pop.jl)
3. `moment_solver.jl` - SDP construction (depends on sparse.jl)
4. `interface.jl` - Orchestration (depends on all above)

**Parallel Work (can be done alongside):**
- `algebra_constructors.jl` (after pop.jl)
- `sos_solver.jl` (after moment_solver.jl)
- `gns.jl` (after sparse.jl)
- `complex_moment_solver.jl` (after moment_solver.jl)

**Final Cleanup:**
- FastPolynomials `utils.jl` legacy removal
- Test migration
- Documentation

## Risk Mitigation

**Risk: Breaking all user code**
- Mitigation: Clear migration guide with before/after examples
- Accept: This is intentional per user request

**Risk: Algebra type mismatch errors**
- Mitigation: Julia's type system catches at compile time
- Benefit: No runtime surprises

**Risk: Performance regression**
- Mitigation: Benchmark before/after each phase
- Monitor: Basis generation, moment matrix construction

**Risk: Test breakage cascade**
- Mitigation: Migrate tests incrementally with each phase
- Validate: Run tests after each file change

## Testing Strategy

### Incremental Re-enablement
- Start with ALL tests in `runtests.jl` commented out
- Re-enable tests as each phase completes (see Phase 5 for order)
- User manually verifies results match expected values

### FastPolynomials Tests (Baseline)
- `make test-FastPoly` should pass throughout refactoring
- These tests use new API and don't depend on NCTSSoS optimization layer

### Numerical Validation
- User will compare results against known solutions
- Tolerance: results within 1e-6 of expected (SDP solver tolerance)
- Known test problems documented in test files

### Performance Benchmarks
- `make bench TARGET=main` after each phase
- Flag any regression > 10%

## Success Criteria

### Phase 1 Complete
- [ ] All core types have algebra parameter
- [ ] pop.jl, sparse.jl, moment_solver.jl refactored
- [ ] Basic optimization works with new types

### Phase 2 Complete
- [ ] All algebra constructors use registry
- [ ] No `@ncpolyvar` usage in NCTSSoS

### Phase 3 Complete
- [ ] Full interface works with new types
- [ ] cs_nctssos end-to-end functional

### Phase 4 Complete
- [ ] Legacy `Variable` removed from FastPolynomials
- [ ] No backward compatibility shims

### Phase 5 Complete
- [ ] All tests pass
- [ ] No performance regression
- [ ] Documentation updated

## Example: Before/After

### Before (Legacy API)
```julia
using NCTSSoS
@ncpolyvar x[1:2] y[1:2] z[1:2]
sys = pauli_algebra(2)
ham = 0.5 * (x[1]*x[2] + y[1]*y[2] + z[1]*z[2])
pop = cpolyopt(ham, sys)
config = SolverConfig(optimizer=Clarabel.Optimizer, order=2)
result = cs_nctssos(pop, config)
```

### After (Registry API)
```julia
using NCTSSoS
sys = pauli_algebra(2)  # Returns registry + σx, σy, σz
σx, σy, σz = sys.σx, sys.σy, sys.σz
ham = 0.5 * (σx[1]*σx[2] + σy[1]*σy[2] + σz[1]*σz[2])
pop = polyopt(ham, sys.registry)
config = SolverConfig(optimizer=Clarabel.Optimizer, order=2)
result = cs_nctssos(pop, config)
```

**Key differences:**
- No `@ncpolyvar` macro
- `pauli_algebra` returns registry directly
- `polyopt` takes registry, not Variable vectors
- Algebra type `PauliAlgebra` propagates through entire pipeline
