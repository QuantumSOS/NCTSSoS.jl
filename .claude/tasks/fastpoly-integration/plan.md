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

#### Step 1.1: Refactor pop.jl
- [ ] Add algebra type parameter to `OptimizationProblem{A, P}`
- [ ] Change `variables::Vector{Variable}` → `registry::VariableRegistry{A}`
- [ ] Remove `comm_gps` (use algebra type for commutation rules)
- [ ] Remove `is_unipotent`, `is_projective` (inferred from algebra type)
- [ ] Update `polyopt()` signature to accept registry

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

#### Step 1.2: Refactor sparse.jl
- [ ] Add algebra type to `CorrelativeSparsity{A, P, M}` and `TermSparsity{A, M}`
- [ ] Replace `cliques::Vector{Vector{Variable}}` → `cliques::Vector{Vector{T}}` (indices)
- [ ] Implement `subregistry(reg, indices)` for clique-local basis generation
- [ ] Replace `get_basis(vars, order)` → `get_ncbasis(sub_registry, order)`
- [ ] Add `variable_indices(poly)` → returns `Vector{T}` of indices
- [ ] Update graph construction with position-based index mapping
- [ ] Handle signed indices as distinct variables (Fermionic/Bosonic)

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

#### Step 1.3: Refactor & Merge moment_solver.jl + complex_moment_solver.jl
- [ ] Create single symbolic `MomentProblem{A, T, M, P}` (no JuMP model)
- [ ] Store constraints as `Vector{Tuple{Symbol, Matrix{P}}}` (:Zero, :PSD, :HPSD)
- [ ] Basis is `Vector{Polynomial}` (simplified polynomials as basis elements)
- [ ] Add `is_complex_problem(A, C)` trait for dispatch
- [ ] Add `solve_moment_problem(mp, optimizer)` interface for direct solving
- [ ] Delete complex_moment_solver.jl after merge

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

### Phase 2: Algebra Constructor Refactoring

#### Step 2.1: Refactor algebra_constructors.jl
- [ ] Replace `@ncpolyvar` with `create_*_variables` for each algebra
- [ ] Remove manual equality constraints (algebra simplification handles Pauli rules, etc.)
- [ ] Remove `is_unipotent`, `is_projective`, `comm_gps` from return type
- [ ] Return simplified NamedTuple: `(registry=..., var1=..., var2=...)`
- [ ] Users add problem-specific constraints via `polyopt(...; eq_constraints=[...])`

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

#### Step 3.3: Refactor gns.jl
- [ ] Change signature: `reconstruct(H, registry::VariableRegistry{A}, H_deg)`
- [ ] Use `get_ncbasis_monomials(registry, deg)` instead of `get_basis`
- [ ] Return `Vector{Matrix}` ordered by registry iteration
- [ ] Replace `monomial(var)` with `Monomial{A}([var_idx])`

### Phase 4: Remove Legacy Code

#### Step 4.1: Clean up utils.jl in FastPolynomials
- [ ] Remove `Variable` struct
- [ ] Remove `@ncpolyvar` macro
- [ ] Remove `get_basis` wrapper (keep `get_ncbasis` only)
- [ ] Remove `variables(poly) → Vector{Variable}`
- [ ] Keep `AbstractPolynomial` as union type (for flexibility)

#### Step 4.2: Update NCTSSoS.jl exports
- [ ] Remove legacy type exports
- [ ] Export new algebra constructors
- [ ] Update `using .FastPolynomials` imports

### Phase 5: Test Migration

#### Step 5.1: Update test/pop.jl
- [ ] Migrate to registry-based API
- [ ] Test all algebra types

#### Step 5.2: Update test/sparse.jl
- [ ] Test new clique decomposition
- [ ] Verify basis generation

#### Step 5.3: Update test/moment_solver.jl
- [ ] Test moment matrix construction
- [ ] Verify numerical results

#### Step 5.4: Update remaining tests
- [ ] test/sos_solver.jl
- [ ] test/interface.jl
- [ ] test/state_poly_opt.jl
- [ ] test/trace_poly_opt.jl
- [ ] test/test_gns.jl
- [ ] test/algebra_constructors.jl

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

### Per-Phase Testing
- After each step, run affected tests
- Fix breakages before proceeding
- Validate numerical equivalence

### End-to-End Validation
- Known problems with known solutions
- Compare results within tolerance (1e-10)
- Verify all algebra types work

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
