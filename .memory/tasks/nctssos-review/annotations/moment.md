# Review: src/optimization/moment.jl

**Lines**: 841
**Purpose**: Moment relaxation construction and direct solving

## Structure

```
MomentProblem{A,T,M,P}           # Symbolic moment problem (Polynomial)
StateMomentProblem{A,ST,T,M,P}   # Symbolic state moment problem (NCStatePolynomial)

_build_constraint_matrix()       # Build polynomial-valued constraint matrix
moment_relax()                   # Main entry: problem → symbolic constraints
solve_moment_problem()           # Instantiate JuMP model and solve
```

## Key Types

### MomentProblem (lines 46-50)
```julia
struct MomentProblem{A,T,M,P}
    objective::P
    constraints::Vector{Tuple{Symbol, Matrix{P}}}  # (cone, matrix)
    total_basis::Vector{M}
end
```
- `constraints` stores symbolic polynomial matrices, not JuMP expressions
- Cone types: `:Zero`, `:PSD`, `:HPSD`

### StateMomentProblem (lines 590-594)
- Parallel structure for NCStatePolynomial
- Stores `block_basis` with each constraint (needed for SOS dualization)

## Algorithm Flow

### `moment_relax()` (lines 134-195)

1. **Compute total basis**: Union of all `basis[i]† * mono * basis[j]` products
2. **Determine cone**: `_is_complex_problem(A)` → PSD or HPSD
3. **Build constraint matrices**: Loop over cliques/term sparsity blocks
   - Moment matrix: `poly = 1`
   - Localizing: `poly = constraint`
4. **Add parity constraints**: For FermionicAlgebra only

### `solve_moment_problem()` (lines 240-258)

Dispatches to:
- `_solve_real_moment_problem()`: Standard PSD cone
- `_solve_complex_moment_problem()`: Hermitian embedding [Re, -Im; Im, Re]

## Quality Assessment

| Aspect | Rating | Notes |
|--------|--------|-------|
| Documentation | ✓✓✓ | Thorough docstrings |
| Correctness | ✓✓✓ | Hermitian embedding correct |
| Type safety | ✓✓ | `Matrix{Any}` in complex solve (noted in comment) |
| Code structure | ✓✓ | Duplication between Polynomial/NCStatePolynomial |

## Key Observations

### Strengths

1. **Symbolic representation**: Clean separation of problem construction from solving
2. **Hermitian embedding** (lines 363-367):
   ```julia
   [Re(H), -Im(H); Im(H), Re(H)] in PSD
   ```
   Standard and correct.

3. **Parity superselection** (lines 483-565):
   - Fermionic algebras require even-parity operators for non-zero expectations
   - Correctly adds zero constraints for odd-parity entries
   - Well-documented rationale

4. **Canonical lookup** (lines 414-415):
   ```julia
   canon_mono = symmetric_canon(expval(mono))
   haskey(monomap, canon_mono) ? ... : zero
   ```
   Handles missing monomials gracefully.

### Concerns

1. **Type instability** (line 347):
   ```julia
   mat_re = Matrix{Any}(undef, dim, dim)
   ```
   Comment acknowledges this - suggests typed matrices would improve performance.

2. **Missing monomials treated as zero** (lines 411, 460):
   Monomials not in basis get zero expectation - could silently miss terms if
   basis generation is incomplete. Should this warn?

3. **Nested comprehension** (lines 142-151):
   Complex nested iteration for total basis computation. Works but hard to read.

4. **Duplicate code** (~200 lines):
   `MomentProblem` and `StateMomentProblem` implementations nearly identical.
   Could abstract with shared helper functions.

## Data Flow

```
PolyOpt + CorrelativeSparsity + TermSparsity
                │
                ▼
       moment_relax()
                │
                ▼
       MomentProblem (symbolic)
           ┌───┴───┐
           ▼       ▼
  solve_moment_problem()   sos_dualize()
           │                    │
           ▼                    ▼
     JuMP Model           SOSProblem
           │                    │
           ▼                    ▼
      optimize!()          optimize!()
           │                    │
           ▼                    ▼
      objective            objective
```

## Code Quality: Good

Well-structured symbolic abstraction. Main issues are type stability and code duplication - both acknowledged in comments.
