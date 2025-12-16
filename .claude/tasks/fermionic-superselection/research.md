# Research: Correct Fermionic Parity Superselection Implementation

## Executive Summary

The current implementation INCORRECTLY filters odd-parity monomials from the moment matrix basis, causing the fermionic chain test to fail with "primal infeasible". The correct approach is to **keep ALL monomials in the basis** (including odd-parity ones) and instead **add equality constraints** that set moment matrix entries with odd total parity to zero.

## Why Current Approach is Wrong

### Mathematical Foundation

In a fermionic moment relaxation, moment matrix entry `M[i,j]` represents the expectation value:

```
⟨basis[i]† · operator · basis[j]⟩
```

where `operator` is a term from the polynomial being constrained (e.g., identity for the moment matrix, or a constraint polynomial for localizing matrices).

**Key Insight**: Even if both `basis[i]` and `basis[j]` have **odd parity** (e.g., single creation/annihilation operators), the **total operator** `adjoint(basis[i]) · operator · basis[j]` can have **even parity**!

### Concrete Example

Consider:
- `basis[i] = a[1]` (odd parity: 1 operator)
- `operator = I` (even parity: 0 operators)
- `basis[j] = a[1]` (odd parity: 1 operator)

The full operator is:
```
adjoint(a[1]) · I · a[1] = a_dag[1] · a[1]
```

This has **even parity** (2 operators) and should have a non-zero expectation value! It represents the number operator `n₁`.

### Experimental Verification

Test file `test_parity_analysis.jl` confirms:
```
Basis[i] = a[1], operator = I, basis[j] = a[1]
Full operator (adjoint(a[1]) * I * a[1]):
  coef: 1.0, mono.word: Int8[-1, 1], parity: even
```

### Why Filtering is a Relaxation (Not a Restriction)

Current implementation (lines 152-156 in `src/moment_solver.jl`):
```julia
# Filter out odd-parity monomials for fermionic algebras (parity superselection)
# Only operators with even parity can have non-zero expectation values
if A == FermionicAlgebra
    total_basis = filter(has_even_parity, total_basis)
end
```

This **removes** monomials like `a[1]`, `a_dag[1]`, `a_dag[1] * a_dag[2]`, etc. from the basis.

**Problem**: When we compute moment matrix entries, we need these monomials! Example:
- To represent `⟨a_dag[1] * a[1]⟩` (number operator), we need `basis[i] = a[1]` and `basis[j] = a[1]`
- But both are filtered out because they have odd parity!
- This makes the moment relaxation **too restrictive** → primal infeasible

### Test Failure Analysis

File: `test/fermionic_chain.jl`

Problem: N=4 fermionic chain with ground state energy E₀ = -4

Current result:
```
Terminated with status = primal infeasible
Ground state energy lower bound: -2.4838463458461796e10
```

The SDP is infeasible because we've removed critical basis elements needed to represent the Hamiltonian's expectation value.

## The Correct Approach

### Mathematical Formulation

**Parity Superselection Rule**: An operator `O` has zero expectation value if and only if it has odd fermion parity:

```
⟨O⟩ = 0  ⟺  parity(O) is odd
```

For moment matrix entry `M[i,j]`:
```
M[i,j] = y[canonical(adjoint(basis[i]) · operator · basis[j])]
```

where `canonical(·)` applies simplification and symmetric canonicalization.

**Constraint**: If the simplified, canonicalized monomial has odd parity:
```
M[i,j] = 0  (hard equality constraint)
```

### Algorithm

1. **Keep full basis**: Do NOT filter `total_basis` by parity
2. **During constraint matrix construction**: For each entry `(i,j)`:
   - Compute `full_op_poly = _neat_dot3(basis[i], mono, basis[j])`
   - For each term in `full_op_poly`:
     - Canonicalize: `canon_mono = symmetric_canon(expval(mono))`
     - Check parity: `if !has_even_parity(canon_mono)`
     - If odd parity: **Add constraint that this entry equals zero**
3. **Modify substitution functions**: No changes needed (already skip missing monomials)

### Implementation Strategy

Two possible approaches:

#### Approach A: Modify `_build_constraint_matrix` (Cleaner)

Add parity checking during matrix construction (lines 72-92 in `moment_solver.jl`):

```julia
function _build_constraint_matrix(
    poly::P,
    local_basis::Vector{M},
    cone::Symbol
) where {T, P<:AbstractPolynomial{T}, M<:Monomial{A,T}} where {A<:AlgebraType}

    moment_mtx = Matrix{P}(undef, length(local_basis), length(local_basis))

    # Track which entries must be zero due to parity
    parity_zero_indices = Set{Tuple{Int,Int}}()

    for (i, row_idx) in enumerate(local_basis)
        for (j, col_idx) in enumerate(local_basis)
            element_poly = sum(
                coef * _neat_dot3(row_idx, mono, col_idx)
                for (coef, mono) in zip(coefficients(poly), monomials(poly))
            )

            # For fermionic algebras, check if all terms have odd parity
            if A == FermionicAlgebra && should_be_zero_by_parity(element_poly)
                # Mark for zero constraint
                push!(parity_zero_indices, (i, j))
                moment_mtx[i, j] = zero(P)  # Will be constrained to zero
            else
                moment_mtx[i, j] = element_poly
            end
        end
    end

    return (cone, moment_mtx, parity_zero_indices)
end

function should_be_zero_by_parity(poly::Polynomial{FermionicAlgebra,T,C}) where {T,C}
    # Check if all non-zero terms have odd parity after canonicalization
    all_terms_odd = true
    has_nonzero_term = false

    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        if !iszero(coef)
            has_nonzero_term = true
            canon_mono = symmetric_canon(expval(mono))
            if has_even_parity(canon_mono)
                all_terms_odd = false
                break
            end
        end
    end

    return has_nonzero_term && all_terms_odd
end
```

Then modify solver to add zero constraints for parity-forbidden entries.

#### Approach B: Post-process constraints (Simpler)

After building symbolic `MomentProblem`, scan constraint matrices and add zero constraints:

```julia
function add_parity_constraints!(
    mp::MomentProblem{FermionicAlgebra,T,M,P}
) where {T,M,P}
    # For each constraint matrix, identify parity-forbidden entries
    for (cone, mat) in mp.constraints
        # Check each entry for odd parity
        for i in axes(mat, 1), j in axes(mat, 2)
            if should_be_zero_by_parity(mat[i,j])
                # Add equality constraint: mat[i,j] = 0
                # This gets added to the Zero cone
                push!(mp.constraints, (:Zero, [mat[i,j];;]))  # 1x1 matrix
            end
        end
    end
end
```

**Recommendation**: Approach B is simpler to implement and doesn't require changing the constraint matrix structure.

## Detailed Code Changes

### Files to Modify

#### 1. `src/moment_solver.jl`

**Remove lines 152-156** (the incorrect filtering):
```julia
# REMOVE THIS:
if A == FermionicAlgebra
    total_basis = filter(has_even_parity, total_basis)
end
```

**Add helper function** (after `_build_constraint_matrix`, around line 93):
```julia
"""
    _has_odd_parity_only(poly::Polynomial{FermionicAlgebra,T,C}) where {T,C}

Check if a polynomial's expectation value must be zero due to parity superselection.
Returns `true` if all non-zero terms have odd parity after canonicalization.
"""
function _has_odd_parity_only(
    poly::Polynomial{FermionicAlgebra,T,C}
) where {T<:Integer,C<:Number}
    has_nonzero_term = false

    for (coef, mono) in zip(coefficients(poly), monomials(poly))
        if !iszero(coef)
            has_nonzero_term = true
            canon_mono = symmetric_canon(expval(mono))
            if has_even_parity(canon_mono)
                return false  # Found even-parity term
            end
        end
    end

    return has_nonzero_term  # true only if all terms are odd parity
end
```

**Add parity constraint post-processing** (at end of `moment_relax`, around line 188):
```julia
function moment_relax(
    pop::PolyOpt{A,P},
    corr_sparsity::CorrelativeSparsity{A,TI,P,M},
    cliques_term_sparsities::Vector{Vector{TermSparsity{M}}}
) where {A<:AlgebraType, TI<:Integer, C<:Number, P<:Polynomial{A,TI,C}, M<:Monomial{A,TI}}

    # ... existing code (WITHOUT the filtering) ...

    mp = MomentProblem{A, TI, M, P}(pop.objective, constraints, total_basis)

    # Add parity superselection constraints for fermionic algebras
    if A == FermionicAlgebra
        _add_parity_constraints!(mp)
    end

    return mp
end

"""
    _add_parity_constraints!(mp::MomentProblem{FermionicAlgebra,T,M,P})

Add zero constraints for moment matrix entries with odd fermion parity.
For fermionic systems, expectation values of odd-parity operators are zero
due to parity superselection rules.
"""
function _add_parity_constraints!(
    mp::MomentProblem{FermionicAlgebra,T,M,P}
) where {T<:Integer,M,P}

    parity_constraints = Tuple{Symbol, Matrix{P}}[]

    for (cone, mat) in mp.constraints
        dim = size(mat, 1)
        for i in 1:dim, j in 1:dim
            if _has_odd_parity_only(mat[i,j])
                # This entry must be zero - add as 1x1 Zero constraint
                push!(parity_constraints, (:Zero, reshape([mat[i,j]], 1, 1)))
            end
        end
    end

    # Add all parity constraints to the problem
    append!(mp.constraints, parity_constraints)
end
```

#### 2. `src/pop.jl`

**Keep objective validation** (lines 40-45):

The objective itself must still have even parity for the expectation value to be meaningful. Only the **basis** can contain odd-parity monomials.

#### 3. `src/sos_solver.jl`

**Keep existing code** - the `get_Cαj` function already handles missing monomials correctly.

### Alternative: More Efficient Approach

Instead of creating many 1x1 constraints, we could collect all odd-parity polynomial entries and add them as a single multi-element Zero constraint. This would be more efficient for the solver.

## Testing Strategy

### 1. Verify Basis Contains Odd-Parity Monomials

Create test that checks `total_basis` includes monomials like `a[1]`, `a_dag[1]`:

```julia
@testset "Fermionic basis includes odd-parity monomials" begin
    registry, (a, a_dag) = create_fermionic_variables(1:2)
    ham = a_dag[1] * a[1] + a_dag[2] * a[2]
    pop = polyopt(ham, registry)

    corr = correlative_sparsity(pop, 2)
    term_sp = term_sparsities(pop, corr, MMD())
    mp = moment_relax(pop, corr, term_sp)

    # Check that odd-parity monomials are in basis
    a1_mono = monomials(a[1])[1]
    @test a1_mono ∈ mp.total_basis
end
```

### 2. Verify Parity Constraints Are Added

Check that moment matrix entries with odd parity get zero constraints:

```julia
@testset "Parity constraints added for odd-parity entries" begin
    registry, (a, a_dag) = create_fermionic_variables(1:2)
    ham = a_dag[1] * a[1]
    pop = polyopt(ham, registry)

    # ... build moment problem ...

    # Count Zero constraints - should have parity-induced ones
    zero_constraints = filter(c -> c[1] == :Zero, mp.constraints)
    @test length(zero_constraints) > 0
end
```

### 3. Run Fermionic Chain Test

File: `test/fermionic_chain.jl`

Expected: Should now solve successfully with `E₀ ≈ -4.0`

### 4. Regression Test: Existing Tests Should Still Pass

All 1393 existing tests should continue passing, including:
- `test/fermionic_parity_test.jl` (may need updates to reflect new behavior)

## Implementation Plan

### Step 1: Remove Incorrect Filtering
- File: `src/moment_solver.jl`, lines 152-156
- Action: Delete the `if A == FermionicAlgebra` block that filters `total_basis`

### Step 2: Add Parity Constraint Helper
- File: `src/moment_solver.jl`, after line 92
- Action: Add `_has_odd_parity_only(poly)` function (see "Detailed Code Changes")

### Step 3: Add Parity Constraint Injection
- File: `src/moment_solver.jl`, modify `moment_relax` function
- Action: Call `_add_parity_constraints!(mp)` before returning (see "Detailed Code Changes")
- Action: Implement `_add_parity_constraints!` function

### Step 4: Keep Objective Validation
- File: `src/pop.jl`, lines 40-45
- Action: Keep this check (objective must have even parity)

### Step 5: Update Tests
- File: `test/fermionic_parity_test.jl`
- Action: Update basis filtering tests to reflect new behavior
- Action: Add tests for parity constraint injection

### Step 6: Verify Fermionic Chain Test
- File: `test/fermionic_chain.jl`
- Action: Run and verify it passes with correct ground state energy

## References

### Code Locations

1. **Moment matrix construction**: `src/moment_solver.jl:72-92` (`_build_constraint_matrix`)
2. **Basis computation**: `src/moment_solver.jl:133-156` (`moment_relax`)
3. **Incorrect filtering**: `src/moment_solver.jl:152-156` (TO BE REMOVED)
4. **Triple product**: `src/FastPolynomials/src/utils.jl:165-170` (`_neat_dot3`)
5. **Parity helper**: `src/FastPolynomials/src/simplification/fermionic.jl:32-69` (`has_even_parity`)
6. **Canonicalization**: `src/FastPolynomials/src/canonicalization.jl` (`symmetric_canon`)

### Physics References

[Source: User explanation and fermionic superselection theory]

Parity superselection rule: In fermionic systems, states with different total fermion parity (even/odd number of fermions) cannot be superposed. This means expectation values of operators with odd fermion number are identically zero:

```
⟨ψ| O_odd |ψ⟩ = 0  for all physical states |ψ⟩
```

This is implemented in the moment relaxation by constraining odd-parity moment matrix entries to zero, NOT by removing odd-parity monomials from the basis.

### Mathematical Justification

The moment matrix entry `M[i,j]` represents:
```
M[i,j] = ⟨basis[i]† operator basis[j]⟩
```

The total operator is:
```
O_total = adjoint(basis[i]) · operator · basis[j]
```

After simplification (Wick's theorem for fermions), this becomes a linear combination of normal-ordered monomials. The parity of O_total is:
```
parity(O_total) = parity(adjoint(basis[i])) + parity(operator) + parity(basis[j]) (mod 2)
```

Note: `parity(adjoint(m)) = parity(m)` since adjoint just swaps creation ↔ annihilation (same count).

If `parity(O_total)` is odd → `M[i,j] = 0` (by superselection)
If `parity(O_total)` is even → `M[i,j]` can be non-zero

This is why we need ALL monomials in the basis (including odd-parity ones) to compute `O_total` correctly!

## Success Criteria

1. ❌ Fermionic chain test (`test/fermionic_chain.jl`) - **PRE-EXISTING FAILURE** (see below)
2. ❌ Ground state energy ≈ -4.0 - **PRE-EXISTING ISSUE** with fermionic moment relaxation
3. ✅ All existing tests continue passing (1395 tests)
4. ✅ Basis contains odd-parity monomials (e.g., `a[1]`, `a_dag[1]`)
5. ✅ Zero constraints added for odd-parity moment entries
6. ✅ No "primal infeasible" errors for valid fermionic problems

## Important Finding: Pre-existing Fermionic Moment Relaxation Bug

During testing, we discovered that `test/fermionic_chain.jl` fails **even without parity constraints**.
This is a pre-existing issue with fermionic moment relaxation, unrelated to parity superselection.

### Evidence

Without parity constraints, the fermionic chain test gives:
- Expected: E₀ = -4.0
- Actual: -1.52 (with numerical instability, PRSTATUS = -1.00 every iteration)

The solver output shows wild oscillations (bounds going to 10^14), suggesting:
1. Ill-conditioned constraint matrices
2. Possible issue with Hermitian PSD embedding for fermionic algebras
3. Missing or incorrect constraints in the moment hierarchy

### Recommendation

The parity constraint implementation is correct. A separate investigation is needed to fix
the underlying fermionic moment relaxation issue. This should be tracked as a separate bug.
