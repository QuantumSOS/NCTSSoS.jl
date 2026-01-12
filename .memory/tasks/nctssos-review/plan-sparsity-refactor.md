# Refactoring Plan: sparsity.jl DRY improvements

**Status**: BLOCKED - Waiting for Monomial hierarchy revamp

## Context
During sparsity.jl review, identified two DRY improvements:
- MAINT-002: `simplified_monomials()` helper to avoid Polynomial wrapping
- MAINT-003: Unify Polynomial/State sparsity functions (~400 lines duplication)

## Prerequisite
Monomial hierarchy revamp must be completed first. The `simplified_monomials()` helper depends on a clean abstraction over:
- `Monomial{A,T}` (SimpleAlgebra)
- `PauliMonomial{T}`
- `PhysicsMonomial{A,T}` (Fermionic/Bosonic)

## Approach
Implement MAINT-002 first (smaller, enables cleaner MAINT-003), then MAINT-003.

---

## Step 1: Add `simplified_monomials()` helper

**File**: `src/simplification/interface.jl` (new file)

**Implementation**:
```julia
"""
    simplified_monomials(m::Monomial) -> Vector{Monomial}

Extract monomials from simplified result, discarding coefficients/phases.
For sparsity computation where only monomial structure matters.
"""
# SimpleAlgebra: simplify returns Monomial directly
simplified_monomials(m::Monomial{A,T}) where {A<:SimpleAlgebra,T} = [simplify(m)]

# PauliAlgebra: extract mono from PauliMonomial
simplified_monomials(m::Monomial{PauliAlgebra,T}) where {T} = [simplify(m).mono]

# Fermionic/Bosonic: extract monomials from PhysicsMonomial
function simplified_monomials(m::Monomial{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T}
    pm = simplify(m)
    iszero(pm) && return Monomial{A,T}[]
    pm.monos
end
```

**Update sparsity.jl** (~6 replacements):
- Line 469: `append!(diagonal_entries, monomials(Polynomial(simplify(diag_mono))))` → `append!(diagonal_entries, simplified_monomials(diag_mono))`
- Line 534-535: `Polynomial(simplify(_neat_dot3(...)))` → use `simplified_monomials` pattern
- Line 592: Similar in `term_sparsity_graph_supp`

---

## Step 2: Unify correlative sparsity structs

**Goal**: Single `CorrelativeSparsity{A,T,P,M}` struct for both Polynomial and NCStatePolynomial

**Changes**:
1. Remove `StateCorrelativeSparsity` struct (line 405-413)
2. Update `CorrelativeSparsity` type params to work for both
3. Update `Base.show` dispatch

---

## Step 3: Unify graph construction

**Create shared helper**:
```julia
function _monomial_var_indices(m::Monomial{A,T}) where {A,T}
    unique([T(abs(idx)) for idx in m.word])
end

function _monomial_var_indices(ncsw::NCStateWord)
    variables(ncsw)
end
```

**Merge** `get_correlative_graph` and `get_state_correlative_graph` into single function.

---

## Step 4: Unify constraint assignment

**Merge** `assign_constraint` and `assign_state_constraint` - already identical after BUG-001 fix.

---

## Step 5: Unify main correlative_sparsity function

**Create dispatch helpers**:
```julia
_get_basis(reg, order, ::Type{<:Polynomial{A,T}}) where {A,T} = get_ncbasis(reg, order)
_get_basis(reg, order, ::Type{<:NCStatePolynomial{C,ST,A,T}}) where {C,ST,A,T} = get_state_basis(reg, order; state_type=ST)
```

---

## Step 6: Unify term sparsity functions

Similar approach for:
- `init_activated_supp`
- `term_sparsities`
- `get_term_sparsity_graph`
- `iterate_term_sparse_supp`
- `term_sparsity_graph_supp`

---

## Verification

1. `make test-minimal` - smoke test
2. `make test-file FILE=test/relaxations/sparsity.jl` - sparsity tests
3. `make test-polynomials` - algebra tests
4. `make test` (if Mosek available) - full solver tests

---

## Estimated Impact

- Before: ~985 lines
- After: ~550-600 lines
- Reduction: ~35-40%

## Critical Files
- `src/simplification/interface.jl` (new)
- `src/optimization/sparsity.jl` (major refactor)
- `src/NCTSSoS.jl` (add include/export)
