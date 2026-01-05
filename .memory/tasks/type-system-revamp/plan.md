# Type System Revamp: NCTSSoS.jl v2

## Goals

1. **Fix comparison bugs** - Type-enforced canonical invariants
2. **API clarity** - Predictable return types from `simplify` and `*`

---

## Key Decisions

| Decision | Choice |
|----------|--------|
| New types | Only `PauliMonomial{T}` (wraps `Term`) |
| Constructor validation | Yes - reject invalid monomials |
| Multiplication returns | Algebra-dependent |
| Keep `Term{M,C}` | Yes - `PauliMonomial` wraps it |
| Migration | Clean break (v2) |
| Moment keys for Pauli | `PauliMonomial` (includes phase) |

---

## Type Changes

### 1. `PauliMonomial{T}` (NEW)

```julia
struct PauliMonomial{T<:Integer} <: AbstractMonomial
    term::Term{Monomial{PauliAlgebra,T}, ComplexF64}
end
```

**Invariants:**
- Underlying monomial has ≤1 operator per site
- Sites sorted ascending
- Phase accumulated from cyclic products

**Constructor:**
```julia
function PauliMonomial(word::Vector{T}) where T
    simplified_term = _simplify_pauli(word)  # Returns Term
    _validate_pauli_canonical!(simplified_term.monomial.word)
    PauliMonomial{T}(simplified_term)
end
```

### 2. `Monomial{A,T}` Constructor Validation

Modify inner constructor to enforce algebra-specific invariants:

| Algebra | Validation |
|---------|------------|
| `PauliAlgebra` | **Reject** if >1 op/site |
| `FermionicAlgebra` | **Reject** if not normal-ordered |
| `BosonicAlgebra` | **Reject** if not normal-ordered |
| `NonCommutativeAlgebra` | **Apply** site-sort |
| `ProjectorAlgebra` | **Apply** P²=P + site-sort |
| `UnipotentAlgebra` | **Apply** U²=I + site-sort |

### 3. Multiplication Return Types

| Algebra | `m1 * m2` returns |
|---------|-------------------|
| NC/Projector/Unipotent | `Monomial{A,T}` |
| Pauli | `PauliMonomial{T}` |
| Fermionic/Bosonic | `Polynomial{A,T,C}` |

### 4. `simplify` Return Types

| Input | Returns |
|-------|---------|
| `Monomial{NC/Proj/Unip}` | `Monomial` (identity - already canonical) |
| `Monomial{Pauli}` | `PauliMonomial` |
| `Monomial{Fermi/Boson}` | `Polynomial` |
| `PauliMonomial` | `PauliMonomial` (identity) |

---

## Implementation Steps

### Phase 1: Add `PauliMonomial` (non-breaking)
- [ ] Create `src/types/pauli_monomial.jl`
- [ ] Define struct wrapping `Term{Monomial{PauliAlgebra,T}, ComplexF64}`
- [ ] Implement constructors with validation
- [ ] Implement `==`, `hash`, `isless`, `degree`
- [ ] Implement iteration protocol `(coef, mono)` for unified processing
- [ ] Implement `*` returning `PauliMonomial`
- [ ] Add display methods
- [ ] Export from module

### Phase 2: Add validation helpers
- [ ] `_validate_pauli_canonical!(word)` - throws if >1 op/site
- [ ] `is_normal_ordered(word)` for Fermi/Boson
- [ ] `_validate_fermionic_canonical!(word)`
- [ ] `_validate_bosonic_canonical!(word)`

### Phase 3: Update `Monomial` constructor (breaking)
- [ ] Add `_validate_monomial_word(::Type{A}, word)` dispatch
- [ ] Modify inner constructor to call validation
- [ ] Add `Monomial_unchecked` for internal use if needed
- [ ] Update tests that create invalid monomials

### Phase 4: Update `simplify` returns (breaking)
- [ ] `simplify(::Monomial{PauliAlgebra})` → `PauliMonomial`
- [ ] `simplify(::Monomial{NC/Proj/Unip})` → identity (already canonical)
- [ ] Update `_collect_simplified_terms!` in polynomial.jl

### Phase 5: Update multiplication dispatch (breaking)
- [ ] `Monomial{Pauli} * Monomial{Pauli}` → `PauliMonomial`
- [ ] Keep Fermi/Boson returning `Polynomial`
- [ ] Update call sites expecting old return types

### Phase 6: Moment matrix integration
- [ ] Update `MomentProblem` key type for Pauli → `PauliMonomial`
- [ ] Update basis construction
- [ ] Update SOS dualization

---

## Critical Files

| File | Changes |
|------|---------|
| `src/types/pauli_monomial.jl` | **NEW** - PauliMonomial definition |
| `src/types/monomial.jl` | Add validation in constructor |
| `src/simplification/pauli.jl` | Return `PauliMonomial` |
| `src/types/polynomial.jl` | Update `_collect_simplified_terms!` |
| `src/optimization/moment.jl` | Key type for Pauli |
| `test/polynomials/simplify.jl` | Update Pauli tests |
| `test/polynomials/monomials.jl` | Test constructor validation |

---

## Design Rationale

### Why wrap `Term` in `PauliMonomial`?
- Preserves `Term` allocation semantics for benchmarking
- Type safety via wrapper without duplicating storage logic
- Can measure overhead of wrapper vs. bare Term

### Why validate in constructor?
- Eliminates class of comparison bugs (non-canonical == checks)
- Type guarantees invariants - no runtime checks needed downstream
- Follows `StateSymbol` pattern which already works well

### Why algebra-dependent `*` return types?
- Reflects algebraic reality: Pauli closed (single term), Fermi/Boson expand
- Type system encodes mathematical properties
- Clearer API than polymorphic `simplify` returns

---

## Migration Notes

- **Version**: v2.0.0 (breaking)
- Old code creating raw Pauli monomials will throw
- Old code expecting `Term` from Pauli simplify needs update
- Provide clear error messages guiding users to new patterns

---

## Open Questions

_Add clarifications here_
