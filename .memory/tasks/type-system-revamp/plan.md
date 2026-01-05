# Type System Revamp: NCTSSoS.jl v2

## Goals

1. **Fix comparison bugs** - Type-enforced canonical invariants
2. **API clarity** - Predictable return types from `simplify` and `*`

---

## Key Decisions

| Decision | Choice |
|----------|--------|
| New types | `PauliMonomial{T}` + `PhysicsMonomial{A,T}` |
| `PauliMonomial` structure | `mono::Monomial{Pauli,T}` + `phase::ComplexF64` |
| `PhysicsMonomial` structure | `coeffs::Vector{Int}` + `monos::Vector{Monomial{A,T}}` |
| Constructor validation | Yes - `Monomial{A,T}` rejects invalid words |
| Multiplication returns | Algebra-dependent (see below) |
| Moment keys for Pauli | `PauliMonomial{T}` (includes phase) |
| State polynomial support | Projector/Unipotent/NC only (drop Pauli/Fermi/Boson) |
| Type hierarchy | Flat `AbstractMonomial{A,T}` supertype |
| Migration | Clean break (v2) |

---

## Type Hierarchy

```
AbstractMonomial{A<:AlgebraType, T<:Integer}
├── Monomial{A,T}           # Bare word, constructor validates
├── PauliMonomial{T}        # Canonical mono + phase ∈ {1,-1,i,-i}
└── PhysicsMonomial{A,T}    # Sum of int×mono (A ∈ {Fermi,Boson})
```

---

## New Types

### 1. `PauliMonomial{T}` (NEW)

```julia
struct PauliMonomial{T<:Integer} <: AbstractMonomial{PauliAlgebra,T}
    mono::Monomial{PauliAlgebra,T}  # canonical form
    phase::ComplexF64               # algebraic phase only: {1, -1, i, -i}
end
```

**Invariants:**
- `mono` has ≤1 operator per site, sites sorted
- `phase` is algebraic (from σ² = I, cyclic products)

**Constructor:**
```julia
function PauliMonomial(word::Vector{T}) where T
    canonical_word, phase = _simplify_pauli_word!(copy(word))
    mono = Monomial{PauliAlgebra,T}(canonical_word)
    PauliMonomial{T}(mono, phase)
end
```

### 2. `PhysicsMonomial{A,T}` (NEW)

```julia
struct PhysicsMonomial{A<:Union{FermionicAlgebra,BosonicAlgebra}, T<:Integer} <: AbstractMonomial{A,T}
    coeffs::Vector{Int}           # integer coefficients
    monos::Vector{Monomial{A,T}}  # normal-ordered monomials
end
```

**Invariants:**
- Each `monos[i]` is normal-ordered (validated by Monomial constructor)
- Represents sum: Σ coeffs[i] × monos[i]
- `coeffs[i]` are integers from (anti)commutation relations

**Constructor:**
```julia
function PhysicsMonomial{A}(word::Vector{T}) where {A,T}
    terms = _simplify_fermionic_word!(copy(word))  # or bosonic
    coeffs = [t[1] for t in terms]
    monos = [Monomial{A,T}(t[2]) for t in terms]
    PhysicsMonomial{A,T}(coeffs, monos)
end
```

### 3. `Monomial{A,T}` (MODIFIED)

Constructor validates canonical form per algebra:

| Algebra | Validation |
|---------|------------|
| `PauliAlgebra` | ≤1 op/site, sites sorted |
| `FermionicAlgebra` | normal-ordered |
| `BosonicAlgebra` | normal-ordered |
| `ProjectorAlgebra` | no P² terms, sites sorted |
| `UnipotentAlgebra` | no U² terms, sites sorted |
| `NonCommutativeAlgebra` | sites sorted |

```julia
function Monomial{A,T}(word::Vector{T}) where {A,T}
    _validate_word!(A, word)  # throws if invalid
    new{A,T}(word)
end
```

---

## Operation Return Types

### Multiplication

| Operation | Returns |
|-----------|---------|
| `Monomial{Pauli} × Monomial{Pauli}` | `PauliMonomial{T}` |
| `Monomial{Fermi} × Monomial{Fermi}` | `PhysicsMonomial{Fermi,T}` |
| `Monomial{Boson} × Monomial{Boson}` | `PhysicsMonomial{Boson,T}` |
| `Monomial{Proj/Unip/NC} × Monomial{...}` | `Monomial{A,T}` |
| `PauliMonomial × PauliMonomial` | `PauliMonomial{T}` (closed) |
| `PhysicsMonomial × PhysicsMonomial` | `PhysicsMonomial{A,T}` (closed) |
| `PauliMonomial × scalar` | `Polynomial{Pauli,T,ComplexF64}` |

### Simplify

| Input | Returns |
|-------|---------|
| `Monomial{A,T}` | `Monomial{A,T}` (identity - already canonical) |
| `PauliMonomial{T}` | `PauliMonomial{T}` (identity) |
| `PhysicsMonomial{A,T}` | `PhysicsMonomial{A,T}` (identity) |

---

## Implementation Steps

### Phase 1: Internal Helpers
- [ ] `_validate_pauli_word!(word)` - check ≤1 op/site, sorted; throw if invalid
- [ ] `_validate_fermionic_word!(word)` - check normal-ordered; throw if invalid
- [ ] `_validate_bosonic_word!(word)` - check normal-ordered; throw if invalid
- [ ] `_validate_projector_word!(word)` - check no P², sorted; throw if invalid
- [ ] `_validate_unipotent_word!(word)` - check no U², sorted; throw if invalid
- [ ] `_validate_nc_word!(word)` - check sorted; throw if invalid
- [ ] Update `_simplify_pauli_word!` to return `(canonical_word, phase)`
- [ ] Update `_simplify_fermionic_word!` to return `Vector{Tuple{Int,Vector{T}}}`
- [ ] Update `_simplify_bosonic_word!` to return `Vector{Tuple{Int,Vector{T}}}`

### Phase 2: Add `AbstractMonomial` hierarchy
- [ ] Define `AbstractMonomial{A,T}` abstract type
- [ ] Make `Monomial{A,T} <: AbstractMonomial{A,T}`
- [ ] Implement shared interface: `degree`, `variable_indices`, `isless`

### Phase 3: Add `PauliMonomial{T}`
- [ ] Create `src/types/pauli_monomial.jl`
- [ ] Define struct with `mono::Monomial{Pauli,T}` + `phase::ComplexF64`
- [ ] Implement constructor calling `_simplify_pauli_word!`
- [ ] Implement `==`, `hash`, `isless` (include phase)
- [ ] Implement `*` returning `PauliMonomial`
- [ ] Implement `adjoint` (conjugate phase)
- [ ] Implement `symmetric_canon` (identity on mono, conjugate phase)
- [ ] Add display methods
- [ ] Export from module

### Phase 4: Add `PhysicsMonomial{A,T}`
- [ ] Create `src/types/physics_monomial.jl`
- [ ] Define struct with `coeffs::Vector{Int}` + `monos::Vector{Monomial{A,T}}`
- [ ] Implement constructor calling `_simplify_*_word!`
- [ ] Implement `==`, `hash`, `isless`
- [ ] Implement `*` (distribute, normal-order, collect)
- [ ] Implement `adjoint`
- [ ] Add display methods
- [ ] Export from module

### Phase 5: Update `Monomial{A,T}` constructor (BREAKING)
- [ ] Add validation dispatch in inner constructor
- [ ] Create `Monomial_unchecked` for internal use
- [ ] Update tests that create invalid monomials

### Phase 6: Update multiplication dispatch (BREAKING)
- [ ] `Monomial{Pauli} * Monomial{Pauli}` → `PauliMonomial`
- [ ] `Monomial{Fermi} * Monomial{Fermi}` → `PhysicsMonomial`
- [ ] `Monomial{Boson} * Monomial{Boson}` → `PhysicsMonomial`
- [ ] Update call sites expecting old return types

### Phase 7: Update basis construction
- [ ] `get_ncbasis` returns `Vector{PauliMonomial{T}}` for Pauli
- [ ] `get_ncbasis` returns `Vector{PhysicsMonomial{A,T}}` for Fermi/Boson
- [ ] Update `CorrelativeSparsity` type parameter M
- [ ] Update `TermSparsity` type parameter M

### Phase 8: Update moment matrix integration
- [ ] `MomentProblem` key type uses wrapper types
- [ ] Update SOS dualization for new key types
- [ ] Verify Hermitian embedding still works

### Phase 9: Restrict state polynomial support
- [ ] Add compile-time check: state polys only for Proj/Unip/NC
- [ ] Remove/deprecate state poly code for Pauli/Fermi/Boson
- [ ] Update tests

### Phase 10: Type conversion
- [ ] `Polynomial(::PauliMonomial)` → `Polynomial{Pauli,T,ComplexF64}`
- [ ] `Polynomial(::PhysicsMonomial)` → `Polynomial{A,T,ComplexF64}` (promote Int)

---

## Critical Files

| File | Changes |
|------|---------|
| `src/types/monomial.jl` | Add `AbstractMonomial`, validation in constructor |
| `src/types/pauli_monomial.jl` | **NEW** - `PauliMonomial{T}` |
| `src/types/physics_monomial.jl` | **NEW** - `PhysicsMonomial{A,T}` |
| `src/simplification/*.jl` | Update return types for `_simplify_*` |
| `src/algorithms/basis.jl` | Return wrapper types from `get_ncbasis` |
| `src/optimization/sparsity.jl` | Update type parameters |
| `src/optimization/moment.jl` | Update key types |
| `src/states/*.jl` | Restrict to simple algebras |
| `test/polynomials/*.jl` | Update for new types |

---

## Design Rationale

### Why separate `mono` + `phase` in `PauliMonomial`?
- Phase is algebraic only {1,-1,i,-i} from simplification
- Scalar coefficients go to `Polynomial`, not absorbed into phase
- Clear separation: algebra structure vs. user coefficients

### Why `PhysicsMonomial` stores sum of monomials?
- Normal-ordering expands: `a₁a₂†a₁† = -a₁†a₂† + a₂†`
- Cannot represent as single monomial
- Integer coefficients preserve exactness from (anti)commutation

### Why validate in `Monomial` constructor?
- Type guarantees invariants - no runtime checks downstream
- Invalid states become unrepresentable
- Follows `StateSymbol` pattern which works well

### Why drop state poly support for Pauli/Fermi/Boson?
- Pauli: complex phases complicate real expectation values
- Fermi/Boson: expansion under normal ordering doesn't fit state poly model
- Simplifies codebase, focuses on working use cases

---

## Migration Notes

- **Version**: v2.0.0 (breaking)
- Old code creating raw Pauli monomials will throw
- Old code expecting `Term` from Pauli operations needs update
- Old code using state polys with Pauli/Fermi/Boson needs migration
- Provide clear error messages guiding users to new patterns
