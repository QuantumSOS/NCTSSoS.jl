# Type System Revamp: NCTSSoS.jl v2

## Goals

1. **Fix comparison bugs** - Type-enforced canonical invariants
2. **API clarity** - Predictable return types from `simplify` and `*`

---

## Key Decisions

| Decision | Choice |
|----------|--------|
| New types | `PauliMonomial{T}` + `PhysicsMonomial{A,T}` |
| `PauliMonomial` structure | `mono::Monomial{Pauli,T}` + `phase_k::UInt8` |
| `PhysicsMonomial` structure | `coeffs::Vector{Int}` + `monos::Vector{Monomial{A,T}}` |
| Constructor validation | Yes - `Monomial{A,T}` rejects invalid words |
| Multiplication returns | Algebra-dependent (see below) |
| Moment keys for Pauli | `PauliMonomial{T}` (includes `phase_k`) |
| JuMP moment variables | Dedup by underlying `Monomial{A,T}` |
| State polynomial support | Projector/Unipotent/NC only (drop Pauli/Fermi/Boson) |
| `Monomial.word` | Internal; copied on construction; treat immutable |
| Pauli phase storage | `phase_k ∈ 0:3` encoding `(im)^phase_k` |
| Type hierarchy | Flat `AbstractMonomial{A,T}` supertype |
| Migration | Clean break (v2) |

---

## Type Hierarchy

```
AbstractMonomial{A<:AlgebraType, T<:Integer}
├── Monomial{A,T}           # Bare word, constructor validates
├── PauliMonomial{T}        # Canonical mono + phase_k ∈ 0:3
└── PhysicsMonomial{A,T}    # Sum of int×mono (A ∈ {Fermi,Boson})
```

---

## New Types

### 1. `PauliMonomial{T}` (NEW)

```julia
struct PauliMonomial{T<:Integer} <: AbstractMonomial{PauliAlgebra,T}
    mono::Monomial{PauliAlgebra,T}  # canonical form
    phase_k::UInt8                  # encodes phase = (im)^phase_k, phase_k ∈ 0:3
end
```

**Invariants:**
- `mono` has ≤1 operator per site, sites sorted
- `phase_k` encodes algebraic phase from σ² = I, cyclic products

**Constructor:**
```julia
function PauliMonomial(word::Vector{T}) where T
    canonical_word, phase_k = _simplify_pauli_word!(copy(word))
    mono = Monomial{PauliAlgebra,T}(canonical_word)
    PauliMonomial{T}(mono, phase_k)
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
| `FermionicAlgebra` | normal-ordered; creators sorted asc; annihilators sorted desc |
| `BosonicAlgebra` | normal-ordered; creators sorted asc; annihilators sorted desc |
| `ProjectorAlgebra` | no P² terms, sites sorted |
| `UnipotentAlgebra` | no U² terms, sites sorted |
| `NonCommutativeAlgebra` | sites sorted |

```julia
function Monomial{A,T}(word::Vector{T}) where {A,T}
    _validate_word!(A, word)  # throws if invalid
    new{A,T}(copy(word))      # internal storage; do not mutate
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

### Phase 1: Internal Helpers ✅
- [x] `_validate_pauli_word!(word)` - check ≤1 op/site, sorted; throw if invalid
- [x] `_validate_fermionic_word!(word)` - check normal-ordered; throw if invalid
- [x] `_validate_bosonic_word!(word)` - check normal-ordered; throw if invalid
- [x] `_validate_projector_word!(word)` - check no P², sorted; throw if invalid
- [x] `_validate_unipotent_word!(word)` - check no U², sorted; throw if invalid
- [x] `_validate_nc_word!(word)` - check sorted; throw if invalid
- [x] Update `_simplify_pauli_word!` to return `(canonical_word, phase_k::UInt8)`
- [x] Update `_simplify_fermionic_word!` to return `Vector{Tuple{Int,Vector{T}}}`
- [x] Update `_simplify_bosonic_word!` to return `Vector{Tuple{Int,Vector{T}}}`

### Phase 2: Add `AbstractMonomial` hierarchy ✅
- [x] Define `AbstractMonomial{A,T}` abstract type
- [x] Make `Monomial{A,T} <: AbstractMonomial{A,T}`
- [x] Implement shared interface: `degree`, `variable_indices`, `isless`

### Phase 3: Add `PauliMonomial{T}` ✅
- [x] Create `src/types/pauli_monomial.jl`
- [x] Define struct with `mono::Monomial{Pauli,T}` + `phase_k::UInt8`
- [x] Implement constructor calling `_simplify_pauli_word!`
- [x] Implement `==`, `hash`, `isless` (include phase_k)
- [x] Implement `*` returning `PauliMonomial`
- [x] Implement `adjoint` (conjugate phase_k)
- [x] Implement `symmetric_canon` (identity on mono, conjugate phase_k)
- [x] Add display methods
- [x] Export from module

### Phase 4: Add `PhysicsMonomial{A,T}` ✅
- [x] Create `src/types/physics_monomial.jl`
- [x] Define struct with `coeffs::Vector{Int}` + `monos::Vector{Monomial{A,T}}`
- [x] Implement constructor calling `_simplify_*_word!`
- [x] Implement `==`, `hash`, `isless`
- [x] Implement `*` (distribute, normal-order, collect)
- [x] Implement `adjoint`
- [x] Add display methods
- [x] Export from module

### Phase 5: Update `Monomial{A,T}` constructor (PARTIAL)
- [x] Add `validate()` function for explicit validation
- [ ] Add validation dispatch in inner constructor (DEFERRED - breaks internal code)
- [ ] Update tests that create invalid monomials

### Phase 6: Update multiplication dispatch (BREAKING) ✅
- [x] `Monomial{Pauli} * Monomial{Pauli}` → `PauliMonomial`
- [x] `Monomial{Fermi} * Monomial{Fermi}` → `PhysicsMonomial`
- [x] `Monomial{Boson} * Monomial{Boson}` → `PhysicsMonomial`
- [ ] Update call sites expecting old return types (DEFERRED to Phase 7-8)

### Phase 7: Update basis construction
- [ ] `get_ncbasis` returns `Vector{PauliMonomial{T}}` for Pauli
- [ ] `get_ncbasis` returns `Vector{PhysicsMonomial{A,T}}` for Fermi/Boson
- [ ] Update `CorrelativeSparsity` type parameter M
- [ ] Update `TermSparsity` type parameter M

### Phase 8: Update moment matrix integration
- [ ] `MomentProblem` key type uses wrapper types
- [ ] Dedup JuMP moment vars by underlying `Monomial{A,T}` (wrapper keys map to scalar/sum of base moments)
- [ ] Update SOS dualization for new key types
- [ ] Verify Hermitian embedding still works

### Phase 9: Restrict state polynomial support
- [ ] Add type parameter constraint: `A <: Union{ProjectorAlgebra, UnipotentAlgebra, NonCommutativeAlgebra}`
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

### Why separate `mono` + `phase_k` in `PauliMonomial`?
- Phase is algebraic only {1,-1,i,-i} from simplification
- Scalar coefficients go to `Polynomial`, not absorbed into phase_k
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
