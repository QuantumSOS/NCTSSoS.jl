# Function Type Specifications for NCTSSoS.jl v2

This document records the input/output type specifications for all functions in the polynomial optimization pipeline, to inform the type system revamp.

---

## 1. Variable Registry Creation

### `create_pauli_variables`
- **File:** `src/types/registry.jl`
- **Input:** `subscripts::AbstractVector`
- **Output:** `(VariableRegistry{PauliAlgebra,T}, (Vector{PauliMonomial{T}}, Vector{PauliMonomial{T}}, Vector{PauliMonomial{T}}))`
- **Notes:** Returns `PauliMonomial{T}` wrapper. Phase stored as `phase_k::UInt8` (encoding `(im)^phase_k`); converted to `ComplexF64` only when forming `Polynomial`/JuMP expressions.

### `create_fermionic_variables`
- **File:** `src/types/registry.jl`
- **Input:** `subscripts::AbstractVector`
- **Output:** `(VariableRegistry{FermionicAlgebra,T}, (Vector{PhysicsMonomial{FermionicAlgebra,T}},))`
- **Notes:** Returns `PhysicsMonomial` - sum of normal-ordered terms with `Int64` coefficients.

### `create_bosonic_variables`
- **File:** `src/types/registry.jl`
- **Input:** `subscripts::AbstractVector`
- **Output:** `(VariableRegistry{BosonicAlgebra,T}, (Vector{PhysicsMonomial{BosonicAlgebra,T}},))`
- **Notes:** Returns `PhysicsMonomial` - sum of normal-ordered terms with `Int64` coefficients (from [a,a†]=1).

### `create_projector_variables`
- **File:** `src/types/registry.jl`
- **Input:** `prefix_subscripts::Vector{Tuple{String, AbstractVector}}`
- **Output:** `(VariableRegistry{ProjectorAlgebra,T}, (Vector{Monomial{ProjectorAlgebra,T}},...))`
- **Notes:** Raw `Monomial` - no wrapper needed (simple algebra).

### `create_unipotent_variables`
- **File:** `src/types/registry.jl`
- **Input:** `prefix_subscripts::Vector{Tuple{String, AbstractVector}}`
- **Output:** `(VariableRegistry{UnipotentAlgebra,T}, (Vector{Monomial{UnipotentAlgebra,T}},...))`
- **Notes:** Raw `Monomial` - no wrapper needed (simple algebra).

### `create_noncommutative_variables`
- **File:** `src/types/registry.jl`
- **Input:** `prefix_subscripts::Vector{Tuple{String, AbstractVector}}`
- **Output:** `(VariableRegistry{NonCommutativeAlgebra,T}, (Vector{Monomial{NonCommutativeAlgebra,T}},...))`
- **Notes:** Raw `Monomial` - no wrapper needed (generic NC).

---

## 2. Core Types

### Type Hierarchy
```
AbstractMonomial{A<:AlgebraType, T<:Integer}
├── Monomial{A,T}           # Bare word (no coefficient)
├── PauliMonomial{T}        # = AbstractMonomial{PauliAlgebra,T} + phase
└── PhysicsMonomial{A,T}    # = AbstractMonomial{A,T} + sum of int×word (A ∈ {Fermi,Boson})
```

### `Monomial{A,T}` (MODIFIED - constructor validates)
- **File:** `src/types/monomial.jl`
- **Type parameters:** `A<:AlgebraType`, `T<:Integer`
- **Fields:** `word::Vector{T}` (internal; do not mutate)
- **Invariants:** **Enforced by constructor** (algebra-dependent):
  - `PauliAlgebra`: ≤1 op/site, sites sorted, no σ² terms
  - `FermionicAlgebra`: normal-ordered; creators sorted asc; annihilators sorted desc (adjoint-closed)
  - `BosonicAlgebra`: normal-ordered; creators sorted asc; annihilators sorted desc
  - `ProjectorAlgebra`: no P² terms, sites sorted
  - `UnipotentAlgebra`: no U² terms, sites sorted
  - `NonCommutativeAlgebra`: sites sorted (commutative per site)
- **Constructor:** Throws error if invariant violated; stores `copy(word)`; word treated immutable by API.

### `Term{M,C}` (existing - unchanged)
- **File:** `src/types/term.jl`
- **Type parameters:** `M<:Monomial`, `C<:Number`
- **Fields:** `coef::C`, `mono::M`
- **Invariants:** None enforced

### `Polynomial{A,T,C}` (existing - unchanged)
- **File:** `src/types/polynomial.jl`
- **Type parameters:** `A<:AlgebraType`, `T<:Integer`, `C<:Number`
- **Fields:** `terms::Vector{Term{Monomial{A,T},C}}`
- **Invariants:** Terms sorted, like terms combined

### `PauliMonomial{T}` (NEW)
- **File:** `src/types/pauli_monomial.jl` (new file)
- **Type parameters:** `T<:Integer`
- **Fields:**
  - `mono::Monomial{PauliAlgebra,T}` - canonical form (≤1 op/site, sorted)
  - `phase_k::UInt8` - accumulated algebraic phase, encoding `(im)^phase_k` with `phase_k ∈ 0:3`
- **Invariants:** Always in canonical form (enforced by constructor)
- **Constructor:** `PauliMonomial(word::Vector{T})` - auto-canonicalizes raw word

### `PhysicsMonomial{A,T}` (NEW)
- **File:** `src/types/physics_monomial.jl` (new file)
- **Type parameters:** `A<:Union{FermionicAlgebra,BosonicAlgebra}`, `T<:Integer`
- **Fields:**
  - `coeffs::Vector{Int}` - integer coefficients (from commutation relations)
  - `monos::Vector{Monomial{A,T}}` - normal-ordered monomials (NOT raw indices!)
- **Invariants:** Each `monos[i]::Monomial{A,T}` is normal-ordered; represents sum Σ coeffs[i] * monos[i]; result is never zero (no full cancellation expected in valid use)
- **Constructor:** `PhysicsMonomial(word::Vector{T})` - normal-orders and collects terms into `Monomial{A,T}` objects

---

## 3. Monomial Operations

### `Base.:(*)` (Monomial × Monomial)
- **Input:** `m1::Monomial{A,T}`, `m2::Monomial{A,T}` (both canonical)
- **Output:** Algebra-dependent:
  - `PauliAlgebra`: `PauliMonomial{T}` (may produce phase, needs re-canonicalization)
  - `FermionicAlgebra`: `PhysicsMonomial{FermionicAlgebra,T}` (expands via anticommutation)
  - `BosonicAlgebra`: `PhysicsMonomial{BosonicAlgebra,T}` (expands via commutation)
  - `Projector/Unipotent/NC`: `Monomial{A,T}` (stays canonical)
- **Notes:** Multiplication of canonical monomials may produce non-canonical result.

### `Base.:(*)` (PauliMonomial × PauliMonomial)
- **Input:** `pm1::PauliMonomial{T}`, `pm2::PauliMonomial{T}`
- **Output:** `PauliMonomial{T}` (closed under multiplication)
- **Notes:** Phases multiply, monomials multiply and re-canonicalize.

### `Base.:(*)` (PhysicsMonomial × PhysicsMonomial)
- **Input:** `pm1::PhysicsMonomial{A,T}`, `pm2::PhysicsMonomial{A,T}`
- **Output:** `PhysicsMonomial{A,T}` (closed under multiplication)
- **Notes:** Distributes over sums, normal-orders each product.

### `Base.:(*)` (PauliMonomial × scalar)
- **Input:** `pm::PauliMonomial{T}`, `c::Number`
- **Output:** `Polynomial{PauliAlgebra,T,ComplexF64}`
- **Notes:** Does NOT absorb into phase; promotes to polynomial.

### `simplify(::Monomial{A})`
- **Input:** `m::Monomial{A,T}` (already canonical by construction)
- **Output:** `m` (identity - already simplified)
- **Notes:** Since Monomial constructor enforces invariants, simplify is identity.

### `Base.adjoint(::Monomial)`
- **Input:** `m::Monomial{A,T}`
- **Output:** `Monomial{A,T}` (reversed word, sign-adjusted for signed algebras)
- **Notes:** For `FermionicAlgebra`, canonical ordering is chosen to make `adjoint(::Monomial)` return canonical `Monomial` without extra sign/term.

### `Base.adjoint(::PauliMonomial)`
- **Input:** `pm::PauliMonomial{T}`
- **Output:** `PauliMonomial{T}`
- **Notes:** Adjoint of mono, conjugate of phase.

---

## 4. Polynomial Operations

### `Polynomial{A,T,C}` construction
- **Default C:** `ComplexF64` for `PauliAlgebra`
- **Storage:** `Vector{Term{Monomial{A,T},C}}` (unchanged)
- **Invariants:** Terms sorted by monomial, like terms combined

### `Base.:(*)` (Polynomial × Polynomial)
- **Input:** `p1::Polynomial{A,T,C}`, `p2::Polynomial{A,T,C}`
- **Output:** `Polynomial{A,T,C}`
- **Notes:** Distributes, multiplies monomials (may need re-canonicalization via wrapper types), collects like terms.

### `Base.:(+)` (Polynomial + PauliMonomial)
- **Input:** `p::Polynomial{PauliAlgebra,T,C}`, `pm::PauliMonomial{T}`
- **Output:** `Polynomial{PauliAlgebra,T,C}`
- **Notes:** Converts `PauliMonomial` to `Term` (`(im)^phase_k` becomes coefficient) and adds.

### `Base.:(+)` (Polynomial + PhysicsMonomial)
- **Input:** `p::Polynomial{A,T,C}`, `pm::PhysicsMonomial{A,T}`
- **Output:** `Polynomial{A,T,C}`
- **Notes:** Converts each term in PhysicsMonomial sum to polynomial Term.

### `simplify(::Polynomial)`
- **Input:** `p::Polynomial{A,T,C}`
- **Output:** `Polynomial{A,T,C}`
- **Notes:** Since constituent monomials are canonical, simplify just combines like terms (identity if already combined).

---

## 5. Canonicalization

### `symmetric_canon(::Monomial{A,T})`
- **Input:** `m::Monomial{A,T}` (canonical)
- **Output:** `Monomial{A,T}`
- **Notes:** For `PauliAlgebra`: identity (single Pauli per site, self-adjoint). For others: min(word, reverse(word)).

### `symmetric_canon(::PauliMonomial{T})`
- **Input:** `pm::PauliMonomial{T}`
- **Output:** `PauliMonomial{T}`
- **Notes:** Same mono (identity on canonical mono), conjugate phase.

### `cyclic_canon(::Monomial{A,T})`
- **Input:** `m::Monomial{A,T}`
- **Output:** `Monomial{A,T}`
- **Notes:** Lexicographically smallest cyclic rotation.

---

## 6. Basis Construction

### `get_ncbasis`
- **Input:** `registry::VariableRegistry{A,T}`, `d::Int`
- **Output:** Algebra-dependent:
  - `PauliAlgebra`: `Vector{PauliMonomial{T}}`
  - `FermionicAlgebra`: `Vector{PhysicsMonomial{FermionicAlgebra,T}}`
  - `BosonicAlgebra`: `Vector{PhysicsMonomial{BosonicAlgebra,T}}`
  - `Projector/Unipotent/NC`: `Vector{Monomial{A,T}}`
- **Notes:** Returns all canonical basis elements up to degree d.

### `get_ncbasis_deg`
- **Input:** `registry::VariableRegistry{A,T}`, `d::Int`
- **Output:** Same as `get_ncbasis` but exactly degree d.
- **Notes:** Subset of `get_ncbasis` for specific degree.

---

## 7. Optimization Setup

### `polyopt`
- **Input:** `objective::Polynomial{A,T,C}`, `registry::VariableRegistry{A,T}`; optional `eq_constraints`, `ineq_constraints`
- **Output:** `PolyOpt{A,T,C}`
- **Notes:** Only accepts `Polynomial` - user must convert wrapper types (PauliMonomial, PhysicsMonomial) to Polynomial first.

---

## 8. Solver Pipeline

### `cs_nctssos`
- **Input:** `pop::PolyOpt{A,T,C}`, `config::SolverConfig`
- **Output:** `PolyOptResult`
- **Notes:** Main solver entry point. Internally uses algebra-appropriate basis types.

### `moment_relax`
- **Input:** `pop::PolyOpt`, `corr_sparsity`, `term_sparsities`
- **Output:** `MomentProblem{A,T,M,P}` where M = basis element type
- **Notes:** Moment matrix keys match basis element type (`PauliMonomial`/`PhysicsMonomial`/`Monomial`). JuMP moment variable deduplication: maintain `Dict{Monomial{A,T}, JuMP.VariableRef}` mapping canonical monomials to JuMP variables; when building moment matrix entries, convert wrapper keys to affine expressions over these base variables (e.g., `PauliMonomial` with phase_k=2 maps to `-1 * base_var`).

### `sos_dualize`
- **Input:** `mp::MomentProblem{A,T,M,P}`
- **Output:** `SOSProblem`
- **Notes:** Dualizes symbolic moment problem to SOS form for JuMP.

---

## 9. Sparsity

### `CorrelativeSparsity{A,T,P,M}`
- **Type parameters:** `A<:AlgebraType`, `T<:Integer`, `P<:Polynomial`, `M<:AbstractMonomial`
- **Fields:**
  - `cliques::Vector{Vector{T}}` - variable index groups
  - `clq_mom_mtx_bases::Vector{Vector{M}}` - basis elements (wrapper types: PauliMonomial for Pauli, Monomial for others)
  - `clq_localizing_mtx_bases::Vector{Vector{Vector{M}}}` - localizing bases
- **Notes:** M matches get_ncbasis return type.

### `TermSparsity{M}`
- **Type parameters:** `M<:AbstractMonomial`
- **Fields:**
  - `block_bases::Vector{Vector{M}}` - block decomposition bases
- **Notes:** M matches CorrelativeSparsity basis element type.

### `correlative_sparsity`
- **Input:** `pop::PolyOpt{A,T,C}`, `order::Int`, `cs_algo`
- **Output:** `CorrelativeSparsity{A,T,Polynomial{A,T,C},M}` where M = algebra-appropriate basis type
- **Notes:** Clique decomposition. Basis elements typed per algebra.

### `term_sparsities`
- **Input:** `activated_supp`, `constraint`, `mom_bases::Vector{M}`, `loc_bases`, `ts_algo`
- **Output:** `Vector{TermSparsity{M}}`
- **Notes:** Block decomposition for term sparsity. M propagates from input bases.

---

## 10. State Polynomial Types

**Supported algebras:** ProjectorAlgebra, UnipotentAlgebra, NonCommutativeAlgebra only.
**NOT supported:** PauliAlgebra, FermionicAlgebra, BosonicAlgebra (complex phases/expansion).
**Enforcement:** Type parameter constraint `A <: Union{ProjectorAlgebra, UnipotentAlgebra, NonCommutativeAlgebra}` on state types.

### `StateSymbol{ST,A,T}`
- **Fields:** `mono::Monomial{A,T}` (canonicalized)
- **Notes:** Stores raw Monomial (not wrappers). Canonicalization depends on ST.

### `StateWord{ST,A,T}`
- **Fields:** `state_syms::Vector{StateSymbol{ST,A,T}}` (sorted)
- **Notes:** Product of expectations ⟨M1⟩⟨M2⟩..., commutative so sorted.

### `NCStateWord{ST,A,T}`
- **Fields:** `sw::StateWord{ST,A,T}`, `nc_word::Monomial{A,T}`
- **Notes:** StateWord × nc_word. nc_word is raw Monomial (simple algebras only).

### `NCStatePolynomial{C,ST,A,T}`
- **Fields:** coefficients and NCStateWords
- **Notes:** Polynomial over NCStateWords with coefficient type C.

---

## 11. Internal Helpers

### Validation functions (NEW - for Monomial constructor)

| Function | Purpose |
|----------|---------|
| `_validate_pauli_word!(word)` | Check ≤1 op/site, sites sorted. Throws on violation. |
| `_validate_fermionic_word!(word)` | Check normal-ordered. Throws on violation. |
| `_validate_bosonic_word!(word)` | Check normal-ordered. Throws on violation. |
| `_validate_projector_word!(word)` | Check no P² terms, sites sorted. Throws on violation. |
| `_validate_unipotent_word!(word)` | Check no U² terms, sites sorted. Throws on violation. |
| `_validate_nc_word!(word)` | Check sites sorted. Throws on violation. |

### Simplification functions (existing - for wrapper constructors)

| Function | Purpose | Returns |
|----------|---------|---------|
| `_simplify_pauli_word!(word)` | Canonicalize Pauli word | `(canonical_word, phase_k::UInt8)` |
| `_simplify_fermionic_word!(word)` | Normal-order fermionic | `Vector{Tuple{Int, Vector{T}}}` (coef, mono pairs) |
| `_simplify_bosonic_word!(word)` | Normal-order bosonic | `Vector{Tuple{Int, Vector{T}}}` (coef, mono pairs) |
| `_simplify_projector_word!(word)` | Apply P²=P, sort | `Vector{T}` (canonical word) |
| `_simplify_unipotent_word!(word)` | Apply U²=I, sort | `Vector{T}` (canonical word) |
| `_simplify_nc_word!(word)` | Sort by site | `Vector{T}` (canonical word) |

### Usage pattern

```julia
# Monomial{Pauli} constructor
function Monomial{PauliAlgebra,T}(word::Vector{T}) where T
    _validate_pauli_word!(word)  # throws if invalid
    new{PauliAlgebra,T}(copy(word))
end

# PauliMonomial constructor
function PauliMonomial(word::Vector{T}) where T
    canonical_word, phase_k = _simplify_pauli_word!(copy(word))
    mono = Monomial{PauliAlgebra,T}(canonical_word)  # validated by simplify
    PauliMonomial{T}(mono, phase_k)
end
```

---

## 12. Type Conversion

### `Polynomial(::PauliMonomial)`
- **Input:** `pm::PauliMonomial{T}`
- **Output:** `Polynomial{PauliAlgebra,T,ComplexF64}`
- **Notes:** Converts PauliMonomial to single-term Polynomial (phase → coefficient).

### `Polynomial(::PhysicsMonomial)`
- **Input:** `pm::PhysicsMonomial{A,T}`
- **Output:** `Polynomial{A,T,ComplexF64}` (Int coeffs promoted)
- **Notes:** Converts PhysicsMonomial sum to Polynomial with ComplexF64 coefficients.

---

## Interview Log

| Date | Function | Decision |
|------|----------|----------|
| 2026-01-05 | create_pauli_variables | Returns PauliMonomial{T}, subscripts::AbstractVector |
| 2026-01-05 | PauliMonomial | `phase_k::UInt8` encoding algebraic phases only (i, -1, -i, 1) |
| 2026-01-05 | PhysicsMonomial | coeffs::Vector{Int}, monos::Vector{Monomial} for sum of normal-ordered terms |
| 2026-01-05 | Monomial{A,T} | Constructor validates canonical form (throws if invalid) |
| 2026-01-05 | simplify(Monomial) | Identity - monomials always canonical |
| 2026-01-05 | PauliMonomial × PauliMonomial | Returns PauliMonomial (closed) |
| 2026-01-05 | PauliMonomial × scalar | Returns Polynomial (not absorbed into phase) |
| 2026-01-05 | get_ncbasis | Returns algebra-appropriate wrapper types |
| 2026-01-05 | polyopt | Accepts Polynomial only (user converts wrappers) |
| 2026-01-05 | MomentProblem keys | Match basis element type |
| 2026-01-05 | PauliMonomial equality | Include phase in ==, isless, hash |
| 2026-01-05 | Type hierarchy | Flat AbstractMonomial{A,T} supertype for all |
| 2026-01-05 | PhysicsMonomial → Poly | Promote Int coeffs to ComplexF64 |
| 2026-01-05 | State polynomials | Only support Projector, Unipotent, NC algebras (drop Pauli, Fermi, Boson) |
| 2026-01-05 | StateSymbol | Stores Monomial{A,T} (not wrappers) |
| 2026-01-05 | NCStateWord.nc_word | Uses `Monomial{A,T}` only (simple algebras only) |
| 2026-01-05 | Sparsity bases | Match get_ncbasis return type (wrapper types) |
| 2026-01-05 | _simplify vs _validate | Separate: _validate checks only (for Monomial), _simplify modifies (for wrappers) |
