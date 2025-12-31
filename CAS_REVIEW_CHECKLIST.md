# NCTSSoS.jl CAS Review Checklist

## 0. Reference Oracle: NCTSSOS (local)

**Local repo**: `/Users/yushengzhao/projects/NCTSSOS`

**Goal (CAS correctness)**: For operations on variables/monomials, `NCTSSoS`
should produce the same structural outcomes as an equivalent construction in
`NCTSSOS` (used only as an oracle for the *free noncommutative* word algebra and
for simple rewrite-style constraints like `"unipotent"`).

### What counts as "the same"
- **Primary comparison object**: expanded *word of variable indices* (tuple/vector
  of generator IDs), plus scalar coefficient when relevant.
- **Do not compare pretty-printing** (`show`/`string`) except as a debugging aid.

### Mapping / equivalences used for comparison

#### Variables
- `NCTSSoS`: variables returned by `create_*_variables` are **degree-1**
  `Monomial`s (there is no separate variable type).
- `NCTSSOS`: variables come from `DynamicPolynomials.@ncpolyvar x[1:n]`.

#### Word (monomial)
- `NCTSSoS`: `Monomial{A,T}([idx...])` stores `.word::Vector{T}` directly.
- `NCTSSOS`: represent a word `w = [i1,i2,...,ik]` as `prod(x[w])` and extract its
  `monomials(...)`/`coefficients(...)` or expand `DP.Monomial` fields like
  `w.vars` and `w.z` (see `NCTSSOS/src/utils.jl:reduce` for expansion logic).

#### Adjoint / star
- `NCTSSoS`: `adjoint(m)` / `m'` (reverse word; reverse+negate for signed types).
- `NCTSSOS`: `star(w::Mono)` in `NCTSSOS/src/utils.jl` (reverse word, respecting
  its exponent blocks).

#### Canonicalization (word-level)
- `NCTSSoS`: `symmetric_canon`, `cyclic_canon`, `cyclic_symmetric_canon`
  (see `src/algorithms/canonicalization.jl`).
- `NCTSSOS`: `_sym_canon`, `_cyclic_canon`
  (see `NCTSSOS/src/utils.jl`).

#### Constraint-like simplification (word-level oracle)
- `NCTSSOS/src/utils.jl:constraint_reduce!` supports `"unipotent"` (removes
  adjacent equal pairs) and "projector-like" idempotency behavior for other
  values (removes adjacent duplicates).

### Caveat: site-based commutation in NCTSSoS
`NCTSSoS` encodes "site" inside unsigned indices and uses site-aware sorting in:
- `simplify(::Monomial{NonCommutativeAlgebra,<:Unsigned})`
- `symmetric_canon(::Vector{<:Unsigned})` and `cyclic_symmetric_canon(::Vector{<:Unsigned})`

`NCTSSOS` does **not** have this site encoding. For direct oracle comparisons:
- Prefer **single-site registries** (one prefix group), so site-aware sorting is
  a no-op, or
- Compare after mapping `NCTSSoS` indices to operator IDs using
  `decode_operator_id` (and ignoring site), when appropriate for the specific
  check.

---

## 1. AlgebraType (Abstract Singleton)

**File**: `src/types/algebra.jl`

**Purpose**: Dispatch tag for algebra-specific simplification rules

### Functions
- [x] `default_coeff_type(::Type{A}) → Type` — Returns default coefficient type (`ComplexF64` for Pauli, `Float64` for others)
- [x] `show(io, ::A)` — Display algebra name

### Concrete Subtypes
- [x] `NonCommutativeAlgebra` — no simplification within site
- [x] `PauliAlgebra` — σ²=I, cyclic products (σxσy=iσz)
- [x] `FermionicAlgebra` — anticommutation {a,a†}=1, a²=0
- [x] `BosonicAlgebra` — commutation [b,b†]=1
- [x] `ProjectorAlgebra` — idempotent P²=P
- [x] `UnipotentAlgebra` — involution U²=I

---

## 2. VariableRegistry{A, T}

**File**: `src/types/registry.jl`

**Purpose**: Bidirectional mapping between indices and symbols

### Core Functions
- [x] `length(reg) → Int` — Number of variables
- [x] `getindex(reg, idx::T) → Symbol` — Get symbol by index
- [x] `getindex(reg, sym::Symbol) → T` — Get index by symbol
- [x] `in(sym, reg) → Bool` — Check if symbol exists
- [x] `symbols(reg) → Vector{Symbol}` — All symbols in sorted order
- [x] `indices(reg) → Vector{T}` — All indices in sorted order
- [x] `subregistry(reg, subset) → VariableRegistry` — Create subset registry
- [x] `show(io, reg)` — Display registry

### Factory Functions (return `(registry, monomials)`)
- [x] `create_pauli_variables(1:N)` — Creates `σx, σy, σz` per site
- [x] `create_fermionic_variables(1:N)` — Creates `a, a†` per mode
- [x] `create_bosonic_variables(1:N)` — Creates `c, c†` per mode
- [x] `create_projector_variables([("P", 1:N)])` — Creates `P` per site
- [x] `create_unipotent_variables([("U", 1:N)])` — Creates `U` per site
- [x] `create_noncommutative_variables([("x", 1:N)])` — Creates `x` per site

---

## 3. Monomial{A, T} <: AbstractMonomial

**File**: `src/types/monomial.jl`

**Purpose**: Word representation of non-commutative operator product

**Review reminder**: In `NCTSSoS`, "variables" are degree-1 monomials produced by
`create_*_variables` (so variable-level behavior is monomial behavior).

### Arithmetic
- [x] `*(m1, m2) → Monomial` — Concatenate words (no simplification)
- [x] `^(m, n::Int) → Monomial` — Repeated multiplication
- [x] `+(m1, m2) → Polynomial` — Create 2-term polynomial
- [x] `-(m1, m2) → Polynomial` — Subtract monomials
- [x] `-(m) → Term` — Negate (coefficient -1)

### Adjoint
- [x] `adjoint(m)` / `m'` → Monomial — Hermitian conjugate (reverse, negate signed)

### Properties
- [x] `degree(m) → Int` — Word length
- [x] `isone(m) → Bool` — Empty word?
- [x] `one(::Type{Monomial{A,T}}) → Monomial` — Identity monomial
- [x] `variable_indices(m) → Set{T}` — Set of indices in word

### Comparison
- [x] `isless(m1, m2) → Bool` — Graded lexicographic order
- [x] `==(m1, m2) → Bool` — Word equality
- [x] `hash(m, h) → UInt` — Hash including algebra type

### Iteration Protocol
- [x] `iterate(m)` — Yields single `(1, m)` pair
- [x] `length(m) → 1` — Single element for iteration
- [x] `eltype(::Type{Monomial{A,T}})` — Returns `Tuple{C, Monomial{A,T}}`

### Utilities
- [x] `expval(m) → Monomial` — Identity (for API compat with StatePolynomial)
- [x] `coeff_type(::Type{...}) → Type` — Default coefficient type
- [x] `show(io, m)` — Display with registry if available

### NCTSSOS oracle notes (free NC only)
- [x] `*(m1,m2)` matches `prod(x[word1])*prod(x[word2])` word concatenation
- [x] `m'` matches `NCTSSOS.star(m)` (reverse of expanded word)
- [x] `degree(m)` matches expanded word length in `NCTSSOS`

---

## 4. Term{M, C}

**File**: `src/types/term.jl`

**Purpose**: Coefficient-monomial pair from simplification

### Arithmetic
- [x] `*(c::Number, t) → Term` — Scale coefficient
- [x] `*(t, c::Number) → Term` — Scale coefficient (commutative)
- [x] `-(t) → Term` — Negate coefficient

### Predicates
- [x] `isone(t) → Bool` — Coefficient=1 and identity monomial?
- [x] `iszero(t) → Bool` — Zero coefficient?

### Identity Elements
- [x] `one(::Type{Term{M,C}}) → Term` — Identity term
- [x] `zero(::Type{Term{M,C}}) → Term` — Zero term

### Comparison
- [x] `==(t1, t2) → Bool` — Same coef and monomial?
- [x] `hash(t, h) → UInt` — Hash

### Iteration Protocol
- [x] `iterate(t)` — Yields single `(coef, mono)` pair
- [x] `length(t) → 1` — Single element
- [x] `eltype(::Type{Term{M,C}})` — Returns `Tuple{C, M}`

### Utilities
- [x] `coeff_type(::Type{...}) → Type` — Returns `C`
- [x] `show(io, t)` — Display

---

## 5. Polynomial{A, T, C} <: AbstractPolynomial{C}

**File**: `src/types/polynomial.jl`

**Purpose**: Sorted, deduplicated sum of terms

### Construction
- [x] `Polynomial(m::Monomial)` — From monomial (coef=1)
- [x] `Polynomial(t::Term)` — From single term
- [x] `Polynomial(terms::Vector{Term})` — From term vector (auto-processes)
- [x] `Polynomial{A,T,C}(c::Number)` — Constant polynomial

### Accessors
- [x] `terms(p) → Vector{Term}` — Get terms
- [x] `coefficients(p) → Vector{C}` — Extract coefficients
- [x] `monomials(p) → Vector{Monomial}` — Extract monomials
- [x] `degree(p) → Int` — Max monomial degree (`-Inf` if zero)
- [x] `variable_indices(p) → Set{T}` — All variable indices

### Predicates
- [x] `iszero(p) → Bool` — No terms?
- [x] `isone(p) → Bool` — Single identity term with coef=1?

### Arithmetic: Polynomial-Polynomial
- [x] `+(p1, p2) → Polynomial` — Add polynomials
- [x] `-(p1, p2) → Polynomial` — Subtract
- [x] `*(p1, p2) → Polynomial` — Multiply (with simplification)

### Arithmetic: Scalar
- [x] `*(c, p) → Polynomial` — Scalar multiplication (left)
- [x] `*(p, c) → Polynomial` — Scalar multiplication (right)
- [x] `/(p, c) → Polynomial` — Scalar division
- [x] `+(p, c) → Polynomial` — Add constant
- [x] `+(c, p) → Polynomial` — Add constant (commutative)
- [x] `-(p, c) → Polynomial` — Subtract constant
- [x] `-(c, p) → Polynomial` — Constant minus polynomial

### Arithmetic: Monomial
- [x] `+(p, m) → Polynomial` — Add monomial to polynomial
- [x] `+(m, p) → Polynomial` — Add monomial to polynomial
- [x] `-(p, m) → Polynomial` — Subtract monomial
- [x] `-(m, p) → Polynomial` — Monomial minus polynomial
- [x] `*(m, p) → Polynomial` — Monomial times polynomial
- [x] `*(p, m) → Polynomial` — Polynomial times monomial
- [x] `*(c, m) → Polynomial` — Scalar times monomial

### Unary Operations
- [x] `-(p) → Polynomial` — Negate all coefficients
- [x] `^(p, n::Int) → Polynomial` — Power by squaring
- [x] `adjoint(p)` / `p'` → Polynomial — Hermitian conjugate

### Identity Elements
- [x] `zero(::Type{Polynomial{...}}) → Polynomial` — Zero polynomial
- [x] `zero(p) → Polynomial` — Zero of same type
- [x] `one(::Type{Polynomial{...}}) → Polynomial` — Identity polynomial
- [x] `one(p) → Polynomial` — Identity of same type
- [x] `copy(p) → Polynomial` — Shallow copy

### Comparison
- [x] `==(p1, p2) → Bool` — Term-by-term equality
- [x] `isless(p1, p2) → Bool` — Graded lex ordering
- [x] `hash(p, h) → UInt` — Hash

### Simplification
- [x] `simplify(p) → Polynomial` — Apply algebra rules to all terms
- [x] `_process_terms(terms, C) → Vector{Term}` — Sort, dedupe, remove zeros
- [x] `_collect_simplified_terms!(result, coef, simplified)` — Helper for simplify

### Type Conversion
- [x] `convert(::Type{Polynomial{A,T,C2}}, p) → Polynomial` — Coefficient type conversion

### Iteration Protocol
- [x] `iterate(p)` — Yields `(coef, mono)` pairs
- [x] `iterate(p, state)` — Continue iteration
- [x] `length(p) → Int` — Number of terms
- [x] `eltype(::Type{Polynomial{A,T,C}})` — Returns `Tuple{C, Monomial{A,T}}`

### Utilities
- [x] `coeff_type(::Type{...}) → Type` — Returns `C`
- [x] `show(io, p)` — Display as sum

---

## 6. Simplification (per Algebra)

### NonCommutativeAlgebra
**File**: `src/simplification/noncommutative.jl`

- [x] `simplify(m::Monomial{NonCommutativeAlgebra}) → Monomial` — Sort by site
- [x] `simplify!(m) → Monomial` — In-place version
- [x] `_simplify_nc_word!(word) → Vector` — In-place word simplification

### PauliAlgebra
**File**: `src/simplification/pauli.jl`

- [x] `simplify(m::Monomial{PauliAlgebra}) → Term` — Apply σ²=I, cyclic products
- [x] `_pauli_site(idx) → Int` — Extract site from index
- [x] `_pauli_type(idx) → Int` — Extract Pauli type (0=X, 1=Y, 2=Z)
- [x] `_pauli_index(site, type) → Int` — Create index from site+type
- [x] `_pauli_product(type1, type2) → (phase, result_type)` — Cyclic product table

### FermionicAlgebra
**File**: `src/simplification/fermionic.jl`

- [x] `simplify(m::Monomial{FermionicAlgebra}) → Polynomial` — Normal order with anticommutation

### BosonicAlgebra
**File**: `src/simplification/bosonic.jl`

- [x] `simplify(m::Monomial{BosonicAlgebra}) → Polynomial` — Normal order with commutation

### ProjectorAlgebra
**File**: `src/simplification/projector.jl`

- [x] `simplify(m::Monomial{ProjectorAlgebra}) → Monomial` — Remove consecutive duplicates (P²=P)

### UnipotentAlgebra
**File**: `src/simplification/unipotent.jl`

- [x] `simplify(m::Monomial{UnipotentAlgebra}) → Monomial` — Remove pairs (U²=I)

---

## 7. Index Encoding Utilities

**File**: `src/types/algebra.jl`

**Purpose**: Site-based index packing for commutation rules

- [x] `site_bits(::Type{T<:Unsigned}) → Int` — Bits for site encoding
- [x] `max_sites(::Type{T}) → Int` — Max sites for type
- [x] `max_operators(::Type{T}) → Int` — Max operators per site
- [x] `encode_index(T, op_id, site) → T` — Pack operator + site
- [x] `decode_site(idx::T) → Int` — Extract site
- [x] `decode_operator_id(idx::T) → Int` — Extract operator ID
- [x] `select_uint_type(n_ops, n_sites) → Type` — Choose smallest UInt

---

## 8. Algorithms

### Basis Generation
**File**: `src/algorithms/basis.jl`

#### What NCTSSoS implements
- [x] `_generate_all_words(indices, d)` — All words of *exact* length `d`
- [x] `get_ncbasis_deg(reg, d)` — For each word, build `Monomial` then `simplify`
  it, and return a `Polynomial` (one polynomial per input word)
- [x] `get_ncbasis(reg, d)` — Concatenate degrees `0:d`

#### NCTSSOS oracle
- [x] `NCTSSOS.get_ncbasis(n, d; ind=..., binary=false)` in
  `NCTSSOS/src/utils.jl` (returns `Vector{Vector{UInt16}}` words)

#### Review checks
- [x] Enumeration matches for free NC, single-site variables
  (`binary=false`, degree-by-degree counts `n^d`)
  - Verified in `test/polynomials/simplify.jl` "NCTSSOS Oracle: get_ncbasis counts"
- [x] Degree-0 basis contains identity word / identity polynomial
  - Verified in `test/polynomials/basis.jl` and oracle tests
- [x] If comparing multi-site `NCTSSoS` bases, document the mapping decision:
  compare encoded indices directly vs compare by `decode_operator_id`
  - Oracle tests use `decode_operator_id` for word extraction

### Canonicalization
**File**: `src/algorithms/canonicalization.jl`

#### What NCTSSoS implements (word-level)
- [x] `symmetric_canon(word)` — `min(word, reverse(word))`
- [x] `cyclic_canon(word)` — Lexicographically minimal cyclic rotation
- [x] `cyclic_symmetric_canon(word)` — `min(cyclic_canon(word), cyclic_canon(reverse(word)))`

#### Monomial + Polynomial interfaces
- [x] `canonicalize(m; cyclic=false|true)` — symmetric vs cyclic+symmetric
- [x] `canonicalize(p; cyclic=false|true)` — canonicalize monomials, then merge
  like terms via polynomial constructor

#### NCTSSOS oracle
- [x] `_sym_canon(word)` in `NCTSSOS/src/utils.jl`
  - Verified in `test/polynomials/simplify.jl` "NCTSSOS Oracle: symmetric_canon (word-level)"
- [x] `_cyclic_canon(word)` in `NCTSSOS/src/utils.jl`
  - Verified in `test/polynomials/simplify.jl` "NCTSSOS Oracle: cyclic_canon (word-level)"
- [x] `min(_cyclic_canon(word), _cyclic_canon(reverse(word)))` corresponds to
  cyclic+symmetric mode
  - Verified in `test/polynomials/simplify.jl` "NCTSSOS Oracle: cyclic_symmetric_canon (word-level)"

#### Caveat: site-aware specializations in NCTSSoS
- [x] `symmetric_canon(::Vector{<:Unsigned})` and
  `cyclic_symmetric_canon(::Vector{<:Unsigned})` sort by `decode_site` before
  comparing; `NCTSSOS` does not model this.
  - Implementation in `src/algorithms/canonicalization.jl` lines 127-142 and 328-340
- [x] For oracle comparisons, prefer single-site words or explicitly state how
  you "forget the site" (e.g., map indices to `decode_operator_id` first).
  - Oracle tests use signed Int16 to avoid site-aware dispatch

---

## Review Notes

Use this section to record findings during review:

### Issues Found
<!-- 
- [ ] Issue description (file:line)
-->

### Suggestions
<!--
- [ ] Suggestion description
-->

### Questions
<!--
- [ ] Question about design/implementation
-->
