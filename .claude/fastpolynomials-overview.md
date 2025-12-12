# FastPolynomials: A Comprehensive Overview

## What is FastPolynomials?

FastPolynomials is a high-performance Julia library for **non-commutative polynomial optimization** within NCTSSoS. It's designed to be the fastest polynomial library for polynomial optimization in Julia, inspired by DynamicPolynomials.jl.

The key innovation: it handles **multiple algebra types** (Pauli spins, fermions, bosons, projectors, etc.) with **zero-overhead compile-time dispatch** by encoding algebra rules in the type system.

---

## Core Components

### 1. Algebra Type System (`algebra_types.jl`)

Six algebra types, each with distinct algebraic rules:

- **`NonCommutativeAlgebra`**: Generic (xy ≠ yx, no simplification)
- **`PauliAlgebra`**: Pauli spins (σᵢ² = I, anticommutation, cyclic products → complex phases)
- **`FermionicAlgebra`**: Fermions ({aᵢ, aⱼ†} = δᵢⱼ anticommutation)
- **`BosonicAlgebra`**: Bosons ([cᵢ, cⱼ†] = δᵢⱼ commutation)
- **`ProjectorAlgebra`**: Projectors (Pᵢ² = Pᵢ idempotency)
- **`UnipotentAlgebra`**: Unipotent (U² = I)

**Design principle:** Singleton types enable compile-time specialization—no runtime overhead.

---

### 2. Variable Registry (`variable_registry.jl`)

**Type:** `VariableRegistry{A<:AlgebraType, T<:Integer}`

**Purpose:** Bidirectional mapping between symbolic variable names and integer indices.

**Key features:**
- Supports contiguous (Vector) and sparse (Dict) index schemes
- Type-safe index encoding (bit-packing for site-based algebras)
- Variable creation functions for each algebra:
  ```julia
  create_pauli_variables(subscripts) → (registry, (σx, σy, σz))
  create_fermionic_variables(subscripts) → (registry, (a, a⁺))
  create_bosonic_variables(subscripts) → (registry, (c, c⁺))
  # ... etc for other algebras
  ```

**Index encoding strategies:**
- **Pauli**: Site-first ordering (σx₁, σy₁, σz₁, σx₂, ...)
- **Fermionic/Bosonic**: Signed integers (positive = annihilation, negative = creation)
- **Others**: Bit-packed site encoding

---

### 3. Monomial (`monomial.jl`)

**Type:** `Monomial{A<:AlgebraType, T<:Integer}`

**Core representation:**
```julia
struct Monomial{A,T}
    word::Vector{T}      # Sequence of variable indices
    hash::UInt64         # Precomputed hash for fast checks
end
```

**Key features:**
- **Word representation**: e.g., [1,3,1,3] = x₁z₁x₁z₁
- **Precomputed hash**: 50-100× faster inequality checks
- **Algebra type in type parameter**: Zero overhead
- **Flexible integer types**: UInt8/16/32/64 or Int8/16/32/64 depending on algebra

**Operations:**
- Multiplication: concatenates words
- Adjoint: reverses word (and negates for signed types)
- Ordering: Degree-first (graded), then lexicographic

---

### 4. Term (`term.jl`)

**Type:** `Term{M<:AbstractMonomial, C<:Number}` (mutable)

**Purpose:** Coefficient-monomial pair—the output of simplification.

```julia
mutable struct Term{M,C}
    coefficient::C      # Float64 or ComplexF64
    monomial::M         # Simplified monomial
end
```

**Coefficient types by algebra:**
- **Pauli**: `ComplexF64` (cyclic products introduce i phases)
- **Others**: `Float64` (real coefficients)

**Mutability rationale:** Enables efficient coefficient accumulation.

---

### 5. Polynomial (`polynomial.jl`)

**Type:** `Polynomial{A<:AlgebraType, T<:Integer, C<:Number}`

**Core representation:**
```julia
struct Polynomial{A,T,C}
    terms::Vector{Term{Monomial{A,T}, C}}  # Sorted, unique, non-zero
end
```

**Invariants (enforced by constructor):**
1. Terms sorted by monomial
2. No duplicate monomials
3. All coefficients non-zero

**Operations:** Full arithmetic (`+`, `-`, `*`, `/`, `^`), adjoint, degree queries, etc.

---

### 6. ComposedMonomial (`composed_monomial.jl`)

**Type:** `ComposedMonomial{Ts<:Tuple}`

**Purpose:** Tensor products across **different** algebra types (e.g., Pauli ⊗ Fermionic).

```julia
struct ComposedMonomial{Ts<:Tuple}
    components::Ts       # Tuple of Monomial{A,T} with different A
    hash::UInt64
end
```

---

## Architecture Overview

### Layer 1: Type Foundation (Compile-Time Configuration)

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         AlgebraType (Singleton)                         │
│  Abstract type defining algebraic rules via compile-time dispatch       │
├─────────┬──────────┬────────────┬────────────┬───────────┬──────────────┤
│ NonComm │  Pauli   │ Fermionic  │  Bosonic   │ Projector │  Unipotent   │
│         │ (σᵢ²=I)  │ ({aᵢ,aⱼ†}) │ ([cᵢ,cⱼ†]) │  (Pᵢ²=Pᵢ) │   (U²=I)     │
└─────────┴──────────┴────────────┴────────────┴───────────┴──────────────┘
```

**Purpose:** The `AlgebraType` singleton types are the foundation of the entire architecture. They enable zero-overhead compile-time dispatch—Julia's compiler can specialize every function for each algebra type without runtime cost. This layer defines the "physics" of the system: what simplification rules apply, what coefficient types are needed (real vs complex), and how indices should be encoded.

**Key Insight:** By making algebra types part of the type parameters (e.g., `Monomial{PauliAlgebra, UInt8}`), we get compile-time guarantees that prevent mixing incompatible algebra types and enable the compiler to optimize algebra-specific operations aggressively.

---

### Layer 2: Variable Management (Symbol ↔ Index Mapping)

```
┌─────────────────────────────────────────────────────────────────────────┐
│                  VariableRegistry{A<:AlgebraType, T<:Integer}           │
│  Bidirectional mapping: Symbol names ⟷ Integer indices                 │
│                                                                          │
│  • Contiguous scheme: Vector-based (dense indices)                      │
│  • Sparse scheme: Dict-based (arbitrary indices)                        │
│  • Algebra-specific encoding (bit-packing, signed indices)              │
└────────────────────────────────┬────────────────────────────────────────┘
                                 │
                    ┌────────────┴────────────┐
                    │   Variable Creation     │
                    │   create_*_variables()  │
                    └─────────────────────────┘
                                 │
                    Returns: (registry, variables)
```

**Purpose:** The `VariableRegistry` bridges the gap between user-facing symbolic names (like `σx[1]`, `a[2]`) and the internal integer representation used by monomials. It's the entry point for users—all polynomial construction begins by creating variables through the registry.

**Key Insight:** Different algebras need different index encoding strategies:
- **Pauli**: Site-first ordering (3 operators per site) for efficient commutation checks
- **Fermionic/Bosonic**: Signed integers where sign distinguishes creation vs annihilation
- **Others**: Bit-packed site indices for multi-site systems

The registry enforces type consistency and provides algebra-specific variable creation functions that return properly structured variables.

---

### Layer 3: Monomial Representation (Sequences of Variables)

```
                    ┌─────────────────────────────────────┐
                    │   AbstractMonomial (Abstract Type)  │
                    └──────────────┬──────────────────────┘
                                   │
              ┌────────────────────┴────────────────────┐
              │                                         │
    ┌─────────▼────────────┐              ┌────────────▼──────────────┐
    │  Monomial{A,T}       │              │  ComposedMonomial{Ts}     │
    │  Single algebra type │              │  Tensor product across    │
    │                      │              │  multiple algebra types   │
    │  word::Vector{T}     │              │                           │
    │  hash::UInt64        │              │  components::Tuple{...}   │
    │                      │              │  hash::UInt64             │
    │  Example:            │              │                           │
    │  [1,3,1,3]           │              │  Example:                 │
    │  = x₁z₁x₁z₁          │              │  (σₓ₁, a₁†a₂)             │
    │                      │              │  = Pauli ⊗ Fermionic      │
    └──────────────────────┘              └───────────────────────────┘
```

**Purpose:** Monomials represent **unsimplified** products of variables. A `Monomial{A,T}` is simply a sequence of variable indices (the "word") along with a precomputed hash for fast comparisons. They carry algebra type information but don't apply simplification rules yet.

**Key Insight:** The separation between single-algebra `Monomial` and multi-algebra `ComposedMonomial` is crucial:
- **Monomial**: Homogeneous—all variables from same algebra (e.g., all Pauli operators). Fast, simple representation.
- **ComposedMonomial**: Heterogeneous—tensor products of monomials from different algebras (e.g., Pauli ⊗ Fermionic for hybrid quantum systems). More complex but enables modeling systems with multiple operator types.

The precomputed hash enables O(1) inequality checks in typical cases—critical for polynomial term deduplication.

---

### Layer 4: Simplification Engine (Algebra-Specific Reduction)

```
                    ┌──────────────────────────────────┐
                    │  Monomial (Unsimplified)         │
                    │  e.g., [σₓ, σₓ, σᵧ] (raw word)   │
                    └────────────┬─────────────────────┘
                                 │
                                 │ simplify(monomial)
                                 │ Dispatch on AlgebraType
                                 ▼
            ┌────────────────────────────────────────────────┐
            │     Simplification Dispatch Layer              │
            │     (Compile-time selection based on A)        │
            ├──────┬─────────┬──────────┬─────────┬──────────┤
            │Pauli │Fermionic│ Bosonic  │Projector│Unipotent │
            │Sort+ │Anticomm │Commute+  │Remove   │Remove    │
            │Cyclic│ +Sign   │Correction│Duplicate│Pairs     │
            └──────┴─────────┴──────────┴─────────┴──────────┘
                                 │
                                 │ Returns: Term or Vector{Term}
                                 ▼
                    ┌──────────────────────────────────┐
                    │  Term{M, C} (Simplified)         │
                    │  coefficient::C × monomial::M    │
                    │  e.g., i × [σz] (simplified)     │
                    └──────────────────────────────────┘
```

**Purpose:** The simplification engine applies algebra-specific reduction rules to transform unsimplified monomials into canonical form. This is where the algebraic "physics" actually happens—Pauli operators get combined, fermionic operators get normal-ordered, etc.

**Key Insight:** Different algebras require different output structures:
- **Linear simplifications** (Pauli, Projector, Unipotent): One monomial → one term (coefficient may change, monomial gets reduced)
- **Non-linear simplifications** (Fermionic, Bosonic): One monomial → multiple terms (e.g., aᵢ†aⱼ → -aⱼaᵢ† + δᵢⱼ for bosons)

The output is always a `Term` or `Vector{Term}`, where each term pairs a coefficient with a simplified monomial. Coefficients can be real or complex depending on the algebra.

---

### Layer 5: Polynomial Representation (Collections of Terms)

```
            ┌─────────────────────────────────────────────┐
            │  Polynomial{A, T, C}                        │
            │                                             │
            │  terms::Vector{Term{Monomial{A,T}, C}}     │
            │                                             │
            │  Invariants:                                │
            │  1. Sorted by monomial (graded lex order)  │
            │  2. No duplicate monomials                  │
            │  3. All coefficients non-zero               │
            │                                             │
            │  Example:                                   │
            │  3.5×[x₁] + 2.1×[x₁,x₂] - 1.0×[x₂,x₁]      │
            └─────────────────────────────────────────────┘
                                 │
                                 │ Arithmetic operations
                                 │ (+, -, *, /, ^, adjoint)
                                 ▼
            ┌─────────────────────────────────────────────┐
            │  Result: New Polynomial (invariants hold)   │
            │  • Terms merged if monomials equal          │
            │  • Zero coefficients removed                │
            │  • Sorted order maintained                  │
            └─────────────────────────────────────────────┘
```

**Purpose:** `Polynomial` is the primary user-facing type—a collection of terms with strong invariants enforced by the constructor. All polynomial arithmetic goes through this type, ensuring mathematical correctness and optimal representation.

**Key Insight:** The three invariants are critical for performance and correctness:
1. **Sorted terms**: Enable binary search, fast equality checks, and predictable iteration
2. **No duplicates**: Ensure unique representation (2x + 3x becomes 5x)
3. **No zeros**: Avoid wasted space and simplify logic (0×x₁ is just removed)

These invariants are maintained across all operations—addition automatically merges like terms, multiplication produces a new sorted polynomial, etc.

---

### Layer 6: Higher-Level Structures (Basis Generation & Optimization)

```
┌────────────────────────────────────────────────────────────────────────┐
│                         High-Level API Layer                           │
├────────────────────────┬──────────────────────┬────────────────────────┤
│   Basis Generation     │   Canonicalization   │  State Polynomials     │
│                        │                      │                        │
│  get_ncbasis(reg, d)   │  symmetric_canon()   │  StatePolynomial{...}  │
│  get_ncbasis_deg(...)  │  cyclic_canon()      │  (quantum states)      │
│                        │  cyclic_sym_canon()  │                        │
│  Returns:              │                      │                        │
│  Vector{Polynomial}    │  Canonical words for │  Specialized type for  │
│                        │  optimization        │  density matrices      │
└────────────────────────┴──────────────────────┴────────────────────────┘
                                     │
                                     ▼
                    ┌──────────────────────────────────┐
                    │  NCTSSoS Optimization Framework  │
                    │  Sum-of-Squares Relaxations      │
                    └──────────────────────────────────┘
```

**Purpose:** These higher-level structures use the core polynomial types to enable optimization workflows. Basis generation creates complete sets of polynomials for moment relaxations, canonicalization exploits symmetry to reduce problem size, and state polynomials specialize the representation for quantum state optimization.

**Key Insight:** This layer is where FastPolynomials connects to the broader NCTSSoS framework:
- **Basis generation**: Creates all monomials up to degree d, simplifies them, and returns as polynomials—needed for SOS moment matrices
- **Canonicalization**: Finds canonical representatives for equivalence classes (e.g., x₁x₂ ≡ x₂x₁ under symmetry)—reduces optimization variables
- **State polynomials**: Special structure for density matrices (Hermitian, trace-normalized)—enables quantum state optimization

---

### Complete Data Flow Example: User Input → Optimized Polynomial

```
1. USER INPUT
   "Create Pauli variables for 2 sites"
   ↓
   create_pauli_variables(1:2)

2. REGISTRY CREATION
   VariableRegistry{PauliAlgebra, UInt8}
   Mapping: σx₁→1, σy₁→2, σz₁→3, σx₂→4, σy₂→5, σz₂→6
   Returns: (registry, (σx, σy, σz))

3. MONOMIAL CONSTRUCTION
   User writes: σx[1] * σy[1] * σx[2]
   Internal: Monomial{PauliAlgebra, UInt8}([1, 2, 4], hash)

4. SIMPLIFICATION
   simplify(monomial) dispatches to Pauli simplifier
   Rules: σx₁σy₁ = iσz₁, then σz₁σx₂ (different sites commute)
   Result: Term(im, Monomial([3, 4], hash))  ← Complex coefficient!

5. POLYNOMIAL CONSTRUCTION
   Polynomial{PauliAlgebra, UInt8, ComplexF64}(
     [Term(im, [3,4])]  ← Single term, invariants satisfied
   )

6. POLYNOMIAL ARITHMETIC
   p1 = im×σz₁σx₂
   p2 = 2.0×σz₁σx₂ + 3.0×σx₁
   p1 + p2 → merges like terms → (2.0+im)×σz₁σx₂ + 3.0×σx₁

7. BASIS GENERATION (for optimization)
   get_ncbasis(registry, 2) →
   All degree ≤2 polynomials: [I, σx₁, σy₁, ..., σx₁σx₂, ...]

8. SOS OPTIMIZATION (in NCTSSoS framework)
   Use basis to construct moment matrix, solve SDP
```

**Purpose:** This end-to-end flow shows how all layers work together to transform user input into optimized mathematical objects ready for sum-of-squares relaxations.

**Key Insight:** Each layer adds value:
- Layer 1 (AlgebraType): Enables type-safe algebra-specific operations
- Layer 2 (Registry): Translates user symbols to efficient integers
- Layer 3 (Monomial): Compact representation of variable products
- Layer 4 (Simplification): Applies algebraic rules to canonical form
- Layer 5 (Polynomial): Maintains invariants for correct arithmetic
- Layer 6 (Basis/Optimization): Connects to SOS framework

The architecture is designed so each layer has clear responsibilities and minimal coupling, enabling maintainability and extensibility while maintaining performance.

---

## Algebra Systems Deep Dive

| Algebra | Rules | Index Type | Coeff Type | Simplification Example |
|---------|-------|------------|------------|------------------------|
| **NonCommutative** | None | Unsigned | Real | No reduction |
| **Pauli** | σᵢ² = I, anticomm, cyclic | Site-first UInt | Complex | σₓσᵧ → iσz |
| **Fermionic** | {aᵢ, aⱼ†} = δᵢⱼ | Signed Int | Real | aᵢaⱼ → -aⱼaᵢ |
| **Bosonic** | [cᵢ, cⱼ†] = δᵢⱼ | Signed Int | Real | cᵢcⱼ† → cⱼ†cᵢ + δᵢⱼ |
| **Projector** | Pᵢ² = Pᵢ | Bit-packed UInt | Real | PᵢPᵢ → Pᵢ |
| **Unipotent** | U² = I | Bit-packed UInt | Real | UᵢUᵢ → I |

---

## Key Files

| File | Purpose |
|------|---------|
| `FastPolynomials.jl` | Module definition, includes, exports |
| `algebra_types.jl` | Singleton algebra types, site encoding |
| `variable_registry.jl` | Symbol ↔ index mapping, variable creation |
| `monomial.jl` | Single-algebra monomial type |
| `term.jl` | Coefficient-monomial pair |
| `polynomial.jl` | Polynomial type, arithmetic |
| `composed_monomial.jl` | Tensor products across algebras |
| `simplification/*.jl` | Algebra-specific simplification algorithms |
| `basis.jl` | Basis generation for optimization |
| `canonicalization.jl` | Symmetric/cyclic canonical forms |
| `state_*.jl` | State polynomial types for quantum states |

---

## Performance Highlights

1. **Precomputed Hashes**: Monomial equality checks are O(1) in typical cases
2. **Type Stability**: All algebra types in type parameters (no runtime dispatch)
3. **Bit-Packing**: Site-based algebras use bit fields for efficient encoding
4. **Minimal Integer Types**: Automatically selects smallest UInt/Int that fits

---

## Main API Surface

**Variable Creation:**
```julia
registry, (σx, σy, σz) = create_pauli_variables(1:3)
registry, (a, a⁺) = create_fermionic_variables(1:4)
```

**Polynomial Operations:**
```julia
degree(p), coefficients(p), monomials(p)
simplify(m), simplify!(m)
adjoint(p), is_symmetric(p)
```

**Basis Generation:**
```julia
basis = get_ncbasis(registry, degree)      # Up to degree d
basis_d = get_ncbasis_deg(registry, d)     # Exactly degree d
```

---

## Summary

FastPolynomials is a **highly specialized, type-safe polynomial library** that achieves high performance through:
- Zero-overhead compile-time dispatch via singleton types
- Algebra-specific simplification producing canonical forms
- Flexible index encoding strategies per algebra
- Strong invariants on polynomial representation
- Registry-based variable management

It forms the computational foundation for NCTSSoS's non-commutative sum-of-squares optimization over quantum operator algebras.
