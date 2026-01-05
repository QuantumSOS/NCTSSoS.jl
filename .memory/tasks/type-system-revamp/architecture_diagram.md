# Type System Architecture Diagram

## 1. Type Hierarchy

```
                            AbstractMonomial{A,T}
                                     │
                 ┌───────────────────┼───────────────────┐
                 │                   │                   │
                 ▼                   ▼                   ▼
         Monomial{A,T}       PauliMonomial{T}    PhysicsMonomial{A,T}
         ┌──────────┐        ┌─────────────┐     ┌────────────────────┐
         │word::    │        │mono::       │     │coeffs::Vector{Int} │
         │Vector{T} │        │Monomial{    │     │monos::Vector{      │
         └──────────┘        │ Pauli,T}    │     │ Monomial{A,T}}     │
              │              │phase::      │     └────────────────────┘
              │              │ComplexF64   │              │
              │              └─────────────┘              │
              │                    │                      │
              │              A = PauliAlgebra       A ∈ {Fermionic,
              │                                          Bosonic}
              │
    ┌─────────┴─────────┬─────────────────┬──────────────┐
    │                   │                 │              │
    ▼                   ▼                 ▼              ▼
Pauli              Fermionic          Bosonic      Projector/
Algebra            Algebra            Algebra      Unipotent/NC
```

## 2. Constructor Flow

```
                         Raw Word (Vector{T})
                                │
                ┌───────────────┴───────────────┐
                │                               │
                ▼                               ▼
    ┌─────────────────────┐         ┌─────────────────────┐
    │  Monomial{A,T}()    │         │  Wrapper Types      │
    │  (direct construct) │         │  (canonicalizing)   │
    └─────────────────────┘         └─────────────────────┘
                │                               │
                ▼                               ▼
    ┌─────────────────────┐         ┌─────────────────────┐
    │ _validate_*_word!() │         │ _simplify_*_word!() │
    │  • check invariants │         │  • canonicalize     │
    │  • throw if invalid │         │  • return canonical │
    └─────────────────────┘         └─────────────────────┘
                │                               │
                ▼                               ▼
    ┌─────────────────────┐         ┌─────────────────────┐
    │ Canonical Monomial  │◄────────│ Wrapper with        │
    │ (invariants hold)   │         │ canonical mono      │
    └─────────────────────┘         └─────────────────────┘


    Example: Pauli

    [1,2,1,3]  ─────────────────────────────────────┐
        │                                           │
        ▼                                           ▼
    Monomial{Pauli}([1,2,1,3])              PauliMonomial([1,2,1,3])
        │                                           │
        ▼                                           ▼
    _validate_pauli_word!                   _simplify_pauli_word!
        │                                           │
        ▼                                           ▼
    ERROR! (σx appears twice)               ([2,3], phase=i)
                                                    │
                                                    ▼
                                            PauliMonomial{T}(
                                              mono=Monomial{Pauli}([2,3]),
                                              phase=im
                                            )
```

## 3. Multiplication Return Types

```
                    Monomial{A} × Monomial{A}
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        ▼                     ▼                     ▼
   A = Pauli            A = Fermi/Boson      A = Proj/Unip/NC
        │                     │                     │
        ▼                     ▼                     ▼
┌───────────────┐    ┌────────────────┐    ┌───────────────┐
│PauliMonomial  │    │PhysicsMonomial │    │Monomial{A,T}  │
│(may have      │    │(sum of normal- │    │(stays         │
│ phase ≠ 1)    │    │ ordered terms) │    │ canonical)    │
└───────────────┘    └────────────────┘    └───────────────┘


    Wrapper × Wrapper (closed operations)

    PauliMonomial × PauliMonomial ──────► PauliMonomial
           │
           ▼
    phases multiply, monos multiply,
    result re-canonicalized


    PhysicsMonomial × PhysicsMonomial ──► PhysicsMonomial
           │
           ▼
    distribute over sums,
    normal-order each product,
    collect like terms
```

## 4. Polynomial Integration

```
                        Polynomial{A,T,C}
                    ┌───────────────────────┐
                    │terms::Vector{Term{    │
                    │  Monomial{A,T}, C}}   │
                    └───────────────────────┘
                              ▲
                              │
            ┌─────────────────┼─────────────────┐
            │                 │                 │
    ┌───────┴───────┐ ┌───────┴───────┐ ┌───────┴───────┐
    │PauliMonomial  │ │PhysicsMonomial│ │  Monomial     │
    │    + scalar   │ │  conversion   │ │   + scalar    │
    └───────────────┘ └───────────────┘ └───────────────┘
            │                 │                 │
            ▼                 ▼                 ▼
    Polynomial{Pauli,   Polynomial{A,T,    Term{Mono,C}
    T,ComplexF64}       ComplexF64}             │
            │                 │                 │
            └─────────────────┴─────────────────┘
                              │
                              ▼
                    polyopt(poly, registry)
                              │
                              ▼
                      PolyOpt{A,T,C}
```

## 5. Optimization Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│                        Variable Creation                             │
│  create_pauli_variables(1:n) ──► (Registry, (σx[], σy[], σz[]))     │
│                                         │                            │
│                                         ▼                            │
│                               Vector{PauliMonomial{T}}               │
└─────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      Build Polynomial                                │
│   objective = Σ cᵢ * (pm₁ * pm₂ * ...)                              │
│                          │                                           │
│                          ▼                                           │
│              Polynomial{PauliAlgebra,T,ComplexF64}                   │
└─────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     polyopt(objective, registry)                     │
│                              │                                       │
│                              ▼                                       │
│                      PolyOpt{Pauli,T,ComplexF64}                     │
└─────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    correlative_sparsity(pop, order, algo)            │
│                              │                                       │
│         ┌────────────────────┴────────────────────┐                 │
│         ▼                                         ▼                 │
│   clq_mom_mtx_bases                    clq_localizing_mtx_bases     │
│   Vector{Vector{PauliMonomial{T}}}     Vector{Vector{Vector{PM}}}   │
│                              │                                       │
│                              ▼                                       │
│              CorrelativeSparsity{Pauli,T,Poly,PauliMonomial}        │
└─────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     term_sparsities(...)                             │
│                              │                                       │
│                              ▼                                       │
│              Vector{TermSparsity{PauliMonomial{T}}}                 │
│                   block_bases: Vector{Vector{PM}}                   │
└─────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     moment_relax(pop, cs, ts)                        │
│                              │                                       │
│                              ▼                                       │
│        MomentProblem{Pauli,T,PauliMonomial{T},Polynomial}           │
│                                                                      │
│   Moment matrix keys: PauliMonomial{T} (includes phase!)            │
└─────────────────────────────────────────────────────────────────────┘
                                          │
                                          ▼
┌─────────────────────────────────────────────────────────────────────┐
│                     sos_dualize(mp)                                  │
│                              │                                       │
│                              ▼                                       │
│                    SOSProblem (JuMP model)                          │
│                              │                                       │
│                              ▼                                       │
│                    optimize!(model)                                  │
│                              │                                       │
│                              ▼                                       │
│                     PolyOptResult                                    │
└─────────────────────────────────────────────────────────────────────┘
```

## 6. State Polynomial Types (Restricted)

```
    Supported Algebras Only:
    ┌─────────────────────────────────────────┐
    │  ProjectorAlgebra                       │
    │  UnipotentAlgebra                       │
    │  NonCommutativeAlgebra                  │
    └─────────────────────────────────────────┘

    NOT Supported (dropped):
    ┌─────────────────────────────────────────┐
    │  PauliAlgebra      (complex phases)     │
    │  FermionicAlgebra  (expansion)          │
    │  BosonicAlgebra    (expansion)          │
    └─────────────────────────────────────────┘


    State Type Hierarchy:

    StateSymbol{ST,A,T}              ⟨M⟩
         │
         ▼
    StateWord{ST,A,T}                ⟨M₁⟩⟨M₂⟩...⟨Mₖ⟩
         │                           (sorted, commutative)
         ▼
    NCStateWord{ST,A,T}              ⟨M₁⟩⟨M₂⟩ × nc_word
    ┌──────────────────┐
    │sw::StateWord     │
    │nc_word::Monomial │◄──── Raw Monomial{A,T} (simple algebras)
    └──────────────────┘
         │
         ▼
    NCStatePolynomial{C,ST,A,T}      Σ cᵢ × NCStateWord
```

## 7. Internal Helper Split

```
    ┌─────────────────────────────────────────────────────────────┐
    │                    _validate_*_word!()                       │
    │                                                              │
    │  Purpose: Check invariants ONLY (no modification)           │
    │  Called by: Monomial{A,T} constructor                       │
    │  Behavior: Throw error if invariant violated                │
    │                                                              │
    │  _validate_pauli_word!(word)     ≤1 op/site, sorted         │
    │  _validate_fermionic_word!(word) normal-ordered             │
    │  _validate_bosonic_word!(word)   normal-ordered             │
    │  _validate_projector_word!(word) no P², sorted              │
    │  _validate_unipotent_word!(word) no U², sorted              │
    │  _validate_nc_word!(word)        sorted by site             │
    └─────────────────────────────────────────────────────────────┘

    ┌─────────────────────────────────────────────────────────────┐
    │                    _simplify_*_word!()                        │
    │                                                              │
    │  Purpose: Canonicalize raw word                             │
    │  Called by: PauliMonomial, PhysicsMonomial constructors     │
    │  Behavior: Modify word in-place, return canonical form      │
    │                                                              │
    │  _simplify_pauli_word!(word)     → (canonical, phase)       │
    │  _simplify_fermionic_word!(word) → [(coef, mono), ...]      │
    │  _simplify_bosonic_word!(word)   → [(coef, mono), ...]      │
    │  _simplify_projector_word!(word) → canonical_word           │
    │  _simplify_unipotent_word!(word) → canonical_word           │
    │  _simplify_nc_word!(word)        → canonical_word           │
    └─────────────────────────────────────────────────────────────┘
```

## 8. Key Design Principles

```
┌─────────────────────────────────────────────────────────────────────┐
│  1. TYPE SAFETY VIA CONSTRUCTION                                    │
│     • Monomial{A,T} is ALWAYS canonical (constructor validates)     │
│     • PauliMonomial{T} wraps canonical mono + phase                 │
│     • PhysicsMonomial{A,T} stores sum of canonical monos            │
│     • Invalid states are UNREPRESENTABLE                            │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│  2. ALGEBRA-DEPENDENT DISPATCH                                       │
│     • Multiplication returns type appropriate to algebra            │
│     • Basis construction returns wrapper types                      │
│     • simplify() on Monomial is identity (already canonical)        │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│  3. PHASE SEMANTICS                                                  │
│     • PauliMonomial.phase: algebraic only {1, -1, i, -i}            │
│     • Scalar multiplication → Polynomial (not absorbed)             │
│     • Equality includes phase (pm1 == pm2 iff mono AND phase equal) │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│  4. CLEAN BOUNDARIES                                                 │
│     • polyopt() accepts Polynomial only (explicit conversion)       │
│     • State polys restricted to simple algebras                     │
│     • _validate vs _simplify clearly separated                      │
└─────────────────────────────────────────────────────────────────────┘
```
