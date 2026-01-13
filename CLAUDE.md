# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

Testing: see `TESTING.md`.

```bash
# Documentation
make servedocs         # Live-reload docs server
make init-docs         # Setup docs environment first
```

## Architecture

### Type System Hierarchy
```
AlgebraType (abstract)
├── MonoidAlgebra        # normal form = single monomial
│   ├── NonCommutativeAlgebra, ProjectorAlgebra, UnipotentAlgebra
├── TwistedGroupAlgebra  # normal form = scalar × monomial
│   └── PauliAlgebra
└── PBWAlgebra           # normal form = sum of monomials
    ├── FermionicAlgebra, BosonicAlgebra
```

### Core Types
- `NormalMonomial{A,T}` - Immutable canonical word
- `Monomial{A,T,C,W}` - User-facing wrapper (coeff + word(s))
- `Polynomial{A,T,C}` - Mutable mapping NormalMonomial → coefficient
- `VariableRegistry{A,T}` - Bidirectional symbol ↔ index mapping

### Type Parameter Convention
Algebra type `A` comes first (enables dispatch), then index type `T`, then coefficient type `C`:
```julia
Polynomial{A,T,C}  # A = algebra, T = index integer, C = coefficient
```

### Key Directories
- `src/types/` - Core type definitions (algebra, registry, monomial, polynomial)
- `src/simplification/` - Algebra-specific simplification rules (one file per algebra)
- `src/optimization/` - SDP relaxation pipeline (problem → sparsity → moment/sos → JuMP)
- `src/states/` - State polynomial types for quantum information
- `test/oracles/` - Reference values from NCTSSOS Python for validation

### Optimization Pipeline
```
polyopt() → cs_nctssos() → compute_sparsity() → moment/sos relaxation → JuMP model
```

## Key Patterns

1. **Immutable normal forms, mutable containers**: `NormalMonomial` is immutable; `Polynomial` is mutable
2. **Algebra-specific simplification**: Each algebra in `simplification/` defines `simplify(::Type{A}, word)`
3. **Registry for symbols**: Variables created via `create_*_variables()` are stored in a registry
4. **Oracle validation**: Problem tests compare against NCTSSOS Python reference implementations

## Solvers

- **CI**: COSMO (free, slower)
- **Local**: Mosek (requires license, faster)
