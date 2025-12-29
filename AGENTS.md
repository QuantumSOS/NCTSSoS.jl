# AGENTS.md - NCTSSoS.jl Agentic Coding Guide

## Build & Test Commands

```bash
# Local development
make init                                      # Precompile dependencies
julia --project -e 'using Pkg; Pkg.test()'     # Run all tests (requires test deps)

# Remote testing (preferred - has Mosek license)
ssh a800 "cd /home/yushengzhao/NCTSSoS.jl && LOCAL_TESTING=true make test"

# Run single test file
julia --project -e 'include("test/pop.jl")'                    # Sparsity tests
julia --project -e 'include("test/heisenberg.jl")'             # Physics model
julia --project -e 'include("test/polynomials/runtests.jl")'   # All polynomial tests
julia --project -e 'include("test/polynomials/simplify.jl")'   # Single polynomial test

# Test subsets
julia --project -e 'include("test/polynomials/runtests.jl")'   # Polynomial algebra only
LOCAL_TESTING=true julia --project -e 'include("test/sos_solver.jl")'  # SOS solver

# Documentation
make servedocs                                 # Live preview docs
make examples                                  # Generate Literate.jl examples
```

## Reference Implementation

The original NCTSSOS package is at `/home/yushengzhao/NCTSSOS` on a800.
Use it to verify results:

```julia
# Original API (NCTSSOS)
using NCTSSOS, DynamicPolynomials
@ncpolyvar x[1:2] y[1:2]
f = x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2]  # CHSH
opt, data = nctssos_first([-f], [x;y], 1, TS=false, partition=2, constraint="unipotent")
# Expected: opt ≈ -2.8284 (quantum bound = 2√2)

# New API (NCTSSoS.jl) - MINIMIZES by default
using NCTSSoS, Clarabel
reg, (U,) = create_unipotent_variables([("U", 1:4)])
CHSH = Polynomial(U[1])*Polynomial(U[3]) + Polynomial(U[1])*Polynomial(U[4]) + 
       Polynomial(U[2])*Polynomial(U[3]) - Polynomial(U[2])*Polynomial(U[4])
prob = polyopt(-CHSH, reg)  # Negate for maximization
result = cs_nctssos(prob, SolverConfig(optimizer=Clarabel.Optimizer, order=1))
# result.objective ≈ -2.8284, so max(CHSH) = -result.objective ≈ 2.8284
```

## Code Style Guidelines

### Formatting
- Use `JuliaFormatter` with `style = "blue"`
- 4-space indentation, no tabs
- Max line length: 92 characters

### Imports
```julia
# Order: stdlib → external packages → local modules
using LinearAlgebra, SparseArrays          # stdlib
using JuMP, CliqueTrees, Graphs            # external
using NCTSSoS: Monomial, Polynomial        # selective imports preferred

# Avoid: using NCTSSoS (pulls everything into namespace)
```

### Types
- Parametric types for flexibility: `Polynomial{A<:AlgebraType, T<:Integer, C<:Number}`
- Type annotations in function signatures: `function foo(x::Vector{T}) where {T<:Number}`
- Singleton types for zero-cost dispatch: `struct PauliAlgebra <: AlgebraType end`

### Naming Conventions
| Element | Convention | Example |
|---------|------------|---------|
| Functions | `snake_case` | `symmetric_canon`, `get_ncbasis` |
| Types | `PascalCase` | `Polynomial`, `VariableRegistry` |
| Constants | `SCREAMING_SNAKE` | `SUBSCRIPT_DIGITS` |
| Internal helpers | `_prefix` | `_process_terms`, `_simplify_pauli` |
| Type parameters | Single uppercase | `A`, `T`, `C` |

### Functions
- Keep under 20 lines; extract helpers
- Single responsibility
- Use multiple dispatch over conditionals on types
- Document with docstrings (triple quotes)

### Error Handling
```julia
@assert length(word) > 0 "Empty word not allowed"           # Invariants
error("Fermionic objective must have even parity")          # User errors
throw(ArgumentError("Invalid algebra type: $A"))            # Specific exceptions
```

## Architecture

### Two-Layer Design
```
┌─────────────────────────────────────────────────────────────┐
│  Optimization Layer (src/optimization/)                      │
│  - moment.jl, sos.jl: SDP relaxation                        │
│  - sparsity.jl: Correlative/term sparsity                   │
│  - interface.jl: cs_nctssos, polyopt                        │
├─────────────────────────────────────────────────────────────┤
│  Polynomial Core (src/types/, src/simplification/)          │
│  - Monomial, Term, Polynomial types                         │
│  - Algebra-specific simplification rules                    │
└─────────────────────────────────────────────────────────────┘
```

### Six Algebra Types
| Type | Rules | Index Type |
|------|-------|------------|
| `PauliAlgebra` | σ²=I, σxσy=iσz | `UInt8+` |
| `FermionicAlgebra` | {a,a†}=1, a²=0 | `Int8+` (signed) |
| `BosonicAlgebra` | [b,b†]=1 | `Int8+` (signed) |
| `ProjectorAlgebra` | P²=P | `UInt8+` |
| `UnipotentAlgebra` | U²=I | `UInt8+` |
| `NonCommutativeAlgebra` | None | `UInt8+` |

### Key APIs
```julia
# Variable creation (returns registry + monomials)
reg, (σx, σy, σz) = create_pauli_variables(1:N)
reg, (a, a†) = create_fermionic_variables(1:N)

# Polynomial construction (auto-simplifies)
H = Polynomial(σx[1]) * Polynomial(σx[2]) + ...

# Optimization (MINIMIZES objective)
prob = polyopt(objective, registry; eq_constraints=[], ineq_constraints=[])
config = SolverConfig(optimizer=Clarabel.Optimizer, order=1)
result = cs_nctssos(prob, config)
# result.objective is the minimum value
```

## Test Structure

```
test/
├── runtests.jl              # Main entry, runs all tests
├── polynomials/             # Core polynomial tests
│   ├── runtests.jl          # Polynomial test suite
│   ├── simplify.jl          # Algebra simplification
│   └── canonicalization.jl  # Symmetric/cyclic canon
├── pop.jl                   # Correlative sparsity
├── sparse.jl                # Term sparsity  
├── moment_solver.jl         # Moment SDP
├── sos_solver.jl            # SOS dualization
├── heisenberg.jl            # Physics: Heisenberg model (LOCAL_TESTING)
├── bell_ineq.jl             # Physics: Bell inequalities (LOCAL_TESTING)
└── fermionic_parity_test.jl # Fermionic superselection
```

Tests requiring Mosek run only with `LOCAL_TESTING=true` environment variable.

## Common Pitfalls

1. **Monomial vs Polynomial simplification**: `Monomial` multiplication doesn't auto-simplify.
   Use `Polynomial(m1) * Polynomial(m2)` or `simplify(m)`.

2. **Variable creation syntax**: `create_pauli_variables(1:N)` not `create_pauli_variables(N)`.

3. **Optimization direction**: `cs_nctssos` MINIMIZES. For max(f): use `-cs_nctssos(polyopt(-f, reg), config).objective`.

4. **Registry required**: `polyopt(obj, registry)` - registry is mandatory.

5. **Coefficient types**: JuMP requires `Float64` or `ComplexF64`, not integers.
