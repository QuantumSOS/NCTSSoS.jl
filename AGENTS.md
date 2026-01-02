# AGENTS.md - NCTSSoS.jl

Julia package for non-commutative polynomial optimization via SDP relaxations with sparsity exploitation.

## Build/Test Commands

```bash
make init                    # Precompile package
make test                    # Full suite with Mosek (--local)
make test-ci                 # CI suite (COSMO, no physics)
make test-polynomials        # Polynomial algebra only
make test-solvers            # SDP solver tests
make test-quality            # Code quality checks

# Single test file (REPL workflow)
julia --project -e 'using Pkg; Pkg.instantiate(); include("test/polynomials/monomials.jl")'

# Multiple test groups
julia --project -e 'using Pkg; Pkg.test(test_args=["--polynomials", "--solvers"])'
```

### Test Structure
```
test/
├── polynomials/     # Core algebra (no solver)
├── quality/         # Aqua, ExplicitImports, Doctest
├── solvers/problems/  # Problem-based tests (CHSH, I3322, Heisenberg)
├── physics/         # Physics models (--local only, requires Mosek)
├── oracles/         # NCTSSOS reference values
│   ├── problems/    # Problem definitions
│   ├── scripts/     # Oracle generation scripts (run with NCTSSOS + Mosek)
│   └── results/     # Oracle values for comparison
└── setup.jl         # Solver config (Mosek vs COSMO)
```

### Oracle Testing Workflow
1. Define problem in `test/oracles/problems/{problem}.jl`
2. Run `test/oracles/scripts/nctssos_{problem}.jl` (server with NCTSSOS + Mosek)
3. Record results in `test/oracles/results/{problem}_oracles.jl`
4. Write comparison test in `test/solvers/problems/{problem}.jl`

## Code Style

### Naming
- Types: `PascalCase` - `Polynomial`, `PolyOptResult`
- Functions: `snake_case` - `variable_indices`, `correlative_sparsity`
- Constants: `SCREAMING_SNAKE` - `CHSH_ORACLES`
- Private helpers: `_` prefix - `_process_terms`
- Mutating functions: `!` suffix - `simplify!`

### Type Parameters
- `A<:AlgebraType`: Algebra type for dispatch
- `T<:Integer`: Index type (UInt16 self-adjoint, Int32 creation/annihilation)
- `C<:Number`: Coefficient type (Float64 or ComplexF64)
- Use `where` clauses: `function foo(p::Polynomial{A,T,C}) where {A,T,C}`

### Imports
```julia
using SparseArrays, LinearAlgebra, JuMP  # Module imports
import CliqueTrees.cliquetree             # Extend specific functions
using NCTSSoS: variable_indices, expval   # Internal access in tests
```

### Type Stability
```julia
# Typed allocations, avoid Any
terms = Term{Monomial{A,T},C}[]  # GOOD
terms = []                       # BAD

# @inline for small hot functions
@inline site_bits(::Type{T}) where {T<:Unsigned} = sizeof(T) * 2
```

### Error Handling
```julia
throw(ArgumentError("Cannot compare monomials of different algebras: $A1 vs $A2"))
throw(DomainError(n, "Polynomial exponent must be non-negative"))
@assert site >= 1 && site <= ms "Site $site out of range"
```

### Core Types
```julia
Monomial{A<:AlgebraType, T<:Integer}  # Word representation
Term{M<:AbstractMonomial, C<:Number}   # Coefficient + monomial
Polynomial{A,T,C}                       # Sum of terms (sorted, deduplicated)
```

**Algebra Types** (singleton structs for dispatch)
- `NonCommutativeAlgebra`: Generic NC, no simplification
- `PauliAlgebra`: σ²=I, cyclic products (ComplexF64)
- `FermionicAlgebra`: Anticommutation, normal ordering (signed indices)
- `BosonicAlgebra`: Commutation, normal ordering (signed indices)
- `ProjectorAlgebra`: P²=P idempotent
- `UnipotentAlgebra`: U²=I involution

**Invariants**: Polynomial terms sorted by monomial, unique, non-zero coefficients

### Testing Patterns
```julia
@testset "Feature" begin
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    oracle = CHSH_ORACLES["CHSH_Dense_d1"]
    
    config = SolverConfig(optimizer=SOLVER, order=1)
    result = cs_nctssos(pop, config)
    
    @test result.objective ≈ oracle.opt atol=1e-6
    @test_throws DomainError m^(-1)
end
```

### Solver Config
```julia
# CI (COSMO)
config = SolverConfig(optimizer=COSMO.Optimizer, order=1)

# Local (Mosek)
config = SolverConfig(optimizer=Mosek.Optimizer, order=2, cs_algo=MF(), ts_algo=MMD())
```

### Quality Checks
- **Aqua.jl**: Ambiguities, unbound args, piracy
- **ExplicitImports.jl**: Import hygiene
- **Doctest**: Docstring examples

```julia
# Ignored non-public accesses
check_all_qualified_accesses_are_public(NCTSSoS, 
    ignore=(:Zeros, :PositiveSemidefiniteConeSquare, :power_by_squaring, :show_default, :HasLength))
```

## Key Files

- `src/NCTSSoS.jl`: Module definition, exports
- `src/types/`: Core types (algebra, monomial, term, polynomial)
- `src/simplification/`: Algebra-specific simplification
- `src/optimization/`: SDP relaxation, sparsity exploitation
- `test/setup.jl`: Solver configuration

## Dependencies

Core: `JuMP`, `LinearAlgebra`, `SparseArrays`, `CliqueTrees`, `ChordalGraph`, `Graphs`
Test: `Clarabel` (default), `COSMO` (CI), `MosekTools` (local)
Quality: `Aqua`, `ExplicitImports`
