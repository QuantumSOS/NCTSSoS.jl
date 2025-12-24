# AGENTS.md - NCTSSoS.jl Agentic Coding Guide

## Build & Test Commands
```bash
make init                                    # Precompile dependencies
LOCAL_TESTING=true make test                 # Run all tests (uses efficient solver)
make test-FastPoly                           # Run FastPolynomials tests only
julia --project -e 'include("test/pop.jl")' # Run single test file
ssh a800 "cd /home/yushengzhao/NCTSSoS.jl && LOCAL_TESTING=true make test"  # Remote test
```

## Code Style (Blue Style + Julia Idioms)
- **Formatting**: Use `JuliaFormatter` with `style = "blue"`
- **Imports**: Group stdlib → external packages → local modules; use `using X: a, b` for selective imports
- **Types**: Prefer parametric types (`Monomial{A,T}`), use type annotations in function signatures
- **Naming**: `snake_case` for functions/variables, `PascalCase` for types, prefix internal helpers with `_`
- **Functions**: Keep functions short (<20 lines), single responsibility, use multiple dispatch
- **Error Handling**: Use `@assert` for invariants, throw specific exceptions with context
- **Tests**: Use `@testset` blocks, test edge cases, group related tests hierarchically
- **Comments**: Minimal—code should be self-documenting; use docstrings for public API

## Architecture Notes
- Two-layer: `FastPolynomials` (polynomial library) + `NCTSSoS` (optimization framework)
- Six algebra types dispatch via `AlgebraType` parameter: `NonCommutative`, `Pauli`, `Fermionic`, `Bosonic`, `Projector`, `Unipotent`
- Polynomial invariants: sorted terms, no duplicates, no zero coefficients
