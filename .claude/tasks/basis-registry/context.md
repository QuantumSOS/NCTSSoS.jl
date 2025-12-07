# Task: basis-registry

## Request
Rewrite `basis.jl` to use `VariableRegistry` for type consistency and correct variable indices.

## Current Branch
`fastpoly-basis-registry` (branched from `fastpolynomial-redesign`)

## STATUS: COMPLETE ✅

## Key Decisions
1. **New API:** `get_ncbasis(registry::VariableRegistry{A,T}, d)` - registry-aware
2. **VariableRegistry{A,T}:** Add algebra type A as type parameter
3. **Return type:** `Vector{Polynomial}` - each element is the simplified form of one input monomial
4. **1-to-1 mapping preserved:** n^d input words → n^d output polynomials

## Final Implementation

### API
```julia
get_ncbasis(registry::VariableRegistry{A,T}, d::Int) -> Vector{Polynomial{A,T,ComplexF64}}
get_ncbasis_deg(registry::VariableRegistry{A,T}, d::Int) -> Vector{Polynomial{A,T,ComplexF64}}
```

### Behavior by Algebra Type
| Algebra | Input Word | Output Polynomial |
|---------|-----------|-------------------|
| NonCommutative | `x₁x₂` | `1.0 * x₁x₂` (1 term) |
| Pauli | `σ₁ˣσ₁ʸ` | `i * σ₁ᶻ` (1 term) |
| Fermionic | `a₁a₁†` | `1 - a₁†a₁` (2 terms) |
| Fermionic | `a₁a₁` | `0` (nilpotent) |
| Bosonic | `c₁c₁†` | `1 + c₁†c₁` (2 terms) |

## Commits
- `f8949af` - feat(fastpoly): add algebra type parameter to VariableRegistry
- `a1c4acf` - feat(fastpoly): update create_*_variables for VariableRegistry{A,T}
- `8587502` - feat(fastpoly): add registry-aware get_ncbasis and get_ncbasis_deg
- `bb2bf62` - test(fastpoly): add tests for registry-based basis API
- `319e10a` - refactor(fastpoly): migrate tests to registry-based basis API
- `effee37` - refactor(fastpoly): change get_ncbasis_* to return Vector{Polynomial}

## Handoff Summary
- **Completed**: Full registry-based basis generation with Vector{Polynomial} output
- **Key Design**: Each polynomial = simplified form of one input monomial (1-to-1 mapping)
- **Tests**: All FastPolynomials tests passing
- **Deferred**: External code updates (sparse.jl, gns.jl, etc.) - separate task
