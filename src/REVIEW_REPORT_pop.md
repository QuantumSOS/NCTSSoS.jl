# Review Report: src/pop.jl

## Summary
The file `src/pop.jl` defines the core data structures for polynomial optimization problems (`PolyOpt` and `StatePolyOpt`). The code is clean, well-documented, and generally follows Julia best practices. However, a significant type instability was identified in the struct definitions regarding the `registry` field.

## Findings

### Medium Severity: Type Instability in Struct Definitions
- **Location**: `PolyOpt` (lines 53-58) and `StatePolyOpt` (lines 243-248).
- **Issue**: The `registry` field is typed as `VariableRegistry{A}`, which is an **abstract type** because `VariableRegistry` takes two type parameters (`A` and `T`), but only `A` is provided.
- **Impact**: Any access to `pop.registry` returns an abstractly typed object, causing dynamic dispatch and boxing, which degrades performance for operations involving variable lookups.
- **Evidence**: `verify_stability.jl` output shows `Body::VARIABLEREGISTRY{NONCOMMUTATIVEALGEBRA}` and `Warning: Field 'registry' is not concrete!`.
- **Recommendation**: Add the index type parameter `T` to the `PolyOpt` and `StatePolyOpt` struct definitions, e.g., `PolyOpt{A, T, P}` or ensure `P` and `registry` share the same `T`.

### Low Severity: Documentation Precision
- **Location**: Docstring for `StatePolyOpt` (line 238).
- **Issue**: The example uses `one(Monomial)`, but `Monomial` usually requires type parameters (e.g., `Monomial{UnipotentAlgebra, UInt8}`).
- **Recommendation**: Update example to use a specific monomial type or `one(typeof(term))` as seen in tests.

## Code Quality
- **Style**: Adheres to Blue Style.
- **Validation**: `polyopt` constructor correctly enforces that coefficient types are not integers (line 104).
- **Constraints**: Constraints are properly deduplicated using `unique!`.

## Testing
- **Coverage**: Thoroughly covered by `test/pop.jl` and `test/state_poly_opt.jl`.
- **Verification**: `polyopt` constructor itself is type stable, but the resulting struct fields are not fully concrete.
