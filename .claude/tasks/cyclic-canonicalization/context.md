# Cyclic Canonicalization for MaxEntangled StateWords

## Task Status: COMPLETE

## Summary

Fixed cyclic-symmetric canonicalization for `MaxEntangled` (trace) StateWords. Trace polynomials require cyclic canonicalization because `tr(ABC) = tr(BCA) = tr(CAB)`.

## Problem

For `StateWord{MaxEntangled}`, monomials inside `tr()` were only being canonicalized using involution (symmetric) canonicalization: `min(m, adjoint(m))`. This meant that cyclically equivalent monomials like `tr(x1*x2*x3)` and `tr(x2*x3*x1)` were treated as different.

### Example of the Bug

```julia
m1 = x[1] * x[2] * x[3]  # word: [1, 2, 3]
m2 = x[2] * x[3] * x[1]  # word: [2, 3, 1]

sw1 = tr(m1)  # Should equal sw2
sw2 = tr(m2)

# Before fix: sw1 != sw2 (BUG)
# After fix:  sw1 == sw2 (CORRECT)
```

## Solution

Added state-type-specific canonicalization via `_state_canon()` function:

```julia
function _state_canon(::Type{Arbitrary}, m::Monomial{A,T})
    # For arbitrary states: <M> = <M†>, symmetric only
    _involution_canon(m)
end

function _state_canon(::Type{MaxEntangled}, m::Monomial{A,T})
    # For trace states: tr(ABC) = tr(BCA) = tr(C†B†A†)
    cyclic_symmetric_canon(m)
end
```

The `StateWord{ST}` constructor now calls `_state_canon(ST, m)` instead of `_involution_canon(m)`.

## Files Modified

- `src/FastPolynomials/src/state_word.jl`:
  - Modified `StateWord{ST}` constructor to use `_state_canon(ST, m)`
  - Added `_state_canon` function with dispatches for `Arbitrary` and `MaxEntangled`

- `test/fastpoly_test/state_word.jl`:
  - Added "Cyclic Canonicalization for Trace" test section
  - Added "Involution (not Cyclic) Canonicalization for Arbitrary" test section

- `test/trace_poly_opt.jl`:
  - Updated comments to accurately describe remaining test limitations

## Remaining Limitations (Unrelated to Cyclic Canonicalization)

### Example 6.1 (ProjectorAlgebra)

The `StatePolyOpt` solver doesn't apply projector-specific simplification rules (`P² = P`) the same way `main` branch's `PolyOpt` with `is_projective=true` does. This results in:
- Different SDP structure (108 vs 81 constraints)
- Different objective value (-0.17 vs -0.047)

**TODO**: Add projector simplification to StatePolyOpt solver.

### Example 6.2.1 (Order Requirement)

At order=2, the relaxation gives -8.0 (not tight). Order=3 gives the tight bound of -4.0. The test uses order=2, so it's correctly skipped.

## Test Results

All tests pass:
- FastPolynomials: 1275 tests pass
- trace_poly_opt.jl: Example 6.2.0 and 6.2.2 pass; 6.1 and 6.2.1 correctly skipped
- state_poly_opt.jl: All tests pass

## Key Insight: Main Branch Approach

On `main` branch, the `get_state_basis()` function in `src/FastPolynomials/src/simplify.jl` (line 236) already handled this:

```julia
canon_algo = ST == MaxEntangled ? cyclic_canonicalize : symmetric_canonicalize
```

The current branch's refactored type system needed this same logic applied at the `StateWord` constructor level.
