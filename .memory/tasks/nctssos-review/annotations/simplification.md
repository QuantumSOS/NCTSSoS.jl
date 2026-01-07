# Simplification Module Annotations

## Files
- `src/simplification/noncommutative.jl` - No rules (site sort only)
- `src/simplification/pauli.jl` - σ²=I, phases
- `src/simplification/fermionic.jl` - Wick's theorem
- `src/simplification/bosonic.jl` - Commutation expansion
- `src/simplification/projector.jl` - P²=P
- `src/simplification/unipotent.jl` - U²=I
- `src/simplification/site_helpers.jl` - Utility functions

---

## Dispatch Architecture

Entry point: `simplify(p::Polynomial)` in `polynomial.jl:1203`

Return types by algebra:
- `Monomial` - NonCommutative, Projector, Unipotent
- `Term` - Pauli (with phase coefficient)
- `Polynomial` - Fermionic, Bosonic (Wick expansion)

---

## pauli.jl

*(To be filled during review)*

---

## fermionic.jl

*(To be filled during review)*

---

## bosonic.jl

*(To be filled during review)*
