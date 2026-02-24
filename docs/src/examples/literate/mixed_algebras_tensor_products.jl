# # [Mixed Systems via Tensor Products](@id mixed-algebras-tensor-products)
#
# `NCTSSoS.jl` supports multiple algebra types (Pauli, fermionic, bosonic, …).
# For **mixed systems**, you can represent tensor-product basis elements using
# [`ComposedMonomial`](@ref), which stores one monomial per factor algebra.
#
# A mixed-system **term** is typically represented as a `(coefficient, ComposedMonomial)` pair.

using NCTSSoS

# Helper: pretty-print any object using a registry for symbolic names.
repr_with_registry(reg, x) = sprint(show, x; context=:registry => reg)

# ## 1) Build factor-algebra monomials

reg_p, (σx, σy, σz) = create_pauli_variables(1:1)
reg_f, (a, a⁺) = create_fermionic_variables(1:1)

# A Pauli product on the same site can introduce a complex phase:
# σx₁ σy₁ = i σz₁
coef_p, mono_p = collect(σx[1] * σy[1])[1]

# A fermionic number operator is already normal-ordered:
mono_f = monomials(a⁺[1] * a[1])[1]

repr_with_registry(reg_p, mono_p), repr_with_registry(reg_f, mono_f)

# ## 2) Form the tensor product

cm = ComposedMonomial((mono_p, mono_f))

# Carry the phase as an external coefficient:
mixed_term = (coef_p, cm)
mixed_term

# `simplify` on a `ComposedMonomial` simplifies each component and returns a list
# of `(coefficient, ComposedMonomial)` pairs (one entry for each expansion term).
#
# Here, both components are already canonical, so the result is a single unit term:
simplify(cm)

