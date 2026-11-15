# # [PBW Algebra Showcase (Fermionic / Bosonic)](@id pbw-algebras-showcase)
#
# This page demonstrates two PBW algebras supported by `NCTSSoS.Next`:
# - `Fermionic()` (CAR / anticommutation)
# - `Bosonic()` (CCR / commutation)
#
# In PBW algebras, multiplying monomials can expand into a **sum of monomials**.
# Variables are polynomials. Their product can contain multiple normal words;
# `Next.terms` iterates the `(coefficient, word)` pairs.

using NCTSSoS.Next

repr_with_registry(reg, x) = sprint(show, x; context=:registry => reg)

# ## 1) Fermionic operators (CAR)
#
# Fermionic operators satisfy the canonical anticommutation relations (CAR):
# ```math
# \{a_i, a_j^\dagger\} = \delta_{ij}, \qquad \{a_i, a_j\} = 0, \qquad \{a_i^\dagger, a_j^\dagger\} = 0.
# ```
#
# Create annihilation operators `a[i]` and creation operators `a⁺[i]`:

reg_f, (a, a⁺) = fermionic_variables(1:2)

# ### CAR example: `a₁ a₁† = 𝟙 - a₁† a₁`

p_car = a[1] * a⁺[1]
repr_with_registry(reg_f, p_car)

# Two terms: identity + normal-ordered number operator:
length(Next.terms(p_car))
@assert length(Next.terms(p_car)) == 2 #src
@assert p_car == 1 - a⁺[1] * a[1] #src

# ### Parity helper
#
# A normal word has even fermionic parity when its degree is even.
# Extract the word from this single-term polynomial:

m_num = last(only(Next.terms(a⁺[1] * a[1])));
iseven(Next.degree(m_num))
@assert iseven(Next.degree(m_num)) #src

# ### Nilpotency shows up as the zero polynomial
#
# Nilpotent products like `a₁ a₁` simplify to the zero `Polynomial`:

iszero(a[1] * a[1])
@assert iszero(a[1] * a[1]) #src

# ## 2) Bosonic operators (CCR)
#
# Bosonic operators satisfy the canonical commutation relations (CCR):
# ```math
# [c_i, c_j^\dagger] = \delta_{ij}, \qquad [c_i, c_j] = 0, \qquad [c_i^\dagger, c_j^\dagger] = 0.
# ```
#
# Create annihilation operators `c[i]` and creation operators `c⁺[i]`:

reg_b, (c, c⁺) = bosonic_variables(1:1)

# ### CCR example: `c₁ c₁† = c₁† c₁ + 𝟙`

p_ccr = c[1] * c⁺[1]
repr_with_registry(reg_b, p_ccr)

# Again, two terms:
length(Next.terms(p_ccr))
@assert length(Next.terms(p_ccr)) == 2 #src
@assert p_ccr == c⁺[1] * c[1] + 1 #src

# Bosonic operators are not nilpotent:
iszero(c[1] * c[1])
@assert !iszero(c[1] * c[1]) #src

# ## Summary
#
# - PBW algebras can introduce extra terms during normal ordering.
# - Fermionic: sign changes + nilpotency.
# - Bosonic: no sign changes, but delta correction terms.
