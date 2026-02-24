```@meta
EditURL = "../literate/pbw_algebras_showcase.jl"
```

# [PBW Algebra Showcase (Fermionic / Bosonic)](@id pbw-algebras-showcase)

This page demonstrates the two `PBWAlgebra` variants supported by `NCTSSoS.jl`:
- [`FermionicAlgebra`](@ref) (CAR / anticommutation)
- [`BosonicAlgebra`](@ref) (CCR / commutation)

In PBW algebras, multiplying monomials can expand into a **sum of monomials**.
In `NCTSSoS.jl`, this means `NormalMonomial * NormalMonomial` returns a [`Polynomial`](@ref)
that may contain multiple terms.

````julia
using NCTSSoS

repr_with_registry(reg, x) = sprint(show, x; context=:registry => reg)
````

````
repr_with_registry (generic function with 1 method)
````

## 1) FermionicAlgebra (CAR)

Fermionic operators satisfy the canonical anticommutation relations (CAR):
```math
\{a_i, a_j^\dagger\} = \delta_{ij}, \qquad \{a_i, a_j\} = 0, \qquad \{a_i^\dagger, a_j^\dagger\} = 0.
```

Create annihilation operators `a[i]` and creation operators `a⁺[i]`:

````julia
reg_f, (a, a⁺) = create_fermionic_variables(1:2)
````

````
(VariableRegistry with 4 variables: a⁺₂, a⁺₁, a₁, a₂, (NCTSSoS.NormalMonomial{NCTSSoS.FermionicAlgebra, Int8}[Int8[1], Int8[2]], NCTSSoS.NormalMonomial{NCTSSoS.FermionicAlgebra, Int8}[Int8[-1], Int8[-2]]))
````

### CAR example: `a₁ a₁† = 𝟙 - a₁† a₁`

````julia
p_car = a[1] * a⁺[1]
repr_with_registry(reg_f, p_car)
````

````
"1 + -a⁺₁a₁"
````

Two terms: identity + normal-ordered number operator:

````julia
length(terms(p_car))
````

````
2
````

### Parity helper

The helper [`has_even_parity`](@ref) checks parity of a **monomial** (not a polynomial):

````julia
m_num = monomials(a⁺[1] * a[1])[1]  # monomial a₁†a₁
has_even_parity(m_num)
````

````
true
````

### Nilpotency shows up as the zero polynomial

Nilpotent products like `a₁ a₁` simplify to the zero `Polynomial`:

````julia
iszero(a[1] * a[1])
````

````
true
````

## 2) BosonicAlgebra (CCR)

Bosonic operators satisfy the canonical commutation relations (CCR):
```math
[c_i, c_j^\dagger] = \delta_{ij}, \qquad [c_i, c_j] = 0, \qquad [c_i^\dagger, c_j^\dagger] = 0.
```

Create annihilation operators `c[i]` and creation operators `c⁺[i]`:

````julia
reg_b, (c, c⁺) = create_bosonic_variables(1:1)
````

````
(VariableRegistry with 2 variables: c⁺₁, c₁, (NCTSSoS.NormalMonomial{NCTSSoS.BosonicAlgebra, Int8}[Int8[1]], NCTSSoS.NormalMonomial{NCTSSoS.BosonicAlgebra, Int8}[Int8[-1]]))
````

### CCR example: `c₁ c₁† = c₁† c₁ + 𝟙`

````julia
p_ccr = c[1] * c⁺[1]
repr_with_registry(reg_b, p_ccr)
````

````
"1 + c⁺₁c₁"
````

Again, two terms:

````julia
length(terms(p_ccr))
````

````
2
````

Bosonic operators are not nilpotent:

````julia
iszero(c[1] * c[1])
````

````
false
````

## Summary

- PBW algebras can introduce extra terms during normal ordering.
- Fermionic: sign changes + nilpotency.
- Bosonic: no sign changes, but delta correction terms.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

