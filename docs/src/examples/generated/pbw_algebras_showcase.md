```@meta
EditURL = "../literate/pbw_algebras_showcase.jl"
```

# [PBW Algebra Showcase (Fermionic / Bosonic)](@id pbw-algebras-showcase)

This page demonstrates two PBW algebras supported by `NCTSSoS`:
- `Fermionic()` (CAR / anticommutation)
- `Bosonic()` (CCR / commutation)

In PBW algebras, multiplying monomials can expand into a **sum of monomials**.
Variables are polynomials. Their product can contain multiple normal words;
`NCTSSoS.terms` iterates the `(coefficient, word)` pairs.

````julia
using NCTSSoS

repr_with_registry(reg, x) = sprint(show, x; context=:registry => reg)
````

````
repr_with_registry (generic function with 1 method)
````

## 1) Fermionic operators (CAR)

Fermionic operators satisfy the canonical anticommutation relations (CAR):
```math
\{a_i, a_j^\dagger\} = \delta_{ij}, \qquad \{a_i, a_j\} = 0, \qquad \{a_i^\dagger, a_j^\dagger\} = 0.
```

Create annihilation operators `a[i]` and creation operators `a⁺[i]`:

````julia
reg_f, (a, a⁺) = fermionic_variables(1:2)
````

````
(Registry(Fermionic, 2 sites, 4 variables), ([a1, a2], [ad1, ad2]))
````

### CAR example: `a₁ a₁† = 𝟙 - a₁† a₁`

````julia
p_car = a[1] * a⁺[1]
repr_with_registry(reg_f, p_car)
````

````
"1.0 + (-1.0) ad1 a1"
````

Two terms: identity + normal-ordered number operator:

````julia
length(NCTSSoS.terms(p_car))
````

````
2
````

### Parity helper

A normal word has even fermionic parity when its degree is even.
Extract the word from this single-term polynomial:

````julia
m_num = last(only(NCTSSoS.terms(a⁺[1] * a[1])));
iseven(NCTSSoS.degree(m_num))
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

## 2) Bosonic operators (CCR)

Bosonic operators satisfy the canonical commutation relations (CCR):
```math
[c_i, c_j^\dagger] = \delta_{ij}, \qquad [c_i, c_j] = 0, \qquad [c_i^\dagger, c_j^\dagger] = 0.
```

Create annihilation operators `c[i]` and creation operators `c⁺[i]`:

````julia
reg_b, (c, c⁺) = bosonic_variables(1:1)
````

````
(Registry(Bosonic, 1 sites, 2 variables), ([b1], [bd1]))
````

### CCR example: `c₁ c₁† = c₁† c₁ + 𝟙`

````julia
p_ccr = c[1] * c⁺[1]
repr_with_registry(reg_b, p_ccr)
````

````
"1.0 + bd1 b1"
````

Again, two terms:

````julia
length(NCTSSoS.terms(p_ccr))
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

