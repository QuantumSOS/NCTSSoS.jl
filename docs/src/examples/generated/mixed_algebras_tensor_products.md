```@meta
EditURL = "../literate/mixed_algebras_tensor_products.jl"
```

# [Mixed Systems via Tensor Products](@id mixed-algebras-tensor-products)

`NCTSSoS.jl` supports multiple algebra types (Pauli, fermionic, bosonic, …).
For **mixed systems**, you can represent tensor-product basis elements using
[`ComposedMonomial`](@ref), which stores one monomial per factor algebra.

A mixed-system **term** is typically represented as a `(coefficient, ComposedMonomial)` pair.

````julia
using NCTSSoS
````

Helper: pretty-print any object using a registry for symbolic names.

````julia
repr_with_registry(reg, x) = sprint(show, x; context=:registry => reg)
````

````
repr_with_registry (generic function with 1 method)
````

## 1) Build factor-algebra monomials

````julia
reg_p, (σx, σy, σz) = create_pauli_variables(1:1)
reg_f, (a, a⁺) = create_fermionic_variables(1:1)
````

````
(VariableRegistry with 2 variables: a⁺₁, a₁, (NCTSSoS.NormalMonomial{NCTSSoS.FermionicAlgebra, Int8}[Int8[1]], NCTSSoS.NormalMonomial{NCTSSoS.FermionicAlgebra, Int8}[Int8[-1]]))
````

A Pauli product on the same site can introduce a complex phase:
σx₁ σy₁ = i σz₁

````julia
coef_p, mono_p = collect(σx[1] * σy[1])[1]
````

````
(0.0 + 1.0im, UInt8[0x03])
````

A fermionic number operator is already normal-ordered:

````julia
mono_f = monomials(a⁺[1] * a[1])[1]

repr_with_registry(reg_p, mono_p), repr_with_registry(reg_f, mono_f)
````

````
("σz₁", "a⁺₁a₁")
````

## 2) Form the tensor product

````julia
cm = ComposedMonomial((mono_p, mono_f))
````

````
ComposedMonomial(UInt8[0x03], Int8[-1, 1])
````

Carry the phase as an external coefficient:

````julia
mixed_term = (coef_p, cm)
mixed_term
````

````
(0.0 + 1.0im) * ComposedMonomial(UInt8[0x03], Int8[-1, 1])
````

`simplify` on a `ComposedMonomial` simplifies each component and returns a list
of `(coefficient, ComposedMonomial)` pairs (one entry for each expansion term).

Here, both components are already canonical, so the result is a single unit term:

````julia
simplify(cm)
````

````
1-element Vector{Tuple{ComplexF64, NCTSSoS.ComposedMonomial{Tuple{NCTSSoS.PauliAlgebra, NCTSSoS.FermionicAlgebra}, Tuple{NCTSSoS.NormalMonomial{NCTSSoS.PauliAlgebra, UInt8}, NCTSSoS.NormalMonomial{NCTSSoS.FermionicAlgebra, Int8}}}}}:
 ComposedMonomial(UInt8[0x03], Int8[-1, 1])
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

