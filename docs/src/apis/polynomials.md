# Fast Polynomials

## Algebra Types

```@autodocs
Modules = [NCTSSoS]
Pages = ["types/algebra.jl"]
```

### Categories

`NCTSSoS.jl` groups algebras by the **normal form of monomials** (i.e. what `simplify(m::Monomial)` returns):

- `MonoidAlgebra`: normal form stays a single monomial (monoid ring viewpoint). [@monoidRing]
- `TwistedGroupAlgebra`: normal form is a scalar/phase times a monomial (twisted group algebra). [@twistedGroupAlgebra]
- `PBWAlgebra`: normal form can expand into a sum of monomials (PBW-type rewriting/normal ordering). [@pbwAlgebraOscar]

## Variables

```@autodocs
Modules = [NCTSSoS]
Pages = ["types/registry.jl"]
```

## Monomials

```@autodocs
Modules = [NCTSSoS]
Pages = ["types/monomial.jl"]
```

## Polynomials

```@autodocs
Modules = [NCTSSoS]
Pages = ["types/polynomial.jl"]
```

## State Polynomial

```@autodocs
Modules = [NCTSSoS]
Pages = ["states/polynomial.jl", "states/word.jl"]
```

## Simplification Interface

```@autodocs
Modules = [NCTSSoS]
Pages = ["simplification/pauli.jl", "simplification/fermionic.jl", "simplification/bosonic.jl", "simplification/projector.jl", "simplification/unipotent.jl", "simplification/noncommutative.jl"]
```

## Utilities

```@autodocs
Modules = [NCTSSoS]
Pages = ["utils.jl"]
```
