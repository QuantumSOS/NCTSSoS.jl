# Polynomials

The modeling API uses a registry to keep operators in one algebra. Build all
polynomials of a problem with that registry. Different registries cannot be
mixed in one problem.

## Algebras and variables

`pauli_variables`, `fermionic_variables`, `bosonic_variables` and
`site_variables` return a registry and operator arrays. Use `NCTSSoS.mixed_variables`
for commuting components with different algebras in one registry.

```@docs
NCTSSoS.Pauli
NCTSSoS.AntiunitaryPauli
NCTSSoS.Bosonic
NCTSSoS.Fermionic
NCTSSoS.Free
NCTSSoS.Projector
NCTSSoS.Unipotent
NCTSSoS.SiteAlgebra
NCTSSoS.pauli_variables
NCTSSoS.fermionic_variables
NCTSSoS.bosonic_variables
NCTSSoS.site_variables
NCTSSoS.mixed_variables
```

## State and trace symbols

```@docs
NCTSSoS.RealState
NCTSSoS.RealTrace
NCTSSoS.state_symbol
NCTSSoS.trace_symbol
```

## Construction and inspection

Use these public calls rather than reading packed storage fields.

```@docs
NCTSSoS.variables
NCTSSoS.polynomial
NCTSSoS.terms
NCTSSoS.word_letters
NCTSSoS.support_letters
NCTSSoS.sites
NCTSSoS.degree
NCTSSoS.operator_name
```

## Basis selection

```@docs
NCTSSoS.basis
NCTSSoS.newton_chip_basis
```
