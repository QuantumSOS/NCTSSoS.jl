# Polynomials

The modeling API uses a registry to keep operators in one algebra. Build all
polynomials of a problem with that registry. Different registries cannot be
mixed in one problem.

## Algebras and variables

`pauli_variables`, `fermionic_variables`, `bosonic_variables` and
`site_variables` return a registry and operator arrays. Use `Next.mixed_variables`
for commuting components with different algebras in one registry.

```@docs
NCTSSoS.Next.Pauli
NCTSSoS.Next.SiteAlgebra
NCTSSoS.Next.bosonic_variables
NCTSSoS.Next.site_variables
NCTSSoS.Next.mixed_variables
```

## State and trace symbols

```@docs
NCTSSoS.Next.RealState
NCTSSoS.Next.RealTrace
NCTSSoS.Next.state_symbol
NCTSSoS.Next.trace_symbol
```

## Construction and inspection

Use these public calls rather than reading packed storage fields.

```@docs
NCTSSoS.Next.variables
NCTSSoS.Next.polynomial
NCTSSoS.Next.terms
NCTSSoS.Next.word_letters
NCTSSoS.Next.support_letters
NCTSSoS.Next.sites
NCTSSoS.Next.degree
NCTSSoS.Next.operator_name
```

## Basis selection

```@docs
NCTSSoS.Next.basis
NCTSSoS.Next.newton_chip_basis
```
