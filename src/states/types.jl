"""
    StateType

Abstract type for quantum state types.

Subtypes determine how expectation values are displayed and interpreted:
- `Arbitrary`: General quantum state (displayed as angle brackets)
- `MaxEntangled`: Maximally entangled/trace state (displayed as "tr")

# Design
StateType is used as a type parameter in StateWord and related types to
distinguish between different state semantics at the type level.

# Examples
```jldoctest
julia> Arbitrary <: StateType
true

julia> MaxEntangled <: StateType
true
```
"""
abstract type StateType end

"""
    Arbitrary <: StateType

Represents an arbitrary quantum state (no special structure assumed).

Expectations with respect to Arbitrary states are displayed using angle brackets:
`<M1><M2>...` for a product of expectations.

# Examples
```jldoctest
julia> Arbitrary()
Arbitrary()

julia> Arbitrary() isa StateType
true
```
"""
struct Arbitrary <: StateType end

"""
    MaxEntangled <: StateType

Represents a maximally entangled state (trace state).

Used for trace optimization problems where expectations are normalized traces.
Displayed as `tr(M1)tr(M2)...` for a product of expectations.

# Examples
```jldoctest
julia> MaxEntangled()
MaxEntangled()

julia> MaxEntangled() isa StateType
true
```
"""
struct MaxEntangled <: StateType end

# Show method for clean output
Base.show(io::IO, s::StateType) = print(io, nameof(typeof(s)), "()")
