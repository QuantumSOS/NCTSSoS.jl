# Pipeline layer: block family, strengthening, symmetry reduction.
# User of contract 1 and provider of contract 2 (`../contracts/conic.jl`).
# Include the files of this layer here.

# Group membership uses the core's canonical Pauli action. Other word
# families and contract-only fixture actions keep their native representation.
_group_action(::Type{W},g) where {W<:AbstractWord} = g
_group_action(::Type{W},g::SitePermutation) where {W<:PauliWord} = CliffordMap(W,length(g.images),g)
_group_action(::Type{W},g::PauliConjugation{W}) where {W<:PauliWord} = CliffordMap(g)

include("real.jl")
include("moments.jl")
include("orbits.jl")
include("finite.jl")
include("lattice.jl")
include("abelian.jl")
include("discovery.jl")
include("linear.jl")
