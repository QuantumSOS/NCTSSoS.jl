# Pipeline layer: block family, strengthening, symmetry reduction.
# User of contract 1 and provider of contract 2 (`../contracts/conic.jl`).
# Include the files of this layer here.

# Group membership uses one action representation for each word family, so
# whole-site and letter/axis maps share composition, equality and hashing.
_group_action(::Type{W},g) where {W<:AbstractWord} = g
_group_action(::Type{W},g::SitePermutation) where {W<:PauliWord} = CliffordMap(W,length(g.images),g)
_group_action(::Type{W},g::PauliConjugation{W}) where {W<:PauliWord} = CliffordMap(g)
_group_action(::Type{W},g::SitePermutation) where {W<:Union{SiteWord,FunctionalWord}} = LetterPermutation(W,g)

include("real.jl")
include("moments.jl")
include("orbits.jl")
include("finite.jl")
include("lattice.jl")
include("abelian.jl")
include("discovery.jl")
include("linear.jl")
include("families.jl")
include("selection.jl")
include("coefficient_sectors.jl")
include("strengthening.jl")
include("hamiltonian.jl")
include("sectors.jl")
include("density_face.jl")
include("fermionic.jl")
include("continuous.jl")
include("continuous_blocks.jl")
include("sector_elimination.jl")
include("description.jl")
include("display.jl")
