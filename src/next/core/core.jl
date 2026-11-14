# Core layer: packed words for each algebra, polynomials, registry.
# Provider of contract 1 (`../contracts/words.jl`). Include the files of this layer here.

include("pauli.jl")
include("fermi.jl")
include("site.jl")
include("polynomial.jl")
include("registry.jl")
include("actions.jl")
include("moments.jl")
