"""
    NCTSSoS.Next

Work area of the redesign (packed words, one pipeline, solver standard form).
Nothing in this module is exported from `NCTSSoS`, and no name here is stable.
The current package does not use this module. At the end of the redesign this
module replaces the current code and the name `Next` goes away.

Layers and their directories:

- `contracts/`: the two interfaces between the layers. Owned by the architect.
- `core/`: packed words for each algebra, polynomials, registry.
- `pipeline/`: block family, strengthening, symmetry reduction.
- `conic/`: moment form, SOS form and solver writers.
- `model/`: the user-facing problem, relaxation and result types.
"""
module Next

using JuMP: JuMP, @constraint, @objective, @variable
import LinearAlgebra

include("contracts/words.jl")
include("contracts/conic.jl")
include("core/core.jl")
include("conic/conic.jl")
include("model/model.jl")
include("pipeline/pipeline.jl")
include("model/detection.jl")
include("precompile.jl")

end
