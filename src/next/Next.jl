"""
    NCTSSoS.Next

Work area of the redesign (packed words, one pipeline, solver standard form).
Nothing in this module is exported from `NCTSSoS`, and no name here is stable.
The current package does not use this module. At the end of the redesign this
module replaces the current code and the name `Next` goes away.

Load the declarative API with `using NCTSSoS.Next`. It can also be loaded with
`using NCTSSoS`; their exports do not clash. `Next.Symmetry` and `Next.variables`
are public, qualified APIs, not exports. For the common symmetry declaration,
pass a tuple or vector of actions and derivations to `Problem(H; symmetry=...)`.
Use `Next.Symmetry(...; moment_only=(...), continuous=(...))` when a unitary
action should constrain moments without closing the block family.

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
include("model/front.jl")
include("model/gns.jl")

public Symmetry, variables

export pauli_variables, fermionic_variables, site_variables,
    Pauli, Fermionic, Free, Projector, Unipotent, SiteAlgebra,
    Problem, Positive, Zero, SitePermutation, CliffordMap, AntiunitaryPauli,
    SpinRotation, spin_rotations, NumberRotation,
    detect_symmetry, Relaxation, Correlative, PhysicsInformed, TermSplit, Contiguous,
    DensityMatrix, InvariantCone, Stationarity, StateOptimality, Sector, FermionicConditions,
    inspect, solve, expectation, gns, evaluate, SCSBackend, MosekBackend, JuMPBackend

include("precompile.jl")

end
