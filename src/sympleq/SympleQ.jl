# =============================================================================
# SympleQ find-side bridge for Pauli Hamiltonian symmetries
# =============================================================================
#
# This subsystem implements the public, reproducible half of SympleQ:
# tableau construction, anticommutation/cycle graph construction,
# colour-preserving automorphisms, symplectic synthesis, best-effort phase
# phase solving/verification, and conversion of accepted Clifford actions into
# the existing `SymmetrySpec` seam.
#
# The exploit side from the SympleQ papers is intentionally not here.

include("tableau.jl")
include("../next/sympleq/graph.jl")
include("../next/sympleq/cycles.jl")
include("../next/sympleq/automorphism.jl")
include("../next/sympleq/symplectic.jl")
include("../next/sympleq/phase.jl")
include("bridge.jl")
