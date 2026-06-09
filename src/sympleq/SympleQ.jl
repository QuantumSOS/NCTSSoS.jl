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
include("graph.jl")
include("cycles.jl")
include("automorphism.jl")
include("symplectic.jl")
include("phase.jl")
include("bridge.jl")
