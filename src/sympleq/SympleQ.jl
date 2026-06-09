# =============================================================================
# SympleQ find-side bridge for Pauli Hamiltonian symmetries
# =============================================================================
#
# This subsystem implements the public, reproducible half of SympleQ:
# tableau construction, anticommutation/cycle graph construction,
# colour-preserving automorphisms, symplectic synthesis, best-effort phase
# verification, and conversion of supported binary Clifford actions into the
# existing `SymmetrySpec` seam via `CliffordSymmetry`.
#
# The exploit side from the SympleQ papers is intentionally not here.

include("tableau.jl")
include("graph.jl")
include("cycles.jl")
include("automorphism.jl")
include("symplectic.jl")
include("phase.jl")
include("bridge.jl")
