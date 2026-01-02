# I_3322 Bell Inequality Oracle Results
# ======================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_i3322.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (maximizes -f)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: I3322_{sparsity}_d{order}
#   Each variant tested at its specified order from I3322_VARIANTS
#
# NOTE on CS_TS convergence:
#   CS_TS gives opt≈-1.0 and does NOT converge to the optimal -0.2509 even at order 6.
#   This is similar to CHSH - the combined sparsity loses cross-clique constraints.
#   Tested in NCTSSoS.jl: order 3-6 all give opt=-1.0 with gap=0.749.

# Generated 2026-01-02 from NCTSSOS on a800 with MosekTools
const I3322_ORACLES = Dict(
    "I3322_Dense_d2" => (opt=-0.25093979763149354, sides=[28], nuniq=154),
    "I3322_CS_d2" => (opt=-0.25907173425147956, sides=[17, 10, 10], nuniq=43),
    "I3322_TS_d3" => (opt=-0.25147090316663245, sides=[13, 12, 12, 12, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2], nuniq=121),
    "I3322_CS_TS_d3" => (opt=-0.999999989225524, sides=[17, 17, 13, 13, 13, 13, 11, 11, 11, 11, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 7, 7, 9, 9, 9, 8, 7], nuniq=43),
    "I3322_CS_TS_d4" => (opt=-0.9999999846837061, sides=[52, 49, 46, 43, 43, 43, 37, 37, 37, 37, 32, 32, 32, 32, 28, 28, 28, 28, 21, 21, 18, 18, 14, 14, 21, 21, 18, 18, 14, 14], nuniq=77),
)
