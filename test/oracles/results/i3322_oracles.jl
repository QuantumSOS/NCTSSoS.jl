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

# TODO: Run nctssos_i3322.jl on a800 to generate actual values
const I3322_ORACLES = Dict(
    "I3322_Dense_d2" => (opt=-0.25, sides=Int[], nuniq=0),
    "I3322_CS_d2" => (opt=-0.25, sides=Int[], nuniq=0),
    "I3322_TS_d3" => (opt=-0.25, sides=Int[], nuniq=0),
    "I3322_CS_TS_d3" => (opt=-0.25, sides=Int[], nuniq=0),
)
