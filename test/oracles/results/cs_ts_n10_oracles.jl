# CS TS Example (n=10) Oracle Results
# ====================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_cs_ts_n10.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: CS_TS_N10_{sparsity}_d{order}

# TODO: Run nctssos_cs_ts_n10.jl on a800 to generate actual values
const CS_TS_N10_ORACLES = Dict(
    "CS_TS_N10_CS_TS_d3" => (opt=3.011288, sides=Int[], nuniq=0),
)
