# CS TS Example (n=10) Oracle Results
# ====================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_cs_ts_n10.jl
#
# Format: (opt, nblocks, nuniq)
#   opt     = optimal value (minimization)
#   nblocks = total number of moment matrix blocks
#   nuniq   = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: CS_TS_N10_{sparsity}_d{order}
#
# NOTE: Full sides array omitted (2982 blocks). Use nblocks for validation.

# Generated 2026-01-02 from NCTSSOS on a800 with MosekTools
const CS_TS_N10_ORACLES = Dict(
    "CS_TS_N10_CS_TS_d3" => (opt=3.011288353315061, nblocks=2982, nuniq=6042),
)
