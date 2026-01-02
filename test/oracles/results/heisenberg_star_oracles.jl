# Heisenberg Star Graph Oracle Results
# =====================================
# Generated from NCTSSOS (local) with MosekTools
# Run: cd ~/projects/NCTSSOS && julia --project path/to/nctssos_heisenberg_star.jl
#
# Format: (opt, sides, nuniq)
#   opt   = optimal value (minimization)
#   sides = vector of moment matrix block sizes
#   nuniq = length(data.ksupp) = unique moment indices (affine constraints)
#
# Naming: HeisenbergStar_{sparsity}_n{sites}_d{order}
#
# NOTE: Expected optimal is -1.0 for star graph at any size (tight bound)

# Generated 2026-01-02 from NCTSSOS on a800 with MosekTools
const HEISENBERG_STAR_ORACLES = Dict(
    "HeisenbergStar_Dense_n10_d1" => (opt=-0.9999999999683749, sides=[46], nuniq=1036),
    "HeisenbergStar_CS_n10_d1" => (opt=-0.999999999385143, sides=[33, 29, 29, 18, 18, 18, 18, 18], nuniq=826),
    "HeisenbergStar_Dense_n8_n8_d1" => (opt=-0.9999999963007001, sides=[29], nuniq=407),
)
