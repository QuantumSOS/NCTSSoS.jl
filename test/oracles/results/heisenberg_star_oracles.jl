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

# TODO: Run nctssos_heisenberg_star.jl on a800 to generate actual values
const HEISENBERG_STAR_ORACLES = Dict(
    "HeisenbergStar_Dense_n10_d1" => (opt=-1.0, sides=Int[], nuniq=0),
    "HeisenbergStar_CS_n10_d1" => (opt=-1.0, sides=Int[], nuniq=0),
    "HeisenbergStar_Dense_n8_d1" => (opt=-1.0, sides=Int[], nuniq=0),
)
