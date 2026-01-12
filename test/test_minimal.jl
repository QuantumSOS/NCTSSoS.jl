# =============================================================================
# test/test_minimal.jl
# =============================================================================
# Minimal test suite covering critical algorithm paths.
# Run: make test-minimal (~25s COSMO, ~10s Mosek)
#
# Coverage:
# ---------
# | # | Test                  | Algorithm Path              | Source File          |
# |---|-----------------------|-----------------------------|----------------------|
# | 1 | Dense baseline        | No sparsity (baseline)      | chsh_simple.jl       |
# | 2 | Correlative sparsity  | MF clique decomposition     | chsh_simple.jl       |
# | 3 | Term sparsity         | MMD block reduction         | nc_example1.jl       |
# | 4 | Combined CS+TS        | Both + constraints          | nc_example2.jl       |
# | 5 | Dualization           | SOS ≈ Moment equivalence    | dualization.jl       |
#
# Why these tests?
# ----------------
# - Cover all major algorithm branches without redundancy
# - Higher orders use same algorithm, just larger matrices
# - Physics tests use same solver path, just different algebras
# =============================================================================

using Test

@testset "Minimal Suite" begin
    # Bell inequality tests (Dense, CS, TS at order=1)
    include("problems/bell_inequalities/chsh_simple.jl")

    # NC polynomial tests (TS, constrained)
    include("problems/nc_polynomial/nc_example1.jl")
    include("problems/nc_polynomial/nc_example2.jl")

    # Dualization equivalence
    include("relaxations/dualization.jl")
end
