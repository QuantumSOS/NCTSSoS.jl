# test/relaxations/runtests.jl
# Tests: Relaxation algorithm components
#
# Coverage:
# - interface.jl: PolyOpt constructors (Polynomial)
# - sos.jl: SOS dualization components (Cαj)
# - sparsity.jl: Term sparsity + result-structure checks (non-correlated)
# - gns.jl: GNS reconstruction
# - dualization.jl: SOS ≈ Moment equivalence

using Test

@testset "Relaxations" begin
    include("interface.jl")
    include("sos.jl")
    include("sparsity.jl")
    include("symmetry.jl")
    include("gns.jl")
    include("gns_pipeline.jl")
    include("dualization.jl")
end
