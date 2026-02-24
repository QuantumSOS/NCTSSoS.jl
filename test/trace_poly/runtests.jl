using Test, NCTSSoS
import NCTSSoS: tr, simplify

@testset "Trace Polynomial" begin
    include("word_api.jl")
    include("polynomial_api.jl")
    include("basis.jl")
    include("chsh_trace.jl")
    include("optimization_paths.jl")
    include("coverage_edges.jl")
end
