isdefined(@__MODULE__, :SOLVER) || include("../../TestUtils.jl")

module PublicRelaxationChecks
using Test, TOML
using NCTSSoS.Next
import ..SOLVER

include("particle_number.jl")
include("sparsity.jl")
include("pipeline_legacy.jl")

function check(;solver=JuMPBackend(SOLVER),sos_solver=nothing,references=false)
    @testset "Public relaxation contracts" begin
        particle_number(solver)
        sparsity(solver)
        PipelineLegacyChecks.check(;solver,sos_solver,references)
    end
end
end

PublicRelaxationChecks.check()
