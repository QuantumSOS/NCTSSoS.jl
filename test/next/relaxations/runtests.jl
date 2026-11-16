isdefined(@__MODULE__, :SOLVER) || include("../../TestUtils.jl")

module PublicRelaxationChecks
using Test, TOML
using NCTSSoS.Next
import ..SOLVER

include("particle_number.jl")
include("sparsity.jl")

function check(;solver=JuMPBackend(SOLVER))
    @testset "Public relaxation contracts" begin
        particle_number(solver)
        sparsity(solver)
    end
end
end

PublicRelaxationChecks.check()
