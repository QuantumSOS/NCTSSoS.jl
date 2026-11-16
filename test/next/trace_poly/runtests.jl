isdefined(@__MODULE__, :SOLVER) || include("../../TestUtils.jl")

module PublicFunctionalChecks
using Test, LinearAlgebra, TOML
using NCTSSoS.Next
import ..SOLVER

include("basis.jl")
include("polynomial_api.jl")
include("../state_poly/polyopt_constructor.jl")
include("../state_poly/arithmetic.jl")
include("../problems/state_polynomial/survey_comparative.jl")

# As in the problem pilot, licensed checks are explicit, not environment-driven.
function check(;solver=JuMPBackend(SOLVER),sos_solver=nothing)
    @testset "Public functional contracts" begin
        basis_contracts()
        polynomial_contracts()
        state_constructor()
        state_arithmetic()
        state_chsh(sos_solver)
        survey_witnesses()
        survey_comparative(solver,sos_solver)
    end
end
end

PublicFunctionalChecks.check()
