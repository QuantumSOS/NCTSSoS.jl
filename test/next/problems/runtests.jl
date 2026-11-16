# Also runs on its own in an environment with the test dependencies installed:
# julia --project=<test-env> test/next/problems/runtests.jl
isdefined(@__MODULE__, :SOLVER) || include("../../TestUtils.jl")

module NextProblemChecks
using Test, LinearAlgebra, TOML
using NCTSSoS.Next
import ..SOLVER

function fixture(file, id)
    data = TOML.parsefile(joinpath(@__DIR__, "../../data/expectations", file * ".toml"))
    return only(filter(c -> c["id"] == id, data["cases"]))["expected"]
end

include("nc_polynomial/nc_example1.jl")
include("nc_polynomial/nc_example2.jl")
include("bell_inequalities/chsh_simple.jl")
include("bell_inequalities/pironio_toy.jl")
include("condensed_matter/ising.jl")
include("fermionic/free_fermion.jl")
include("condensed_matter/harmonic_oscillator.jl")

# SOS checks are opt-in: call check(; sos_solver=MosekBackend(form=:sos))
# after loading Mosek and checking its licence. The CI baseline needs no licence.
function check(; solver=JuMPBackend(SOLVER), sos_solver=nothing)
    @testset "Next: reviewed problem fixtures" begin
        nc_example1(solver, sos_solver)
        nc_example2(solver, sos_solver)
        chsh_simple(solver)
        pironio_toy(solver, sos_solver)
        ising(solver)
        free_fermion(solver)
        harmonic_oscillator(solver)
    end
end
end

NextProblemChecks.check()
