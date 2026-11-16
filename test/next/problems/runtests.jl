# Also runs on its own in an environment with the test dependencies installed:
# julia --project=<test-env> test/next/problems/runtests.jl
isdefined(@__MODULE__, :SOLVER) || include("../../TestUtils.jl")

module NextProblemChecks
using Test, LinearAlgebra, TOML
using NCTSSoS
import Clarabel, JuMP
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
include("fermionic/fermionic.jl")
include("fermionic/fermionic_chain.jl")
include("fermionic/xy_model.jl")
include("condensed_matter/bose_hubbard.jl")
include("benchmarks.jl")
include("condensed_matter/heisenberg.jl")
include("fermionic_legacy.jl")
include("bell_inequalities/higher_order.jl")
include("trace_literature.jl")
include("local_references.jl")

# SOS checks are opt-in: call check(; sos_solver=MosekBackend(form=:sos))
# after loading Mosek and checking its licence. The CI baseline needs no licence.
function check(; solver=JuMPBackend(SOLVER), sos_solver=nothing,
    benchmark_solver=BENCHMARK_SOLVER,references=false)
    @testset "NCTSSoS: reviewed problem fixtures" begin
        nc_example1(solver, sos_solver)
        nc_example2(solver, sos_solver)
        chsh_simple(solver)
        pironio_toy(solver, sos_solver)
        ising(solver)
        free_fermion(solver)
        harmonic_oscillator(solver)
        fermionic_parity(solver)
        kitaev_chain(solver)
        xx_field(solver)
        xy_model(solver)
        bose_hubbard(solver)
        benchmark_contracts(benchmark_solver, sos_solver)
        heisenberg_star(solver,sos_solver;references)
        heisenberg_spin(solver)
        heisenberg_structure()
        FermionicLegacyChecks.check(;solver,sos_solver,references)
        higher_bell(sos_solver)
        trace_literature(sos_solver)
        LocalReferenceChecks.check(;references,mosek=sos_solver isa MosekBackend)
    end
end
end

NextProblemChecks.check()
