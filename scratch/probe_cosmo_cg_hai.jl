#!/usr/bin/env julia
using Pkg
const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)

function ensure_pkg(name)
    try
        Base.require(Main, Symbol(name))
        println(name, " import: ok")
    catch err
        println(name, " import failed; installing...")
        showerror(stdout, err); println()
        Pkg.add(name)
        Base.require(Main, Symbol(name))
        println(name, " import after install: ok")
    end
end

ensure_pkg("IterativeSolvers")
ensure_pkg("LinearMaps")
import COSMO
println("CGIndirectKKTSolver defined = ", isdefined(COSMO, :CGIndirectKKTSolver))
isdefined(COSMO, :CGIndirectKKTSolver) || error("COSMO.CGIndirectKKTSolver is still unavailable")

using JuMP
model = Model(optimizer_with_attributes(
    COSMO.Optimizer,
    "verbose" => true,
    "kkt_solver" => COSMO.CGIndirectKKTSolver,
    "max_iter" => 50,
))
@variable(model, x >= 1)
@objective(model, Min, x)
optimize!(model)
println("termination_status = ", termination_status(model))
println("objective_value = ", objective_value(model))
