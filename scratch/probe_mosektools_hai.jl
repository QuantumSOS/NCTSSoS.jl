#!/usr/bin/env julia

using Pkg
const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)

try
    @eval import MosekTools
    println("MosekTools import: ok")
catch err
    println("MosekTools import failed; installing MosekTools...")
    showerror(stdout, err)
    println()
    Pkg.add(PackageSpec(name="MosekTools", version="0.15.9"))
    @eval import MosekTools
    println("MosekTools import after install: ok")
end

using JuMP
model = Model(MosekTools.Optimizer)
@variable(model, x >= 1)
@objective(model, Min, x)
optimize!(model)
println("termination_status = ", termination_status(model))
println("primal_status = ", primal_status(model))
println("objective_value = ", objective_value(model))
