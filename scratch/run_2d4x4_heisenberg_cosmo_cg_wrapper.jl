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
        showerror(stdout, err)
        println()
        Pkg.add(name)
        Base.require(Main, Symbol(name))
        println(name, " import after install: ok")
    end
end

ensure_pkg("IterativeSolvers")
ensure_pkg("LinearMaps")
ENV["SOLVER"] = "cosmo_cg"
include(joinpath(@__DIR__, "run_2d4x4_heisenberg_clifford_ed.jl"))
