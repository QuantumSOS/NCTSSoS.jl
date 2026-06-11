#!/usr/bin/env julia

# Remote wrapper: install/enable MosekTools in the active root environment, then
# run the 4x4 Heisenberg Clifford experiment with SOLVER=mosek.
# This is intentionally scratch-only; root Project.toml is not changed locally
# unless this wrapper is run locally.

using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)

function ensure_mosektools!()
    try
        @eval import MosekTools
        println("MosekTools import: ok")
        return nothing
    catch err
        println("MosekTools import failed; installing MosekTools into remote root project...")
        showerror(stdout, err)
        println()
        Pkg.add(PackageSpec(name="MosekTools", version="0.15.9"))
        @eval import MosekTools
        println("MosekTools import after install: ok")
        return nothing
    end
end

ensure_mosektools!()
include(joinpath(@__DIR__, "run_2d4x4_heisenberg_clifford_ed.jl"))
