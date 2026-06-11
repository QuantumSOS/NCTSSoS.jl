#!/usr/bin/env julia

# Remote wrapper: instantiate the root project, ensure MosekTools exists in the
# active remote environment, then run the 4x4 order-3 distance-limited Clifford
# experiment. This is scratch-only; easy-ssh pushes the local project first, so
# any Pkg.add here mutates only the remote copy.

using Pkg

const ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT)
Pkg.instantiate()

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
include(joinpath(@__DIR__, "run_2d4x4_heisenberg_order3_distance3_clifford.jl"))
