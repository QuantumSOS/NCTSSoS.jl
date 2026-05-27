#!/usr/bin/env julia
include(joinpath(@__DIR__, "driver.jl"))
main(vcat(["--system=h2", "--norb=4", "--nelec-per-cell=2"], ARGS))
