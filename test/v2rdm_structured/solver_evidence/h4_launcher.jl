#!/usr/bin/env julia
include(joinpath(@__DIR__, "driver.jl"))
main(vcat(["--system=h4", "--norb=8", "--nelec-per-cell=4", "--include-one-d=true"], ARGS))
