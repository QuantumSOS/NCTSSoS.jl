# Regenerate only the pauli_clifford_symmetry example markdown (requires Mosek).
import Pkg
Pkg.activate("docs")
using Literate

config = Dict("execute" => true, "flavor" => Literate.DocumenterFlavor())
Literate.markdown(
    "docs/src/examples/literate/pauli_clifford_symmetry.jl",
    "docs/src/examples/generated";
    config=config,
    name="pauli_clifford_symmetry",
)
println("REGEN-OK")
