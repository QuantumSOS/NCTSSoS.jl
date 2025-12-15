module BenchMonomials
using BenchmarkTools
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: neat_dot, star

const SUITE = BenchmarkGroup()

# Create variables using new typed API
const reg_nc, (x,) = create_noncommutative_variables([("x", 1:10)])
const reg_proj, (xp,) = create_projector_variables([("x", 1:10)])
const reg_uni, (xu,) = create_unipotent_variables([("x", 1:10)])

SUITE["Basis Creation"] = @benchmarkable get_ncbasis(x, 3)

SUITE["Neat dot"] = @benchmarkable neat_dot(a, b) setup = (
    a = x[6] * x[8] * x[7] * x[1] * x[2];
    b = x[5] * x[1] * x[3] * x[7] * x[4]
)

SUITE["Star"] = @benchmarkable star(a) setup = (a = x[6] * x[8] * x[7] * x[1] * x[2])

SUITE["Multiplication"] = @benchmarkable a * b setup = (
    a = x[6] * x[8] * x[7] * x[1] * x[2];
    b = x[5] * x[1] * x[3] * x[7] * x[4]
)

SUITE["Neat dot with multiplication"] = @benchmarkable neat_dot(a, b * c) setup = (
    a = x[6] * x[8] * x[7];
    b = x[5] * x[1] * x[3];
    c = x[3] * x[4] * x[5]
)

# Simplification benchmarks using typed algebras
SUITE["Simplification NonCommutative"] = @benchmarkable simplify(a) setup = (
    a = x[6] * x[3] * x[6] * x[1] * x[2] * x[5] * x[4]
)

SUITE["Simplification Projector"] = @benchmarkable simplify(a) setup = (
    a = xp[6] * xp[3] * xp[6] * xp[1] * xp[2] * xp[5] * xp[4]
)

SUITE["Simplification Unipotent"] = @benchmarkable simplify(a) setup = (
    a = xu[6] * xu[3] * xu[6] * xu[1] * xu[2] * xu[5] * xu[4]
)

end

BenchMonomials.SUITE
