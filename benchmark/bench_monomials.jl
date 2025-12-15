module BenchMonomials
using BenchmarkTools
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: monomial, monomials, neat_dot, star, SimplifyAlgorithm

const SUITE = BenchmarkGroup()

reg, (x,) = create_noncommutative_variables([("x", 1:10)])
var_vec = [x[1], x[2], x[2], x[1], x[3]]
z_vec = [10, 20, 2, 0, 3]

SUITE["Monomial Creation"] = @benchmarkable monomial(var_vec, z_vec)

SUITE["Basis Creation"] = @benchmarkable monomials(x, Val(3))

SUITE["Compare different degree"] = @benchmarkable cmp(a, b) setup = (a = monomial([x[1], x[2]], [2, 1]); b = monomial([x[1], x[2]], [1, 1])) # TODO: each 1.5ns

SUITE["Compare same degree"] = @benchmarkable cmp(a, b) setup = (a = monomial([x[1], x[2], x[1]], [2, 1, 1]); b = monomial([x[1], x[2]], [3, 1])) # TODO: reach 2.38ns

SUITE["Neat dot"] = @benchmarkable neat_dot(a, b) setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [8, 4, 3, 2, 9, 8, 10, 6, 8, 9]))

SUITE["Neat dot small total degree"] = @benchmarkable neat_dot(a, b) setup = (
    a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [1, 2, 1, 1, 2, 1, 2, 1, 1, 2]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [2, 1, 1, 2, 1, 2, 2, 2, 1, 1])
)

SUITE["Star"] = @benchmarkable star(a) setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]))

SUITE["Multiplication"] = @benchmarkable a * b setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [8, 5, 3, 6, 5, 10, 2, 8, 10, 7]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [8, 4, 3, 2, 9, 8, 10, 6, 8, 9]))

SUITE["Multiplication small total degree"] = @benchmarkable a * b setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [1, 1, 2, 1, 2, 2, 1, 2, 2, 1]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [2, 1, 1, 2, 1, 1, 2, 1, 2, 1]))

SUITE["Neat dot with multiplication"] = @benchmarkable neat_dot(a, b * c) setup = (a = monomial(x[[6, 8, 7, 1, 2, 5, 8, 3, 4, 2]], [1, 1, 2, 1, 2, 2, 1, 2, 2, 1]); b = monomial(x[[5, 1, 3, 7, 4, 8, 7, 6, 3, 9]], [2, 1, 1, 2, 1, 1, 2, 1, 2, 1]); c = monomial(x[[3, 4, 5, 7]], [2, 1, 2, 1]))

SUITE["_comm!"] = @benchmarkable _comm!(mono, comm_gps) setup = (
    mono = monomial([x[1], x[4], x[2], x[5], x[6]], [1, 2, 3, 4, 5, 6]);
    comm_gps = Dict(zip([x[1:3]; x[4:6]], [fill(1, 3); fill(2, 3)]))
)


SUITE["Simplification vanilla"] = @benchmarkable simplify!(a, sa_vanilla) setup = (
    a = monomial(x[[6, 3, 6, 1, 2, 5, 4]], [1, 1, 1, 2, 1, 3, 2]);
    sa_vanilla = SimplifyAlgorithm(comm_gps=[x], is_unipotent=false, is_projective=false)
)

SUITE["Simplification cm_gp"] = @benchmarkable simplify!(a, sa_cmgp) setup = (
    a = monomial(x[[6, 3, 6, 1, 2, 5, 4]], [1, 1, 1, 2, 1, 3, 2]);
    sa_cmgp = SimplifyAlgorithm(comm_gps=[x[1:3], x[4:6]], is_unipotent=false, is_projective=false)
)

SUITE["Simplification projective"] = @benchmarkable simplify!(a, sa_proj) setup = (
    a = monomial(x[[6, 3, 6, 1, 2, 5, 4]], [1, 1, 1, 2, 1, 3, 2]);
    sa_proj = SimplifyAlgorithm(comm_gps=[x[1:3], x[4:6]], is_unipotent=false, is_projective=true)
)

SUITE["Simplification unipotent"] = @benchmarkable simplify!(a, sa_uni) setup = (
    a = monomial(x[[6, 3, 6, 1, 2, 5, 4]], [1, 1, 1, 2, 1, 3, 2]);
    sa_uni = SimplifyAlgorithm(comm_gps=[x[1:3], x[4:6]], is_unipotent=true, is_projective=false)
)

SUITE["Simplification projective"] = @benchmarkable simplify!(a, sa_proj) setup = (
    a = monomial(x[[6, 3, 6, 1, 2, 5, 4]], [1, 1, 1, 2, 1, 3, 2]);
    sa_proj = SimplifyAlgorithm(comm_gps=[x[1:3], x[4:6]], is_unipotent=false, is_projective=true)
)

end

BenchMonomials.SUITE
