module BenchPolynomials
using BenchmarkTools
using NCTSSoS.FastPolynomials

const SUITE = BenchmarkGroup()

function poly_create(x, n)
    f = 0.0
    for i = 1:n
        jset = max(1, i-5):min(n, i+1)
        jset = setdiff(jset, i)
        f += (2x[i] + 5 * x[i]^3 + 1)^2
        f -= sum([
            4x[i] * x[j] +
            10x[i]^3 * x[j] +
            2x[j] +
            4x[i] * x[j]^2 +
            10x[i]^3 * x[j]^2 +
            2x[j]^2 for j in jset
        ])
        f += sum([
            x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset
        ])
    end
end

n = 5

reg, (x,) = create_noncommutative_variables([("x", 1:n)])

SUITE["Polynomial Creation"] = @benchmarkable poly_create(x, n)

p1 = x[1]^2 + 2x[2]*x[3] + 3x[3]*x[1] + 4x[4]^100*x[5]^2

SUITE["Polynomial get variables"] = @benchmarkable variables(p1)
end

BenchPolynomials.SUITE
