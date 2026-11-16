# Frozen boundary outputs of the legacy-backed detector, before its type port.
# Compare actions on every X/Y/Z letter, rather than private Clifford storage.
module DetectorReference
using Test, Random, TOML
import NCTSSoS: Next as N

function cases()
    out = Pair{String,Any}[]
    for n in (2, 3, 4), kind in (:heisenberg, :ising)
        _, (x,y,z) = N.pauli_variables(n)
        h = kind == :heisenberg ?
            sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:n-1) :
            -sum(z[i]*z[i+1] for i in 1:n-1)-0.7sum(x)
        push!(out, "$kind-open-$n" => N.Problem(h))
    end
    r = N.Registry(N.Pauli(),3,[:left=>(1,0),:right=>(1,2),:other=>(2,0),:axis=>(2,2)])
    v = N.variables(r)
    x,z,y,w = v[:left],v[:right],v[:other],v[:axis]
    push!(out, "unnamed-axes" => N.Problem(x*y+z*w-x*z*y*w))
    rng = MersenneTwister(414)
    _, ops = N.pauli_variables(4)
    h = sum(randn(rng)*ops[a][i]*ops[b][j] for i in 1:4 for j in i+1:4 for a in 1:3 for b in 1:3)+
        sum(randn(rng)*ops[a][i] for i in 1:4 for a in 1:3)
    push!(out, "asymmetric-seed414" => N.Problem(h))
    push!(out, "zero" => N.Problem(zero(h)))
    push!(out, "one" => N.Problem(one(h)))
    return out
end

function signature(problem, generators)
    W = eltype(problem.objective.data.words)
    return [join((begin
        phase, image = N.act(g,N.pauli_letter(W,i,a))
        c = N.scale(phase,1+0im)
        "$(real(c)),$(imag(c)):$(join(N.legacy_letters(image),','))"
    end for i in 1:problem.objective.registry.nsites for a in 0:2), ";") for g in generators]
end

function records()
    return [begin
        s = N.detect_symmetry(problem)
        Dict("id" => id, "expected" => Dict(
            "block" => signature(problem,s.generators),
            "moment_only" => signature(problem,s.moment_only)))
    end for (id,problem) in cases()]
end

function check(path=joinpath(@__DIR__,"../../data/expectations/next_detector.toml"))
    fixture = TOML.parsefile(path)
    @testset "Legacy-backed detector signatures" begin
        actual = records()
        @test length(actual) == length(fixture["cases"]) == 10
        for (got,expected) in zip(actual,fixture["cases"])
            @test got["id"] == expected["id"]
            @test got["expected"]["block"] == expected["expected"]["block"]
            @test got["expected"]["moment_only"] == expected["expected"]["moment_only"]
        end
    end
end
end
