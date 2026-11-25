# julia --project prototype/pauli/test_square_seed.jl
module SquareSeedChecks

import NCTSSoS as NC
using Test
include("PauliWorkflowRedesign.jl")
using .PauliWorkflowRedesign
const WF = PauliWorkflowRedesign

word(p) = last(only(NC.terms(p)))

@testset "Domain-complete square stationarity and shared RDM coordinates" begin
    _, (x,y,z) = NC.pauli_variables(4)
    edges = [(1,2),(1,3),(2,4),(3,4)]
    H = (3//8)*sum(a[i]*a[j] for a in (x,y,z) for (i,j) in edges) + 2
    W = typeof(word(x[1]))
    for rdm in (false,true)
        seed = build_square_xxx_seed(H,2; rdm, stability=true)
        expected_rows = [one(W); word.([x;y;z]);
            [word(a[i]*b[j]) for (i,j) in edges for a in (x,y,z) for b in (x,y,z)]]
        @test seed.rows == sort!(expected_rows)
        words = NC.basis_rows(W,NC.Pauli(),1:4,4)
        domain = rdm ? Set(words) : Set(last(NC.mul(u,v)) for u in seed.rows for v in seed.rows)
        @test seed.domain == domain
        exhaustive = Dict(w=>im*(H*u-u*H) for w in words for u in (WF._wf_word(H,w),))
        filter!(kv -> !iszero(last(kv)) && all(t -> last(t) in domain,NC.terms(last(kv))),exhaustive)
        @test allunique(seed.stationarity_words)
        @test Dict(zip(seed.stationarity_words,seed.strengthening.equalities)) == exhaustive
        @test any(w -> w ∉ seed.rows,seed.stationarity_words)
        @test isempty(seed.strengthening.rows)
        @test seed.plaquettes == (rdm ? [[1,2,3,4]] : Vector{Int}[])

        # The complex forms of the Gram and density cones must use the same
        # variables. Stationarity and stability must not add any variables.
        block = NC.RowBlock(seed.rows,[(1.0+0im,one(W))],NC.HPSD,0)
        objective = [(ComplexF64(c),w) for (c,w) in NC.terms(H)]
        program,table,_ = NC.family_program([block],objective)
        gram_side, old_nvars = only(filter(b -> b.kind==NC.HPSD,program.blocks)).side, program.nvars
        extended = WF._append_strengthening!(program,table,seed.strengthening)
        @test [b.side for b in extended.blocks if b.kind==NC.HPSD] == (rdm ? [gram_side,16] : [gram_side])
        @test old_nvars == extended.nvars == length(domain)
        u = x[1]
        @test first(seed.strengthening.inequalities) == u*H*u-H
        @test length(seed.strengthening.inequalities) == 12
    end
end

@testset "Square torus seed and extensible motif orbits" begin
    L = 4
    _, (x,y,z) = NC.pauli_variables(L^2)
    east = [i % L == 0 ? i-L+1 : i+1 for i in 1:L^2]
    north = mod1.((1:L^2) .+ L,L^2)
    H = sum(a[i]*(a[east[i]]+a[north[i]]) for a in (x,y,z) for i in 1:L^2)
    seed = build_square_xxx_seed(H,L)
    @test length(seed.rows) == 337
    @test length(seed.plaquettes) == 16
    @test [1,4,13,16] in seed.plaquettes
    @test word(y[4]*z[1]) in seed.rows
    @test word(z[1]*z[6]) ∉ seed.rows
    @test all(p -> all(t -> last(t) in seed.domain,NC.terms(p)),seed.strengthening.equalities)
    @test isempty(seed.strengthening.inequalities)
    @test word(x[1]*y[6]) in seed.stationarity_words
    u = x[1]*x[2]*x[6]
    covered = [w in seed.domain for (_,w) in NC.terms(H*u-u*H)]
    @test any(covered) && !all(covered)
    @test word(u) ∉ seed.stationarity_words

    corner = ((0,0),(1,0),(0,1))
    added = build_square_xxx_seed(H,L; motifs=(corner,((1,1),(1,2),(2,1))), rdm=false)
    @test length(added.rows) == 337 + 4*L^2*27 # full D4 orbit, no duplicate added motif
    @test issubset(seed.rows,added.rows)
    @test issubset(seed.domain,added.domain)
    @test issubset(seed.stationarity_words,added.stationarity_words)
    @test word(y[4]*x[1]*y[13]) in added.rows # reflected corner crossing both seams
    @test word(x[1]*x[2]*x[3]) ∉ added.rows
    @test isempty(added.strengthening.psd_blocks)
    @test_throws ArgumentError build_square_xxx_seed(H,L; motifs=(((0,0),(1,1)),))
    @test_throws ArgumentError build_square_xxx_seed(H+x[1]*x[2],L)
end

end
