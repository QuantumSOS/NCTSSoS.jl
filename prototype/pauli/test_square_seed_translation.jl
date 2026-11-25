# julia --project prototype/pauli/test_square_seed_translation.jl [L=4]
module SquareSeedTranslationChecks

import NCTSSoS as NC
using Test
include("PauliWorkflowRedesign.jl")
using .PauliWorkflowRedesign

function translate_polynomial(p::NC.Polynomial{W,C}, shift) where {W,C}
    terms = NC.PolynomialBuffer{W,C}()
    for (c,u) in NC.terms(p)
        a,w = NC.act(shift,u)
        NC.push_term!(terms,NC.scale(a,c),w)
    end
    return NC._polynomial(p.registry,terms)
end

function translation_oracle(L)
    coords = [(a,b) for b in 0:L-1 for a in 0:L-1]
    translations = [NC.SitePermutation([mod(a+i,L)+L*mod(b+j,L)+1 for (a,b) in coords]) for (i,j) in coords]
    to_origin = [translations[mod(-a,L)+L*mod(-b,L)+1] for (a,b) in coords]
    return translations, to_origin
end

function normalized_image(w, to_origin)
    best, action = w, first(to_origin)
    found = false
    for i in NC.sites(w)
        s = to_origin[i]
        v = last(NC.act(s,w))
        if !found || v < best
            best, action = v, s
            found = true
        end
    end
    return best, action
end

function assert_plaquette_orbit(quotient, dense, translations)
    origin = only(quotient.plaquettes)
    matrix = only(quotient.strengthening.psd_blocks).matrix
    # Wrapped faces have a different sorted ket order. Check the operator
    # entries after the corresponding permutation of the four qubit factors.
    old = Dict(Tuple(sites)=>block.matrix for (sites,block) in zip(dense.plaquettes,dense.strengthening.psd_blocks))
    @test all(translations) do t
        image = [NC._image(t,s) for s in origin]
        region = sort(image)
        permutation = [findfirst(==(s),image) for s in region]
        kets = [1+sum(((j>>(4-k))&1)<<(4-i) for (i,k) in enumerate(permutation)) for j in 0:15]
        old[Tuple(region)][kets,kets] == map(p -> translate_polynomial(p,t),matrix)
    end
end

L = isempty(ARGS) ? 4 : parse(Int,only(ARGS))
n = L^2
_, axes = NC.pauli_variables(n)
east = [i % L == 0 ? i-L+1 : i+1 for i in 1:n]
north = mod1.((1:n).+L,n)
H = (-3//8)*sum(a[i]*(a[east[i]]+a[north[i]]) for a in axes for i in 1:n) + 2
translations, to_origin = translation_oracle(L)

for motifs in ((),(((0,0),(1,0),(0,1)),))
    @testset "$(L)×$(L) translation quotient (motifs=$(length(motifs)))" begin
        dense = build_square_xxx_seed(H,L; motifs,stability=true)
        quotient = build_square_xxx_seed(H,L; motifs,stability=true,translation_quotient=true)
        @test quotient.rows == dense.rows
        normalized = Set(first(normalized_image(w,to_origin)) for w in dense.domain)
        @test Set(quotient.domain) == normalized
        @test all(in(quotient.domain),dense.domain)
        @test length(quotient.domain) < length(dense.domain)
        @test length(quotient.strengthening.inequalities) == 3
        stable = quotient.strengthening.inequalities
        @test Set(translate_polynomial(p,t) for p in stable for t in translations) == Set(dense.strengthening.inequalities)
        assert_plaquette_orbit(quotient,dense,translations)

        by_source = Dict(zip(quotient.stationarity_words,quotient.strengthening.equalities))
        @test allunique(quotient.stationarity_words)
        @test Set(quotient.stationarity_words) == Set(first(normalized_image(w,to_origin)) for w in dense.stationarity_words)
        @test all(zip(dense.stationarity_words,dense.strengthening.equalities)) do (w,p)
            rep,t = normalized_image(w,to_origin)
            get(by_source,rep,nothing) == translate_polynomial(p,t)
        end

        # Membership is defined for physical words, not packed padding sites.
        W = eltype(quotient.rows)
        @test NC.pauli_letter(W,n+1,0) ∉ quotient.domain
    end
end

end
