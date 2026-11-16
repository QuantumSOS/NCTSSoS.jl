# Exact action on unnormalized occupation states; no cutoff or polynomial
# rewrite is used in this oracle. Creation has coefficient 1, annihilation n.
function ladder_action(word, occupation)
    state, amplitude = collect(occupation), 1
    for op in reverse(word)
        mode = abs(op)
        if op > 0
            state[mode] == 0 && return (Tuple(state), 0)
            amplitude *= state[mode]
            state[mode] -= 1
        else
            state[mode] += 1
        end
    end
    return (Tuple(state), amplitude)
end

function simplify_oracles()
    @testset "Complete bosonic contractions" begin
        for (modes, alphabet, occupations, max_length) in
            ((1,[1,-1],((0,),(1,),(4,)),6), (2,[1,-1,2,-2],((0,2),(3,1)),4))
            registry,(b,bd) = bosonic_variables(modes)
            ops = Dict([i=>b[i] for i in 1:modes]...,[-i=>bd[i] for i in 1:modes]...)
            labels = Dict(NCTSSoS.operator_name(registry,last(only(NCTSSoS.terms(p))))=>i for (i,p) in ops)
            for len in 0:max_length, letters in Iterators.product(ntuple(_->alphabet,len)...)
                word = collect(letters)
                p = foldl(*,(ops[i] for i in word);init=one(b[1]))
                for occupation in occupations
                    state,amplitude = ladder_action(word,occupation)
                    expected = amplitude==0 ? Dict{typeof(state),ComplexF64}() : Dict(state=>ComplexF64(amplitude))
                    actual = Dict{typeof(state),ComplexF64}()
                    for (coefficient,term) in NCTSSoS.terms(p)
                        normal = [labels[NCTSSoS.operator_name(registry,l)] for l in NCTSSoS.word_letters(term)]
                        target,value = ladder_action(normal,occupation)
                        value==0 && continue
                        actual[target] = get(actual,target,0)+coefficient*value
                    end
                    @test actual == expected
                end
            end
        end
    end
    @testset "Fixed monoid reductions" begin
        # Basic squares and one cascading pair are already checked in model
        # and core tests. These longer fixed cases came from the old oracle.
        for (rule,cases) in (
            (Unipotent(), [([1,1,2,3,3,2],Int[]),([1,2,1,2],[1,2,1,2]),([1,1,1],[1])]),
            (Projector(), [([1,1,1],[1]),([1,1,2,2,1],[1,2,1]),([1,2,1],[1,2,1])]),
        )
            _,(x,) = site_variables(rule,1,3)
            for (input,expected) in cases
                @test foldl(*,(x[i] for i in input);init=one(x[1])) == foldl(*,(x[i] for i in expected);init=one(x[1]))
            end
            _,(x,y) = site_variables(rule,2,2)
            expected = rule isa Unipotent ? one(x[1]) : x[1]*y[1]
            @test y[1]*x[1]*y[1]*x[1] == expected
        end
    end
    @testset "Independent normal-word counts" begin
        # Same counts as the old NCTSSOS tables, derived without the package:
        # n^k Free words; n*(n-1)^(k-1) reduced projector/involution words.
        # Degree zero and n=2/order=2 are already in model/inspection.jl.
        for rule in (Free(),Projector(),Unipotent()), n in 1:4, order in 1:3
            n==2 && order==2 && continue
            registry,_ = site_variables(rule,1,n)
            expected = rule isa Free ? sum(n^k for k in 0:order) :
                1+sum(n*(n-1)^(k-1) for k in 1:order)
            @test length(NCTSSoS.basis(registry,order)) == expected
        end
        # One selected operator at site 1, two at site 2. Restrict the public
        # basis iterator to that alphabet; no custom/internal registry is used.
        for rule in (Projector(),Unipotent())
            registry,_ = site_variables(rule,2,2)
            selected = Set((:s1_1,:s2_1,:s2_2))
            rows = filter(collect(NCTSSoS.basis(registry,3))) do p
                all(l->NCTSSoS.operator_name(registry,l) in selected,
                    NCTSSoS.support_letters(last(only(NCTSSoS.terms(p)))))
            end
            @test length(rows) == 12
        end
    end
end
