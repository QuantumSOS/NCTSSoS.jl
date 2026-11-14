# Conformance checks for contract 1 (`src/next/contracts/words.jl`).
# The tests of each word type call `check_word_contract(W, samples)`.
using Test
using NCTSSoS

const Nx = NCTSSoS.Next

"A test sink: a dictionary from word to coefficient."
struct TermDict{W,C}
    terms::Dict{W,C}
end
TermDict{W,C}() where {W,C} = TermDict{W,C}(Dict{W,C}())
function Nx.push_term!(s::TermDict{W,C}, c::Number, w::W) where {W,C}
    s.terms[w] = get(s.terms, w, zero(C)) + c
    return s
end
clean(s::TermDict) = filter(kv -> abs(kv.second) > 1e-12, s.terms)

"A test sink that counts terms and does not allocate."
mutable struct CountSink
    count::Int
end
Nx.push_term!(s::CountSink, c::Number, w) = (s.count += 1; s)

product_terms(a::W, b::W, ::Type{C}) where {W,C} = Nx.mul!(TermDict{W,C}(), a, b, one(C))
function times(p::TermDict{W,C}, b::W) where {W,C}
    out = TermDict{W,C}()
    for (w, c) in p.terms
        Nx.mul!(out, w, b, c)
    end
    return out
end
function times(a::W, p::TermDict{W,C}) where {W,C}
    out = TermDict{W,C}()
    for (w, c) in p.terms
        Nx.mul!(out, a, w, c)
    end
    return out
end
function star_terms(p::TermDict{W,C}) where {W,C}
    out = TermDict{W,C}()
    for (w, c) in p.terms
        ph, s = Nx.star(w)
        Nx.push_term!(out, Nx.scale(ph, conj(c)), s)
    end
    return out
end

function _mul_allocations(a, b, c, n)
    sink = CountSink(0)
    Nx.mul!(sink, a, b, c)
    return @allocated for _ in 1:n
        Nx.mul!(sink, a, b, c)
    end
end
function _scalar_allocations(a, b, n)
    acc = 0
    return @allocated for _ in 1:n
        acc += Nx.degree(Nx.star(a)[2]) + (isless(a, b) ? 1 : 0) + (hash(a) == hash(b) ? 1 : 0) + Nx.degree(Nx.moment_key(b)[1])
    end
end

"""
    check_word_contract(W, samples; C = ComplexF64, triples = 200)

Test that the word type `W` satisfies contract 1 on `samples`. `C` is the
coefficient type for products. Associativity is tested on the first `triples`
triples of samples.
"""
function check_word_contract(::Type{W}, samples::Vector{W}; C::Type=ComplexF64, triples::Int=200) where {W<:Nx.AbstractWord}
    @testset "word contract: $W" begin
        @test isbitstype(W)
        e = one(W)
        @test isone(e)
        @test Nx.degree(e) == 0
        @test isempty(collect(Nx.sites(e)))
        @test Nx.ProductStyle(W) isa Nx.ProductStyle

        @testset "inspection, adjoint and moment key" begin
            for a in samples
                @test Nx.degree(a) >= 0
                s = collect(Nx.sites(a))
                @test eltype(s) == Int && issorted(s) && allunique(s)
                p, b = Nx.star(a)
                q, back = Nx.star(b)
                @test back == a && isone(conj(p) * q)
                @test Nx.degree(b) == Nx.degree(a)
                key, _, _ = Nx.moment_key(a)
                @test !isless(a, key)
                @test Nx.moment_key(key) == (key, one(Nx.Phase), false)
                @test Nx.moment_reality(key) in (:real, :imaginary, :complex)
                @test clean(product_terms(e, a, C)) == Dict(a => one(C))
                @test clean(product_terms(a, e, C)) == Dict(a => one(C))
            end
        end

        @testset "order and hash" begin
            for a in samples, b in samples
                @test (a == b) + isless(a, b) + isless(b, a) == 1
                a == b && @test hash(a) == hash(b)
                Nx.degree(a) < Nx.degree(b) && @test isless(a, b)
            end
            sorted = sort(samples)
            @test issorted(sorted)
        end

        @testset "product: associative, and star reverses it" begin
            n = 0
            for a in samples, b in samples, c in samples
                (n += 1) > triples && break
                ab = product_terms(a, b, C)
                @test clean(times(ab, c)) == clean(times(a, product_terms(b, c, C)))
                sb, sa = Nx.star(b), Nx.star(a)
                reversed = Nx.mul!(TermDict{W,C}(), sb[2], sa[2], Nx.scale(sb[1] * sa[1], one(C)))
                @test clean(star_terms(ab)) == clean(reversed)
            end
        end

        @testset "no allocation" begin
            a, b = samples[end], samples[max(1, end - 1)]
            _mul_allocations(a, b, one(C), 2)
            _scalar_allocations(a, b, 2)
            @test _mul_allocations(a, b, one(C), 1000) == 0
            @test _scalar_allocations(a, b, 1000) == 0
        end
    end
    return nothing
end
