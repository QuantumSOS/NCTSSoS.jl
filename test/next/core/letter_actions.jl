module LetterActionChecks
using Test, Random, LinearAlgebra, NCTSSoS
import NCTSSoS: Next as N

image(g,w)=last(N.act(g,w))
product(a,b)=last(N.mul(a,b))
function exercise!(out,sink,g,w)
    for i in eachindex(out)
        out[i]=image(g,w)
        empty!(sink); N.act!(sink,g,w,2.0+0im)
    end
end
function allocations(g,w)
    out=fill(w,4); sink=N.PolynomialBuffer{typeof(w),ComplexF64}(1)
    exercise!(out,sink,g,w)
    bytes=@allocated exercise!(out,sink,g,w)
    @test all(==(image(g,w)),out)
    @test length(sink)==1 && (sink.coefficients[1],sink.words[1])==(2.0+0im,image(g,w))
    bytes
end
function mixed(rng,W)
    a,b,c,d,e=[N.site_letter(W,rand(rng,1:2),rand(rng,1:3)) for _ in 1:5]
    if rand(rng,Bool)
        long=first(N.moment_key(product(product(a,b),c)))
        return product(product(long,first(N.moment_key(d))),e)
    end
    s=first(N.moment_key(product(a,b)))
    product(product(s,s),c)
end
function matrix_word(w::N.FunctionalWord{F,A,B,S},ops,psi) where {F,A,B,S}
    id=Matrix{Float64}(I,4,4)
    segment(lo,hi)=foldl((m,i)->m*ops[3*(N._fsite(N._ftoken(w,i),S)-1)+Int(N._ftoken(w,i)>>S)],lo:hi;init=id)
    phi(m)=F===N.RealTrace ? LinearAlgebra.tr(m)/4 : dot(psi,m*psi)
    value=segment(1,N._ftail(w)); i=N._ftail(w)+1
    while i<=N.degree(w)
        j=N._fstop(w,i); value=value*phi(segment(i,j)); i=j+1
    end
    value
end

function check()
    rng=MersenneTwister(442)
    @testset "Letter-selected complete mixed rows" begin
        for F in (N.RealState,N.RealTrace),A in (N.Free,N.Projector,N.Unipotent)
            algebra=N.FunctionalAlgebra(F(),N.SiteAlgebra(A(),3))
            W=N.word_type(algebra,2,4)
            letters=[N.site_letter(W,s,l) for s in 1:2 for l in 1:3]
            full=N.basis_rows(W,algebra,1:2,2)
            for mask in 0:63
                selected=letters[[!iszero(mask&(1<<(i-1))) for i in 1:6]]
                expected=filter(w->all(l->l in selected,N.support_letters(w)),full)
                actual=N.basis_rows(W,algebra,reverse(selected),2)
                @test actual==expected
                @test allunique(actual)
            end
            selected=letters[[1,3,4]]
            full=N.basis_rows(W,algebra,1:2,3)
            @test N.basis_rows(W,algebra,selected,3)==filter(w->all(l->l in selected,N.support_letters(w)),full)
            @test N.basis_rows(W,algebra,selected,0)==[one(W)]
            @test_throws r"distinct" N.basis_rows(W,algebra,[letters[1],letters[1]],2)
            for bad in (one(W),first(N.moment_key(letters[1])),product(letters[1],letters[2]))
                @test_throws r"operator letters" N.basis_rows(W,algebra,[bad],2)
            end
            @test_throws r"operator letters" N.basis_rows(W,algebra,Any[letters[1],:not_a_letter],2)
            @test_throws r"nonnegative" N.basis_rows(W,algebra,letters,-1)
        end
        # A huge ambient alphabet must not be enumerated for a small letter set.
        algebra=N.FunctionalAlgebra(N.RealState(),N.SiteAlgebra(N.Free(),1000))
        W=N.word_type(algebra,1000,4); a=N.site_letter(W,997,999)
        @test length(N.basis_rows(W,algebra,[a],2))==7
    end

    @testset "Letter permutation laws and independent matrices" begin
        for F in (N.RealState,N.RealTrace),A in (N.Free,N.Projector,N.Unipotent)
            r,(a,b)=N.site_variables(A(),2,3;functional=F())
            vars=vcat(a,b); images=[5,6,4,1,2,3]
            g=N.LetterPermutation([vars[i]=>vars[images[i]] for i in 1:6])
            h=N.LetterPermutation(a[1]=>a[2],a[2]=>a[1])
            gh=N.compose(g,h); hg=N.compose(h,g); inverse=inv(g); id=one(g)
            @test gh!=hg
            @test N.compose(g,inverse)==id
            @test N.compose(inverse,g)==id
            @test hash(id)==hash(N.LetterPermutation(a[1]=>a[1]))
            @test N.LetterPermutation(reverse([vars[i]=>vars[images[i]] for i in 1:6]))==g
            @test inv(inverse)==g && N.isone(id)
            matrices=Matrix{Float64}[]
            for i in 1:6
                v=randn(rng,2); v/=norm(v); p=v*v'
                localop=A===N.Free ? [i/7 0.2;0.2 -i/9] : A===N.Projector ? p : Matrix{Float64}(I,2,2)-2p
                push!(matrices,i<=3 ? kron(localop,Matrix{Float64}(I,2,2)) : kron(Matrix{Float64}(I,2,2),localop))
            end
            psi=randn(rng,4); psi/=norm(psi)
            for K in (1,2,4,8)
                W=N.FunctionalWord{F,A,4,2,K}
                site=N.SitePermutation([2,1]); lifted=N.LetterPermutation(W,site)
                named=N.LetterPermutation([a[i]=>b[i] for i in 1:3]...,[b[i]=>a[i] for i in 1:3]...)
                @test lifted==named && hash(lifted)==hash(named)
                for _ in 1:30
                    w,v=mixed(rng,W),mixed(rng,W)
                    @test image(id,w)==w
                    @test image(gh,w)==image(g,image(h,w))
                    @test image(inverse,image(g,w))==w
                    @test image(g,product(w,v))==product(image(g,w),image(g,v))
                    @test image(g,last(N.star(w)))==last(N.star(image(g,w)))
                    @test image(g,first(N.moment_key(w)))==first(N.moment_key(image(g,w)))
                    @test image(lifted,w)==image(site,w)
                    @test matrix_word(image(gh,w),matrices,psi) ≈ matrix_word(w,matrices[[6,5,4,1,2,3]],psi) atol=1e-10
                end
                w=mixed(rng,W)
                @test allocations(g,w)==0
                @test allocations(site,w)==0
                # Cross token/limb boundaries with repeated scalar factors.
                symbol=first(N.moment_key(N.site_letter(W,1,1)))
                for _ in 1:N._fcap(W)-N.degree(w)-1; w=product(w,symbol); end
                @test image(inverse,image(g,w))==w
                @test allocations(g,w)==0
            end
        end
        for A in (N.Free,N.Projector,N.Unipotent)
            r,(a,b)=N.site_variables(A(),2,3)
            g=N.LetterPermutation(a[1]=>a[3],a[3]=>a[1])
            exchange=N.LetterPermutation([a[i]=>b[4-i] for i in 1:3]...,[b[i]=>a[4-i] for i in 1:3]...)
            combined=N.compose(g,exchange)
            for d in (2,8,24)
                W=N.word_type(r,d)
                x=N.site_letter(W,1,1); y=N.site_letter(W,2,3)
                w=product(x,y)
                @test image(g,w)==product(N.site_letter(W,1,3),y)
                @test image(inv(g),image(g,w))==w
                @test allocations(g,w)==0
                @test allocations(exchange,w)==0
                for _ in 1:10
                    randomword()=foldl((w,_)->product(w,N.site_letter(W,rand(rng,1:2),rand(rng,1:3))),1:d÷2;init=one(W))
                    u,v=randomword(),randomword()
                    @test image(combined,u)==image(g,image(exchange,u))
                    @test image(inv(combined),image(combined,u))==u
                    @test image(combined,product(u,v))==product(image(combined,u),image(combined,v))
                    @test image(combined,last(N.star(u)))==last(N.star(image(combined,u)))
                end
            end
        end
    end
end
end
LetterActionChecks.check()

module SignedLetterChecks
using Test, Random, LinearAlgebra, NCTSSoS
import NCTSSoS: Next as N
import ..LetterActionChecks: mixed, matrix_word, product

function action(g,pw)
    p,w=pw; q,v=N.act(g,w); (p*q,v)
end
function multiply(pu,qv)
    p,u=pu; q,v=qv; r,w=N.mul(u,v); (p*q*r,w)
end
function adjoint_word(pw)
    p,w=pw; q,v=N.star(w); (conj(p)*q,v)
end
function exercise!(out,sink,g,w)
    for i in eachindex(out)
        out[i]=N.act(g,w)
        empty!(sink); N.act!(sink,g,w,2.0+0im)
    end
end
function check()
    rng=MersenneTwister(443)
    @testset "Signed letter laws and matrix substitution" begin
        for F in (nothing,N.RealState(),N.RealTrace()),A in (N.Free,N.Unipotent)
            r,(a,b)=N.site_variables(A(),2,3;functional=F)
            variables=vcat(a,b); images=[5,6,4,1,2,3]; signs=[-1,1,-1,1,-1,1]
            pairs=[variables[i]=>signs[i]*variables[images[i]] for i in 1:6]
            g=N.LetterPermutation(pairs); h=N.LetterPermutation(a[1]=>-a[2],a[2]=>a[1])
            gh=N.compose(g,h); inverse=inv(g)
            @test gh!=N.compose(h,g)
            @test N.compose(g,inverse)==one(g)==N.compose(inverse,g)
            @test inv(inverse)==g
            @test hash(g)==hash(N.LetterPermutation(reverse(pairs)))
            matrices=Matrix{Float64}[]
            for i in 1:6
                v=randn(rng,2); v/=norm(v)
                op=A===N.Free ? [i/7 .2;.2 -i/9] : Matrix{Float64}(I,2,2)-2v*v'
                push!(matrices,i<=3 ? kron(op,Matrix{Float64}(I,2,2)) : kron(Matrix{Float64}(I,2,2),op))
            end
            psi=randn(rng,4); psi/=norm(psi)
            mapped=[signs[i]*matrices[images[i]] for i in 1:6]
            ordinary(w,ops)=foldl((m,l)->m*ops[3*(N._letter_location(l)[1]-1)+N._letter_location(l)[2]],N.word_letters(w);init=Matrix{Float64}(I,4,4))
            value(w,ops)=isnothing(F) ? ordinary(w,ops) : matrix_word(w,ops,psi)
            widths=isnothing(F) ? (N.word_type(r,24),) : Tuple(N.FunctionalWord{typeof(F),A,4,2,K} for K in (1,2,4,8))
            for W in widths
                sample()=isnothing(F) ? foldl((w,_)->product(w,N.site_letter(W,rand(rng,1:2),rand(rng,1:3))),1:5;init=one(W)) : mixed(rng,W)
                for _ in 1:30
                    w,v=sample(),sample(); pw=(one(N.Phase),w); pv=(one(N.Phase),v)
                    @test action(one(g),pw)==pw
                    @test action(inverse,action(g,pw))==pw
                    @test action(gh,pw)==action(g,action(h,pw))
                    @test action(g,multiply(pw,pv))==multiply(action(g,pw),action(g,pv))
                    @test action(g,adjoint_word(pw))==adjoint_word(action(g,pw))
                    phase,image=N.act(g,w)
                    @test N.scale(phase,1.)*value(image,matrices) ≈ value(w,mapped) atol=1e-10
                    if !isnothing(F)
                        @test N.act(g,first(N.moment_key(w)))==(phase,first(N.moment_key(image)))
                    end
                end
                w=sample(); out=fill(N.act(g,w),4); sink=N.PolynomialBuffer{W,ComplexF64}(1)
                exercise!(out,sink,g,w)
                @test (@allocated exercise!(out,sink,g,w))==0
                phase,image=N.act(g,w)
                @test all(==((phase,image)),out)
                @test sink.coefficients[1]==N.scale(phase,2.0+0im) && sink.words[1]==image
            end
        end
    end
end
end
SignedLetterChecks.check()
