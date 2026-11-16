function arithmetic_contracts()
    @testset "Public polynomial values and input preservation" begin
        r,(vars,)=site_variables(Free(),1,4)
        x,y,z,w=vars
        # Compare coefficients by named words, never by storage order.
        terms(p)=Dict(Tuple(NCTSSoS.operator_name(r,l) for l in NCTSSoS.word_letters(word))=>c
            for (c,word) in NCTSSoS.terms(p))
        p=NCTSSoS.polynomial(r,(x*y,2z,-3x))
        @test terms(p)==Dict((:s1_1,:s1_2)=>1,(:s1_3,)=>2,(:s1_1,)=>-3)
        @test terms(5one(x))==Dict(()=>5)
        @test terms(x*y*z)==Dict((:s1_1,:s1_2,:s1_3)=>1)
        @test terms(3x*y)==Dict((:s1_1,:s1_2)=>3)
        @test NCTSSoS.degree(x+2y*z*w)==3
        @test NCTSSoS.degree(zero(p))==0 # Approved replacement for the old -Inf convention.
        @test isone(one(p)) && NCTSSoS.degree(one(p))==0
        @test isempty(NCTSSoS.terms(zero(p)))
        @test terms(-p)==Dict((:s1_1,:s1_2)=>-1,(:s1_3,)=>-2,(:s1_1,)=>3)
        @test terms((x+2y)+(3y+4z))==Dict((:s1_1,)=>1,(:s1_2,)=>5,(:s1_3,)=>4)
        @test terms(3x-(x+2y))==Dict((:s1_1,)=>2,(:s1_2,)=>-2)
        @test terms(sum([x,x*y,z]))==Dict((:s1_1,)=>1,(:s1_1,:s1_2)=>1,(:s1_3,)=>1)
        @test sort([NCTSSoS.degree(word) for (_,word) in NCTSSoS.terms(x+x^2)])==[1,2]
        @test terms(3*(2x*y))==terms((2x*y)*3)==Dict((:s1_1,:s1_2)=>6)
        @test terms(Float32(2)*(x+2y))==Dict((:s1_1,)=>2,(:s1_2,)=>4)
        @test terms((6x)/2)==Dict((:s1_1,)=>3)
        @test_throws ArgumentError p/0
        @test terms((2x)*(3y))==Dict((:s1_1,:s1_2)=>6)
        @test terms((x+2y)*(3x))==Dict((:s1_1,:s1_1)=>3,(:s1_2,:s1_1)=>6)
        @test terms(2x+3)==terms(3+2x)==Dict((:s1_1,)=>2,()=>3)
        @test terms(2x-3)==Dict((:s1_1,)=>2,()=>-3)
        @test terms(5-2x)==Dict((:s1_1,)=>-2,()=>5)
        @test terms(3x-y)==Dict((:s1_1,)=>3,(:s1_2,)=>-1)
        @test terms(y-3x)==Dict((:s1_1,)=>-3,(:s1_2,)=>1)
        @test p*one(p)==one(p)*p==p
        @test iszero(p*zero(p)) && iszero(zero(p)*p) && iszero(0p)
        @test (x*y+2z+3)'==y*x+2z+3
        a,b=(1+2im)*x*y+3z,2y*z+5
        saved=deepcopy.(collect.(NCTSSoS.terms.((a,b))))
        a+b; a*b; a'
        @test collect.(NCTSSoS.terms.((a,b)))==saved
        @test (a*b)'==b'*a'
        support(p)=Set(NCTSSoS.operator_name(r,l) for (_,word) in NCTSSoS.terms(p)
            for l in NCTSSoS.support_letters(word))
        @test support(x*y+y*z)==Set((:s1_1,:s1_2,:s1_3))
        @test isempty(support(zero(p)))
        @test support(x^3)==Set((:s1_1,))
        _,(left,right)=site_variables(Free(),2,3)
        @test left[1]*(right[2]+2right[3])==(right[2]+2right[3])*left[1]
    end
    @testset "Public algebra factories and adjoints" begin
        for factory in (pauli_variables,fermionic_variables,bosonic_variables)
            _,ops=factory(3)
            @test all(length(a)==3 for a in ops)
        end
        for rule in (Free(),Projector(),Unipotent())
            _,ops=site_variables(rule,2,3)
            @test length(ops)==2 && all(length(a)==3 for a in ops)
        end
        _,(x,y,z)=pauli_variables(2)
        @test (2x[1])*(3x[1])==6one(x[1])
        @test NCTSSoS.degree(x[1]^2)==0
        @test ((1+2im)*x[1])'==(1-2im)*x[1]
        @test 3x[1]*x[2]==x[1]*x[2]*3
        for factory in (fermionic_variables,bosonic_variables)
            r,(a,ad)=factory(3)
            @test NCTSSoS.degree(ad[1]*a[1]*ad[2]*a[2])==4
            @test (ad[1]*a[2]+2ad[2]*a[1])'==ad[2]*a[1]+2ad[1]*a[2]
            names=Set(NCTSSoS.operator_name(r,l) for (_,word) in NCTSSoS.terms(ad[3]*ad[1]*a[2])
                for l in NCTSSoS.support_letters(word))
            @test names==Set(factory===fermionic_variables ? (:ad3,:ad1,:a2) : (:bd3,:bd1,:b2))
        end
    end
end
