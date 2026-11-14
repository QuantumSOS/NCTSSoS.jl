struct TestFermiMap <: Nx.AbstractAction
    permutation::Vector{Int}
    phase::Int
    anti::Bool
end
Nx.is_antiunitary(g::TestFermiMap)=g.anti
Base.:(==)(g::TestFermiMap,h::TestFermiMap)=(g.permutation,g.phase,g.anti)==(h.permutation,h.phase,h.anti)
Base.hash(g::TestFermiMap,h::UInt)=hash((g.permutation,g.phase,g.anti),h)
Base.one(g::TestFermiMap)=TestFermiMap(collect(eachindex(g.permutation)),0,false)
Base.isone(g::TestFermiMap)=g==one(g)
Nx.compose(g::TestFermiMap,h::TestFermiMap)=TestFermiMap(g.permutation[h.permutation],mod(g.phase+(g.anti ? -h.phase : h.phase),4),xor(g.anti,h.anti))
function Base.inv(g::TestFermiMap)
    p=invperm(g.permutation)
    TestFermiMap(p,mod(g.anti ? g.phase : -g.phase,4),g.anti)
end
function Nx.act(g::TestFermiMap,w::Nx.FermiWord)
    letters=Nx.legacy_letters(w)
    image=[sign(l)*g.permutation[abs(l)] for l in letters]
    inversions=sum(image[i]>image[j] for i in eachindex(image) for j in (i+1):length(image);init=0)
    phase=Nx.Phase(2inversions+g.phase*sum(sign,letters;init=0))
    phase,Nx.word_from_legacy(typeof(w),sort!(image))
end
function coordinate_value(m,y)
    (m.re==0 ? 0im : m.real_coefficient*y[m.re])+(m.im==0 ? 0im : m.imag_coefficient*y[m.im])
end

@testset "Complex orbit moments and CAR sums" begin
    W=Nx.word_type(Nx.Fermionic(),3,2)
    a=Nx.fermi_letter(W,1); adag=Nx.fermi_letter(W,-1)
    for phase in 0:3
        anti=TestFermiMap([1,2,3],phase,true)
        table=Nx.OrbitMoments(W,(anti,))
        m=Nx.orbit_coordinates!(table,a); md=Nx.orbit_coordinates!(table,adag)
        y=collect(1.:length(table.labels))
        z=coordinate_value(m,y)
        @test coordinate_value(md,y)==conj(z)
        @test Nx.scale(Nx.Phase(phase),z)==conj(z)
        @test !iszero(z) # an antiunitary sign does not force a complex moment to zero
    end
    table=Nx.OrbitMoments(W,(TestFermiMap([1,2,3],1,false),))
    @test coordinate_value(Nx.orbit_coordinates!(table,a),[1.])==0
    rows=[one(W);[Nx.fermi_letter(W,i) for i in [-3,-2,-1,1,2,3]]]
    for flux in (0.,0.37)
        c=-cis(flux)
        objective=Tuple{ComplexF64,W}[]
        for i in 1:3
            j=mod1(i+1,3)
            push!(objective,(c,Nx.word_from_legacy(W,[-i,j])))
            push!(objective,(conj(c),Nx.word_from_legacy(W,[-j,i])))
        end
        dense,_=dense_program(rows,objective)
        p,table,_=Nx.cyclic_program(rows,objective,Nx.SitePermutation([2,3,1]),3)
        @test table.representatives isa Nx.MomentTable{W}
        plain,_=Nx.plain_program(rows,objective)
        expected=sum(min(0.,-2cos(flux+2pi*k/3)) for k in 0:2)
        @test solve_reference(p)[1]≈solve_reference(dense)[1] atol=1e-6
        @test solve_reference(p)[1]≈solve_reference(plain)[1] atol=1e-6
        @test solve_reference(p)[1]≈expected atol=1e-6
        if iszero(flux)
            group=[TestFermiMap([2,3,1],0,false),TestFermiMap([3,2,1],0,false)]
            generic,_,_=Nx.finite_program(rows,objective,group)
            anti,_,_=Nx.finite_program(rows,objective,group;moment_generators=(TestFermiMap([1,2,3],1,true),))
            @test all(b.kind==Nx.PSD for b in anti.blocks)
            @test [b.side for b in anti.blocks]==[b.side for b in generic.blocks]
            @test solve_reference(anti)[1]≈expected atol=1e-6
        end
    end
    g=Nx.compose(test_cnot(),test_conjugate(2))
    @test_throws r"No same-side real basis" Nx.finite_program(toy_words(2,1),[],[one(g)];moment_generators=(g,))
end
