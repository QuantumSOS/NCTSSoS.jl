# Bridge only for testing: production actions are owned by the core layer.
function Nx.act(g::TestPauliMap,w::Nx.PauliWord)
    toy=Nx.word_from_legacy(ToyPauli,Nx.legacy_letters(w))
    p,image=Nx.act(g,toy)
    p,Nx.word_from_legacy(typeof(w),Nx.legacy_letters(image))
end
function core_action(g::TestPauliMap,::Type{W}) where {W}
    images=[Nx.act(g,toy_letter(i,a)) for a in (0,2) for i in 1:(length(g.images)÷3)]
    Nx.CliffordMap([Nx.word_from_legacy(W,Nx.legacy_letters(w)) for (_,w) in images],
        [p.k==2 for (p,_) in images];antiunitary=Nx.is_antiunitary(g))
end
@testset "Pipeline accepts packed words through contract 1" begin
    n=4; W=Nx.word_type(Nx.Pauli(),n,2)
    convertword(w)=Nx.word_from_legacy(W,Nx.legacy_letters(w))
    rows=toy_words(n,1); objective=ring_objective(n,:heisenberg)
    p,_,_=Nx.cyclic_program(rows,objective,test_shift(n),n)
    packed,_,_=Nx.cyclic_program(convertword.(rows),[(c,convertword(w)) for (c,w) in objective],Nx.SitePermutation([2,3,4,1]),n)
    @test packed.nvars==p.nvars
    @test packed.forms.ids==p.forms.ids
    @test packed.forms.coefs≈p.forms.coefs
    @test packed.objective_coefs==p.objective_coefs
    packed_rows=Nx.basis_rows(W,Nx.Pauli(),1:n,1)
    packed_objective=[(c,convertword(w)) for (c,w) in objective]
    generic,_,_=Nx.finite_program(packed_rows,packed_objective,
        [Nx.SitePermutation([2,3,4,1]),Nx.SitePermutation([4,3,2,1])])
    plain,_=Nx.plain_program(packed_rows,packed_objective)
    @test solve_reference(generic)[1]≈solve_reference(plain)[1] atol=1e-6
    real_generic,_,_=Nx.finite_program(packed_rows,packed_objective,
        [Nx.SitePermutation([2,3,4,1]),Nx.SitePermutation([4,3,2,1])];
        moment_generators=(Nx.AntiunitaryPauli(one(W)),))
    @test [b.side for b in real_generic.blocks]==[b.side for b in generic.blocks]
    @test all(b.kind==Nx.PSD for b in real_generic.blocks)
    @test solve_reference(real_generic)[1]≈solve_reference(plain)[1] atol=1e-6
    atoms=[[Nx.pauli_letter(W,i,a) for a in 0:2] for i in 1:n]
    discovered,_,_=Nx.generators_program(packed_rows,packed_objective,
        [Nx.SitePermutation([2,3,4,1]),Nx.SitePermutation([4,3,2,1])],atoms)
    @test solve_reference(discovered)[1]≈solve_reference(plain)[1] atol=1e-6
    mixed=Nx.AbstractAction[Nx.SitePermutation([2,3,4,1]),
        core_action(test_map(n;axes=(1,2,0)),W),
        core_action(test_map(n;axes=(1,0,2),signs=(1,-1,1)),W),Nx.SitePermutation([4,3,2,1])]
    mixed_program,_,_=Nx.generators_program(packed_rows,packed_objective,
        [mixed;Nx.AntiunitaryPauli(one(W))],atoms)
    @test all(b.kind==Nx.PSD for b in mixed_program.blocks)
    @test solve_reference(mixed_program)[1]≈solve_reference(plain)[1] atol=1e-6
    finite_mixed,_,_=Nx.finite_program(packed_rows,packed_objective,mixed;
        moment_generators=(Nx.AntiunitaryPauli(one(W)),))
    @test solve_reference(finite_mixed)[1]≈solve_reference(plain)[1] atol=1e-6
end
