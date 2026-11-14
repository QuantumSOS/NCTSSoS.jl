# Bridge only for testing: production actions are owned by the core layer.
function Nx.act(g::TestPauliMap,w::Nx.PauliWord)
    toy=Nx.word_from_legacy(ToyPauli,Nx.legacy_letters(w))
    p,image=Nx.act(g,toy)
    p,Nx.word_from_legacy(typeof(w),Nx.legacy_letters(image))
end
@testset "Pipeline accepts packed words through contract 1" begin
    n=4; W=Nx.word_type(Nx.Pauli(),n,2)
    convertword(w)=Nx.word_from_legacy(W,Nx.legacy_letters(w))
    rows=toy_words(n,1); objective=ring_objective(n,:heisenberg)
    p,_,_=Nx.cyclic_program(rows,objective,test_shift(n),n)
    packed,_,_=Nx.cyclic_program(convertword.(rows),[(c,convertword(w)) for (c,w) in objective],test_shift(n),n)
    @test packed.nvars==p.nvars
    @test packed.forms.ids==p.forms.ids
    @test packed.forms.coefs≈p.forms.coefs
    @test packed.objective_coefs==p.objective_coefs
end
