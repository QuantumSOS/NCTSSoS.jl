@testset "Verified translation discovery and label guards" begin
    for (nx,ny) in ((2,3),(3,3),(6,6))
        tx,ty,point=torus_actions(nx,ny)
        atoms=[[toy_letter(i,a) for a in 0:2] for i in 1:(nx*ny)]
        structure=Nx.discover_lattice([point;tx;ty],atoms)
        @test length(structure.translations.elements)==nx*ny
        @test isempty(structure.moment_generators)
    end
    n=5; shift=test_shift(n); reflection=test_reflect(n)
    atoms=[[toy_letter(i,a) for a in 0:2] for i in 1:n]
    # Neither supplied generator is a translation: their product is one.
    structure=Nx.discover_lattice([reflection,Nx.compose(shift,reflection)],atoms)
    @test length(structure.translations.elements)==5
    @test_throws r"bounded generator/product search" Nx.discover_lattice([reflection],atoms)
    @test Nx.character_sector_closure(structure.translations,structure.point_generators,[1])==[1]
    @test_throws r"physical state charge" Nx.character_sector_closure(structure.translations,structure.point_generators,[2])
    @test length(Nx.character_sector_closure(structure.translations,structure.point_generators,[2];union=true))==2
    # A guard computes a label orbit, not particle-number equations or a
    # modified program. Actual physical state selection is strengthening.
    rows=toy_words(n,1); objective=ring_objective(n,:heisenberg)
    p,_,_=Nx.generators_program(rows,objective,[reflection,Nx.compose(shift,reflection)],atoms)
    reference,_=dense_program(rows,objective)
    @test solve_reference(p)[1]≈solve_reference(reference)[1] atol=1e-6
    # Finite spin-charge labels: a quarter turn has characters ±i, and a
    # spin flip exchanges them. They are row labels, not state charges.
    rotation=test_map(2;axes=(1,0,2),signs=(1,-1,1))
    flip=test_map(2;signs=(1,-1,-1))
    labels=Nx.TranslationGroup([rotation],[4])
    @test_throws ArgumentError Nx.character_sector_closure(labels,[flip],[2])
    @test Nx.character_sector_closure(labels,[flip],[2];union=true)==[2,4]
    rows=toy_words(2,1); objective=ring_objective(2,:heisenberg)
    labelled,_,_=Nx.abelian_lattice_program(rows,objective,labels,[flip])
    dense,_=dense_program(rows,objective)
    @test solve_reference(labelled)[1]≈solve_reference(dense)[1] atol=1e-6
    @test_throws r"physical state charge" Nx.abelian_lattice_program(rows,objective,labels,[flip];selected_characters=[2])
    selected,_,report=Nx.abelian_lattice_program(rows,objective,labels,[flip];selected_characters=[2],union_characters=true)
    @test report.selected_characters==[2,4]
    @test length(selected.blocks)<length(labelled.blocks)
    @test all(b.kind==Nx.HPSD for b in selected.blocks) # no state-sector equations
    # The maximally mixed state is still feasible, although it is not in a
    # fixed total-spin sector. Row selection must not constrain that state.
    mixed=zeros(selected.nvars); mixed[1]=1
    @test minimum(eigmin(Hermitian(A)) for A in evaluated_blocks(selected,mixed))>=-1e-9
end
