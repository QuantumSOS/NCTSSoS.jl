function test_cnot()
    x1,y1,z1 = (toy_letter(1,a) for a in 0:2)
    x2,y2,z2 = (toy_letter(2,a) for a in 0:2)
    ws = [Nx.mul(x1,x2)[2], Nx.mul(y1,x2)[2], z1, x2, Nx.mul(z1,y2)[2], Nx.mul(z1,z2)[2]]
    TestPauliMap([(one(Nx.Phase),w) for w in ws],false)
end
function toy_matrix(w,n)
    matrices = (ComplexF64[0 1;1 0], ComplexF64[0 -im;im 0], ComplexF64[1 0;0 -1])
    letters = Dict((l-1)÷3+1 => (l-1)%3+1 for l in Nx.legacy_letters(w))
    foldl(kron, [haskey(letters,i) ? matrices[letters[i]] : Matrix{ComplexF64}(I,2,2) for i in 1:n])
end
@testset "Degree-changing Clifford closure" begin
    g = test_cnot()
    U = ComplexF64[1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0]
    for w in toy_words(2,2)
        p,v = Nx.act(g,w)
        @test Nx.scale(p,1.0+0im) * toy_matrix(v,2) == U*toy_matrix(w,2)*U'
    end
    x1,z1,x2,z2 = toy_letter(1,0),toy_letter(1,2),toy_letter(2,0),toy_letter(2,2)
    rows = [one(ToyPauli),x1,z1,x2,z2]
    objective = [(1.,x1),(1.,Nx.mul(x1,x2)[2]),(1.,z2),(1.,Nx.mul(z1,z2)[2]),(.5,x2),(.3,z1)]
    @test_throws ArgumentError Nx.finite_program(rows,objective,[g])
    p,table,_ = Nx.cyclic_program(rows,objective,one(g),1;moment_generators=(g,))
    @test only(p.blocks).side == length(rows)
    @test Nx.orbit_moment!(table,x1) == Nx.orbit_moment!(table,Nx.mul(x1,x2)[2])
    # Independent reference: retain the original block and add the CNOT
    # relations explicitly as linear equalities, instead of quotienting ids.
    dense,dt = dense_program(rows,objective)
    plain = solve_reference(dense)[1]
    for w in toy_words(2,2)
        phase,v = Nx.act(g,w)
        iw = moment_ids!(dt,w)[1]; iv = moment_ids!(dt,v)[1]
        Nx.push_term!(dense.forms,1.,iw)
        Nx.push_term!(dense.forms,-Nx.scale(phase,1.),iv)
        Nx.close_form!(dense.forms); Nx.push_block!(dense,Nx.ZERO,1)
    end
    dense.nvars=length(dt.labels)
    reduced = solve_reference(p)[1]
    @test reduced ≈ solve_reference(dense)[1] atol=1e-6
    @test reduced > plain + 0.01
end
