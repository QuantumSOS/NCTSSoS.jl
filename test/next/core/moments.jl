function lookup_allocations(t,words)
    out = Vector{Tuple{Int32,Int32}}(undef,length(words))
    function fill_ids!()
        for i in eachindex(words); out[i]=Nx.moment_ids!(t,words[i]); end
    end
    fill_ids!()
    return @allocated fill_ids!()
end

@testset "Growing polynomial buffers and canonical guard" begin
    W = Nx.PauliWord{1}
    a,b,out,scratch = (Nx.PolynomialBuffer{W,ComplexF64}() for _ in 1:4)
    for i in 1:100
        Nx.push_term!(a,1.0,Nx.pauli_letter(W,mod1(i,20),i%3))
    end
    Nx.push_term!(b,2.0,one(W))
    @test length(a)==100
    @test_throws ArgumentError Nx.add!(out,a,b)
    Nx.canonicalize!(a)
    @test_throws ArgumentError Nx.add!(out,a,b)
    Nx.canonicalize!(b)
    Nx.add!(out,a,b)
    @test buffer_dict(out)==merge(buffer_dict(a),Dict(one(W)=>2.0+0im))
    Nx.add!(scratch,out,b)
    @test buffer_dict(scratch)[one(W)]==4.0
    Nx.mul!(out,a,b,scratch)
    @test buffer_dict(out)==Dict(w=>2c for (w,c) in buffer_dict(a))
    copied = Nx.PolynomialBuffer{W,ComplexF64}()
    copyto!(copied,a)
    Nx.add!(out,copied,b)
    Nx.push_term!(copied,1.0,one(W))
    @test_throws ArgumentError Nx.add!(out,copied,b)
    Nx.star!(out,copied)
    Nx.add!(scratch,out,b)
    @test buffer_dict(scratch)[one(W)]==3.0
end

@testset "Representative-only moment table" begin
    for (algebra,W) in ((Nx.Pauli(),Nx.PauliWord{1}),(Nx.Fermionic(),Nx.FermiWord{1}))
        words = Nx.basis_rows(W,algebra,1:5,3)
        t = Nx.MomentTable{W}(0)
        expected = Dict{W,Tuple{Int32,Int32}}()
        for w in words
            expected[w] = Nx.moment_ids!(t,w)
        end
        @test length(t)==length(words)
        @test all(Nx.moment_ids(t,w)==id for (w,id) in expected)
        Nx.reserve!(t,3length(words))
        @test all(Nx.moment_ids(t,w)==id for (w,id) in expected)
        @test lookup_allocations(t,words)==0
        @test Nx.moment_ids(t,one(W))==(1,0)
        @test Nx.nvars(t)==length(t.labels)
    end
    W = Nx.FermiWord{1}
    t = Nx.MomentTable{W}()
    a,ad = Nx.fermi_letter(W,1),Nx.fermi_letter(W,-1)
    Nx.moment_ids!(t,a)
    @test Nx.moment_ids(t,ad)==(0,0) # no implicit adjoint alias
    @test length(t)==2
    for (reality,expected) in ((:real,(2,0)),(:imaginary,(0,2)),(:complex,(2,3)))
        forced = Nx.MomentTable{W}()
        @test Nx.moment_ids!(forced,a,reality)==expected
        @test Nx.moment_ids!(forced,a,reality)==expected
        for other in (:real,:imaginary,:complex)
            other==reality && continue
            @test_throws ArgumentError Nx.moment_ids!(forced,a,other)
        end
    end
    @test_throws ArgumentError Nx.moment_ids!(t,ad,:zero)
    @test_throws ArgumentError Nx.moment_ids!(t,one(W),:imaginary)
end

@testset "Plain programs equal the dense reference" begin
    W = Nx.PauliWord{1}
    sparse_rows = Nx.basis_rows(W,Nx.Pauli(),1:3,1)
    uncovered = Nx.mul(Nx.mul(Nx.pauli_letter(W,1,0),Nx.pauli_letter(W,2,0))[2],Nx.pauli_letter(W,3,0))[2]
    error = try
        Nx.plain_program(sparse_rows,[(1.0,uncovered)])
        nothing
    catch e
        e
    end
    @test error isa ArgumentError
    @test occursin("degree 2",sprint(showerror,error))
    @test occursin(string(uncovered),sprint(showerror,error))
    covered,_ = Nx.plain_program(Nx.basis_rows(W,Nx.Pauli(),1:3,2),[(1.0,uncovered)])
    @test length(covered.objective_ids)==1
    for (algebra,W) in ((Nx.Pauli(),Nx.PauliWord{2}),(Nx.Fermionic(),Nx.FermiWord{2}),
            ((Nx.SiteAlgebra(rule,2),Nx.word_type(Nx.SiteAlgebra(rule,2),3,4)) for rule in (Nx.Free(),Nx.Projector(),Nx.Unipotent()))...)
        rows = Nx.basis_rows(W,algebra,1:3,2)
        @test issorted(rows) && allunique(rows)
        @test all(Nx.degree(w)<=2 for w in rows)
        objective = [(0.7,rows[2]),(-1.2,rows[end])]
        p,t = Nx.plain_program(rows,objective)
        ref,rt = dense_program(rows,objective)
        @test p.nvars==ref.nvars
        @test p.forms.offsets==ref.forms.offsets
        @test p.forms.ids==ref.forms.ids
        @test p.forms.coefs==ref.forms.coefs
        @test Dict(zip(p.objective_ids,p.objective_coefs))==
            Dict(id=>c for (id,c) in zip(ref.objective_ids,ref.objective_coefs) if !iszero(c))
        @test t.labels==rt.labels
    end
end
