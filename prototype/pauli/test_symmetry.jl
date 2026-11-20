# Standalone: julia --project prototype/pauli/test_symmetry.jl
module WorkflowSymmetryChecks
import NCTSSoS as NC
using LinearAlgebra
using Test
include("symmetry.jl")

const XYZ = (ComplexF64[0 1;1 0],ComplexF64[0 -im;im 0],ComplexF64[1 0;0 -1])
word(p) = last(only(NC.terms(p)))

# Independent tensor-product reference. Obtain the letters from the public
# polynomial variables, not the multiplication or action code under test.
function matrices(vars)
    n = length(first(vars))
    operators = Dict(word(vars[a][i])=>(i,a) for i in 1:n for a in 1:3)
    function dense(w)
        axes = zeros(Int,n)
        for u in NC.word_letters(w)
            i,a = operators[u]
            axes[i] = a
        end
        reduce(kron,(a==0 ? Matrix{ComplexF64}(I,2,2) : XYZ[a] for a in axes))
    end
    dense
end

function gram(rows,dense,rho)
    ops = dense.(rows)
    Hermitian([tr(rho*A'*B) for A in ops,B in ops])
end

function coordinates(table,dense,rho)
    [l.part==0x02 ? imag(tr(rho*dense(l.word))) : real(tr(rho*dense(l.word))) for l in table.labels]
end

function evaluation(program,y)
    map(program.blocks) do block
        A = zeros(ComplexF64,block.side,block.side)
        f = block.first_form
        for j in axes(A,2),i in j:block.side
            A[i,j] = sum((program.forms.coefs[k]*y[program.forms.ids[k]] for k in NC.term_range(program.forms,f));init=0.0)
            A[j,i] = conj(A[i,j]); f+=1
        end
        Hermitian(A)
    end
end

dense_polynomial(p,dense) = sum(c*dense(w) for (c,w) in NC.terms(p))
function family_for(H,groups;order=2)
    W = eltype(H.data.words)
    [NC.RowBlock(NC.basis_rows(W,H.registry.algebra,group,order),[(1.0+0im,one(W))],NC.HPSD,0) for group in groups]
end

function assemble(H,family;spatial_mode=:none,lattice_shape=nothing,realify=true,opts...)
    symmetry = workflow_symmetry(H;family,spatial_mode,lattice_shape,opts...)
    workflow_symmetry_program(family,collect(NC.terms(H)),symmetry;spatial_mode,lattice_shape,realify)
end

@testset "Pauli workflow symmetry" begin
    @testset "Discover spatial orbits without geometry or site order" begin
        # Both sides of the term-incidence size cutoff. Swap two site labels
        # so an index shift or reversal is not a symmetry.
        for (n,edges) in ((60,[(i,mod1(i+1,60)) for i in 1:60]),
                ((L^2,[(i,j) for i in 1:L^2 for j in
                    (L*((i-1)÷L)+mod1(i+1,L),mod1(i+L,L^2))]) for L in (4,6))...)
            _,vars = NC.pauli_variables(n)
            label = collect(1:n); label[2],label[7] = label[7],label[2]
            H = sum(a[label[i]]*a[label[j]] for a in vars for (i,j) in edges)
            discover(h) = workflow_symmetry(h;parity=false,antiunitary=false,spin=false,axis_permutations=false).generators
            x = first(vars)
            function orbit(gs)
                reached = Set([word(x[1])]); pending = collect(reached)
                while !isempty(pending)
                    u = pop!(pending)
                    for g in gs
                        v = last(NC.act(g,u))
                        if v ∉ reached
                            push!(reached,v); push!(pending,v)
                        end
                    end
                end
                reached
            end
            @test orbit(discover(H)) == Set(word.(x))
            @test orbit(discover(H+x[1])) == Set([word(x[1])])
        end
    end

    @testset "Heisenberg m=0, all J and shared singlet moments" begin
        r,vars = NC.pauli_variables(4)
        x,y,z = vars
        H = sum(j*v[j]*v[j+1] for v in vars for j in 1:3)
        chi = x[1]*(y[2]*z[3]-z[2]*y[3]) + y[1]*(z[2]*x[3]-x[2]*z[3]) + z[1]*(x[2]*y[3]-y[2]*x[3])
        dense = matrices(vars)
        for objective in (H,H+chi/3)
            family = family_for(objective,[1:4])
            p,table,report = assemble(objective,family)
            @test report.route==:spin_j_m0
            @test report.labels==[[0,1,2]]
            @test [b.side÷(i in report.doubled ? 2 : 1) for (i,b) in enumerate(p.blocks)]==[7,10,6]
            objective==H && @test isempty(report.doubled)
            dotspin(i,j) = sum(v[i]*v[j] for v in vars)
            pairings = [dotspin(1,2)*dotspin(3,4),dotspin(1,3)*dotspin(2,4),dotspin(1,4)*dotspin(2,3)]
            quartic = dense_polynomial(sum(c*p for (c,p) in zip((1,2,4),pairings)),dense)
            rho = Hermitian((Matrix{ComplexF64}(I,16,16)+dense_polynomial(objective,dense)/32+quartic/128)/16)
            vals = coordinates(table,dense,rho)
            blocks = evaluation(p,vals)
            # Independent dense Gram spectrum, not a second use of the sparse
            # column constructor. Repeat each multiplicity block 2J+1 times.
            lifted = sort!(vcat((repeat(eigvals(A)[1:(J+1 in report.doubled ? 2 : 1):end],2J+1)
                for (J,A) in zip(0:2,blocks))...))
            reference = eigvals(gram(only(family).rows,dense,rho))
            @test lifted ≈ reference atol=2e-12
            @test sum(p.objective_coefs .* vals[p.objective_ids]) ≈ real(tr(rho*dense_polynomial(objective,dense))) atol=2e-12
            abc = NC.moment_terms(table,word(x[1]*y[2]*z[3]))
            @test !isnothing(abc)
            @test sum((c*vals[id] for (id,c) in abc);init=0.0) ≈ (objective==H ? 0.0 : 1/96) atol=1e-12
            xxzz = NC.moment_terms(table,word(x[1]*x[2]*z[3]*z[4]))
            @test !isnothing(xxzz)
            d1,d2,d3 = [tr(rho*dense_polynomial(p,dense)) for p in pairings]
            @test sum((c*vals[id] for (id,c) in xxzz);init=0.0) ≈ (4d1-d2-d3)/30 atol=1e-12
        end

        family = family_for(H,[1:3,2:4])
        p,_,report = assemble(H,family)
        @test report.labels==[[0,1,2],[0,1,2]]
        @test [b.side for b in p.blocks]==[4,6,3,4,6,3]

        symmetric = sum(v[j]*v[j+1] for v in vars for j in 1:3)
        p,_,report = assemble(symmetric,family;spatial_mode=:clique_orbits)
        @test report.route==:spin_j_m0
        @test report.labels==[[0,1,2]]
        @test [b.side for b in p.blocks]==[4,6,3]
    end

    @testset "Asymmetric axial fields and sector preservation" begin
        r,vars = NC.pauli_variables(2); x,y,z = vars
        H = x[1]*x[2]+y[1]*y[2]+2z[1]*z[2]+z[1]+3z[2]
        dense = matrices(vars)
        family = family_for(H,[1:2])
        p,table,report = assemble(H,family)
        @test report.route==:axial_charge
        @test report.labels==[collect(-2:2)] # K does not reverse Z charge
        rho = (Matrix{ComplexF64}(I,4,4)+dense_polynomial(H,dense)/12)/4
        vals = coordinates(table,dense,rho)
        @test sort!(vcat(eigvals.(evaluation(p,vals))...)) ≈ eigvals(gram(only(family).rows,dense,rho)) atol=2e-12

        J = sum(v[1]*v[2] for v in vars)
        K = NC.CliffordMap(NC.AntiunitaryPauli(one(typeof(word(x[1])))))
        theta = NC.CliffordMap(NC.AntiunitaryPauli(word(y[1]*y[2])))
        @test last(NC.act(K,word(z[1])))==word(z[1]) && first(NC.act(K,word(z[1]))).k==0
        @test first(NC.act(theta,word(z[1]))).k==2
        target = NC.Sector(z[1]+z[2],2)
        s = workflow_symmetry(J;target,observable=z[1])
        @test NC.is_antiunitary(last(s.moment_only)) && last(s.moment_only)==K
        @test all(g->first(NC.act(g,word(z[1]))).k==0 && last(NC.act(g,word(z[1])))==word(z[1]),
            (s.generators...,s.moment_only...))
        _,_,unpaired = assemble(J,family;target)
        _,_,paired = assemble(J,family;target=NC.Sector(z[1]+z[2],0))
        @test unpaired.labels==[collect(-2:2)] && paired.labels==[collect(0:2)]
    end

    @testset "Axial pairing on a proper subset" begin
        r,vars=NC.pauli_variables(5); x,y,z=vars
        H=x[2]*x[4]+z[2]*z[4]+2y[2]*y[4]
        G=H+y[2]-2y[4]
        family=family_for(H,[[2,4]])
        _,_,paired=assemble(H,family)
        p,t,unpaired=assemble(G,family)
        @test paired.labels==[collect(0:2)] # K reverses the Y charge
        @test unpaired.labels==[collect(-2:2)]
        dense=matrices(vars)
        rho=(Matrix{ComplexF64}(I,32,32)+dense_polynomial(G,dense)/12)/32
        blocks=evaluation(p,coordinates(t,dense,rho))
        spectra=[eigvals(A)[1:(i in unpaired.doubled ? 2 : 1):end] for (i,A) in enumerate(blocks)]
        @test sort!(vcat(spectra...)) ≈ eigvals(gram(only(family).rows,dense,rho)) atol=2e-12

        W=typeof(word(x[1]))
        flip=NC.CliffordMap(NC.PauliConjugation(word(prod(x)),false))
        axial=Set(word.(y[[2,4]]))
        swap=_symmetry_action(W,5,[1,4,3,2,5])
        outside=_symmetry_action(W,5,[1,3,2,5,4])
        other_axis=NC.compose(flip,first(_workflow_axes(W,5)))
        @test _workflow_reverses_charge(NC.compose(flip,swap),axial)
        @test !_workflow_reverses_charge(NC.compose(flip,outside),axial)
        @test !_workflow_reverses_charge(other_axis,axial)
    end

    @testset "Full commutator kernel and selected rows" begin
        r,vars=NC.pauli_variables(2); x,y,z=vars
        # Orthogonal relative frame: the three exact conserved generators
        # have different axes on sites 1 and 2.
        H = x[1]*y[2]+y[1]*z[2]+z[1]*x[2]
        dense = matrices(vars); hd = dense_polynomial(H,dense)
        letters = word.([x[1],y[1],z[1],x[2],y[2],z[2]])
        ops=dense.(letters)
        C=hcat([vec(A*hd-hd*A) for A in ops]...)
        sym=workflow_symmetry(H;spatial_mode=:none)
        @test length(sym.continuous)==size(C,2)-rank(C)==3
        for d in sym.continuous
            A=sum(c*dense(w) for (w,c) in d.columns)
            @test norm(A*hd-hd*A)<1e-12
        end
        # A projected 6×6 ad_H would incorrectly find six generators for
        # this all-two-body anisotropic Hamiltonian.
        G = x[1]*x[2]+2y[1]*y[2]+3z[1]*z[2]
        @test isempty(workflow_symmetry(G;spatial_mode=:none).continuous)
        p,_,report = assemble(G,family_for(G,[1:2];order=1))
        @test report.route==:pauli_grading
        @test sort([b.side for b in p.blocks])==[1,2,2,2]
        _,_,isotropic=assemble(sum(v[1]*v[2] for v in vars),family_for(G,[1:2];order=1);spin=false)
        @test length(only(isotropic.labels))==2
        F = [sort([one(typeof(word(x[1]))),word(x[1]),word(y[2])])]
        s = workflow_symmetry(H;family=F)
        @test all(g->Set(last(NC.act(g,w)) for w in only(F))==Set(only(F)),(s.generators...,s.moment_only...))
        @test all(d->all(w->all(in(only(F)),keys(NC._derivative(d,w))),only(F)),s.continuous)

        none = workflow_symmetry(G;spatial_mode=:none,parity=false,antiunitary=false,spin=false,axis_permutations=false)
        @test all(isempty,(none.generators,none.moment_only,none.continuous))
        p,_,rep = workflow_symmetry_program(family_for(G,[1:2];order=1),collect(NC.terms(G)),none;realify=false)
        @test only(p.blocks).kind==NC.HPSD && only(p.blocks).side==7 && !hasproperty(rep,:doubled)
    end

    @testset "Labelled pair generators, not a lattice assumption" begin
        # Petersen graph: every automorphism is induced by permuting the
        # five symbols in these 2-subsets, so the group has order 120 and no
        # element of order 10. A regular abelian group of order 10 is cyclic.
        pairs = [(i,j) for i in 1:5 for j in i+1:5]
        edges = [(i,j) for i in 1:10 for j in i+1:10 if isdisjoint(pairs[i],pairs[j])]
        a,b = first(edges)
        options = (;spin=false, parity=false, antiunitary=false, axis_permutations=false)
        for labels in (1:10,[3,8,1,10,2,7,5,9,4,6])
            _,(x,y,z) = NC.pauli_variables(10)
            H = sum(x[labels[i]]*x[labels[j]] for (i,j) in edges)
            generators = workflow_symmetry(H;options...).generators
            @test length(generators) < 10
            @test length(NC.OrbitGroup(generators)) == 120
            letters = [word.([x[i],y[i],z[i]]) for i in 1:10]
            @test isnothing(NC.discover_lattice(collect(generators),letters;required=false))
            bond = x[labels[a]]*x[labels[b]]
            @test length(NC.OrbitGroup(workflow_symmetry(H+bond;options...).generators)) == 8
            @test length(NC.OrbitGroup(workflow_symmetry(H;observable=z[labels[a]],options...).generators)) == 12
            row = word(y[labels[b]])
            sym = workflow_symmetry(H;observable=z[labels[a]],family=[[one(typeof(row)),row]],options...)
            @test length(NC.OrbitGroup(sym.generators)) == 4
            @test all(g -> last(NC.act(g,row)) == row,sym.generators)
        end

        # The incident multisets do not tell these orientations apart. A
        # reflection would preserve the unlabelled cycle but changes X_i Y_j.
        _,(x,y,_) = NC.pauli_variables(7)
        labels = [2,5,7,1,3,6,4]
        H = sum(x[labels[i]]*y[labels[mod1(i+1,7)]] for i in 1:7)
        @test length(NC.OrbitGroup(workflow_symmetry(H;options...).generators)) == 7
    end

    @testset "DFT on a ring and a rectangular torus" begin
        for n in (4,5)
            r,vars=NC.pauli_variables(n); x,y,z=vars
            if n==4 # 2×2 torus with unequal horizontal and vertical couplings.
                H = x[1]*x[2]+x[3]*x[4]+2z[1]*z[3]+2z[2]*z[4]
                copies = [1,1,1,1]
            else
                H = sum(x[j]*x[mod1(j+1,n)]+2y[j]*y[mod1(j+1,n)]+3z[j]*z[mod1(j+1,n)] for j in 1:n)
                copies = [1,2,2]
            end
            F=family_for(H,[1:n];order=1)
            p,t,rep=assemble(H,F;spatial_mode=:lattice)
            @test rep.route==:lattice
            @test [b.side for b in p.blocks]==[4;fill(3,length(copies)-1)]
            n==4 && @test rep.reduction.translation_order==4
            dense=matrices(vars)
            rho=(Matrix{ComplexF64}(I,2^n,2^n)+dense_polynomial(H,dense)/32)/2^n
            spectra=eigvals.(evaluation(p,coordinates(t,dense,rho)))
            @test sort!(vcat((repeat(λ,k) for (λ,k) in zip(spectra,copies))...)) ≈ eigvals(gram(only(F).rows,dense,rho)) atol=2e-12
            @test_throws ArgumentError assemble(H+z[1],F;spatial_mode=:lattice)
        end

        r,vars=NC.pauli_variables(4)
        H=sum(v[j]*v[mod1(j+1,4)] for v in vars for j in 1:4)
        F=family_for(H,[1:4])
        p,t,rep=assemble(H,F;spatial_mode=:lattice)
        @test rep.route==:lattice && hasproperty(rep,:continuous)
        dense=matrices(vars)
        rho=(Matrix{ComplexF64}(I,16,16)+dense_polynomial(H,dense)/32)/16
        reduced=evaluation(p,coordinates(t,dense,rho))
        reference=gram(only(F).rows,dense,rho)
        @test minimum(eigmin,reduced) ≈ eigmin(reference) atol=2e-12
        @test maximum(eigmax,reduced) ≈ eigmax(reference) atol=2e-12
    end

    @testset "Explicit translations above the graph cap" begin
        nx,ny = shape = (10,8)
        n=nx*ny
        r,(x,y,z)=NC.pauli_variables(n)
        H=sum(x[i+(j-1)*nx]*x[mod1(i+1,nx)+(j-1)*nx] +
            2z[i+(j-1)*nx]*z[i+(mod1(j+1,ny)-1)*nx] for i in 1:nx for j in 1:ny)
        @test n+3length(collect(NC.terms(H)))>512
        F=family_for(H,[1:n];order=1)
        options=(;spatial_mode=:lattice,spin=false,parity=false,axis_permutations=false,antiunitary=false)
        _,_,automatic=assemble(H,F;options...)
        @test isnothing(automatic.lattice_shape) && automatic.reduction.translation_order==n
        p,_,rep=assemble(H,F;lattice_shape=shape,realify=false,options...)
        @test rep.lattice_shape==shape && rep.reduction.translation_order==n
        @test [b.side for b in p.blocks]==[4;fill(3,n-1)]
        @test_throws ArgumentError workflow_symmetry(H;lattice_shape=(5,8),options...)
        @test_throws ArgumentError workflow_symmetry(H;lattice_shape=(8,10),options...)
        @test_throws ArgumentError workflow_symmetry(H+x[1]*x[2];lattice_shape=shape,options...)
        @test_throws ArgumentError workflow_symmetry(H+x[1]*y[2]-x[1]*x[2];lattice_shape=shape,options...)
        @test_throws ArgumentError workflow_symmetry(H;observable=z[1],lattice_shape=shape,options...)
        @test_throws ArgumentError workflow_symmetry(H;family=family_for(H,[1:(n-1)];order=1),lattice_shape=shape,options...)
    end
end
end
