# The swap/projector equalities are operator identities, not scalar constraints.
function heisenberg_pairs(n; projectors=false)
    pairs=[(i,j) for i in 1:n for j in i+1:n]
    _,(v,)=site_variables(projectors ? Projector() : Unipotent(),1,length(pairs))
    index=Dict(p=>k for (k,p) in enumerate(pairs))
    op(i,j)=v[index[minmax(i,j)]]
    equalities=unique([projectors ?
        op(i,j)*op(j,k)+op(j,k)*op(i,j)-(op(i,j)+op(j,k)-op(i,k))/2 :
        op(i,j)*op(j,k)+op(j,k)*op(i,j)-op(i,j)-op(j,k)-op(i,k)+1
        for i in 1:n for j in 1:n for k in 1:n if i!=j && j!=k && i!=k])
    h=projectors ? -sum(2op(i,mod1(i+1,n))+op(i,mod1(i+2,n)) for i in 1:n) :
        sum(op(1,j) for j in 2:n)
    Problem(h;constraints=Zero.(equalities))
end

function heisenberg_star(solver,sos_solver;references=false)
    @testset "Heisenberg star and Majumdar–Ghosh" begin
        for n in (8,10)
            if !references
                @test_skip "Heisenberg star n=$n needs references=true (over 20 seconds)"
                continue
            end
            problem=heisenberg_pairs(n)
            result=solve(problem,Relaxation(order=1);solver)
            @test result.bound ≈ -1 atol=1e-6
            if n==10
                sparse=solve(problem,Relaxation(order=1,grouping=Correlative(algorithm=MF()));solver)
                @test sparse.bound ≈ -1 atol=1e-6
            elseif isnothing(sos_solver)
                @test_skip "Heisenberg star n=8 SOS needs a licensed backend"
            else
                @test solve(problem,Relaxation(order=1);solver=sos_solver).bound ≈ -1 atol=1e-6
            end
        end
        @test solve(heisenberg_pairs(6;projectors=true),Relaxation(order=1);solver).bound ≈ -9 atol=1e-4
    end
end

function heisenberg_spin(solver)
    @testset "Spin invariance and singlet restriction are separate" begin
        for n in (2,4)
            registry,(x,y,z)=pauli_variables(n)
            h=n==2 ? (x[1]*x[2]+y[1]*y[2]+z[1]*z[2])/4 :
                sum(op[i]*op[mod1(i+1,n)] for op in (x,y,z) for i in 1:n)/4
            s2=(sum(x)^2+sum(y)^2+sum(z)^2)/4
            spatial=n==2 ? (SitePermutation([2,1]),) :
                (SitePermutation([2,3,4,1]),SitePermutation([4,3,2,1]))
            plain=solve(Problem(h),Relaxation(order=1);solver)
            # Independent spin matrices give the exact unrestricted minimum.
            eye=Matrix{ComplexF64}(I,2,2)
            matrices=(ComplexF64[0 1;1 0],ComplexF64[0 -im;im 0],ComplexF64[1 0;0 -1])
            localop(m,i)=reduce(kron,[j==i ? m : eye for j in 1:n])
            edges=n==2 ? [(1,2)] : [(i,mod1(i+1,n)) for i in 1:n]
            exact=eigmin(Hermitian(sum(localop(m,i)*localop(m,j) for m in matrices for (i,j) in edges)/4))
            @test exact ≈ (n==2 ? -0.75 : -2) atol=1e-12
            @test plain.bound <= exact+1e-6
            n==2 && @test plain.bound ≈ -0.75 atol=1e-6
            for (label,symmetry,strengthen) in (
                (:spatial,spatial,()),
                (:rotation,(SpinRotation(registry,:z),),()),
                (:singlet,(),(Sector(s2,0),)),
                (:combined,(spatial...,SpinRotation(registry,:z)),(Sector(s2,0),)))
                built=Next.build(Problem(h;symmetry),Relaxation(order=1;strengthen))
                result=solve(built;solver)
                @test result.bound ≈ plain.bound atol=1e-6
                expected=fixture("next_heisenberg","n$(n)_$(label)")
                @test Next.block_sides(built)==expected["sides"]
                @test built.sizes.moments==expected["moments"]
            end
        end
    end
end

function heisenberg_structure()
    @testset "Heisenberg support-split sizes at 4d29c6e" begin
        n=16
        registry,(x,y,z)=pauli_variables(n)
        h=sum(op[i]*op[mod1(i+1,n)] for op in (x,y,z) for i in 1:n)/4
        spatial=(SitePermutation([2:n;1]),SitePermutation(reverse(1:n)))
        # Conjugation by the product of Z: X and Y change sign, Z does not.
        sign=CliffordMap([last(only(Next.terms(p))) for p in [x;z]],[trues(n);falses(n)])
        for (id,order,row_rule,generators,rows) in (
            ("n16_order2",2,nothing,spatial,1129),
            ("n16_contiguous4",4,Contiguous(periodic=true),(spatial...,sign),1921))
            problem=Problem(h;symmetry=(generators...,SpinRotation(registry,:z)))
            built=Next.build(problem,Relaxation(;order,row_rule))
            expected=fixture("next_heisenberg",id)
            @test built.sizes.rows==rows
            @test Next.block_sides(built)==expected["sides"]
            @test built.sizes.moments==expected["moments"]
        end
    end
end
