module FermionicLegacyChecks
using Test, LinearAlgebra, TOML
using NCTSSoS.Next

function fixture(file,id)
    data=TOML.parsefile(joinpath(@__DIR__,"../../data/expectations",file*".toml"))
    only(filter(c->c["id"]==id,data["cases"]))["expected"]
end

function matrices(n)
    [reduce(kron,[j<i ? [1 0;0 -1] : j==i ? [0 1;0 0] : [1 0;0 1] for j in 1:n]) for i in 1:n]
end

function hubbard(a,ad,n,periodic)
    bonds=periodic ? [(i,mod1(i+1,n)) for i in 1:n] : [(i,i+1) for i in 1:n-1]
    -sum(ad[i+s]*a[j+s]+ad[j+s]*a[i+s] for (i,j) in bonds for s in (0,n))+
        4sum(ad[i]*a[i]*ad[i+n]*a[i+n] for i in 1:n)
end

function hubbard_checks(solver;references=false)
    @testset "Hubbard physical sectors" begin
        for n in (2,4)
            n==4 && !references && continue
            r,(a,ad)=fermionic_variables(2n)
            h=hubbard(a,ad,n,true)
            m=matrices(2n); matrix=hubbard(m,adjoint.(m),n,true)
            exact=eigmin(Hermitian(matrix))
            mask=(1<<n)-1
            keep=[s+1 for s in 0:(1<<(2n))-1 if count_ones(s&mask)==n÷2 && count_ones(s>>n)==n÷2]
            canonical=eigmin(Hermitian(matrix[keep,keep]))
            @test exact ≈ (n==2 ? -2.4721359549995796 : -3.4185507188738447) atol=1e-12
            @test canonical ≈ (n==2 ? exact : -2.102748483462074) atol=1e-12
            n==4 && @test exact<canonical
            plain=solve(Problem(h),Relaxation(order=2);solver)
            @test plain.bound ≈ fixture("hubbard","periodic_n$(n)_u4_order2")["objective"] atol=(n==2 ? 1e-6 : 5e-6)
            @test plain.bound<=exact+1e-6
            up=sum(ad[i]*a[i] for i in 1:n)
            down=sum(ad[i]*a[i] for i in n+1:2n)
            relaxation=Relaxation(order=2,strengthen=(Sector(up,n÷2),Sector(down,n÷2)))
            fixed=solve(Problem(h),relaxation;solver)
            if n==4
                @test fixed.bound ≈ fixture("hubbard","periodic_n4_u4_order2_canonical")["objective"] atol=1e-6
                @test fixed.bound<=canonical+1e-6
                @test fixed.bound > -3.5
            else
                symmetric=solve(Problem(h;symmetry=(SitePermutation([2,1,4,3]),NumberRotation(r))),relaxation;solver)
                @test symmetric.bound ≈ fixed.bound atol=1e-5
            end
        end
        r,(a,ad)=fermionic_variables(4)
        h=hubbard(a,ad,2,false)
        m=matrices(4); matrix=hubbard(m,adjoint.(m),2,false)
        keep=[s+1 for s in 0:15 if count_ones(s)==2]
        exact=eigmin(Hermitian(matrix[keep,keep]))
        @test exact ≈ 2-sqrt(8) atol=1e-12
        relaxation=Relaxation(order=2,strengthen=Sector(sum(ad[i]*a[i] for i in 1:4),2))
        plain=solve(Problem(h),relaxation;solver)
        reduced=solve(Problem(h;symmetry=(NumberRotation(r),)),relaxation;solver)
        @test plain.bound ≈ fixture("hubbard","open_n2_u4_dimer_order2")["objective"] atol=5e-6
        @test plain.bound ≈ exact atol=5e-6
        @test reduced.bound ≈ plain.bound atol=5e-6
        @test reduced.bound ≈ exact atol=5e-6
    end
end

function molecule(a,ad,h,v)
    result=sum(h[p,q]*ad[p+s]*a[q+s] for p in 1:2,q in 1:2,s in (0,2))
    for p in 1:2,q in 1:2,r in 1:2,t in 1:2,s in (0,2),u in (0,2)
        result+=0.5v[p,q,r,t]*ad[p+s]*ad[q+u]*a[t+u]*a[r+s]
    end
    result
end

function molecular_checks(solver)
    @testset "Plain molecular number sectors" begin
        _,(a,ad)=fermionic_variables(4)
        h=[-1.252477303982 0;0 -0.475934275367]
        v=zeros(2,2,2,2)
        v[1,1,1,1]=0.674493166046; v[2,2,2,2]=0.69739794957
        v[1,2,1,2]=v[2,1,2,1]=0.663472101055
        v[1,2,2,1]=v[2,1,1,2]=v[1,1,2,2]=v[2,2,1,1]=0.181287518306
        ham=molecule(a,ad,h,v)
        m=matrices(4); matrix=molecule(m,adjoint.(m),h,v)
        keep=[s+1 for s in 0:15 if count_ones(s)==2]
        exact=eigmin(Hermitian(matrix[keep,keep]))
        @test exact ≈ fixture("fermionic_symmetry","h2_sto3g_n2_spin_adapted")["objective"] atol=1e-10
        number=sum(ad[i]*a[i] for i in 1:4)
        result=solve(Problem(ham),Relaxation(order=2,strengthen=Sector(number,2));solver)
        @test result.bound ≈ exact atol=5e-6

        h=[-7.9386904776328775 0.27373827553230023;0.27373827553230023 -1.7707398712114135]
        ham=sum(h[p,q]*ad[p+s]*a[q+s] for p in 1:2,q in 1:2,s in (0,2))
        rows=[one(ham);[ad[p+s]*a[q+s] for s in (0,2) for p in 1:2 for q in 1:2]]
        result=solve(Problem(ham),Relaxation(;rows,strengthen=Sector(number,2));solver)
        @test result.bound ≈ 2eigmin(Symmetric(h)) atol=5e-6
        @test result.bound ≈ fixture("fermionic_symmetry","be_1s_2s_one_body_n2_order1")["objective"] atol=5e-6
    end
end

function check(;solver,sos_solver=nothing,references=false)
    for backend in (isnothing(sos_solver) ? (solver,) : (solver,sos_solver))
        hubbard_checks(backend;references)
        molecular_checks(backend)
    end
end
end
