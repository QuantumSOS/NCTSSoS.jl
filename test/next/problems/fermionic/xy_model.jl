function xy_model(solver)
    @testset "XY boundary sectors" begin
        n=4
        eye=Matrix{ComplexF64}(I,2,2)
        xmat,ymat=ComplexF64[0 1;1 0],ComplexF64[0 -im;im 0]
        localop(m,i)=reduce(kron,[j==i ? m : eye for j in 1:n])
        matrix=sum((localop(xmat,i)*localop(xmat,mod1(i+1,n))+
            localop(ymat,i)*localop(ymat,mod1(i+1,n)))/4 for i in 1:n)
        @test eigmin(Hermitian(matrix)) ≈ -sqrt(2) atol=1e-12
        _,(x,y,_)=pauli_variables(n)
        h=sum((x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)])/4 for i in 1:n)
        built=NCTSSoS.build(Problem(h),Relaxation(order=2))
        @test solve(built;solver).bound ≈ fixture("xy_model","pbc_pauli_n4_order2")["objective"] atol=1e-6
        structure=fixture("next_manybody","xy_pauli_n4_order2")
        @test NCTSSoS.block_sides(built) == structure["sides"]
        @test built.sizes.moments == structure["moments"]
        _,(a,ad)=fermionic_variables(n)
        for (sign,id,exact) in ((-1,"apbc_fermionic_n4_order1",-sqrt(2)),(1,"pbc_fermionic_n4_order1",-1))
            h=(sum(ad[i]*a[i+1]+ad[i+1]*a[i] for i in 1:n-1)+sign*(ad[n]*a[1]+ad[1]*a[n]))/2
            built=NCTSSoS.build(Problem(h),Relaxation(order=1))
            result=solve(built;solver)
            @test result.bound ≈ exact atol=1e-6
            @test result.bound ≈ fixture("xy_model",id)["objective"] atol=1e-6
            structure=fixture("next_manybody","fermionic_n4_order1")
            @test NCTSSoS.block_sides(built) == structure["sides"]
            @test built.sizes.moments == structure["moments"]
        end
    end
end
