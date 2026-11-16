function xx_field(solver)
    @testset "XX chain with field: different program" begin
        n = 4
        xmat,zmat = [0. 1;1 0],[1. 0;0 -1]
        eye = Matrix{Float64}(I,2,2)
        localop(m,i) = reduce(kron,[j==i ? m : eye for j in 1:n])
        matrix = -sum(localop(xmat,i)*localop(xmat,i+1) for i in 1:n-1)-sum(localop(zmat,i) for i in 1:n)
        exact = eigmin(Hermitian(matrix))
        _,(x,_,z) = pauli_variables(n)
        h = -sum(x[i]*x[i+1] for i in 1:n-1)-sum(z)
        built = Next.build(Problem(h),Relaxation(order=1))
        result = solve(built;solver)
        oracle = fixture("next_manybody","xx_field_n4_order1")
        # The old -7 is the trivial coefficient bound, not the new program's
        # optimum. Cause in the old program not established.
        @test result.bound <= exact+1e-6
        @test result.bound >= -7-1e-6
        @test result.bound ≈ oracle["objective"] atol=1e-6
        @test Next.block_sides(built) == oracle["sides"]
        @test result.sizes.moments == oracle["moments"]
    end
end

# Explicit Jordan–Wigner matrices, independent of package multiplication.
function chain_matrices(n)
    z = ComplexF64[1 0;0 -1]
    lower = ComplexF64[0 0;1 0]
    eye = Matrix{ComplexF64}(I,2,2)
    return [reduce(kron,[j<i ? z : j==i ? lower : eye for j in 1:n]) for i in 1:n]
end

function kitaev_chain(solver)
    @testset "Kitaev chain reviewed phases" begin
        matrices = chain_matrices(4)
        energy(a,ad,mu,delta) = -mu*sum(ad.*a)-sum(ad[i]*a[i+1]+ad[i+1]*a[i] for i in 1:3)+
            delta*sum(a[i]*a[i+1]+ad[i+1]*ad[i] for i in 1:3)
        _,(a,ad) = fermionic_variables(4)
        for (mu,delta,id) in ((0,1,"sweet_spot_n4_order1"),(3,1,"trivial_phase_mu3_n4_order1"),
                (0,0.5,"asymmetric_pairing_delta0p5_n4_order1"),(0,0,"pure_hopping_n4_order1"))
            oracle = fixture("kitaev_chain",id)
            exact = eigmin(Hermitian(energy(matrices,adjoint.(matrices),mu,delta)))
            built = Next.build(Problem(energy(a,ad,mu,delta)),Relaxation(order=1))
            result = solve(built;solver)
            @test result.bound ≈ exact atol=1e-6
            @test result.bound ≈ oracle["objective"] atol=1e-6
            structure = fixture("next_manybody","fermionic_n4_order1")
            @test Next.block_sides(built) == structure["sides"]
            @test built.sizes.moments == structure["moments"]
        end
        # Independently retain the Majorana-to-Pauli sweet-spot identity.
        gamma = [isodd(k) ? matrices[(k+1)÷2]+matrices[(k+1)÷2]' :
            im*(matrices[k÷2]'-matrices[k÷2]) for k in 1:8]
        majorana = sum(im*gamma[2i]*gamma[2i+1] for i in 1:3)
        eye = Matrix{ComplexF64}(I,2,2); xmat=ComplexF64[0 1;1 0]
        xs = [reduce(kron,[j==i ? xmat : eye for j in 1:4]) for i in 1:4]
        @test isapprox(majorana,sum(xs[i]*xs[i+1] for i in 1:3);atol=1e-12,rtol=0)
        @test eigmin(Hermitian(majorana)) ≈ -3 atol=1e-12
        _,(x,_,_) = pauli_variables(4)
        built = Next.build(Problem(sum(x[i]*x[i+1] for i in 1:3)),Relaxation(order=1))
        result = solve(built;solver)
        @test result.bound ≈ fixture("fermionic_chain","kitaev_majorana_n4_order1")["objective"] atol=1e-6
        structure = fixture("next_manybody","majorana_pauli_n4_order1")
        @test Next.block_sides(built) == structure["sides"]
        @test built.sizes.moments == structure["moments"]
    end
end
