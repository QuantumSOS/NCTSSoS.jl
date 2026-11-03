using Test, NCTSSoS, JuMP, LinearAlgebra

@testset "seniority-zero paired constraints" begin
    @test_throws ArgumentError build_paired_moment_data(ones(2,3); pairs=1)
    @test_throws ArgumentError build_paired_moment_data([1.0 2; 0 3]; pairs=1)
    @test_throws ArgumentError build_paired_moment_data(ones(2,2); pairs=3)
    @test_throws ArgumentError build_paired_moment_data(ones(2,2); pairs=1, conditions=:T2)

    function energy(J, pairs, conditions)
        linear = build_paired_moment_data(J; pairs, conditions)
        M = typeof(first(values(linear.key_to_monomial)))
        T = eltype(linear.identity)
        P = Polynomial{FermionicAlgebra,T,Float64}
        mp = NCTSSoS._moment_problem_with_linear(FermionicAlgebra,T,M,P,
            zero(P),Tuple{Symbol,Matrix{P}}[],M[],0,linear)
        model = NCTSSoS.sos_dualize(mp).model
        set_optimizer(model,SOLVER)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end

    # With one pair the exact Hamiltonian is J itself. Asymmetric diagonals and
    # mixed hopping signs expose sign, double-counting, and number mistakes.
    J = [-0.7 0.2 -0.4; 0.2 0.9 0.1; -0.4 0.1 1.3]
    for conditions in (:PQG,:PQGT)
        @test energy(J,0,conditions) ≈ 0 atol=2e-6
        @test energy(J,1,conditions) ≈ eigmin(Symmetric(J)) atol=2e-6
        @test energy(J,3,conditions) ≈ LinearAlgebra.tr(J) atol=2e-6
    end

    # Four levels / two pairs: independently assemble the six occupation states.
    L, pairs = 4, 2
    J = Matrix(Diagonal([0.25,0.5,0.75,1.0])) - 0.2ones(L,L)
    states = [s for s in 0:15 if count_ones(s)==pairs]
    H = zeros(6,6)
    for (col,s) in enumerate(states), j in 1:L
        s & (1<<(j-1)) == 0 && continue
        H[col,col] += J[j,j]
        for i in 1:L
            s & (1<<(i-1)) != 0 && continue
            target = xor(s,(1<<(i-1)) | (1<<(j-1)))
            H[findfirst(==(target),states),col] += J[i,j]
        end
    end
    exact = eigmin(Symmetric(H))
    pqg = energy(J,pairs,:PQG)
    pqgt = energy(J,pairs,:PQGT)
    @test pqg <= exact + 2e-6
    @test pqgt >= pqg - 2e-6
    @test pqgt ≈ exact atol=2e-6
end
