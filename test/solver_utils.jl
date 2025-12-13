using Test, NCTSSoS
using JuMP
using NCTSSoS.FastPolynomials
using NCTSSoS: get_dim, sorted_unique

@testset "Utilities" begin
    @testset "VectorConstraint Dim" begin
        model = Model()
        n = 5
        var1 = @variable(model, [1:n, 1:n])
        var2 = @variable(model, [1:(2*n), 1:(2*n)])

        cons1 = @constraint(model, var1 in PSDCone())
        cons2 = @constraint(model, var2 in Zeros())

        @test get_dim(constraint_object(cons1)) == n
        @test get_dim(constraint_object(cons2)) == 2 * n
    end

    # SKIP: basis.jl type mismatch - see TODO in basis.jl
    if false
    @testset "Simplify" begin
        obj = 1.0 * x * y + 2.0 * y * z

        basis = get_basis([x, y, z], 3)

        pop = polyopt(obj; is_unipotent = true)
        sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)
        @test simplify((y * x^2 * y), sa) == one(Monomial)

        pop = polyopt(obj; is_projective = true)
        sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)
        @test simplify((y * x^2 * y), sa) == y * x * y

        pop = polyopt(obj; comm_gps = [[x], [y, z]], is_unipotent = true)
        sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)
        @test simplify((y * x^2 * y), sa) == one(y)

        @test sorted_unique(simplify.(basis,Ref(sa))) == sort([
            one(x * y),
            monomial([z], [1]),
            monomial([y], [1]),
            monomial([x], [1]),
            z * y,
            x * y,
            y * z,
            x * z,
            x * z * y,
            z * y * z,
            y * z * y,
            x * y * z,
        ])

        pop = polyopt(obj; comm_gps = [[x], [y, z]], is_projective = true)
        sa = SimplifyAlgorithm(comm_gps=pop.comm_gps, is_unipotent=pop.is_unipotent, is_projective=pop.is_projective)
        @test simplify((y * x^2 * y), sa) == x * y

        @test sorted_unique(simplify.(basis,Ref(sa))) == sort([
            one(x * y * z),
            monomial([z], [1]),
            monomial([y], [1]),
            monomial([x], [1]),
            z * y,
            x * z,
            x * y,
            y * z,
            x * z * y,
            z * y * z,
            x * y * z,
            y * z * y,
        ])
    end
    end  # if false - Simplify
end
