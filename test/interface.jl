using Test, NCTSSoS

using Clarabel
if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

if haskey(ENV, "LOCAL_TESTING") 
    @testset "1D Transverse Field Ising Model" begin
        N = 3
        @ncpolyvar x[1:N] y[1:N] z[1:N]

        J = 1.0
        h = 2.0
        for (periodic, true_ans) in zip((true, false), (-1.0175918, -1.0104160))
            ham = sum(-complex(J / 4) * z[i] * z[mod1(i + 1, N)] for i in 1:(periodic ? N : N - 1)) + sum(-h / 2 * x[i] for i in 1:N)

            eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

            pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

            solver_config = SolverConfig(optimizer=SOLVER, order=2)

            res = cs_nctssos(pop, solver_config)
            @test res.objective / N ≈ true_ans atol = 1e-6
        end
    end
end

# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "Naive Example" begin
    N = 1
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    ham = sum(ComplexF64(1/2) * op[1] for op in [x,y,z])

    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    solver_config = SolverConfig(optimizer=SOLVER, order=1)

    @test_throws ErrorException cs_nctssos(pop, solver_config; dualize=false)
    res = cs_nctssos(pop, solver_config)
    @test res.objective ≈ -0.8660254037844387 atol = 1e-6
end
end  # if false - Naive Example

# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "Naive Example 2" begin
    N = 1
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    ham = one(ComplexF64) * x[1] * y[1] + one(ComplexF64) * y[1] * x[1] 

    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    solver_config = SolverConfig(optimizer=SOLVER, order=3)

    res = cs_nctssos(pop, solver_config)
    @test res.objective ≈ -0.0 atol = 1e-6
end
end  # if false - Naive Example 2

if haskey(ENV, "LOCAL_TESTING")
    @testset "1D Heisenberg Chain" begin
        N = 6
        @ncpolyvar x[1:N] y[1:N] z[1:N]

        ham = sum(ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2)

        res = cs_nctssos(pop, solver_config)

        @test res.objective / N ≈ -0.467129 atol = 1e-6
    end

    @testset "Example" begin
        @ncpolyvar x[1:3]
        @ncpolyvar y[1:3]
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
            x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]  # objective function

        pop = polyopt(-f; comm_gps=[x, y], is_projective=true)

        for (cs_algo, ts_algo, ans) in zip([NoElimination(), MF(), MF()],
            [NoElimination(), MMD(), MaximalElimination()],
            [-0.2508755573198166, -0.9999999892255513, -0.2512780696727863])
            solver_config = SolverConfig(optimizer=SOLVER; order=3, cs_algo=cs_algo, ts_algo=ts_algo)
            result = cs_nctssos(pop, solver_config)
            @test isapprox(result.objective, ans; atol=1e-5)
        end
    end
end

# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "Majumdar Gosh Model" begin
    num_sites = 6
    J1_interactions =
        unique!([tuple(sort([i, mod1(i + 1, num_sites)])...) for i = 1:num_sites])
    J2_interactions =
        unique!([tuple(sort([i, mod1(i + 2, num_sites)])...) for i = 1:num_sites])

    J1 = 2.0
    J2 = 1.0

    true_ans = -num_sites / 4 * 6

    ij2idx_dict = Dict(
        zip(
            [(i, j) for i in 1:num_sites, j in 1:num_sites if j > i],
            1:(num_sites*(num_sites-1)÷2),
        ),
    )
    @ncpolyvar hij[1:(num_sites*(num_sites-1)÷2)]

    objective = (
        sum([J1 * hij[ij2idx_dict[(i, j)]] for (i, j) in J1_interactions]) + sum([J2 * hij[ij2idx_dict[(i, j)]] for (i, j) in J2_interactions])
    )

    gs = unique!([
        (
            hij[ij2idx_dict[tuple(sort([i, j])...)]] *
            hij[ij2idx_dict[tuple(sort([j, k])...)]] +
            hij[ij2idx_dict[tuple(sort([j, k])...)]] *
            hij[ij2idx_dict[tuple(sort([i, j])...)]] -
            0.5 * (
                hij[ij2idx_dict[tuple(sort([i, j])...)]] +
                hij[ij2idx_dict[tuple(sort([j, k])...)]] -
                hij[ij2idx_dict[tuple(sort([i, k])...)]]
            )
        ) for i in 1:num_sites, j in 1:num_sites, k in 1:num_sites if
        (i != j && j != k && i != k)
    ])

    pop = polyopt(
        -objective;
        eq_constraints=gs,
        is_projective=true,
    )

    solver_config = SolverConfig(optimizer=SOLVER; order=1)

    result = cs_nctssos(pop, solver_config)

    @test isapprox(result.objective, true_ans; atol=1e-4)
end
end  # if false - Majumdar Gosh Model

# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "Problem Creation Interface" begin
    n = 2
    @ncpolyvar x[1:n]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0
    # change struct name
    pop = polyopt(f; ineq_constraints=[g], eq_constraints=[h1])

    solver_config = SolverConfig(
        optimizer=SOLVER;
        order=2,
        cs_algo=MF(),
        ts_algo=MMD(),
    )

    result = cs_nctssos(pop, solver_config)
    @test isapprox(result.objective, -1.0; atol=1e-4)

    result_higher = cs_nctssos_higher(pop, result, solver_config)
    @test isapprox(result.objective, result_higher.objective; atol=1e-4)
end
end  # if false - Problem Creation Interface


# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "README Example Unconstrained" begin
    @ncpolyvar x[1:3]
    f =
        1.0 +
        x[1]^4 +
        x[2]^4 +
        x[3]^4 +
        x[1] * x[2] +
        x[2] * x[1] +
        x[2] * x[3] +
        x[3] * x[2]

    pop = polyopt(f)

    solver_config_dense = SolverConfig(optimizer=SOLVER)

    result_dense = cs_nctssos(pop, solver_config_dense)

    result_cs =
        cs_nctssos(pop, SolverConfig(optimizer=SOLVER; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(
        pop,
        SolverConfig(optimizer=SOLVER; cs_algo=MF(), ts_algo=MMD()),
    )

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)
end
end  # if false - README Example Unconstrained

# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "README Example Constrained" begin
    @ncpolyvar x[1:2]
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h1 = x[1] * x[2] + x[2] * x[1] - 2.0

    pop = polyopt(f; ineq_constraints=[g], eq_constraints=[h1])

    result_dense = cs_nctssos(pop, SolverConfig(optimizer=SOLVER))

    result_cs =
        cs_nctssos(pop, SolverConfig(optimizer=SOLVER; cs_algo=MF()))

    @test isapprox(result_dense.objective, result_cs.objective, atol=1e-4)

    result_cs_ts = cs_nctssos(
        pop,
        SolverConfig(optimizer=SOLVER; cs_algo=MF(), ts_algo=MMD()),
    )

    @test isapprox(result_cs.objective, result_cs_ts.objective, atol=1e-4)

    result_cs_ts_higher = cs_nctssos_higher(
        pop,
        result_cs_ts,
        SolverConfig(optimizer=SOLVER; cs_algo=MF(), ts_algo=MMD()),
    )

    @test isapprox(result_dense.objective, result_cs_ts_higher.objective, atol=1e-4)
end
end  # if false - README Example Constrained
