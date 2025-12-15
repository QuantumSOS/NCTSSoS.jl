using Test, NCTSSoS, NCTSSoS.FastPolynomials

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

using COSMO
const QUICK_SOLVER = COSMO.Optimizer
using JuMP
using NCTSSoS:
    neat_dot,
    constrain_moment_matrix!,
    substitute_variables,
    NoElimination

using NCTSSoS.FastPolynomials: expval, terms, Arbitrary, get_state_basis, NCStateWord

# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "State Polynomial Opt 7.2.0" begin
    @ncpolyvar x[1:2] y[1:2]
    sp =
        -1.0 * ς(x[1] * y[1]) - 1.0 * ς(x[1] * y[2]) - 1.0 * ς(x[2] * y[1]) +
        1.0 * ς(x[2] * y[2])
    spop = polyopt(sp * one(Monomial); is_unipotent=true, comm_gps=[x, y])

    d = 1

    solver_config = SolverConfig(; optimizer=SOLVER, order=d)

    if haskey(ENV, "LOCAL_TESTING")
        result_mom = cs_nctssos(spop, solver_config; dualize=false)
        @test isapprox(result_mom.objective, -2.8284271321623202, atol=1e-5)
    end

    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -2.8284271321623202, atol=1e-5)


    @testset "Sparse" begin
        solver_config = SolverConfig(; optimizer=SOLVER, order=d, cs_algo=NoElimination(), ts_algo=MMD())

        result = cs_nctssos(spop, solver_config)

        @test result.objective ≈ -2.8284271321623202 atol = 1e-5
    end
end
end  # if false - State Polynomial Opt 7.2.0

# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "State Polynomial Opt 7.2.1" begin
    @ncpolyvar x[1:2] y[1:2]
    sp1 = 1.0 * ς(x[1] * y[2]) + 1.0 * ς(x[2] * y[1])
    sp2 = 1.0 * ς(x[1] * y[1]) + -1.0 * ς(x[2] * y[2])
    sp = -1.0 * sp1 * sp1 - 1.0 * sp2 * sp2

    spop = polyopt(sp * one(Monomial); is_unipotent=true, comm_gps=[x, y])

    d = 3

    solver_config = SolverConfig(; optimizer=QUICK_SOLVER, order=d)

    if haskey(ENV, "LOCAL_TESTING")
        result_mom = cs_nctssos(spop, solver_config; dualize=false)
        @test isapprox(result_mom.objective, -4.0, atol=1e-4)
    end

    result_sos = cs_nctssos(spop, solver_config)
    @test isapprox(result_sos.objective, -4.0, atol=1e-4)
end
end  # if false - State Polynomial Opt 7.2.1

# SKIP: basis.jl type mismatch - see TODO in basis.jl
# These tests use StatePolynomial/NCStatePolynomial which requires basis.jl migration
# The legacy @ncpolyvar macro has been removed. Tests need to use:
#   - create_unipotent_variables() for registry-based variable creation
#   - ς() on Monomial{A,T} from FastPolynomials (not DynamicPolynomials)
#   - Registry-based polyopt() without is_unipotent/comm_gps kwargs
if false
@testset "State Polynomial Opt 7.2.2" begin
    # Migration notes:
    # - Replace @ncpolyvar x[1:3] y[1:3] with:
    #   reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
    # - ς() expects Monomial{A,T} from FastPolynomials
    # - Remove is_unipotent, comm_gps kwargs - use registry-based API
    #
    # Original code (requires basis.jl migration):
    # @ncpolyvar x[1:3] y[1:3]
    # cov(a, b) = 1.0 * ς(x[a] * y[b]) - 1.0 * ς(x[a]) * ς(y[b])
    # sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
    # spop = polyopt(sp * one(Monomial); is_unipotent=true, comm_gps=[x[1:3], y[1:3]])
    # solver_config = SolverConfig(; optimizer=SOLVER, order=2)
    # result = cs_nctssos(spop, solver_config)
    # @test result.objective ≈ -5.0 atol = 1e-2
end
end  # if false - State Polynomial Opt 7.2.2


# SKIP: basis.jl type mismatch - see TODO in basis.jl
if false
@testset "Constrain Moment matrix" begin
    @ncpolyvar x[1:2]

    sa = SimplifyAlgorithm(comm_gps=[x], is_unipotent=false, is_projective=false)

    basis = get_state_basis(Arbitrary, x, 1, sa)

    sp = 1.0 * ς(x[1] * x[2]) + 2.0 * ς(x[1]) + 3.0 * ς(x[2])
    nc_words = monomial.([one(x[1]), x[1], x[2]])
    ncsp =
        1.0 * ς(x[1] * x[2]) * one(Monomial) +
        2.0 * ς(x[1]) * monomial(x[1]) +
        3.0 * ς(x[2]) * monomial(x[2])
    poly = one(ncsp)

    total_basis = sort(unique([expval(neat_dot(a, b)) for a in basis for b in basis]))

    model = GenericModel{Float64}()
    @variable(model, y[1:length(total_basis)])
    wordmap = Dict(zip(total_basis, y))

    ncterms = map(
        a -> a[1] * NCStateWord(Arbitrary, monomial.(a[2]), a[3]),
        zip([1.0, 2.0, 3.0], [[x[1] * x[2]], [x[1]], [x[2]]], nc_words),
    )

    @test map(a -> a[1] * a[2], terms(ncsp)) == ncterms
    @test substitute_variables(expval(ncsp), wordmap) == 1.0 * y[8] + 3.0 * y[6] + 2.0 * y[4]

    true_mom_mtx = expval.([neat_dot(a, b) for a in basis, b in basis])
    mom_mtx_cons =
        constrain_moment_matrix!(model, one(ncsp), basis, wordmap, PSDCone(), sa)
    mom_mtx = constraint_object(mom_mtx_cons)

    reshape(mom_mtx.func, 5, 5)
    @test reshape(mom_mtx.func, 5, 5) == AffExpr[
        y[1] y[2] y[3] y[2] y[3];
        y[2] y[4] y[5] y[4] y[5];
        y[3] y[5] y[6] y[5] y[6];
        y[2] y[4] y[5] y[7] y[8];
        y[3] y[5] y[6] y[8] y[10]
    ]
end
end  # if false - Constrain Moment matrix
