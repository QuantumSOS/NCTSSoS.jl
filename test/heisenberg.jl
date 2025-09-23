using NCTSSoS, NCTSSoS.FastPolynomials, Test
using JuMP

if haskey(ENV, "LOCAL_TESTING") 
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

@testset "J1 J2 Model" begin
    T = ComplexF64
	N = 4
    J1 = 1.0
    J2s = 0.1:0.1:2.1
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    energy_lower_bounds = zeros(length(J2s))

    for (idx, J2) in enumerate(J2s)
        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

        pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)

        res = cs_nctssos_higher(pop, res, solver_config)
        energy_lower_bounds[idx] = res.objective / N
    end
    for val in energy_lower_bounds
        println(val, ",")
    end
end



@testset "XXX Model" begin
    T = ComplexF64
    N = 6
    J1 = 1.0
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    res = cs_nctssos(pop, solver_config)

    @test res.objective / N ≈ -0.467129 atol = 1e-6
end


@testset "J1 J2 Model" begin
    T = ComplexF64
	N = 6
    J1 = 1.0
    J2 = 0.2
    @ncpolyvar x[1:N] y[1:N] z[1:N]

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

    eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

    pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

    solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

    res = cs_nctssos(pop, solver_config)

    res = cs_nctssos_higher(pop, res, solver_config)

    @test res.objective / N ≈ -0.4270083225302217 atol = 1e-6
end

# @testset "2D Model" begin
#     T = ComplexF64
#     Nx = 3
#     Ny = 3
#     N = Nx * Ny
#     J1 = 1.0
#     J2 = 0.0
#     @ncpolyvar x[1:N] y[1:N] z[1:N]

#     LI = LinearIndices((1:Nx, 1:Ny))

#     ham = sum(T(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + T(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + T(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + T(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [x, y, z] for i in 1:Nx for j in 1:Ny)

#     eq_cons = reduce(vcat, [[x[i] * y[i] - im * z[i], y[i] * x[i] + im * z[i], y[i] * z[i] - im * x[i], z[i] * y[i] + im * x[i], z[i] * x[i] - im * y[i], x[i] * z[i] + im * y[i]] for i in 1:N])

#     pop = cpolyopt(ham; eq_constraints=eq_cons, comm_gps=[[x[i], y[i], z[i]] for i in 1:N], is_unipotent=true)

#     solver_config = SolverConfig(optimizer=SOLVER, order=3, cs_algo=MF(), ts_algo=MMD())

#     res = cs_nctssos(pop, solver_config)

#     res = cs_nctssos_higher(pop, res, solver_config)
#     res = cs_nctssos_higher(pop, res, solver_config)

#     res = cs_nctssos_higher(pop, res, solver_config) # -4.390300714054776 = -0.4878111904505307 * N
#     res = cs_nctssos_higher(pop, res, solver_config) # -4.381164563801521 = -0.48679606264461345 

#     @test res.objective / N ≈ -0.44100019443650207 atol = 1e-6
# end


