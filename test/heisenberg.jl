using NCTSSoS, NCTSSoS.FastPolynomials, Test
using JuMP

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

# Heisenberg model tests use PauliAlgebra via the pauli_algebra() helper.
# PauliAlgebra automatically handles simplification rules (sigma^2 = I, cyclic products).

@testset "J1 J2 Model (N=4)" begin
    T = ComplexF64
    N = 4
    J1 = 1.0
    J2s = 0.1:0.1:2.1

    # Create Pauli algebra system
    sys = pauli_algebra(N)
    x, y, z = sys.sx, sys.sy, sys.sz

    energy_lower_bounds = zeros(length(J2s))

    for (idx, J2) in enumerate(J2s)
        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        # Use new registry-based API - Pauli simplification is automatic
        pop = polyopt(ham, sys.registry)

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

    # Create Pauli algebra system
    sys = pauli_algebra(N)
    x, y, z = sys.sx, sys.sy, sys.sz

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

    # Use new registry-based API - Pauli simplification is automatic
    pop = polyopt(ham, sys.registry)

    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    res = cs_nctssos(pop, solver_config)

    @test res.objective / N ≈ -0.467129 atol = 1e-6
end


@testset "J1 J2 Model (N=6)" begin
    T = ComplexF64
    N = 6
    J1 = 1.0
    J2 = 0.2

    # Create Pauli algebra system
    sys = pauli_algebra(N)
    x, y, z = sys.sx, sys.sy, sys.sz

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

    # Use new registry-based API - Pauli simplification is automatic
    pop = polyopt(ham, sys.registry)

    solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

    res = cs_nctssos(pop, solver_config)

    res = cs_nctssos_higher(pop, res, solver_config)

    @test res.objective / N ≈ -0.4270083225302217 atol = 1e-6
end

# 2D Model commented out - too slow for regular testing
# @testset "2D Model" begin
#     T = ComplexF64
#     Nx = 3
#     Ny = 3
#     N = Nx * Ny
#     J1 = 1.0
#     J2 = 0.0
#
#     sys = pauli_algebra(N)
#     x, y, z = sys.sx, sys.sy, sys.sz
#
#     LI = LinearIndices((1:Nx, 1:Ny))
#
#     ham = sum(T(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + T(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + T(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + T(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [x, y, z] for i in 1:Nx for j in 1:Ny)
#
#     pop = polyopt(ham, sys.registry)
#
#     solver_config = SolverConfig(optimizer=SOLVER, order=3, cs_algo=MF(), ts_algo=MMD())
#
#     res = cs_nctssos(pop, solver_config)
#
#     res = cs_nctssos_higher(pop, res, solver_config)
#     res = cs_nctssos_higher(pop, res, solver_config)
#
#     res = cs_nctssos_higher(pop, res, solver_config) # -4.390300714054776 = -0.4878111904505307 * N
#     res = cs_nctssos_higher(pop, res, solver_config) # -4.381164563801521 = -0.48679606264461345
#
#     @test res.objective / N ≈ -0.44100019443650207 atol = 1e-6
# end
