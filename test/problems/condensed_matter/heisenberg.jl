# Heisenberg Model Tests
# =======================
# Tests ground state energy bounds for Heisenberg spin chains.

using Test, NCTSSoS

# Load solver configuration if running standalone
@isdefined(SOLVER) || include(joinpath(dirname(@__FILE__), "..", "..", "standalone_setup.jl"))

@testset "XXX Model (N=6)" begin
    T = ComplexF64
    N = 6
    J1 = 1.0

    registry, (x, y, z) = create_pauli_variables(1:N)

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] for op in [x, y, z] for i in 1:N)

    pop = polyopt(ham, registry)

    solver_config = SolverConfig(optimizer=SOLVER, order=2)

    res = cs_nctssos(pop, solver_config)

    @test res.objective / N ≈ -0.467129 atol = 1e-6
end

@testset "J1 J2 Model (N=4)" begin
    T = ComplexF64
    N = 4
    J1 = 1.0
    J2s = 0.1:0.1:2.1

    registry, (x, y, z) = create_pauli_variables(1:N)

    energy_lower_bounds = zeros(length(J2s))

    for (idx, J2) in enumerate(J2s)
        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        pop = polyopt(ham, registry)

        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)

        res = cs_nctssos_higher(pop, res, solver_config)
        energy_lower_bounds[idx] = res.objective / N
    end

    # Just verify we get reasonable bounds
    @test all(e -> e < 0, energy_lower_bounds)
end

@testset "J1 J2 Model (N=6)" begin
    T = ComplexF64
    N = 6
    J1 = 1.0
    J2 = 0.2

    registry, (x, y, z) = create_pauli_variables(1:N)

    ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

    pop = polyopt(ham, registry)

    solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

    res = cs_nctssos(pop, solver_config)

    res = cs_nctssos_higher(pop, res, solver_config)

    @test res.objective / N ≈ -0.4270083225302217 atol = 1e-6
end

@testset "2D Model (3x3)" begin
    T = ComplexF64
    Nx = 3
    Ny = 3
    N = Nx * Ny
    J1 = 1.0
    J2 = 0.0

    registry, (x, y, z) = create_pauli_variables(1:N)

    LI = LinearIndices((1:Nx, 1:Ny))

    ham = sum(T(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(i, mod1(j + 1, Ny))]] + T(J1 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), j)]] + T(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j + 1, Ny))]] + T(J2 / 4) * op[LI[CartesianIndex(i, j)]] * op[LI[CartesianIndex(mod1(i + 1, Nx), mod1(j - 1, Ny))]] for op in [x, y, z] for i in 1:Nx for j in 1:Ny)

    pop = polyopt(ham, registry)

    solver_config = SolverConfig(optimizer=SOLVER, order=3, cs_algo=MF(), ts_algo=MMD())

    res = cs_nctssos(pop, solver_config)

    res = cs_nctssos_higher(pop, res, solver_config)
    res = cs_nctssos_higher(pop, res, solver_config)
    res = cs_nctssos_higher(pop, res, solver_config)
    @test_broken res.objective ≈ -4.390300714054776 atol = 1e-6
    res = cs_nctssos_higher(pop, res, solver_config)
    @test_broken res.objective ≈ -4.381164563801521 atol = 1e-6
    @test_broken res.objective / N ≈ -0.44100019443650207 atol = 1e-6
end
