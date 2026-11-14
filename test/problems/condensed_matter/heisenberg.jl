# Heisenberg Model Tests
# Tests ground state energy bounds for Heisenberg spin chains.

using Test, NCTSSoS, JuMP

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

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

    for J2 in J2s
        ham = sum(T(J1 / 4) * op[i] * op[mod1(i + 1, N)] + T(J2 / 4) * op[i] * op[mod1(i + 2, N)] for op in [x, y, z] for i in 1:N)

        pop = polyopt(ham, registry)

        solver_config = SolverConfig(optimizer=SOLVER, order=2, ts_algo=MMD())

        res = cs_nctssos(pop, solver_config)

        higher = cs_nctssos_higher(pop, res, solver_config)
        # Couple opposite-site pairs: the ground state is either two
        # triplets coupled to spin zero, or two singlets. J2 bonds count twice.
        exact_energy = min(-2J1 + J2, -3J2)
        @test is_solved_and_feasible(higher.model)
        @test higher.objective >= res.objective - 1e-6
        @test higher.objective ≈ exact_energy atol = 1e-6
    end
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

    higher = cs_nctssos_higher(pop, res, solver_config)

    # Independent diagonalization of the 64×64 spin-1/2 Hamiltonian gives
    # E0 = -2.562049935181331, attained after one sparse higher step.
    exact_per_site = -0.4270083225302217
    @test is_solved_and_feasible(res.model)
    @test is_solved_and_feasible(higher.model)
    @test res.objective <= higher.objective + 1e-6
    @test higher.objective / N ≈ exact_per_site atol = 1e-6
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

    # Six periodic rows/columns are three-spin triangles. Each triangle has
    # H_triangle = (S_total² - 9/4)/2 >= -3/4, giving an SOS bound of -4.5.
    # Independent diagonalization of the 512×512 spin matrix gives
    # E0 = -3.969001749928519 = N * -0.44100019443650207.
    # A sparse relaxation need not reach E0; test both bound directions and
    # monotonicity, rather than old iteration-specific solver outputs.
    exact_energy = N * -0.44100019443650207
    @test is_solved_and_feasible(res.model)
    @test -4.5 - 1e-6 <= res.objective <= exact_energy + 1e-6
    for _ in 1:4
        higher = cs_nctssos_higher(pop, res, solver_config)
        @test is_solved_and_feasible(higher.model)
        @test res.objective <= higher.objective + 1e-6
        @test -4.5 - 1e-6 <= higher.objective <= exact_energy + 1e-6
        res = higher
    end
end
