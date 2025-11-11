using NCTSSoS
using JuMP
using MosekTools
using Profile

const SOLVER = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)

function chsh_problem(order::Int = 6)
    @ncpolyvar x[1:2]
    @ncpolyvar y[1:2]

    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    pop = polyopt(f; comm_gps = [x, y], is_unipotent = true)

    solver_config = SolverConfig(optimizer = SOLVER; order = order)

    result = cs_nctssos(pop, solver_config; dualize=true)

    return result
end

println("Warming up...")
# Warm-up run
chsh_problem()

println("\nStarting profiling...")
# Profile the actual run
Profile.clear()
@profile chsh_problem(12)

Profile.print( mincount=10000)

