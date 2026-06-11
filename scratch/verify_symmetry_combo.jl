# Round 2: verify that manual + auto generators combine, and what the union buys.
using NCTSSoS, COSMO

const MOI = NCTSSoS.MOI
const SOLVER = MOI.OptimizerWithAttributes(COSMO.Optimizer, MOI.Silent() => true)

function run_case(name, pop, basis, spec)
    config = SolverConfig(
        optimizer    = SOLVER,
        moment_basis = basis,
        cs_algo      = NoElimination(),
        ts_algo      = NoElimination(),
        symmetry     = spec,
    )
    res = cs_nctssos(pop, config)
    rep = res.symmetry
    println("== ", name, " ==")
    println("  objective          = ", res.objective)
    println("  group_order        = ", rep.group_order)
    println("  psd_block_sizes    = ", rep.psd_block_sizes)
    println("  invariant_moments  = ", rep.invariant_moment_count)
    println("  unique moments     = ", res.n_unique_moment_matrix_elements)
    return res
end

# 2-site
registry, (sx, sy, sz) = create_pauli_variables(1:2)
H = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (sx, sy, sz))
pop = polyopt(H, registry)
basis = [one(sx[1]), sx[1], sx[2], sy[1], sy[2], sz[1], sz[2]]

swap = CliffordSymmetry(:SWAP, 1, 2)
auto = sympleq_symmetry_spec(H)
combo = SymmetrySpec(auto.clifford_generators..., swap)
run_case("2-site auto + manual SWAP", pop, basis, combo)

# 4-site
N = 4
registry4, (sx4, sy4, sz4) = create_pauli_variables(1:N)
H4 = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (sx4, sy4, sz4) for i in 1:N
)
pop4 = polyopt(H4, registry4)
basis4 = [one(sx4[1]); sx4; sy4; sz4]

M4 = typeof(sx4[1])
timgs = Dict{M4,Tuple{Int,M4}}()
rimgs = Dict{M4,Tuple{Int,M4}}()
for i in 1:N
    for op in (sx4, sy4, sz4)
        timgs[op[i]] = (1, op[mod1(i + 1, N)])
        rimgs[op[i]] = (1, op[N + 1 - i])
    end
end
translation4 = CliffordSymmetry(timgs)
reflection4 = CliffordSymmetry(rimgs)

auto4 = sympleq_symmetry_spec(H4)
combo4 = SymmetrySpec(auto4.clifford_generators..., translation4, reflection4)
run_case("4-site auto + translation + reflection", pop4, basis4, combo4)

println("ALL CASES DONE")
