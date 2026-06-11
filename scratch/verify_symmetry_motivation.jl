# Verification script for the pauli_clifford_symmetry docs rework.
# Runs license-free (COSMO). Prints group orders, PSD block sizes,
# invariant moment counts, and objectives for every claim the new docs make.

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
    if rep === nothing
        println("  moment_matrix_sizes = ", res.moment_matrix_sizes)
        println("  unique moments      = ", res.n_unique_moment_matrix_elements)
    else
        println("  group_order        = ", rep.group_order)
        println("  psd_block_sizes    = ", rep.psd_block_sizes)
        println("  invariant_moments  = ", rep.invariant_moment_count)
        println("  basis half/full    = ", (rep.basis_half_size, rep.basis_full_size))
        println("  unique moments     = ", res.n_unique_moment_matrix_elements)
    end
    return res
end

# ---------------------------------------------------------------- 2-site
println("################ 2-site Heisenberg ################")
registry, (sx, sy, sz) = create_pauli_variables(1:2)
H = sum(ComplexF64(1 / 4) * op[1] * op[2] for op in (sx, sy, sz))
pop = polyopt(H, registry)
basis = [one(sx[1]), sx[1], sx[2], sy[1], sy[2], sz[1], sz[2]]

run_case("dense", pop, basis, nothing)

swap = CliffordSymmetry(:SWAP, 1, 2)
run_case("manual SWAP", pop, basis, SymmetrySpec(swap))

# global Hadamard: X<->Z, Y->-Y on BOTH sites (axis exchange)
M = typeof(sx[1])
gH = CliffordSymmetry(Dict{M,Tuple{Int,M}}(
    sx[1] => (1, sz[1]), sz[1] => (1, sx[1]), sy[1] => (-1, sy[1]),
    sx[2] => (1, sz[2]), sz[2] => (1, sx[2]), sy[2] => (-1, sy[2]),
))
# global S: X->Y, Y->-X on BOTH sites (axis rotation about z)
gS = CliffordSymmetry(Dict{M,Tuple{Int,M}}(
    sx[1] => (1, sy[1]), sy[1] => (-1, sx[1]),
    sx[2] => (1, sy[2]), sy[2] => (-1, sx[2]),
))

run_case("manual SWAP + globalH", pop, basis, SymmetrySpec(swap, gH))
run_case("manual SWAP + globalH + globalS", pop, basis, SymmetrySpec(swap, gH, gS))
run_case("manual globalH + globalS (no swap)", pop, basis, SymmetrySpec(gH, gS))

auto = sympleq_symmetry_spec(H)
println("auto generators found: ", length(auto.clifford_generators))
run_case("auto sympleq", pop, basis, auto)

# ---------------------------------------------------------------- 4-site
println()
println("################ 4-site Heisenberg ring ################")
N = 4
registry4, (sx4, sy4, sz4) = create_pauli_variables(1:N)
H4 = sum(
    ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
    for op in (sx4, sy4, sz4) for i in 1:N
)
pop4 = polyopt(H4, registry4)
basis4 = [one(sx4[1]); sx4; sy4; sz4]

run_case("dense 4-site", pop4, basis4, nothing)

# manual lattice translation i -> i+1 (mod 4)
M4 = typeof(sx4[1])
timgs = Dict{M4,Tuple{Int,M4}}()
for i in 1:N
    j = mod1(i + 1, N)
    for op in (sx4, sy4, sz4)
        timgs[op[i]] = (1, op[j])
    end
end
translation4 = CliffordSymmetry(timgs)
run_case("manual translation", pop4, basis4, SymmetrySpec(translation4))

# manual reflection i -> N+1-i
rimgs = Dict{M4,Tuple{Int,M4}}()
for i in 1:N
    j = N + 1 - i
    for op in (sx4, sy4, sz4)
        rimgs[op[i]] = (1, op[j])
    end
end
reflection4 = CliffordSymmetry(rimgs)
run_case("manual translation + reflection", pop4, basis4,
    SymmetrySpec(translation4, reflection4))

# global S on all 4 sites — is the axis rotation invariant here too?
simgs = Dict{M4,Tuple{Int,M4}}()
for i in 1:N
    simgs[sx4[i]] = (1, sy4[i])
    simgs[sy4[i]] = (-1, sx4[i])
end
gS4 = CliffordSymmetry(simgs)
try
    run_case("manual translation + reflection + globalS", pop4, basis4,
        SymmetrySpec(translation4, reflection4, gS4))
catch err
    println("  translation+reflection+globalS FAILED: ", sprint(showerror, err))
end

auto4 = sympleq_symmetry_spec(H4)
println("auto generators found (4-site): ", length(auto4.clifford_generators))
run_case("auto sympleq 4-site", pop4, basis4, auto4)

# combine auto with the manual global S, in case SympleQ misses it
try
    combo = SymmetrySpec(auto4.clifford_generators..., gS4)
    run_case("auto + globalS 4-site", pop4, basis4, combo)
catch err
    println("  auto+globalS FAILED: ", sprint(showerror, err))
end

println()
println("ALL CASES DONE")
