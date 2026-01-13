# NCTSSOS Oracle Script: CHSH Bell Inequality (All Formulations)
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_chsh.jl
#
# Problem: CHSH Bell inequality with unipotent constraint (U²=I)
# Variables: x[1:2] (Alice), y[1:2] (Bob)
# Objective: f = x₁y₁ + x₁y₂ + x₂y₁ - x₂y₂
# Expected optimal: 2√2 ≈ 2.8284 (Tsirelson bound)
# Sign convention: Both NCTSSOS and NCTSSoS minimize, so opt ≈ -2.8284
#
# This script generates oracles for THREE formulations:
# 1. Basic NC polynomial (nctssos_first / cs_nctssos_first)
# 2. State polynomial (pstateopt_first) - uses ς() operator
# 3. Trace polynomial (ptraceopt_first) - uses tr() operator

include("oracle_utils.jl")

# PART 1: Basic NC Polynomial Formulation
# Uses nctssos_first / cs_nctssos_first with @ncpolyvar

const CHSH_NC_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1),
    (name="CS", cs="MF", ts=false, order=1),
    (name="TS", cs=false, ts="MD", order=1),
    (name="CS_TS", cs="MF", ts="MD", order=2),
]

const CHSH_NCTSSOS_PARAMS = (
    partition = 2,
    constraint = "unipotent",
)

print_header("CHSH Bell Inequality - NC Polynomial")

@ncpolyvar x[1:2] y[1:2]
obj = -(x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2])
vars = [x; y]

println("Objective: -(x₁y₁ + x₁y₂ + x₂y₁ - x₂y₂)")
println("Constraint: $(CHSH_NCTSSOS_PARAMS.constraint)")
println("Expected: -2√2 ≈ -2.8284")
println()

results_nc = map(CHSH_NC_VARIANTS) do v
    key = "CHSH_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), CS=$(v.cs), TS=$(v.ts))")

    if v.cs == false
        opt, data = nctssos_first([obj], vars, v.order;
            TS=v.ts, partition=CHSH_NCTSSOS_PARAMS.partition,
            constraint=CHSH_NCTSSOS_PARAMS.constraint)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first([obj], vars, v.order;
            TS=v.ts, CS=v.cs, partition=CHSH_NCTSSOS_PARAMS.partition,
            constraint=CHSH_NCTSSOS_PARAMS.constraint)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

# PART 2: State Polynomial Formulation
# Uses pstateopt_first with support/coefficient format
# State polynomial: sp = -ς(x₁y₁) - ς(x₁y₂) - ς(x₂y₁) + ς(x₂y₂)
# Variable mapping: x₁=1, x₂=2, y₁=3, y₂=4

const CHSH_STATE_VARIANTS = [
    (name="State_Dense", ts=false, order=1),
    (name="State_TS", ts="MD", order=1),
]

print_header("CHSH Bell Inequality - State Polynomial")

# Support format: [[vars]] for ς(product of vars)
# ς(x₁y₁) = [[1,3]], ς(x₁y₂) = [[1,4]], etc.
supp_state = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]
coe_state = [-1.0, -1.0, -1.0, 1.0]

println("Objective: -ς(x₁y₁) - ς(x₁y₂) - ς(x₂y₁) + ς(x₂y₂)")
println("Variables: n=4, vargroup=[2;2] (group sizes: Alice=2, Bob=2)")
println("Constraint: unipotent")
println("Expected: -2√2 ≈ -2.8284")
println()

results_state = map(CHSH_STATE_VARIANTS) do v
    key = "CHSH_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")

    opt, data = pstateopt_first(supp_state, coe_state, 4, v.order;
        vargroup=[2, 2],  # Group sizes: 2 vars in party 1, 2 vars in party 2
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# PART 3: Trace Polynomial Formulation
# Uses ptraceopt_first with support/coefficient format
# Trace polynomial: p = -tr(x₁y₁) - tr(x₁y₂) - tr(x₂y₁) + tr(x₂y₂)
# Variable mapping: x₁=1, x₂=2, y₁=3, y₂=4

const CHSH_TRACE_VARIANTS = [
    (name="Trace_Dense", ts=false, order=1),
    (name="Trace_TS", ts="MD", order=1),
]

print_header("CHSH Bell Inequality - Trace Polynomial")

# Support format: [[vars]] for tr(product of vars)
# tr(x₁y₁) = [[1,3]], tr(x₁y₂) = [[1,4]], etc.
supp_trace = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]
coe_trace = [-1.0, -1.0, -1.0, 1.0]

println("Objective: -tr(x₁y₁) - tr(x₁y₂) - tr(x₂y₁) + tr(x₂y₂)")
println("Variables: n=4")
println("Constraint: unipotent")
println("Expected: -2√2 ≈ -2.8284")
println()

results_trace = map(CHSH_TRACE_VARIANTS) do v
    key = "CHSH_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")

    opt, data = ptraceopt_first(supp_trace, coe_trace, 4, v.order;
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# Summary
all_results = vcat(results_nc, results_state, results_trace)
print_summary("CHSH", all_results)
