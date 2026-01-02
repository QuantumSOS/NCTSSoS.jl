# NCTSSOS Oracle Script: CHSH Bell Inequality
# =============================================
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_chsh.jl
#
# Problem: CHSH Bell inequality with unipotent constraint (U²=I)
# Variables: x[1:2] (Alice), y[1:2] (Bob)
# Objective: f = x₁y₁ + x₁y₂ + x₂y₁ - x₂y₂
# Expected optimal: 2√2 ≈ 2.8284 (Tsirelson bound)
# Sign convention: Both NCTSSOS and NCTSSoS minimize, so opt ≈ -2.8284

include("oracle_utils.jl")

# Sparsity variants
const CHSH_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=1),
    (name="CS", cs="MF", ts=false, order=1),
    (name="TS", cs=false, ts="MD", order=1),
    (name="CS_TS", cs="MF", ts="MD", order=2),
]

# NCTSSOS parameters
const CHSH_NCTSSOS_PARAMS = (
    partition = 2,
    constraint = "unipotent",
)

print_header("CHSH Bell Inequality")

# Problem setup
@ncpolyvar x[1:2] y[1:2]
obj = -(x[1]*y[1] + x[1]*y[2] + x[2]*y[1] - x[2]*y[2])
vars = [x; y]

println("Objective: -(x₁y₁ + x₁y₂ + x₂y₁ - x₂y₂)")
println("Constraint: $(CHSH_NCTSSOS_PARAMS.constraint)")
println("Expected: -2√2 ≈ -2.8284")
println()

results = map(CHSH_VARIANTS) do v
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

print_summary("CHSH", results)
