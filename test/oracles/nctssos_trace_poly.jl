# NCTSSOS Oracle Script: Trace Polynomial Optimization
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_trace_poly.jl
#
# Problems from Section 6 of the NCTSSOS paper.
# Trace polynomials use the tr() operator.
# API: ptraceopt_first for pure trace polynomial optimization
#      Support format: [[w1], [w2]] means tr(w1)·tr(w2)

include("oracle_utils.jl")

# Example 6.1: Projector Algebra (Toy Example)
# Variables: x[1:3] (projector, P²=P)
# Objective: p = tr(x₁x₂x₃) + tr(x₁x₂)·tr(x₃)
# Expected: ≈ -0.0467 (order 2), ≈ -0.03125 (order 3)

const TRACE_6_1_VARIANTS = [
    (name="Dense", ts=false, order=2),
    (name="Dense", ts=false, order=3),
    (name="TS", ts="MD", order=2),
    (name="TS", ts="MD", order=3),
]

print_header("Trace Polynomial 6.1 (Projector Algebra)")

# tr(x₁x₂x₃) + tr(x₁x₂)·tr(x₃)
supp_6_1 = [
    [[1,2,3]],     # tr(x₁x₂x₃)
    [[1,2], [3]],  # tr(x₁x₂)·tr(x₃)
]
coe_6_1 = [1.0, 1.0]

println("# Variables: x[1:3] (projector, P²=P)")
println("# Objective: tr(x₁x₂x₃) + tr(x₁x₂)·tr(x₃)")
println("# Expected: ≈ -0.0467 (d=2), ≈ -0.03125 (d=3)")
println()

results_6_1 = map(TRACE_6_1_VARIANTS) do v
    key = "Trace_6_1_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = ptraceopt_first(supp_6_1, coe_6_1, 3, v.order;
        TS=v.ts == false ? false : v.ts,
        constraint="projection")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# Example 6.2.0: CHSH Trace Polynomial
# Variables: x[1:2], y[1:2] (unipotent, U²=I)
# Variable mapping: x1=1, x2=2, y1=3, y2=4
# Objective: p = -tr(x₁y₁) - tr(x₁y₂) - tr(x₂y₁) + tr(x₂y₂)
# Expected: -2√2 ≈ -2.8284

const TRACE_6_2_0_VARIANTS = [
    (name="Dense", ts=false, order=1),
    (name="TS", ts="MD", order=1),
]

print_header("Trace Polynomial 6.2.0 (CHSH)")

supp_6_2_0 = [
    [[1,3]],  # tr(x₁y₁)
    [[1,4]],  # tr(x₁y₂)
    [[2,3]],  # tr(x₂y₁)
    [[2,4]],  # tr(x₂y₂)
]
coe_6_2_0 = [-1.0, -1.0, -1.0, 1.0]

println("# Variables: x[1:2]={1,2}, y[1:2]={3,4}")
println("# Objective: -tr(x₁y₁) - tr(x₁y₂) - tr(x₂y₁) + tr(x₂y₂)")
println("# Expected: -2√2 ≈ -2.8284")
println()

results_6_2_0 = map(TRACE_6_2_0_VARIANTS) do v
    key = "Trace_6_2_0_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = ptraceopt_first(supp_6_2_0, coe_6_2_0, 4, v.order;
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# Example 6.2.1: Squared Trace Expressions
# Variables: x[1:2], y[1:2] (unipotent)
# Variable mapping: x1=1, x2=2, y1=3, y2=4
# Objective: -(tr(x₁y₂) + tr(x₂y₁))² - (tr(x₁y₁) - tr(x₂y₂))²
# Expected: -4.0

const TRACE_6_2_1_VARIANTS = [
    (name="Dense", ts=false, order=2),
    (name="TS", ts="MD", order=2),
]

print_header("Trace Polynomial 6.2.1 (Squared Traces)")

# Same structure as state 7.2.1 but with tr() instead of ς()
supp_6_2_1 = [
    [[1,4], [1,4]],  # tr(x₁y₂)²
    [[2,3], [2,3]],  # tr(x₂y₁)²
    [[1,4], [2,3]],  # tr(x₁y₂)tr(x₂y₁)
    [[1,3], [1,3]],  # tr(x₁y₁)²
    [[2,4], [2,4]],  # tr(x₂y₂)²
    [[1,3], [2,4]],  # tr(x₁y₁)tr(x₂y₂)
]
coe_6_2_1 = [-1.0, -1.0, -2.0, -1.0, -1.0, 2.0]

println("# Variables: x[1:2]={1,2}, y[1:2]={3,4}")
println("# Objective: -(tr(x₁y₂) + tr(x₂y₁))² - (tr(x₁y₁) - tr(x₂y₂))²")
println("# Expected: -4.0")
println()

results_6_2_1 = map(TRACE_6_2_1_VARIANTS) do v
    key = "Trace_6_2_1_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = ptraceopt_first(supp_6_2_1, coe_6_2_1, 4, v.order;
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# Example 6.2.2: Covariance Trace Polynomial
# Variables: x[1:3], y[1:3] (unipotent)
# Variable mapping: x1=1, x2=2, x3=3, y1=4, y2=5, y3=6
# cov(a,b) = tr(xₐyᵦ) - tr(xₐ)tr(yᵦ)
# Objective: -cov(1,1) - cov(1,2) - cov(1,3) - cov(2,1) - cov(2,2) + cov(2,3) - cov(3,1) + cov(3,2)
# Expected: -5.0

const TRACE_6_2_2_VARIANTS = [
    (name="Dense", ts=false, order=2),
    (name="TS", ts="MD", order=2),
]

print_header("Trace Polynomial 6.2.2 (Covariance)")

# Same structure as state 7.2.2
supp_6_2_2 = [
    [[1,4]],     # tr(x₁y₁)
    [[1], [4]],  # tr(x₁)tr(y₁)
    [[1,5]],     # tr(x₁y₂)
    [[1], [5]],  # tr(x₁)tr(y₂)
    [[1,6]],     # tr(x₁y₃)
    [[1], [6]],  # tr(x₁)tr(y₃)
    [[2,4]],     # tr(x₂y₁)
    [[2], [4]],  # tr(x₂)tr(y₁)
    [[2,5]],     # tr(x₂y₂)
    [[2], [5]],  # tr(x₂)tr(y₂)
    [[2,6]],     # tr(x₂y₃) with negative cov
    [[2], [6]],  # tr(x₂)tr(y₃)
    [[3,4]],     # tr(x₃y₁)
    [[3], [4]],  # tr(x₃)tr(y₁)
    [[3,5]],     # tr(x₃y₂) with negative cov
    [[3], [5]],  # tr(x₃)tr(y₂)
]
coe_6_2_2 = [
    -1.0, 1.0,   # cov(1,1)+
    -1.0, 1.0,   # cov(1,2)+
    -1.0, 1.0,   # cov(1,3)+
    -1.0, 1.0,   # cov(2,1)+
    -1.0, 1.0,   # cov(2,2)+
     1.0, -1.0,  # cov(2,3)-
    -1.0, 1.0,   # cov(3,1)+
     1.0, -1.0,  # cov(3,2)-
]

println("# Variables: x[1:3]={1,2,3}, y[1:3]={4,5,6}")
println("# Objective: Σ±cov(i,j) where cov(a,b) = tr(xₐyᵦ) - tr(xₐ)tr(yᵦ)")
println("# Expected: -5.0")
println()

results_6_2_2 = map(TRACE_6_2_2_VARIANTS) do v
    key = "Trace_6_2_2_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = ptraceopt_first(supp_6_2_2, coe_6_2_2, 6, v.order;
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# Summary
all_results = vcat(results_6_1, results_6_2_0, results_6_2_1, results_6_2_2)
print_summary("TRACE_POLY", all_results)
