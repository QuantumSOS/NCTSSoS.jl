# NCTSSOS Oracle Script: State Polynomial Optimization (Extended)
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_state_poly_extended.jl
#
# Problems from Section 7.2 of the NCTSSOS paper.
# State polynomials with product terms (mixword support format).
# API: pstateopt_first with multi-word support entries [[w1], [w2]] for ς(w1)ς(w2)

include("oracle_utils.jl")

# 7.2.1: Squared Expectations
# Variables: x[1:2], y[1:2] (unipotent, U²=I)
# Variable mapping: x1=1, x2=2, y1=3, y2=4
# Objective: sp = -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²
#          = -ς(x₁y₂)² - ς(x₂y₁)² - 2ς(x₁y₂)ς(x₂y₁) - ς(x₁y₁)² - ς(x₂y₂)² + 2ς(x₁y₁)ς(x₂y₂)
# Expected: -4.0 (at order=3)
#
# Support format: [[word1], [word2]] means ς(word1)·ς(word2)
# [[word]] means ς(word)² when repeated, or just ς(word) if single entry

const STATE_7_2_1_VARIANTS = [
    (name="Dense", ts=false, order=3),
    (name="TS", ts="MD", order=3),
]

print_header("State Polynomial 7.2.1 (Squared Expectations)")

# sp = -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²
# Expand: -ς(x₁y₂)² - ς(x₂y₁)² - 2ς(x₁y₂)ς(x₂y₁) - ς(x₁y₁)² - ς(x₂y₂)² + 2ς(x₁y₁)ς(x₂y₂)
# Variable mapping: x1=1, x2=2, y1=3, y2=4
# x₁y₂ = [1,4], x₂y₁ = [2,3], x₁y₁ = [1,3], x₂y₂ = [2,4]
supp_7_2_1 = [
    [[1,4], [1,4]],  # ς(x₁y₂)²
    [[2,3], [2,3]],  # ς(x₂y₁)²
    [[1,4], [2,3]],  # ς(x₁y₂)ς(x₂y₁)
    [[1,3], [1,3]],  # ς(x₁y₁)²
    [[2,4], [2,4]],  # ς(x₂y₂)²
    [[1,3], [2,4]],  # ς(x₁y₁)ς(x₂y₂)
]
coe_7_2_1 = [-1.0, -1.0, -2.0, -1.0, -1.0, 2.0]

println("# Variables: x[1:2]={1,2}, y[1:2]={3,4}")
println("# Objective: -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²")
println("# Expected: -4.0")
println()

results_7_2_1 = map(STATE_7_2_1_VARIANTS) do v
    key = "State_7_2_1_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = pstateopt_first(supp_7_2_1, coe_7_2_1, 4, v.order;
        vargroup=[2, 2],
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# 7.2.2: Covariance Expression
# Variables: x[1:3], y[1:3] (unipotent)
# Variable mapping: x1=1, x2=2, x3=3, y1=4, y2=5, y3=6
# cov(a,b) = ς(xₐyᵦ) - ς(xₐ)ς(yᵦ)
# Objective: sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
# Expected: -5.0

const STATE_7_2_2_VARIANTS = [
    (name="Dense", ts=false, order=2),
    (name="TS", ts="MD", order=2),
]

print_header("State Polynomial 7.2.2 (Covariance)")

# Build support for cov(a,b) = ς(xₐyᵦ) - ς(xₐ)ς(yᵦ)
# Each cov term has 2 supports: [[xa, yb]] for ς(xₐyᵦ) and [[xa], [yb]] for ς(xₐ)ς(yᵦ)
# Signs: cov(1,1)+, cov(1,2)+, cov(1,3)+, cov(2,1)+, cov(2,2)+, cov(2,3)-, cov(3,1)+, cov(3,2)-
supp_7_2_2 = [
    [[1,4]],     # ς(x₁y₁)
    [[1], [4]],  # ς(x₁)ς(y₁)
    [[1,5]],     # ς(x₁y₂)
    [[1], [5]],  # ς(x₁)ς(y₂)
    [[1,6]],     # ς(x₁y₃)
    [[1], [6]],  # ς(x₁)ς(y₃)
    [[2,4]],     # ς(x₂y₁)
    [[2], [4]],  # ς(x₂)ς(y₁)
    [[2,5]],     # ς(x₂y₂)
    [[2], [5]],  # ς(x₂)ς(y₂)
    [[2,6]],     # ς(x₂y₃) with negative cov
    [[2], [6]],  # ς(x₂)ς(y₃)
    [[3,4]],     # ς(x₃y₁)
    [[3], [4]],  # ς(x₃)ς(y₁)
    [[3,5]],     # ς(x₃y₂) with negative cov
    [[3], [5]],  # ς(x₃)ς(y₂)
]
# For cov+: coef of ς(xy) is -1, coef of ς(x)ς(y) is +1 (minimize, so -cov)
# For cov-: coef of ς(xy) is +1, coef of ς(x)ς(y) is -1
# We minimize, so objective = -sp = -Σ±cov = Σ∓cov
# cov(1,1)+ -> -ς(x₁y₁) + ς(x₁)ς(y₁)
# cov(2,3)- -> +ς(x₂y₃) - ς(x₂)ς(y₃)
coe_7_2_2 = [
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
println("# Objective: Σ±cov(i,j) where cov(a,b) = ς(xₐyᵦ) - ς(xₐ)ς(yᵦ)")
println("# Expected: -5.0")
println()

results_7_2_2 = map(STATE_7_2_2_VARIANTS) do v
    key = "State_7_2_2_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = pstateopt_first(supp_7_2_2, coe_7_2_2, 6, v.order;
        vargroup=[3, 3],
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# 7.2.3: Mixed State Polynomial
# Variables: x[1:2], y[1:2] (unipotent)
# Variable mapping: x1=1, x2=2, y1=3, y2=4
# From NCTSSOS stateopt.jl example
# Objective has: linear terms, mixed monomials, products of expectations, squared expectations
# Expected: -3.5114802

const STATE_7_2_3_VARIANTS = [
    (name="Dense", ts=false, order=2),
    (name="TS", ts="MD", order=2),
]

print_header("State Polynomial 7.2.3 (Mixed)")

# From NCTSSOS stateopt.jl Example 7.2.3:
# supp = [[[2]], [[3]], [[4]], [[1;3]], [[2;3]], [[1;4]], [[2;4]], [[1], [3]],
# [[2], [3]], [[2], [4]], [[1], [1]], [[4], [4]]]
# coe = -[1; 1; 1; -1; 1; 1; 1; -1; -1; -1; -1; -1]
supp_7_2_3 = [
    [[2]],        # ς(x₂)
    [[3]],        # ς(y₁)
    [[4]],        # ς(y₂)
    [[1,3]],      # ς(x₁y₁)
    [[2,3]],      # ς(x₂y₁)
    [[1,4]],      # ς(x₁y₂)
    [[2,4]],      # ς(x₂y₂)
    [[1], [3]],   # ς(x₁)ς(y₁)
    [[2], [3]],   # ς(x₂)ς(y₁)
    [[2], [4]],   # ς(x₂)ς(y₂)
    [[1], [1]],   # ς(x₁)²
    [[4], [4]],   # ς(y₂)²
]
coe_7_2_3 = -[1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

println("# Variables: x[1:2]={1,2}, y[1:2]={3,4}")
println("# Objective: Mixed linear, product, and squared expectations")
println("# Expected: -3.5114802")
println()

results_7_2_3 = map(STATE_7_2_3_VARIANTS) do v
    key = "State_7_2_3_$(v.name)_d$(v.order)"
    println("# $(v.name) (order=$(v.order), TS=$(v.ts))")
    
    opt, data = pstateopt_first(supp_7_2_3, coe_7_2_3, 4, v.order;
        vargroup=[2, 2],
        TS=v.ts == false ? false : v.ts,
        constraint="unipotent")
    result = extract_oracle(key, opt, data; use_cs=false)
    print_oracle(result)
    println()
    result
end

# Summary
all_results = vcat(results_7_2_1, results_7_2_2, results_7_2_3)
print_summary("STATE_POLY_EXTENDED", all_results)
