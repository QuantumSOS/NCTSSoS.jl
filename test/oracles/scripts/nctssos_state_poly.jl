# NCTSSOS Oracle Script: State Polynomial Optimization
# =====================================================
# Run on a800 with MosekTools:
#   cd /home/yushengzhao/NCTSSOS && julia --project path/to/nctssos_state_poly.jl
#
# Output: Oracle values for state polynomial problems (Section 7.2)
#
# API: stateopt_first for state polynomial optimization
#      Uses coefficient/support format for state expectations

using NCTSSOS, DynamicPolynomials
using MosekTools

println("=" ^ 70)
println("NCTSSOS Oracle: State Polynomial Optimization")
println("=" ^ 70)
println()

# Store results
results = Dict{String, NamedTuple}()

# Helper to extract oracle info
function extract_stateopt(name, opt, data)
    sides = [size(M, 1) for M in data.moment]
    nuniq = length(data.ksupp)
    results[name] = (opt=opt, sides=sides, nuniq=nuniq)
    println("\"$name\" => (opt=$opt, sides=$sides, nuniq=$nuniq),")
end

function extract_cs_stateopt(name, opt, data)
    sides = [size(M, 1) for clique in data.moment for M in clique]
    nuniq = length(data.ksupp)
    results[name] = (opt=opt, sides=sides, nuniq=nuniq)
    println("\"$name\" => (opt=$opt, sides=$sides, nuniq=$nuniq),")
end

# =============================================================================
# 7.2.0: CHSH State Polynomial
# =============================================================================
println("\n# 7.2.0: CHSH State Polynomial")
println("# Expected: -2√2 ≈ -2.8284")
println()

@ncpolyvar x[1:2] y[1:2]

# Dense (order=1)
begin
    coe = [-1.0, -1.0, -1.0, 1.0]
    supp = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]
    
    println("# Dense (order=1)")
    opt, data = stateopt_first(coe, supp, [x;y], 1; 
                               TS=false, partition=2, constraint="unipotent")
    extract_stateopt("State_7_2_0_Dense_d1", opt, data)
end

# Term sparsity (order=1)
begin
    coe = [-1.0, -1.0, -1.0, 1.0]
    supp = [[[1,3]], [[1,4]], [[2,3]], [[2,4]]]
    
    println("# TS (order=1)")
    opt, data = stateopt_first(coe, supp, [x;y], 1; 
                               TS="MD", partition=2, constraint="unipotent")
    extract_stateopt("State_7_2_0_TS_d1", opt, data)
end

# =============================================================================
# 7.2.1: Squared Expectations
# =============================================================================
println("\n# 7.2.1: Squared Expectations")
println("# Expected: -4.0 (at order=3)")
println()

@ncpolyvar x[1:2] y[1:2]

# Dense (order=3)
begin
    # sp = -(ς(x₁y₂) + ς(x₂y₁))² - (ς(x₁y₁) - ς(x₂y₂))²
    # Expanded: -<14><14> - 2<14><23> - <23><23> - <13><13> + 2<13><24> - <24><24>
    coe = [-1.0, -2.0, -1.0, -1.0, 2.0, -1.0]
    supp = [[[1,4], [1,4]], [[1,4], [2,3]], [[2,3], [2,3]], 
            [[1,3], [1,3]], [[1,3], [2,4]], [[2,4], [2,4]]]
    
    println("# Dense (order=3)")
    opt, data = stateopt_first(coe, supp, [x;y], 3; 
                               TS=false, partition=2, constraint="unipotent")
    extract_stateopt("State_7_2_1_Dense_d3", opt, data)
end

# =============================================================================
# 7.2.2: Covariance Expression
# =============================================================================
println("\n# 7.2.2: Covariance Expression")
println("# Expected: -5.0")
println()

@ncpolyvar x[1:3] y[1:3]

# Dense (order=2)
begin
    # cov(i,j) = <x[i]*y[j]> - <x[i]><y[j]>
    # sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
    coe = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0,  # cov(1,1), cov(1,2), cov(1,3)
           1.0, -1.0, 1.0, -1.0, -1.0, 1.0,  # cov(2,1), cov(2,2), -cov(2,3)
           1.0, -1.0, -1.0, 1.0]             # cov(3,1), -cov(3,2)
    supp = [[[1,4]], [[1], [4]], [[1,5]], [[1], [5]], [[1,6]], [[1], [6]],
            [[2,4]], [[2], [4]], [[2,5]], [[2], [5]], [[2,6]], [[2], [6]],
            [[3,4]], [[3], [4]], [[3,5]], [[3], [5]]]
    
    println("# Dense (order=2)")
    opt, data = stateopt_first(coe, supp, [x;y], 2; 
                               TS=false, partition=3, constraint="unipotent")
    extract_stateopt("State_7_2_2_Dense_d2", opt, data)
end

# =============================================================================
# 7.2.3: Mixed State Polynomial
# =============================================================================
println("\n# 7.2.3: Mixed State Polynomial")
println("# Expected: -3.5114802")
println()

@ncpolyvar x[1:2] y[1:2]

# Dense (order=2)
begin
    # From NCTSSOS: coe = -[1; 1; 1; -1; 1; 1; 1; -1; -1; -1; -1; -1]
    coe = [-1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    supp = [[[2]], [[3]], [[4]], [[1,3]], [[2,3]], [[1,4]], [[2,4]], 
            [[1], [3]], [[2], [3]], [[2], [4]], [[1], [1]], [[4], [4]]]
    
    println("# Dense (order=2)")
    opt, data = stateopt_first(coe, supp, [x;y], 2; 
                               TS=false, partition=2, constraint="unipotent")
    extract_stateopt("State_7_2_3_Dense_d2", opt, data)
end

# Term sparsity (order=2)
begin
    coe = [-1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    supp = [[[2]], [[3]], [[4]], [[1,3]], [[2,3]], [[1,4]], [[2,4]], 
            [[1], [3]], [[2], [3]], [[2], [4]], [[1], [1]], [[4], [4]]]
    
    println("# TS (order=2)")
    opt, data = stateopt_first(coe, supp, [x;y], 2; 
                               TS="MD", partition=2, constraint="unipotent")
    extract_stateopt("State_7_2_3_TS_d2", opt, data)
end

# =============================================================================
# Summary
# =============================================================================
println()
println("=" ^ 70)
println("STATE_POLY ORACLE SUMMARY")
println("=" ^ 70)
println()
println("const STATE_POLY_ORACLES = Dict(")
for (name, res) in sort(collect(results), by=first)
    println("    \"$name\" => (opt=$(res.opt), sides=$(res.sides), nuniq=$(res.nuniq)),")
end
println(")")
