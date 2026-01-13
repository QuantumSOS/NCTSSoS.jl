# NCTSSOS Oracle Script: NC Polynomial Benchmarks
# Run on server with NCTSSOS + MosekTools:
#   cd ~/NCTSSOS && julia --project path/to/nctssos_benchmarks.jl
#
# Benchmarks included:
#   1. Broyden Banded (n=4, degree=6) - Global min = 0
#   2. Broyden Tridiagonal (n=6, degree=4) - Global min = 0
#   3. Chained Singular (n=8, degree=4) - Global min = 0
#   4. Chained Wood (n=8, degree=4) - Global min ≈ 1
#
# Note: Rosenbrock already has oracles in nctssos_rosenbrock.jl

include("oracle_utils.jl")

# Broyden Banded (n=4, d=3)
# f = n + Σᵢ₌₁ⁿ [ 4xᵢ + 4xᵢ² + 10xᵢ³ + 20xᵢ⁴ + 25xᵢ⁶
#                 + Σⱼ∈Jᵢ (-2xⱼ - 2xⱼ² - 4xᵢxⱼ - 4xᵢxⱼ² - 10xᵢ³xⱼ - 10xᵢ³xⱼ²)
#                 + Σⱼ,ₖ∈Jᵢ (xⱼxₖ + 2xⱼ²xₖ + xⱼ²xₖ²) ]
# where Jᵢ = {max(1,i-5), ..., min(n,i+1)} \ {i}
# Global minimum: 0 at origin
# Degree: 6, so d=3 is required

const BROYDEN_BANDED_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=3, n=4),
    (name="CS", cs="MF", ts=false, order=3, n=4),
    (name="CS_TS", cs="MF", ts="MD", order=3, n=4),
]

print_header("Broyden Banded (n=4)")

println("# Degree: 6, Order: 3")
println("# Expected: 0.0")
println()

results_broyden_banded = map(BROYDEN_BANDED_VARIANTS) do v
    n = v.n
    @ncpolyvar x[1:n]

    # Build Broyden Banded objective
    obj = n * 1.0
    for i in 1:n
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)

        # Single variable terms
        obj += 4.0 * x[i] + 4.0 * x[i]^2 + 10.0 * x[i]^3 + 20.0 * x[i]^4 + 25.0 * x[i]^6

        # Cross terms with j in Jᵢ
        for j in jset
            obj += -2.0 * x[j] - 2.0 * x[j]^2 - 4.0 * x[i] * x[j] -
                   4.0 * x[i] * x[j]^2 - 10.0 * x[i]^3 * x[j] - 10.0 * x[i]^3 * x[j]^2
        end

        # Double sum terms
        for j in jset
            for k in jset
                obj += x[j] * x[k] + 2.0 * x[j]^2 * x[k] + x[j]^2 * x[k]^2
            end
        end
    end

    key = "BroydenBanded_$(v.name)_n$(v.n)_d$(v.order)"
    println("# $(v.name) (n=$(v.n), order=$(v.order), CS=$(v.cs), TS=$(v.ts))")

    if v.cs == false
        opt, data = nctssos_first([obj], x, v.order; TS=v.ts)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first([obj], x, v.order; TS=v.ts, CS=v.cs)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

# Broyden Tridiagonal (n=6, d=2)
# f = n + Σᵢ₌₁ⁿ gᵢ(x)
# Global minimum: 0 at origin
# Degree: 4, so d=2 is sufficient

const BROYDEN_TRIDIAGONAL_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2, n=6),
    (name="CS", cs="MF", ts=false, order=2, n=6),
    (name="CS_TS", cs="MF", ts="MD", order=2, n=6),
]

print_header("Broyden Tridiagonal (n=6)")

println("# Degree: 4, Order: 2")
println("# Expected: 0.0")
println()

results_broyden_tridiagonal = map(BROYDEN_TRIDIAGONAL_VARIANTS) do v
    n = v.n
    @ncpolyvar x[1:n]

    # Build Broyden Tridiagonal objective
    obj = n * 1.0

    # First variable (i=1)
    obj += 5.0 * x[1]^2 + 4.0 * x[1]^4 + 4.0 * x[2]^2 - 12.0 * x[1]^3 -
           12.0 * x[1] * x[2] + 6.0 * x[1] + 8.0 * x[1]^2 * x[2] - 4.0 * x[2]

    # Middle variables (i=2,...,n-1)
    for i in 2:(n-1)
        obj += 5.0 * x[i]^2 + 4.0 * x[i]^4 + x[i-1]^2 + 4.0 * x[i+1]^2 -
               12.0 * x[i]^3 - 6.0 * x[i-1] * x[i] - 12.0 * x[i] * x[i+1] +
               6.0 * x[i] + 4.0 * x[i-1] * x[i]^2 + 8.0 * x[i]^2 * x[i+1] +
               4.0 * x[i-1] * x[i+1] - 2.0 * x[i-1] - 4.0 * x[i+1]
    end

    # Last variable (i=n)
    obj += 5.0 * x[n]^2 + 4.0 * x[n]^4 + x[n-1]^2 - 12.0 * x[n]^3 -
           6.0 * x[n-1] * x[n] + 6.0 * x[n] + 4.0 * x[n-1] * x[n]^2 - 2.0 * x[n-1]

    key = "BroydenTridiagonal_$(v.name)_n$(v.n)_d$(v.order)"
    println("# $(v.name) (n=$(v.n), order=$(v.order), CS=$(v.cs), TS=$(v.ts))")

    if v.cs == false
        opt, data = nctssos_first([obj], x, v.order; TS=v.ts)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first([obj], x, v.order; TS=v.ts, CS=v.cs)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

# Chained Singular (n=8, d=2)
# f = Σᵢ₌₁,₃,₅,...,ₙ₋₃ [ (xᵢ + 10xᵢ₊₁)² + 5(xᵢ₊₂ - xᵢ₊₃)²
#                        + (xᵢ₊₁ - 2xᵢ₊₂)⁴ + 10(xᵢ - 10xᵢ₊₃)⁴ ]
# Global minimum: 0 at origin
# Degree: 4, so d=2 is sufficient

const CHAINED_SINGULAR_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2, n=8),
    (name="CS", cs="MF", ts=false, order=2, n=8),
    (name="CS_TS", cs="MF", ts="MD", order=2, n=8),
]

print_header("Chained Singular (n=8)")

println("# Degree: 4, Order: 2")
println("# Expected: 0.0")
println()

results_chained_singular = map(CHAINED_SINGULAR_VARIANTS) do v
    n = v.n
    @ncpolyvar x[1:n]

    # NC expansion helpers
    # (a + b)² = a² + ab + ba + b²
    nc_square(a, b) = a^2 + a * b + b * a + b^2

    # (a + b)⁴ = ((a+b)²)²
    function nc_fourth_power(a, b)
        ab = a + b
        ab2 = ab * ab
        return ab2 * ab2
    end

    obj = 0.0 * x[1]  # Start with polynomial type

    for i in 1:2:(n-3)
        # (xᵢ + 10xᵢ₊₁)²
        obj += nc_square(x[i], 10.0 * x[i+1])

        # 5(xᵢ₊₂ - xᵢ₊₃)²
        obj += 5.0 * nc_square(x[i+2], -1.0 * x[i+3])

        # (xᵢ₊₁ - 2xᵢ₊₂)⁴
        obj += nc_fourth_power(x[i+1], -2.0 * x[i+2])

        # 10(xᵢ - 10xᵢ₊₃)⁴
        obj += 10.0 * nc_fourth_power(x[i], -10.0 * x[i+3])
    end

    key = "ChainedSingular_$(v.name)_n$(v.n)_d$(v.order)"
    println("# $(v.name) (n=$(v.n), order=$(v.order), CS=$(v.cs), TS=$(v.ts))")

    if v.cs == false
        opt, data = nctssos_first([obj], x, v.order; TS=v.ts)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first([obj], x, v.order; TS=v.ts, CS=v.cs)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

# Chained Wood (n=8, d=2)
# f = (21n - 41) + Σᵢ₌₁,₃,₅,...,ₙ₋₃ gᵢ(x)
# Global minimum: ≈1
# Degree: 4, so d=2 is sufficient

const CHAINED_WOOD_VARIANTS = [
    (name="Dense", cs=false, ts=false, order=2, n=8),
    (name="CS", cs="MF", ts=false, order=2, n=8),
    (name="CS_TS", cs="MF", ts="MD", order=2, n=8),
]

print_header("Chained Wood (n=8)")

println("# Degree: 4, Order: 2")
println("# Expected: ≈1.0")
println()

results_chained_wood = map(CHAINED_WOOD_VARIANTS) do v
    n = v.n
    @ncpolyvar x[1:n]

    # Constant term
    obj = Float64(21 * n - 41)

    for i in 1:2:(n-3)
        obj += (-2.0) * x[i] + x[i]^2 + 100.0 * x[i]^4 - 200.0 * x[i]^2 * x[i+1] +
               (-40.0) * x[i+1] + 110.1 * x[i+1]^2 + 19.8 * x[i+1] * x[i+3] +
               (-2.0) * x[i+2] + x[i+2]^2 + 90.0 * x[i+2]^4 - 180.0 * x[i+2]^2 * x[i+3] +
               (-40.0) * x[i+3] + 100.1 * x[i+3]^2
    end

    key = "ChainedWood_$(v.name)_n$(v.n)_d$(v.order)"
    println("# $(v.name) (n=$(v.n), order=$(v.order), CS=$(v.cs), TS=$(v.ts))")

    if v.cs == false
        opt, data = nctssos_first([obj], x, v.order; TS=v.ts)
        result = extract_oracle(key, opt, data; use_cs=false)
    else
        opt, data = cs_nctssos_first([obj], x, v.order; TS=v.ts, CS=v.cs)
        result = extract_oracle(key, opt, data; use_cs=true)
    end
    print_oracle(result)
    println()
    result
end

# Summary
all_results = vcat(
    results_broyden_banded,
    results_broyden_tridiagonal,
    results_chained_singular,
    results_chained_wood
)
print_summary("BENCHMARKS", all_results)
