# SOS (Sum of Squares) Dualization Tests
# Tests specific to SOS dualization components:
#   - Cαj coefficient extraction
#
# Note: Problem-specific SOS tests are in problems/ subdirectory.

using Test, NCTSSoS, JuMP
using SparseArrays
using NCTSSoS: get_Cαj

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

@testset "SOS Components" begin
    @testset "Cαj" begin
        # Test get_Cαj with new API: Vector{Monomial}, Matrix{Polynomial}
        registry, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Create a properly typed polynomial matrix
        poly_matrix = [
            1.0*x[1] 1.0*x[2]
            1.0*x[2] 1.0*x[1]+1.0*x[3]
        ]

        # Basis: the monomials themselves (sorted)
        basis_monomials = sort([one(x[1]), x[1], x[2], x[3]])

        C_α_js = get_Cαj(basis_monomials, poly_matrix)

        # Just verify structure - we have 5 non-zero entries
        @test length(C_α_js) == 5
        @test all(v -> v == 1.0, values(C_α_js))
    end

    @testset "Cαj complex" begin
        registry, (x,) = create_noncommutative_variables([("x", 1:2)])
        basis_polys = get_ncbasis(registry, 2)
        # For NC algebra, each basis poly is single-term; extract the monomial
        basis_monomials = [monomials(p)[1] for p in basis_polys]

        # Create properly typed polynomial matrix
        localizing_mtx = [
            1.0*x[1]-1.0*one(x[1])   1.0*x[2]^2-1.0*one(x[1])
            1.0*x[2]^2-1.0*one(x[1])   1.0*x[2]^3
        ]
        C_α_js = get_Cαj(basis_monomials, localizing_mtx)

        @test !isempty(C_α_js)
        @test all(v -> v == 1.0 || v == -1.0, values(C_α_js))
    end
end
