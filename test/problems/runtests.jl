# Problem-based tests runner.

using Test

@testset "Problems" begin
    # Bell inequalities (CHSH, I_3322)
    @testset "Bell Inequalities" begin
        include("bell_inequalities/runtests.jl")
    end

    # NC polynomial examples
    @testset "NC Polynomial" begin
        include("nc_polynomial/runtests.jl")
    end

    # State polynomial examples
    @testset "State Polynomial" begin
        include("state_polynomial/runtests.jl")
    end

    # Trace polynomial examples
    @testset "Trace Polynomial" begin
        include("trace_polynomial/runtests.jl")
    end

    # Benchmarks (Rosenbrock, etc.)
    @testset "Benchmarks" begin
        include("benchmarks/runtests.jl")
    end

    # Tests requiring --local (Mosek)
    if USE_LOCAL
        # Condensed matter physics
        @testset "Condensed Matter" begin
            include("condensed_matter/runtests.jl")
        end

        # Quantum networks
        @testset "Quantum Networks" begin
            include("quantum_networks/runtests.jl")
        end

        # Fermionic systems
        @testset "Fermionic" begin
            include("fermionic/runtests.jl")
        end
    end
end
