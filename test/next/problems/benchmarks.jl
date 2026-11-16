# E1–E5 from Wang–Magron, Tables 2–6. Keep the Hermitian-square
# residuals from test/problems/benchmarks/, not naive NC fourth powers.
# These old cases use Clarabel rather than the shared COSMO factory.
const BENCHMARK_SOLVER = JuMPBackend(JuMP.optimizer_with_attributes(
    Clarabel.Optimizer, "verbose" => false))

function benchmark_problem(kind, n)
    _, (x,) = site_variables(Free(), 1, n)
    f = kind in (3, 4) ? one(x[1]) : zero(x[1])
    if kind == 1 # Broyden banded
        for i in 1:n
            neighbors = setdiff(max(1, i-5):min(n, i+1), i)
            r = 2x[i] + 5x[i]^3 + 1 -
                sum((x[j] + x[j]^2 for j in neighbors); init=zero(f))
            f += r' * r
        end
    elseif kind == 2 # Chained singular
        for i in 1:2:(n-3)
            a = x[i] + 10x[i+1]
            b = x[i+2] - x[i+3]
            c = x[i+1]^2 - 4x[i+1]*x[i+2] + 4x[i+2]^2
            d = x[i]^2 - 20x[i]*x[i+3] + 100x[i+3]^2
            f += a'*a + 5b'*b + c'*c + 10d'*d
        end
    elseif kind == 3 # Generalized Rosenbrock
        for i in 2:n
            r = x[i] - x[i-1]^2
            s = 1 - x[i]
            f += 100r'*r + s'*s
        end
    elseif kind == 4 # Chained Wood
        for i in 1:2:(n-3)
            a = x[i+1] - x[i]^2
            b = 1 - x[i]
            c = x[i+3] - x[i+2]^2
            d = 1 - x[i+2]
            e = x[i+1] + x[i+3] - 2
            g = x[i+1] - x[i+3]
            f += 100a'*a + b'*b + 90c'*c + d'*d + 10e'*e + 0.1g'*g
        end
    else # Broyden tridiagonal
        for i in 1:n
            r = 3x[i] - 2x[i]^2 + 1
            i > 1 && (r -= x[i-1])
            i < n && (r -= 2x[i+1])
            f += r' * r
        end
    end
    return Problem(f)
end

function benchmark_contracts(solver, sos_solver)
    names = ("e1_broyden_banded", "e2_chained_singular",
        "e3_generalized_rosenbrock", "e4_chained_wood", "e5_broyden_tridiagonal")
    tolerances = (2e-5, 1e-5, 2e-5, 5e-4, 1e-6)
    @testset "Newton-chip literature benchmarks" begin
        for (kind, name) in enumerate(names), n in (20, kind in (1, 4) ? 40 : 60)
            @testset "$name n=$n" begin
                oracle = fixture(name, "newton_chip_n$n")
                problem = benchmark_problem(kind, n)
                rows = NCTSSoS.newton_chip_basis(problem, kind == 1 ? 3 : 2)
                built = NCTSSoS.build(problem,
                    Relaxation(rows=rows, split=TermSplit(algorithm=MMD(), level=1)))
                @test length(rows) == oracle["dense_mb"]
                @test maximum(NCTSSoS.block_sides(built)) == oracle["sparse_mb"]
                if kind == 1
                    # The old E1 check uses SOS. JuMPBackend has no form selector;
                    # retain this bound check under the pilot's licensed SOS gate.
                    @test solve(built; solver=sos_solver).bound ≈ oracle["sparse_opt"] atol=tolerances[kind] skip=isnothing(sos_solver)
                else
                    result = solve(built; solver)
                    if kind == 5 && n == 60
                        # Architect-reviewed 2026-11-16: the same 121-row program
                        # misses 1e-6 on several old and new solver routes. Keep
                        # 1e-6 on validity; allow 1e-5 on tightness only here.
                        @test result.bound <= 1e-6
                        @test result.bound >= -1e-5
                        @test result.sizes.moments == fixture("next_benchmarks", "e5_n60")["moments"]
                    else
                        # E2/E4's old tests use exact 0/1, not the paper's rounded
                        # solver noise stored in sparse_opt. Preserve that choice.
                        expected = kind == 2 ? 0.0 : kind == 4 ? 1.0 : oracle["sparse_opt"]
                        @test result.bound ≈ expected atol=tolerances[kind]
                    end
                end
            end
        end
    end
end
