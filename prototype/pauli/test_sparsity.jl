using Test: @test, @testset, @test_throws
import NCTSSoS as NC

# Include beside the prototype's test runner; pass its module and solver.
function check_workflow_sparsity(workflow; solver)
    support = workflow.workflow_support
    split = workflow.workflow_term_sparsity
    function bound(p)
        result = NC.solve_standard(NC.moment_standard_form(p), solver)
        isnothing(result.bound) && error("Sparsity check requires an optimum; got $(result.status): $(result.raw_status)")
        return result.bound
    end
    function cone!(p, kind, side, forms; label=0)
        for form in forms
            for (id, c) in form
                NC.push_term!(p.forms, c, id)
            end
            NC.close_form!(p.forms)
        end
        NC.push_block!(p, kind, side, label)
    end

    @testset "Pauli workflow sparsity" begin
        @testset "3-spin path: support bound and graph powers" begin
            _, (x, y, z) = NC.pauli_variables(3)
            H = sum(a[i]*a[i+1] for a in (x, y, z), i in 1:2)
            reflection = NC.SitePermutation([3, 2, 1])
            problem = NC.Problem(H; symmetry=(reflection,))
            full = support(H, 2)
            edges = support(H, 2; support_mode=:correlative_cliques, power_graph_radius=1,
                symmetry_generators=(reflection,))
            powered = support(H, 2; support_mode=:correlative_cliques, power_graph_radius=2)
            @test full.row_counts == powered.row_counts == [37]
            @test edges.regions == [[1, 2], [2, 3]]
            @test edges.row_counts == [16, 16]
            @test (full.relation, edges.relation) == (:exact, :relaxation)
            @test_throws ArgumentError support(H, 2; max_rows=36)
            dense = NC.build(problem, NC.Relaxation(family=full.family))
            sparse = NC.build(problem, NC.Relaxation(family=edges.family))
            @test length(sparse.sizes.closure.block_orbits) == 1
            @test bound(dense.program) ≈ -4 atol=2e-6
            @test bound(sparse.program) ≈ -6 atol=2e-6
            exact = split(dense.program)
            fixed = split(dense.program; ts_algo=:MMD, ts_step=1)
            @test exact.report.zero_lift_audit
            @test bound(exact.program) ≈ -4 atol=2e-6
            @test bound(fixed.program) <= -4 + 2e-6

            strengthened = NC.build(problem, NC.Relaxation(family=full.family,
                strengthen=(NC.StateOptimality(), NC.DensityMatrix(1:3), NC.Stationarity())))
            guarded = split(strengthened.program)
            retained = filter(b -> b.label != 0, strengthened.program.blocks)
            @test count(b -> b.kind in (NC.PSD, NC.HPSD), retained) >= 2
            @test any(b -> b.kind == NC.ZERO, retained)
            @test all(id -> id in guarded.report.active_ids,
                (strengthened.program.forms.ids[t] for b in retained for f in NC.form_range(b)
                    for t in NC.term_range(strengthened.program.forms, f)))
        end

        @testset "Periodic lattice stencils are spin-closed" begin
            _, (x, _, _) = NC.pauli_variables(7)
            chain = support(x[1], 3; support_mode=:contiguous_lattice, two_body_range_r=2)
            @test chain.row_counts == [1 + 7*3 + 14*9 + 7*27]
            @test [1, 6, 7] in chain.regions
            @test [1, 3, 7] ∉ chain.regions
            @test [1, 6] in chain.regions
            @test [1, 4] ∉ chain.regions
            _, (v, _, _) = NC.pauli_variables(12)
            grid = support(v[1], 4; support_mode=:contiguous_lattice, lattice_shape=(4, 3), two_body_range_r=1)
            @test grid.row_counts == [1 + 12*3 + 24*9 + 48*27 + 12*81]
            @test [1, 4, 9, 12] in grid.regions
            @test [1, 2, 3, 4] ∉ grid.regions
            @test [1, 6] ∉ grid.regions
            extended = support(v[1], 2; support_mode=:contiguous_lattice, lattice_shape=(4, 3), two_body_range_r=2)
            @test [1, 6] in extended.regions

            _, (u, _, w) = NC.pauli_variables(4)
            H = sum(u[i]*u[mod1(i+1, 4)] + 0.17*w[i] for i in 1:4)
            translation = NC.SitePermutation([2, 3, 4, 1])
            problem = NC.Problem(H; symmetry=(translation,))
            lattice = support(H, 2; support_mode=:contiguous_lattice, two_body_range_r=1,
                symmetry_generators=(translation,))
            dense = NC.build(problem, NC.Relaxation(family=support(H, 2).family))
            sparse = NC.build(problem, NC.Relaxation(family=lattice.family))
            @test lattice.row_counts == [49]
            @test sparse.sizes.closure.route == :lattice
            @test bound(sparse.program) <= bound(dense.program) + 2e-6
        end

        @testset "Cross-cone closure and cancellation" begin
            c(id) = [(id, 1)]
            p = NC.MomentProgram{Float64}(9)
            push!(p.objective_ids, 2); push!(p.objective_coefs, -1)
            # Seed 2 plus the equality seed 5 induce 3 inside the first block;
            # then 3 plus inequality seed 6 induce 4 in the second block.
            cone!(p, NC.PSD, 3, [c(1), c(2), c(3), c(1), c(5), c(1)])
            cone!(p, NC.PSD, 3, [c(1), c(3), c(4), c(1), c(6), c(1)])
            cancel = [(2, 1e16), (2, 1.0), (2, -1e16), (2, -1.0), (9, 1.0)]
            cone!(p, NC.PSD, 2, [c(1), cancel, c(1)])
            # PSD stand-ins for registered KMS and density cones are immutable.
            cone!(p, NC.PSD, 2, [c(1), c(7), c(1)]; label=3)
            cone!(p, NC.PSD, 2, [c(1), c(8), c(1)]; label=4)
            cone!(p, NC.ZERO, 1, [c(5)]; label=1)
            cone!(p, NC.NONNEGATIVE, 1, [c(6)]; label=2)
            exact = split(p)
            @test exact.report.steps == 3
            @test exact.report.active_ids == collect(1:8)
            @test exact.report.cliques[3] == [[1], [2]]
            @test [(b.kind, b.side, b.label) for b in exact.program.blocks[end-3:end]] ==
                [(b.kind, b.side, b.label) for b in p.blocks[end-3:end]]
            @test split(p; protected_ids=(9,)).report.cliques[3] == [[1, 2]]
            @test_throws ArgumentError split(p; split_blocks=[4])
        end

        @testset "Exact sparse PSD versus principal-clique tests" begin
            unit, zero_form = [(1, 1)], Tuple{Int,Int}[]
            path = NC.MomentProgram{ComplexF64}(2)
            push!(path.objective_ids, 2); push!(path.objective_coefs, -1)
            # Eigenvalues are 1 and 1 ± sqrt(5)t. Each edge test alone gives t ≤ 1/2.
            cone!(path, NC.HPSD, 3, [unit, [(2, im)], zero_form, unit, [(2, 2)], unit])
            exact = split(path; ts_algo=:AggregateDualChordal)
            relaxed = split(path; ts_algo=:MMD, ts_step=1)
            @test exact.report.relation == :exact
            @test relaxed.report.relation == :relaxation
            @test length(exact.report.auxiliary_ids) == 8
            @test bound(path) ≈ -inv(sqrt(5)) atol=2e-7
            @test bound(exact.program) ≈ -inv(sqrt(5)) atol=2e-7
            @test bound(relaxed.program) ≈ -0.5 atol=2e-7

            cycle = NC.MomentProgram{Float64}(2)
            push!(cycle.objective_ids, 2); push!(cycle.objective_coefs, -1)
            edge = [(2, 1)]
            cone!(cycle, NC.PSD, 4, [unit, edge, zero_form, edge, unit, edge, zero_form, unit, edge, unit])
            exact = split(cycle; ts_algo=:AggregateDualChordal)
            merged = split(cycle; ts_algo=:AggregateDualChordal, merge_threshold=4)
            relaxed = split(cycle; ts_algo=:MMD, ts_step=1)
            @test sort(length.(exact.report.cliques[1])) == [3, 3]
            @test length(merged.report.cliques[1]) == 1
            # The even cycle has eigenvalue 1-2t, while each completed principal
            # clique is a 3-vertex path. A free fill entry would change the exact bound.
            @test bound(exact.program) ≈ -0.5 atol=2e-7
            @test bound(merged.program) ≈ -0.5 atol=2e-7
            @test bound(relaxed.program) ≈ -inv(sqrt(2)) atol=2e-7
        end
    end
end
