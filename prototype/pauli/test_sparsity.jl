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

        @testset "Hamiltonian-derived pairs, connected products and coverage" begin
            _, (x, y, z) = NC.pauli_variables(5)
            H = x[1]*x[3] + x[3]*x[4] + y[3]*y[5]
            words(p) = Set(last(only(NC.terms(q))) for q in only(p.family))
            word(p) = last(only(NC.terms(p)))
            pair_sites(p) = Set(Tuple(r) for r in p.regions if length(r) == 2)
            near = support(H, 2; support_mode=:graph_distance, power_graph_radius=1, max_rows=43)
            wide = support(H, 2; support_mode=:graph_distance)
            local2 = support(H, 2; support_mode=:local_products)
            @test pair_sites(near) == Set([(1,3), (3,4), (3,5)])
            @test pair_sites(wide) == Set([(1,3), (3,4), (3,5), (1,4), (1,5), (4,5)])
            @test near.generation == (considered=46, retained=43, dense_rows=106, products=0)
            @test wide.row_counts == [70]
            @test Set(w for w in words(local2) if NC.degree(w) == 2) ==
                Set(word.([x[1]*x[3], x[3]*x[4], y[3]*y[5], x[1]*x[4]]))
            @test local2.generation == (considered=26, retained=20, dense_rows=106, products=9)
            @test word(x[1]*z[3]*y[5]) in words(support(H, 3; support_mode=:local_products))
            @test word(x[1]*y[3]*x[4]*y[5]) ∈ words(support(H, 4; support_mode=:local_products, local_product_depth=3))
            @test word(x[1]*y[3]*x[4]*y[5]) ∉ words(support(H, 4; support_mode=:local_products))
            @test_throws ArgumentError support(H, 2; support_mode=:graph_distance, power_graph_radius=1, max_rows=42)

            # A chordless cycle with a tail and an isolated site, with
            # unequal and mixed-axis interactions rather than XXX bonds.
            _, (a,b,c) = NC.pauli_variables(6)
            irregular = a[4]*b[2] - 2b[2]*c[5] + 0.5a[5]*c[1] +
                3a[1]*b[4] - 0.25c[5]*b[3] + 0.125a[4]
            motifs(h, order; kwargs...) = support(h, order;
                support_mode=:graph_motifs,power_graph_radius=1,kwargs...)
            m4 = motifs(irregular, 4)
            pairs = [[1,4],[1,5],[2,4],[2,5],[3,5]]
            triples = [[1,2,4],[1,2,5],[1,3,5],[1,4,5],[2,3,5],[2,4,5]]
            @test m4.regions == sort!([[i] for i in 1:6] ∪ pairs ∪ triples ∪ [[1,2,4,5]])
            @test m4.row_counts == [1 + 6*3 + 5*9 + 6*27 + 81]
            @test motifs(irregular, 3).row_counts == [only(m4.row_counts) - 81]
            @test [1,2,4,5] ∉ motifs(irregular + a[1]*b[2], 4).regions
            @test_throws ArgumentError motifs(irregular, 4;max_rows=only(m4.row_counts)-1)
            @test_throws ArgumentError motifs(irregular, 5)

            obs2 = b[2]*c[6]
            m2 = motifs(irregular, 2;observable=obs2)
            @test word(obs2) in words(m2) && [2,6] ∉ m2.regions
            plan = workflow.compile_relaxation(irregular;order=2,observable=obs2,
                sparsity=workflow.PauliSparsitySpec(support_mode=:graph_motifs,power_graph_radius=1,
                    ts_algo=:None,rewrite_complete_local_rdm=false))
            @test plan.reports.support.mode === :graph_motifs
            @test plan.reports.support.row_counts == m2.row_counts
            @test :graph_motifs in plan.provenance.active_relaxations
            @test support(irregular,2).row_counts == [154] # The default is still full.

            # Off-stencil two-body and higher-than-order moments. A relabeling
            # changes the natural sorted-site cut, but not the family of cuts.
            obs = z[2]*y[5] + x[1]*z[2]*y[4]*z[5]
            perm = [2, 4, 1, 5, 3]
            action = NC.SitePermutation(perm)
            Hp = x[perm[1]]*x[perm[3]] + x[perm[3]]*x[perm[4]] + y[perm[3]]*y[perm[5]]
            obsp = z[perm[2]]*y[perm[5]] + x[perm[1]]*z[perm[2]]*y[perm[4]]*z[perm[5]]
            for mode in (:graph_distance, :local_products, :graph_motifs)
                s = support(H, 2; support_mode=mode, observable=obs)
                t = support(Hp, 2; support_mode=mode, observable=obsp)
                selected = words(s)
                @test Set(last(NC.act(action, w)) for w in selected) == words(t)
                @test s.relation == :relaxation && pair_sites(s) == pair_sites(support(H, 2; support_mode=mode))
                @test all(w -> any(v -> last(NC.mul(v,w)) in selected, selected), last.(collect(NC.terms(obs))))
                @test_throws ArgumentError support(H, 2; support_mode=mode, observable=prod(x))
            end

            s = support(H, 2; support_mode=:local_products, observable=z[2]*y[5])
            rows = [sort!(collect(words(s)))]
            W = eltype(only(rows))
            gram = [NC.RowBlock(only(rows), [(1.0+0im,one(W))], NC.HPSD, 0)]
            objective = [(ComplexF64(c),w) for (c,w) in NC.terms(H)]
            p, table, _ = workflow.workflow_symmetry_program(gram, objective, NC.FamilySymmetry();
                spatial_mode=:none, realify=false)
            @test length(p.forms) == 21*22÷2 # One full triangle on the selected rows.
            rdm_word = word(y[1]*y[2]*y[5])
            @test isnothing(NC.moment_terms(table, rdm_word))
            extra = workflow.workflow_strengthening(H, 2; moment_rows=rows,
                commutator_linear_order=1, rdm_positivity_k=3, rdm_regions=[[1,2,5]])
            p = workflow._append_strengthening!(p, table, extra)
            observed = first.(NC.moment_terms(table, word(z[2]*y[5])))
            rdmid = only(NC.moment_terms(table, rdm_word))[1]
            guarded = split(p; protected_ids=observed)
            @test guarded.report.zero_lift_audit && all(in(guarded.report.active_ids), (rdmid, observed...))
            @test any(b -> b.label > 0 && b.kind == NC.ZERO, guarded.program.blocks)
            @test any(b -> b.label > 0 && b.kind == NC.HPSD && b.side == 8, guarded.program.blocks)
        end

        @testset "Selected rows give the exact small-spin ground energy" begin
            _, (x, y, z) = NC.pauli_variables(3)
            H = sum(a[i]*a[i+1] for a in (x,y,z), i in 1:2)
            # Two unit Heisenberg bonds: S_13=1 and total spin 1/2 give -4.
            full = support(H, 2)
            near = support(H, 2; support_mode=:graph_distance, power_graph_radius=1)
            wide = support(H, 2; support_mode=:graph_distance)
            local2 = support(H, 2; support_mode=:local_products)
            @test near.row_counts == [28]
            @test local2.row_counts == [19]
            @test wide.family == full.family && wide.relation == :exact
            for (s, rotations) in ((near, 3), (local2, 0))
                family = [last.(only.(NC.terms.(f))) for f in s.family]
                symmetry = workflow.workflow_symmetry(H; family, spatial_mode=:none)
                @test !isempty(symmetry.generators) # Discrete global axis actions survived.
                @test length(symmetry.continuous) == rotations
                p = NC.build(NC.Problem(H), NC.Relaxation(family=s.family)).program
                @test bound(p) ≈ -4 atol=2e-6
            end
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
