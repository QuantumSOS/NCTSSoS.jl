using Test: @test, @testset
import NCTSSoS as NC

function check_workflow_local_density(workflow; solver)
    @testset "Audited complete local density replacement" begin
        registry, (x, y, z) = NC.pauli_variables(3)
        H = y[1]*x[2] + z[2]*x[3] + 0.2*y[3] + 0.4*x[1]*x[3]
        W = NC.word_type(registry, 2)
        rows(sites) = [NC._polynomial(registry, w) for w in NC.basis_rows(W, NC.Pauli(), sites, 2)]
        a, b = rows(1:2), rows(2:3)
        # A third overlapping block uses only part of the available algebra.
        cross = [one(x[1]), x[1], x[3]]
        wrong_phase = [a[1:end-1]; -a[end]]
        duplicate = [a[1:end-1]; a[end-1]]
        problem = NC.Problem(H; constraints=(NC.ExpectationZero(z[2]), NC.ExpectationPositive(x[2] + 0.5)))
        candidates = workflow.workflow_local_density_candidates
        @test isempty(candidates(H, [a[1:end-1], wrong_phase, duplicate]))
        for right in (b, b[1:end-1])
            family = [a, right, cross]
            audited = candidates(H, family)
            @test [c.family_index for c in audited] == (length(right) == 16 ? [1, 2] : [1])
            built = NC.build(problem, NC.Relaxation(; family, strengthen=(NC.DensityMatrix([1, 3]),)))
            rewritten = workflow.workflow_local_density(built, audited)
            short = rewritten.program
            @test length(rewritten.report.applied) == length(audited)
            @test isempty(rewritten.report.skipped)
            @test [b.side for b in short.blocks[1:3]] == (length(right) == 16 ? [4, 4, 3] : [4, 15, 3])
            retained(p) = [(b.kind, b.side, b.label, [
                [(p.forms.ids[t], p.forms.coefs[t]) for t in NC.term_range(p.forms, f)]
                for f in NC.form_range(b)]) for b in p.blocks[3:end]]
            @test retained(short) == retained(built.program)
            @test short.nvars == built.program.nvars
            @test Set(short.forms.ids) == Set(built.program.forms.ids)
            @test short.objective_ids == built.program.objective_ids
            @test short.objective_coefs == built.program.objective_coefs
            term = workflow.workflow_term_sparsity(short)
            @test all(i -> !haskey(term.report.cliques, i),
                findall(b -> b.label < 0, short.blocks))
            original = NC.solve_standard(NC.moment_standard_form(built.program), solver)
            replacement = NC.solve_standard(NC.moment_standard_form(short), solver)
            @test original.bound !== nothing && replacement.bound !== nothing
            @test original.bound ≈ replacement.bound atol=2e-6
        end
    end
end
