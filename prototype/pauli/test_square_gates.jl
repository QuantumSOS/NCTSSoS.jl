# Standalone, with NCTSSoS and MosekTools on LOAD_PATH:
#   julia --project=MY_SOLVER_ENV prototype/pauli/test_square_gates.jl
using LinearAlgebra: Hermitian, eigmin, kron
using Random: MersenneTwister
using Test
import NCTSSoS as NC
using NCTSSoS: JuMP
using MosekTools: MosekTools

include("PauliWorkflowRedesign.jl")
const Workflow = PauliWorkflowRedesign
const MOI = JuMP.MOI

function square_gate_solve(program; tolerance=1e-11)
    model,_ = NC.to_jump(NC.moment_standard_form(program))
    JuMP.set_optimizer(model,MosekTools.Optimizer)
    JuMP.set_silent(model)
    for parameter in ("PFEAS","DFEAS","REL_GAP")
        JuMP.set_optimizer_attribute(model,"MSK_DPAR_INTPNT_CO_TOL_"*parameter,tolerance)
    end
    JuMP.set_optimizer_attribute(model,"MSK_IPAR_NUM_THREADS",1)
    JuMP.optimize!(model)
    status = JuMP.termination_status(model)
    if status==MOI.OPTIMAL
        value = JuMP.objective_value(model)
        @test abs(value-JuMP.dual_objective_value(model)) <= 1e-9
        return value
    end
    status==MOI.DUAL_INFEASIBLE && return -Inf
    error("Square gate SDP did not solve: $status ($(JuMP.raw_status(model))).")
end

function square_gate_energy_floor!(p, floor)
    for (id,c) in zip(p.objective_ids,p.objective_coefs)
        NC.push_term!(p.forms,c,id)
    end
    NC.push_term!(p.forms,-floor,1)
    NC.close_form!(p.forms)
    NC.push_block!(p,NC.NONNEGATIVE,1,1)
    p
end

function square_gate_bonds(L, steps)
    site(x,y) = mod(x,L)+1+L*mod(y,L)
    sort!(unique(minmax(site(x,y),site(x+dx,y+dy)) for x in 0:L-1 for y in 0:L-1 for (dx,dy) in steps))
end

function square_gate_model(L)
    _,xyz = NC.pauli_variables(L^2)
    edges = square_gate_bonds(L,((1,0),(0,1)))
    diagonals = square_gate_bonds(L,((1,1),(-1,1)))
    dot_bond(i,j) = sum(a[i]*a[j] for a in xyz)
    H = sum(dot_bond(i,j) for (i,j) in edges)+(1//2)*sum(dot_bond(i,j) for (i,j) in diagonals)
    (;xyz,H,edges,diagonals)
end

function square_gate_ed(edges,diagonals)
    unit = ComplexF64[1 0;0 1]
    xyz = (ComplexF64[0 1;1 0],ComplexF64[0 -im;im 0],ComplexF64[1 0;0 -1])
    bond(i,j) = sum(reduce(kron,(k in (i,j) ? a : unit for k in 1:4)) for a in xyz)
    eigmin(Hermitian(sum(bond(i,j) for (i,j) in edges)+0.5sum(bond(i,j) for (i,j) in diagonals)))
end

word(p) = last(only(NC.terms(p)))

function square_gate_gram_replay(program, table, report)
    recovered = report.gram_rows()
    @test length.(recovered)==[b.side for b in program.blocks]
    worst = 0.0
    for (block,rows) in zip(program.blocks,recovered), j in eachindex(rows), i in j:length(rows)
        f = block.first_form+NC.tri_offset(block.side,i,j)
        expected = Dict(program.forms.ids[t]=>program.forms.coefs[t] for t in NC.term_range(program.forms,f))
        replay = Dict{Int32,ComplexF64}()
        for (c,w) in NC.terms(rows[i]'*rows[j])
            for (id,a) in NC.moment_terms(table,w)
                replay[id] = get(replay,id,0im)+c*a
            end
        end
        for id in union(keys(expected),keys(replay))
            actual = get(replay,id,0im)
            block.kind==NC.PSD && (actual=real(actual))
            worst = max(worst,abs(actual-get(expected,id,0im)))
        end
    end
    @test worst <= 1e-11
    worst
end

@testset "Adaptive square gates" begin
    @testset "G1: complex vs real K gauge, including cross parity" begin
        _,(x,y,z) = NC.pauli_variables(4)
        rows = word.([one(x[1]),z[1],z[2],z[1]*z[2],y[1]])
        rng = MersenneTwister(293)
        errors = Float64[]
        for _ in 1:3
            a,b,c = rand(rng,-8:8,3)./8
            H = a*z[1]+b*z[2]+c*z[1]*z[2]
            complex,_,_ = Workflow.workflow_square_symmetry_program(H,rows,2;spin=false,spatial=false,real_gauge=false)
            real,table,_ = Workflow.workflow_square_symmetry_program(H,rows,2;spin=false,spatial=false)
            @test only(complex.blocks).kind == NC.HPSD
            block = only(real.blocks)
            @test block.kind == NC.PSD && block.side == length(rows)
            # In the K-fixed basis, Z₁†(iY₁)=X₁. Pick this nonzero
            # cross-parity matrix entry regardless of the sorted word order.
            order = sort(rows)
            j,i = sort([findfirst(==(word(z[1])),order),findfirst(==(word(y[1])),order)])
            f = block.first_form+NC.tri_offset(block.side,i,j)
            id = only(NC.moment_terms(table,word(x[1])))[1]
            @test [(real.forms.ids[t],real.forms.coefs[t]) for t in NC.term_range(real.forms,f)] == [(id,1)]
            copt,ropt = square_gate_solve.((complex,real);tolerance=1e-12)
            exact = minimum(a*s+b*t+c*s*t for s in (-1,1) for t in (-1,1))
            @test copt ≈ exact atol=1e-11 rtol=0
            @test ropt ≈ exact atol=1e-11 rtol=0
            push!(errors,abs(copt-ropt))
        end
        println("G1 absolute solver differences: ",errors)
        @test maximum(errors) <= 1e-12
    end

    @testset "G2: four actual SDP rows, independent ED, and ablations" begin
        (;H,xyz,edges,diagonals) = square_gate_model(2)
        @test length(edges)==4 && length(diagonals)==2
        d1,d2 = [sum(a[i]*a[j] for a in xyz) for (i,j) in diagonals]
        ident = one(H)
        total_spin = [sum(a) for a in xyz]
        rows = [ident,d1+d2,d1-d2,total_spin...]
        # The spin triplet has three copies of the same square; the other
        # two squares use the even and odd combinations of diagonal bonds.
        # Match Pauli coefficients exactly, independently of MOSEK's tolerance.
        exact(p) = Workflow._wf_exact(p)
        @test exact(H)+7 == (1//2)*sum(exact(s)^2 for s in total_spin) +
            (1//16)*(exact(2ident-d1-d2)^2 + exact(d1-d2)^2)
        p,table,report = Workflow.workflow_square_symmetry_program(H,rows,2;recover_gram_rows=true)
        @test report.row_dimension==6
        @test Set((s.J,s.momentum,s.side) for s in report.sectors) == Set([(0,(0,0),2),(0,(1,1),1),(1,(0,0),1)])
        @test sort([b.side for b in p.blocks])==[1,1,2]
        @test all(b->b.kind==NC.PSD,p.blocks)
        replay = square_gate_gram_replay(p,table,report)
        ed = square_gate_ed(edges,diagonals)
        full = square_gate_solve(p)
        @test ed ≈ -7 atol=1e-12
        @test full ≈ ed atol=1e-10
        no_triplet = square_gate_solve(first(Workflow.workflow_square_symmetry_program(H,rows[1:3],2)))
        # Without B2g the unconstrained four-body moment can grow
        # quadratically with the diagonal sum. An improving linear ray
        # need not exist; test a bounded restriction of this ablation.
        no_b2g_program = first(Workflow.workflow_square_symmetry_program(H,rows[[1,2,4,5,6]],2))
        no_b2g = square_gate_solve(square_gate_energy_floor!(no_b2g_program,-8))
        println("G2 SDP=$full ED=$ed no_J1=$no_triplet no_B2g_with_floor=$no_b2g max_gram_entry_error=$replay")
        @test no_triplet < full-1e-3
        @test no_b2g ≈ -8 atol=1e-9 rtol=0

        diagonal_axes = [a[i]*a[j] for (i,j) in diagonals for a in xyz]
        proposed = word.([ident;reduce(vcat,xyz);diagonal_axes])
        @test length(proposed)==19
        @test_throws ArgumentError Workflow.workflow_square_symmetry_program(H,proposed,2)
        # Actual spin closure of XX, YY, ZZ includes symmetric off-axis
        # pairs (spin 2), but does not need an antisymmetric diagonal pair.
        sym_axes = [xyz[a][i]*xyz[b][j]+xyz[b][i]*xyz[a][j]
            for (i,j) in diagonals for a in 1:2 for b in a+1:3]
        closed = [ident;reduce(vcat,xyz);diagonal_axes;sym_axes]
        q,_,qreport = Workflow.workflow_square_symmetry_program(H,closed,2)
        @test qreport.row_dimension==25
        @test any(s->s.J==2,qreport.sectors)
        qpruned,qaudit = Workflow.workflow_term_sparsity(q;prune_constant_singletons=true)
        @test qaudit.zero_lift_audit
        @test sum(b.side for b in qpruned.blocks)>4
        @test square_gate_solve(qpruned) ≈ ed atol=1e-10
    end

    @testset "G3: exact zero lift and conservative singleton pruning" begin
        (;H,xyz) = square_gate_model(3)
        rows = word.([one(H);reduce(vcat,xyz)])
        dense,_,_ = Workflow.workflow_square_symmetry_program(H,rows,3;spin=false,spatial=false,real_gauge=false)
        symmetric,table,report = Workflow.workflow_square_symmetry_program(H,rows,3;recover_gram_rows=true)
        @test any(s->s.J==1 && s.momentum!=(0,0),report.sectors)
        @test sum((2s.J+1)*s.star_size*s.irrep_copies*s.side for s in report.sectors)==length(rows)
        replay = square_gate_gram_replay(symmetric,table,report)
        pruned,closure = Workflow.workflow_term_sparsity(symmetric;prune_constant_singletons=true)
        @test closure.zero_lift_audit
        @test !isempty(closure.pruned_singletons)
        @test length(pruned.blocks)<length(symmetric.blocks)
        baseline, reduced, sparse = square_gate_solve.((dense,symmetric,pruned))
        println("G3 3×3 dense=$baseline symmetry=$reduced pruned=$sparse removed=$(closure.pruned_singletons) max_gram_entry_error=$replay")
        @test abs(baseline-reduced)<=1e-9
        @test abs(reduced-sparse)<=1e-9

        # The 1–2 edge is a bridge, pinned through the equality. The 3rd
        # row is isolated but its nonconstant square controls the objective.
        p = NC.MomentProgram{Float64}(3)
        for entries in ([(1,1)],[(3,1)],[],[(1,1)],[],[(2,1)])
            for (id,c) in entries
                NC.push_term!(p.forms,c,id)
            end
            NC.close_form!(p.forms)
        end
        NC.push_block!(p,NC.PSD,3)
        NC.push_term!(p.forms,1,3)
        NC.push_term!(p.forms,-0.5,1)
        NC.close_form!(p.forms)
        NC.push_block!(p,NC.ZERO,1,1)
        push!(p.objective_ids,2)
        push!(p.objective_coefs,1)
        q,audit = Workflow.workflow_term_sparsity(p;prune_constant_singletons=true)
        @test Set(Tuple.(audit.cliques[1]))==Set([(1,2),(3,)])
        @test isempty(audit.pruned_singletons)
        @test square_gate_solve(q) ≈ 0 atol=1e-10
    end
end
