# Standalone G4 check. Set PAULI_PRICING_MOSEK=1 with MosekTools and a license
# available, or run with Clarabel in the Julia load path.
using Test
using LinearAlgebra: BLAS, Hermitian, eigmin
import NCTSSoS as NC
include("PauliWorkflowRedesign.jl")
include("../../test/Expectations.jl")
using .PauliWorkflowRedesign
using .TestExpectations: expectations_oracle
const Workflow = PauliWorkflowRedesign

if get(ENV,"PAULI_PRICING_MOSEK","0") == "1"
    import MosekTools
    const optimizer = NC.JuMP.optimizer_with_attributes(MosekTools.Optimizer,
        "MSK_IPAR_LOG"=>0, "MSK_IPAR_NUM_THREADS"=>1)
else
    import Clarabel
    const optimizer = NC.JuMP.optimizer_with_attributes(Clarabel.Optimizer,
        "verbose"=>false, "chordal_decomposition_enable"=>false, "max_threads"=>1)
end
BLAS.set_num_threads(1)

function check_g4(L)
    n = L^2
    site(x,y) = mod(x,L) + 1 + L*mod(y,L)
    bonds(directions) = Set(minmax(site(x,y),site(x+dx,y+dy))
        for x in 0:L-1 for y in 0:L-1 for (dx,dy) in directions)
    spins = last(NC.pauli_variables(n))
    Hamiltonian(edges) = sum(a[i]*a[j] for (i,j) in edges for a in spins)
    H = Hamiltonian(bonds(((1,0),(0,1)))) + Hamiltonian(bonds(((1,1),(1,-1))))/2
    W = eltype(H.data.words)
    seed = [one(W); [last(only(NC.terms(a[i]))) for a in spins for i in 1:n]]
    shape = (L,L)
    sym = PauliSymmetrySpec(spatial_mode=:circulant_dft, rotated_lie_detection=false,
        cartan_u1_weyl=false, su2_mz0_tree=false)
    sparsity = PauliSparsitySpec(support_mode=:contiguous_lattice, rewrite_complete_local_rdm=false)
    fixed_symmetry = Workflow.workflow_symmetry(H; lattice_shape=shape, Workflow._symmetry_options(sym)...)
    objective = [(ComplexF64(c),w) for (c,w) in NC.terms(H)]
    states = NamedTuple[]
    function solve_motifs(motifs)
        if isempty(motifs)
            plan = compile_relaxation(H; order=1, symmetry=sym, sparsity, lattice_shape=shape)
            state = square_pricing_state(solve_ground_state(plan,H; optimizer))
        else
            # Minimal external seed fixture; the same compiler stages and exact
            # symmetry/term policy as the initial workflow, with explicit rows.
            rows = Set(seed)
            for motif in motifs, x in 0:L-1, y in 0:L-1, sx in (-1,1), sy in (-1,1), swap in (true,false)
                sites = sort!([site(x+sx*(swap ? j : i),y+sy*(swap ? i : j)) for (i,j) in motif])
                union!(rows,(w for w in NC.basis_rows(W,NC.Pauli(),sites,length(sites)) if NC.degree(w)==length(sites)))
            end
            rows = sort!(collect(rows))
            family = [NC.RowBlock(rows,[(1.0+0im,one(W))],NC.HPSD,0)]
            p, table, _ = Workflow.workflow_symmetry_program(family,objective,fixed_symmetry;
                spatial_mode=sym.spatial_mode, lattice_shape=shape, realify=false)
            program = Workflow.workflow_term_sparsity(Workflow._real_workflow_program(p)).program
            NC.refresh_coverage!(table,program)
            form = NC.moment_standard_form(program)
            A = form.A
            constraint_ids = BitSet([1; [j+1 for j in axes(A,2) if A.colptr[j+1] > A.colptr[j]]])
            sol = NC.solve_standard(form,NC.JuMPBackend(optimizer))
            state = (; rows, table, moments=sol.values, constraint_ids, sol.bound, sol.status)
        end
        push!(states,state)
        return state
    end
    priced = price_square(H,L,solve_motifs; batch_size=2, max_rounds=2, max_diameter=1,
        max_rows=L==2 ? 180 : 350, witness_tol=2e-7, bound_tol=2e-6)
    history = priced.history
    for r in history
        kinds = [s.cut.kind for s in r.scores if s.cut !== nothing]
        println("G4 L=$L round=$(r.round) rows=$(r.row_count) bound=$(r.bound) " *
            "candidates=$(length(r.scores)) range=$(count(==(:range),kinds)) schur=$(count(==(:schur),kinds)) " *
            "selected=$(r.selected) solve_s=$(round(r.solve_seconds; digits=3)) price_s=$(round(r.pricing_seconds; digits=3))")
    end
    println("G4 L=$L stop=$(priced.stop_reason)")
    @test length(history) >= 2
    @test !isempty(first(history).selected)
    upper = expectations_oracle("expectations/square_pauli_pricing.toml","j1j2_l$L").opt
    @test all(r -> r.status == NC.JuMP.MOI.OPTIMAL && r.bound <= upper+2e-6, history)
    @test all(r -> r.delta >= -2e-6, history[2:end])
    state = first(states)
    entry(a,b) = begin
        phase, w = NC.mul(a,b)
        terms = NC.moment_terms(state.table,w)
        @assert terms !== nothing && all(id == 1 || id in state.constraint_ids for (id,c) in terms if !iszero(c))
        NC.scale(phase,NC.moment_value(state.table,w,state.moments))
    end
    cuts = filter(s -> s.cut !== nothing,first(history).scores)
    # Compare the anchored quotient with an independent fully translated shell.
    # A large local pool isolates candidate quotienting on this small fixture.
    domain = Set(NC._gns_domain(state.table))
    callback(_) = (; state..., domain)
    orbit = price_square(H,L,callback; pricing_mode=:orbit, local_radius=2,
        local_row_limit=length(state.rows), max_diameter=1, max_rounds=0, witness_tol=2e-7)
    orbit_round = only(orbit.history)
    reference = Dict(s.word=>s.cut for s in first(history).scores)
    expected_reps = Set(Workflow._pricing_word_rep(w,L) for w in keys(reference))
    @test Set(s.word for s in orbit_round.scores) == expected_reps
    @test orbit_round.pricing.partial && orbit_round.pricing.products < first(history).pricing.products
    orbit_cuts = filter(s -> s.cut !== nothing,orbit_round.scores)
    @test !isempty(orbit_cuts)
    @test all(s -> reference[s.word] !== nothing, orbit_cuts)
    for s in Iterators.flatten((cuts,orbit_cuts))
        rows = state.rows[s.cut.indices]
        base = Hermitian([entry(a,b) for a in rows,b in rows])
        @test eigmin(base) >= -2e-7
        extended = [rows; s.word]
        @test eigmin(Hermitian([entry(a,b) for a in extended,b in extended])) < -2e-7
    end
    capped = price_square(H,L,callback; pricing_mode=:orbit, local_radius=0,
        local_row_limit=1, max_diameter=1, max_rounds=0, witness_tol=2e-7)
    @test only(capped.history).pricing.largest_pool == 1
    @test capped.stop_reason === :partial_no_violation
end

@testset "Square adaptive pricing" begin
    @testset "Only complete known minors can cut" begin
        _, (x,y,z) = NC.pauli_variables(4)
        word(p) = last(only(NC.terms(p)))
        H = z[1]*z[2]
        p = word(H)
        W = typeof(p)
        rows = [one(W),word(x[1]),word(x[2])]
        table = NC.MomentTable{W}()
        # M[2,3] is allocated but initially unmeasured. Its chosen 0 value
        # would falsely make the full minor indefinite, unlike both known minors.
        words = word.([x[1],x[2],x[1]*x[2],H,y[1]*z[2],z[1]*y[2]])
        for w in words
            NC.moment_ids!(table,w)
        end
        moments = [0.,0.,0.,0.,0.9,0.8]
        constraint_ids = BitSet([1,2,3,5,6,7])
        state = (; rows, table, moments, constraint_ids, bound=0.0, status=:fixture)
        score() = only(filter(s -> s.word == p, only(price_square(H,2,_ -> state; max_rounds=0).history).scores))
        missing = score()
        @test (missing.known_cross,missing.minors,missing.cut) == (3,2,nothing)
        push!(constraint_ids,4)
        schur = score().cut
        @test schur.kind === :schur && schur.schur ≈ 0.45
        moments[3] = 1
        range = score().cut
        @test range.kind === :range && range.range_residual ≈ 0.1/sqrt(2)
        @test range.schur == Inf && range.lambda_min < 0
    end
    @testset "G4: $(L)×$(L)" for L in (2,4)
        check_g4(L)
    end
end
