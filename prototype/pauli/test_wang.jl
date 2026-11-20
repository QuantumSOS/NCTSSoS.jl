# Correctness controls, not benchmarks. Run in the locked Mosek environment.
module WangWorkflowChecks
using LinearAlgebra, Test
import NCTSSoS as NC
import MosekTools
include("PauliWorkflowRedesign.jl")
using .PauliWorkflowRedesign
using .NC: JuMP

const AXES = (ComplexF64[0 1;1 0], ComplexF64[0 -im;im 0], ComplexF64[1 0;0 -1])
const ID2 = Matrix{ComplexF64}(I,2,2)
const SOLVER = JuMP.optimizer_with_attributes(MosekTools.Optimizer,
    JuMP.MOI.Silent()=>true,"MSK_IPAR_NUM_THREADS"=>1)

function word_matrix(w,n)
    localops = fill(ID2,n)
    for i in NC.legacy_letters(w)
        localops[(i-1)÷3+1] = AXES[mod1(i,3)]
    end
    return reduce(kron,localops)
end

function check_bound(plan,H,value; atol=2e-6)
    r = solve_ground_state(plan,H;optimizer=SOLVER)
    @test r.status == JuMP.MOI.OPTIMAL
    @test r.bound ≈ value atol=atol
    return r.bound
end

function evaluate_blocks(plan,y,label)
    p = plan.program
    [Symmetric([sum((p.forms.coefs[t]*y[p.forms.ids[t]] for t in NC.term_range(p.forms,
        b.first_form + NC.tri_offset(b.side,max(i,j),min(i,j))));init=0.0)
        for i in 1:b.side,j in 1:b.side]) for b in p.blocks if b.label == label]
end

function ring(label)
    n = length(label)
    _,v = NC.pauli_variables(n)
    H = sum(a[label[i]]*a[label[mod1(i+1,n)]] for a in v for i in 1:n)/(4n)
    return H,v
end

@testset "Wang compiler" begin
    BLAS.set_num_threads(1)
    order = 2
    spec(mode;radius=2) = PauliSparsitySpec(support_mode=mode,power_graph_radius=radius,
        ts_algo=:None,rewrite_complete_local_rdm=false)

    @testset "Hand-checkable ring, relabeling, and a broken generator" begin
        for label in (1:4,[3,1,4,2]),mode in (:full,:graph_distance,:local_products)
            H,_ = ring(label)
            plan = compile_relaxation(H;order,wang_matched=true,sparsity=spec(mode))
            @test plan.reports.symmetry.reduction.translation_orders == [4]
            if mode === :local_products
                # Fewer rows; it need not reproduce the full SDP's optimum.
                r = solve_ground_state(plan,H;optimizer=SOLVER)
                @test r.status == JuMP.MOI.OPTIMAL
                @test -3/4-2e-6 <= r.bound <= -1/2+2e-6
                @test length(only(plan.reports.original_rows)) < 67
            else
                check_bound(plan,H,-1/2)
                @test length(only(plan.reports.original_rows)) == 67
            end
        end
        H,(x,_,z) = ring([3,1,4,2])
        defect = H + z[1]/16
        observable = x[1]*z[3]
        broken = compile_relaxation(defect;order,wang_matched=true,sparsity=spec(:graph_distance),observable)
        dense = compile_relaxation(defect;order,symmetry=PauliSymmetrySpec(spatial_mode=:none,
            time_reversal_gauge=false,anti_sign_realify=false,klein4_sign=false,xyz_permutation=false,
            rotated_lie_detection=false,cartan_u1_weyl=false,su2_mz0_tree=false),sparsity=spec(:full))
        reference = solve_ground_state(dense,defect;optimizer=SOLVER)
        @test reference.status == JuMP.MOI.OPTIMAL
        check_bound(broken,defect,reference.bound)
        @test broken.reports.symmetry.route != :lattice
        @test !isnothing(NC.moment_terms(broken.table,last(only(NC.terms(observable)))))
        @test_throws ArgumentError compile_relaxation(defect;order,wang_matched=true,
            strengthening=PauliStrengtheningSpec(commutator_linear_order=3))
    end

    @testset "Shared stationarity, state optimality and density moments" begin
        n = 8
        H,_ = ring([7,1,4,6,2,3,8,5])
        mat = Hermitian(sum(c*word_matrix(w,n) for (c,w) in NC.terms(H)))
        energies,states = eigen(mat)
        ψ = states[:,1]
        strength = PauliStrengtheningSpec(commutator_linear_order=7,kms_psd_order=2,rdm_positivity_k=8)
        plan = compile_relaxation(H;order,wang_matched=true,sparsity=spec(:graph_distance;radius=4),
            strengthening=strength,wang_lso_span=4)
        check_bound(plan,H,first(energies))
        wang = plan.reports.analysis.wang
        @test plan.reports.analysis.active_equality_rows == wang.lso.rows > 0
        @test wang.pso.rows == length(only(plan.reports.original_rows))-1 == 276
        @test wang.rdm.sides == [binomial(n,k) for k in 0:n÷2]
        y = [real(dot(ψ,word_matrix(l.word,n)*ψ)) for l in plan.table.labels]
        @test all(A -> abs(only(A)) < 2e-10,evaluate_blocks(plan,y,-3))
        @test minimum(eigmin.(evaluate_blocks(plan,y,-1))) >= -2e-10

        # The uniform mixture of six aligned spin product states lies at the
        # top of the XXX spectrum. Its RDM is physical but violates PSO.
        high = map(plan.table.labels) do l
            ids = NC.legacy_letters(l.word)
            isempty(ids) ? 1.0 : iseven(length(ids)) && allequal(mod.(ids,3)) ? 1/3 : 0.0
        end
        @test minimum(eigmin.(evaluate_blocks(plan,high,-1))) < -1e-3
        @test minimum(eigmin.(evaluate_blocks(plan,high,-2))) >= -2e-10

        # Compare RDM entries (not only positivity) to the independent ED ket,
        # permuted from registry bit order into the discovered region order.
        region = wang.rdm.region
        for (weight,A) in enumerate(evaluate_blocks(plan,y,-2))
            bits = filter(b -> count_ones(b)==weight-1,0:(1<<n)-1)
            amps = [ψ[1+sum(((b>>(n-k))&1)<<(n-i) for (k,i) in enumerate(region))] for b in bits]
            @test Matrix(A) ≈ real(amps*amps') atol=2e-10
        end
    end

    @testset "Relabeled square (analytic ferromagnetic control)" begin
        L = 4
        n = L^2
        label = [3,1,6,7,2,8,5,4,12,16,10,13,11,9,15,14]
        _,v = NC.pauli_variables(n)
        site(x,y) = label[mod1(x,L)+L*(mod1(y,L)-1)]
        H = -sum(a[site(x,y)]*(a[site(x+1,y)]+a[site(x,y+1)])
            for a in v for x in 1:L for y in 1:L)/(4n)
        # The lower floor −1/4 per bond is attained by the all-up state.
        plan = compile_relaxation(H;order,wang_matched=true,sparsity=spec(:graph_distance;radius=4))
        @test plan.reports.symmetry.reduction.translation_orders == [4,4]
        @test length(only(plan.reports.original_rows)) == 1129
        @test plan.reports.analysis.wang.layout.shape === :square
        check_bound(plan,H,-1/2)
    end
end
end
