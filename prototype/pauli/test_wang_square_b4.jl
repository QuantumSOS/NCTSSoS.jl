# Build-only; no optimizer or license. Run with julia --project --threads=4.
module WangSquareB4Checks
using Test, LinearAlgebra
import NCTSSoS as NC
include("PauliWorkflowRedesign.jl")
import .PauliWorkflowRedesign as WF
include("benchmark/cases.jl")

function model(L)
    _,axes = NC.pauli_variables(L^2)
    sum(a[i]*a[j] for a in axes for (i,j) in bonds("square",L))/(4L^2), axes
end

word(p) = last(only(NC.terms(p)))

@testset "Periodic radius-four box (support only)" begin
    H,(x,y,_) = model(16)
    support = WF.workflow_support(H,4;support_mode=:wang_square_b4,max_rows=155137)
    words = Set(word(p) for p in only(support.family))
    pair(i,j) = word(x[i]*y[j])
    @test length(words) == 155137
    @test pair(1,5) in words && pair(1,69) in words  # (4,0), (4,4)
    @test pair(1,6) ∉ words                          # (5,0)
    @test pair(1,205) in words                       # (-4,-4), across limbs
end

@testset "Same square B4 SDP with compact product moments" begin
    BLAS.set_num_threads(1)
    H,(x,y,z) = model(6)
    sparsity = WF.PauliSparsitySpec(support_mode=:wang_square_b4,ts_algo=:None,rewrite_complete_local_rdm=false)
    strengthening = WF.PauliStrengtheningSpec(commutator_linear_order=5,kms_psd_order=2,rdm_positivity_k=8)
    options = (;order=4,sparsity,strengthening,wang_matched=true,wang_lso_span=12,moi_form=:dualized)
    @test_throws ArgumentError WF.compile_relaxation(H+x[1]*x[2];options...,square_local_moments=true)
    @test_throws ArgumentError WF.compile_relaxation(H;options...,square_local_moments=true,observable=z[1])
    reference = WF.compile_relaxation(H;options...)
    compact = WF.compile_relaxation(H;options...,square_local_moments=true)
    @test reference.reports.original_rows == compact.reports.original_rows
    @test length(only(compact.reports.original_rows)) == 14527
    wang = compact.reports.analysis.wang
    @test (wang.pso.rows,wang.rdm.sides,wang.lso.rows) == (5778,[1,8,28,56,70],70)
    a,b = reference.standard,compact.standard
    @test (a.A.colptr,a.A.rowval,a.A.nzval,a.b,a.c,a.offset,a.cones) ==
          (b.A.colptr,b.A.rowval,b.A.nzval,b.b,b.c,b.offset,b.cones)
    @test compact.table.covered == reference.table.covered
    # Table keys may differ. Check that *physical moments* name exactly
    # the same columns, including newly chosen representative words.
    for (lhs,rhs) in ((reference.table,compact.table),(compact.table,reference.table))
        @test all(enumerate(lhs.labels)) do (id,label)
            NC.moment_terms(rhs,label.word) == ((id,1),)
        end
    end
    @test NC.moment_terms(compact.table,word(x[1]*y[2]*z[3])) == ()
    @test isnothing(NC.moment_terms(compact.table,word(prod(x))))
end
end
