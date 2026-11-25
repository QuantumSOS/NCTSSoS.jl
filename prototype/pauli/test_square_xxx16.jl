# julia --project --threads=16,0 prototype/pauli/test_square_xxx16.jl [--full]
using Test, SparseArrays, SHA
using LinearAlgebra: BLAS
import NCTSSoS
include("PauliWorkflowRedesign.jl")
using .PauliWorkflowRedesign
include("benchmark/square_moment.jl")

function square_model(L)
    _, axes = NCTSSoS.pauli_variables(L^2)
    east(i) = i % L == 0 ? i-L+1 : i+1
    north(i) = mod1(i+L,L^2)
    H = sum(a[i]*(a[east(i)]+a[north(i)]) for i in 1:L^2 for a in axes)/(4L^2)
    return H,axes
end

@testset "Independent acceptance gates" begin
    moi = NCTSSoS.JuMP.MOI
    status = (;status=moi.OPTIMAL,primal=moi.FEASIBLE_POINT,dual=moi.FEASIBLE_POINT)
    # KKT at x=z=0, and two controlled changes that break different gates.
    form = NCTSSoS.MomentStandardForm(sparse([1],[1],[1.0]),[0.0],[0.0],0.0,
        [NCTSSoS.DataCone(:psd,1,1,1,1)])
    results(x) = square_point_diagnostics!(Dict{String,Any}(),form,[x],[0.0])
    r = results(2e-8)
    @test r["original_primal_psd_violation"] == 2e-8
    @test square_acceptance!(r,status;strict=false)
    @test !square_acceptance!(r,status;strict=true)
    @test r["status_accepted"] && !r["strict_original_1e-8"]
    r = results(1e-8)
    @test square_acceptance!(r,status;strict=true)
    @test !square_acceptance!(r,(;status...,status=moi.SLOW_PROGRESS);strict=true)
    @test !r["status_accepted"] && r["strict_original_1e-8"]
    delete!(r,"original_balance_inf")
    @test !square_acceptance!(r,status;strict=true)
end

@testset "Witness is opt-in and requires coverage" begin
    H,_ = square_model(4)
    plan = compile_minimal_square(H,4)
    @test count(c -> c.kind == :psd && c.side == 16,plan.standard.cones) == 1
    @test_throws ArgumentError compile_square_xxx_witness(H,16)
    @test_throws ArgumentError compile_minimal_square(H,4;physical_witness3x3=true)
    @test_throws ArgumentError compile_minimal_square(H,4;physical_witness3x3=true,recover_folded_gram=true)
end

for L in ("--full" in ARGS ? (4,16) : (4,))
    @testset "Projected physical witness on $(L)×$(L)" begin
        H,axes = square_model(L)
        # Reproduce the historical decomposition. A different BLAS reduction
        # order can legitimately retain an extra Hermitian sector.
        previous_threads = BLAS.get_num_threads()
        plan = try
            L == 16 && BLAS.set_num_threads(16)
            compile_square_xxx_witness(H,L)
        finally
            BLAS.set_num_threads(previous_threads)
        end
        f = plan.standard
        if L == 16 # Independent Tencent form, not a structural extrapolation.
            @test length(plan.rows) == 56_065
            @test size(f.A) == (28_878,12_143)
            @test length(f.cones) == 4_152
            @test nnz(f.A) == 129_897_180
            @test nnz(f.A[end,:]) == 158
        end
        @test (f.cones[end].kind,f.cones[end].side,f.cones[end].first,f.cones[end].length) ==
            (:psd,1,length(f.b),1)
        cols,values = findnz(f.A[end,:])

        path = joinpath(@__DIR__,"square_witness3x3.tsv")
        @test bytes2hex(open(sha256,path)) == "5f816cefdb7cdee4015083dd91a1fa3497df5622c09eb96cb399422f0f72bde7"

        # Build the physical polynomial before projecting it, using the public
        # Pauli arithmetic and moment lookup instead of the row builder.
        tile = [1+i+L*j for j in 0:2 for i in 0:2]
        dotpair(i,j) = sum(a[tile[i]]*a[tile[j]] for a in axes)
        physical = open(path) do io
            readline(io)
            sum(eachline(io)) do line
                fields = parse.(Int,split(line))
                pairs = Iterators.partition(fields[1:end-1],2)
                fields[end]*prod(dotpair(i,j) for (i,j) in pairs if i != 0)
            end
        end
        denom = big(68_719_476_736)
        expected = Dict{Int,Rational{BigInt}}(1 => 7_607_443_585//denom)
        for (coef,w) in NCTSSoS.terms(physical), (id,c) in NCTSSoS.moment_terms(plan.table,w)
            expected[id] = get(expected,id,0) + BigInt(real(coef))*real(c)/denom
        end
        constant = pop!(expected,1)
        filter!(p -> !iszero(last(p)),expected)
        @test Dict(cols .+ 1 .=> -values) == Dict(id => Float64(c) for (id,c) in expected)
        err = sum(abs(expected[id+1]+Rational{BigInt}(v)) for (id,v) in zip(cols,values))
        @test f.b[end] == Float64(constant+err,RoundUp)
    end
end
