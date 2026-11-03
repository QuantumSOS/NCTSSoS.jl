using Test, NCTSSoS, LinearAlgebra

@testset "structured real RDM finite symmetry" begin
    h = Dict(0 => [0.0 -2.0; -2.0 0.0])
    v = zeros(2,2,2,2); v[1,1,1,1] = v[2,2,2,2] = 4.0
    options = (; nk=1,norb=2,nelec_per_cell=2,blocking=:spin,
        real_moments=true,include_one_d=true,three_index=:T1T2prime)
    original = build_pqg_moment_data(h,Dict((0,0,0,0)=>v); options...)
    spatial = [FermionicModePermutation(1=>2,2=>1,3=>4,4=>3)]
    flip = [FermionicModePermutation(1=>3,2=>4,3=>1,4=>2)]
    result = NCTSSoS.reduce_v2rdm_symmetry(original,spatial; moment_generators=flip)
    @test length(result.linear.moments) < length(original.moments)
    @test maximum(b.size for b in result.linear.psd_blocks_lin) < maximum(b.size for b in original.psd_blocks_lin)
    y = Dict(k=>sin(i) for (i,k) in enumerate(result.linear.moments))
    full = Dict(k=>(s==0 ? 0.0 : s*y[r]) for (k,(s,r)) in result.reconstruction)
    evaluate(form,values) = sum((c*values[k] for (k,c) in form); init=0.0)
    @test evaluate(original.objective_lin,full) ≈ evaluate(result.linear.objective_lin,y)
    next = 1
    for (block,bases) in zip(original.psd_blocks_lin,result.transforms)
        A = [evaluate(f,full) for f in block.entries]
        rebuilt = zeros(block.size,block.size)
        for U in bases
            B = [evaluate(f,y) for f in result.linear.psd_blocks_lin[next].entries]
            rebuilt += U'*B*U
            next += 1
        end
        @test rebuilt ≈ A atol=1e-10
    end
    broken_h = Dict(0 => [0.3 -2.0; -2.0 0.0])
    broken = build_pqg_moment_data(broken_h,Dict((0,0,0,0)=>v); options...)
    @test_throws ArgumentError NCTSSoS.reduce_v2rdm_symmetry(broken,spatial)
end
