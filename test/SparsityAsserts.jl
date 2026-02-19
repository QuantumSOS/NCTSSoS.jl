using Test
using NCTSSoS

# Shared test helpers for solver-free structural validation of sparsity results.
# Safe to include multiple times (guards against method overwrite warnings).

if !isdefined(@__MODULE__, :flatten_sizes)
    flatten_sizes(sizes) = reduce(vcat, sizes; init=Int[])
end

if !isdefined(@__MODULE__, :moment_block_sizes)
    moment_block_sizes(sparsity::NCTSSoS.SparsityResult) =
        map(sparsity.cliques_term_sparsities) do ts
            length.(ts[1].block_bases)
        end
end

if !isdefined(@__MODULE__, :flatten_moment_block_sizes)
    flatten_moment_block_sizes(sparsity::NCTSSoS.SparsityResult) = flatten_sizes(moment_block_sizes(sparsity))
end

if !isdefined(@__MODULE__, :max_moment_block_size)
    function max_moment_block_size(sparsity::NCTSSoS.SparsityResult)
        sizes = flatten_moment_block_sizes(sparsity)
        isempty(sizes) ? 0 : maximum(sizes)
    end
end

if !isdefined(@__MODULE__, :assert_term_sparsity_alignment)
    function assert_term_sparsity_alignment(sparsity::NCTSSoS.SparsityResult)
        cs = sparsity.corr_sparsity
        nclq = length(cs.cliques)

        @test length(sparsity.initial_activated_supps) == nclq
        @test length(sparsity.cliques_term_sparsities) == nclq
        @test length(cs.clq_cons) == nclq
        @test length(cs.clq_mom_mtx_bases) == nclq
        @test length(cs.clq_localizing_mtx_bases) == nclq

        for i in 1:nclq
            # moment matrix TS + one TS per localizing matrix (per assigned constraint)
            @test length(sparsity.cliques_term_sparsities[i]) == 1 + length(cs.clq_cons[i])

            # Moment matrix blocks cover the full clique basis
            mom_basis = cs.clq_mom_mtx_bases[i]
            mom_blocks = sparsity.cliques_term_sparsities[i][1].block_bases
            union_blocks = reduce(vcat, mom_blocks; init=eltype(mom_basis)[])
            @test Set(union_blocks) == Set(mom_basis)

            # Localizing blocks cover each localizing basis
            local_bases = cs.clq_localizing_mtx_bases[i]
            local_ts = sparsity.cliques_term_sparsities[i][2:end]
            @test length(local_bases) == length(local_ts)

            for (basis, ts) in zip(local_bases, local_ts)
                blocks = ts.block_bases
                union_local = reduce(vcat, blocks; init=eltype(basis)[])
                @test Set(union_local) == Set(basis)
            end
        end

        return nothing
    end
end

