using Test, NCTSSoS

flatten_moment_block_sizes(sparsity::NCTSSoS.SparsityResult) =
    reduce(vcat, [length.(ts_vec[1].block_bases) for ts_vec in sparsity.cliques_term_sparsities]; init=Int[])

function _chsh_unipotent_pop()
    reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
    f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
    return polyopt(-f, reg)
end

@testset "Term Sparsity (Structural, Solver-Free)" begin
    @testset "CHSH (unipotent, order=1): blocks + nuniq" begin
        pop = _chsh_unipotent_pop()

        cases = [
            ("Dense", NoElimination(), NoElimination(), [5], 11),
            ("CS-only (MF)", MF(), NoElimination(), [4, 4], 10),
            ("TS-only (MMD)", NoElimination(), MMD(), [3, 3, 1], 6),
        ]

        for (name, cs_algo, ts_algo, expected_blocks, expected_nuniq) in cases
            @testset "$name" begin
                config = SolverConfig(optimizer=nothing, order=1, cs_algo=cs_algo, ts_algo=ts_algo)
                sparsity = compute_sparsity(pop, config)

                blocks = flatten_moment_block_sizes(sparsity)
                @test sort(blocks) == sort(expected_blocks)

                mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
                @test mp.n_unique_moment_matrix_elements == expected_nuniq
            end
        end
    end

    @testset "iterate_term_sparse_supp: connected components vs dense" begin
        registry, (x,) = create_noncommutative_variables([("x", 1:2)])
        T = eltype(indices(registry))
        M = NormalMonomial{NonCommutativeAlgebra,T}
        P = Polynomial{NonCommutativeAlgebra,T,Float64}

        x_idx = [x[i].word[1] for i in 1:2]
        basis = M[one(M), M([x_idx[1]]), M([x_idx[2]])]
        activated = M[one(M), M([x_idx[1]])]

        ts_cc = NCTSSoS.iterate_term_sparse_supp(activated, one(P), basis, MaximalElimination())
        @test sort(length.(ts_cc.block_bases)) == [1, 2]

        ts_dense = NCTSSoS.iterate_term_sparse_supp(activated, one(P), basis, NoElimination())
        @test length.(ts_dense.block_bases) == [3]
    end
end
