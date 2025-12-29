# Sparsity Verification Tests
# ============================
# Tests that sparsity exploitation produces correct results
# compared to dense (non-sparse) computation.

using Test, NCTSSoS

@testset "Bell Inequalities with Sparsity" begin
    @testset "CHSH - Dense vs Sparse" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        
        f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]
        pop = polyopt(-f, reg)
        
        # Dense (no sparsity)
        config_dense = SolverConfig(optimizer=SOLVER, order=1, 
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # Correlative sparsity only (MF)
        config_cs = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_cs = cs_nctssos(pop, config_cs)
        
        # Term sparsity only (MMD)
        config_ts = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=NoElimination(), ts_algo=MMD())
        result_ts = cs_nctssos(pop, config_ts)
        
        expected = -2.8284
        
        @test result_dense.objective ≈ expected atol=1e-3
        @test result_cs.objective ≈ expected atol=1e-3
        @test result_ts.objective ≈ expected atol=1e-3
        
        @test result_dense.objective ≈ result_cs.objective atol=1e-4
        @test result_dense.objective ≈ result_ts.objective atol=1e-4
    end
    
    @testset "I_3322 - Dense vs Sparse" begin
        reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + 
            x[2] * (y[1] + y[2] - y[3]) + 
            x[3] * (y[1] - y[2]) - 
            x[1] - 2.0 * y[1] - y[2]
        
        pop = polyopt(-f, reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # With MF correlative sparsity
        config_mf = SolverConfig(optimizer=SOLVER, order=2,
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_mf = cs_nctssos(pop, config_mf)
        
        expected = -0.25
        
        @test result_dense.objective ≈ expected atol=1e-2
        @test result_mf.objective ≈ expected atol=1e-2
        @test result_mf.objective <= result_dense.objective + 1e-2
    end
end

@testset "Trace Polynomials with Sparsity" begin
    @testset "CHSH Trace - Dense vs Sparse" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]
        
        p = -1.0 * tr(x[1] * y[1]) - tr(x[1] * y[2]) - tr(x[2] * y[1]) + tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=1,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(tpop, config_dense)
        
        # With MaximalElimination term sparsity
        config_ts = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=NoElimination(), ts_algo=MaximalElimination())
        result_ts = cs_nctssos(tpop, config_ts)
        
        expected = -2.8284
        
        @test result_dense.objective ≈ expected atol=1e-3
        @test result_ts.objective ≈ expected atol=1e-3
        @test result_dense.objective ≈ result_ts.objective atol=1e-4
    end
    
    @testset "Covariance Trace - Dense vs Sparse" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:6)])
        x = vars[1:3]
        y = vars[4:6]
        
        cov(i, j) = tr(x[i] * y[j]) - tr(x[i]) * tr(y[j])
        p = -1.0 * (cov(1, 1) + cov(1, 2) + cov(1, 3) + cov(2, 1) + cov(2, 2) - cov(2, 3) + cov(3, 1) - cov(3, 2))
        tpop = polyopt(p * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(tpop, config_dense)
        
        # With MF + MMD
        config_sparse = SolverConfig(optimizer=SOLVER, order=2,
                                     cs_algo=MF(), ts_algo=MMD())
        result_sparse = cs_nctssos(tpop, config_sparse)
        
        expected = -5.0
        
        @test result_dense.objective ≈ expected atol=1e-3
        @test result_sparse.objective ≈ expected atol=1e-3
        @test result_dense.objective ≈ result_sparse.objective atol=1e-3
    end
end

@testset "State Polynomials with Sparsity" begin
    @testset "CHSH State - Dense vs Sparse" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
        
        sp = -ς(x[1] * y[1]) - ς(x[1] * y[2]) - ς(x[2] * y[1]) + ς(x[2] * y[2])
        spop = polyopt(sp * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=1,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(spop, config_dense)
        
        # With MMD term sparsity
        config_ts = SolverConfig(optimizer=SOLVER, order=1,
                                 cs_algo=NoElimination(), ts_algo=MMD())
        result_ts = cs_nctssos(spop, config_ts)
        
        expected = -2.8284
        
        @test result_dense.objective ≈ expected atol=1e-3
        @test result_ts.objective ≈ expected atol=1e-3
        @test result_dense.objective ≈ result_ts.objective atol=1e-4
    end
    
    @testset "Covariance State - Dense vs Sparse" begin
        reg, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])
        
        cov(a, b) = 1.0 * ς(x[a] * y[b]) - ς(x[a]) * ς(y[b])
        sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)
        spop = polyopt(sp * one(typeof(x[1])), reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(spop, config_dense)
        
        # With MF + MMD
        config_sparse = SolverConfig(optimizer=SOLVER, order=2,
                                     cs_algo=MF(), ts_algo=MMD())
        result_sparse = cs_nctssos(spop, config_sparse)
        
        expected = -5.0
        
        @test result_dense.objective ≈ expected atol=1e-2
        @test result_sparse.objective ≈ expected atol=1e-2
        @test result_dense.objective ≈ result_sparse.objective atol=1e-2
    end
end

@testset "Constrained POP with Sparsity" begin
    @testset "Ball Constraint - Dense vs Sparse" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:2)])
        
        f = 2.0 - x[1]^2 + x[1]*x[2]^2*x[1] - x[2]^2 + x[1]*x[2]*x[1]*x[2] + x[2]*x[1]*x[2]*x[1] +
            x[1]^3*x[2] + x[2]*x[1]^3 + x[1]*x[2]^3 + x[2]^3*x[1]
        
        g1 = 1.0 - x[1]^2
        g2 = 1.0 - x[2]^2
        
        pop = polyopt(f, reg; ineq_constraints=[g1, g2])
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # With MMD term sparsity
        config_ts = SolverConfig(optimizer=SOLVER, order=2,
                                 cs_algo=NoElimination(), ts_algo=MMD())
        result_ts = cs_nctssos(pop, config_ts)
        
        # Both should produce negative lower bounds for this minimization problem
        @test result_dense.objective < 0
        @test result_ts.objective < 0
        @test result_dense.objective < 5.0
        @test result_ts.objective < 5.0
    end
    
    @testset "Rosenbrock - Dense vs Sparse" begin
        n = 6
        reg, (x,) = create_noncommutative_variables([("x", 1:n)])
        
        f = Float64(n) * one(typeof(x[1]))
        for i = 2:n
            f = f + 100.0 * x[i-1]^4 - 200.0 * x[i-1]^2 * x[i] - 2.0 * x[i] + 101.0 * x[i]^2
        end
        
        pop = polyopt(f, reg)
        
        # Dense
        config_dense = SolverConfig(optimizer=SOLVER, order=2,
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        
        # With MF correlative sparsity
        config_cs = SolverConfig(optimizer=SOLVER, order=2,
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_cs = cs_nctssos(pop, config_cs)
        
        # With both sparsities
        config_both = SolverConfig(optimizer=SOLVER, order=2,
                                   cs_algo=MF(), ts_algo=MMD())
        result_both = cs_nctssos(pop, config_both)
        
        # All should give similar results (lower bounds)
        @test result_cs.objective <= result_dense.objective + 1e-2
        @test result_both.objective <= result_dense.objective + 1e-2
    end
end

@testset "Sparsity Algorithm Variants" begin
    @testset "Correlative Sparsity Algorithms" begin
        reg, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
        
        f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + 
            x[2] * (y[1] + y[2] - y[3]) + 
            x[3] * (y[1] - y[2]) - 
            x[1] - 2.0 * y[1] - y[2]
        
        pop = polyopt(-f, reg)
        
        expected = -0.25
        
        # NoElimination (dense) - should give correct result
        config_dense = SolverConfig(optimizer=SOLVER, order=2, 
                                    cs_algo=NoElimination(), ts_algo=NoElimination())
        result_dense = cs_nctssos(pop, config_dense)
        @test result_dense.objective ≈ expected atol=1e-2
        
        # MF - should give correct result (good clique tree)
        config_mf = SolverConfig(optimizer=SOLVER, order=2, 
                                 cs_algo=MF(), ts_algo=NoElimination())
        result_mf = cs_nctssos(pop, config_mf)
        @test result_mf.objective ≈ expected atol=1e-2
        
        # AsIsElimination can give looser bounds
        config_asis = SolverConfig(optimizer=SOLVER, order=2, 
                                   cs_algo=AsIsElimination(), ts_algo=NoElimination())
        result_asis = cs_nctssos(pop, config_asis)
        @test result_asis.objective < 0
    end
    
    @testset "Term Sparsity Algorithms" begin
        reg, (vars,) = create_unipotent_variables([("v", 1:4)])
        x = vars[1:2]
        y = vars[3:4]
        
        p = -1.0 * tr(x[1] * y[1]) - tr(x[1] * y[2]) - tr(x[2] * y[1]) + tr(x[2] * y[2])
        tpop = polyopt(p * one(typeof(x[1])), reg)
        
        expected = -2.8284
        
        for (name, algo) in [
            ("NoElimination", NoElimination()),
            ("MMD", MMD()),
            ("MaximalElimination", MaximalElimination())
        ]
            config = SolverConfig(optimizer=SOLVER, order=1, cs_algo=NoElimination(), ts_algo=algo)
            result = cs_nctssos(tpop, config)
            @test result.objective ≈ expected atol=1e-3
        end
    end
end
