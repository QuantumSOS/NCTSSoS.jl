@testset "Physics-informed Pauli sparsity" begin
    reg, (x, y, z) = create_pauli_variables(1:4)
    config = SolverConfig(optimizer=SOLVER, order=2)
    blocks(s) = [b for ts in s.cliques_term_sparsities for b in ts[1].block_bases]
    nested(a, b) = all(old -> any(new -> issubset(old, new), blocks(b)), blocks(a))

    # A single four-site interaction is not six pair interactions.
    many = polyopt(2.5 * z[1]*z[2]*z[3]*z[4] - 0.25*x[2] + 7.0, reg)
    interactions = hamiltonian_interactions(many)
    @test Set(zip(interactions.supports, interactions.weights)) ==
        Set([([1, 2, 3, 4], 2.5), ([2], 0.25)])
    seed = compute_sparsity(many, config, PhysicsSparsity(max_region_sites=1, max_block_size=3))
    @test maximum(length, blocks(seed)) <= 3
    @test all(b -> all(w -> degree(w) <= 2, b), blocks(seed))
    seed_result = cs_nctssos(many, config;
        physics=PhysicsSparsity(max_region_sites=1, max_block_size=3), dualize=false)
    @test isapprox(seed_result.objective, 7.0 - 2.5 - 0.25; atol=2e-5)

    # Couplings, not just topology, change the first affordable merge.
    uniform = polyopt(1.0*z[1]*z[2] + z[2]*z[3] + z[3]*z[4], reg)
    defect = polyopt(1.0*z[1]*z[2] + z[2]*z[3] + 10z[3]*z[4], reg)
    u = compute_sparsity(uniform, config, PhysicsSparsity(max_region_sites=2))
    d = compute_sparsity(defect, config, PhysicsSparsity(max_region_sites=2))
    first_pair(s) = first(b for b in blocks(s) if length(b) == 16)
    @test x[1] in first_pair(u)
    @test x[3] in first_pair(d)
    @test !(x[1] in first_pair(d))

    # Independent region and row caps preserve earlier principal matrices.
    for strategy in (:cut, :radius)
        previous = compute_sparsity(defect, config,
            PhysicsSparsity(; max_region_sites=1, max_block_size=3, strategy))
        for (sites, rows) in [(2, 3), (2, 8), (3, 8), (4, 8), (4, 67)]
            current = compute_sparsity(defect, config,
                PhysicsSparsity(; max_region_sites=sites, max_block_size=rows, strategy))
            @test nested(previous, current)
            @test maximum(length, blocks(current)) <= rows
            previous = current
        end
        @test only(blocks(previous)) == get_ncbasis(reg, 2)
        disconnected = polyopt(1.0*x[1] + z[4], reg)
        endpoint = compute_sparsity(disconnected, config,
            PhysicsSparsity(; max_region_sites=4, max_block_size=67, strategy))
        @test only(blocks(endpoint)) == get_ncbasis(reg, 2)
    end

    # Frustration: each bond alone permits -1, but all three cannot be -1.
    # The exact minimum of the classical triangle is -1, not the seed's -3.
    reg3, (_, _, zz) = create_pauli_variables(1:3)
    triangle = polyopt(1.0*(zz[1]*zz[2] + zz[2]*zz[3] + zz[1]*zz[3]), reg3)
    early = cs_nctssos(triangle, config;
        physics=PhysicsSparsity(max_region_sites=1), dualize=false)
    final = cs_nctssos(triangle, config;
        physics=PhysicsSparsity(max_region_sites=3, max_block_size=37), dualize=false)
    dense = cs_nctssos(triangle, config; dualize=false)
    @test isapprox(early.objective, -3; atol=2e-5)
    @test isapprox(final.objective, -1; atol=2e-5)
    @test isapprox(final.objective, dense.objective; atol=2e-5)
    dual = cs_nctssos(triangle, config;
        physics=PhysicsSparsity(max_region_sites=3, max_block_size=37))
    @test isapprox(dual.objective, final.objective; atol=2e-5)

    # Uniform distance ties must enter together; a strong bond breaks the tie.
    weighted_triangle = polyopt(triangle.objective + 9.0*zz[1]*zz[3], reg3)
    radius = PhysicsSparsity(max_region_sites=2, strategy=:radius)
    @test maximum(length, blocks(compute_sparsity(triangle, config, radius))) == 4
    @test maximum(length, blocks(compute_sparsity(weighted_triangle, config, radius))) == 16
    radius_result = cs_nctssos(triangle, config;
        physics=PhysicsSparsity(max_region_sites=3, max_block_size=37, strategy=:radius), dualize=false)
    @test isapprox(radius_result.objective, dense.objective; atol=2e-5)

    @test_throws ArgumentError PhysicsSparsity(max_block_size=2)
    @test_throws ArgumentError compute_sparsity(many,
        SolverConfig(optimizer=SOLVER, order=1), PhysicsSparsity())
    constrained = polyopt(uniform.objective, reg; eq_constraints=[1.0*x[1]])
    @test_throws ArgumentError compute_sparsity(constrained, config, PhysicsSparsity())
end

@testset "Graph-selected joint Pauli basis" begin
    reg, ops = create_pauli_variables(1:10)
    pop = polyopt(heisenberg_chain_hamiltonian(ops),reg)
    adjacent = physics_moment_basis(pop,4)
    @test adjacent == pauli_contiguous_chain_basis(reg,4)
    @test length(adjacent) == 1201
    ranged = physics_moment_basis(pop,4;pair_distance=4)
    @test length(ranged) == 1471
    @test issubset(adjacent,ranged)
    x,y,z = ops
    @test only(monomials(x[1]*y[5])) in ranged
    @test !(only(monomials(x[1]*y[6])) in ranged)
    @test only(monomials(z[9]*x[10]*y[1]*z[2])) in ranged
    @test !(only(monomials(x[1]*y[3]*z[5])) in ranged)
    @test length(physics_moment_basis(pop,4;pair_distance=5)) == 1516
    @test_throws ArgumentError physics_moment_basis(pop,4;pair_distance=4,max_rows=1470)
    @test_throws ArgumentError physics_moment_basis(pop,5)
    # Relabeling the graph must relabel the basis, not change its selection.
    perm = [2,5,1,9,7,3,10,4,8,6]
    g = pauli_site_permutation(perm)
    relabeled = polyopt(NCTSSoS._act_polynomial(g,pop.objective),reg)
    @test physics_moment_basis(relabeled,4;pair_distance=4) ==
        sort!([last(NCTSSoS._act_monomial(g,m)) for m in ranged])

    # A 4x4 periodic square also has winding four-cycles. The graph-only
    # selector intentionally retains these, beyond the paper's plaquettes.
    square_reg, square_ops = create_pauli_variables(1:16)
    site(i,j) = mod1(i,4)+4*(mod1(j,4)-1)
    h = sum(sum(axis[site(i,j)]*(axis[site(i+1,j)]+axis[site(i,j+1)])
        for axis in square_ops) for i in 1:4,j in 1:4)/4
    square = physics_moment_basis(polyopt(h,square_reg),4;pair_distance=4)
    @test count(m->degree(m)==2,square) == 120*9
    @test count(m->degree(m)==3,square) == 96*27
    @test count(m->degree(m)==4,square) == 24*81
    @test length(square) == 5665

    # The public moment_basis interface retains joint tests: it must not turn
    # these rows back into independent two-site matrices.
    small_reg, small_ops = create_pauli_variables(1:3)
    small = polyopt(heisenberg_chain_hamiltonian(small_ops;periodic=false),small_reg)
    selected = physics_moment_basis(small,2;pair_distance=2)
    result = cs_nctssos(small,SolverConfig(optimizer=SOLVER,moment_basis=selected);dualize=false)
    @test result.objective ≈ -1 atol=2e-5
end

@testset "Geometry-ball joint Pauli basis" begin
    reg, (x, y, z) = create_pauli_variables(1:5)
    graph = path_graph(5)
    local_rows = physics_moment_basis(reg, graph; radius=0)
    rows = physics_moment_basis(reg, graph; radius=1, max_rows=160)
    @test length(local_rows) == 16
    @test length(rows) == 1 + 5*3 + 7*9 + 3*27
    @test issubset(local_rows, rows)
    # The middle ball contains disconnected endpoints, unlike cliques of G.
    @test only(monomials(x[1]*y[3])) in rows
    @test only(monomials(x[1]*y[2]*z[3])) in rows
    @test !(only(monomials(x[1]*y[4])) in rows)
    @test !(only(monomials(x[1]*y[3]*z[5])) in rows)
    @test physics_moment_basis(reg, graph; radius=2) == get_ncbasis(reg, 5)
    @test_throws ArgumentError physics_moment_basis(reg, graph; max_rows=159)
    @test_throws ArgumentError physics_moment_basis(reg, graph; radius=2, max_rows=1023)
    @test_throws ArgumentError physics_moment_basis(reg, graph; radius=-1)
    @test_throws ArgumentError physics_moment_basis(reg, path_graph(4))
    @test_throws ArgumentError physics_moment_basis(reg, SimpleDiGraph(graph))

    perm = [3, 1, 5, 2, 4]
    relabeled = SimpleGraph(5)
    for edge in edges(graph)
        add_edge!(relabeled, perm[src(edge)], perm[dst(edge)])
    end
    action = pauli_site_permutation(perm)
    @test physics_moment_basis(reg, relabeled) ==
        sort!([last(NCTSSoS._act_monomial(action, m)) for m in rows])

    disconnected = SimpleGraph(5)
    add_edge!(disconnected, 1, 2)
    add_edge!(disconnected, 3, 4)
    @test length(physics_moment_basis(reg, disconnected; radius=10)) == 34

    ring_reg, (rx, ry, rz) = create_pauli_variables(1:6)
    ring_rows = physics_moment_basis(ring_reg, cycle_graph(6))
    @test length(ring_rows) == 289
    @test only(monomials(rx[6]*ry[1]*rz[2])) in ring_rows
    # {1,3,5} is a clique in G^2, but is in no radius-one ball of C6.
    @test !(only(monomials(rx[1]*ry[3]*rz[5])) in ring_rows)
    translation = pauli_site_permutation([2, 3, 4, 5, 6, 1])
    @test sort!([last(NCTSSoS._act_monomial(translation, m)) for m in ring_rows]) == ring_rows

    # Cross-ball tests recover the CHSH bound, not four independent -1 terms.
    small_reg, (sx, _, sz) = create_pauli_variables(1:2)
    pop = polyopt(1.0*(sx[1]*sx[2] + sx[1]*sz[2] + sz[1]*sx[2] - sz[1]*sz[2]), small_reg)
    selected = physics_moment_basis(small_reg, path_graph(2); radius=0)
    result = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, moment_basis=selected); dualize=false)
    @test result.objective ≈ -2sqrt(2) atol=2e-5
end
