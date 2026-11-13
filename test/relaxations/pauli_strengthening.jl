using LinearAlgebra

@testset "Pauli subsystem and perturbative positivity" begin
    reg, (x,y,z) = create_pauli_variables(1:2)
    config = SolverConfig(optimizer=SOLVER,order=1)
    h = -sum(axis[1]*axis[2] for axis in (x,y,z))/4
    pop = polyopt(h,reg)
    base = pauli_strengthened_relaxation(pop,config)
    rdm = pauli_strengthened_relaxation(pop,config;subsystems=[[1,2],[1]])
    @test rdm.moment_problem.total_basis == base.moment_problem.total_basis
    @test !rdm.ground_state_only

    # Independently evaluate the density matrices on a complex physical state.
    # The imaginary amplitudes detect a transposed density / wrong Y sign.
    C = Complex{Rational{BigInt}}
    eye = Matrix{C}(I,2,2)
    paulis = map(p -> C.(p), ([0 1;1 0], [0 -im;im 0], [1 0;0 -1]))
    operators = Dict(axis[i].word[1] => (i == 1 ? kron(p,eye) : kron(eye,p))
        for (axis,p) in zip((x,y,z),paulis) for i in 1:2)
    operator(m) = prod((operators[i] for i in m.word);init=Matrix{C}(I,4,4))
    psi = C[1,im,2,-im]
    rho = psi*psi'/7
    expectation(p, state=rho) = sum(c*LinearAlgebra.tr(state*operator(m)) for (c,m) in terms(p);init=zero(C))
    @test expectation.(rdm.extra.constraints[1][2]) == rho
    @test expectation.(rdm.extra.constraints[2][2]) ==
        [sum(rho[2i-2+k,2j-2+k] for k in 1:2) for i in 1:2,j in 1:2]
    @test sum(diag(rdm.extra.constraints[1][2])) == one(h)

    anisotropic = polyopt(0.5*x[1]*x[2]+0.25*y[1]*z[2]+0.125*z[1]-0.375*x[2]-0.0625*y[2],reg)
    q = [x[1],y[1],z[2]]
    full_config = SolverConfig(optimizer=SOLVER,moment_basis=get_ncbasis(reg,1))
    pert = pauli_strengthened_relaxation(anisotropic,full_config;perturbations=q)
    hm = sum(C(c)*operator(m) for (c,m) in terms(anisotropic.objective))
    qm = operator.(pert.extra.qs)
    expected = [LinearAlgebra.tr(rho*(a*hm*b'-(hm*a*b'+a*b'*hm)/2)) for a in qm,b in qm]
    @test expectation.(only(pert.extra.constraints)[2]) == expected
    @test expected == expected'
    @test pert.ground_state_only
    @test pert.moment_problem.total_basis == base.moment_problem.total_basis

    # Perturbative positivity must exclude an excited state, not be mislabeled
    # as a condition satisfied by every density matrix.
    one_reg, (ox,oy,oz) = create_pauli_variables(1:1)
    single = pauli_strengthened_relaxation(polyopt(1.0*oz[1],one_reg),config;
        perturbations=[ox[1]],stationarity=false)
    @test only(only(single.extra.constraints)[2]) == -2oz[1]

    # Stationarity is independent of the optional perturbative positive cone.
    stationary = pauli_strengthened_relaxation(polyopt(1.0*oz[1],one_reg),config;
        stationarity_operators=[oy[1]],moment_domain=:extend)
    @test isempty(stationary.extra.constraints)
    @test only(stationary.extra.equalities) == -2ox[1]
    @test stationary.ground_state_only
    stationary_rows = only(stationary.original.positive_rows)
    v = [isone(m) || m == oz[1] ? 1 : 0 for m in stationary_rows]
    checked = verify_pauli_strengthening_certificate(stationary,-1,[v*v'/2],[0])
    @test checked.identity_exact && checked.certified && checked.ground_state_only

    three_reg, (tx,ty,tz) = create_pauli_variables(1:3)
    incomplete = polyopt(1.0*tz[1]*tz[3],three_reg)
    @test_throws ArgumentError pauli_strengthened_relaxation(incomplete,config;subsystems=[[1,2,3]])
    @test_throws ArgumentError pauli_strengthened_relaxation(incomplete,config;perturbations=[tx[1],ty[2]])
    @test_throws ArgumentError pauli_strengthened_relaxation(pop,config;subsystems=[[1,1]])
    @test_throws ArgumentError pauli_strengthened_relaxation(pop,config;subsystems=[[1,2]],max_block_size=3)

    # Exact certificate: H + I/4 is the singlet projector. The normalized
    # subsystem matrix is essential to this factor and to overlap consistency.
    witness = zeros(C,4,4)
    witness[2:3,2:3] = C[1 -1;-1 1]/2
    grams = [zeros(C,7,7),witness,zeros(C,2,2)]
    certificate = verify_pauli_strengthening_certificate(rdm,-1//4,grams)
    @test certificate.certified && certificate.identity_exact
    @test certificate.certified_bound == -1//4
    @test !verify_pauli_strengthening_certificate(rdm,-1//4,
        [zeros(C,7,7),-witness,zeros(C,2,2)]).certified
    corrected = verify_pauli_strengthening_certificate(rdm,-1//4+1//100,grams)
    @test corrected.certified_bound == -1//4
    factors = [zeros(C,7,1),reshape(C[0,1,-1,0],4,1),zeros(C,2,1)]
    @test verify_pauli_strengthening_certificate(rdm,-1//4,factors;factors=true).certified

    # A certificate using the ground-state-only cone: Z + 2I =
    # (I+Z)^2 + (1/2)*(-2Z). Its positive multipliers must be checked too.
    single_rows = only(single.original.positive_rows)
    vector = [isone(m) || m == oz[1] ? 1 : 0 for m in single_rows]
    ground_certificate = verify_pauli_strengthening_certificate(single,-2,
        [vector*vector',fill(1//2,1,1)])
    @test ground_certificate.certified && ground_certificate.identity_exact
    @test ground_certificate.ground_state_only

    complex_result = pauli_strengthened_relaxation(polyopt(1.0*(ox[1]+oy[1]+oz[1]),one_reg),
        config;subsystems=[[1]])
    complex_certificate = verify_pauli_strengthening_certificate(complex_result,-3,
        [zeros(C,4,4),C[4 1-im;1+im 2]])
    @test complex_certificate.certified && complex_certificate.identity_exact
    complex_reduced = pauli_strengthened_relaxation(polyopt(1.0*(ox[1]+oy[1]+oz[1]),one_reg),
        config;subsystems=[[1]],clifford_generators=[pauli_site_permutation([1])])
    complex_grams = [zeros(C,size(complex_reduced.reduction.moment_problem.constraints[b.constraint_index][2]))
        for b in complex_reduced.reduction.blocks]
    qcomplex = C[4 1-im;1+im 2]
    push!(complex_grams,[real.(qcomplex) -imag.(qcomplex);imag.(qcomplex) real.(qcomplex)]/2)
    @test only(complex_reduced.extra_embedded)
    complex_checked = verify_pauli_strengthening_certificate(complex_reduced,-3,complex_grams)
    @test complex_checked.certified && complex_checked.identity_exact

    # On this ferromagnetic pair, the order-one moment relaxation gives -3/4;
    # density positivity restores the physical lower bound -1/4.
    bounds = Float64[]
    for result in (base,rdm)
        model, _ = build_jump_model(result.moment_problem)
        set_optimizer(model,SOLVER)
        set_silent(model)
        optimize!(model)
        @test termination_status(model) == JuMP.MOI.OPTIMAL
        push!(bounds,objective_value(model))
    end
    @test bounds[1] ≈ -0.75 atol=2e-5
    @test bounds[2] ≈ -0.25 atol=2e-5

    # Both full spin rotations and a site swap preserve this entire family.
    vectors = [[axis[i] for axis in (x,y,z)] for i in 1:2]
    actions = [SuppliedDerivation(reg,Dict(v[i].word[1] =>
        (i == a ? Polynomial(v[b]) : i == b ? -Polynomial(v[a]) : zero(Polynomial(v[i])))
        for v in vectors for i in 1:3)) for a in 1:3 for b in a+1:3]
    swap = [pauli_site_permutation([2,1])]
    reduced = pauli_strengthened_relaxation(pop,config;subsystems=[[1,2]],
        perturbations=[x;y;z],actions,clifford_generators=swap)
    model, _ = build_jump_model(reduced.moment_problem)
    set_optimizer(model,SOLVER)
    set_silent(model)
    optimize!(model)
    @test termination_status(model) == JuMP.MOI.OPTIMAL
    @test objective_value(model) ≈ -0.25 atol=2e-5
    @test reduced.moment_domain == :unchanged
    reduced_grams = [zeros(C,size(reduced.reduction.moment_problem.constraints[b.constraint_index][2]))
        for b in reduced.reduction.blocks]
    for (i,(_,mat)) in enumerate(reduced.extra.constraints)
        d = size(mat,1)
        gram = zeros(C,reduced.extra_embedded[i] ? 2d : d,reduced.extra_embedded[i] ? 2d : d)
        i == 1 && (gram[1:d,1:d] = witness)
        push!(reduced_grams,gram)
    end
    reduced_certificate = verify_pauli_strengthening_certificate(reduced,-1//4,
        reduced_grams,zeros(Int,length(reduced.extra.equalities)))
    @test reduced_certificate.certified && reduced_certificate.identity_exact
    @test reduced_certificate.symmetry_averaged
    @test_throws ArgumentError pauli_strengthened_relaxation(pop,config;
        subsystems=[[1]],clifford_generators=swap)
    @test_throws ArgumentError pauli_strengthened_relaxation(pop,config;
        perturbations=[x[1]],actions)
end

@testset "Extended shared Pauli moments" begin
    reg,ops = create_pauli_variables(1:3)
    pop = polyopt(heisenberg_chain_hamiltonian(ops),reg)
    config = SolverConfig(optimizer=SOLVER,order=1)
    base = pauli_strengthened_relaxation(pop,config)
    @test_throws ArgumentError pauli_strengthened_relaxation(pop,config;subsystems=[[1,2,3]])
    full = pauli_strengthened_relaxation(pop,config;
        subsystems=[[1,2,3],[1,3]],moment_domain=:extend)
    @test only(full.original.positive_rows) == only(base.original.positive_rows)
    @test length(full.moment_problem.total_basis) == 64
    @test allunique(full.moment_problem.total_basis)
    @test full.moment_domain == :extend
    @test !full.ground_state_only
    @test_throws ArgumentError pauli_strengthened_relaxation(pop,config;moment_domain=:typo)

    # Independent physical matrices fix signs, normalization, and shared overlap
    # variables. The asymmetric complex state detects lost off-diagonal moments.
    paulis = ([0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    eye = Matrix{ComplexF64}(I,2,2)
    operators = Dict(axis[i].word[1] => reduce(kron,[j == i ? p : eye for j in 1:3])
        for (axis,p) in zip(ops,paulis) for i in 1:3)
    operator(m) = prod((operators[i] for i in m.word);init=Matrix{ComplexF64}(I,8,8))
    psi = ComplexF64[1,im,2,-im,3,1+im,-2,2im]
    rho = psi*psi'/sum(abs2,psi)
    expectation(p) = sum(c*LinearAlgebra.tr(rho*operator(m)) for (c,m) in terms(p);init=0im)
    @test expectation.(full.extra.constraints[1][2]) ≈ rho atol=1e-14
    @test expectation.(full.extra.constraints[2][2]) ≈
        [sum(rho[4fld(a,2)+2s+mod(a,2)+1,4fld(b,2)+2s+mod(b,2)+1] for s in 0:1)
            for a in 0:3, b in 0:3] atol=1e-14

    vectors = [[axis[i] for axis in ops] for i in 1:3]
    actions = [SuppliedDerivation(reg,Dict(v[i].word[1] =>
        (i == a ? Polynomial(v[b]) : i == b ? -Polynomial(v[a]) : zero(Polynomial(v[i])))
        for v in vectors for i in 1:3)) for a in 1:3 for b in a+1:3]
    reduced = pauli_strengthened_relaxation(pop,config;
        subsystems=[[1,2,3]],moment_domain=:extend,actions)
    for result in (full,reduced)
        model,_ = build_jump_model(result.moment_problem)
        set_optimizer(model,SOLVER)
        set_silent(model)
        optimize!(model)
        @test termination_status(model) == JuMP.MOI.OPTIMAL
        @test objective_value(model) ≈ -0.75 atol=2e-5
    end
end

@testset "Short chain orbits and magnetization density blocks" begin
    reg,ops = create_pauli_variables(1:4)
    h = heisenberg_chain_hamiltonian(ops)
    pop = polyopt(h,reg)
    # Distance two includes antipodal pairs of orbit length two, alongside
    # single-site and adjacent-pair rows of orbit length four.
    basis = physics_moment_basis(pop,2;pair_distance=2)
    paulis = ([0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    eye = Matrix{ComplexF64}(I,2,2)
    operators = Dict(axis[i].word[1] => reduce(kron,[j == i ? p : eye for j in 1:4])
        for (axis,p) in zip(ops,paulis) for i in 1:4)
    operator(m) = prod((operators[i] for i in m.word);init=Matrix{ComplexF64}(I,16,16))
    # An anisotropic thermal state is invariant under the discrete symmetries
    # but not all spin rotations; this exercises more than the optimum alone.
    hm = sum(weight*operators[axis[i].word[1]]*operators[axis[mod1(i+1,4)].word[1]]
        for (axis,weight) in zip(ops,(0.7,1.3,0.2)) for i in 1:4)
    rho = exp(-hm)
    rho /= LinearAlgebra.tr(rho)
    expectation(p) = sum(c*LinearAlgebra.tr(rho*operator(m)) for (c,m) in terms(p);init=0im)
    main,report = pauli_translation_invariant_moment_relaxation(pop,ops,2;basis)
    spectrum = Float64[]
    for (label,(_,mat)) in zip(report.block_labels,main.constraints)
        values = eigvals(Hermitian(expectation.(mat)))
        append!(spectrum,values)
        label.momentum == 1 && append!(spectrum,values)
    end
    dense = [LinearAlgebra.tr(rho*(operator(a)'*operator(b))) for a in basis,b in basis]
    @test sort(spectrum) ≈ eigvals(Hermitian(dense)) atol=1e-10

    # Use a smaller main basis so four-site density positivity actually adds
    # variables. All magnetization sectors are retained, not a singlet sector.
    rows = physics_moment_basis(pop,1)
    xyz = only(monomials(ops[1][1]*ops[2][2]*ops[3][3]))
    bounds = Float64[]
    for sign in (1,-1), mode in (:none,:magnetization)
        result = pauli_chain_strengthened_relaxation(polyopt(sign*h,reg),ops;
            basis=rows,subsystems=[collect(1:4)],stationarity_operators=[xyz],
            moment_domain=:extend,subsystem_symmetry=mode)
        @test isempty(result.extra.qs)
        @test !isempty(result.extra.equalities)
        @test result.ground_state_only
        @test length(result.original.positive_rows[1]) == length(rows)
        if mode == :magnetization
            @test sort([size(mat,1) for (_,mat) in result.extra.constraints]) == [1,1,4,4,6]
        end
        model,_ = build_jump_model(result.moment_problem)
        set_optimizer(model,SOLVER)
        set_silent(model)
        optimize!(model)
        @test termination_status(model) == JuMP.MOI.OPTIMAL
        target = sign == 1 ? -2.0 : -1.0
        @test objective_value(model) ≈ target atol=2e-5
        push!(bounds,objective_value(model))
    end
    @test bounds[1] ≈ bounds[2] atol=2e-5
    @test bounds[3] ≈ bounds[4] atol=2e-5
    @test_throws ArgumentError pauli_chain_strengthened_relaxation(pop,ops;
        basis=rows,subsystems=[collect(1:4)],moment_domain=:extend,
        subsystem_symmetry=:magnetization,max_block_size=5)
end

@testset "Historical chain stationarity selection" begin
    reg,ops = create_pauli_variables(1:10)
    pop = polyopt(heisenberg_chain_hamiltonian(ops),reg)
    basis = physics_moment_basis(pop,4;pair_distance=4)
    selected = pauli_chain_stationarity_operators(pop,ops;basis,subsystems=[collect(1:8)])
    # Counts from executing QMBCertify b18830a's historical selector.
    @test selected.candidate_count == 182
    @test length(selected.operators) == 50
    @test selected.domain_size == 641
    @test (selected.unsupported,selected.zero_equations,selected.duplicates) == (0,24,108)
    @test only(monomials(ops[1][1]*ops[2][2]*ops[3][4])) in selected.operators
    @test !(only(monomials(ops[1][1]*ops[2][2]*ops[3][3])) in selected.operators)
    @test allunique(selected.operators)
    @test length(selected.equations) == length(selected.operators)
    @test issorted(selected.moments) && allunique(selected.moments)
    @test length(selected.moments) == selected.domain_size
    @test all(!iszero(eq) && all(m in selected.moments for m in monomials(eq))
        for eq in selected.equations)
    # Independent Pauli multiplication fixes the sign and factor in -i[H,Q],
    # not just which candidates survive. H has bond coefficient 1/4 here.
    q = only(monomials(ops[1][1]*ops[2][2]*ops[3][4]))
    M = eltype(first(ops))
    T = eltype(first(first(ops)).word)
    expected_equation = Polynomial([(c,M(T.(word))) for (c,word) in (
        (-1//2,[1,4,8,11]), (1//2,[1,4,8,26]), (1//2,[1,4,8,29]),
        (-1//2,[1,5,7,14]), (1//2,[1,5,10,14]), (-1//2,[1,5,11,13]),
        (1//2,[1,7]), (-1//2,[1,10]))])
    @test selected.equations[findfirst(==(q),selected.operators)] == expected_equation

    # First-order rows support only same-axis pairs after the global sign
    # projection. Axis permutations identify XX/YY/ZZ; translations and mirror
    # identify distances d and 8-d. This enumerates the domain independently.
    small_reg,small_ops = create_pauli_variables(1:8)
    small_pop = polyopt(heisenberg_chain_hamiltonian(small_ops),small_reg)
    x = first(small_ops)
    M = eltype(x)
    expected = [one(M); [only(monomials(x[1]*x[j])) for j in 2:5]]
    closed = pauli_chain_stationarity_operators(small_pop,small_ops;
        basis=physics_moment_basis(small_pop,1))
    @test Set(closed.moments) == Set(expected)
    # This asymmetric row set is NOT translation closed. Orbit-pair pruning
    # without a closure check would invent pairs absent from the original basis.
    incomplete = pauli_chain_stationarity_operators(small_pop,small_ops;
        basis=[one(M),x[1],x[2],x[4]])
    @test Set(incomplete.moments) == Set(expected[1:4])
    @test_throws ArgumentError pauli_chain_stationarity_operators(
        polyopt(small_pop.objective+x[1],small_reg),small_ops;
        basis=physics_moment_basis(small_pop,1))
end

@testset "Complete spin and lattice reduction" begin
    using Graphs: cycle_graph
    reg,ops = create_pauli_variables(1:5)
    basis = physics_moment_basis(reg,cycle_graph(5);radius=1)
    patches = [sort(mod1.(i:i+3,5)) for i in 1:5]
    qs = [reduce(vcat,ops); [only(monomials(a[i]*b[mod1(i+1,5)]))
        for i in 1:5 for a in ops for b in ops]]
    h = heisenberg_chain_hamiltonian(ops)
    paulis = ([0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    matrix = sum(reduce(kron,[j in (i,mod1(i+1,5)) ? p : Matrix{ComplexF64}(I,2,2)
        for j in 1:5]) for i in 1:5 for p in paulis)/4
    # Odd rings cannot be singlets; ferromagnetic ground states have maximal
    # spin. Both controls catch an accidental total-spin-sector restriction.
    for sign in (1,-1), compressed in (false,true)
        result = pauli_chain_strengthened_relaxation(polyopt(sign*h,reg),ops;
            basis,subsystems=patches,perturbations=qs,compress=compressed)
        @test result.continuous_spin
        @test result.coefficient_error <= 1e-9
        @test length(result.spin.words) < length(result.spin.index)
        model,_ = build_jump_model(result.moment_problem)
        set_optimizer(model,SOLVER)
        set_silent(model)
        optimize!(model)
        @test termination_status(model) == JuMP.MOI.OPTIMAL
        @test objective_value(model) ≈ eigmin(Hermitian(sign*matrix)) atol=2e-5
        if sign == 1 && compressed
            # Arbitrary invariant moment assignments, not just the optimum,
            # must annihilate every infinitesimal spin rotation exactly.
            coordinates = collect(1:length(result.spin.words))
            values = result.spin.kernel*coordinates
            expectation(p) = sum(c*values[result.spin.index[m]]
                for (c,m) in terms(result.spin.discrete(p));init=0//1)
            @test all(iszero(expectation(differentiate(a,m)))
                for a in result.spin.actions for m in result.original.words)
        end
    end
    @test_throws ArgumentError pauli_chain_strengthened_relaxation(
        polyopt(h+0.1*first(ops)[1],reg),ops;basis)
    @test_throws ArgumentError pauli_chain_strengthened_relaxation(
        polyopt(h,reg),ops;basis,subsystems=[[1,2,3,4]])
end
