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
