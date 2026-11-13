using Test, NCTSSoS, LinearAlgebra

@testset "Supplied continuous generators: exact validation only" begin
    # SO(5): two free vectors, not two commuting sites or involutions.
    reg, (w,) = create_noncommutative_variables([("w", 1:10)])
    x, y = w[1:5], w[6:10]
    @test all(a * b != b * a for a in x, b in y)
    h = -sum(x[i] * y[i] + y[i] * x[i] for i in 1:5) / 2
    balls = [1 - sum(a * a for a in x), 1 - sum(a * a for a in y)]
    pop = polyopt(h, reg; ineq_constraints=balls)
    actions = [SuppliedDerivation(reg, Dict(v.word[1] =>
        (i == a || i == a+5 ? Polynomial(w[i+ b-a]) :
         i == b || i == b+5 ? -Polynomial(w[i-b+a]) : zero(h))
        for (i, v) in enumerate(w))) for a in 1:5 for b in a+1:5]
    bases = [vcat([one(w[1])], w), [one(w[1])], [one(w[1])]]
    report = validate_continuous_symmetry(pop, actions; positivity_bases=bases,
        compact_metric=Matrix{Int}(I, 10, 10))
    @test report.validation == :exact
    @test report.compactness == :verified
    @test report.reduction == :not_performed
    @test report.relation_to_input == :not_applicable
    @test report.physical_target == :unchanged
    @test report.positivity_blocks == 3
    @test all(iszero(differentiate(a, h)) for a in actions)
    @test_throws ArgumentError validate_continuous_symmetry(
        polyopt(h, reg; ineq_constraints=[1 - x[1]^2]), actions; positivity_bases=bases)
    # A closed union does not make each positivity block closed.
    @test_throws ArgumentError validate_continuous_symmetry(pop, actions;
        positivity_bases=[[one(w[1]), x[1]], w[2:end]])
    @test_throws ArgumentError validate_continuous_symmetry(pop, actions;
        positivity_bases=[[one(w[1]), x[1]]])

    # Nonorthogonal column-image representation K=P^-1 J P, metric P'P.
    # The objective is kept zero here to isolate metric and relation validation.
    p = Matrix{Rational{BigInt}}(I, 10, 10)
    p[1, 2] = p[6, 7] = 1
    j = zeros(Rational{BigInt}, 10, 10)
    j[2, 1] = j[7, 6] = 1
    j[1, 2] = j[6, 7] = -1
    k = inv(p) * j * p
    transformed = SuppliedDerivation(reg, Dict(w[c].word[1] =>
        sum(Polynomial((k[r, c], w[r])) for r in 1:10) for c in 1:10))
    zero_pop = polyopt(zero(h), reg)
    @test validate_continuous_symmetry(zero_pop, [transformed];
        positivity_bases=[w], compact_metric=p' * p).compactness == :verified
    @test_throws ArgumentError validate_continuous_symmetry(zero_pop, [transformed];
        positivity_bases=[w], compact_metric=Matrix{Int}(I, 10, 10))
    @test_throws ArgumentError validate_continuous_symmetry(zero_pop, [transformed];
        positivity_bases=[w], compact_metric=zeros(Int, 10, 10))

    # Equality spans may rotate; inequalities cannot use the same span rule.
    @test validate_continuous_symmetry(polyopt(h, reg; eq_constraints=Polynomial.(w)),
        actions; positivity_bases=bases).validation == :exact
    for epsilon in (1, 1e-12, big(1)//10^12)
        breaking = Polynomial((epsilon, x[1]))
        @test_throws ArgumentError validate_continuous_symmetry(polyopt(breaking, reg),
            actions; positivity_bases=bases)
        @test_throws ArgumentError validate_continuous_symmetry(polyopt(h, reg; eq_constraints=[breaking]),
            actions; positivity_bases=bases)
        @test_throws ArgumentError validate_continuous_symmetry(polyopt(h, reg; moment_eq_constraints=[breaking]),
            actions; positivity_bases=bases)
        @test_throws ArgumentError validate_continuous_symmetry(polyopt(h, reg; ineq_constraints=[breaking]),
            actions; positivity_bases=bases)
    end

    # Polynomial images really run through Leibniz, without a compact claim.
    free, (v,) = create_noncommutative_variables([("v", 1:2)])
    nonlinear = SuppliedDerivation(free, Dict(v[1].word[1] => v[2]^2,
        v[2].word[1] => zero(Polynomial(v[2]))))
    @test differentiate(nonlinear, v[1]^2) == v[2]^2 * v[1] + v[1] * v[2]^2
    nonlinear_pop = polyopt(v[2]^2, free)
    @test validate_continuous_symmetry(nonlinear_pop, [nonlinear];
        positivity_bases=[[one(v[1]), v[1], v[2]^2]]).compactness == :unverified
    @test_throws ArgumentError validate_continuous_symmetry(nonlinear_pop, [nonlinear];
        positivity_bases=[[one(v[1]), v[1], v[2]^2]], compact_metric=[1 0; 0 1])
    scaling = SuppliedDerivation(free, Dict(a.word[1] => a for a in v))
    @test validate_continuous_symmetry(polyopt(zero(Polynomial(v[1])), free), [scaling];
        positivity_bases=[v]).compactness == :unverified
    @test_throws ArgumentError validate_continuous_symmetry(polyopt(zero(Polynomial(v[1])), free),
        [scaling]; positivity_bases=[v], compact_metric=[1 0; 0 1])

    # Differentiate raw involution/idempotency relations, not their reduced zero.
    for create in (create_unipotent_variables, create_projector_variables)
        r, (u,) = create([("u", 1:1)])
        bad = SuppliedDerivation(r, Dict(u[1].word[1] => u[1]))
        @test_throws ArgumentError validate_continuous_symmetry(polyopt(0.0*u[1], r),
            [bad]; positivity_bases=[[one(u[1]), u[1]]])
    end
    # Cross-site commutation is a defining relation too.
    r, (a, b) = create_noncommutative_variables([("a", 1:2), ("b", 1:1)])
    bad = SuppliedDerivation(r, Dict(a[1].word[1] => zero(Polynomial(a[1])),
        a[2].word[1] => zero(Polynomial(a[1])), b[1].word[1] => a[1]))
    @test_throws ArgumentError validate_continuous_symmetry(polyopt(0.0*a[1], r),
        [bad]; positivity_bases=[vcat(a,b)])

    # Pauli rotation guards and the same-order singlet/non-singlet inputs.
    r, (sx, sy, sz) = create_pauli_variables(1:2)
    letters = vcat(sx, sy, sz)
    rotations = [SuppliedDerivation(r, Dict(v.word[1] => image for (v,image) in
        zip(letters, images))) for images in
        (vcat(zero.(Polynomial.(sx)), Polynomial.(sz), -Polynomial.(sy)),
         vcat(-Polynomial.(sz), zero.(Polynomial.(sy)), Polynomial.(sx)),
         vcat(Polynomial.(sy), -Polynomial.(sx), zero.(Polynomial.(sz))))]
    spin_h = sx[1]*sx[2] + sy[1]*sy[2] + sz[1]*sz[2]
    spin_rows = [vcat([one(sx[1])], letters)]
    for (couplings, expected) in (((1,2,3),0), ((1,1,2),1), ((1,1,1),3),
        ((1,1+big(1)//10^12,2),0))
        # Counts only these supplied rotations, not automatic Lie detection.
        xyz = sum(Polynomial((c, only(monomials(pair)))) for (c,pair) in
            zip(couplings, (sx[1]*sx[2], sy[1]*sy[2], sz[1]*sz[2])))
        @test count(a -> iszero(differentiate(a, xyz)), rotations) == expected
    end
    for problem in (polyopt(spin_h,r), polyopt(-spin_h,r),
        polyopt(-spin_h,r;eq_constraints=[spin_h-1]))
        @test validate_continuous_symmetry(problem, rotations;
            positivity_bases=spin_rows, compact_metric=Matrix{Int}(I,6,6)).validation == :exact
    end
    @test_throws ArgumentError validate_continuous_symmetry(polyopt(-(sx[1]+2sy[1]+3sz[1]),r),
        rotations; positivity_bases=spin_rows)
    @test_throws ArgumentError validate_continuous_symmetry(polyopt(spin_h,r;eq_constraints=[sz[1]*1e-12]),
        rotations; positivity_bases=spin_rows)

    # Canonical anticommutation/commutation relations and adjoint compatibility.
    for create in (create_fermionic_variables, create_bosonic_variables)
        r, (a, ad) = create(1:2)
        phase = SuppliedDerivation(r, Dict(v.word[1] => image for (v,image) in
            zip(vcat(a,ad), vcat(im .* a, -im .* ad))))
        hopping = ad[1]*a[2] + ad[2]*a[1]
        basis = [vcat([one(a[1])], a, ad)]
        @test validate_continuous_symmetry(polyopt(hopping,r), [phase];
            positivity_bases=basis, compact_metric=Matrix{Int}(I,4,4)).compactness == :verified
        @test_throws ArgumentError validate_continuous_symmetry(
            polyopt(hopping + a[1]*a[2] + ad[2]*ad[1],r), [phase]; positivity_bases=basis)
        bad = SuppliedDerivation(r, Dict(v.word[1] => im*v for v in vcat(a,ad)))
        @test_throws ArgumentError validate_continuous_symmetry(polyopt(0.0*hopping,r),
            [bad]; positivity_bases=basis)
        # This one respects adjoints but violates the defining relation's 1.
        scaling = SuppliedDerivation(r, Dict(v.word[1] => v for v in vcat(a,ad)))
        @test_throws ArgumentError validate_continuous_symmetry(polyopt(0.0*hopping,r),
            [scaling]; positivity_bases=basis)
        if create === create_bosonic_variables
            @test differentiate(phase, a[1]^20) == 20im * a[1]^20
            @test_throws ArgumentError differentiate(phase, a[1]^21)
        end
    end
end

@testset "Exact structured compact metrics" begin
    reg,(v,) = create_noncommutative_variables([("v",1:2)])
    x,y = v
    pop = polyopt(0.0*x,reg)
    action = SuppliedDerivation(reg,Dict(x.word[1]=>2y,y.word[1]=>-Polynomial(x)))
    # Unequal weights matter: K is skew for diag(2,1), but not for identity.
    for metric in ([2 0;0 1], Diagonal([2,1]), [2//7 0;0 1//7])
        @test validate_continuous_symmetry(pop,[action];positivity_bases=[v],
            compact_metric=metric).compactness == :verified
    end
    epsilon = big(1)//big(10)^40
    for metric in ([1 0;0 1], [0 0;0 0], [-2 0;0 -1], [2 0;0 0],
        [2 epsilon;epsilon 1], [2 epsilon;0 1])
        @test_throws ArgumentError validate_continuous_symmetry(pop,[action];
            positivity_bases=[v],compact_metric=metric)
    end
    @test_throws ArgumentError validate_continuous_symmetry(polyopt(epsilon*x,reg),[action];
        positivity_bases=[v],compact_metric=[2 0;0 1])
    # The dense fallback must still use conjugate transpose and certify strict
    # positivity, not merely nonnegative eigenvalues or a positive diagonal.
    rotation = SuppliedDerivation(reg,Dict(x.word[1]=>Polynomial(y),y.word[1]=>-Polynomial(x)))
    @test validate_continuous_symmetry(pop,[rotation];positivity_bases=[v],
        compact_metric=[2 im;-im 2]).compactness == :verified
    @test_throws ArgumentError validate_continuous_symmetry(pop,[rotation];
        positivity_bases=[v],compact_metric=[1 im;-im 1])
end
