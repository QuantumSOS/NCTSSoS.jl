using Test, NCTSSoS, LinearAlgebra, SparseArrays
using JuMP, COSMO
using NCTSSoS: invariant_moment_map, reconstruct_continuous_moments, continuous_moment_coefficients
using NCTSSoS: continuous_moment_relaxation

function supplied_vector_rotations(reg, vectors)
    n = length(first(vectors))
    letters = reduce(vcat, vectors)
    z = zero(Polynomial(first(letters)))
    return [SuppliedDerivation(reg, Dict(v[i].word[1] =>
        (i == a ? Polynomial(v[b]) : i == b ? -Polynomial(v[a]) : z)
        for v in vectors for i in 1:n)) for a in 1:n for b in a+1:n]
end

@testset "Shared exact invariant moment map" begin
    reg, (x,y,z) = create_pauli_variables(1:2)
    actions = supplied_vector_rotations(reg, [[x[i],y[i],z[i]] for i in 1:2])
    words = get_ncbasis(reg, 2)
    map = invariant_moment_map(actions, words)
    @test size(map.reconstruction) == (16,2)
    @test issparse(map.reconstruction)
    h = x[1]*x[2] + y[1]*y[2] + z[1]*z[2]
    @test continuous_moment_coefficients(map, h) == [0,3]
    for c in (-1//1, 1//3)
        moments = reconstruct_continuous_moments(map, [1,c])
        @test moments[one(x[1])] == 1
        @test all(moments[v] == 0 for v in vcat(x,y,z))
        for (a, b) in Iterators.product((x,y,z), (x,y,z))
            @test moments[only(monomials(a[1]*b[2]))] == (a == b ? c : 0)
        end
    end
    @test_throws ArgumentError invariant_moment_map(actions, [one(x[1]), x[1]])

    free, (w,) = create_noncommutative_variables([("w",1:10)])
    x,y = w[1:5],w[6:10]
    rotations = supplied_vector_rotations(free, [x,y])
    map = invariant_moment_map(rotations, get_ncbasis(free,2))
    @test size(map.reconstruction) == (111,5)
    # b and c stay independent until selecting the real moment section.
    moments = reconstruct_continuous_moments(map, [1,2,3,5,7])
    @test moments[only(monomials(x[5]*y[5]))] != moments[only(monomials(y[5]*x[5]))]
    for a in rotations, word in map.words
        @test all(iszero, continuous_moment_coefficients(map, differentiate(a, word)))
    end
end

@testset "Actual continuous cone assembly" begin
    reg, (x,y,z) = create_pauli_variables(1:2)
    actions = supplied_vector_rotations(reg, [[x[i],y[i],z[i]] for i in 1:2])
    h = x[1]*x[2] + y[1]*y[2] + z[1]*z[2]
    config = SolverConfig(optimizer=nothing, order=1)
    reduced = continuous_moment_relaxation(polyopt(h,reg),config,actions; compact_metric=Matrix{Int}(I,6,6))
    @test [b.size for b in reduced.moment_problem.linear.psd_blocks_lin] == [1,2]
    @test length(reduced.moment_problem.linear.moments) == 2
    @test reduced.relation_to_input == :equivalent
    @test reduced.physical_target == :unchanged
    @test [b.dimension for b in reduced.blocks] == [1,3]
    exact_grams = [zeros(Rational{BigInt},1,1), fill(big(3)//2,2,2)]
    certificate = verify_continuous_certificate(reduced,-3,exact_grams)
    @test certificate.identity_exact && certificate.grams_psd_exact
    @test certificate.certified_bound == -3
    wrong = verify_continuous_certificate(reduced,-3,[zeros(Int,1,1),ones(Int,2,2)])
    @test !wrong.identity_exact
    @test wrong.certified_bound < -3
    @test !verify_continuous_certificate(reduced,-3,[zeros(Int,1,1),[1 2;2 1]]).certified
    for c in (-1//1,1//3)
        moments = reconstruct_continuous_moments(reduced,[1,c])
        rebuilt = reconstruct_continuous_matrices(reduced,[1,c])
        for (idx,mat) in rebuilt
            direct = map(p -> sum(a*moments[m] for (a,m) in terms(p);init=0),reduced.original.constraints[idx][2])
            @test mat == direct
            @test minimum(eigvals(Hermitian(ComplexF64.(mat)))) > -1e-12
        end
    end
    invariant_only = continuous_moment_relaxation(polyopt(h,reg),config,actions;
        compact_metric=Matrix{Int}(I,6,6), compress=false)
    @test [b.size for b in invariant_only.moment_problem.linear.psd_blocks_lin] == [7]

    # Independent identity: h+3 = sum((sigma_i[1]+sigma_i[2])^2)/2.
    # Traceless single-site diagonal noise averages to zero under rotations.
    rows = only(invariant_only.original.positive_rows)
    known = zeros(Rational{BigInt},length(rows),length(rows))
    for axis in (x,y,z)
        v = [m in axis ? 1 : 0 for m in rows]
        known += (big(1)//2)*(v*v')
    end
    noise = zeros(Int,size(known))
    ix,iy = findfirst(==(x[1]),rows),findfirst(==(y[1]),rows)
    noise[ix,ix],noise[iy,iy] = 1,-1
    for scale in (1,10)
        gram = known+scale*noise
        @test gram[iy,iy] < 0 # A negative diagonal rules out positive semidefiniteness.
        averaged = verify_continuous_certificate(invariant_only,-3,[gram])
        @test averaged.identity_exact && averaged.grams_psd_exact
        @test averaged.certified_bound == -3
        @test only(values(averaged.grams)) == known
    end
    @test !verify_continuous_certificate(invariant_only,-3,[known-2I]).certified

    free, (w,) = create_noncommutative_variables([("w",1:10)])
    x,y = w[1:5],w[6:10]
    actions = supplied_vector_rotations(free,[x,y])
    h = -sum(x[i]*y[i]+y[i]*x[i] for i in 1:5)/2
    pop = polyopt(h,free; ineq_constraints=[1-sum(t^2 for t in x),1-sum(t^2 for t in y)])
    reduced = continuous_moment_relaxation(pop,config,actions; compact_metric=Matrix{Int}(I,10,10))
    @test [b.size for b in reduced.moment_problem.linear.psd_blocks_lin] == [1,2,1,1]
    @test size(reduced.moment_map.reconstruction) == (111,5)
    @test size(reduced.section,2) == 4
    @test [b.dimension for b in reduced.blocks] == [1,5,1,1]
    certificate = verify_continuous_certificate(reduced,-1,
        [zeros(Int,1,1),(big(5)//2)*[1 -1;-1 1],fill(1//2,1,1),fill(1//2,1,1)])
    @test certificate.identity_exact && certificate.grams_psd_exact
    @test certificate.certified_bound == -1
    @test verify_continuous_certificate(reduced,-1+1e-12,
        [zeros(Int,1,1),(big(5)//2)*[1 -1;-1 1],fill(1//2,1,1),fill(1//2,1,1)]).certified_bound === nothing
    moments = reconstruct_continuous_moments(reduced,[1,1//5,1//5,1//5])
    matrices = reconstruct_continuous_matrices(reduced,[1,1//5,1//5,1//5])
    @test matrices[1][2:end,2:end] == kron(fill(1//5,2,2),Matrix{Int}(I,5,5))
    for (idx,mat) in matrices
        @test mat == map(p -> sum(a*moments[m] for (a,m) in terms(p);init=0),reduced.original.constraints[idx][2])
    end
    solver = optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8,"max_iter"=>50000)
    for dualize in (false,true)
        solved = NCTSSoS.solve_sdp(reduced.moment_problem,solver;dualize)
        @test solved.objective ≈ -1 atol=2e-5
    end
end

@testset "Supplied search and same-relaxation solves" begin
    solver = optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8,"max_iter"=>50000)
    config = SolverConfig(optimizer=solver,order=1)
    reg,(x,y,z) = create_pauli_variables(1:2)
    candidates = supplied_vector_rotations(reg,[[x[i],y[i],z[i]] for i in 1:2])
    rows = [vcat([one(x[1])],x,y,z)]
    metric = Matrix{Int}(I,6,6)
    for (couplings,dimension) in (((1,2,3),0),((1,1,2),1),((1,1,1),3),((1,1+big(1)//10^12,2),0))
        h = sum(Polynomial((big(c)//1,only(monomials(pair)))) for (c,pair) in
            zip(couplings,(x[1]*x[2],y[1]*y[2],z[1]*z[2])))
        search = detect_continuous_symmetries(polyopt(h,reg),candidates;positivity_bases=rows,compact_metric=metric)
        @test search.dimension == dimension
        @test search.search_class == :real_span_of_supplied_compact_generators
    end
    h = x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
    duplicate_search = detect_continuous_symmetries(polyopt(h,reg),[candidates;candidates[1]];
        positivity_bases=rows,compact_metric=metric)
    @test duplicate_search.dimension == 3
    for (pop,expected) in ((polyopt(h,reg),-3),(polyopt(-h,reg),-3),(polyopt(-h,reg;eq_constraints=[h-1]),-1))
        plain = cs_nctssos(pop,config;dualize=false)
        result = continuous_nctssos(pop,config,candidates;compact_metric=metric)
        @test plain.objective ≈ expected atol=2e-5
        @test result.objective ≈ expected atol=2e-5
        dual = NCTSSoS.sos_dualize(result.reduction.moment_problem)
        set_optimizer(dual.model,solver)
        set_silent(dual.model)
        optimize!(dual.model)
        @test objective_value(dual.model) ≈ expected atol=2e-5
        checked = verify_continuous_certificate(result.reduction,objective_value(dual.model),
            [value.(q) for q in dual.model[:psd_duals]],value.(dual.model[:zero_duals]))
        @test Float64(checked.residual_coefficient_bound) < 2e-5
        if expected == -1
            exact = verify_continuous_certificate(result.reduction,-1,
                [zeros(Int,b.size,b.size) for b in result.reduction.moment_problem.linear.psd_blocks_lin],[-1])
            @test exact.identity_exact && exact.certified_bound == -1
            @test exact.equality_coefficients == [-1]
        end
    end
    for epsilon in (1.0,1e-12)
        @test_throws ArgumentError continuous_moment_relaxation(polyopt(-epsilon*(x[1]+2y[1]+3z[1]),reg),config,candidates;compact_metric=metric)
        @test_throws ArgumentError continuous_moment_relaxation(polyopt(h,reg;eq_constraints=[epsilon*z[1]]),config,candidates;compact_metric=metric)
    end
    @test_throws ArgumentError continuous_moment_relaxation(polyopt(h,reg),
        SolverConfig(optimizer=nothing,moment_basis=[one(x[1]),x[1],x[2]]),candidates;compact_metric=metric)

    reg,(a,ad) = create_fermionic_variables(1:2)
    phases = [SuppliedDerivation(reg,Dict(v.word[1] =>
        (abs(v.word[1]) == abs(a[mode].word[1]) ? (v.word[1] > 0 ? im*v : -im*v) : zero(Polynomial(v)))
        for v in vcat(a,ad))) for mode in 1:2]
    rows = [vcat([one(a[1])],a,ad)]
    hopping = ad[1]*a[2]+ad[2]*a[1]
    pairing = a[1]*a[2]+ad[2]*ad[1]
    for (h,dimension) in ((hopping,1),(hopping+pairing,0),(pairing,1))
        found = detect_continuous_symmetries(polyopt(h,reg),phases;positivity_bases=rows,compact_metric=Matrix{Int}(I,4,4))
        @test found.dimension == dimension
    end
    common_phase = detect_continuous_symmetries(polyopt(hopping,reg),phases;
        positivity_bases=rows,compact_metric=Matrix{Int}(I,4,4)).generators
    for dualize in (false,true)
        solved = continuous_nctssos(polyopt(hopping,reg),config,common_phase;
            compact_metric=Matrix{Int}(I,4,4),dualize)
        @test solved.objective ≈ -1 atol=2e-5
        @test sum(b.dimension*b.multiplicity for b in solved.reduction.blocks) == 5
        @test any(b.embedded for b in solved.reduction.blocks)
        if dualize
            @test Float64(solved.certificate.residual_coefficient_bound) < 2e-5
        else
            @test size(solved.matrices[1]) == (5,5)
            for (idx,mat) in solved.matrices
                direct = map(p -> sum(c*solved.moments[m] for (c,m) in terms(p);init=0),
                    solved.reduction.original.constraints[idx][2])
                @test mat ≈ direct atol=1e-12
            end
        end
    end
end

@testset "Nonorthogonal SO(5) and unbounded control" begin
    reg,(w,) = create_noncommutative_variables([("w",1:10)])
    metric = Matrix{Rational{BigInt}}(I,10,10)
    frame = copy(metric)
    frame[1,2] = frame[6,7] = 1
    metric = frame'*frame
    physical = inv(frame)'*w
    x,y = physical[1:5],physical[6:10]
    actions = SuppliedDerivation[]
    for a in 1:5, b in a+1:5
        rotation = zeros(Rational{BigInt},10,10)
        for offset in (0,5)
            rotation[offset+b,offset+a] = 1
            rotation[offset+a,offset+b] = -1
        end
        k = inv(frame)*rotation*frame
        push!(actions,SuppliedDerivation(reg,Dict(w[j].word[1] =>
            sum(k[i,j]*w[i] for i in 1:10) for j in 1:10)))
    end
    h = -sum(x[i]*y[i]+y[i]*x[i] for i in 1:5)/2
    balls = [1-sum(t^2 for t in x),1-sum(t^2 for t in y)]
    solver = optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,"eps_abs"=>1e-8,"eps_rel"=>1e-8,"max_iter"=>50000)
    config = SolverConfig(optimizer=solver,order=1)
    solved = continuous_nctssos(polyopt(h,reg;ineq_constraints=balls),config,actions;compact_metric=metric)
    @test solved.objective ≈ -1 atol=2e-5
    @test [b.size for b in solved.reduction.moment_problem.linear.psd_blocks_lin] == [1,2,1,1]
    @test any(b.metric != Matrix{Int}(I,b.dimension,b.dimension) for b in solved.reduction.blocks)
    for (idx,mat) in solved.matrices
        @test mat ≈ map(p -> sum(c*solved.moments[m] for (c,m) in terms(p);init=0),
            solved.reduction.original.constraints[idx][2]) atol=1e-12
    end
    # Independently known identity: h+1 = sum((x-y)^2)/2 + (ball_x+ball_y)/2.
    known = zeros(Rational{BigInt},11,11)
    eye5 = Matrix{Rational{BigInt}}(I,5,5)
    known[2:end,2:end] = inv(frame)*[eye5 -eye5;-eye5 eye5]*inv(frame)'/2
    @test solved.reduction.original.positive_rows[1] == [one(first(w));w]
    exact_grams = Matrix{Complex{Rational{BigInt}}}[]
    for idx in solved.reduction.original.positive_indices
        blocks = filter(b -> b.original_index == idx,solved.reduction.blocks)
        inverse = inv(hcat([b.columns for b in blocks]...))
        q = inverse*(idx == 1 ? known : fill(1//2,1,1))*inverse'
        offset = 0
        for block in blocks
            d,m = block.dimension,block.multiplicity
            push!(exact_grams,[LinearAlgebra.tr(q[offset+(i-1)*d+1:offset+i*d,offset+(j-1)*d+1:offset+j*d]*block.metric)
                for i in 1:m,j in 1:m])
            offset += d*m
        end
    end
    exact = verify_continuous_certificate(solved.reduction,-1,exact_grams)
    @test exact.identity_exact && exact.grams_psd_exact && exact.certified_bound == -1
    @test exact.grams[1] == known
    moment_only = continuous_moment_relaxation(polyopt(h,reg;ineq_constraints=balls),
        config,actions;compact_metric=metric,compress=false)
    averaged = verify_continuous_certificate(moment_only,-1,
        [known,fill(1//2,1,1),fill(1//2,1,1)])
    @test averaged.identity_exact && averaged.grams_psd_exact && averaged.certified_bound == -1
    @test averaged.grams[1] == known
    dual = continuous_nctssos(polyopt(h,reg;ineq_constraints=balls),config,actions;compact_metric=metric,dualize=true)
    @test dual.objective ≈ -1 atol=2e-5
    @test Float64(dual.certificate.residual_coefficient_bound) < 2e-5
    @test_throws ArgumentError continuous_moment_relaxation(polyopt(h,reg;ineq_constraints=[balls;1-w[1]]),
        config,actions;compact_metric=metric)

    rotations = supplied_vector_rotations(reg,[w[1:5],w[6:10]])
    free_h = -sum(w[i]*w[i+5]+w[i+5]*w[i] for i in 1:5)/2
    unbounded = continuous_moment_relaxation(polyopt(free_h,reg),config,rotations;
        compact_metric=Matrix{Int}(I,10,10))
    for scale in (1,100)
        coordinates = [1,scale,scale,scale]
        moments = reconstruct_continuous_moments(unbounded,coordinates)
        @test sum(c*moments[m] for (c,m) in terms(free_h)) == -5scale
        @test minimum(eigvals(Hermitian(ComplexF64.(reconstruct_continuous_matrices(unbounded,coordinates)[1])))) > -1e-10
    end
    model,_ = NCTSSoS.build_jump_model(unbounded.moment_problem)
    set_optimizer(model,solver)
    optimize!(model)
    @test termination_status(model) in (JuMP.MOI.DUAL_INFEASIBLE,JuMP.MOI.ALMOST_DUAL_INFEASIBLE)

    # This compact action needs sqrt(2) for its complex irreducible basis.
    # Exact moment reduction still works; unavailable cone splitting must reject.
    reg,(v,) = create_noncommutative_variables([("v",1:2)])
    action = SuppliedDerivation(reg,Dict(v[1].word[1]=>-2v[2],v[2].word[1]=>Polynomial(v[1])))
    pop = polyopt(v[1]^2+2v[2]^2,reg)
    @test_throws ArgumentError continuous_moment_relaxation(pop,config,[action];compact_metric=[2 0;0 1])
    partial = continuous_moment_relaxation(pop,config,[action];compact_metric=[2 0;0 1],compress=false)
    @test partial.relation_to_input == :equivalent
    @test [b.size for b in partial.moment_problem.linear.psd_blocks_lin] == [3]
    # Gram averaging must still work without rational irreducible splitting.
    # The second input is indefinite, but its averaged Gram is diag(0,1,2).
    for raw in (diagm([0,4//3,4//3]),diagm([0,-2,8]))
        certificate = verify_continuous_certificate(partial,0,[raw])
        @test certificate.identity_exact && certificate.certified_bound == 0
    end

    fixed = SuppliedDerivation(reg,Dict(t.word[1]=>zero(Polynomial(t)) for t in v))
    custom = SolverConfig(optimizer=solver,moment_basis=[one(v[1]),v[1]])
    with_equality = continuous_moment_relaxation(polyopt(v[1]^2,reg;moment_eq_constraints=[v[1]*v[2]]),
        custom,[fixed];compact_metric=[1 0;0 1])
    @test only(monomials(v[2]*v[1])) in with_equality.moment_map.words
    @test verify_continuous_certificate(with_equality,0,[[0 0;0 1]],
        zeros(Int,length(with_equality.equalities))).identity_exact

    reg,(v,) = create_noncommutative_variables([("v",1:4)])
    action = SuppliedDerivation(reg,Dict(v[i].word[1]=>p for (i,p) in enumerate(
        [Polynomial(v[3]),Polynomial(v[3]),-v[1]-v[2],zero(Polynomial(v[4]))])))
    custom = SolverConfig(optimizer=solver,moment_basis=[one(v[4]),v[4]])
    unused = continuous_nctssos(polyopt(v[4]^2,reg;eq_constraints=[v[1]+v[2],Polynomial(v[3])]),
        custom,[action];compact_metric=Matrix{Int}(I,4,4))
    @test unused.objective ≈ 0 atol=2e-5
    @test unused.moments[v[1]] == unused.moments[v[2]] == 0
end

@testset "Nonorthogonal scalar/vector complement" begin
    reg,(v,) = create_noncommutative_variables([("v",1:4)])
    # Mix a fixed letter into a rotating vector. Euclidean subtraction, instead
    # of the invariant-metric complement, fails to span this representation.
    frame = Rational{BigInt}[1 1 1//2 0;0 1 0 0;0 1//3 1 0;0 0 0 1]
    actions = SuppliedDerivation[]
    for (a,b) in ((2,3),(2,4),(3,4))
        k = zeros(Rational{BigInt},4,4)
        k[b,a] = 1
        k[a,b] = -1
        transformed = inv(frame)*k*frame
        push!(actions,SuppliedDerivation(reg,Dict(v[j].word[1]=>
            Polynomial([(transformed[i,j],v[i]) for i in 1:4]) for j in 1:4)))
    end
    physical_metric = inv(frame)*inv(frame)'
    objective = Polynomial([(physical_metric[i,j],only(monomials(v[i]*v[j]))) for i in 1:4 for j in 1:4])
    pop = polyopt(objective,reg)
    reduced = continuous_moment_relaxation(pop,SolverConfig(optimizer=nothing,order=1),actions;
        compact_metric=frame'*frame)
    @test sort([(b.dimension,b.multiplicity) for b in reduced.blocks]) == [(1,2),(3,1)]
    coordinates = zeros(Int,size(reduced.section,2))
    coordinates[1] = 1
    expected = zeros(Int,5,5)
    expected[1,1] = 1
    @test only(values(reconstruct_continuous_matrices(reduced,coordinates))) == expected
end

@testset "Order-two equivalence and combined equality certificates" begin
    solver = optimizer_with_attributes(COSMO.Optimizer,"verbose"=>false,
        "eps_abs"=>1e-8,"eps_rel"=>1e-8,"max_iter"=>50000)
    config = SolverConfig(optimizer=solver,order=2)
    reg,(x,y,z) = create_pauli_variables(1:2)
    actions = supplied_vector_rotations(reg,[[x[i],y[i],z[i]] for i in 1:2])
    h = x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
    pop = polyopt(h,reg)
    plain = cs_nctssos(pop,config;dualize=false)
    moment_only = continuous_nctssos(pop,config,actions;compact_metric=Matrix{Int}(I,6,6),compress=false)
    compressed = continuous_nctssos(pop,config,actions;compact_metric=Matrix{Int}(I,6,6))
    @test plain.objective ≈ -3 atol=2e-5
    @test moment_only.objective ≈ plain.objective atol=2e-5
    @test compressed.objective ≈ plain.objective atol=2e-5
    @test [b.size for b in moment_only.reduction.moment_problem.linear.psd_blocks_lin] == [32]
    @test [b.size for b in compressed.reduction.moment_problem.linear.psd_blocks_lin] == [2,6,1]
    @test compressed.reduction.section == moment_only.reduction.section
    @test reconstruct_continuous_matrices(compressed.reduction,compressed.coordinates) ==
        reconstruct_continuous_matrices(moment_only.reduction,compressed.coordinates)

    reg,(v,) = create_noncommutative_variables([("v",1:2)])
    x,y = v
    actions = supplied_vector_rotations(reg,[v])
    g = x^2+y^2-1
    commutator = x*y-y*x
    pop = polyopt(commutator^2,reg;eq_constraints=[g])
    solved = continuous_nctssos(pop,config,actions;compact_metric=[1 0;0 1],dualize=true)
    plain = cs_nctssos(pop,config;dualize=false)
    @test plain.objective ≈ -1 atol=2e-5
    @test solved.objective ≈ plain.objective atol=2e-5
    reduction = solved.reduction
    @test count(lift -> count(!iszero,lift)>1,reduction.equalities) == 2
    @test Float64(solved.certificate.residual_coefficient_bound) < 2e-5

    # Independent exact identity: c²+1 = (x²-y²)²+(xy+yx)²+g²-2xgx-2ygy.
    # Build its known Gram in the original rows, then express it in the returned
    # representation basis. No solver rounding or tolerance supplies this proof.
    rows = only(reduction.original.positive_rows)
    known = zeros(Rational{BigInt},length(rows),length(rows))
    for p in (x^2-y^2,x*y+y*x,g)
        coefficients = Dict(m=>Rational{BigInt}(c) for (c,m) in terms(p))
        vector = [get(coefficients,m,big(0)//1) for m in rows]
        known += vector*vector'
    end
    inverse = inv(hcat([b.columns for b in reduction.blocks]...))
    transformed = inverse*known*inverse'
    grams = []
    offset = 0
    for block in reduction.blocks
        d,m = block.dimension,block.multiplicity
        gram = [LinearAlgebra.tr(transformed[offset+(i-1)*d+1:offset+i*d,
            offset+(j-1)*d+1:offset+j*d]*block.metric) for i in 1:m,j in 1:m]
        push!(grams,block.embedded ? [real.(gram) -imag.(gram);imag.(gram) real.(gram)]/2 : gram)
        offset += d*m
    end
    exact = verify_continuous_certificate(reduction,-1,grams,[0,-2,-2])
    @test exact.identity_exact && exact.grams_psd_exact && exact.certified_bound == -1
    @test only(values(exact.grams)) == known
    lifted_equality = sum(c*p for (c,p) in zip(exact.equality_coefficients,reduction.original.equalities))
    @test lifted_equality == -2*(x*g*x+y*g*y)
    # Unequal multipliers require the averaged lift; the symmetric pair above
    # would also pass if reconstruction incorrectly used unaveraged tests.
    asymmetric = verify_continuous_certificate(reduction,-1,grams,[0,-4,0])
    @test asymmetric.identity_exact && asymmetric.grams_psd_exact && asymmetric.certified_bound == -1
    @test asymmetric.equality_coefficients == exact.equality_coefficients
    @test !verify_continuous_certificate(reduction,-1,grams,[0,-2,0]).identity_exact
end

@testset "Closed physics blocks and combined Pauli symmetry" begin
    solver = optimizer_with_attributes(COSMO.Optimizer, "verbose"=>false,
        "eps_abs"=>1e-8, "eps_rel"=>1e-8, "max_iter"=>50000)
    config = SolverConfig(optimizer=solver, order=2)
    reg, (x,y,z) = create_pauli_variables(1:3)
    actions = supplied_vector_rotations(reg, [[x[i],y[i],z[i]] for i in 1:3])
    metric = Matrix{Int}(I,9,9)
    bond(i,j) = x[i]*x[j] + y[i]*y[j] + z[i]*z[j]
    pop = polyopt(bond(1,2)+bond(2,3), reg)
    reflection = [pauli_site_permutation([3,2,1])]
    for (sites, expected) in ((2,-6), (3,-4))
        physics = PhysicsSparsity(strategy=:connected, max_region_sites=sites, max_block_size=37)
        continuous = continuous_nctssos(pop, config, actions; compact_metric=metric, physics)
        combined = continuous_nctssos(pop, config, actions;
            compact_metric=metric, physics, clifford_generators=reflection)
        sparsity = combined.reduction.sparsity
        plain = NCTSSoS.solve_sdp(NCTSSoS.moment_relax(pop, sparsity.corr_sparsity,
            sparsity.cliques_term_sparsities), solver; dualize=false)
        # Independent three-spin identity gives -4; separate singlet edges give -6.
        @test plain.objective ≈ expected atol=2e-5
        @test continuous.objective ≈ plain.objective atol=2e-5
        @test combined.objective ≈ plain.objective atol=2e-5
        @test combined.reduction.original.positive_rows == continuous.reduction.original.positive_rows
        for (idx, mat) in combined.matrices
            direct = map(p -> sum(c*combined.moments[m] for (c,m) in terms(p); init=0),
                combined.reduction.original.constraints[idx][2])
            @test mat ≈ direct atol=1e-12
        end
        for (m, value) in combined.moments
            sign, image = NCTSSoS._act_monomial(only(reflection), m)
            @test value ≈ sign*combined.moments[image] atol=1e-12
        end
    end

    # A non-tight control prevents accidentally restricting to the ground state
    # or to a total-spin sector: degree one gives -4.5, whereas the energy is -3.
    triangle = polyopt(bond(1,2)+bond(2,3)+bond(1,3), reg)
    order1 = SolverConfig(optimizer=solver, order=1)
    cycle = [pauli_site_permutation([2,3,1])]
    zero_action = [SuppliedDerivation(reg, Dict(i=>zero(triangle.objective) for i in indices(reg)))]
    plain = cs_nctssos(triangle, order1; dualize=false)
    @test plain.objective ≈ -4.5 atol=2e-5
    for generators in (zero_action, actions)
        combined = continuous_nctssos(triangle, order1, generators;
            compact_metric=metric, clifford_generators=cycle)
        @test combined.objective ≈ plain.objective atol=2e-5
    end

    # Spatial averaging must move multipliers between different edge blocks.
    # Three times one bond certificate averages to H+9, exactly, not H+3.
    physics = PhysicsSparsity(strategy=:connected, max_region_sites=2, max_block_size=16)
    reduction = continuous_moment_relaxation(triangle, config, actions;
        compact_metric=metric, physics, clifford_generators=cycle, compress=false)
    grams = []
    selected = false
    for block in reduction.blocks
        rows = reduction.original.positive_rows[findfirst(==(block.original_index), reduction.original.positive_indices)]
        q = zeros(Rational{BigInt}, length(rows), length(rows))
        if !selected && all(m in rows for m in (x[1],x[2],y[1],y[2],z[1],z[2]))
            for axis in (x,y,z)
                v = [m == axis[1] || m == axis[2] ? 1 : 0 for m in rows]
                q += (big(3)//2)*(v*v')
            end
            selected = true
        end
        push!(grams, block.embedded ? [q zero(q); zero(q) q]/2 : q)
    end
    @test selected
    certificate = verify_continuous_certificate(reduction, -9, grams)
    @test certificate.identity_exact && certificate.grams_psd_exact
    @test certificate.certified_bound == -9
    @test all(any(!iszero, q) for q in values(certificate.grams))

    # A global axis rotation normalizes (but does not commute with) all three
    # spin generators. A single-axis subgroup is not preserved by that rotation.
    axes = pauli_global_axis_rotation_generators((x,y,z))
    axis_reduction = continuous_moment_relaxation(triangle, order1, actions;
        compact_metric=metric, clifford_generators=axes)
    @test axis_reduction.section == continuous_moment_relaxation(triangle, order1, actions;
        compact_metric=metric).section
    @test_throws ArgumentError continuous_moment_relaxation(triangle, order1, actions[1:1];
        compact_metric=metric, clifford_generators=axes)
    asymmetric = polyopt(bond(1,2)+2bond(2,3), reg)
    @test_throws ArgumentError continuous_moment_relaxation(asymmetric, order1, actions;
        compact_metric=metric, clifford_generators=reflection)
    @test_throws ArgumentError continuous_moment_relaxation(polyopt(triangle.objective,reg;
        eq_constraints=[bond(1,2)-1]), config, actions; compact_metric=metric, clifford_generators=cycle)

    # Signed transport between blocks is essential: averaging (1+y₁)²
    # with (1-y₂)² gives H+2 for H=y₁-y₂, not 2+y₁+y₂.
    signed_reg, (sx,sy,sz) = create_pauli_variables(1:2)
    signed_pop = polyopt(sy[1]-sy[2], signed_reg)
    around_y = supplied_vector_rotations(signed_reg, [[sx[i],sy[i],sz[i]] for i in 1:2])[2:2]
    signed_swap = CliffordSymmetry(Dict(
        sx[1]=>(1,sz[2]), sz[1]=>(1,sx[2]), sy[1]=>(-1,sy[2]),
        sx[2]=>(1,sz[1]), sz[2]=>(1,sx[1]), sy[2]=>(-1,sy[1])); nqubits=2)
    signed_reduction = continuous_moment_relaxation(signed_pop, order1, around_y;
        compact_metric=Matrix{Int}(I,6,6),
        physics=PhysicsSparsity(strategy=:connected,max_region_sites=1,max_block_size=4),
        clifford_generators=[signed_swap], compress=false)
    signed_grams = []
    for block in signed_reduction.blocks
        pos = findfirst(==(block.original_index),signed_reduction.original.positive_indices)
        rows = signed_reduction.original.positive_rows[pos]
        v = [sy[1] in rows && (isone(m) || m == sy[1]) ? 1 : 0 for m in rows]
        q = v*v'
        push!(signed_grams, block.embedded ? [q zero(q);zero(q) q]/2 : q)
    end
    signed_certificate = verify_continuous_certificate(signed_reduction,-2,signed_grams)
    @test signed_certificate.identity_exact && signed_certificate.grams_psd_exact
    @test signed_certificate.certified_bound == -2

    # Degree is not an ignored request, and closure is charged to the row cap.
    r4, ops4 = create_pauli_variables(1:4)
    p4 = polyopt(heisenberg_chain_hamiltonian(ops4; periodic=false), r4)
    policy = PhysicsSparsity(strategy=:connected, max_region_sites=4, max_block_size=256)
    b2 = only(compute_sparsity(p4, SolverConfig(optimizer=nothing,order=2), policy).corr_sparsity.clq_mom_mtx_bases)
    b4 = only(compute_sparsity(p4, SolverConfig(optimizer=nothing,order=4), policy).corr_sparsity.clq_mom_mtx_bases)
    @test length(b2) == 67
    @test length(b4) == 256
    @test issubset(b2, b4)
    @test maximum(degree, b4) == 4
    @test_throws ArgumentError compute_sparsity(p4, SolverConfig(optimizer=nothing,order=4),
        PhysicsSparsity(strategy=:connected,max_region_sites=4,max_block_size=67))
    @test_throws ArgumentError compute_sparsity(pop, order1,
        PhysicsSparsity(max_region_sites=1,max_block_size=3); actions)
end
