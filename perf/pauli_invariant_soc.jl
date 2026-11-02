# Run one isolated case/method, e.g.:
# julia --project=<environment with NCTSSoS, MosekTools, JuMP, Graphs>
#     perf/pauli_invariant_soc.jl rich3 combined result.toml
# The report records cold cost separately, then three warmed independent builds.
using NCTSSoS, JuMP, MosekTools, LinearAlgebra, Graphs, TOML, Statistics

const SOC_OPTIMIZER = optimizer_with_attributes(MosekTools.Optimizer,
    "MSK_IPAR_NUM_THREADS" => 1,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-9,
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-9,
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-9)

function soc_case(name)
    k = name == "rich4" ? 4 : 3
    reg, (x, y, z) = create_pauli_variables(1:k)
    if name == "supplied"
        h = -(x[1]*y[2] + x[1]*y[3] + x[2]*x[3] + y[2]*y[3] + z[2]*z[3])
    elseif name in ("rich3", "rich4")
        # First draw: Random.seed!(382); rand(-5:5, 27) ./ 5, Julia 1.12.7.
        weights = [2,-5,5,-5,4,-4,-5,-4,-3,3,-3,-1,-3,-5,-2,-3,-4,-4,4,2,4,0,4,2,-5,4,-2] ./ 5
        reg3, _ = create_pauli_variables(1:3)
        pairs = filter(m -> NCTSSoS.degree(m) == 2, get_ncbasis(reg3, 2))
        h = sum(c*m for (c, m) in zip(weights, pairs))
        name == "rich4" && (h += -0.7*(x[1]*y[4] + x[2]*x[4] + y[2]*y[4] + z[2]*z[4]) + 0.23z[1] + 0.17x[4])
    else
        error("Unknown case $name")
    end
    return polyopt(h, reg), k
end

function soc_label(mono, k)
    letters = fill('I', k)
    for letter in mono.word
        letters[(letter-1) ÷ 3 + 1] = "XYZ"[(letter-1) % 3 + 1]
    end
    return String(letters)
end

function soc_exact_energy(h, k)
    matrices = Dict('I' => ComplexF64[1 0; 0 1], 'X' => ComplexF64[0 1; 1 0],
        'Y' => ComplexF64[0 -im; im 0], 'Z' => ComplexF64[1 0; 0 -1])
    H = zeros(ComplexF64, 2^k, 2^k)
    for (c, m) in terms(h)
        H += c * reduce(kron, [matrices[letter] for letter in soc_label(m, k)])
    end
    return eigmin(Hermitian(H))
end

function soc_anticommutation!(model, mp)
    monos = [mp.linear.key_to_monomial[key] for key in mp.linear.moments]
    graph = SimpleGraph(length(monos))
    for i in eachindex(monos), j in 1:(i-1)
        NCTSSoS._pauli_anticommutes(monos[i], monos[j]) && add_edge!(graph, i, j)
    end
    cliques = filter(c -> length(c) > 1, maximal_cliques(graph))
    for clique in cliques
        @constraint(model, [1.0; model[:y_re][clique]] in SecondOrderCone())
    end
    return cliques
end

function soc_trial(pop, method)
    method in ("bare", "anti", "soc", "combined", "one_soc", "soc4", "extra1", "extra3", "extra5") || error("Unknown method $method")
    extra = startswith(method, "extra") ? parse(Int, method[6:end]) : 0
    cliques = Vector{Int}[]
    setup = @timed begin
        basis = get_ncbasis(pop.registry, 1)
        extra > 0 && append!(basis, filter(m -> NCTSSoS.degree(m) == 2, get_ncbasis(pop.registry, 2))[1:extra])
        config = SolverConfig(optimizer=SOC_OPTIMIZER, moment_basis=basis)
        sparsity = compute_sparsity(pop, config)
        mp = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
        model, _ = build_jump_model(mp)
        set_optimizer(model, SOC_OPTIMIZER)
        set_silent(model)
        method in ("anti", "combined") && (cliques = soc_anticommutation!(model, mp))
        clusters = NCTSSoS._soc_clusters(pop, sparsity, 64)
        method == "soc4" && filter!(cluster -> length(cluster) == 4, clusters)
    end
    rounds = method in ("soc", "combined", "soc4") ? 3 : method == "one_soc" ? 1 : 0
    result = separate_pauli_invariant_socs!(model, mp, clusters;
        max_rounds=rounds, batch_size=method == "one_soc" ? 1 : 2)
    row = Dict("construction_s" => setup.time, "construction_bytes" => setup.bytes,
        "separation_s" => result.separation_seconds, "solve_s" => result.solve_seconds,
        "primal" => result.objective, "dual" => result.bound,
        "status" => string(termination_status(model)), "variables" => num_variables(model),
        "moments" => length(mp.linear.moments), "basis_words" => length(basis),
        "anticommutation_cones" => length(cliques), "invariant_cones" => length(result.records),
        "objectives" => result.objectives)
    return (; row, result, mp, basis, cliques)
end

# Exact-arithmetic postprocessing of numerical duals into a quantum lower bound.
# Every physical Pauli expectation is in [-1,1], the identity is 1, and imaginary
# parts vanish. PSD multipliers B*B' are exactly positive semidefinite; SOC
# multipliers satisfy t² >= ||v||² exactly. Residual coefficients cost their l1
# norm. This certifies a lower bound even if solver dual feasibility is imperfect.
soc_rat(x) = rationalize(BigInt, x; tol=0)
soc_small_rat(x) = round(BigInt, x * 10^9) // big(10)^9

function soc_certificate(model, mp)
    residual = Dict(v => big(0)//big(1) for v in all_variables(model))
    constant = big(0)//big(1)
    function accumulate!(expr, weight)
        if expr isa Number
            constant += weight * soc_rat(expr)
        else
            constant += weight * soc_rat(JuMP.constant(expr))
            for (coef, variable) in linear_terms(expr)
                residual[variable] += weight * soc_rat(coef)
            end
        end
    end
    accumulate!(objective_function(model), big(1)//big(1))
    for cref in all_constraints(model; include_variable_in_set_constraints=false)
        object = constraint_object(cref)
        f, set = object.func, object.set
        d = dual(cref)
        if set isa JuMP.MOI.EqualTo
            q = soc_small_rat(d)
            accumulate!(f, -q)
            constant += q * soc_rat(set.value)
        elseif set isa JuMP.MOI.SecondOrderCone
            q = soc_small_rat.(d)
            q[1] = ceil(BigInt, norm(Float64.(q[2:end])) * 10^9 + 1) // big(10)^9
            @assert q[1] >= 0 && q[1]^2 >= sum(abs2, q[2:end])
            # Replace irrational radii sqrt(3), sqrt(7) by rational upper
            # bounds, checked by squaring, hence use a slightly weaker cone.
            radius = JuMP.constant(f[1])
            radius2 = round(Int, radius^2)
            upper = ceil(BigInt, radius * 10^9 + 1) // big(10)^9
            @assert upper^2 >= radius2
            for i in eachindex(f)
                accumulate!(f[i], -q[i])
            end
            constant -= q[1] * (upper - soc_rat(radius))
        elseif set isa JuMP.MOI.PositiveSemidefiniteConeSquare
            eigen = LinearAlgebra.eigen(Symmetric((Matrix(d) + transpose(Matrix(d))) / 2))
            B = soc_small_rat.(eigen.vectors * Diagonal(sqrt.(max.(0, eigen.values))))
            Q = B * transpose(B)
            F = reshape(f, set.side_dimension, set.side_dimension)
            for i in axes(F, 1), j in axes(F, 2)
                accumulate!(F[i,j], -Q[i,j])
            end
        else
            error("Unsupported certificate cone $(typeof(set))")
        end
    end
    y = model[:y_re]
    identity = mp.linear.moment_index[mp.linear.identity]
    residual_l1 = sum(abs(residual[y[i]]) for i in eachindex(y) if i != identity)
    lower = constant + residual[y[identity]] - residual_l1
    return Dict("lower" => Float64(lower), "lower_rational" => string(lower),
        "affine_constant" => string(constant), "identity_index" => identity,
        "real_residuals" => string.([residual[v] for v in y]),
        "residual_l1" => Float64(residual_l1), "residual_l1_rational" => string(residual_l1))
end

function soc_benchmark(case, method, output)
    BLAS.set_num_threads(1)
    pop, k = soc_case(case)
    cold = @elapsed soc_trial(pop, method)
    trials = [soc_trial(pop, method) for _ in 1:3]
    last = trials[end]
    certificate = soc_certificate(last.result.model, last.mp)
    exact = soc_exact_energy(pop.objective, k)
    @assert all(t -> abs(t.row["primal"] - t.row["dual"]) < 1e-6, trials)
    @assert certificate["lower"] <= last.result.objective + 1e-8
    @assert certificate["lower"] > last.result.bound - 1e-4 (certificate=certificate["lower"], dual=last.result.bound)
    @assert certificate["lower"] <= exact + 1e-8
    data = Dict("case" => case, "method" => method, "julia" => string(VERSION),
        "solver" => solver_name(last.result.model), "cpu" => Sys.cpu_info()[1].model,
        "solver_version" => JuMP.MOI.get(backend(last.result.model), JuMP.MOI.SolverVersion()),
        "mosektools_version" => string(pkgversion(MosekTools)),
        "threads" => 1, "tolerance" => 1e-9, "cold_s" => cold,
        "peak_rss_bytes" => Sys.maxrss(), "exact_energy" => exact,
        "certified_lower" => certificate["lower"], "certificate_rational" => certificate["lower_rational"],
        "certificate_details" => certificate,
        "hamiltonian" => Dict(soc_label(m,k) => real(c) for (c,m) in terms(pop.objective)),
        "basis" => soc_label.(last.basis, k), "samples" => [t.row for t in trials],
        "moments" => [soc_label(last.mp.linear.key_to_monomial[key],k) for key in last.mp.linear.moments],
        "anticommutation_cliques" => last.cliques,
        "constraints" => [Dict("cluster"=>r.cluster, "r"=>r.r, "labels"=>r.labels,
            "constant"=>r.constant, "constant_squared"=>r.constant_squared,
            "round"=>r.round, "violation"=>r.violation) for r in last.result.records])
    open(output, "w") do io
        TOML.print(io, data; sorted=true)
    end
    println((case=case, method=method, certified_lower=certificate["lower"], exact=exact,
        construction_s=median(t.row["construction_s"] for t in trials),
        solve_s=median(t.row["solve_s"] for t in trials), peak_rss_bytes=Sys.maxrss()))
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 3 || error("Expected CASE METHOD OUTPUT.toml")
    soc_benchmark(ARGS...)
end
