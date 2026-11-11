# Small geometry-only controls, not a paper-scale or fair timing benchmark.
# Run with the local package, Graphs, JuMP, and licensed MosekTools available.
using NCTSSoS, Graphs, JuMP, MosekTools, LinearAlgebra, Test
import NCTSSoS: degree
const MOI = JuMP.MOI

function main()
    solver = optimizer_with_attributes(MosekTools.Optimizer,
        "MSK_IPAR_NUM_THREADS"=>1, "MSK_IPAR_LOG"=>0,
        "MSK_DPAR_OPTIMIZER_MAX_TIME"=>120.0,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>1e-8,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>1e-8,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>1e-8)
    license = Model(solver)
    set_silent(license)
    @variable(license, v >= 1)
    @objective(license, Min, v)
    optimize!(license)
    @assert termination_status(license) == MOI.OPTIMAL

    println("case,j2,variant,rows,lower,exact,relative_gap,reduced_violation,termination,primal_status,dual_status,accepted,seconds")
    paulis = ([0 1;1 0], [0 -im;im 0], [1 0;0 -1])
    for (label, n) in [("ring6",6), ("square_plaquette4",4)]
        reg, ops = create_pauli_variables(1:n)
        geometry = cycle_graph(n)
        nn = [(src(e),dst(e)) for e in edges(geometry)]
        # Four sites are a single open square plaquette, NOT the paper's
        # larger periodic square. Its next neighbors are the two diagonals.
        nnn = unique([minmax(i,mod1(i+2,n)) for i in 1:n])
        nn_h = sum(sum(axis[i]*axis[j] for axis in ops) for (i,j) in nn)/4
        selected = physics_moment_basis(reg, geometry; radius=1, max_rows=1000)
        paper = physics_moment_basis(polyopt(nn_h,reg),4; pair_distance=1)
        spatial = SymmetrySpec([pauli_site_permutation([2:n;1]),
            pauli_site_permutation(reverse(1:n))]; offblock_check=:full)
        for j2 in (0.0,0.3)
            h = nn_h + j2*sum(sum(axis[i]*axis[j] for axis in ops) for (i,j) in nnn)/4
            pop = polyopt(h,reg)
            matrix = zeros(ComplexF64,2^n,2^n)
            for (pairs, coupling) in ((nn,1.0),(nnn,j2)), (i,j) in pairs, p in paulis
                matrix += coupling/4 * reduce(kron,
                    [k in (i,j) ? p : Matrix{ComplexF64}(I,2,2) for k in 1:n])
            end
            exact = eigmin(Hermitian(matrix))
            variants = label == "ring6" ? (:geometry,:paper_family,:occurrence) :
                (:geometry,:geometry_plain,:paper_family,:occurrence)
            geometry_value = NaN
            for variant in variants
                started = time()
                basis = variant == :paper_family ? paper : selected
                if label == "ring6" && variant != :occurrence
                    mp, _ = pauli_translation_invariant_moment_relaxation(pop,ops,maximum(degree,basis);
                        basis,axis_permutation_symmetry=true)
                    model = NCTSSoS.sos_dualize(mp).model
                    set_optimizer(model,solver)
                    set_silent(model)
                    optimize!(model)
                else
                    config = variant == :occurrence ?
                        SolverConfig(optimizer=solver,order=3,cs_algo=MF(),ts_algo=MMD()) :
                        SolverConfig(optimizer=solver,moment_basis=basis,
                            symmetry=variant == :geometry_plain ? nothing : spatial)
                    result = cs_nctssos(pop,config)
                    model = result.model
                end
                primal, dual = objective_value(model), dual_objective_value(model)
                lower = objective_sense(model) == MOI.MAX_SENSE ? primal : dual
                gap = abs(primal-dual)/max(1.0,abs(primal),abs(dual))
                violation = maximum(values(primal_feasibility_report(model;atol=0.0));init=0.0)
                accepted = termination_status(model) == MOI.OPTIMAL &&
                    primal_status(model) == MOI.FEASIBLE_POINT && dual_status(model) == MOI.FEASIBLE_POINT &&
                    gap <= 1e-7 && violation <= 1e-7 && lower <= exact+1e-7
                println(join((label,j2,variant,variant == :occurrence ? "generated" : length(basis),
                    lower,exact,gap,violation,termination_status(model),primal_status(model),
                    dual_status(model),accepted,time()-started),','))
                flush(stdout)
                @test accepted
                if variant == :geometry
                    geometry_value = lower
                elseif variant == :geometry_plain
                    @test lower ≈ geometry_value atol=1e-7
                end
            end
        end
    end
end

@testset "Geometry-only controls" begin
    main()
end
