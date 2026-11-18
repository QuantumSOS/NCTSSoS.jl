# Run with Clarabel available, as for test_workflow.jl.
module PauliStationarityChecks

import NCTSSoS as NC
using Clarabel
using Test
include("PauliWorkflowRedesign.jl")
using .PauliWorkflowRedesign
include("benchmark/cases.jl")

@testset "Stationarity on the chain-8 RDM stencils" begin
    n = 8
    _, (x,y,z) = NC.pauli_variables(n)
    H = (1/(4n))*sum(x[i]*x[j] + y[i]*y[j] + z[i]*z[j] for (i,j) in bonds("chain",n))
    regions = [sort([mod1(i+j,n) for j in 0:3]) for i in 1:n]
    plan = compile_relaxation(H; order=3,moi_form=:auto,
        symmetry=PauliSymmetrySpec(spatial_mode=:circulant_dft),
        sparsity=PauliSparsitySpec(support_mode=:contiguous_lattice,two_body_range_r=4,
            cs_algo=:NoElimination,ts_algo=:Maximal,ts_step=0),
        strengthening=PauliStrengtheningSpec(commutator_linear_order=3,rdm_positivity_k=4),
        rdm_regions=regions)
    @test plan.reports.symmetry.route === :lattice
    @test plan.reports.support.row_counts == [493] # no enlargement of the moment basis
    @test plan.reports.analysis.rdm.regions == n

    # Many emitted rows have the same content after symmetry. Count exact
    # independent constraints on the compiler's reduced moment coordinates.
    forms = plan.program.forms
    Q = Rational{BigInt}
    pivots = Dict{Int,Dict{Int,Q}}()
    for b in plan.program.blocks
        b.kind === NC.ZERO && b.label > 0 || continue
        for f in b.first_form:b.first_form+b.side-1
            row = Dict{Int,Q}()
            for t in NC.term_range(forms,f)
                id = Int(forms.ids[t])
                row[id] = get(row,id,zero(Q)) + Q(forms.coefs[t])
            end
            filter!(pair -> !iszero(last(pair)),row)
            NC._insert_generator_equation!(pivots,row;constant=false)
        end
    end
    @test length(pivots) == 1

    # Use the SOS form of this very same plan for the local Clarabel control.
    # The benchmark still uses :auto with Mosek; the default moment form can
    # stall at reduced accuracy with Clarabel on this rank-deficient model.
    model, _ = NC.to_jump(NC.sos_standard_form(plan.standard))
    NC.JuMP.set_optimizer(model,Clarabel.Optimizer)
    NC.JuMP.set_silent(model)
    NC.JuMP.set_optimizer_attribute(model,"chordal_decomposition_enable",false)
    NC.JuMP.optimize!(model)
    @test NC.JuMP.termination_status(model) == NC.JuMP.MOI.OPTIMAL
    @test NC.JuMP.primal_status(model) == NC.JuMP.MOI.FEASIBLE_POINT
    @test NC.JuMP.dual_status(model) == NC.JuMP.MOI.FEASIBLE_POINT
    @test NC.JuMP.objective_value(model) ≈ exact_energy("chain",n) atol=1e-6
end

end
