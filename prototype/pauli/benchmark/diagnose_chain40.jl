# From the reviewed checkout (using its existing CVM benchmark environment):
# JULIA_LOAD_PATH="@:$PWD:@stdlib" julia --startup-file=no --threads=1,0 \
#   --project="$HOME/pauli-bench-env" prototype/pauli/benchmark/diagnose_chain40.jl [raw|compiled]
using LinearAlgebra, SparseArrays, JuMP, MosekTools
import Mosek
include("cases.jl")
include("prototype_adapter.jl")

BLAS.set_num_threads(1)
Threads.nthreads(:default) == 1 && Threads.nthreads(:interactive) == 0 || error("use --threads=1,0")
mode = isempty(ARGS) ? "compiled" : only(ARGS)
mode in ("raw","compiled") || error("expected raw or compiled")

const DIAGNOSTIC_TOL = 1e-8
const DIAGNOSTIC_PARAMETERS = (
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => DIAGNOSTIC_TOL,
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => DIAGNOSTIC_TOL,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => DIAGNOSTIC_TOL,
    "MSK_IPAR_NUM_THREADS" => 1)

function diagnose(form, name)
    model,_ = NCTSSoS.to_jump(form)
    set_optimizer(model, optimizer_with_attributes(Mosek.Optimizer, DIAGNOSTIC_PARAMETERS...))
    for (key,value) in DIAGNOSTIC_PARAMETERS
        get_optimizer_attribute(model,key) == value || error("unexpected $key")
    end
    println((;name,A_size=size(form.A),A_nnz=nnz(form.A)))
    flush(stdout)
    optimize!(model)
    term, primal, dual = termination_status(model), primal_status(model), dual_status(model)
    accepted_status = term == MOI.OPTIMAL && primal == dual == MOI.FEASIBLE_POINT
    p,d = accepted_status ? (objective_value(model),dual_objective_value(model)) : (NaN,NaN)
    gap = abs(p-d)/max(1,abs(p),abs(d))
    accepted = accepted_status && isfinite(p) && isfinite(d) && gap <= DIAGNOSTIC_TOL
    println((;name,term,primal,dual,raw=raw_status(model),accepted,p,d,gap,solver_seconds=solve_time(model)))
    return accepted
end

L = 40
_,(x,y,z) = NCTSSoS.pauli_variables(L)
H = (1/(4L))*sum(x[i]*x[j]+y[i]*y[j]+z[i]*z[j] for (i,j) in bonds("chain",L))
plan = compile_relaxation(H; order=2,
    symmetry=PauliSymmetrySpec(spatial_mode=:circulant_dft),
    sparsity=PauliSparsitySpec(support_mode=:contiguous_lattice,two_body_range_r=4,
        cs_algo=:NoElimination,ts_algo=:Maximal,ts_step=0))
@assert plan.reports.symmetry.route === :lattice && plan.reports.term_sparsity.zero_lift_audit
moment = mode == "raw" ? NCTSSoS.moment_standard_form(plan.program) : plan.standard
println((;mode,auto_form=plan.moi_form,mosek=Mosek.getversion(),julia=VERSION,
    empty_moment_columns=count(j->all(iszero,view(moment.A.nzval,nzrange(moment.A,j))),axes(moment.A,2))))
accepted_moment = diagnose(moment,:dualized)
accepted_sos = diagnose(NCTSSoS.sos_standard_form(moment),:primalized)
exit(accepted_moment && accepted_sos ? 0 : 1)
