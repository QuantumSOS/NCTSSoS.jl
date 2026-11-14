# Tests of the conic layer of NCTSSoS.Next. Include the test files of this layer here.
using Test, NCTSSoS, JuMP, LinearAlgebra
isdefined(@__MODULE__,:Nx) || include("../contract_checks.jl")
isdefined(@__MODULE__,:dense_program) || include("../reference.jl")
isdefined(@__MODULE__,:SOLVER) || include("../../TestUtils.jl")

function conic_bound(p)
    model,x = Nx.to_jump(p)
    set_optimizer(model,SOLVER)
    optimize!(model)
    @test termination_status(model) in (JuMP.MOI.OPTIMAL,JuMP.MOI.ALMOST_OPTIMAL)
    return objective_value(model)
end
function check_forms(p;expected=nothing,tolerance=2e-5)
    value = conic_bound(p)
    moment = conic_bound(Nx.moment_standard_form(p))
    sos = conic_bound(Nx.sos_standard_form(p))
    @test moment ≈ value atol=tolerance
    @test sos ≈ value atol=tolerance
    isnothing(expected) || @test value ≈ expected atol=tolerance
    @info "Standard-form comparison" reference=value moment sos expected
    return value
end
function add_form!(p,terms)
    for (id,c) in terms; Nx.push_term!(p.forms,c,id); end
    Nx.close_form!(p.forms)
end

@testset "Standard forms: all cones, constants and complex entries" begin
    p = Nx.MomentProgram{ComplexF64}(6)
    # Variables x,y,z,u,t. The optimum is x=2,y=3,z=13/4,u=9/4,t=√13.
    add_form!(p,[(2,1),(1,-2),(3,im),(1,-3im)])
    Nx.push_block!(p,Nx.ZERO,1)
    add_form!(p,[(4,1)])
    Nx.push_block!(p,Nx.NONNEGATIVE,1)
    for id in (6,2,3); add_form!(p,[(id,1)]); end
    Nx.push_block!(p,Nx.SOC,3)
    for terms in ([(5,1)],[(3,1)],[(1,4)]); add_form!(p,terms); end
    Nx.push_block!(p,Nx.PSD,2)
    for terms in ([(4,1)],[(2,1),(3,im)],[(1,4)]); add_form!(p,terms); end
    Nx.push_block!(p,Nx.HPSD,2,77)
    append!(p.objective_ids,Int32[1,4,5,6]); append!(p.objective_coefs,[5.,.7,.3,1.])
    standard = Nx.moment_standard_form(p)
    @test standard.b-standard.A*[2.,3.,13/4,9/4,sqrt(13)] ≈
        [0,0,13/4,sqrt(13),2,3,9/4,3sqrt(2),4,13/4,2sqrt(2),3sqrt(2),4]
    dual = Nx.sos_standard_form(standard)
    @test dual.A==transpose(standard.A)
    @test dual.b == -standard.c && dual.c == -standard.b
    @test first(dual.cones).kind==:free
    @test last(standard.cones).label==last(dual.cones).label==77
    check_forms(p;expected=5+.7*13/4+.3*9/4+sqrt(13))
    W = Nx.PauliWord{1}
    constant,_ = Nx.plain_program([one(W)],[(2.5,one(W))])
    @test size(Nx.moment_standard_form(constant).A)==(1,0)
    check_forms(constant;expected=2.5)
end

@testset "Plain route and both forms equal cs_nctssos" begin
    reg,(x,y,z) = create_pauli_variables(1:4)
    ham = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3)
    old = cs_nctssos(polyopt(ham,reg),SolverConfig(optimizer=SOLVER,order=1);dualize=false).objective
    W = Nx.PauliWord{2}
    rows = Nx.basis_rows(W,Nx.Pauli(),63:66,1)
    objective = [(1.,Nx.mul(Nx.pauli_letter(W,i,a),Nx.pauli_letter(W,i+1,a))[2]) for i in 63:65 for a in 0:2]
    p,_ = Nx.plain_program(rows,objective)
    ref,_ = dense_program(rows,objective)
    @test conic_bound(ref) ≈ old atol=2e-5
    check_forms(p;expected=old)

    reg,(a,ad) = create_fermionic_variables(1:3)
    ham = -(ad[1]*a[2]+ad[2]*a[1])-.7*(ad[2]*a[3]+ad[3]*a[2])+1.3*ad[1]*a[1]*ad[2]*a[2]+.4*ad[3]*a[3]
    old = cs_nctssos(polyopt(ham,reg),SolverConfig(optimizer=SOLVER,order=2);dualize=false).objective
    W = Nx.FermiWord{1}
    rows = Nx.basis_rows(W,Nx.Fermionic(),1:3,2)
    @test length(rows)==22 # 1 + 6 letters + choose(6,2) normal words
    objective = [(c,Nx.word_from_legacy(W,m.word)) for (c,m) in terms(ham)]
    p,_ = Nx.plain_program(rows,objective)
    ref,_ = dense_program(rows,objective)
    @test conic_bound(ref) ≈ old atol=2e-4
    check_forms(p;expected=old,tolerance=2e-4)

    reg,(x,y) = create_unipotent_variables([("x",1:2),("y",1:2)])
    ham = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    old = cs_nctssos(polyopt(ham,reg),SolverConfig(optimizer=SOLVER,order=1);dualize=false).objective
    algebra = Nx.SiteAlgebra(Nx.Unipotent(),2)
    W = Nx.word_type(algebra,2,2)
    rows = Nx.basis_rows(W,algebra,1:2,1)
    objective = [(-s,Nx.mul(Nx.site_letter(W,1,i),Nx.site_letter(W,2,j))[2]) for (i,j,s) in ((1,1,1.),(1,2,1.),(2,1,1.),(2,2,-1.))]
    real_program,_ = Nx.plain_program(rows,objective;C=Float64)
    @test only(real_program.blocks).kind==Nx.PSD
    @test real_program.nvars==11
    @test conic_bound(real_program) ≈ old atol=2e-5
    @test_throws r"real objective coefficients" Nx.plain_program(rows,[(im,rows[2])];C=Float64)
    p,_ = Nx.plain_program(rows,objective)
    check_forms(p;expected=old)
    @test old ≈ -2sqrt(2) atol=2e-5
end
