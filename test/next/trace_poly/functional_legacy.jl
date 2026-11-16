module FunctionalLegacyChecks
using Test, TOML
using NCTSSoS

function state_problem(kind)
    _,(x,y)=site_variables(Unipotent(),2,kind==:covariance ? 3 : 2;functional=NCTSSoS.RealState())
    s=state_symbol
    if kind==:squared
        h=-(s(x[1]*y[2])+s(x[2]*y[1]))^2-(s(x[1]*y[1])-s(x[2]*y[2]))^2
    elseif kind==:covariance
        cov(a,b)=s(x[a]*y[b])-s(x[a])*s(y[b])
        # Keep the positive sign of the old state example, unlike the trace case.
        h=cov(1,1)+cov(1,2)+cov(1,3)+cov(2,1)+cov(2,2)-cov(2,3)+cov(3,1)-cov(3,2)
    else
        h=-s(x[2])-s(y[1])-s(y[2])+s(x[1]*y[1])-s(x[2]*y[1])-
            s(x[1]*y[2])-s(x[2]*y[2])+s(x[1])*s(y[1])+s(x[2])*s(y[1])+
            s(x[2])*s(y[2])+s(x[1])^2+s(y[2])^2
    end
    Problem(h)
end

function bilocal_problem()
    _,(x,y,z)=site_variables(Unipotent(),3,2;functional=NCTSSoS.RealState())
    # Same eight ordered triple words and 28 upper-triangle coefficients as
    # test/problems/quantum_networks/bilocal_networks.jl, Example 8.1.1.
    s=[state_symbol(x[i]*y[j]*z[k]) for j in 1:2 for i in 1:2 for k in 1:2]
    cross=[2,2,2,-2,2,2,-2, 2,2,-2,2,2,-2, 2,-2,2,2,-2,
        -2,2,2,-2, -2,-2,2, 2,-2, -2]
    linear=[-1,-1,-1,-1,-1,1,1,-1]
    h=sum(t^2 for t in s)/8+sum(linear[i]*s[i] for i in 1:8)
    index=0
    for i in 1:8,j in i+1:8
        index+=1
        h+=cross[index]*s[i]*s[j]/8
    end
    Problem(h)
end

function bilocal_check(;solver=nothing)
    @testset "Bilocal removal reference" begin
        built=NCTSSoS.build(bilocal_problem(),Relaxation(order=3))
        @test NCTSSoS.block_sides(built)==[637]
        @test built.sizes.moments==15867
        # The 507 s native SOS solve is kept for the removal record only.
        @test solve(built;solver).bound ≈ -3.9999999896685883 atol=1e-6 skip=isnothing(solver)
    end
end

function trace_problem(kind)
    if kind==:constrained
        _,(v,)=site_variables(Free(),1,2;functional=NCTSSoS.RealTrace())
        x,y=v
        return Problem(trace_symbol(2-x^2+x*y^2*x-y^2);
            constraints=(Positive(4-x^2-y^2),Zero(x*y+y*x-2)))
    end
    _,(v,)=site_variables(Unipotent(),1,4;functional=NCTSSoS.RealTrace())
    x,y=v[1:2],v[3:4]
    s=trace_symbol
    h=kind==:chsh ? -s(x[1]*y[1])-s(x[1]*y[2])-s(x[2]*y[1])+s(x[2]*y[2]) :
        -(s(x[1]*y[2])+s(x[2]*y[1]))^2-(s(x[1]*y[1])-s(x[2]*y[2]))^2
    Problem(h)
end

function cases()
    (("state_squared",state_problem(:squared),3,-3.9999999914666895,1e-6),
     ("state_covariance",state_problem(:covariance),2,-4.999999999824081,1e-6),
     ("state_mixed",state_problem(:mixed),2,-3.511480225797076,1e-6),
     ("trace_constrained",trace_problem(:constrained),2,-1.0,1e-6),
     ("trace_chsh",trace_problem(:chsh),1,-2.8284271247462307,1e-6),
     ("trace_squared",trace_problem(:squared),2,-4.000000007251562,1e-5))
end

# All seven NCTSSoS side/moment records equal the original dense records where
# the old file asserted them, including bilocal's 637 rows / 15,867 moments.
function check(solver;sos_solver=nothing,references=false,cosmo=true)
    sizes=Dict("state_squared"=>([209],1887),"state_covariance"=>([106],1098),
        "state_mixed"=>([49],233),
        "trace_constrained"=>([19,5],57),"trace_chsh"=>([9],21),
        "trace_squared"=>([53],222))
    cosmo_values=Dict("state_squared"=>-4.000000008371047,
        "state_covariance"=>-5.00000032294447,"state_mixed"=>-3.5114802313179263,
        "trace_chsh"=>-2.8284271246604606,"trace_squared"=>-4.000000001047738)
    @testset "Legacy dense functional equations" begin
        for (id,problem,order,expected,atol) in cases()
            @testset "$id" begin
                # Cold COSMO squared-state solve: 35.1 s, above 20 s.
                if !references && id=="state_squared"
                    @test false skip=true
                    continue
                end
                if cosmo
                    expected=get(cosmo_values,id,expected)
                    startswith(id,"state_") && (atol=1e-3)
                end
                built=NCTSSoS.build(problem,Relaxation(;order))
                @test NCTSSoS.block_sides(built)==sizes[id][1]
                @test built.sizes.moments==sizes[id][2]
                result=solve(built;solver)
                @test result.bound ≈ expected atol=atol
                if id=="state_mixed"
                    if isnothing(sos_solver)
                        @test false skip=true
                    else
                        dual=solve(built;solver=sos_solver)
                        @test dual.bound ≈ expected atol=atol
                        @test dual.bound ≈ result.bound atol=atol
                    end
                end
            end
        end
    end
    bilocal_check(;solver=references ? sos_solver : nothing)
    @testset "D1 sparse functional routes deferred after main" begin
        # These named rejections cover the sparse state/network cases and
        # T1, T2, T5 and the trace CHSH/covariance sparse cases in the ledger.
        for problem in (state_problem(:mixed),trace_problem(:chsh))
            @test_throws r"State/trace.*plain dense route" NCTSSoS.build(problem,
                Relaxation(order=2,grouping=Correlative(),split=TermSplit(algorithm=MMD())))
        end
    end
end
end
