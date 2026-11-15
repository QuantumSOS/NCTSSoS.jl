if get(ENV,"PROJECTION_ONLY","false")=="true"
    path=joinpath(@__DIR__,"validate.jl")
    include_string(Main,first(split(read(path,String),"@testset \"Small-system")),path)
else
    include("validate.jl")
end

@testset "Sparse coefficient projection" begin
    reg,(x,y,z)=create_pauli_variables(1:2)
    rows=get_ncbasis(reg,2); C=collect(eachindex(rows))
    H=x[1]+2y[2]+3z[1]*x[2]
    coords=E.moment_coordinates(H,rows,[C],CliffordSymmetry[])
    model=Model(); @variable(model,t[1:coords.nvars])
    U=[complex(sin(3i+j),cos(i+2j)) for i in 1:5,j in eachindex(rows)]
    got=only(E.projected_forms(rows,C,[U],coords,t))
    M=[begin
        phase,w=E.product(a,b);k,s=coords.coordinates[w];phase*s*t[k]
    end for a in rows,b in rows]
    expected=U*M*U'
    @test maximum(abs(coefficient(got[i,j]-expected[i,j],v))
        for i in 1:5,j in 1:5,v in t)<1e-10

    reg,ops=create_pauli_variables(1:4)
    H=heisenberg_chain_hamiltonian(ops;coupling=1/4,periodic=true)
    gs=sympleq_symmetry_spec(H).clifford_generators
    rows=E.close_rows(get_ncbasis(reg,2),gs)
    groups=[collect(eachindex(rows))]
    coords=E.moment_coordinates(H,rows,groups,gs)
    actions=E.actions_on(rows,gs)
    plans=[(o.representative,E.decompose_signed(o.stabilizers,length(o.representative)))
        for o in E.block_orbits(groups,actions,length(rows))]
    bounds=Float64[]
    for numeric_projection in (false,true)
        built=E.assemble(H,rows,groups,gs,OPT;dual_form=true,
            coordinates=coords,block_plans=plans,numeric_projection)
        optimize!(built.model)
        @test termination_status(built.model)==MOI.OPTIMAL
        push!(bounds,objective_value(built.model))
        values=built.moment_values()
        M=[begin
            phase,w=E.product(a,b);k,s=coords.coordinates[w]
            k==0 ? 0.0im : phase*s*values[k]
        end for a in rows,b in rows]
        @test eigmin(Hermitian(M))>=-2e-6
    end
    @test abs(bounds[1]-bounds[2])<1e-7
    @test abs(last(bounds)+2)<1e-7
    println("PROJECTION_BOUNDS ",bounds)
end
println("PROJECTION_CHECK_PASSED")
