# Public contracts from test/correlated_sparsity/core_pipeline_numeric.jl
# and e9_chained_singular_polydisc.jl. Literature fixtures stay unchanged.
function correlated_chain(x)
    X,Y,Z=x
    return X^2-X*Y-Y*X+3Y^2-2X*Y*X+2X*Y^2*X-Y*Z-Z*Y+6Z^2+
        9Y^2*Z+9Z*Y^2-54Z*Y*Z+142Z*Y^2*Z
end

function correlated_polyball()
    _,(x,)=site_variables(Free(),1,4)
    X1,X2,X3,X4=x
    f1 =
        4.0 - X1 + 3.0 * X2 - 3.0 * X3 - 3.0 * X1^2 - 7.0 * X1 * X2 + 6.0 * X1 * X3 -
        X2 * X1 - 5.0 * X3 * X1 + 5.0 * X3 * X2 - 5.0 * X1^3 - 3.0 * X1^2 * X3 +
        4.0 * X1 * X2 * X1 - 6.0 * X1 * X2 * X3 + 7.0 * X1 * X3 * X1 +
        2.0 * X1 * X3 * X2 - X1 * X3^2 - X2 * X1^2 + 3.0 * X2 * X1 * X2 -
        X2 * X1 * X3 - 2.0 * X2^3 - 5.0 * X2^2 * X3 - 4.0 * X2 * X3^2 -
        5.0 * X3 * X1^2 + 7.0 * X3 * X1 * X2 + 6.0 * X3 * X2 * X1 -
        4.0 * X3 * X2^2 - X3^2 * X1 - 2.0 * X3^2 * X2 + 7.0 * X3^3
    f2 =
        -1.0 + 6.0 * X2 + 5.0 * X3 + 3.0 * X4 - 5.0 * X2^2 + 2.0 * X2 * X3 +
        4.0 * X2 * X4 - 4.0 * X3 * X2 + X3^2 - X3 * X4 + X4 * X2 - X4 * X3 +
        2.0 * X4^2 - 7.0 * X2^3 + 4.0 * X2 * X3^2 + 5.0 * X2 * X3 * X4 -
        7.0 * X2 * X4 * X3 - 7.0 * X2 * X4^2 + X3 * X2^2 + 6.0 * X3 * X2 * X3 -
        6.0 * X3 * X2 * X4 - 3.0 * X3^2 * X2 - 7.0 * X3^2 * X4 + 6.0 * X3 * X4 * X2 -
        3.0 * X3 * X4 * X3 - 7.0 * X3 * X4^2 + 3.0 * X4 * X2^2 - 7.0 * X4 * X2 * X3 -
        X4 * X2 * X4 - 5.0 * X4 * X3^2 + 7.0 * X4 * X3 * X4 + 6.0 * X4^2 * X2 -
        4.0 * X4^3
    return Problem(f1+f1'+f2+f2';constraints=(
        Positive(1-X1^2-X2^2-X3^2),Positive(1-X2^2-X3^2-X4^2)))
end

function block_letter_positions(built,count)
    registry=built.problem.objective.registry
    letters=collect(Next.basis(registry,1))[2:end]
    position=Dict(Next.operator_name(registry,last(only(Next.terms(p))))=>i for (i,p) in enumerate(letters))
    return sort!([sort!(unique!([position[Next.operator_name(registry,l)] for w in block for l in Next.word_letters(w)]))
        for block in collect(Next.block_rows(built))[1:count]])
end

function correlated_polydisc(n)
    _,(x,)=site_variables(Free(),1,n)
    f=zero(x[1])
    for i in 1:2:n-3
        a=x[i]+10x[i+1]; b=x[i+2]-x[i+3]
        c=x[i+1]^2-4x[i+1]*x[i+2]+4x[i+2]^2
        d=x[i]^2-20x[i]*x[i+3]+100x[i+3]^2
        f+=a'*a+5b'*b+c'*c+10d'*d
    end
    return Problem(f;constraints=Positive.([1 .- x.^2;x .- 1/3]))
end

function correlated_structure()
    @testset "Bipartite letter groups" begin
        _,(a,)=site_variables(Free(),1,6)
        x,y=a[1:3],a[4:6]
        f=-(x[1]*(y[1]+y[2]+y[3])+x[2]*(y[1]+y[2]-y[3])+
            x[3]*(y[1]-y[2])-x[1]-2y[1]-y[2])
        @test_throws ArgumentError Problem(f)
        problem=Problem((f+f')/2)
        for (key,algorithm) in (("mf",MF()),("as_is_elimination",AsIsElimination()))
            built=Next.build(problem,Relaxation(order=3,grouping=Correlative(;algorithm)))
            counts=[length(unique([l for w in block for l in Next.word_letters(w)])) for block in Next.block_rows(built)]
            @test maximum(counts)==fixture("correlated_structure","correlative_sparsity_example2_order3")["max_clique_by_algo"][key]
        end
    end
    @testset "Correlative localizers and overlapping groups" begin
        _,(x,)=site_variables(Free(),1,2)
        problem=Problem(2-x[1]^2+x[1]*x[2]^2*x[1]-x[2]^2;
            constraints=(Positive(4-x[1]^2-x[2]^2),Zero(x[1]*x[2]+x[2]*x[1]-2)))
        for algorithm in (NoElimination(),MF(),AsIsElimination())
            built=Next.build(problem,Relaxation(order=2,grouping=Correlative(;algorithm)))
            @test length.(collect(Next.block_rows(built)))==[7,3,3]
            polyball=Next.build(correlated_polyball(),Relaxation(order=2,grouping=Correlative(;algorithm)))
            groups=algorithm isa NoElimination ? [[1,2,3,4]] : [[1,2,3],[2,3,4]]
            @test block_letter_positions(polyball,length(groups))==groups
        end
    end
end

function correlated_numeric(solver,sos_solver)
    @testset "Correlative and term-sparse chain" begin
        _,(x,)=site_variables(Free(),1,3)
        problem=Problem(correlated_chain(x);constraints=Positive.([1 .- x.^2;x .- 1/3]))
        cs=Next.build(problem,Relaxation(order=3,grouping=Correlative(algorithm=MF())))
        @test Next.block_sides(cs)==[15,15,fill(7,8)...]
        @test cs.sizes.moments==149
        @test solve(cs;solver).bound ≈ fixture("correlated_pipeline","CS_d3")["objective"] atol=1e-5
        @test solve(cs;solver=sos_solver).bound ≈ fixture("correlated_pipeline","CS_d3")["objective"] atol=1e-5 skip=isnothing(sos_solver)
        higher=Next.build(problem,Relaxation(order=3,split=TermSplit(algorithm=MMD(),level=2)))
        # The old counter excludes localizer-only moments (123); both full
        # programs contain 288 scalar moment variables.
        @test higher.sizes.moments==288
        @test solve(higher;solver).bound ≈ fixture("correlated_pipeline","TS_d3")["objective"] atol=1e-5
        first=Next.build(problem,Relaxation(order=3,split=TermSplit(algorithm=MMD())))
        @test solve(first;solver=sos_solver).bound ≈ fixture("correlated_pipeline","TS_d3")["objective"] atol=1e-4 skip=isnothing(sos_solver)
    end
    @testset "Chordal fill reaches the exact minimum" begin
        _,(x,)=site_variables(Free(),1,3)
        X,Y,Z=x
        problem=Problem(correlated_chain(x))
        rows=[one(X),X,Y,Z,Y*X,Y*Z]
        for (level,sides) in ((1,[3,3,3,3]),(2,[4,4,4]))
            built=Next.build(problem,Relaxation(;rows,split=TermSplit(algorithm=MMD(),level=level)))
            @test sort(Next.block_sides(built))==sides
            bound=solve(built;solver).bound
            if level==1
                @test bound < -0.003
            else
                @test bound ≈ 0 atol=1e-6
            end
        end
    end
end

function correlated_sos(sos_solver)
    @testset "E8 polyball hierarchy" begin
        bounds=Float64[]
        for (id,relaxation) in (
            ("E8_sparse_d2",Relaxation(order=2,grouping=Correlative(algorithm=MF()))),
            ("E8_sparse_d3",Relaxation(order=3,grouping=Correlative(algorithm=MF()))),
            ("E8_dense_d2",Relaxation(order=2)))
            if isnothing(sos_solver)
                @test false skip=true
                continue
            end
            bound=solve(correlated_polyball(),relaxation;solver=sos_solver).bound
            @test bound ≈ fixture("correlated_pipeline",id)["objective"] atol=2e-6
            push!(bounds,bound)
        end
        if !isempty(bounds)
            @test bounds[1] < bounds[2]
            @test bounds[2] ≈ bounds[3] atol=2e-6
        end
    end
    @testset "E9 constrained chained singular" begin
        bounds=Float64[]
        for (n,sparse) in ((4,false),(4,true),(8,true),(12,true))
            problem=correlated_polydisc(n)
            relaxation=Relaxation(order=2,grouping=sparse ? Correlative(algorithm=MF()) : nothing)
            built=Next.build(problem,relaxation)
            id="E9_$(sparse ? "sparse" : "dense")_n$(n)_d2"
            oracle=fixture("correlated_pipeline",id)
            @test Next.block_sides(built)==[oracle["sides"];fill(sparse ? 4 : 5,sparse ? 6(n-2) : 2n)]
            @test built.sizes.moments==oracle["n_unique_moment_matrix_elements"]
            expected=fixture("correlated_structure","correlative_sparsity_e9_polydisc_n$(n)_order2")
            @test block_letter_positions(built,length(oracle["sides"]))==expected[sparse ? "mf" : "no_elimination"]["cliques"]
            if isnothing(sos_solver)
                @test false skip=true
            else
                bound=solve(built;solver=sos_solver).bound
                @test bound ≈ oracle["objective"] atol=1e-2
                push!(bounds,bound)
            end
        end
        isempty(bounds) || @test bounds[1] ≈ bounds[2] atol=1e-3
    end
end
