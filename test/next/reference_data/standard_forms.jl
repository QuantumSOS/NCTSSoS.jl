module StandardFormReference
using Test, SparseArrays, SHA, TOML, JuMP, COSMO
import NCTSSoS as N
include("combination_inputs.jl")
const C = CombinationMatrix

# Preserve the reviewer's byte encoding, including removal of tautological
# zero-cone rows. The input snapshot is fixed independently of the live matrix.
fingerprint(p,r) = fingerprint(N.moment_standard_form(N.build(p,r).program))

function fingerprint(f; comparison="exact")
    A = f.A
    rownz = zeros(Int,size(A,1))
    for (row,v) in zip(A.rowval,A.nzval)
        iszero(v) || (rownz[row] += 1)
    end
    keep = trues(size(A,1)); dropped = 0
    cones = Tuple{Symbol,Int}[]
    for c in f.cones
        rows = c.first:c.first+c.length-1
        if c.kind == :zero
            for row in rows
                (rownz[row] == 0 && iszero(f.b[row])) && (keep[row] = false; dropped += 1)
            end
        end
        n = count(@view keep[rows])
        n > 0 && push!(cones,(c.kind,n))
    end
    B = A[keep,:]
    if comparison == "invariant"
        sides = sort!([(c.kind,c.side) for c in f.cones])
        return string("rows ",size(B,1)," tautologies ",dropped,
            " variables ",size(B,2)," cones ",sides)
    end
    io = IOBuffer(); rows,cols,vals = findnz(B)
    if comparison == "rounded"
        # Numerical support and values, not the BLAS-specific low bits.
        roundvalue(v) = abs(v) <= 1e-11 ? 0.0 : round(v;digits=9) + 0.0
        support = abs.(vals) .> 1e-11
        write(io,Int64.(rows[support])); write(io,Int64.(cols[support]))
        write(io,roundvalue.(vals[support]))
        write(io,roundvalue.(f.b[keep])); write(io,roundvalue.(f.c)); write(io,roundvalue(f.offset))
        for c in cones
            print(io,c)
        end
        return string(bytes2hex(sha256(take!(io)))[1:16]," rows ",size(B,1),
            " cols ",size(B,2)," tautologies ",dropped)
    end
    comparison == "exact" || error("Unknown comparison: $comparison")
    write(io,Int64.(rows)); write(io,Int64.(cols)); write(io,vals)
    write(io,f.b[keep]); write(io,f.c); write(io,f.offset)
    for c in cones
        print(io,c)
    end
    string(bytes2hex(sha256(take!(io)))[1:16]," rows ",size(B,1)," tautologies ",dropped)
end

attempt(f) = try f() catch e; "ERROR " * first(sprint(showerror,e),60) end

function records(record_case=(i,p,r)->fingerprint(p,r))
    fixtures = Dict{Any,Any}()
    fx(problem,n) = get!(() -> C.fixture(problem,n),fixtures,(problem,n))
    out = String[]
    for (i,c) in enumerate(C.smoke_cases())
        f = fx(c.problem,c.n)
        value = attempt() do
            record_case(length(out)+1,N.Problem(f.H;constraints=f.constraints,
                symmetry=C.declaration(f,c.finite,c.continuous)),
                C.relaxation(f,c.family,C.strengthening(f,c.strengthening)))
        end
        push!(out,"smoke $i $value")
    end
    for problem in (:heisenberg,:xxz,:ising), n in (4,6)
        f = fx(problem,n)
        t = N.SitePermutation([2:n;1]); r = N.SitePermutation([1;n:-1:2])
        K = N.AntiunitaryPauli(one(eltype(f.H.data.words)))
        flip = N.CliffordMap(N.PauliConjugation(first(prod(f.operators[1:n]).data.words),false))
        for (label,gens) in (("K",()),("r+K",(r,)),("t+r+K",(t,r)),
                ("flip+K",(flip,)),("t+r+flip+K",(t,r,flip))),
            (rl,relax) in (("o1",N.Relaxation(order=1)),("o2",N.Relaxation(order=2)),
                ("o2 grouped split",N.Relaxation(order=2,grouping=N.Correlative(),split=N.TermSplit())),
                ("o1 density",N.Relaxation(order=1,strengthen=(N.DensityMatrix(1:3),))))
            problem == :ising && occursin("flip",label) && continue
            value = attempt(() -> record_case(length(out)+1,N.Problem(f.H;symmetry=N.Symmetry(gens...;moment_only=(K,))),relax))
            push!(out,"anti $problem n$n $label $rl $value")
        end
    end
    return out
end

function check(path=joinpath(@__DIR__,"../../data/expectations/next_standard_forms.toml");
    reference_solver::Union{Nothing,N.MosekBackend}=nothing)
    expected = TOML.parsefile(path)["cases"]
    # The full suite has already compiled this writer. All twenty small solves
    # together take under a second warm. The licensed Mosek suite passes its
    # backend explicitly to retain the stricter reference check on these cases.
    backend = isnothing(reference_solver) ? N.JuMPBackend(optimizer_with_attributes(COSMO.Optimizer,
        "verbose"=>false,"decompose"=>false,"eps_abs"=>1e-8,
        "eps_rel"=>1e-8,"max_iter"=>50_000)) : reference_solver
    atol = isnothing(reference_solver) ? 1e-5 : 1e-6
    @testset "Frozen standard forms and rejected combinations" begin
        actual = records() do i,p,r
            want = expected[i]
            comparison = get(want,"comparison","exact")
            f = N.moment_standard_form(N.build(p,r).program)
            if comparison == "invariant"
                @testset "$(want["id"])" begin
                    solution = N.solve_standard(f,backend)
                    @test solution.status == MOI.OPTIMAL
                    # First-order COSMO on degenerate density-matrix blocks:
                    # measured error 4.5e-7 across four kernels, 2.74e-6 in CI.
                    # Licensed Mosek retains atol=1e-6 on the same twenty forms.
                    @test solution.bound ≈ want["expected"]["objective"] atol=atol rtol=0
                end
            end
            fingerprint(f;comparison)
        end
        @test length(actual) == length(expected) == 242
        for (got,want) in zip(actual,expected)
            @test got == want["expected"]["record"]
        end
    end
end
end
