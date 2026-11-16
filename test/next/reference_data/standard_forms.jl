module StandardFormReference
using Test, SparseArrays, SHA, TOML
import NCTSSoS: Next as N
include("combination_inputs.jl")
const C = CombinationMatrix

# Preserve the reviewer's byte encoding, including removal of tautological
# zero-cone rows. The input snapshot is fixed independently of the live matrix.
function fingerprint(p,r)
    f = N.moment_standard_form(N.build(p,r).program)
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
    B = A[keep,:]; io = IOBuffer(); rows,cols,vals = findnz(B)
    write(io,Int64.(rows)); write(io,Int64.(cols)); write(io,vals)
    write(io,f.b[keep]); write(io,f.c); write(io,f.offset)
    for c in cones
        print(io,c)
    end
    string(bytes2hex(sha256(take!(io)))[1:16]," rows ",size(B,1)," tautologies ",dropped)
end

attempt(f) = try f() catch e; "ERROR " * first(sprint(showerror,e),60) end

function records()
    fixtures = Dict{Any,Any}()
    fx(problem,n) = get!(() -> C.fixture(problem,n),fixtures,(problem,n))
    out = String[]
    for (i,c) in enumerate(C.smoke_cases())
        f = fx(c.problem,c.n)
        value = attempt() do
            fingerprint(N.Problem(f.H;constraints=f.constraints,
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
            value = attempt(() -> fingerprint(N.Problem(f.H;symmetry=N.Symmetry(gens...;moment_only=(K,))),relax))
            push!(out,"anti $problem n$n $label $rl $value")
        end
    end
    return out
end

function check(path=joinpath(@__DIR__,"../../data/expectations/next_standard_forms.toml"))
    expected = TOML.parsefile(path)["cases"]
    @testset "Frozen standard forms and rejected combinations" begin
        actual = records()
        @test length(actual) == length(expected) == 242
        for (got,want) in zip(actual,expected)
            @test got == want["expected"]["record"]
        end
    end
end
end
