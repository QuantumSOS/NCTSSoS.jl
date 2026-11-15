# Standalone acceptance product; intentionally not included by Pkg.test.
# julia --project test/next/pipeline/run_combination_matrix.jl --full results.tsv
# Add --solve to check the selected bounds/queries (requires SCS in the environment).
include("combination_matrix.jl")
const C=CombinationMatrix
"--solve" in ARGS && (@eval import SCS)

function main(args)
    cases="--full" in args ? C.combinations() : C.smoke_cases()
    solved=Set(C.solved_cases())
    "--solve" in args && !("--full" in args) && (cases=unique([cases;C.solved_cases()]))
    problemarg=findfirst(a->startswith(a,"--problem="),args)
    isnothing(problemarg) || (cases=filter(c->c.problem==Symbol(split(args[problemarg],'=')[2]),cases))
    finitearg=findfirst(a->startswith(a,"--finite="),args)
    isnothing(finitearg) || (cases=filter(c->c.finite==Symbol(split(args[finitearg],'=')[2]),cases))
    countarg=findfirst(a->startswith(a,"--limit="),args)
    isnothing(countarg) || (cases=cases[1:min(length(cases),parse(Int,split(args[countarg],'=')[2]))])
    paths=filter(a->!startswith(a,"--"),args)
    path=isempty(paths) ? "combination_results.tsv" : only(paths)
    solver=if "--solve" in args
        C.Nx.SCSBackend(verbose=0,eps_abs=1e-8,eps_rel=1e-8)
    else
        nothing
    end
    fixtures=Dict((c.problem,c.n)=>C.fixture(c.problem,c.n) for c in unique(c->(c.problem,c.n),cases))
    counts=Dict("program"=>0,"ArgumentError"=>0,"DEFECT"=>0)
    unaccepted=0
    started=time()
    open(path,"w") do io
        C.write_table(io,[])
        for (i,c) in enumerate(cases)
            selected=c in solved
            r=C.run(c;f=fixtures[(c.problem,c.n)],solver=selected ? solver : nothing)
            C.write_table(io,[r];header=false)
            counts[r.outcome]+=1
            r.solve_status in ("not_selected","not_checked","OPTIMAL/OPTIMAL") || (unaccepted+=1)
            r.outcome=="DEFECT" && println(stderr,"DEFECT at $i: $c\n$(r.message)")
            if i%100==0 || i==length(cases)
                flush(io)
                println(stderr,"$i/$(length(cases)) $(round(time()-started;digits=2)) s $counts")
                flush(stderr)
            end
        end
    end
    print(read(path,String))
    counts["DEFECT"]==0 || error("Combination matrix has $(counts["DEFECT"]) defects; see $path")
    unaccepted==0 || error("$unaccepted selected solves have no accepted bound; they are not passes. See $path")
end

main(ARGS)
