using NCTSSoS, LinearAlgebra
include("functional-sparse-models.jl")
const M=SparseFunctionalModels
const N=NCTSSoS.Next
BLAS.set_num_threads(1)
for (name,model,order) in (("trace-chsh",M.chsh(true),1),("bilocal",M.bilocal(),3),("t2-n20",M.t2(20),2))
    for sparse in (false,true)
        cs=sparse && startswith(name,"t2") ? MF() : NoElimination()
        ts=sparse ? MMD() : NoElimination()
        rel=N.BlockRelaxation(;order,grouping=sparse && startswith(name,"t2") ? N.Correlative() : nothing,
            split=sparse ? N.TermSplit(algorithm=MMD()) : nothing)
        for route in (:old,:new)
            GC.gc()
            data=@timed (route==:old ? M.old_build(model.old,order,cs,ts) : N.build(model.new,rel))
            result=data.value
            sides=route==:old ? [size(m,1) for (_,m) in result.moment.constraints] : N.block_sides(result)
            println((;name,sparse,route,seconds=data.time,bytes=data.bytes,blocks=length(sides),maxside=maximum(sides)))
            flush(stdout)
        end
    end
end
