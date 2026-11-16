# No solver or license needed. Run with the package environment:
# julia --project perf/next/bosonic.jl [cost|sector|narrow|fixed]
# `fixed` replaces only word_type in this disposable process for the width A/B.
using NCTSSoS, LinearAlgebra, Test
const N=NCTSSoS.Next

function samples(label,call;count=3)
    value=call()
    for _ in 1:count
        GC.gc()
        t=@timed call()
        println("SAMPLE ",label," seconds=",t.time," bytes=",t.bytes)
    end
    return value
end

function hamiltonian(b,bd,edges)
    number=bd.*b
    return -sum(bd[i]*b[j]+bd[j]*b[i] for (i,j) in edges)+sum(n*n-n for n in number)
end
function next_model(n,edges)
    _,(b,bd)=N.bosonic_variables(n)
    return N.Problem(hamiltonian(b,bd,edges))
end
function sector_cost()
    n=8; edges=[(i,mod1(i+1,n)) for i in 1:n]
    model=samples("next_sector_model",()->begin
        _,(b,bd)=N.bosonic_variables(n)
        (N.Problem(hamiltonian(b,bd,edges)),N.Relaxation(order=2,strengthen=N.Sector(sum(bd.*b),n)))
    end;count=5)
    built=samples("next_sector_assembly",()->N.build(model...);count=5)
    form=samples("next_sector_standard",()->N.moment_standard_form(built.program);count=5)
    println("SECTOR equations_next=",sum(b.side for b in built.program.blocks if b.kind==N.ZERO),
        " next_sizes=",built.sizes," standard=",size(form.A))
end
function kernel(s,a,b,repeats)
    total=0
    for _ in 1:repeats
        empty!(s); N.mul!(s,a,b,1.0)
        total+=length(s)
    end
    return total
end
function cost()
    edges=[(i,mod1(i+1,8)) for i in 1:8]
    p=samples("next_model",()->next_model(8,edges))
    built=samples("next_assembly",()->N.build(p,N.Relaxation(order=2)))
    form=samples("next_standard",()->N.moment_standard_form(built.program))
    println("DIMENSIONS ",built.sizes," standard=",size(form.A))
    W=N.word_type(N.Bosonic(),8,8)
    a=N.word_from_legacy(W,[1,1,8]); b=N.word_from_legacy(W,[-8,-1,-1])
    sink=N.PolynomialBuffer{W,Float64}(6)
    N.mul!(sink,a,b,1.0)
    @test kernel(sink,a,b,1)==6
    @test (@allocated kernel(sink,a,b,1))==0
    samples("next_product_10000",()->kernel(sink,a,b,10000))
end
function width(label)
    edges=[(4y+x,4y+x+1) for y in 0:3 for x in 1:3]
    append!(edges,[(4y+x,4(y+1)+x) for y in 0:2 for x in 1:4])
    p=next_model(16,edges)
    result=samples(label,()->N.build(p,N.Relaxation(order=2)))
    println("WIDTH ",label," word=",N.word_type(N.Bosonic(),16,4),
        " bytes=",sizeof(N.word_type(N.Bosonic(),16,4))," sizes=",result.sizes)
end

BLAS.set_num_threads(1)
println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," Julia_threads=",Threads.nthreads()," BLAS_threads=",BLAS.get_num_threads())
mode=isempty(ARGS) ? "cost" : only(ARGS)
if mode=="fixed"
    @eval N function word_type(::Bosonic,n::Integer,d::Integer)
        n>=0 && 0<=d<=255 || throw(ArgumentError("Bosonic capacity exceeded."))
        return BosonWord{8,max(1,cld(n,8))}
    end
end
mode=="cost" ? cost() : mode=="sector" ? sector_cost() : width(mode)
