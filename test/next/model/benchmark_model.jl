# Build-only benchmark: polynomial/registry construction is outside the timer;
# checked Problem/Relaxation construction is inside. No optimizer is attached or run.
include("../core/benchmark_build.jl")

model_program(H,order) = Nx.build(Nx.Problem(H),Nx.Relaxation(order=order))
model_standard(H,order) = Nx.moment_standard_form(model_program(H,order).program)
model_jump(H,order) = Nx.to_jump(model_standard(H,order))
function benchmark_model(name,H,order)
    built = model_program(H,order)
    W = eltype(built.table.keys)
    println(name," rows=",built.sizes.rows," keys=",built.sizes.moments," word_bytes=",sizeof(W))
    measure_build("$name model checked rows+program",model_program,H,order)
    measure_build("$name model checked through standard",model_standard,H,order)
    measure_build("$name model checked through JuMP",model_jump,H,order)
end
function model_benchmarks()
    println("Julia=",VERSION," CPU=",Sys.CPU_NAME," threads=",Threads.nthreads()," samples=7 GC outside timer")
    _,(x,y,z) = Nx.pauli_variables(10)
    H = sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:9)
    benchmark_model("Heisenberg10/d2",H,2)
    _,(a,ad) = Nx.fermionic_variables(4)
    H = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    benchmark_model("Hubbard2/d2",H,2)
    _,(x,y) = Nx.site_variables(Nx.Unipotent(),2,2)
    H = -(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2])
    benchmark_model("CHSH/d2",H,2)
end
abspath(PROGRAM_FILE)==(@__FILE__) && model_benchmarks()
