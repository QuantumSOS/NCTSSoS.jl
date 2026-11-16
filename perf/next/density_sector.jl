# Exact density-face acceptance and kernel cost. No solver for `kernel`.
# julia --project perf/next/density_sector.jl kernel
# julia --project perf/next/density_sector.jl check  # Mosek and SCS required
using NCTSSoS, LinearAlgebra, SparseArrays, Random, Test
import LinearAlgebra: tr
const N=NCTSSoS
BLAS.set_num_threads(1)

function kernel_cost(n)
    paulis=([0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    spins=[sum(foldl(kron,(s==j ? P : Matrix{ComplexF64}(I,2,2) for s in 1:n)) for j in 1:n)/2 for P in paulis]
    S2=Complex{Rational{BigInt}}.(sum(S*S for S in spins))
    operators=[1=>S2]
    face=N._density_kernel(operators)
    @test size(face.basis,2)==binomial(n,n÷2)-binomial(n,n÷2-1)
    @test all(iszero,sparse(S2)*sparse(face.basis))
    for _ in 1:3
        GC.gc(); t=@timed N._density_kernel(operators)
        println("KERNEL sites=",n," side=",2^n," nullity=",size(t.value.basis,2)," seconds=",t.time," bytes=",t.bytes)
    end
end

function physical_matrix(p,n)
    paulis=([1 0;0 1],[0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    word(w)=begin
        axes=Dict((l-1)÷3+1=>(l-1)%3+1 for l in N.legacy_letters(w))
        foldl(kron,(ComplexF64.(paulis[get(axes,s,0)+1]) for s in 1:n))
    end
    sum(c*word(w) for (c,w) in N.terms(p)),word
end

function check_sector(n,complex_axis)
    registry,(x,y,z)=N.pauli_variables(n)
    Q=complex_axis ? 3y[1]+4z[1] : ((sum(x))^2+(sum(y))^2+(sum(z))^2)/4
    q=complex_axis ? 5 : 0
    H=complex_axis ? z[2]+0.2y[1]*y[2] : sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    kinds=(N.DensityMatrix(1:n),N.Sector(Q,q))
    p=N.Problem(H); r=N.Relaxation(order=1,strengthen=kinds)
    built=N.build(p,r)
    W=eltype(built.table.representatives.keys)
    sectors=N._checked_sectors(p,kinds,N.Symmetry(),W)
    extra=N._pauli_extra(first(kinds),W,1)
    face=N._density_face(extra,sectors)
    B=ComplexF64.(face.basis)
    U=Matrix(qr(B).Q)[:,1:size(B,2)]
    M,word=physical_matrix(Q,n); M-=q*I
    @test norm(M*B)<1e-10
    e=eigen(Hermitian(M)); V=e.vectors[:,findall(v->abs(v)<1e-10,e.values)]
    @test size(V,2)==size(B,2)
    density=deepcopy(built.program)
    filter!(b->b.label==1 && b.kind in (N.PSD,N.HPSD),density.blocks)
    for seed in 1:5
        A=randn(MersenneTwister(seed),ComplexF64,size(V,2),size(V,2))
        rho=V*A*A'*V'; rho/=tr(rho)
        # Orthogonalize only for the physical spectrum comparison. The
        # production cone uses the exact nonorthogonal kernel basis.
        @test eigvals(Hermitian(U'*rho*U))≈last(eigvals(Hermitian(rho)),size(V,2)) atol=1e-12
        coords=[label.part==2 ? imag(tr(rho*word(label.word))) : real(tr(rho*word(label.word))) for label in built.table.labels]
        values=CombinationMatrix.cone_values(built.program,coords)
        @test minimum(values.spectra)>-1e-10
        @test maximum(abs,values.equalities)<1e-10
        @test CombinationMatrix.cone_values(density,coords).spectra≈eigvals(Hermitian(B'*rho*B)) atol=1e-10
    end
    # Keep the original density cone in the comparator, not the new helper.
    old=N.build(p,N.Relaxation(order=1,strengthen=(last(kinds),)))
    N._append_extra!(old.program,old.table,extra)
    old.program.nvars=length(old.table.labels)
    previous=NaN
    for (name,program) in (("old",old.program),("face",built.program))
        result=N.solve_standard(N.moment_standard_form(program),N.MosekBackend(MSK_IPAR_NUM_THREADS=1))
        println("COMPARE n=",n," complex=",complex_axis," route=",name," status=",result.status," bound=",result.bound)
        @test result.bound!==nothing
        name=="old" ? (previous=result.bound) : (@test abs(result.bound-previous)<1e-6)
    end
end

function check_controls()
    C=CombinationMatrix
    f=C.fixture(:ising,4)
    r=N.Relaxation(order=1,grouping=N.Correlative(),split=N.TermSplit(),strengthen=(N.DensityMatrix(1:4),N.Sector(f.Q,0)))
    p=N.Problem(f.H;symmetry=C.declaration(f,:antiunitary,:none))
    b=N.build(p,r)
    face=only(b.sizes.closure.strengthening.density_faces)
    @test (face.before,face.after)==(16,6)
    for backend in (N.MosekBackend(MSK_IPAR_NUM_THREADS=1),N.SCSBackend(verbose=0,eps_abs=1e-8,eps_rel=1e-8))
        result=N.solve_standard(N.moment_standard_form(b.program),backend)
        println("ISING solver=",typeof(backend)," status=",result.status," bound=",result.bound)
        @test result.bound!==nothing && abs(result.bound)<1e-6
    end
    registry,(x,y,z)=N.pauli_variables(2)
    g=N.clifford_cnot(N.word_type(registry,2),2,1,2)
    p=N.Problem(z[1]+x[2];symmetry=(g,))
    b=N.build(p,N.Relaxation(order=1,strengthen=(N.DensityMatrix([1]),N.Sector(z[1],1))))
    @test only(b.sizes.closure.strengthening.density_faces).after==1
    result=N.solve_standard(N.moment_standard_form(b.program),N.MosekBackend(MSK_IPAR_NUM_THREADS=1))
    println("CNOT status=",result.status," bound=",result.bound)
    @test result.bound≈0 atol=1e-6
end

if abspath(PROGRAM_FILE)==@__FILE__
    println("ENV Julia=",VERSION," CPU=",Sys.CPU_NAME," Julia_threads=",Threads.nthreads()," BLAS_threads=",BLAS.get_num_threads())
    mode=isempty(ARGS) ? "kernel" : only(ARGS)
    if mode=="kernel"
        for n in (4,6,8); kernel_cost(n); end
    elseif mode=="check"
        @eval import Mosek, SCS
        include(joinpath(@__DIR__,"../../test/next/pipeline/combination_matrix.jl"))
        @testset "Density sector acceptance" begin
            check_sector(2,true)
            check_sector(4,false)
            check_controls()
        end
    else
        error("Use kernel or check")
    end
end
