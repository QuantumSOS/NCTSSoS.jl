# Run independently of the front-door compiler: julia --project prototype/pauli/test_strengthening.jl
module PauliStrengtheningChecks

import NCTSSoS as NC
import LinearAlgebra
using Test
include("strengthening.jl")
include("certificates.jl")

function dense(p)
    letters = ([0 1;1 0], [0 -im;im 0], [1 0;0 -1], [1 0;0 1])
    n = p.registry.nsites
    out = zeros(Complex{Rational{BigInt}}, 2^n, 2^n)
    for (c,w) in NC.terms(p)
        op = fill(4,n)
        for l in NC.legacy_letters(w)
            pos,a = divrem(l-1,3)
            op[pos+1] = a+1
        end
        out += c*reduce(kron, letters[a] for a in op)
    end
    return out
end
expect(p, ψ) = LinearAlgebra.dot(ψ,dense(p)*ψ)/sum(abs2,ψ)
word(p) = last(only(NC.terms(p)))

@testset "One-sided conserved sectors" begin
    _, (x,y,z) = NC.pauli_variables(2)
    H = z[1] + 2z[2]
    sector = workflow_strengthening(H,1;target_mz_sector=2,rdm_positivity_k=2)
    ψ = [1,0,0,0]
    @test all(p -> iszero(expect(p,ψ)), sector.equalities)
    @test expect(x[1]*(z[1]+z[2]-2)*x[1],ψ) == -2 # the incorrect sandwich excludes ψ
    @test any(p -> !iszero(expect(p,[0,0,0,1])), sector.equalities)
    @test [size(b.matrix,1) for b in sector.psd_blocks] == [1,2,1]
    @test [sum(diag[diagidx,diagidx] for diagidx in axes(diag,1))
           for diag in map(b -> expect.(b.matrix,Ref(ψ)),sector.psd_blocks)] == [1,0,0]
    @test_throws ArgumentError workflow_strengthening(H,1;target_spin_j=0)
    @test_throws ArgumentError workflow_strengthening(H,1;target_mz_sector=1)
    @test_throws ArgumentError workflow_strengthening(H,1;target_mz_sector=2,energy_target=:global_ground)
    heis = x[1]*x[2]+y[1]*y[2]+z[1]*z[2]
    singlet = workflow_strengthening(heis,1;target_spin_j=0)
    @test all(p -> iszero(expect(p,[0,1,-1,0])), singlet.equalities)
    @test any(p -> !iszero(expect(p,[0,1,1,0])), singlet.equalities)
    @test_throws ArgumentError workflow_strengthening(heis,1;target_spin_j=0,target_mz_sector=2)
    _, (_,_,singlez) = NC.pauli_variables(1)
    @test all(iszero,workflow_strengthening(singlez[1],1;target_spin_j=1//2).target.sectors)
end

@testset "KMS ground-space domain and stationarity" begin
    _, (x,y,z) = NC.pauli_variables(2)
    H = z[1] + 2z[2]
    kms = workflow_strengthening(H,1;commutator_linear_order=1,kms_psd_order=1,
                                 kms_rows=[x[1],y[1]],kms_frustration_guard=false)
    M = only(kms.psd_blocks).matrix
    @test expect.(M, Ref([0,0,0,1])) == [2 -2im; 2im 2]
    @test expect(M[1,1],[1,0,0,0]) == -2
    @test all(p -> iszero(expect(p,[0,0,0,1])), kms.equalities)
    @test any(p -> !iszero(expect(p,[1,1,0,0])), kms.equalities)
    @test_throws ArgumentError workflow_strengthening(H,2;target_mz_sector=0,kms_psd_order=2)
    @test_throws ArgumentError workflow_strengthening(H,2;target_mz_sector=0,kms_psd_order=2,kms_rows=x)
    exchange = x[1]*x[2]+y[1]*y[2]
    sector = workflow_strengthening(H,2;target_mz_sector=0,kms_psd_order=2,kms_rows=[exchange])
    @test expect(only(only(sector.psd_blocks).matrix),[0,1,0,0]) == 8
    near = workflow_strengthening(z[1]+z[2]*(1//10^9),1;kms_psd_order=1)
    @test size(only(near.psd_blocks).matrix) == (2,2)
end

@testset "Local density and support" begin
    _, (x,y,z) = NC.pauli_variables(2)
    H = z[1]+2y[2]
    wf = workflow_strengthening(H,1;lightcone_depth=2,rdm_regions=[(2,)],
                               rg_coarse_grain_depth=1,rg_maps=((depth=1,sites=(2,),isometry=reshape([0,im],2,1)),))
    @test Set(wf.rows) == Set([word(one(H)),word(z[1]),word(y[2]),word(z[1]*y[2])])
    ket = [1,3im,0,0]
    M, R = (b.matrix for b in wf.psd_blocks)
    @test expect.(M,Ref(ket)) == [1 -3im;3im 9].//10
    @test expect(only(R),ket) == 9//10
    @test_throws ArgumentError workflow_strengthening(H,1;rg_coarse_grain_depth=1)
    @test_throws ArgumentError workflow_strengthening(H,1;rg_coarse_grain_depth=1,
        rg_maps=((depth=1,sites=(2,),isometry=reshape([1,1],2,1)),))
    @test workflow_density_rewrite(H,word.([one(H),x[2],y[2],z[2]]),(2,)).matrix == M
    @test_throws ArgumentError workflow_density_rewrite(H,word.([one(H),x[2],y[2]]),(2,))

    # For a degree-six basis on 72 sites enumeration would be infeasible.
    _, (xl,yl,zl) = NC.pauli_variables(72)
    Ham = zl[65]+2zl[66]
    local_rows = word.([one(Ham),xl[65],yl[65],zl[65]])
    lattice = workflow_strengthening(Ham,6;moment_rows=[local_rows],
        commutator_linear_order=6,kms_psd_order=6,rdm_regions=[(65,66)])
    @test length(lattice.equalities) == 2
    @test [size(b.matrix) for b in lattice.psd_blocks] == [(2,2),(4,4)]
    @test all(w -> all(in((65,66)),NC.sites(w)),
        (w for b in lattice.psd_blocks for p in b.matrix for (_,w) in NC.terms(p)))
end

@testset "Verified rational bounds and asymmetric states" begin
    _, (x,y,z) = NC.pauli_variables(2)
    H = -z[1]
    gram = (rows=[x[1],y[1]], gram=[1 im;-im 1].*(1//2))
    cert = certify_rational_sos(H,-99//100;grams=(gram,))
    @test (cert.status,cert.bound,cert.residual_l1) == (:verified_rational,-1,1//100)
    @test_throws ArgumentError certify_rational_sos(H,-1;
        grams=((rows=gram.rows,gram=[1 1;1 1-1//10^12]),))
    @test_throws ArgumentError certify_rational_sos(H,-1;
        grams=((rows=gram.rows,gram=[0 im;-im 1]),))
    window = certify_energy_window(H,cert,[2,2im,0,0])
    @test (window.lower,window.upper) == (-1,-1)
    @test_throws ArgumentError workflow_strengthening(H+z[2],1;energy_window=window)
    @test_throws ArgumentError workflow_strengthening(H,1;energy_window=window,energy_target=:energy_window,commutator_linear_order=1)
    @test_throws ArgumentError workflow_strengthening(H,1;energy_window=window,energy_target=:energy_window,kms_psd_order=1)
    @test_throws ArgumentError workflow_strengthening(H,1;energy_upper_bound=-1.0)
    success = workflow_strengthening(H,1;energy_window=window,variance_center=-1,variance_upper_bound=0)
    @test length(success.inequalities) == 3
    @test all(p -> iszero(expect(p,[2,2im,0,0])), success.inequalities)
    @test_throws ArgumentError workflow_strengthening(H,1;energy_window=window,variance_center=0,variance_upper_bound=0)

    interval = certify_observable_interval(H,y[2],window;
        lower=(level=-1,squares=((one(H)+y[2])*(1//2),(one(H)+y[2])*(1//2))),
        upper=(level=-1,squares=((one(H)-y[2])*(1//2),(one(H)-y[2])*(1//2))))
    @test (interval.lower,interval.upper,interval.status) == (-1,1,:verified_rational)
    @test expect(y[2],[1,-im,0,0]) == interval.lower
    @test expect(y[2],[1,im,0,0]) == interval.upper # complex conjugation changes O, not H
    energy_observable = certify_observable_interval(H,H,window;
        lower=(level=-1,window_multipliers=(1,0)), upper=(level=1,window_multipliers=(0,1)))
    @test (energy_observable.lower,energy_observable.upper) == (-1,-1)
    @test_throws ArgumentError certify_observable_interval(H,H,window;
        lower=(level=-1,window_multipliers=(-1,0)), upper=(level=1,))

    _, (_,_,sz) = NC.pauli_variables(1)
    Hsec = sz[1] # selected +Z sector is an excited global eigenspace
    sc = certify_rational_sos(Hsec,1;target_mz_sector=1,sector_multipliers=(one(Hsec),))
    @test sc.bound == 1
    @test_throws ArgumentError certify_energy_window(Hsec,sc,[0,1])
    @test_throws ArgumentError certify_energy_window(Hsec,sc,[1,0];target_mz_sector=nothing)
    ws = certify_energy_window(Hsec,sc,[1,0])
    @test ws.upper == 1
    @test_throws ArgumentError workflow_strengthening(Hsec,1;energy_window=ws)
end

end
