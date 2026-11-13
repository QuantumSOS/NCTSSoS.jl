using Test, NCTSSoS, LinearAlgebra

@testset "structured real RDM finite symmetry" begin
    h = Dict(0 => [0.0 -2.0; -2.0 0.0])
    v = zeros(2,2,2,2); v[1,1,1,1] = v[2,2,2,2] = 4.0
    options = (; nk=1,norb=2,nelec_per_cell=2,blocking=:spin,
        real_moments=true,include_one_d=true,three_index=:T1T2prime)
    original = build_pqg_moment_data(h,Dict((0,0,0,0)=>v); options...)
    spatial = [FermionicModePermutation(1=>2,2=>1,3=>4,4=>3)]
    flip = [FermionicModePermutation(1=>3,2=>4,3=>1,4=>2)]
    result = NCTSSoS.reduce_v2rdm_symmetry(original,spatial; moment_generators=flip)
    @test length(result.linear.moments) < length(original.moments)
    @test maximum(b.size for b in result.linear.psd_blocks_lin) < maximum(b.size for b in original.psd_blocks_lin)
    y = Dict(k=>sin(i) for (i,k) in enumerate(result.linear.moments))
    full = Dict(k=>(s==0 ? 0.0 : s*y[r]) for (k,(s,r)) in result.reconstruction)
    evaluate(form,values) = sum((c*values[k] for (k,c) in form); init=0.0)
    @test evaluate(original.objective_lin,full) ≈ evaluate(result.linear.objective_lin,y)
    next = 1
    for (block,bases) in zip(original.psd_blocks_lin,result.transforms)
        A = [evaluate(f,full) for f in block.entries]
        rebuilt = zeros(block.size,block.size)
        for U in bases
            B = [evaluate(f,y) for f in result.linear.psd_blocks_lin[next].entries]
            rebuilt += U'*B*U
            next += 1
        end
        @test rebuilt ≈ A atol=1e-10
    end
    broken_h = Dict(0 => [0.3 -2.0; -2.0 0.0])
    broken = build_pqg_moment_data(broken_h,Dict((0,0,0,0)=>v); options...)
    @test_throws ArgumentError NCTSSoS.reduce_v2rdm_symmetry(broken,spatial)
end

@testset "PQG through the polynomial symmetry and sparsity interface" begin
    for nk in (2, 3)
        # Include three momenta: modulo two cannot distinguish pair momentum
        # from momentum transfer. No structured assembler supplies this basis.
        reg, (a, ad) = create_fermionic_variables(1:2nk)
        pairs = [(i,j) for i in 1:2nk for j in i+1:2nk]
        basis = [one(a[1]);
            [only(monomials(a[i]*a[j])) for (i,j) in pairs];
            [only(monomials(ad[i]*ad[j])) for (i,j) in pairs];
            [only(monomials(ad[i]*a[j])) for i in 1:2nk for j in 1:2nk]]
        layout = FermionicModeLayout(
            Dict(i => (i+1)÷2 for i in 1:2nk);
            spin2_of=Dict(i => (isodd(i) ? 1 : -1) for i in 1:2nk),
            irrep_of=Dict(k+1 => k for k in 0:nk-1),
            irrep_table=AbelianIrrepTable(0;
                multiply=(x,y)->mod(x+y,nk), dual=x->mod(-x,nk)),
        )
        symmetry = SymmetrySpec(sector=FermionicSectorSpec(
            mode_layout=layout, split_number=true, split_spin=true, split_irrep=true))
        h = Dict(k => fill(-1.0+0.4k,1,1) for k in 0:nk-1)
        eri = Dict((k,l,k,l) => fill(ComplexF64(0.3+0.1k+0.2l),1,1,1,1)
            for k in 0:nk-1 for l in 0:nk-1)
        if nk==2
            eri[(0,0,1,1)] = fill(0.8im,1,1,1,1)
            eri[(1,1,0,0)] = fill(-0.8im,1,1,1,1)
        end
        ham = sum(h[k][1,1]/nk * ad[2k+s]*a[2k+s] for k in 0:nk-1 for s in 1:2)
        for ((k,l,r,t),v) in eri, s in 1:2, u in 1:2
            ham += v[1,1,1,1]/(2nk^2) * ad[2k+s]*ad[2l+u]*a[2t+u]*a[2r+s]
        end
        number = sum(ad[i]*a[i] for i in 1:2nk)
        pop = polyopt(ham, reg; moment_eq_constraints=[number-nk])
        cfg = SolverConfig(optimizer=SOLVER, moment_basis=basis, symmetry=symmetry)
        sparsity = compute_sparsity(pop,cfg)
        mp, report = NCTSSoS.moment_relax_symmetric(pop,
            sparsity.corr_sparsity,sparsity.cliques_term_sparsities,symmetry)
        neutral = NCTSSoS._fermionic_sector_label(one(a[1]),symmetry.sector)
        @test all(m->NCTSSoS._fermionic_sector_label(m,symmetry.sector)==neutral,mp.total_basis)
        reference = build_pqg_moment_data(h,eri; nk,norb=1,nelec_per_cell=1,
            blocking=:spin,spin_resolved_trace=false,singlet_s2=false)
        @test sort(report.psd_block_sizes) == sort([b.size for b in reference.psd_blocks_lin])

        # Compare every matrix polynomial after matching row order, not just
        # block counts or an optimal value on a symmetric example.
        full_basis = only(sparsity.cliques_term_sparsities[1][1].block_bases)
        matrices = [matrix for (cone,matrix) in mp.constraints if cone==:HPSD]
        for (label,matrix) in zip(report.block_labels,matrices)
            rows = filter(row->NCTSSoS._fermionic_sector_label(row,symmetry.sector)==label,full_basis)
            ref = only(filter(b->Set(b.meta.row_labels)==Set(rows),reference.psd_blocks_lin))
            order = [findfirst(==(row),ref.meta.row_labels) for row in rows]
            for i in eachindex(rows), j in eachindex(rows)
                expected = sum((c*reference.key_to_monomial[k]
                    for (k,c) in ref.entries[order[i],order[j]]);init=zero(ham))
                @test matrix[i,j] == expected
            end
        end

        nk==2 || continue
        # In the two-electron occupation basis, the paired states have these
        # energies and hopping amplitude. All four other states have -0.6875.
        expected_energy = eigmin(Hermitian([-0.925 0.2im; -0.2im -0.45]))
        for (cs,ts) in ((NoElimination(),NoElimination()),
                        (NoElimination(),MMD()), (MMD(),MMD()))
            sparse_cfg = SolverConfig(optimizer=SOLVER,moment_basis=basis,
                cs_algo=cs,ts_algo=ts,symmetry=symmetry)
            reduced = cs_nctssos(pop,sparse_cfg;dualize=false)
            plain = cs_nctssos(pop,SolverConfig(optimizer=SOLVER,moment_basis=basis,
                cs_algo=cs,ts_algo=ts);dualize=false)
            @test reduced.objective ≈ plain.objective atol=2e-6
            # This two-electron fixture is exact even with sparse blocks.
            # Arbitrary Hamiltonians are not promised dense/sparse equality.
            @test reduced.objective ≈ expected_energy atol=2e-6
            @test all(label->label isa FermionicSectorLabel,reduced.symmetry.block_labels)
        end
    end

    # Exercise genuinely separate orbital cliques, not just MMD on a complete
    # orbital graph. Number-sector phases preserve each clique individually.
    reg, (a,ad) = create_fermionic_variables(1:4)
    ham = -(ad[1]*a[2]+ad[2]*a[1]) - 2im*ad[3]*a[4] +
        2im*ad[4]*a[3] + 0.3*ad[3]*a[3]
    pop = polyopt(ham,reg)
    symmetry = SymmetrySpec(sector=FermionicSectorSpec(split_number=true))
    cfg = SolverConfig(optimizer=SOLVER,order=1,cs_algo=MMD(),ts_algo=MMD(),symmetry=symmetry)
    result = cs_nctssos(pop,cfg;dualize=false)
    @test length(result.sparsity.corr_sparsity.cliques)==2
    @test result.objective ≈ -1+eigmin(Hermitian([0.3 -2im; 2im 0.0])) atol=2e-6
    @test maximum(result.symmetry.psd_block_sizes)==2
    @test result.symmetry.basis_half_size==9
    @test cs_nctssos(pop,cfg).objective ≈ result.objective atol=2e-6
    # The sparse path must still reject an objective that changes charge.
    @test_throws ArgumentError cs_nctssos(polyopt(ham+a[1]+ad[1],reg),cfg;dualize=false)

    # Charged rows still determine the truncation degree: filtering them must
    # not drop the scalar number equality when identity is the only neutral row.
    reg, (a,ad) = create_fermionic_variables(1:2)
    number = sum(ad[i]*a[i] for i in 1:2)
    fixed = polyopt(number,reg;moment_eq_constraints=[number-1])
    cfg = SolverConfig(optimizer=SOLVER,moment_basis=[one(a[1]);a],symmetry=symmetry)
    @test cs_nctssos(fixed,cfg;dualize=false).objective ≈ 1 atol=2e-6
end
