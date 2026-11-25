"""
    compile_minimal_square(H, L; motifs=(), pair_shells=0, stationarity=true, rdm=true, stability=false)

Build the periodic square XXX relaxation in the input Hamiltonian's units.
The translation quotient represents all equivalent plaquette cones by one.
Return the original Pauli rows, final shared moment coordinates, and explicit
moment-minimization standard form (not an energy certificate).
"""
function compile_minimal_square(H, L; motifs=(), pair_shells=0, stationarity=true, rdm=true,
        rdm_shape=(2,2), additional_rdm_offsets=(), stability=false, recover_gram_rows=false, recover_folded_gram=false)
    seed = build_square_xxx_seed(H,L; motifs, pair_shells, stationarity, rdm, rdm_shape, additional_rdm_offsets, stability, translation_quotient=true)
    program,table,symmetry = workflow_square_symmetry_program(H,seed.rows,L;
        seed.strengthening, recover_gram_rows, recover_folded_gram)
    real_program = _real_workflow_program(program)
    pruned,pruning = workflow_term_sparsity(real_program; ts_algo=:Maximal, prune_constant_singletons=true)
    NC.refresh_coverage!(table,pruned)
    standard = NC.moment_standard_form(pruned)
    A = standard.A
    constraint_ids = BitSet((1,))
    for j in axes(A,2)
        any(!iszero(A.nzval[i]) for i in A.colptr[j]:A.colptr[j+1]-1) && push!(constraint_ids,j+1)
    end
    lifts = recover_gram_rows || recover_folded_gram ?
        ()->square_lifts(H,seed,program,real_program,pruned,pruning,symmetry;form=standard) : nothing
    (;seed.rows,seed.domain,table,constraint_ids,standard,program=pruned,symmetry,pruning,lifts)
end
