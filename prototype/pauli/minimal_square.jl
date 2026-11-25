"""
    compile_minimal_square(H, L; motifs=(), pair_shells=0, stationarity=true, rdm=true, stability=false)

Build the periodic square XXX relaxation in the input Hamiltonian's units.
The translation quotient represents all equivalent plaquette cones by one.
Return the original Pauli rows, final shared moment coordinates, and explicit
moment-minimization standard form (not an energy certificate).

`physical_witness3x3=true` appends the certified local inequality to `standard`
only, with no new coordinates. This requires coverage of its Pauli observables
and is incompatible with Gram lift recovery. `program` remains the pre-cut program.
"""
function compile_minimal_square(H, L; motifs=(), pair_shells=0, stationarity=true, rdm=true,
        rdm_shape=(2,2), additional_rdm_offsets=(), stability=false, recover_gram_rows=false,
        recover_folded_gram=false, physical_witness3x3=false)
    if physical_witness3x3
        L >= 3 || throw(ArgumentError("The 3×3 physical witness requires L ≥ 3."))
        !(recover_gram_rows || recover_folded_gram) ||
            throw(ArgumentError("The 3×3 physical witness does not support Gram lift recovery."))
    end
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
    if physical_witness3x3
        standard = _append_square_witness3x3(standard,table,constraint_ids,L,eltype(seed.rows))
    end
    (;seed.rows,seed.domain,table,constraint_ids,standard,program=pruned,symmetry,pruning,lifts)
end

"""
    compile_square_xxx_witness(H, L)

Compile the periodic nearest-neighbor XXX recipe used for the 16×16 run,
in `H`'s units. The compiled `standard` includes the certified
3×3 physical inequality. No optimizer, saved point or coordinate map is used.
"""
function compile_square_xxx_witness(H, L)
    motifs = (((0,0),(1,0),(0,1)), ((0,0),(1,0),(2,0)))
    additional_rdm_offsets = (((0,0),(1,0),(2,0),(3,0)), ((0,0),(-1,0),(1,0),(0,-1),(0,1)))
    return compile_minimal_square(H,L; motifs, pair_shells=2, rdm_shape=(2,3),
        additional_rdm_offsets, physical_witness3x3=true)
end
