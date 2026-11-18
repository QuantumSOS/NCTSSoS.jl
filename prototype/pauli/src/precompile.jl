using PrecompileTools: @compile_workload

# Exercise symbolic compilation, without a solver, measured input, or cached plan.
@compile_workload let
    _, axes = NC.pauli_variables(3)
    h = sum(a[1]*a[2] + 2a[2]*a[3] for a in axes)
    compile_relaxation(h; order=2, preset=:ExactNPA, moi_form=:auto)
end

# Load the periodic sparse assembly path with a different Hamiltonian from the
# comparison cases. PauliWord{1} and PauliWord{2} need separate machine code.
@compile_workload for n in (5,65)
    _, (x,y,z) = NC.pauli_variables(n)
    h = sum(x[i]*x[mod1(i+1,n)] + y[i]*y[mod1(i+1,n)] + 2z[i]*z[mod1(i+1,n)] for i in 1:n)
    compile_relaxation(h; order=2, moi_form=:dualized,
        symmetry=PauliSymmetrySpec(spatial_mode=:circulant_dft),
        sparsity=PauliSparsitySpec(support_mode=:contiguous_lattice,two_body_range_r=1))
end
