# # [Dense GNS Reconstruction](@id dense-gns-construction)
#
# In this page, we reconstruct dense GNS models from two inputs:
# 1. a dense Hankel (moment) matrix,
# 2. a direct `cs_nctssos(...; dualize=false)` solve result.
#
# A Hankel matrix stores moments `H[i, j] = <b_i^* b_j>` on a monomial basis.
# Here "dense" means one full moment block rather than clique-sparse data. The
# result-based path currently supports the package's `MonoidAlgebra` families
# (`NonCommutativeAlgebra`, `ProjectorAlgebra`, and `UnipotentAlgebra`) when the
# solve keeps a single dense moment block.

# ## Setup

using NCTSSoS
using COSMO
using JuMP
using LinearAlgebra

repr_with_registry(reg, mono) = sprint(show, mono; context=:registry => reg);

# ## 1) Reconstruct from a dense Hankel matrix
#
# This `7x7` fixture is the rounded flat Hankel extension from Example 2.7 of
# Klep, Povh, and Volcic, *Minimizer Extraction in Polynomial Optimization Is
# Robust* (SIAM J. Optim. 2018; preprint:
# <https://optimization-online.org/wp-content/uploads/2017/10/6249.pdf>). The
# same matrix is kept locally as `example_27_fixture()` in
# `test/relaxations/gns.jl`. The paper prints the extracted operators in an
# orthonormal quotient-space basis; below we show both the raw quotient-basis
# operators returned by `gns_reconstruct` and the orthonormalized matrices that
# can be compared directly with Example 2.7. The degree-2 dense basis is
# `(1, theta[1], theta[2], theta[1]^2, theta[1] * theta[2], theta[2] * theta[1],
# theta[2]^2)`.

reg, (theta,) = create_noncommutative_variables([("theta", 1:2)]);
basis_strings = repr_with_registry.(Ref(reg), get_ncbasis(reg, 2));
# `basis_strings` lists the degree-2 dense basis used to index the Hankel matrix.

basis_strings

#-

hankel = [
    1.00001 0.499907 0.500102 1.0483 -0.5483 -0.5483 1.0484;
    0.499907 1.0483 -0.548283 1.0627 -0.0144 -0.6090 0.0606;
    0.500102 -0.548283 1.04827 -0.0144 -0.5340 0.0606 0.9878;
    1.0483 1.0627 -0.0144 1.4622 -0.3995 -0.8006 0.7863;
    -0.5483 -0.0144 -0.5340 -0.3995 0.3852 0.1917 -0.7256;
    -0.5483 -0.6090 0.0606 -0.8006 0.1917 0.4411 -0.3804;
    1.0484 0.0606 0.9878 0.7863 -0.7256 -0.3804 1.3682;
]
# `hankel` is the dense moment matrix on `basis_strings`.

model = gns_reconstruct(hankel, reg, 2; atol=0.1, rtol=1e-6);
# `model` is a `GNSModel`; `model.report` stores flatness and moment diagnostics.

model.report

#-

@assert model.report.flat
@assert model.report.rank_full == 2
@assert model.report.rank_quotient == 2
@assert model.report.max_moment_error < 1e-3

quotient_basis_strings = repr_with_registry.(Ref(reg), model.quotient_basis);
# `quotient_basis_strings` are the basis vectors kept in the quotient space.

quotient_basis_strings

#-

model.cyclic_vector
# `model.cyclic_vector` is the canonical cyclic vector fixed by the identity basis element.

theta1_idx = only(variable_indices(theta[1]));
theta2_idx = only(variable_indices(theta[2]));

operator_action = (
    theta_1=round.(model.operators[theta1_idx] * model.cyclic_vector, digits=4),
    theta_2=round.(model.operators[theta2_idx] * model.cyclic_vector, digits=4),
);
# `operator_action` shows how each recovered generator acts on the cyclic vector.

operator_action

#-

quotient_operator_matrices = (
    theta_1=round.(model.operators[theta1_idx], digits=4),
    theta_2=round.(model.operators[theta2_idx], digits=4),
);
# `quotient_operator_matrices` are the recovered generator matrices in
# `model.quotient_basis` coordinates. This basis is generally not orthonormal, so
# these are not the symmetric matrices printed in Example 2.7.

quotient_operator_matrices

#-

gram_factor = cholesky(Symmetric(model.gram)).U
paper_sign = Diagonal([1.0, -1.0])

paper_operator_matrices = (
    theta_1=round.(paper_sign * (gram_factor * model.operators[theta1_idx] / gram_factor) * paper_sign, digits=4),
    theta_2=round.(paper_sign * (gram_factor * model.operators[theta2_idx] / gram_factor) * paper_sign, digits=4),
);
# `paper_operator_matrices` moves to an orthonormal basis of the quotient space
# (`gram_factor' * gram_factor == model.gram`) and then fixes the arbitrary sign
# of the second basis vector to match the convention used in Example 2.7.

paper_operator_matrices

#-

@assert isapprox(paper_operator_matrices.theta_1, [0.5019 -0.8931; -0.8931 0.1727]; atol=3e-3)
@assert isapprox(paper_operator_matrices.theta_2, [0.4981 0.8939; 0.8939 0.0825]; atol=3e-3)
# These are the rounded Example 2.7 matrices after accounting for basis
# normalization and the sign choice of one orthonormal basis vector.

# ## 2) Solve, extract dense moments, and reconstruct
#
# The dense result path starts from `cs_nctssos(...; dualize=false)`. The solve
# result carries one dense moment block, which `dense_moment_solution` turns into
# a reusable container for GNS reconstruction.

solver = optimizer_with_attributes(
    COSMO.Optimizer,
    "verbose" => false,
    "eps_abs" => 1e-7,
    "eps_rel" => 1e-7,
);

reg_u, (u,) = create_unipotent_variables([("u", 1:1)]);
pop = polyopt(-1.0 * u[1], reg_u)
# `pop` is a one-variable polynomial optimization problem over a unipotent generator.

result = cs_nctssos(pop, SolverConfig(optimizer=solver, order=1); dualize=false)
# `result` is a `PolyOptResult`; dense moment data is attached on the supported
# dense, non-dualized path.

result.objective

#-

sol = dense_moment_solution(result);
# `sol` is a `DenseMomentSolution`; it stores the dense basis, moments, and Hankel matrix.

sol_summary = (order=sol.order, basis=repr_with_registry.(Ref(reg_u), sol.basis));
sol_summary

#-

model_direct = gns_reconstruct(result; atol=1e-7, rtol=1e-7);
# `model_direct` is the dense `GNSModel` reconstructed straight from the solve result.

model_direct.report

#-

@assert model_direct.report.flat
@assert model_direct.quotient_basis == [one(u[1])]
@assert model_direct.cyclic_vector == [1.0]

u_idx = only(variable_indices(u[1]));
u_operator = round.(model_direct.operators[u_idx], digits=6);
u_operator
# This `1x1` operator is the recovered action of `u[1]` in the quotient space.

verify_gns(model_direct, sol)
# `verify_gns` recomputes the maximum moment mismatch from the reconstructed model.

# ## Limits of the dense path
#
# Dense result-based reconstruction currently assumes:
# - `dualize=false`,
# - one dense moment block,
# - a `MonoidAlgebra` registry.
#
# If those conditions are not met, `dense_moment_solution(result)` throws an
# `ArgumentError`.

# ## Next steps
#
# - Compare with the other Literate workflows under `docs/src/examples/literate/`.
# - See the user-interface API pages for `polyopt`, `SolverConfig`, and `cs_nctssos`.
