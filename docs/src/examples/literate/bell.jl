# # [Bell inequalities](@id bell-inequalities)

# Bell inequalities are mathematical expressions that test whether the predictions of quantum mechanics can be explained by local hidden variable theories. They were first introduced by John Stewart Bell in 1964 and have since become fundamental tools in quantum information theory and quantum foundations.
# A Bell inequality is typically expressed as a linear combination of expectation values of observables, with bounds that differ between classical and quantum theories. In the classical case, these inequalities must be satisfied if the system can be described by local hidden variables. However, quantum mechanics can violate these inequalities, demonstrating the non-local nature of quantum correlations.
# The general form of a Bell inequality can be written as:
#
# ```math
# \sum_{i,j} c_{ij} \langle A_i B_j \rangle \leq C
# ```
#
# where $A_i$ and $B_j$ are observables measured by two parties (traditionally called Alice and Bob), $c_{ij}$ are real coefficients, and $C$ is the classical bound. Quantum mechanics can violate this inequality, with the maximum violation known as the quantum bound.

# ## Linear Bell inequalities

# ### CHSH inequality
# The most famous Bell inequality is the CHSH (Clauser-Horne-Shimony-Holt) inequality, which involves two parties, each measuring two observables. For unipotent (squared to 1) observables $A_1, A_2$ measured by Alice and $B_1, B_2$ measured by Bob. We define the objective function as:
#
# ```math
# f(A_1, A_2, B_1, B_2) = \langle A_1B_1 \rangle + \langle A_1B_2 \rangle + \langle A_2B_1 \rangle - \langle A_2B_2 \rangle
# ```
#
# The CHSH inequality is then given by $$f(A_1, A_2, B_1, B_2) \leq 2,$$ which must be satisfied by any local hidden variable theory. However, quantum mechanics can violate this inequality up to the value $2\sqrt{2}$, known as Tsirelson's bound. This violation demonstrates that quantum mechanics cannot be described by any local hidden variable theory.
# The CHSH inequality is particularly important because it is the simplest non-trivial Bell inequality and has been experimentally verified numerous times, providing strong evidence for the non-local nature of quantum mechanics.

# An upper bound on the maximal quantum violation of the CHSH inequality can be computed using the following code:

using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)

# Create unipotent variables (operators that square to identity)
# x = (A_1, A_2), y = (B_1, B_2)
registry, (x, y) = create_unipotent_variables([("x", 1:2), ("y", 1:2)])
f = 1.0 * x[1] * y[1] + x[1] * y[2] + x[2] * y[1] - x[2] * y[2]  # objective function

# The registry encodes all algebra constraints automatically!
pop = polyopt(f, registry)

solver_config = SolverConfig(optimizer=SILENT_MOSEK,  # solver backend
    order=1                    # relaxation order
)
result = cs_nctssos(pop, solver_config)
@show result.objective  # upper bound on the maximal quantum violation

# The resulting upper bound is very close to the theoretical exact value $2\sqrt{2} \approx 2.8284271247461903$ (accurate up to 7 decimals!!).

# The `create_unipotent_variables` function creates variables with the UnipotentAlgebra type,
# which automatically encodes that these operators square to the identity (U^2 = I).
# This is appropriate for Pauli operators and other self-inverse observables.
#
# The registry groups variables by their labels ("x" and "y"), and variables with different
# labels automatically commute with each other. This replaces the old `comm_gps` parameter.
#
# In the solver configuration, we use [Mosek](https://www.mosek.com/) as the SDP solver backend.

# ### $I_{3322}$ inequality

# The $I_{3322}$ inequality is a more complex inequality that involves two parties, each measuring three observables. Let $A_1, A_2, A_3$ be the projective (squared to itself) observables measured by Alice and $B_1, B_2, B_3$ be the projective observables measured by Bob. We define the objective function as [pal2010maximal](@cite):
#
# ```math
# f(A_1, A_2, A_3, B_1, B_2, B_3) = \langle A_1(B_1+B_2+B_3) \rangle + \langle A_2(B_1+B_2-B_3) \rangle\\
# + \langle A_3(B_1-B_2) \rangle
# - \langle A_1 \rangle - 2\langle B_1 \rangle - \langle B_2 \rangle.
# ```
#
# In classical mechanics, the inequality $f(A_1, A_2, A_3, B_1, B_2, B_3) \leq 0$ must be satisfied. However, quantum mechanics can violate this inequality up to the value $0.25$. This violation demonstrates that quantum mechanics cannot be described by any local hidden variable theory.

# An upper bound on the maximal quantum violation of the $I_{3322}$ inequality can be computed using the following code:

using NCTSSoS, MosekTools

# Create projector variables (operators that square to themselves: P^2 = P)
registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]

pop = polyopt(-f, registry)

solver_config = SolverConfig(optimizer=SILENT_MOSEK, order=2)

result = cs_nctssos(pop, solver_config)
@show result.objective

# The `create_projector_variables` function creates variables with the ProjectorAlgebra type,
# which automatically encodes that these operators are idempotent (P^2 = P).
# This is appropriate for projection operators like $|0\rangle\langle 0|$ and $|1\rangle\langle 1|$.
#
# The resulting upper bound is close to the theoretically exact value $0.25$. By increasing the relaxation order, this upper bound could be further improved.

# #### Reducing the SDP Size by exploiting sparsity

# To reach the theoretically exact value $0.25$, one may increase the relaxation order [magronSparsePolynomialOptimization2023](@cite).

using NCTSSoS, MosekTools

registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]

pop = polyopt(-f, registry)

solver_config = SolverConfig(optimizer=SILENT_MOSEK, order=3)

@time result = cs_nctssos(pop, solver_config)
@show result.objective
# Indeed, by increasing the relaxation order to 3, we tighten the upper bound toward the theoretical value $0.25$ [pal2010maximal](@cite).
#
# However, keep increasing the order leads to large-scale SDPs that are computationally expensive. To reduce the SDP size, we may exploit the sparsity of the problem [magronSparsePolynomialOptimization2023](@cite). There are two types of sparsities:
#
# 1. **Correlative Sparsity**: exploiting the fact that few variable products appear in the objective function. Therefore, we could break down the objective function into smaller parts, each involving fewer variables. This reduces the matrix size and the number of constraints of the SDP, making it more tractable.
#
# 2. **Term Sparsity**: exploiting the fact that few monomials appear in the objective function. By identifying and removing unnecessary monomials, we can further reduce the matrix size and the number of constraints of the SDP.

# To take advantage of these sparsity patterns:

using NCTSSoS, MosekTools

registry, (x, y) = create_projector_variables([("x", 1:3), ("y", 1:3)])
f = 1.0 * x[1] * (y[1] + y[2] + y[3]) + x[2] * (y[1] + y[2] - y[3]) +
    x[3] * (y[1] - y[2]) - x[1] - 2 * y[1] - y[2]

pop = polyopt(-f, registry)

solver_config = SolverConfig(optimizer=SILENT_MOSEK, order=6, cs_algo=MF())

@time result = cs_nctssos(pop, solver_config)
@show result.objective

# Using almost half of the time, we are able to improve the $7$-th digit of the upper bound!

# ## Nonlinear Bell Inequalities

# Nonlinear Bell inequalities are extensions of the standard linear Bell inequalities. Instead of being linear combinations of expectation values, they involve polynomial functions of these expectation values. These inequalities arise naturally when considering more complex scenarios, such as multi-party settings or when the parties can perform sequences of measurements.
#
# The significance of nonlinear Bell inequalities in quantum information lies in their ability to detect non-locality in situations where linear inequalities might fail. They can provide tighter bounds on classical correlations and reveal quantum non-locality in a broader range of experimental setups. Furthermore, studying nonlinear Bell inequalities helps in understanding the structure of quantum correlations and the boundary between classical and quantum physics more deeply. They are also relevant in the context of quantum cryptography and communication complexity, where understanding the limits of classical and quantum correlations is crucial.

# ### Covariance Bell Inequality

# The covariance Bell inequality is a nonlinear Bell inequality that involves the covariance of measurements. It can be expressed as:
#
# ```math
# \text{Cov}(A, B) = \langle A B \rangle - \langle A \rangle \langle B \rangle,
# ```
#
# where $A$ and $B$ are observables measured by two parties. The covariance Bell inequality is nonlinear because it involves the product of expectation values of two observables.
#
# Let us define the objective function as:
# ```math
# f(A_1,A_2,A_3, B_1,B_2,B_3) = \text{Cov}(A_1, B_1) + \text{Cov}(A_1, B_2) + \text{Cov}(A_1,B_3)  + \\ \text{Cov}(A_2, B_1) + \text{Cov}(A_2, B_2) - \text{Cov}(A_2, B_3) + \text{Cov}(A_3, B_1) - \text{Cov}(A_3,B_2).
# ```
#
# It was shown that $f(A_1,A_2,A_3,B_1,B_2,B_3) \leq \frac{9}{2}$ in classical models, while it attains the quantum violation $5$ with a maximally entangled state in a spatial quantum model [pozsgay2017Covariance](@cite).
#
# An *open question* is: what is the maximal quantum violation that the covariance Bell inequality can attain in spatial quantum models. We can tackle this question using state polynomial optimization [klep2024State](@cite).

using NCTSSoS, MosekTools

# Create unipotent variables for Alice's and Bob's observables
registry, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])

# Identity monomial for converting StatePolynomial to NCStatePolynomial
const ID = one(NormalMonomial{UnipotentAlgebra,UInt8})

# covariance function
cov(a, b) = 1.0 * ς(x[a] * y[b]) * ID - 1.0 * ς(x[a]) * ς(y[b]) * ID

# objective function
sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)

spop = polyopt(sp, registry)

solver_config = SolverConfig(
    optimizer=SILENT_MOSEK,          # solver backend
    order=2                             # relaxation order
)

result = cs_nctssos(spop, solver_config)
@show result.objective

# !!! note "Typing Unicodes"
#     You can type the unicode characters in the code by using `\varsigma` and pressing `Tab` to get the unicode character `ς`.
#
# The resulting upper bound is very close to the previously known best value $5$ (accurate up to 3 decimals!!).

# We can use sparsity to further improve the bound.

using NCTSSoS, MosekTools

registry, (x, y) = create_unipotent_variables([("x", 1:3), ("y", 1:3)])

# Identity monomial for converting StatePolynomial to NCStatePolynomial
const ID = one(NormalMonomial{UnipotentAlgebra,UInt8})

# covariance function
cov(a, b) = 1.0 * ς(x[a] * y[b]) * ID - 1.0 * ς(x[a]) * ς(y[b]) * ID

# objective function
sp = cov(1,1) + cov(1,2) + cov(1,3) + cov(2,1) + cov(2,2) - cov(2,3) + cov(3,1) - cov(3,2)

spop = polyopt(sp, registry)

solver_config = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=3,
    ts_algo=MF()
)

result = cs_nctssos(spop, solver_config)

result_higher = cs_nctssos_higher(spop, result, solver_config)
@show result_higher.objective

# This is accurate to high precision relative to the literature value $5$ [pozsgay2017Covariance](@cite).