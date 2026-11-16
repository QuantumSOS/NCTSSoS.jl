# # [Certifying Ground State Properties](@id certify-property)
#
# Ground-state properties reveal phases and correlations in quantum many-body
# systems [wang2024Certifying](@cite). Exact diagonalization gives a reference
# for small systems, but its cost grows exponentially with the number of spins.
#
# A variational upper bound on the energy can instead restrict a moment
# relaxation. Minimize and maximize an observable over those moments to bound
# its ground-state expectation. Here we check the method against exact
# diagonalization on three spins.
#
# The model is the periodic transverse-field Ising chain, with spin operators
# ``S^a=\sigma^a/2``:
#
# ```math
# H=\frac{J}{4}\sum_{j=1}^{N}Z_jZ_{j+1}
#   +\frac{h}{2}\sum_{j=1}^{N}X_j,\qquad Z_{N+1}=Z_1.
# ```
#
# We bound the energy and the correlation ``\langle Z_1Z_2\rangle/4``.
# The SDP inequalities hold in exact arithmetic. This page uses floating-point
# solvers and checks the results at stated tolerances; it is not a proof with
# rigorous rounding-error bounds.

using Yao, LinearAlgebra
using NCTSSoS, Mosek

N = 3;
J = 1.0;
xs = collect(0.1:0.2:2.0);
ground_state_energy_upper_bounds = Float64[];
s1s2values = Float64[];

# ## Exact diagonalization reference
#
# Construct the Hamiltonian as an `8 × 8` matrix. Store the energy per site
# and the correlation in its lowest eigenvector. These are numerical exact-
# diagonalization references, not DMRG data.

for h in xs
    ham = sum(J / 4 * kron(N, i => Z, mod1(i + 1, N) => Z) for i in 1:N) +
          sum(h / 2 * put(N, i => X) for i in 1:N)
    F = eigen(Hermitian(Matrix(ham)))
    ground_state = F.vectors[:, 1]
    s1s2 = Matrix(kron(N, 1 => Z, 2 => Z)) / 4
    push!(ground_state_energy_upper_bounds, F.values[1] / N)
    push!(s1s2values, real(dot(ground_state, s1s2 * ground_state)))
end

# Check three independently reproduced reference values:

@assert isapprox(s1s2values[1], -0.081304860191; atol=1e-12)
@assert isapprox(s1s2values[5], -0.045516505972; atol=1e-12)
@assert isapprox(s1s2values[10], -0.026857181875; atol=1e-12)

# ## Energy lower bounds from moment relaxations
#
# `pauli_variables` applies the Pauli product rules and commutation between
# sites. No extra algebra constraints are needed. Use a dense order-2
# relaxation and the direct Mosek backend.

registry, (x, y, z) = pauli_variables(N);
solver = MosekBackend();
relaxation = Relaxation(order=2);
ground_state_energy_lower_bounds = Float64[];

for h in xs
    ham = sum(J / 4 * z[j] * z[mod1(j + 1, N)] + h / 2 * x[j] for j in 1:N)
    result = solve(Problem(ham), relaxation; solver)
    push!(ground_state_energy_lower_bounds, result.bound / N)
end

# Columns below are `h`, SDP energy per site, and exact energy per site.
# They agree within `1e-6` for this small model.

hcat(xs, ground_state_energy_lower_bounds, ground_state_energy_upper_bounds)

@assert maximum(abs.(ground_state_energy_lower_bounds - ground_state_energy_upper_bounds)) < 1e-6 #src

# ## Bound a ground-state correlation with one scalar condition
#
# If ``E_{\rm upper}\ge E_0``, every ground state satisfies
# ``\langle E_{\rm upper}-H\rangle\ge0``. Impose this on both the lower-
# and upper-bound problems with `ExpectationPositive(E_upper - ham)`.
# `solve` minimizes, so minimize `-A` for an upper bound on ``\langle A\rangle``.
#
# This is **one scalar expectation inequality**, not the operator inequality
# ``E_{\rm upper}I-H\succeq0``. `Positive(E_upper - ham)` would constrain
# the whole operator and would exclude the excited spectrum. It is not the
# right call for this energy condition.
#
# `ExpectationZero(p)` similarly means ``\langle p\rangle=0``. Both scalar
# constraint types require a Hermitian polynomial. They currently support
# plain `order` or `rows` relaxations only, without symmetry, grouping,
# term splitting or strengthening.

A = z[1] * z[2] / 4;
lower_bounds = Float64[];
upper_bounds = Float64[];

for (i, h) in enumerate(xs)
    ham = sum(J / 4 * z[j] * z[mod1(j + 1, N)] + h / 2 * x[j] for j in 1:N)
    E_upper = ground_state_energy_upper_bounds[i] * N
    constraints = (ExpectationPositive(E_upper - ham),)
    lower = solve(Problem(A; constraints), relaxation; solver)
    upper = solve(Problem(-A; constraints), relaxation; solver)
    push!(lower_bounds, lower.bound)
    push!(upper_bounds, -upper.bound)
end

# Columns are `h`, lower bound, exact correlation, and upper bound.
# Both endpoints agree with the reference within `1e-4`; the enclosure check
# allows `1e-6` for numerical error. No reference value is substituted into
# a computed bound.

hcat(xs, lower_bounds, s1s2values, upper_bounds)

@assert all(lower_bounds .<= s1s2values .+ 1e-6) #src
@assert all(upper_bounds .>= s1s2values .- 1e-6) #src
@assert maximum(abs.(lower_bounds - s1s2values)) < 1e-4 #src
@assert maximum(abs.(upper_bounds - s1s2values)) < 1e-4 #src

# Earlier versions of this page added the scalar energy block after the
# solver's linear cache was built. That constraint was omitted from the
# solve, so the intervals stayed near [-0.25, 0.25]. Those values are not
# a ground-state reference. The new intervals include the energy condition
# and are compared with exact diagonalization above.

# ## Visualize the correlation bounds
#
# The three curves nearly overlap at this scale. The printed table retains
# the interval endpoints so the small differences can be inspected.

using CairoMakie

f = Figure(size=(800, 600));
ax = Axis(f[1, 1],
    xlabel="Transverse field strength (h)",
    ylabel=L"\langle S^z_1 S^z_2 \rangle",
    title="Ground-state correlation bounds");
scatterlines!(ax, xs, upper_bounds, color=:red, label="Upper bound (SDP)", linewidth=2, markersize=8);
scatterlines!(ax, xs, lower_bounds, color=:green, label="Lower bound (SDP)", linewidth=2, markersize=8);
scatterlines!(ax, xs, s1s2values, color=:blue, label="Exact value (ED)", linewidth=2, markersize=8);
axislegend(ax, position=:rb);
ax.xgridvisible = true;
ax.ygridvisible = true;
ax.xgridwidth = 0.5;
ax.ygridwidth = 0.5;
ax.xgridcolor = (:gray, 0.2);
ax.ygridcolor = (:gray, 0.2);
f

# ## Summary
#
# - A valid energy upper bound restricts the possible ground-state moments.
# - Declare that restriction with `ExpectationPositive`, not an operator
#   localizing constraint or a custom solver helper.
# - Minimize an observable and its negative to get two interval endpoints.
# - Check observables separately from energy: an accurate energy bound alone
#   does not imply accurate correlations.
# - A looser variational energy upper bound is valid, but can give a wider
#   correlation interval. Degenerate ground states can also have different
#   correlations; the interval must cover all of them.
