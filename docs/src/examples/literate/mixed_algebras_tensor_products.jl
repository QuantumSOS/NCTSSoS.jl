# # [Mixed Systems via Tensor Products](@id mixed-algebras-tensor-products)
#
# `NCTSSoS.jl` supports multiple algebra types (Pauli, fermionic, bosonic, …).
# For **mixed systems**, declare named tensor components in one registry with
# `Next.mixed_variables`. Products between components commute. Their
# expectations are joint moments; they do not factor into separate expectations.
#
# A mixed-system term has a coefficient and a word. `Next.terms` gives these
# pairs without requiring access to the polynomial's storage.

using NCTSSoS.Next, Mosek, LinearAlgebra
import NCTSSoS.Next as Next

# ## 1) Build factor-algebra monomials

registry, parts = Next.mixed_variables(spin=(Pauli(), 1), matter=(Fermionic(), 1));
σx, σy, σz = parts.spin;
a, ad = parts.matter;

# A Pauli product on the same site can introduce a complex phase:
# σx₁ σy₁ = i σz₁

spin_product = σx[1] * σy[1]
@assert spin_product == im*σz[1] #src

# A fermionic number operator is already normal-ordered:

number = ad[1] * a[1]
@assert number^2 == number #src

# ## 2) Form the tensor product
#
# Multiply the two polynomials. The registry keeps their component names,
# and multiplication applies the relations of each component automatically.

mixed_term = spin_product * number
@assert mixed_term == im*σz[1]*number #src
@assert σz[1]*number == number*σz[1] #src

# The phase is carried by the coefficient. Inspect the term with the public API:

# Pass the registry to `sprint` to display the word with symbolic names.

coefficient, word = only(Next.terms(mixed_term));
(coefficient, sprint(show, word; context=:registry => registry))
@assert coefficient == im #src

# Both factors are already canonical. Removing the phase leaves one unit term:

canonical = σz[1] * number;
[(c, sprint(show, w; context=:registry => registry)) for (c, w) in Next.terms(canonical)]
@assert first(only(Next.terms(canonical))) == 1 #src
@assert last(only(Next.terms(canonical))) == word #src

# ## 3) Check a finite mixed model with matrices
#
# Take `H = n/2 + Z + 2Xn`. In the empty fermion sector, the spin Hamiltonian
# is `Z`, with minimum `-1`. In the occupied sector it is `1/2 + Z + 2X`,
# with minimum `1/2-sqrt(5)`. Thus the exact ground energy is `1/2-sqrt(5)`.
# The matrices below act on spin ⊗ fermion and give a separate reference.

X = [0 1; 1 0];
Y = Complex{Int}[0 -im; im 0];
Z = [1 0; 0 -1];
I2 = Matrix{Int}(I, 2, 2);
A = [0 1; 0 0];
Nm = A' * A;
@assert kron(X*Y, Nm) == im*kron(Z, Nm) #src
Hmatrix = kron(I2, Nm)/2 + kron(Z, I2) + 2*kron(X, Nm);
exact = 0.5-sqrt(5);
eigvals(Hermitian(Hmatrix))
@assert isapprox(eigmin(Hermitian(Hmatrix)), exact; atol=1e-12) #src

# Use the same mixed registry for the objective and its moments. The order-2
# relaxation is checked against the exact reference with absolute tolerance `1e-6`.
# Moment accuracy is separate from bound accuracy. We request tighter Mosek
# tolerances to check both the energy and individual moments at `1e-6`.

H = number/2 + σz[1] + 2*σx[1]*number;
solver = MosekBackend(MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-11,
    MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-11, MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-11);
result = solve(Problem(H), Relaxation(order=2); solver);
println("SDP bound: ", result.bound)
println("Exact minimum: ", exact)
println("bound_kind: ", result.bound_kind)
println("valid_for: ", result.valid_for)
@assert isapprox(result.bound, exact; atol=1e-6) #src

# In the occupied ground state, `⟨n⟩ = 1` and `⟨Xn⟩ = -2/sqrt(5)`.
# Query the product as one polynomial; do not replace it with a product of moments.

(expectation(result, number), expectation(result, σx[1]*number))
@assert isapprox(expectation(result, number), 1; atol=1e-6) #src
@assert isapprox(expectation(result, σx[1]*number), -2/sqrt(5); atol=1e-6) #src

# Mixed registries currently support plain order-based relaxations, custom
# rows and `expectation`. Symmetry, grouping, term splitting, strengthening,
# state/trace components and `gns` are not yet supported. Put mutually
# anticommuting fermionic modes in one component; separate components commute.
