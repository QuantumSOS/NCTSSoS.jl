using NCTSSoS
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_basis, _neat_dot3, expval, simplify, star
using JuMP
using MosekTools
using LinearAlgebra

# Setup
N = 4
sys = pauli_algebra(N)
X_nc, Y_nc, Z_nc = sys.variables
vars = [X_nc; Y_nc; Z_nc]
sa = sys.simplify_algo

# Get basis
basis_unsorted = get_basis(Polynomial{ComplexF64}, vars, 2, sa)
basis = sort(basis_unsorted)
n_basis = length(basis)

println("Basis size: ", n_basis)
println("First 10 basis elements: ", basis[1:10])

# Build Hamiltonian
ham = sum(
    one(ComplexF64) * X_nc[i] * X_nc[mod1(i+1, N)] +
    Y_nc[i] * Y_nc[mod1(i+1, N)] +
    Z_nc[i] * Z_nc[mod1(i+1, N)]
    for i in 1:N
)

println("\nHamiltonian: ", ham)

# Create and solve problem
pop = cpolyopt(ham, sys)
SOLVER = optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0)
solver_config = SolverConfig(optimizer=SOLVER, order=2)

println("\nSolving...")
result = cs_nctssos(pop, solver_config)
model = result.model

println("Objective value: ", objective_value(model))

# Extract dual
cons = all_constraints(model, include_variable_in_set_constraints=true)
println("\nNumber of constraints: ", length(cons))
println("First constraint type: ", typeof(cons[1]))

X_matrix = value.(dual(cons[1]))
println("X_matrix size: ", size(X_matrix))

# Convert to complex Hermitian
n_total = size(X_matrix, 1)
n_half = div(n_total, 2)

X1 = X_matrix[1:n_half, 1:n_half]
X2 = X_matrix[(n_half+1):end, (n_half+1):end]
X3 = X_matrix[(n_half+1):end, 1:n_half]

H_R = (X1 + X2) / 2
H_I = (X3 - X_matrix[1:n_half, (n_half+1):end]') / 2
H_nctssos = H_R + im * H_I

println("H_nctssos size: ", size(H_nctssos))

# Now let's understand what the basis should be
# by looking at what NCTSSoS actually computes
println("\n" * "="^60)
println("INVESTIGATING BASIS CONSTRUCTION")
println("="^60)

# The moment matrix is indexed by expval(basis[i]† * basis[j])
# Let's compute what these are
println("\nComputing expval(basis[i]† * basis[j]) for first few elements:")
for i in 1:5
    for j in 1:5
        b_i_star = star(basis[i])
        product = _neat_dot3(b_i_star, one(basis[1]), basis[j])
        simplified = simplify(expval(product), sa)
        println("  expval(basis[$i]† * basis[$j]) = expval($(b_i_star) * $(basis[j])) = ", simplified)
    end
end

# Check if this matches our basis
println("\nChecking if expval products give us back basis elements:")
for i in 1:min(10, n_basis)
    for j in 1:min(10, n_basis)
        b_i_star = star(basis[i])
        product = _neat_dot3(b_i_star, one(basis[1]), basis[j])
        simplified = simplify(expval(product), sa)

        # Is this in our basis?
        idx = findfirst(b -> b == simplified, basis)
        if isnothing(idx)
            println("  WARNING: expval(basis[$i]† * basis[$j]) = $simplified NOT in basis!")
        end
    end
end
