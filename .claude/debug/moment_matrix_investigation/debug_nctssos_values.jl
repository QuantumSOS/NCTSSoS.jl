using NCTSSoS
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_basis
using JuMP
using MosekTools
using Yao
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

# Solve SDP
ham = sum(
    one(ComplexF64) * X_nc[i] * X_nc[mod1(i+1, N)] +
    Y_nc[i] * Y_nc[mod1(i+1, N)] +
    Z_nc[i] * Z_nc[mod1(i+1, N)]
    for i in 1:N
)

pop = cpolyopt(ham, sys)
SOLVER = optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0)
solver_config = SolverConfig(optimizer=SOLVER, order=2)
result = cs_nctssos(pop, solver_config)
model = result.model

# Extract H_nctssos
cons = all_constraints(model, include_variable_in_set_constraints=true)
X_matrix = value.(dual(cons[1]))

n_total = size(X_matrix, 1)
n_half = div(n_total, 2)

X1 = X_matrix[1:n_half, 1:n_half]
X2 = X_matrix[(n_half+1):end, (n_half+1):end]
X3 = X_matrix[(n_half+1):end, 1:n_half]

H_R = (X1 + X2) / 2
H_I = (X3 - X_matrix[1:n_half, (n_half+1):end]') / 2
H_nctssos = H_R + im * H_I

# Get exact ground state for comparison
H_yao = sum(
    put(N, i=>X) * put(N, mod1(i+1, N)=>X) +
    put(N, i=>Y) * put(N, mod1(i+1, N)=>Y) +
    put(N, i=>Z) * put(N, mod1(i+1, N)=>Z)
    for i in 1:N
)

H_matrix = Matrix(mat(H_yao))
eigendata = eigen(Hermitian(H_matrix))
E0_exact = minimum(real.(eigendata.values))
ground_state_idx = argmin(real.(eigendata.values))
psi_exact = eigendata.vectors[:, ground_state_idx]

# Compute direct
x1x2_op = put(N, 1=>X) * put(N, 2=>X)
x1x2_direct = psi_exact' * Matrix(mat(x1x2_op)) * psi_exact

# Find indices
x1x2_idx = findfirst(b -> string(b) == "x₁¹x₂¹", basis)

println("="^60)
println("COMPARING H_nctssos VALUES")
println("="^60)
println()
println("Direct Yao computation: ⟨X₁X₂⟩ = ", real(x1x2_direct))
println()
println("H_nctssos[2,3] (should be ⟨X₁ * X₂⟩ = ⟨X₁X₂⟩)")
println("  = ", real(H_nctssos[2,3]))
println("  Match? ", isapprox(H_nctssos[2,3], x1x2_direct, atol=1e-6))
println()
println("H_nctssos[1,$x1x2_idx] (should be ⟨1 * X₁X₂⟩ = ⟨X₁X₂⟩)")
println("  = ", real(H_nctssos[1,x1x2_idx]))
println("  Match? ", isapprox(H_nctssos[1,x1x2_idx], x1x2_direct, atol=1e-6))
println()

# Check if they're consistent with each other
println("Are H_nctssos[2,3] and H_nctssos[1,$x1x2_idx] equal?")
println("  ", isapprox(H_nctssos[2,3], H_nctssos[1,x1x2_idx], atol=1e-10))
println()

# Let's look for where the actual value appears
println("="^60)
println("SEARCHING FOR WHERE -0.666... APPEARS IN H_nctssos")
println("="^60)
println()

target = -0.666
for i in 1:length(basis)
    for j in 1:length(basis)
        if abs(real(H_nctssos[i,j]) - target) < 0.01
            println("H_nctssos[$i,$j] ($(basis[i]), $(basis[j])) = ", real(H_nctssos[i,j]))
        end
    end
end
println()

# Check normalized values
println("="^60)
println("CHECKING NORMALIZATION")
println("="^60)
println()
println("H_nctssos[1,1] (should be 1): ", real(H_nctssos[1,1]))
println("Trace of H_nctssos: ", real(tr(H_nctssos)))
println()

# Print first 5x5 block of both matrices
println("="^60)
println("COMPARING FIRST 5x5 BLOCK")
println("="^60)
println()
println("Basis indices 1-5: ", basis[1:5])
println()

# Compute H_manual for these indices
function subscript_to_int(s::AbstractString)
    result = ""
    for c in s
        if c in ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉']
            result *= string(Int(c) - Int('₀'))
        end
    end
    return parse(Int, result)
end

function monomial_to_yao_operator(mono, N)
    if isempty(mono.vars)
        return put(N, 1=>I2)
    end
    op = put(N, 1=>I2)
    first = true
    for (var, exp) in zip(mono.vars, mono.z)
        var_str = string(var.name)
        op_char = var_str[1]
        idx = subscript_to_int(var_str[2:end])
        pauli_gate = if op_char == 'x'
            X
        elseif op_char == 'y'
            Y
        elseif op_char == 'z'
            Z
        else
            error("Unknown")
        end
        for _ in 1:exp
            if first
                op = put(N, idx=>pauli_gate)
                first = false
            else
                op = op * put(N, idx=>pauli_gate)
            end
        end
    end
    return op
end

using NCTSSoS.FastPolynomials: star

H_manual_block = zeros(ComplexF64, 5, 5)
for i in 1:5
    for j in 1:5
        mi_mono = basis[i]
        mj_mono = basis[j]
        mi_dag_op = monomial_to_yao_operator(star(mi_mono), N)
        mj_op = monomial_to_yao_operator(mj_mono, N)
        combined_op = mi_dag_op * mj_op
        op_matrix = Matrix(mat(combined_op))
        H_manual_block[i, j] = psi_exact' * op_matrix * psi_exact
    end
end

println("H_manual block:")
display(round.(real.(H_manual_block), digits=4))
println()
println()

println("H_nctssos block:")
display(round.(real.(H_nctssos[1:5, 1:5]), digits=4))
println()
println()

println("Difference:")
display(round.(real.(H_manual_block - H_nctssos[1:5, 1:5]), digits=4))
