using NCTSSoS
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: get_basis, star, simplify
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

println("="^60)
println("UNDERSTANDING MOMENT MATRIX INDEXING")
println("="^60)

# The moment matrix M[i,j] should represent ⟨basis[i]†, basis[j]⟩
# which in quantum mechanics is ⟨ψ| basis[i]† basis[j] |ψ⟩

println("\nFor Pauli operators, σ† = σ (self-adjoint)")
println("So basis[i]† = basis[i] for single Pauli operators")
println()

# Check a few examples
println("Example: M[2,3] should be ⟨basis[2]†, basis[3]⟩")
println("  basis[2] = ", basis[2])
println("  basis[3] = ", basis[3])
println("  basis[2]† = ", star(basis[2]))
println("  basis[2]† * basis[3] = ", star(basis[2]) * basis[3])
println("  After simplification: ", simplify(star(basis[2]) * basis[3], sa))
println()

# Find X₁X₂ in basis
x1x2_idx = findfirst(b -> string(b) == "x₁¹x₂¹", basis)
println("X₁X₂ is at index ", x1x2_idx, " in basis")
println()

# So M[2,3] represents ⟨X₁X₂⟩
# And M[1, x1x2_idx] also represents ⟨X₁X₂⟩ (since basis[1] = 1)
println("Key insight:")
println("  M[2,3] represents ⟨X₁ * X₂⟩ = ⟨X₁X₂⟩")
println("  M[1,$x1x2_idx] represents ⟨1 * X₁X₂⟩ = ⟨X₁X₂⟩")
println("  These should be the SAME value!")
println()

# Now check with Yao
println("="^60)
println("COMPUTING WITH YAO")
println("="^60)

# Build Hamiltonian
H_yao = sum(
    put(N, i=>X) * put(N, mod1(i+1, N)=>X) +
    put(N, i=>Y) * put(N, mod1(i+1, N)=>Y) +
    put(N, i=>Z) * put(N, mod1(i+1, N)=>Z)
    for i in 1:N
)

# Get ground state
H_matrix = Matrix(mat(H_yao))
eigendata = eigen(Hermitian(H_matrix))
E0_exact = minimum(real.(eigendata.values))
ground_state_idx = argmin(real.(eigendata.values))
psi_exact = eigendata.vectors[:, ground_state_idx]

println("Ground state energy: ", E0_exact)

# Compute ⟨X₁X₂⟩ directly
x1x2_op = put(N, 1=>X) * put(N, 2=>X)
x1x2_direct = psi_exact' * Matrix(mat(x1x2_op)) * psi_exact
println("⟨X₁X₂⟩ (direct): ", real(x1x2_direct))

# Compute ⟨X₁ * X₂⟩ as ⟨X₁⟩⟨X₂⟩? No, that's not right for non-commuting operators
# Actually ⟨X₁ * X₂⟩ = ⟨X₁X₂⟩ for adjacent products
println()

println("="^60)
println("COMPUTING MANUAL MOMENT MATRIX ELEMENTS")
println("="^60)

# Helper to convert subscripts
function subscript_to_int(s::AbstractString)
    result = ""
    for c in s
        if c == '₀'
            result *= '0'
        elseif c == '₁'
            result *= '1'
        elseif c == '₂'
            result *= '2'
        elseif c == '₃'
            result *= '3'
        elseif c == '₄'
            result *= '4'
        elseif c == '₅'
            result *= '5'
        elseif c == '₆'
            result *= '6'
        elseif c == '₇'
            result *= '7'
        elseif c == '₈'
            result *= '8'
        elseif c == '₉'
            result *= '9'
        else
            continue
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
            error("Unknown operator type: $op_char")
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

# Compute H_manual[2,3]
i = 2  # X₁
j = 3  # X₂
mi_mono = basis[i]
mj_mono = basis[j]

mi_dag_op = monomial_to_yao_operator(star(mi_mono), N)
mj_op = monomial_to_yao_operator(mj_mono, N)

combined_op = mi_dag_op * mj_op
op_matrix = Matrix(mat(combined_op))

H_manual_2_3 = psi_exact' * op_matrix * psi_exact

println("H_manual[2,3] = ⟨basis[2]† * basis[3]⟩ = ⟨X₁ * X₂⟩")
println("  = ", H_manual_2_3)
println("  Direct ⟨X₁X₂⟩ = ", x1x2_direct)
println("  Match? ", isapprox(H_manual_2_3, x1x2_direct, atol=1e-10))
println()

# Compute H_manual[1, x1x2_idx]
i = 1  # identity
j = x1x2_idx  # X₁X₂
mi_mono = basis[i]
mj_mono = basis[j]

mi_dag_op = monomial_to_yao_operator(star(mi_mono), N)
mj_op = monomial_to_yao_operator(mj_mono, N)

combined_op = mi_dag_op * mj_op
op_matrix = Matrix(mat(combined_op))

H_manual_1_x1x2 = psi_exact' * op_matrix * psi_exact

println("H_manual[1,$x1x2_idx] = ⟨basis[1]† * basis[$x1x2_idx]⟩ = ⟨1 * X₁X₂⟩")
println("  = ", H_manual_1_x1x2)
println("  Direct ⟨X₁X₂⟩ = ", x1x2_direct)
println("  Match? ", isapprox(H_manual_1_x1x2, x1x2_direct, atol=1e-10))
println()

println("="^60)
println("CONCLUSION")
println("="^60)
println("Both H_manual[2,3] and H_manual[1,$x1x2_idx] compute ⟨X₁X₂⟩")
println("They should be equal: ", isapprox(H_manual_2_3, H_manual_1_x1x2, atol=1e-10))
println()
println("The moment matrix satisfies:")
println("  M[i,j] = ⟨ψ| basis[i]† * basis[j] |ψ⟩")
println()
println("This means different matrix elements can represent the same")
println("expectation value if basis[i]† * basis[j] simplifies to the")
println("same monomial.")
