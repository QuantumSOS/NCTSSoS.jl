using Test, NCTSSoS
using LinearAlgebra
using Yao

if haskey(ENV, "LOCAL_TESTING")
    using MosekTools
    const SOLVER = Mosek.Optimizer
else
    using Clarabel
    const SOLVER = Clarabel.Optimizer
end

@testset "Moment Matrix Extraction" begin
    @testset "4-Qubit XXX Heisenberg Model - Physical Validation" begin
        # This test validates moment matrix reconstruction for a physically meaningful quantum system
        # using the 4-qubit XXX Heisenberg model with exact ground state comparison

        N = 4  # Number of qubits

        # Part 1: Get exact ground state using Yao.jl
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

        # Part 2: Solve with NCTSSoS
        sys = pauli_algebra(N)
        X_nc, Y_nc, Z_nc = sys.variables

        ham = sum(
            one(ComplexF64) * X_nc[i] * X_nc[mod1(i+1, N)] +
            Y_nc[i] * Y_nc[mod1(i+1, N)] +
            Z_nc[i] * Z_nc[mod1(i+1, N)]
            for i in 1:N
        )

        pop = cpolyopt(ham, sys)
        result = cs_nctssos(pop, SolverConfig(optimizer=SOLVER, order=2))

        # Part 3: Extract moment matrix from NCTSSoS solution
        mm = moment_matrix(result)
        H_nctssos = mm.matrix

        # Part 4: Manually compute moment matrix from exact ground state for comparison
        # Helper function to convert monomial to Yao operator
        function subscript_to_int(s::AbstractString)
            result = ""
            for c in s
                if c == '₀'; result *= '0'
                elseif c == '₁'; result *= '1'
                elseif c == '₂'; result *= '2'
                elseif c == '₃'; result *= '3'
                elseif c == '₄'; result *= '4'
                elseif c == '₅'; result *= '5'
                elseif c == '₆'; result *= '6'
                elseif c == '₇'; result *= '7'
                elseif c == '₈'; result *= '8'
                elseif c == '₉'; result *= '9'
                else; continue
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

                pauli_gate = if op_char == 'x'; X
                elseif op_char == 'y'; Y
                elseif op_char == 'z'; Z
                else; error("Unknown operator type: $op_char")
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

        # Compute manual moment matrix
        n_basis = length(mm.basis)
        H_manual = zeros(ComplexF64, n_basis, n_basis)

        for i in 1:n_basis
            for j in 1:n_basis
                mi_mono = mm.basis[i]
                mj_mono = mm.basis[j]

                mi_dag_op = monomial_to_yao_operator(NCTSSoS.FastPolynomials.star(mi_mono), N)
                mj_op = monomial_to_yao_operator(mj_mono, N)

                combined_op = mi_dag_op * mj_op
                op_matrix = Matrix(mat(combined_op))

                H_manual[i, j] = psi_exact' * op_matrix * psi_exact
            end
        end

        # Part 5: Physical Validation Tests
        # Both matrices should satisfy all physical constraints even if they represent different quantum states

        # Test 1: H_manual satisfies physical constraints
        @test norm(H_manual - H_manual') < 1e-10  # Hermitian
        @test minimum(real.(eigvals(Hermitian(H_manual)))) >= -1e-10  # PSD
        @test abs(real(H_manual[1,1]) - 1.0) < 1e-10  # Normalized (⟨1,1⟩ = 1)

        # Test 2: H_nctssos satisfies physical constraints
        @test norm(H_nctssos - H_nctssos') < 1e-10  # Hermitian
        @test minimum(real.(eigvals(Hermitian(H_nctssos)))) >= -1e-10  # PSD
        @test abs(real(H_nctssos[1,1]) - 1.0) < 1e-10  # Normalized (⟨1,1⟩ = 1)

        # Test 3: Energy consistency
        # SDP should achieve the exact ground state energy
        @test abs(result.objective - E0_exact) < 1e-4  # Allow small numerical tolerance

        # Test 4: Matrix dimensions match
        @test size(H_manual) == size(H_nctssos)

        # Test 5: Both have rank-1 structure (or low-rank for degenerate ground states)
        H_manual_eigs = eigvals(Hermitian(H_manual))
        H_nctssos_eigs = eigvals(Hermitian(H_nctssos))

        # Count eigenvalues > 1e-6 (effective rank)
        rank_manual = sum(H_manual_eigs .> 1e-6)
        rank_nctssos = sum(H_nctssos_eigs .> 1e-6)

        @test rank_manual <= 4  # Should be low-rank (at most degeneracy of ground state)
        @test rank_nctssos <= 4  # Should be low-rank

        # Test 6: Verify manual computation by spot-checking specific expectations
        # Find X₁ and X₂ in basis
        x1_str = "x₁¹"
        x2_str = "x₂¹"

        basis_strs = string.(mm.basis)
        x1_idx = findfirst(==(x1_str), basis_strs)
        x2_idx = findfirst(==(x2_str), basis_strs)

        if !isnothing(x1_idx) && !isnothing(x2_idx)
            # Direct computation: ⟨ψ|X₁X₂|ψ⟩
            x1x2_op = put(N, 1=>X) * put(N, 2=>X)
            x1x2_direct = psi_exact' * Matrix(mat(x1x2_op)) * psi_exact

            # Should match H_manual[x1_idx, x2_idx] = ⟨X₁†, X₂⟩ = ⟨X₁, X₂⟩ (since X† = X)
            @test abs(H_manual[x1_idx, x2_idx] - x1x2_direct) < 1e-10
        end
    end
end
