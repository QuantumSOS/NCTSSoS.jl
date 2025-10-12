using Test, NCTSSoS
using NCTSSoS.FastPolynomials
using NCTSSoS: get_basis, neat_dot
using LinearAlgebra

@testset "Pauli Operator Reconstruction" begin
    @testset "2-Qubit |00⟩ State - Pauli X, Y, Z" begin
        # Declare non-commuting variables for Pauli operators on 2 sites
        # x[i], y[i], z[i] represent Pauli X, Y, Z on site i
        @ncpolyvar x[1:2] y[1:2] z[1:2]

        # Combine all operators into a single variable list
        vars_site1 = [x[1], y[1], z[1]]
        vars_site2 = [x[2], y[2], z[2]]
        all_vars = [vars_site1; vars_site2]

        # Test reconstruction for site 1 only (simpler case)
        @testset "Single Site - Site 1" begin
            d = 2
            basis = get_basis(vars_site1, d)

            # For generic normalized state |ψ⟩ = (0.6|0⟩ + 0.8i|1⟩)
            # This state is not an eigenstate of any Pauli operator
            # It has full support on the 2D Hilbert space
            # ⟨ψ|I|ψ⟩ = 1 (normalized)
            # ⟨ψ|X²|ψ⟩ = 1, ⟨ψ|Y²|ψ⟩ = 1, ⟨ψ|Z²|ψ⟩ = 1 (Pauli operators square to identity)
            # This gives a full-rank Hankel matrix for proper GNS reconstruction

            n = length(basis)
            H = zeros(ComplexF64, n, n)

            # Helper function to compute expectation value for |ψ⟩ state
            # Pauli matrices: X = [0 1; 1 0], Y = [0 -i; i 0], Z = [1 0; 0 -1]
            # |ψ⟩ = (0.6|0⟩ + 0.8i|1⟩)

            function expval_psi(mono::Monomial)
                if isone(mono)
                    return 1.0 + 0.0im
                end

                # Compute the matrix representation by multiplying Pauli matrices
                mat = Matrix{ComplexF64}(I, 2, 2)
                for (var, exp) in zip(mono.vars, mono.z)
                    pauli_mat = if var == x[1]
                        ComplexF64[0 1; 1 0]  # X
                    elseif var == y[1]
                        ComplexF64[0 -im; im 0]  # Y
                    elseif var == z[1]
                        ComplexF64[1 0; 0 -1]  # Z
                    else
                        error("Unknown variable: $var")
                    end

                    for _ in 1:exp
                        mat = mat * pauli_mat
                    end
                end

                # Compute ⟨ψ|mat|ψ⟩ for normalized state
                # Use a generic normalized state: |ψ⟩ = (0.6|0⟩ + 0.8i|1⟩)
                # Check: |0.6|² + |0.8i|² = 0.36 + 0.64 = 1 ✓
                ket = ComplexF64[0.6; 0.8im]
                return dot(ket, mat * ket)
            end


            # Fill the Hankel matrix
            # H[i,j] = ⟨basis[i]†, basis[j]⟩ = ⟨basis[i]† * basis[j]⟩
            for i in 1:n
                for j in 1:n
                    # For Hermitian operators, basis[i]† = basis[i]
                    mono_i = basis[i]
                    mono_j = basis[j]

                    # Compute the product monomial
                    product = neat_dot(mono_i, mono_j)

                    H[i, j] = expval_psi(product)
                end
            end

            display(H)

            # Test that H is Hermitian (for real quantum states)
            @test H ≈ H' atol=1e-10

            # Reconstruct the operators
            # Using hankel_deg = 1 to get 2x2 representations (for 2-level system)
            # H_deg=2: full basis [1, x, y, z, x², xy, xz, y², yz, z²]
            # hankel_deg=1: use degree-1 block [1, x, y, z] for rank determination
            X_recon, Y_recon, Z_recon = reconstruct(H, vars_site1, d, 1; rtol=1e-8)

            X_recon_red = map(a -> (real(a) < 1e-6 ? zero(a) : real(a)) + im * (imag(a) < 1e-6 ? zero(a) : imag(a)), X_recon)
            Y_recon_red = map(a -> (real(a) < 1e-6 ? zero(a) : real(a)) + im * (imag(a) < 1e-6 ? zero(a) : imag(a)), Y_recon)
            Z_recon_red = map(a -> (real(a) < 1e-6 ? zero(a) : real(a)) + im * (imag(a) < 1e-6 ? zero(a) : imag(a)), Z_recon)

            X_recon_red^2

            # Expected Pauli matrices in the computational basis
            # Note: The reconstructed matrices may be in a different basis
            # but should have the same eigenvalues
            pauli_X = ComplexF64[0 1; 1 0]
            pauli_Y = ComplexF64[0 -im; im 0]
            pauli_Z = ComplexF64[1 0; 0 -1]

            # Check eigenvalues (basis-independent property)
            @test sort(real(eigvals(X_recon))) ≈ sort(real(eigvals(pauli_X))) atol=1e-6
            @test sort(real(eigvals(Y_recon))) ≈ sort(real(eigvals(pauli_Y))) atol=1e-6
            @test sort(real(eigvals(Z_recon))) ≈ sort(real(eigvals(pauli_Z))) atol=1e-6

            # Check commutation relations: [X,Y] = 2iZ, etc.
            # In the reconstructed basis, these should still hold
            commutator_XY = X_recon * Y_recon - Y_recon * X_recon
            commutator_YZ = Y_recon * Z_recon - Z_recon * Y_recon
            commutator_ZX = Z_recon * X_recon - X_recon * Z_recon

            # Commutators should be proportional to each other (structure constants)
            # [X,Y] should be proportional to Z_recon, etc.
            # For Pauli matrices: [X,Y] = 2iZ
            @test rank(commutator_XY) > 0  # Non-zero commutator
            @test rank(commutator_YZ) > 0
            @test rank(commutator_ZX) > 0

            # Check anti-commutation: {X,Y} = XY + YX should have specific properties
            anticommutator_XY = X_recon * Y_recon + Y_recon * X_recon
            # For Pauli matrices, anti-commutators should be zero
            @test norm(anticommutator_XY) < 1e-6

            # Check that X², Y², Z² ≈ I (up to scalar)
            X2 = X_recon * X_recon
            Y2 = Y_recon * Y_recon
            Z2 = Z_recon * Z_recon

            # All should be proportional to identity
            @test norm(X2 - tr(X2)/2 * I(2)) < 1e-6
            @test norm(Y2 - tr(Y2)/2 * I(2)) < 1e-6
            @test norm(Z2 - tr(Z2)/2 * I(2)) < 1e-6
        end

        @testset "Two Sites - |ψψ⟩ Product State" begin
            # For 2-qubit system in |ψψ⟩ = |ψ⟩⊗|ψ⟩ state where |ψ⟩ = (0.6|0⟩ + 0.8i|1⟩)
            # This is more complex as we need tensor products
            # Focus on site 1 operators (site 2 is in |ψ⟩)

            # Generate basis for site 1 only, degree 1: [1, x₁, y₁, z₁]
            basis = get_basis(vars_site1, 1)
            n = length(basis)

            # Hankel matrix for |ψψ⟩ state
            # When acting with site 1 operators, site 2 is traced out
            # This gives the same as single-site case with |ψ⟩
            H = zeros(ComplexF64, n, n)

            function expval_psi_site1(mono::Monomial)
                if isone(mono)
                    return 1.0 + 0.0im
                end

                mat = Matrix{ComplexF64}(I, 2, 2)
                for (var, exp) in zip(mono.vars, mono.exps)
                    if var in vars_site1
                        pauli_mat = if var == x[1]
                            ComplexF64[0 1; 1 0]
                        elseif var == y[1]
                            ComplexF64[0 -im; im 0]
                        elseif var == z[1]
                            ComplexF64[1 0; 0 -1]
                        else
                            error("Unknown variable: $var")
                        end

                        for _ in 1:exp
                            mat = mat * pauli_mat
                        end
                    end
                end

                ket_psi = ComplexF64[0.6; 0.8im]
                return dot(ket_psi, mat * ket_psi)
            end

            for i in 1:n
                for j in 1:n
                    product = basis[i] * basis[j]
                    H[i, j] = expval_psi_site1(product)
                end
            end

            # Reconstruct site 1 operators
            X1_recon, Y1_recon, Z1_recon = reconstruct(H, vars_site1, 1, 1; rtol=1e-10)

            # Verify they satisfy Pauli algebra
            @test sort(real(eigvals(X1_recon))) ≈ [-1.0, 1.0] atol=1e-6
            @test sort(real(eigvals(Z1_recon))) ≈ [-1.0, 1.0] atol=1e-6

            # Anti-commutation
            @test norm(X1_recon * Y1_recon + Y1_recon * X1_recon) < 1e-6
            @test norm(Y1_recon * Z1_recon + Z1_recon * Y1_recon) < 1e-6
            @test norm(Z1_recon * X1_recon + X1_recon * Z1_recon) < 1e-6
        end
    end
end
