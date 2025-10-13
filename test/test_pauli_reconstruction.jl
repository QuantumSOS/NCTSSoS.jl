using Test, NCTSSoS
using NCTSSoS.FastPolynomials
using NCTSSoS: get_basis, neat_dot
using LinearAlgebra
using LinearAlgebra: tr

function expval_mixed(mono::Monomial)
    if isone(mono)
        return 1.0 + 0.0im
    end

    # Compute the matrix representation by multiplying Pauli matrices
    mat = Matrix{ComplexF64}(I, 2, 2)
    for (var, exp) in zip(mono.vars, mono.z)
        pauli_mat = if var == x
            ComplexF64[0 1; 1 0]  # X
        elseif var == y
            ComplexF64[0 -im; im 0]  # Y
        elseif var == z
            ComplexF64[1 0; 0 -1]  # Z
        else
            error("Unknown variable: $var")
        end

        for _ in 1:exp
            mat = mat * pauli_mat
        end
    end

    # Compute Tr(ρ * mat) for mixed state
    return tr(ρ * mat)
end

@testset "Pauli Operator Reconstruction" begin
    # Declare non-commuting variables for Pauli operators on 2 sites
    # x[i], y[i], z[i] represent Pauli X, Y, Z on site i
    @ncpolyvar x y z

    # Combine all operators into a single variable list
    vars_site1 = [x, y, z]

    # Test reconstruction for site 1 only (simpler case)
    @testset "Single Site - Site 1" begin
        d = 3
        basis = get_basis(vars_site1, d)

        ρ = zero_state * zero_state' * 0.5 + one_state * one_state' * 0.5

        n = length(basis)
        H = zeros(ComplexF64, n, n)

        # Fill the Hankel matrix
        # H[i,j] = ⟨basis[i]†, basis[j]⟩ = ⟨basis[i]† * basis[j]⟩
        for i in 1:n
            for j in 1:n
                # For Hermitian operators, basis[i]† = basis[i]
                mono_i = basis[i]
                mono_j = basis[j]

                # Compute the product monomial
                product = neat_dot(mono_i, mono_j)

                H[i, j] = expval_mixed(product)
            end
        end

        # Test that H is Hermitian (for real quantum states)
        @test H ≈ H' atol = 1e-10

        # Reconstruct the operators
        # Using hankel_deg = 1 to get 2x2 representations (for 2-level system)
        # H_deg=2: full basis [1, x, y, z, x², xy, xz, y², yz, z²]
        # hankel_deg=1: use degree-1 block [1, x, y, z] for rank determination

        X_recon, Y_recon, Z_recon = reconstruct(H, vars_site1, d, 1, 4)

        X = ComplexF64[0 1; 1 0]  # Pauli X
        Y = ComplexF64[0 -im; im 0]  # Pauli Y
        Z = ComplexF64[1 0; 0 -1]  #


        display(round.(X_recon))
        X_recon^2
        display(round.(Y_recon))
        Y_recon^2
        display(round.(Z_recon))
        Z_recon^2

        # Comprehensive Pauli Algebra Relations Tests
        # ============================================
        # The Pauli matrices must satisfy the following relations:
        # 1. X² = Y² = Z² = I (squares equal identity, up to scalar)
        # 2. [X,Y] = XY - YX = 2iZ (commutator relations)
        # 3. {X,Y} = XY + YX = 0 (anti-commutator relations)

        @testset "Pauli Algebra: Squares" begin
            # X², Y², Z² should all be proportional to identity
            X2 = X_recon * X_recon
            Y2 = Y_recon * Y_recon
            Z2 = Z_recon * Z_recon

            # All should be proportional to identity
            @test norm(X2 -  I(4)) < 1e-6
            @test norm(Y2 -  I(4)) < 1e-6
            @test norm(Z2 -  I(4)) < 1e-6
        end

        @testset "Pauli Algebra: Anti-commutators" begin
            # {X,Y} = XY + YX = 0
            # {Y,Z} = YZ + ZY = 0
            # {Z,X} = ZX + XZ = 0
            anticomm_XY = X_recon * Y_recon + Y_recon * X_recon
            anticomm_YZ = Y_recon * Z_recon + Z_recon * Y_recon
            anticomm_ZX = Z_recon * X_recon + X_recon * Z_recon

            @test norm(anticomm_XY) < 1e-6
            @test norm(anticomm_YZ) < 1e-6
            @test norm(anticomm_ZX) < 1e-6
        end

        @testset "Pauli Algebra: Commutators" begin
            # [X,Y] = XY - YX should be proportional to iZ
            # [Y,Z] = YZ - ZY should be proportional to iX
            # [Z,X] = ZX - XZ should be proportional to iY
            comm_XY = X_recon * Y_recon - Y_recon * X_recon
            comm_YZ = Y_recon * Z_recon - Z_recon * Y_recon
            comm_ZX = Z_recon * X_recon - X_recon * Z_recon

            @test norm(comm_XY - 2im * Z_recon) < 1e-5
            @test norm(comm_YZ - 2im * X_recon) < 1e-5
            @test norm(comm_ZX - 2im * Y_recon) < 1e-5
        end
    end
end
