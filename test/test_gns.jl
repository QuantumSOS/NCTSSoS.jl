using Test, NCTSSoS
using NCTSSoS.FastPolynomials
using NCTSSoS.FastPolynomials: neat_dot
using LinearAlgebra
using NCTSSoS: get_basis, neat_dot

@testset "GNS Construction" begin
    @testset "Hankel dictionary utilities" begin
        @ncpolyvar x

        basis = get_basis([x], 2)  # [1, x, x^2]
        H = [
            1.0  0.2  0.7;
            0.2  0.8  0.3;
            0.7  0.3  0.6
        ]

        dict = NCTSSoS.hankel_entries_dict(H, basis)

        key_x = neat_dot(basis[1], basis[2])
        key_x2_first = neat_dot(basis[1], basis[3])
        key_x2_second = neat_dot(basis[2], basis[2])

        @test dict[key_x] == H[1, 2]
        @test dict[key_x2_first] == H[1, 3]
        @test dict[key_x2_second] == H[1, 3]

        K = NCTSSoS.construct_localizing_matrix(dict, x, basis)

        @test size(K) == (3, 3)
        @test K[1, 1] == H[1, 2]
        @test K[1, 2] == H[1, 3]
        @test K[2, 1] == H[1, 3]
        @test K[2, 2] == H[2, 3]
    end
    
    @testset "Localizing Matrix Construction" begin
        @ncpolyvar x y

        # Test the localizing matrix construction directly
        basis = get_basis([x, y], 2)  # [1, x, y, xy, x^2, yx, y^2]

        # We want to reorder to: [1, x, y, x^2, xy, yx, y^2]
        # Original indices:      1   2   3   4    5    6    7
        # Target indices:        1   2   3   5    4    6    7
        perm = [1, 2, 3, 5, 4, 6, 7]

        n = length(basis)
        swap_matrix = zeros(Float64, n, n)
        for i in 1:n
            swap_matrix[i, perm[i]] = 1.0
        end

        H = swap_matrix * [
            1.0000 0.5000 0.5001 1.0483 −0.5483 −0.5483 1.0484;
            0.5000 1.0483 −0.5483 1.0627 −0.0144 −0.6090 0.0606;
            0.5001 −0.5483 1.0484 −0.0144 −0.5340 0.0606 0.9878;
            1.0483 1.0627 −0.0144 1.4622 −0.3995 −0.8006 0.7863;
            −0.5483 −0.0144 −0.5340 −0.3995 0.3852 0.1917 −0.7256;
            −0.5483 −0.6090 0.0606 −0.8006 0.1917 0.4411 −0.3804;
            1.0484 0.0606 0.9878 0.7863 −0.7256 −0.3804 1.3682
        ] * swap_matrix

        hankel_dict = NCTSSoS.hankel_entries_dict(H, basis)
        localizing_basis = filter(m -> degree(m) <= 1, basis)

        K_x = NCTSSoS.construct_localizing_matrix(hankel_dict, x, localizing_basis)

        @test K_x ≈ [
            0.5000 1.0483 −0.5483;
            1.0483 1.0627 −0.0144;
            −0.5483 −0.0144 -0.5340
        ] atol = 1e-4

        K_y = NCTSSoS.construct_localizing_matrix(hankel_dict, y, localizing_basis)

        @test K_y ≈ [
            0.5001 -0.5483 1.0484;
            -0.5483 -0.6090 0.0606;
            1.0484 0.0606 0.9878
        ] atol = 1e-4
    end

@testset "GNS Reconstruction Tests" begin

    @testset "Dimension Mismatch" begin
        @ncpolyvar x y

        # Hankel matrix size doesn't match basis size
        H = [1.0 0.0; 0.0 1.0]  # 2x2 matrix
        # But basis for degree 1 with [x,y] has 3 elements: [1, x, y]

        @test_throws ArgumentError reconstruct(H, [x, y], 1, 1, 2)
    end

    @testset "Example 2.7" begin
        @ncpolyvar x y

        basis = get_basis([x, y], 2)  # [1, x, y, xy, x^2, yx, y^2]

        perm = [1, 2, 3, 5, 4, 6, 7]

        n = length(basis)
        swap_matrix = zeros(Float64, n, n)
        for i in 1:n
            swap_matrix[i, perm[i]] = 1.0
        end

        H = swap_matrix * [
            1.00001 0.499907 0.500102 1.0483 −0.5483 −0.5483 1.0484;
            0.499907 1.0483 −0.548283 1.0627 −0.0144 −0.6090 0.0606;
            0.500102 −0.548283 1.04827 −0.0144 −0.5340 0.0606 0.9878;
            1.0483 1.0627 −0.0144 1.4622 −0.3995 −0.8006 0.7863;
            −0.5483 −0.0144 −0.5340 −0.3995 0.3852 0.1917 −0.7256;
            −0.5483 −0.6090 0.0606 −0.8006 0.1917 0.4411 −0.3804;
            1.0484 0.0606 0.9878 0.7863 −0.7256 −0.3804 1.3682
        ] * swap_matrix

        X_mat, Y_mat = reconstruct(H, [x, y], 2, 1, 2)

        display(X_mat)

        @test X_mat ≈ [
            0.1727 −0.8931;
            −0.8931 0.5019
        ] atol = 1e-3

        @test Y_mat ≈ [
            0.0825 0.8939;
            0.8939 0.4981
        ] atol = 1e-3
    end
end  # GNS Reconstruction Tests
end  # GNS Construction