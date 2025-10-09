using Test, NCTSSoS
using NCTSSoS.FastPolynomials
using LinearAlgebra
using NCTSSoS: get_basis

@testset "GNS Reconstruction Tests" begin
    
    @testset "Simple 2x2 Case" begin
        @ncpolyvar x y
        
        # Create a simple Hankel matrix for basis [1, x, y] (degree 1)
        # This represents a 2x2 matrix representation where x and y commute
        H = [1.0  0.0  0.0;   # <1,1>  <1,x>  <1,y>
             0.0  1.0  0.0;   # <x,1>  <x,x>  <x,y>  
             0.0  0.0  1.0]   # <y,1>  <y,x>  <y,y>
        
        matrices = reconstruct(H, [x, y], 1)
        
        @test length(matrices) == 2
        @test size(matrices[1]) == (3, 3)  # rank should be 3
        @test size(matrices[2]) == (3, 3)
        
        # Test that matrices are real
        @test all(isreal.(matrices[1]))
        @test all(isreal.(matrices[2]))
    end
    
    @testset "Rank Deficient Case" begin
        @ncpolyvar x
        
        # Create a rank-1 Hankel matrix for basis [1, x] (degree 1)  
        H = [1.0  0.5;   # Rank-1 case where second variable is 0.5 * first
             0.5  0.25]
        
        matrices = reconstruct(H, [x], 1)
        
        @test length(matrices) == 1
        @test size(matrices[1]) == (1, 1)  # rank should be 1
    end
    
    @testset "Zero Hankel Matrix" begin
        @ncpolyvar x
        
        # Zero Hankel matrix should throw an error
        H = [0.0 0.0; 0.0 0.0]
        
        @test_throws ArgumentError reconstruct(H, [x], 1)
    end
    
    @testset "Dimension Mismatch" begin
        @ncpolyvar x y
        
        # Hankel matrix size doesn't match basis size
        H = [1.0 0.0; 0.0 1.0]  # 2x2 matrix
        # But basis for degree 1 with [x,y] has 3 elements: [1, x, y]
        
        @test_throws ArgumentError reconstruct(H, [x, y], 1)
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
		# Example usage: coeffs_in_new_order = swap_matrix * coeffs_in_original_order

		swap_matrix
        
        # Simple identity-like Hankel matrix
        H = swap_matrix*[1.0000 0.5000 0.5001 1.0483 −0.5483 −0.5483 1.0484; 
			0.5000 1.0483 −0.5483 1.0627 −0.0144 −0.6090 0.0606;
			0.5001 −0.5483 1.0484 −0.0144 −0.5340 0.0606 0.9878;
			1.0483 1.0627 −0.0144 1.4622 −0.3995 −0.8006 0.7863;
			−0.5483 −0.0144 −0.5340 −0.3995 0.3852 0.1917 −0.7256; 
			−0.5483 −0.6090 0.0606 −0.8006 0.1917 0.4411 −0.3804;
			1.0484 0.0606 0.9878 0.7863 −0.7256 −0.3804 1.3682] *swap_matrix
        
        # Test localizing matrix for x
        K_x = NCTSSoS.construct_localizing_matrix(H, x, basis, 2)
        
        @test K_x == [0.5000 1.0483 −0.5483; 
					  1.0483 1.0627 −0.0144; 
					  −0.5483 −0.0144 -0.5340] 
        
        # Test localizing matrix for y  
        K_y = NCTSSoS.construct_localizing_matrix(H, y, basis, 2)
        
        @test K_y == [0.5001 -0.5483 1.0484; -0.5483 -0.6090 0.0606; 1.0484 0.0606 0.9878] 
    end
    
    @testset "Higher Degree Case" begin
        @ncpolyvar x
        
        # Create Hankel matrix for degree 2: basis [1, x, x^2]
        H = [1.0  0.0  0.0;   # <1,1>   <1,x>   <1,x²>
             0.0  1.0  0.0;   # <x,1>   <x,x>   <x,x²>
             0.0  0.0  1.0]   # <x²,1>  <x²,x>  <x²,x²>
        
        matrices = reconstruct(H, [x], 2)
        
        @test length(matrices) == 1
        @test size(matrices[1]) == (3, 3)  # rank should be 3
    end

	@testset "Example 2.7" begin
		@ncpolyvar x y 

        basis = get_basis([x,y],1)

        # sort(map(b -> x * b, basis_order1))
        # sort(map(b -> y * b, basis_order1))

        H = [1.0000 0.5000 0.5001 
             0.5000 1.0483 −0.5483 
             0.5001 −0.5483 1.0484]

		perm = [1, 3, 2]

		n = length(basis)
		swap_matrix = zeros(Float64, n, n)
		for i in 1:n
			swap_matrix[i, perm[i]] = 1.0
		end

        K1 = swap_matrix * [ 0.5000 1.0483 −0.5483;
            1.0483 1.0627 −0.0144;
            −0.5483 −0.0144 −0.5340] * swap_matrix

        K2 = [
            0.5001 −0.5483 1.0484;
            −0.5483 −0.6090 0.0606; 
            1.0484 0.0606 0.9878
        ]


		X_mat, Y_mat = reconstruct(H,[x,y],1;rtol=1e-4)
		X_mat, Y_mat = reconstruct(K1,[x,y],1;rtol=1e-4)
		X_mat, Y_mat = reconstruct(K2,[x,y],1;rtol=1e-4)

        X_mat

		X_mat
		Y_mat
		@test X_mat = [0.5019 −0.8931; −0.8931 0.1727]
		@test Y_mat = [0.4981 0.8939; 0.8939 0.0825]
	end
end