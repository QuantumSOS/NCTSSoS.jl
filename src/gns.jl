function reconstruct()

end

using LinearAlgebra

H = [1.0000 0.5000 0.5001;0.5000 1.0483 −0.5483;0.5001 −0.5483 1.0484]

U, S, Vt = svd(H)

U * Diagonal(S) * U' 

U

S

nsp = nullspace(H,atol=1e-4)
nsp

U * U'


# rank of H determines the dimension of X_i
X_i = sqrt(Diagonal(S))^-1 * transpose(U) * K_i * U * sqrt(Diagonal(S))^-1





# K_i = 
function localizing_matrix(hankel, cur_var, basis,deg)
	# cur_var is the variable to be localized
	# you will need to take out all entries in hankel matrix 
	# that contains basis[i] * cur_var * basis[j] <= deg
end