# Independent representations of raw operator strings. No package multiplication,
# adjoint, normal ordering, or basis enumeration computes an expected value.
using LinearAlgebra

const ORACLE_PAULI = (Complex{Int}[1 0;0 1], Complex{Int}[0 1;1 0],
    Complex{Int}[0 -im;im 0], Complex{Int}[1 0;0 -1])

# Factorized 2x2 matrices keep the original random 192-site Pauli cases feasible.
function pauli_normal(raw)
    matrices = Dict{Int,Matrix{Complex{Int}}}()
    for l in raw
        site,axis = divrem(l-1,3)
        matrices[site] = get(matrices,site,ORACLE_PAULI[1])*ORACLE_PAULI[axis+2]
    end
    coefficient = 1+0im
    word = Int[]
    for site in sort!(collect(keys(matrices)))
        matrix = matrices[site]
        axis = findfirst(P->!iszero(LinearAlgebra.tr(P*matrix)),ORACLE_PAULI)
        coefficient *= Complex{Int}(LinearAlgebra.tr(ORACLE_PAULI[axis]*matrix)/2)
        axis != 1 && push!(word,3site+axis-1)
    end
    coefficient,word
end

# Jordan-Wigner action on every occupation column, with active modes compressed
# in their original order. Inactive modes are spectators, not truncated modes.
function fermi_matrix(raw,modes)
    dim = 1<<length(modes)
    matrix = zeros(Int,dim,dim)
    for col in 0:dim-1
        row,coefficient = col,1
        for l in reverse(raw)
            i = findfirst(==(abs(l)),modes)-1
            occupied = (row>>i)&1
            if occupied == (l<0 ? 1 : 0)
                coefficient = 0
                break
            end
            coefficient *= isodd(count_ones(row&((1<<i)-1))) ? -1 : 1
            row ⊻= 1<<i
        end
        matrix[row+1,col+1] = coefficient
    end
    matrix
end

# Creation is multiplication by x and annihilation is d/dx. This unnormalized
# occupation representation is exact and has no finite-Fock cutoff error.
function boson_column(raw,modes,column)
    row = collect(column)
    coefficient = big(1)
    for l in reverse(raw)
        i = findfirst(==(abs(l)),modes)
        if l>0
            coefficient *= row[i]
            iszero(coefficient) && return coefficient,Tuple(row)
            row[i] -= 1
        else
            row[i] += 1
        end
    end
    coefficient,Tuple(row)
end

struct OperatorOracle{A}
    algebra::A
    modes::Vector{Int}
    matrices::Dict{Int,Matrix{Float64}}
end
function OperatorOracle(algebra,raws)
    raw = reduce(vcat,raws;init=Int[])
    modes = sort!(unique!(abs.(raw)))
    matrices = Dict{Int,Matrix{Float64}}()
    if algebra isa Nx.SiteAlgebra
        sites = sort!(unique(l&255 for l in raw))
        rng = MersenneTwister(62019)
        id = Matrix{Float64}(I,3,3)
        for l in sort!(unique(raw))
            m = randn(rng,3,3)
            if algebra isa Nx.SiteAlgebra{Nx.Free}
                m = (m+m')/8
            else
                q = Matrix(qr(m).Q)
                u = q[:,1:(isodd(l>>8) ? 1 : 2)]
                m = u*u'
                algebra isa Nx.SiteAlgebra{Nx.Unipotent} && (m=2m-id)
            end
            matrices[l] = reduce(kron,[site==(l&255) ? m : id for site in sites])
        end
    end
    OperatorOracle(algebra,modes,matrices)
end

function operator_matrix(oracle::OperatorOracle,raw)
    oracle.algebra isa Nx.Fermionic && return fermi_matrix(raw,oracle.modes)
    if oracle.algebra isa Nx.Pauli
        sites=sort!(unique((l-1)÷3 for l in oracle.modes))
        locals=[copy(ORACLE_PAULI[1]) for _ in sites]
        for l in raw
            site,axis=divrem(l-1,3)
            locals[findfirst(==(site),sites)] *= ORACLE_PAULI[axis+2]
        end
        return reduce(kron,locals;init=ones(Complex{Int},1,1))
    elseif oracle.algebra isa Nx.Bosonic
        # Mixed-word tests have degree <=4. Five input occupations distinguish
        # their coefficient polynomials; retain every output up to occupation 8.
        @assert length(raw)<=4
        columns=collect(Iterators.product(ntuple(_->0:4,length(oracle.modes))...))
        matrix=zeros(BigInt,9^length(oracle.modes),length(columns))
        for (j,column) in enumerate(columns)
            c,row=boson_column(raw,oracle.modes,column)
            i=1+sum(row[k]*9^(k-1) for k in eachindex(row);init=0)
            matrix[i,j]=c
        end
        return matrix
    end
    dim = isempty(oracle.matrices) ? 1 : size(first(values(oracle.matrices)),1)
    foldl((m,l)->m*oracle.matrices[l],raw;init=Matrix{Float64}(I,dim,dim))
end

function operator_equal(oracle::OperatorOracle,left,right)
    difference = [left; [(-c,w) for (c,w) in right]]
    if oracle.algebra isa Nx.Pauli
        terms = Dict{Vector{Int},ComplexF64}()
        for (c,w) in difference
            phase,word = pauli_normal(w)
            terms[word] = get(terms,word,0im)+c*phase
        end
        return all(iszero,values(terms))
    elseif oracle.algebra isa Nx.Bosonic
        # Each normally ordered coefficient is a polynomial in occupations of
        # degree <= the input degree. Evaluate degree+1 points on every mode.
        degree = maximum((length(w) for (_,w) in difference);init=0)
        for column in Iterators.product(ntuple(_->0:degree,length(oracle.modes))...)
            terms = Dict{Tuple,Complex{BigInt}}()
            for (c,w) in difference
                factor,row = boson_column(w,oracle.modes,column)
                terms[row] = get(terms,row,0)+Complex{BigInt}(c)*factor
            end
            all(iszero,values(terms)) || return false
        end
        return true
    end
    zero_matrix = zero(operator_matrix(oracle,Int[]))
    actual = sum((c*operator_matrix(oracle,w) for (c,w) in left);init=zero_matrix)
    expected = sum((c*operator_matrix(oracle,w) for (c,w) in right);init=zero_matrix)
    isapprox(actual,expected;atol=1e-10,rtol=1e-10)
end

oracle_terms(s::TermDict) = [(c,raw) for (raw,c) in bridge_terms(s)]
raw_adjoint(algebra,raw) = algebra isa Union{Nx.Fermionic,Nx.Bosonic} ? -reverse(raw) : reverse(raw)
