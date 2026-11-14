using SparseArrays: SparseMatrixCSC, sparse

public MomentStandardForm

"A cone in consecutive coordinates. PSD/HPSD coordinates use the Frobenius-isometric svec convention."
struct DataCone
    kind::Symbol
    side::Int
    first::Int
    length::Int
    label::Int32
end

"min c'x + offset, subject to b - A*x in the product cone; x is free."
struct MomentStandardForm
    A::SparseMatrixCSC{Float64,Int}
    b::Vector{Float64}
    c::Vector{Float64}
    offset::Float64
    cones::Vector{DataCone}
end

"max c'z + offset, subject to A*z = b and z in the product cone."
struct SOSStandardForm
    A::SparseMatrixCSC{Float64,Int}
    b::Vector{Float64}
    c::Vector{Float64}
    offset::Float64
    cones::Vector{DataCone}
end

"Lower moment form to sparse real data. HPSD uses lower-column order, re then im off diagonal, both scaled by √2."
function moment_standard_form(p::MomentProgram{C}) where {C}
    validate(p)
    I,J,V = Int[],Int[],Float64[]
    b = Float64[]
    cones = DataCone[]
    function add_component(k,imaginary,scale)
        push!(b,0.0)
        row = length(b)
        for t in term_range(p.forms,k)
            coef = scale*(imaginary ? imag(p.forms.coefs[t]) : real(p.forms.coefs[t]))
            iszero(coef) && continue
            id = p.forms.ids[t]
            if id==1
                b[row] += coef
            else
                push!(I,row); push!(J,Int(id)-1); push!(V,-Float64(coef))
            end
        end
    end
    for block in p.blocks
        first = length(b)+1
        if block.kind in (PSD,HPSD)
            for j in 1:block.side,i in j:block.side
                k = block.first_form+tri_offset(block.side,i,j)
                scale = i==j ? 1.0 : sqrt(2.0)
                add_component(k,false,scale)
                block.kind==HPSD && i!=j && add_component(k,true,scale)
            end
        else
            for k in form_range(block)
                add_component(k,false,1.0)
                block.kind==ZERO && C <: Complex && add_component(k,true,1.0)
            end
        end
        kind = block.kind==ZERO ? :zero : block.kind==NONNEGATIVE ? :nonnegative : block.kind==SOC ? :soc : block.kind==PSD ? :psd : :hpsd
        push!(cones,DataCone(kind,block.side,first,length(b)-first+1,block.label))
    end
    c,offset = zeros(Float64,nvars(p)-1),0.0
    for (id,coef) in zip(p.objective_ids,p.objective_coefs)
        id==1 ? (offset+=coef) : (c[id-1]+=coef)
    end
    return MomentStandardForm(sparse(I,J,V,length(b),length(c)),b,c,offset,cones)
end

function sos_standard_form(p::MomentStandardForm)
    cones = [DataCone(k.kind==:zero ? :free : k.kind,k.side,k.first,k.length,k.label) for k in p.cones]
    return SOSStandardForm(sparse(transpose(p.A)),-p.c,-p.b,p.offset,cones)
end
sos_standard_form(p::MomentProgram) = sos_standard_form(moment_standard_form(p))

function _jump_cones!(model,coordinates,cones)
    for cone in cones
        first,last = cone.first,cone.first+cone.length-1
        if cone.kind==:free
            continue
        elseif cone.kind==:zero
            for i in first:last; @constraint(model,coordinates[i]==0); end
        elseif cone.kind==:nonnegative
            for i in first:last; @constraint(model,coordinates[i]>=0); end
        elseif cone.kind==:soc
            @constraint(model,coordinates[first:last] in JuMP.SecondOrderCone())
        else
            n = cone.side
            if cone.kind==:psd
                M = Matrix{JuMP.AffExpr}(undef,n,n)
                k = first
                for j in 1:n,i in j:n
                    M[i,j] = coordinates[k]/(i==j ? 1.0 : sqrt(2.0))
                    M[j,i] = M[i,j]
                    k += 1
                end
                @constraint(model,LinearAlgebra.Symmetric(M) in JuMP.PSDCone())
            elseif cone.kind==:hpsd
                M = Matrix{JuMP.GenericAffExpr{ComplexF64,JuMP.VariableRef}}(undef,n,n)
                k = first
                for j in 1:n,i in j:n
                    if i==j
                        M[i,j] = coordinates[k]; k+=1
                    else
                        M[i,j] = (coordinates[k]+im*coordinates[k+1])/sqrt(2.0); k+=2
                        M[j,i] = conj(M[i,j])
                    end
                end
                @constraint(model,LinearAlgebra.Hermitian(M,:L) in JuMP.HermitianPSDCone())
            else
                throw(ArgumentError("Unknown standard-form cone $(cone.kind)."))
            end
        end
    end
    return model
end

function to_jump(p::MomentStandardForm)
    model = JuMP.Model()
    @variable(model,x[1:length(p.c)])
    coordinates = p.b-p.A*x
    _jump_cones!(model,coordinates,p.cones)
    @objective(model,Min,p.offset+sum(p.c[i]*x[i] for i in eachindex(p.c)))
    return model,x
end
function to_jump(p::SOSStandardForm)
    model = JuMP.Model()
    @variable(model,z[1:length(p.c)])
    @constraint(model,p.A*z .== p.b)
    _jump_cones!(model,z,p.cones)
    @objective(model,Max,p.offset+sum(p.c[i]*z[i] for i in eachindex(p.c)))
    return model,z
end
