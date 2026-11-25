# The checked-in TSV is the exact physical operator from Tencent cut-recovery.
# Its sites are 1-based, row-major on the 3×3 tile. Each key denotes a product
# of disjoint Pauli dots (XX + YY + ZZ), without spin-1/2 factors.
function _append_square_witness3x3(form, table, allowed, L, ::Type{W}) where {W}
    sites = [1+x+L*y for y in 0:2 for x in 0:2]
    dots = [(i,j) => ntuple(3) do a
        last(NC.mul(NC.pauli_letter(W,sites[i],a-1), NC.pauli_letter(W,sites[j],a-1)))
    end for i in 1:8 for j in i+1:9] |> Dict

    row = open(joinpath(@__DIR__,"square_witness3x3.tsv")) do io
        header = match(r"^# scale=(\d+) floor=(-?\d+) denominator=(\d+)$",readline(io))
        scale, floor, denominator = parse.(BigInt,header.captures)
        @assert scale == denominator
        row = Dict(1 => -floor//scale)
        for line in eachline(io)
            fields = parse.(Int64,split(line))
            key = fields[1:end-1]
            iszero(last(key)) && resize!(key,2)  # Single dot is padded with 0,0.
            edges = [dots[(key[i],key[i+1])] for i in 1:2:length(key)]
            weight = last(fields)//scale
            coordinates = Dict{Int,Rational{BigInt}}()
            for factors in Iterators.product(edges...)
                w = foldl((a,b)->last(NC.mul(a,b)),factors)
                terms = _pricing_existing_terms(table,w,true)
                isnothing(terms) && throw(ArgumentError("3×3 physical witness observable $key is not represented."))
                for (id,c) in terms
                    isreal(c) || error("Non-real Hermitian moment in 3×3 physical witness.")
                    coordinates[id] = get(coordinates,id,0) + real(c)
                end
            end
            for (id,c) in coordinates
                iszero(c) && continue
                (id == 1 || (id in allowed && 1 <= id-1 <= length(form.c))) ||
                    throw(ArgumentError("3×3 physical witness observable $key lacks constrained coordinates."))
                row[id] = get(row,id,0) + weight*c
            end
        end
        row
    end

    constant = pop!(row,1,0//big(1))
    ids = sort!([id for (id,c) in row if !iszero(c)])
    values = Float64[row[id] for id in ids]
    # Free coordinates are Pauli expectations with modulus ≤ 1 on physical
    # states. Widen the half-space for coefficient rounding, then round b up.
    err = sum((abs(row[id]-Rational{BigInt}(v)) for (id,v) in zip(ids,values));init=0//big(1))
    b = Float64(constant+err,RoundUp)
    extra = sparse(ones(Int,length(ids)),ids.-1,-values,1,length(form.c))
    cone = NC.DataCone(:psd,1,length(form.b)+1,1,maximum(c.label for c in form.cones)+1)
    return NC.MomentStandardForm(vcat(form.A,extra),[form.b;b],form.c,form.offset,[form.cones;cone])
end
