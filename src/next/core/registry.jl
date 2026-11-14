"Names belong to the model layer; kernels receive words and never look up names."
struct Registry{A}
    algebra::A
    nsites::Int
    names::Dict{Symbol,Tuple{Int,Int}}
    function Registry(algebra::A, nsites::Integer, entries) where {A}
        nsites >= 0 || throw(ArgumentError("Site count must be nonnegative."))
        names = Dict{Symbol,Tuple{Int,Int}}()
        locations = Set{Tuple{Int,Int}}()
        for (name, (site, label)) in entries
            name isa Symbol || throw(ArgumentError("Registry names must be Symbols."))
            1 <= site <= nsites || throw(ArgumentError("$name has site $site outside 1:$nsites."))
            _check_label(algebra, label)
            haskey(names, name) && throw(ArgumentError("Duplicate operator name $name."))
            location = (Int(site), Int(label))
            location in locations && throw(ArgumentError("Two names refer to the same site and operator label $location."))
            names[name] = location
            push!(locations, location)
        end
        new{A}(algebra, Int(nsites), names)
    end
end
_check_label(::Pauli, label) = 0 <= label <= 2 || throw(ArgumentError("Pauli labels are 0 (X), 1 (Y), or 2 (Z)."))
_check_label(::Fermionic, label) = label in (-1,1) || throw(ArgumentError("Fermionic labels are -1 (creation) or 1 (annihilation)."))
_check_label(a::SiteAlgebra, label) = 1 <= label <= a.nletters || throw(ArgumentError("Site operator label must be in 1:$(a.nletters)."))
Base.length(r::Registry) = length(r.names)
Base.getindex(r::Registry, name::Symbol) = get(r.names, name) do
    throw(ArgumentError("Unknown operator name $name."))
end
word_type(r::Registry, max_letters::Integer) = word_type(r.algebra, r.nsites, max_letters)
function word(r::Registry{Pauli}, name::Symbol, ::Type{W}) where {W<:PauliWord}
    site, label = r[name]
    return pauli_letter(W, site, label)
end
function word(r::Registry{SiteAlgebra{A}}, name::Symbol, ::Type{W}) where {A,W<:SiteWord{A}}
    site, label = r[name]
    return site_letter(W, site, label)
end
function word(r::Registry{Fermionic}, name::Symbol, ::Type{W}) where {W<:FermiWord}
    site, label = r[name]
    return fermi_letter(W, site * label)
end
