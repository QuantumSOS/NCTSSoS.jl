# Keep option displays at the layer that owns their constructors and defaults.
import ..NCTSSoS: MMD

_show_option_algorithm(io,algorithm) = show(IOContext(io,:module=>parentmodule(typeof(algorithm))),algorithm)
function _show_option_algorithm(io,algorithm::MMD)
    print(io,"MMD(")
    algorithm.delta==MMD().delta || print(io,"delta=",algorithm.delta)
    print(io,")")
end

function Base.show(io::IO,option::Correlative)
    print(io,"Correlative(")
    if !isequal(option.algorithm,Correlative().algorithm)
        print(io,"algorithm=")
        _show_option_algorithm(io,option.algorithm)
    end
    print(io,")")
end

function Base.show(io::IO,option::TermSplit)
    print(io,"TermSplit(")
    custom_algorithm=!isequal(option.algorithm,TermSplit().algorithm)
    if custom_algorithm
        print(io,"algorithm=")
        _show_option_algorithm(io,option.algorithm)
    end
    if option.level!=TermSplit().level
        custom_algorithm && print(io,", ")
        print(io,"level=",option.level)
    end
    print(io,")")
end

function Base.show(io::IO,option::PhysicsInformed)
    print(io,"PhysicsInformed(")
    option.level==PhysicsInformed().level || print(io,"level=",option.level)
    print(io,")")
end

function Base.show(io::IO,option::Contiguous)
    print(io,"Contiguous(")
    option.periodic==Contiguous().periodic || print(io,"periodic=",option.periodic)
    print(io,")")
end

function _show_option_sites(io,sites)
    if last(sites)-first(sites)+1==length(sites)
        show(io,first(sites):last(sites))
    else
        show(io,sites)
    end
end

function Base.show(io::IO,kind::DensityMatrix)
    print(io,"DensityMatrix(")
    _show_option_sites(io,kind.sites)
    print(io,")")
end

function Base.show(io::IO,kind::InvariantCone)
    print(io,"InvariantCone(")
    _show_option_sites(io,kind.sites)
    print(io,", ")
    show(io,kind.r)
    print(io,")")
end

Base.show(io::IO,kind::Union{Stationarity,StateOptimality}) = print(io,nameof(typeof(kind)),"()")

function Base.show(io::IO,kind::Sector)
    print(io,"Sector(")
    show(io,kind.operator)
    print(io,", ")
    show(io,kind.eigenvalue)
    print(io,")")
end

function Base.show(io::IO,kind::FermionicConditions)
    print(io,"FermionicConditions(")
    _show_option_sites(io,kind.modes)
    print(io,"; conditions=")
    show(io,Tuple(kind.conditions))
    print(io,")")
end

Base.show(io::IO,::MIME"text/plain",option::Union{Correlative,TermSplit,PhysicsInformed,Contiguous,StrengtheningKind}) = show(io,option)
