_check_label(a::FunctionalAlgebra,label) = _check_label(a.base,label)
function word(r::Registry{<:FunctionalAlgebra},name::Symbol,::Type{W}) where {W<:FunctionalWord}
    site,label=r[name]
    site_letter(W,site,label)
end
function _model_word_type(r::Registry{<:FunctionalAlgebra})
    b=max(1,ndigits(r.nsites;base=2))+ndigits(r.algebra.base.nletters;base=2)
    word_type(r,128÷b)
end

"""
    state_symbol(p)

Symbolic expectation in the `Next.RealState()` chosen on the variable registry.
Linear in p; closes its operator tail and leaves existing scalar factors intact.
The result remains a Polynomial usable directly in Problem. This constructs a
symbol, not a numerical expectation; read a solved outer moment with expectation.
"""
state_symbol(p::Polynomial) = _functional_symbol(RealState,p)
"Normalized symbolic trace in the registry's Next.RealTrace() domain; not an unnormalized matrix trace."
trace_symbol(p::Polynomial) = _functional_symbol(RealTrace,p)
function _functional_symbol(::Type{F},p::Polynomial{W,C}) where {F,W,C}
    p.registry.algebra isa FunctionalAlgebra{F} || throw(ArgumentError("$(F===RealState ? "state_symbol" : "trace_symbol") requires site_variables(...; functional=Next.$(nameof(F))()); use the matching functional from the same registry. Pauli, fermionic, bosonic and mixed domains are unsupported."))
    out=PolynomialBuffer{W,C}(length(p.data))
    for (c,w) in terms(p); push_term!(out,c,first(moment_key(w))); end
    _polynomial(p.registry,out)
end

function _show_functional_segment(io,w::FunctionalWord{F,A,B,S},r,i,j) where {F,A,B,S}
    for k in i:j
        k>i && print(io," ")
        token=_ftoken(w,k); location=(_fsite(token,S),Int(token>>S))
        name=isnothing(r) ? nothing : findfirst(==(location),r.names)
        isnothing(name) ? print(io,"operator(site=",location[1],", label=",location[2],")") : print(io,name)
    end
end
function Base.show(io::IO,w::FunctionalWord{F},r::Union{Registry,Nothing}) where {F}
    isone(w) && return print(io,"1")
    i=_ftail(w)+1; printed=false
    while i<=degree(w)
        printed && print(io," ")
        print(io,F===RealState ? "⟨" : "tr(")
        j=_fstop(w,i); _show_functional_segment(io,w,r,i,j)
        print(io,F===RealState ? "⟩" : ")")
        printed=true; i=j+1
    end
    if _ftail(w)>0
        printed && print(io," ")
        _show_functional_segment(io,w,r,1,_ftail(w))
    end
end
