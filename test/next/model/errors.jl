# A short, executable transcript of five invalid user inputs.
using NCTSSoS
const Nx = NCTSSoS.Next
_,(X,Y,Z) = Nx.pauli_variables(3)
_,(other,_,_) = Nx.pauli_variables(3)
for (input,call) in (
    ("Problem(X[1]*Y[1])",()->Nx.Problem(X[1]*Y[1])),
    ("Problem(X[1]; constraints=[Positive(other[1])])",()->Nx.Problem(X[1];constraints=[Nx.Positive(other[1])])),
    ("Relaxation(order=-1)",()->Nx.Relaxation(order=-1)),
    ("Relaxation(rows=[X[1]])",()->Nx.Relaxation(rows=[X[1]])),
    ("inspect(Problem(X[1]*X[2]*X[3]), Relaxation(order=1))",()->Nx.inspect(Nx.Problem(X[1]*X[2]*X[3]),Nx.Relaxation(order=1))))
    println("> ",input)
    try
        call()
        error("Expected an input rejection")
    catch e
        e isa ArgumentError || rethrow()
        showerror(stdout,e); println()
    end
end
