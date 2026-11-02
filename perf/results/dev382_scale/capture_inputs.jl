# Audit helper, run after timing: reuse the measured generator, without a solve.
# A negative repeat endpoint makes the harness's 0:REPEATS loop empty.
using SHA
root = normpath(joinpath(@__DIR__, "../../.."))
empty!(ARGS)
scratch = tempname()
append!(ARGS, ["antiunitary", "inhomogeneous", "22", "-1", scratch])
include(joinpath(root, "perf/dev382_scale.jl"))
rm(scratch)
problem = construct()
registry, ops = create_pauli_variables(1:N)
basis = pauli_contiguous_chain_basis(ops, 2; periodic=true)
@assert length(basis) == problem.basis_words == 265
function label(m)
    letters = fill('I', N)
    for idx in m.word
        letters[(idx-1) ÷ 3 + 1] = "XYZ"[(idx-1) % 3 + 1]
    end
    return String(letters)
end
h = problem.reduction.original.objective
pop = polyopt(h, registry)
rejection = try
    pauli_translation_invariant_moment_relaxation(pop, ops, 2;
        basis, sign_symmetry=false, reflection=false)
    error("Unexpected acceptance of nonuniform objective")
catch err
    @assert err isa ArgumentError
    message = sprint(showerror, err)
    @assert occursin("one-site translation-invariant objective", message)
    message
end
record = Dict("n"=>N, "seed"=>"not applicable: deterministic coefficients, no RNG",
    "generator"=>"perf/dev382_scale.jl:construct",
    "generator_sha256"=>bytes2hex(sha256(read(joinpath(root,"perf/dev382_scale.jl")))),
    "hamiltonian"=>Dict(label(m)=>real(c) for (c,m) in terms(h)),
    "basis"=>label.(basis), "specialized_rejection"=>rejection)
open(joinpath(@__DIR__,"input-inhomogeneous-22.toml"),"w") do io
    TOML.print(io,record;sorted=true)
end
println("Frozen $(length(record["hamiltonian"])) coefficients and $(length(basis)) basis words")
println(rejection)
