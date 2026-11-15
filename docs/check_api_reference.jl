using NCTSSoS
import NCTSSoS.Next as Next

# Run with: julia --project=docs docs/check_api_reference.jl
# Every public name needs a docstring and exactly one reference home.
reference = joinpath(@__DIR__, "src", "apis")
rendered = Symbol[]
for file in sort(readdir(reference; join=true))
    endswith(file, ".md") || continue
    text = read(file, String)
    @assert !occursin("@autodocs", text) "Use explicit public @docs entries in $file"
    for block in eachmatch(r"```@docs\n(.*?)```"s, text)
        for line in split(strip(block.captures[1]), '\n')
            entry = match(r"^NCTSSoS\.Next(?:\.([A-Za-z_][A-Za-z_0-9!]*))?$", strip(line))
            @assert !isnothing(entry) "Not a Next reference entry: $line"
            name = isnothing(entry.captures[1]) ? :Next : Symbol(entry.captures[1])
            @assert Base.ispublic(Next, name) "Non-public reference entry: $name"
            @assert isdefined(Next, name) && Base.Docs.hasdoc(Next, name) "No docstring: $name"
            push!(rendered, name)
        end
    end
end
public_names = sort(filter(name -> Base.ispublic(Next, name), names(Next)))
no_docstring = filter(name -> !isdefined(Next, name) || !Base.Docs.hasdoc(Next, name), public_names)
missing = setdiff(public_names, rendered)
@assert allunique(rendered) "A public name has more than one reference home"
println("Public names: ", length(public_names))
println("Rendered: ", length(rendered))
println("No docstring (", length(no_docstring), "): ", join(no_docstring, ", "))
println("Missing documented names: ", missing)
@assert isempty(no_docstring) "Add docstrings for the public names listed above"
@assert isempty(missing) "Add the missing public docstrings to the reference"
println("PASS: every public name has a docstring and one reference home")
