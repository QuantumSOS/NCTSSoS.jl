import Pkg; Pkg.activate(@__DIR__);

using Literate

# Process Literate.jl examples
const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples", "literate")
const GENERATED_DIR = joinpath(@__DIR__, "src", "examples", "generated")

# Create generated directory if it doesn't exist
isdir(GENERATED_DIR) || mkdir(GENERATED_DIR)

println("Generating markdown files from Literate.jl examples...")
println("Source directory: $EXAMPLES_DIR")
println("Output directory: $GENERATED_DIR")

# Generate fenced `julia` code blocks (NOT Documenter `@example` blocks), so the
# docs build won't try to execute examples (CI has no Mosek license).
config = Dict(
    "execute" => false,
    "flavor" => Literate.CommonMarkFlavor(),
)

# Convert Literate.jl scripts to markdown
for file in readdir(EXAMPLES_DIR)
    if endswith(file, ".jl")
        name = splitext(file)[1]
        input = joinpath(EXAMPLES_DIR, file)

        println("Processing: $file")

        # Convert Literate.jl script to markdown
        Literate.markdown(input, GENERATED_DIR; config=config, name=name)

        println("  Generated: $(name).md")
    end
end

println("\nMarkdown generation complete!")
