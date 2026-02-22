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

# Execute examples locally (requires Mosek license) and commit the generated
# markdown with outputs.  CI serves the pre-generated files directly.
config = Dict(
    "execute" => true,
    "flavor" => Literate.DocumenterFlavor(),
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

# Update regen stamp for CI gate (no MOSEK needed; hashes tracked inputs only).
include(joinpath(@__DIR__, "examples_stamp.jl"))
ExamplesStamp.write_stamp!()
