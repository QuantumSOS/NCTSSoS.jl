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

# Convert Literate.jl scripts to markdown
for file in readdir(EXAMPLES_DIR)
    if endswith(file, ".jl")
        name = splitext(file)[1]
        input = joinpath(EXAMPLES_DIR, file)
        output = joinpath(GENERATED_DIR, name * ".md")

        # Skip if output file already exists
        if isfile(output)
            println("Skipping: $file (already exists)")
            continue
        end

        println("Processing: $file")

        # Configure to generate regular Julia code blocks instead of @example blocks
        config = Dict("````julia" => "````", "execute" => true)

        # Convert Literate.jl script to markdown
        Literate.markdown(input, GENERATED_DIR; config=config, name=name)

        println("  Generated: $(name).md")
    end
end

println("\nMarkdown generation complete!")