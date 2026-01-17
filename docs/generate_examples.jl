import Pkg; Pkg.activate(@__DIR__);

using Literate
using SHA

# Process Literate.jl examples
const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples", "literate")
const GENERATED_DIR = joinpath(@__DIR__, "src", "examples", "generated")
const CHECKSUM_FILE = joinpath(GENERATED_DIR, ".literate.sha256")

# Create generated directory if it doesn't exist
isdir(GENERATED_DIR) || mkdir(GENERATED_DIR)

println("Generating markdown files from Literate.jl examples...")
println("Source directory: $EXAMPLES_DIR")
println("Output directory: $GENERATED_DIR")

function literate_checksum(files)
    buffer = IOBuffer()
    for file in files
        write(buffer, file)
        write(buffer, UInt8(0))
        write(buffer, read(file))
        write(buffer, UInt8(0))
    end
    return bytes2hex(SHA.sha256(take!(buffer)))
end

# Generate fenced `julia` code blocks (NOT Documenter `@example` blocks), so the
# docs build won't try to execute examples (CI has no Mosek license).
config = Dict(
    "execute" => true,
    "flavor" => Literate.CommonMarkFlavor(),
)

# Convert Literate.jl scripts to markdown
example_files = sort(filter(endswith(".jl"), readdir(EXAMPLES_DIR; join=true)))
checksum = literate_checksum(example_files)
if isfile(CHECKSUM_FILE)
    previous = strip(read(CHECKSUM_FILE, String))
    if previous == checksum
        println("No changes detected in Literate examples; skipping regeneration.")
        exit(0)
    end
end

for input in example_files
    name = splitext(basename(input))[1]
    println("Processing: $(basename(input))")

    # Convert Literate.jl script to markdown
    Literate.markdown(input, GENERATED_DIR; config=config, name=name)

    println("  Generated: $(name).md")
end

write(CHECKSUM_FILE, checksum * "\n")

println("\nMarkdown generation complete!")
