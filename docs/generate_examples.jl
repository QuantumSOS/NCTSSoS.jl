import Pkg; Pkg.activate(@__DIR__);

using Literate
using SHA

# Process Literate.jl examples
const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples", "literate")
const GENERATED_DIR = joinpath(@__DIR__, "src", "examples", "generated")

# Marker prefix for checksum line
const CHECKSUM_MARKER = "<!-- nctssos-literate-source:"

# Create generated directory if it doesn't exist
isdir(GENERATED_DIR) || mkdir(GENERATED_DIR)

"""
    compute_source_hash(filepath::String, config::Dict) -> String

Compute SHA256 hash of source file contents combined with config dictionary.
Returns a hex string.
"""
function compute_source_hash(filepath::String, config::Dict)
    content = read(filepath, String)
    # Include config in hash to detect config changes
    # Sort by keys (strings) to ensure deterministic order
    config_str = string(sort(collect(config); by=first))
    combined = content * config_str
    return bytes2hex(sha256(combined))
end

"""
    extract_stored_hash(filepath::String) -> Union{String, Nothing}

Extract the stored hash from a generated markdown file.
Returns Nothing if no hash marker found or file doesn't exist.
"""
function extract_stored_hash(filepath::String)
    !isfile(filepath) && return nothing

    # Read first 5 lines to find marker
    lines = String[]
    open(filepath, "r") do io
        for _ in 1:5
            eof(io) && break
            push!(lines, readline(io))
        end
    end

    for line in lines
        if startswith(line, CHECKSUM_MARKER)
            # Extract hash from: <!-- nctssos-literate-source: filename.jl sha256: <hash> -->
            m = match(r"sha256:\s*([a-f0-9]+)", line)
            m !== nothing && return m.captures[1]
        end
    end
    return nothing
end

"""
    prepend_checksum_marker(filepath::String, source_name::String, hash::String)

Prepend the checksum marker to a generated markdown file.
"""
function prepend_checksum_marker(filepath::String, source_name::String, hash::String)
    content = read(filepath, String)
    marker_line = "$CHECKSUM_MARKER $source_name sha256: $hash -->\n\n"
    write(filepath, marker_line * content)
end

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
        output = joinpath(GENERATED_DIR, name * ".md")

        # Compute hash of source file + config
        current_hash = compute_source_hash(input, config)

        # Check if regeneration is needed
        stored_hash = extract_stored_hash(output)
        if stored_hash == current_hash
            println("Skipping: $file (unchanged, hash: $(current_hash[1:8])...)")
            continue
        end

        if stored_hash === nothing
            println("Processing: $file (new or missing hash)")
        else
            println("Processing: $file (changed, old: $(stored_hash[1:8])..., new: $(current_hash[1:8])...)")
        end

        # Convert Literate.jl script to markdown
        Literate.markdown(input, GENERATED_DIR; config=config, name=name)

        # Prepend checksum marker to generated file
        prepend_checksum_marker(output, file, current_hash)

        println("  Generated: $(name).md")
    end
end

println("\nMarkdown generation complete!")
