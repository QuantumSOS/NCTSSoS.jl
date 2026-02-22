module ExamplesStamp

using SHA
using TOML

const FORMAT_VERSION = 1
const STAMP_REL_PATH = joinpath("docs", "examples_stamp.toml")
const SEP = UInt8[0]

repo_root() = normpath(@__DIR__, "..")
stamp_path(root::AbstractString = repo_root()) = joinpath(root, STAMP_REL_PATH)

function _git_ls_files(root::AbstractString, pathspecs::Vector{String})::Vector{String}
    cmd = `git -C $root ls-files -z -- $(pathspecs...)`
    raw = read(cmd)
    files = split(String(raw), '\0'; keepempty=false)
    sort!(files)
    return files
end

function _sha256_files(root::AbstractString, relpaths::Vector{String})::String
    ctx = SHA.SHA256_CTX()
    for rel in relpaths
        SHA.update!(ctx, codeunits(rel))
        SHA.update!(ctx, SEP)
        SHA.update!(ctx, read(joinpath(root, rel)))
        SHA.update!(ctx, SEP)
    end
    return bytes2hex(SHA.digest!(ctx))
end

function expected_stamp(root::AbstractString = repo_root())::Dict{String, Any}
    src_files = _git_ls_files(root, ["src"])
    literate_files = _git_ls_files(root, ["docs/src/examples/literate"])
    generator_files = _git_ls_files(root, ["docs/generate_examples.jl"])

    return Dict(
        "format_version" => FORMAT_VERSION,
        "src_hash" => _sha256_files(root, src_files),
        "literate_hash" => _sha256_files(root, literate_files),
        "generator_hash" => _sha256_files(root, generator_files),
    )
end

function write_stamp!(root::AbstractString = repo_root(); path::AbstractString = stamp_path(root))::Nothing
    data = expected_stamp(root)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
    return nothing
end

function verify_stamp(root::AbstractString = repo_root(); path::AbstractString = stamp_path(root))::Nothing
    if !isfile(path)
        println(stderr, "ERROR: Missing examples regen stamp: $(relpath(path, root))")
        println(stderr, "Fix: run `make examples` locally (requires Mosek), then commit:")
        println(stderr, "  - docs/src/examples/generated/*.md")
        println(stderr, "  - docs/examples_stamp.toml")
        exit(1)
    end

    expected = expected_stamp(root)
    actual = TOML.parsefile(path)

    mismatches = String[]
    for (k, v_expected) in expected
        v_actual = get(actual, k, nothing)
        if v_actual != v_expected
            push!(mismatches, k)
        end
    end

    if !isempty(mismatches)
        println(stderr, "ERROR: Stale examples regen stamp: $(relpath(path, root))")
        for k in sort(mismatches)
            v_actual = get(actual, k, "<missing>")
            println(stderr, "  $(k): expected $(expected[k]) got $(v_actual)")
        end
        println(stderr, "Fix: run `make examples` locally (requires Mosek), then commit:")
        println(stderr, "  - docs/src/examples/generated/*.md")
        println(stderr, "  - docs/examples_stamp.toml")
        exit(1)
    end

    return nothing
end

function main(args::Vector{String} = ARGS)::Nothing
    verify_stamp()
    return nothing
end

end # module ExamplesStamp

if abspath(PROGRAM_FILE) == @__FILE__
    ExamplesStamp.main()
end
