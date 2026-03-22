module TestExpectations

using TOML
using NCTSSoS

export expectations_path,
       expectations_load,
       expectations_case,
       expectations_oracle

const TEST_DATA_DIR = joinpath(pkgdir(NCTSSoS), "test", "data")

expectations_path(parts::AbstractString...) = joinpath(TEST_DATA_DIR, parts...)

function expectations_load(relpath::AbstractString)
    path = expectations_path(relpath)
    return TOML.parsefile(path)
end

function expectations_case(data, id::AbstractString)
    haskey(data, "cases") || error("Missing key `cases` in expectations TOML.")
    for case in data["cases"]
        case["id"] == id && return case
    end
    error("Case id not found in expectations TOML: $(repr(id))")
end

function expectations_oracle(relpath::AbstractString, id::AbstractString)
    data = expectations_load(relpath)
    case = expectations_case(data, id)
    haskey(case, "expected") || error("Missing key `expected` for case $(repr(id)).")
    expected = case["expected"]

    haskey(expected, "objective") || error("Missing key `expected.objective` for case $(repr(id)).")
    oracle = (opt=expected["objective"],)

    if haskey(expected, "sides")
        oracle = merge(oracle, (sides=expected["sides"],))
    end

    if haskey(expected, "nuniq")
        oracle = merge(oracle, (nuniq=expected["nuniq"],))
    elseif haskey(expected, "n_unique_moment_matrix_elements")
        oracle = merge(oracle, (nuniq=expected["n_unique_moment_matrix_elements"],))
    end

    return oracle
end

end
