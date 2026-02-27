module TestExpectations

using JSON3
using NCTSSoS

export expectations_path,
       expectations_load,
       expectations_case,
       expectations_oracle

const TEST_DATA_DIR = joinpath(pkgdir(NCTSSoS), "test", "data")

expectations_path(parts::AbstractString...) = joinpath(TEST_DATA_DIR, parts...)

function expectations_load(relpath::AbstractString)
    path = expectations_path(relpath)
    return JSON3.read(read(path, String))
end

function expectations_case(data, id::AbstractString)
    haskey(data, "cases") || error("Missing key `cases` in expectations JSON.")
    for case in data["cases"]
        String(case["id"]) == id && return case
    end
    error("Case id not found in expectations JSON: $(repr(id))")
end

function expectations_oracle(relpath::AbstractString, id::AbstractString)
    data = expectations_load(relpath)
    case = expectations_case(data, id)
    haskey(case, "expected") || error("Missing key `expected` for case $(repr(id)).")
    expected = case["expected"]

    haskey(expected, "objective") || error("Missing key `expected.objective` for case $(repr(id)).")
    oracle = (opt=Float64(expected["objective"]),)

    if haskey(expected, "sides")
        oracle = merge(oracle, (sides=[Int(x) for x in expected["sides"]],))
    end

    if haskey(expected, "nuniq")
        oracle = merge(oracle, (nuniq=Int(expected["nuniq"]),))
    elseif haskey(expected, "n_unique_moment_matrix_elements")
        oracle = merge(oracle, (nuniq=Int(expected["n_unique_moment_matrix_elements"]),))
    end

    return oracle
end

end

