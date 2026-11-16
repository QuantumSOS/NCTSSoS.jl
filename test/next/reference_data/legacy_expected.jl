using TOML
const LEGACY_EXPECTED = Dict(case["id"] => case["expected"]
    for file in ("next_legacy.toml","next_legacy_extra.toml")
    for case in TOML.parsefile(joinpath(@__DIR__,"../../data/expectations",file))["cases"])
