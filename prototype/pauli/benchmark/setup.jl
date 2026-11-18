using Pkg

length(ARGS) == 1 || error("usage: julia --startup-file=no setup.jl NEW_ENV_DIR")
env = abspath(only(ARGS))
ispath(env) && error("the benchmark environment directory must not exist")
mkpath(env)
repo = normpath(joinpath(@__DIR__, "../../.."))
reference = joinpath(env, "QMBCertify")
rev = "0b152dfcd9bcde3403aa4d967642ce5261c80712"
run(`git init -q $reference`)
run(`git -C $reference fetch --depth=1 https://github.com/wangjie212/QMBCertify $rev`)
run(`git -C $reference checkout --detach FETCH_HEAD`)
Pkg.activate(env)
Pkg.develop([PackageSpec(path=repo), PackageSpec(path=reference)])
Pkg.add(["JuMP", "Mosek", "MosekTools"])
Pkg.instantiate()
