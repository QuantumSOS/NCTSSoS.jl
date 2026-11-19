using Pkg

length(ARGS) == 2 || error("usage: julia --startup-file=no setup.jl NEW_ENV_DIR TRANSFER_LOCKS_DIR")
env, locks = abspath.(ARGS)
ispath(env) && error("the benchmark environment directory must not exist")
repo = normpath(joinpath(@__DIR__, "../../.."))
VERSION == v"1.12.7" || error("use the pinned Julia 1.12.7")
for (case, rev) in (("chain", "0b152dfcd9bcde3403aa4d967642ce5261c80712"),
                    ("square", "9a3ae67ff2bb6da1c04b72318ecbfcec71db5495"))
    project = joinpath(env, case)
    mkpath(project)
    for name in ("Project.toml", "Manifest.toml")
        cp(joinpath(locks,case,name),joinpath(project,name))
    end
    symlink(repo,joinpath(project,"NCTSSoS"))
    symlink(joinpath(repo,"prototype/pauli"),joinpath(project,"PauliWorkflowRedesign"))
    reference = joinpath(project,"QMBCertify")
    run(`git init -q $reference`)
    run(`git -C $reference fetch --depth=1 https://github.com/wangjie212/QMBCertify $rev`)
    run(`git -C $reference checkout --detach $rev`)
    Pkg.activate(project)
    Pkg.instantiate()
    for name in ("Project.toml", "Manifest.toml")
        read(joinpath(locks,case,name)) == read(joinpath(project,name)) || error("Pkg changed the pinned $case lockfile")
    end
end
