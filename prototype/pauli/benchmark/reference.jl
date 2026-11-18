import QMBCertify

const REFERENCE_REV = "0b152dfcd9bcde3403aa4d967642ce5261c80712"

function install_reference_probe()
    dir = pkgdir(QMBCertify)
    readchomp(`git -C $dir rev-parse HEAD`) == REFERENCE_REV || error("wrong QMBCertify revision")
    isempty(readchomp(`git -C $dir status --porcelain --untracked-files=no`)) ||
        error("QMBCertify checkout has tracked changes")
    path = joinpath(dir, "src", "bound_gsp.jl")
    source = first(split(read(path, String), "\nfunction arrange("; limit=2))
    for (old, new) in (
        "function GSB(" => "function benchmark_GSB(checkpoint, ",
        "                mons = filter_mons(mons, tsupp, h, H_coe, L, lattice=\"chain\")" => """
            mons = filter_mons(mons, tsupp, h, H_coe, L, lattice="chain")
            checkpoint(:stationarity, model, length(mons))""",
        "            posepsd8!(model, cons, tsupp, L, lattice=lattice)" => """
            checkpoint(:before_rdm, model, nothing)
            posepsd8!(model, cons, tsupp, L, lattice=lattice)
            checkpoint(:after_rdm, model, nothing)""",
        "    time = @elapsed begin\n    optimize!(model)\n    end" => """
            checkpoint(:built, model, nothing)
            time = @elapsed optimize!(model)
            checkpoint(:solved, model, time)""",
        "    objv = objective_value(model)" => """
            status == MOI.OPTIMAL && primal_status(model) == MOI.FEASIBLE_POINT &&
                dual_status(model) == MOI.FEASIBLE_POINT || return nothing, nothing
            objv = objective_value(model)""",
    )
        length(findall(old, source)) == 1 || error("reference probe anchor drift: $old")
        source = replace(source, old => new)
    end
    Base.include_string(QMBCertify, source, path)
end

function reference_options(lattice, L, extra, strengthened)
    !strengthened || (lattice == "chain" && L >= 8) || error("rdm=8 requires a ring with at least 8 sites")
    coefficient = lattice == "chain" ? 3/4 : 3/2
    return (; H_supp=[[1, 4]], H_coe=[coefficient], lattice, extra, SU2_symmetry=false,
        lso=strengthened, lol=strengthened ? 4 : L, pso=0, rdm=strengthened ? 8 : false,
        energy=Float64[], writetofile=false,
        three_type=[1, 1], J2=0, correlation=false, Gram=false, QUIET=true)
end

function reference_run(checkpoint, L, order, settings, options)
    return Base.invokelatest(QMBCertify.benchmark_GSB, checkpoint, options.H_supp, options.H_coe, L, order;
        options..., mosek_setting=settings)
end
