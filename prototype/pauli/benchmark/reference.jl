import QMBCertify

const REFERENCE_REV = "0b152dfcd9bcde3403aa4d967642ce5261c80712"

function check_reference_checkout()
    dir = pkgdir(QMBCertify)
    readchomp(`git -C $dir rev-parse HEAD`) == REFERENCE_REV || error("wrong QMBCertify revision")
    isempty(readchomp(`git -C $dir status --porcelain --untracked-files=no`)) ||
        error("QMBCertify checkout has tracked changes")
end

function reference_options(lattice, L, extra, strengthened)
    !strengthened || (lattice == "chain" && L >= 8) || error("rdm=8 requires a ring with at least 8 sites")
    coefficient = lattice == "chain" ? 3/4 : 3/2
    return (; H_supp=[[1, 4]], H_coe=[coefficient], lattice, extra, SU2_symmetry=false,
        lso=strengthened, lol=strengthened ? 4 : L, pso=0, rdm=strengthened ? 8 : false,
        energy=Float64[], writetofile=false,
        three_type=[1, 1], J2=0, correlation=false, Gram=false, QUIET=false)
end
