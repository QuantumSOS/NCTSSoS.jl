import QMBCertify

const REFERENCE_REVS = Dict("chain" => "0b152dfcd9bcde3403aa4d967642ce5261c80712",
    "square" => "9a3ae67ff2bb6da1c04b72318ecbfcec71db5495")

function check_reference_checkout(lattice)
    dir = pkgdir(QMBCertify)
    readchomp(`git -C $dir rev-parse HEAD`) == REFERENCE_REVS[lattice] || error("wrong QMBCertify revision")
    isempty(readchomp(`git -C $dir status --porcelain`)) || error("QMBCertify checkout is not clean")
end

function reference_options(lattice, L, extra, strengthened)
    chain = lattice == "chain"
    !strengthened || (chain ? L : L^2) >= 8 || error("rdm=8 requires at least 8 sites")
    coefficient = lattice == "chain" ? 3/4 : 3/2
    return (; H_supp=[[1, 4]], H_coe=[coefficient], lattice, extra, SU2_symmetry=false,
        lso=strengthened, lol=strengthened ? 4 : L, pso=strengthened && chain ? 2 : 0, rdm=strengthened ? 8 : false,
        energy=Float64[], writetofile=false,
        three_type=[1, 1], J2=0, correlation=false, Gram=false, QUIET=false)
end
