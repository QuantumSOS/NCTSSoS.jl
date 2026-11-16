# Fresh-process cases from the architect's first_call_probe.jl, plus CHSH.
# Usage: julia --project=<env with NCTSSoS, SCS> perf/next/first_call.jl <case>
t_load = @elapsed begin
    using NCTSSoS, SCS
end
const N = NCTSSoS
case = only(ARGS)
backend = N.SCSBackend(verbose=0, eps_abs=1e-7, eps_rel=1e-7)
function inputs(case)
    if case == "xxz6"                       # close to the workload (4 sites there)
        n = 6; _, (x, y, z) = N.pauli_variables(n)
        h = sum(x[i] * x[mod1(i + 1, n)] + y[i] * y[mod1(i + 1, n)] + 0.7z[i] * z[mod1(i + 1, n)] for i in 1:n)
        return N.Problem(h), N.BlockRelaxation(order=2, grouping=N.Correlative(), symmetry=N.FamilySymmetry([N.SitePermutation([2:n; 1])]))
    elseif case == "xxz12split"             # other size, term split and a reflection: not in the workload
        n = 12; _, (x, y, z) = N.pauli_variables(n)
        h = sum(x[i] * x[mod1(i + 1, n)] + y[i] * y[mod1(i + 1, n)] + 0.7z[i] * z[mod1(i + 1, n)] for i in 1:n)
        return N.Problem(h), N.BlockRelaxation(order=2, grouping=N.Correlative(), split=N.TermSplit(),
            symmetry=N.FamilySymmetry([N.SitePermutation([2:n; 1]), N.SitePermutation([mod1(2 - i, n) for i in 1:n])]))
    elseif case == "xxz12front"             # the same problem in the front-door spelling (pull request 417)
        n = 12; _, (x, y, z) = N.pauli_variables(n)
        h = sum(x[i] * x[mod1(i + 1, n)] + y[i] * y[mod1(i + 1, n)] + 0.7z[i] * z[mod1(i + 1, n)] for i in 1:n)
        return N.Problem(h; symmetry=N.Symmetry(N.SitePermutation([2:n; 1]), N.SitePermutation([mod1(2 - i, n) for i in 1:n]))),
            N.Relaxation(order=2, grouping=N.Correlative(), split=N.TermSplit())
    elseif case == "hubbard"                # fermionic plain route: not in the workload
        _, (a, ad) = N.fermionic_variables(4)
        h = -(ad[1] * a[3] + ad[3] * a[1] + ad[2] * a[4] + ad[4] * a[2]) + 2 * (ad[1] * a[1] * ad[2] * a[2] + ad[3] * a[3] * ad[4] * a[4])
        return N.Problem(h), N.Relaxation(order=2)
    elseif case == "ball"                   # free algebra with a localizing constraint: not in the workload
        _, ops = N.site_variables(N.Free(), 1, 2); a, b = ops[1][1], ops[1][2]
        return N.Problem(a * b + b * a; constraints=(N.Positive(1 - a^2 - b^2),)), N.Relaxation(order=2)
    elseif case == "CHSH"
        _,(a,b)=N.site_variables(N.Unipotent(),2,2)
        return N.Problem(-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2])),N.Relaxation(order=1)
    end
    error("unknown case")
end
t_inputs = @elapsed args = inputs(case)
t_call = @elapsed r = N.solve(args...; solver=backend)
t_second = @elapsed N.solve(args...; solver=backend)
println("PROBE case=", case, " load=", round(t_load; digits=3), " inputs=", round(t_inputs; digits=3),
    " first_call=", round(t_call; digits=3), " second_call=", round(t_second; digits=4), " bound=", r.bound, " status=", r.status)
