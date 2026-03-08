function print_check(label, ok; expected=nothing, actual=nothing, detail=nothing)
    println(ok ? "[PASS] $label" : "[FAIL] $label")
    expected === nothing || println("  expected: ", repr("text/plain", expected))
    actual === nothing || println("  actual:   ", repr("text/plain", actual))
    detail === nothing || println("  detail:   ", detail)
    ok || error("Check failed: $label")
    return ok
end

function check_equal(label, actual, expected)
    print_check(label, actual == expected; expected=expected, actual=actual)
    return nothing
end

function check_isapprox(label, actual, expected; kwargs...)
    detail = join(["$(key)=$(value)" for (key, value) in pairs(kwargs)], ", ")
    print_check(
        label,
        isapprox(actual, expected; kwargs...);
        expected=expected,
        actual=actual,
        detail=isempty(detail) ? nothing : detail,
    )
    return nothing
end

function report_comparison(label, ok; expected=nothing, actual=nothing, detail=nothing)
    println(ok ? "[MATCH] $label" : "[MISMATCH] $label")
    expected === nothing || println("  expected: ", repr("text/plain", expected))
    actual === nothing || println("  actual:   ", repr("text/plain", actual))
    detail === nothing || println("  detail:   ", detail)
    return ok
end

function compare_equal(label, actual, expected)
    report_comparison(label, actual == expected; expected=expected, actual=actual)
    return nothing
end

function compare_isapprox(label, actual, expected; kwargs...)
    detail = join(["$(key)=$(value)" for (key, value) in pairs(kwargs)], ", ")
    report_comparison(
        label,
        isapprox(actual, expected; kwargs...);
        expected=expected,
        actual=actual,
        detail=isempty(detail) ? nothing : detail,
    )
    return nothing
end

function edge_pairs(graph)
    sort([
        (
            min(NCTSSoS.Graphs.src(edge), NCTSSoS.Graphs.dst(edge)),
            max(NCTSSoS.Graphs.src(edge), NCTSSoS.Graphs.dst(edge)),
        )
        for edge in NCTSSoS.Graphs.edges(graph)
    ])
end

function labeled_edges(edges, labels)
    [(labels[i], labels[j]) for (i, j) in edges]
end

function chordal_completion_edges(graph, ts_algo)
    completed = deepcopy(graph)
    for block in NCTSSoS.clique_decomp(graph, ts_algo)
        NCTSSoS.add_clique!(completed, block)
    end
    return edge_pairs(completed)
end

function draw_term_graph(
    labels,
    coords,
    solid_edges;
    dashed_edges=Tuple{Int,Int}[],
    title,
    markersize=34,
    fontsize=18,
)
    fig = Figure(size=(520, 320))
    ax = Axis(fig[1, 1], title=title)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()

    xs = first.(coords)
    ys = last.(coords)

    for (i, j) in solid_edges
        lines!(ax, [xs[i], xs[j]], [ys[i], ys[j]]; color=:black, linewidth=2)
    end
    for (i, j) in dashed_edges
        lines!(ax, [xs[i], xs[j]], [ys[i], ys[j]]; color=:black, linewidth=2, linestyle=:dash)
    end

    scatter!(ax, xs, ys; color=:white, strokecolor=:black, strokewidth=2, markersize=markersize)
    text!(ax, xs, ys; text=labels, align=(:center, :center), fontsize=fontsize)

    xpad = 0.6
    ypad = 0.6
    xlims!(ax, minimum(xs) - xpad, maximum(xs) + xpad)
    ylims!(ax, minimum(ys) - ypad, maximum(ys) + ypad)

    return fig
end
