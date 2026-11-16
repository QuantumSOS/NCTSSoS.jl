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

# Match row words by public operator names, independent of packed word width.
word_key(registry, word) = Tuple(NCTSSoS.operator_name(registry, l) for l in NCTSSoS.word_letters(word))

function completed_blocks(built, rows, registry)
    indices = Dict(word_key(registry, last(only(NCTSSoS.terms(p)))) => i for (i, p) in enumerate(rows))
    [sort([indices[word_key(registry, w)] for w in block]) for block in NCTSSoS.block_rows(built)]
end

function labeled_edges(edges, labels)
    [(labels[i], labels[j]) for (i, j) in edges]
end

function completed_edges(built, rows, registry)
    edges = Set{Tuple{Int,Int}}()
    for block in completed_blocks(built, rows, registry)
        for i in eachindex(block), j in (i+1):length(block)
            push!(edges, (block[i], block[j]))
        end
    end
    return sort!(collect(edges))
end

function draw_term_graph(
    labels,
    coords,
    solid_edges;
    dashed_edges=Tuple{Int,Int}[],
    title,
    markersize=46,
    fontsize=16,
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
