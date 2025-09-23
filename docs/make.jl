using Documenter
using NCTSSoS
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

makedocs(;
    sitename="NCTSSoS.jl",
    pages=Any[
        "Home"=>"index.md",
        "Quick Start"=>"quick_start.md",
        "Manual"=>Any["Polynomials"=>"manual/polynomials.md", "Polynomial Optimization"=>"manual/polynomial_optimization.md", "Sparsities"=>"manual/sparsities.md", "SDP Relaxation"=>"manual/sdp_relaxation.md", "Optimizers"=>"manual/optimizers.md"],
        "Examples"=>Any["Bell inequalities"=>"examples/generated/bell.md", "Trace Polynomial"=>"examples/generated/trace_poly.md",
        "Ground State Energy"=>"examples/generated/ground_state_energy.md",
        "Certifying Ground State Property"=>"examples/generated/certify_ground_state_property.md"
            ],
        "References"=>"reference.md",
        "APIs"=>["User interface" => "apis/interface.md", "Polynomials" => "apis/polynomials.md", "Sparsities" => "apis/sparsities.md", "SDP Relaxation" => "apis/relaxations.md"]
    ],
    plugins=[bib],
    # modules=[NCTSSoS, NCTSSoS.FastPolynomials],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true", assets=String["assets/citations.css"], size_threshold=10^6),
)

deploydocs(; repo="github.com/QuantumSOS/NCTSSoS.jl.git", devbranch="main")
