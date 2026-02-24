using Documenter
using NCTSSoS
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

makedocs(;
    sitename="NCTSSoS.jl",
    doctest=false,
    pages=Any[
        "Home"=>"index.md",
        "Quick Start"=>"quick_start.md",
        "Manual"=>Any["Polynomials"=>"manual/polynomials.md", "Monomials"=>"manual/monomials.md", "Polynomial Optimization"=>"manual/polynomial_optimization.md", "Sparsities"=>"manual/sparsities.md", "SDP Relaxation"=>"manual/sdp_relaxation.md", "Optimizers"=>"manual/optimizers.md"],
        "Examples"=>Any[
            "Bell inequalities"=>"examples/generated/bell.md",
            "Trace Polynomial"=>"examples/generated/trace_poly.md",
            "Monoid Algebra Showcase"=>"examples/generated/monoid_algebras_showcase.md",
            "PBW Algebra Showcase"=>"examples/generated/pbw_algebras_showcase.md",
            "Mixed Systems (Tensor Products)"=>"examples/generated/mixed_algebras_tensor_products.md",
            "Ground State Energy"=>"examples/generated/ground_state_energy.md",
            "Pauli Algebra Interface"=>"examples/generated/pauli_algebra_interface.md",
            "Certifying Ground State Property"=>"examples/generated/certify_ground_state_property.md",
            # "GNS Construction for Pauli Operators"=>"examples/generated/pauli_gns_construction.md"  # TODO: fix bug
        ],
        "References"=>"reference.md",
        "APIs"=>["User interface" => "apis/interface.md", "Polynomials" => "apis/polynomials.md", "Sparsities" => "apis/sparsities.md", "SDP Relaxation" => "apis/relaxations.md"]
    ],
    plugins=[bib],
    # modules=[NCTSSoS],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true", assets=String["assets/citations.css"], size_threshold=10^6),
)

deploydocs(; repo="github.com/QuantumSOS/NCTSSoS.jl.git", devbranch="main")
