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
        "Manual"=>Any[
            "Polynomials"=>"manual/polynomials.md",
            "Monomials"=>"manual/monomials.md",
            "Polynomial Optimization"=>"manual/polynomial_optimization.md",
            "Sparsities"=>"manual/sparsities.md",
            "SDP Relaxation"=>"manual/sdp_relaxation.md",
            "Optimizers"=>"manual/optimizers.md"
        ],
        "Examples"=>Any[
            "Bell inequalities"=>"examples/generated/bell.md",
            "CHSH GNS Reconstruction"=>"examples/generated/chsh_gns_reconstruction.md",
            "Trace Polynomial"=>"examples/generated/trace_poly.md",
            "Newton Chip Method"=>"examples/generated/newton_chip_method.md",
            "Stabilization vs. Exactness"=>"examples/generated/sparsity_convergence.md",
            "Monoid Algebra Showcase"=>"examples/generated/monoid_algebras_showcase.md",
            "PBW Algebra Showcase"=>"examples/generated/pbw_algebras_showcase.md",
            "Mixed Systems (Tensor Products)"=>"examples/generated/mixed_algebras_tensor_products.md",
            "Ground State Energy"=>"examples/generated/ground_state_energy.md",
            "Fermionic Ground State"=>"examples/generated/fermionic_ground_state.md",
            "Kitaev Chain"=>"examples/generated/kitaev_chain.md",
            "Hubbard Model"=>"examples/generated/hubbard_model.md",
            "Bosonic Ground State"=>"examples/generated/bosonic_ground_state.md",
            "Pauli Algebra Interface"=>"examples/generated/pauli_algebra_interface.md",
            "Certifying Ground State Property"=>"examples/generated/certify_ground_state_property.md",
            "GNS Optimizer Extraction"=>"examples/generated/gns_optimizer_extraction.md",
            "GNS Construction Guide"=>"examples/generated/gns_construction_guide.md",
            "GNS Construction for Pauli Operators"=>"examples/generated/pauli_gns_construction.md",
            "Periodic V2RDM (H₄ Chain)"=>"examples/generated/periodic_v2rdm_h4.md",
            "H4-My-Attempt"=>"examples/generated/h4_my_attempt.md",
            "Periodic V2RDM: Why Correlative Sparsity Misses the k-Blocks"=>"examples/generated/periodic_v2rdm_correlative_sparsity.md",
            "Periodic V2RDM: ²D Symmetry Blocks vs. Term-Sparsity Blocks"=>"examples/generated/periodic_v2rdm_block_structure.md",
            "Periodic V2RDM: Building the Paper's PQG Relaxation"=>"examples/generated/periodic_v2rdm_pqg_construction.md",
            "Periodic V2RDM: PQG Basis, Trace Constraints, and COSMO"=>"examples/generated/periodic_v2rdm_pqg_cosmo.md",
            "Periodic V2RDM: D-Only Bridge, Spin Resolution, and Solver Probes"=>"examples/generated/periodic_v2rdm_native_jump.md"
        ],
        "References"=>"reference.md",
        "APIs"=>["User interface" => "apis/interface.md", "Polynomials" => "apis/polynomials.md", "Sparsities" => "apis/sparsities.md", "SDP Relaxation" => "apis/relaxations.md"]
    ],
    plugins=[bib],
    # modules=[NCTSSoS],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true", assets=String["assets/citations.css"], size_threshold=10^6),
)

deploydocs(; repo="github.com/QuantumSOS/NCTSSoS.jl.git", devbranch="main")
