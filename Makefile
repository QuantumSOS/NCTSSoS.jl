JL = julia --project

default: init test

# =============================================================================
# Setup
# =============================================================================

init:
	$(JL) -e 'using Pkg; Pkg.precompile()'

init-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.develop(path="."); Pkg.precompile()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile()'

update-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.update(); Pkg.precompile()'

# =============================================================================
# Testing
# =============================================================================
# Test structure:
#   polynomials/  - Core algebra (no solver)
#   quality/      - Code quality checks
#   solvers/      - SDP solver integration
#   physics/      - Physics models (LOCAL_TESTING only)
#
# Solver config:
#   Default: COSMO (open-source)
#   LOCAL_TESTING=true: Mosek (required for physics tests)
#
# Run single file: julia --project -e 'include("test/solvers/moment.jl")'
# =============================================================================

# Full test suite with Mosek
test:
	LOCAL_TESTING=true $(JL) -e 'using Pkg; Pkg.test()'

# Quick test with COSMO (CI-compatible, skips physics tests)
test-quick:
	$(JL) -e 'using Pkg; Pkg.test()'

# Polynomial tests only (fast, no solver needed)
test-polynomials:
	$(JL) -e 'using NCTSSoS, Test; include("test/polynomials/runtests.jl")'

# =============================================================================
# Documentation
# =============================================================================

servedocs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs=["docs/src/assets","docs/src/generated"])'

examples:
	$(JL) docs/generate_examples.jl

# =============================================================================
# Cleanup
# =============================================================================

clean:
	rm -rf docs/build
	find . -name "*.cov" -delete

.PHONY: init init-docs update update-docs test test-quick test-polynomials servedocs examples clean
