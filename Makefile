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
#   physics/      - Physics models (--local only, requires Mosek)
#
# Flags:
#   --local        Use Mosek solver (required for physics tests)
#   --polynomials  Run polynomial algebra tests
#   --quality      Run code quality checks
#   --solvers      Run SDP solver tests
#   --physics      Run physics model tests
#
# Usage:
#   make test              - Full suite with Mosek
#   make test-ci           - CI suite (no physics, uses COSMO)
#   make test-polynomials  - Just polynomial algebra
#   make test-solvers      - Just SDP solver tests
#   make test-physics      - Just physics models (needs Mosek)
#   make test-quality      - Just code quality checks
#
# Direct Pkg.test usage:
#   Pkg.test("NCTSSoS"; test_args=["--polynomials"])
#   Pkg.test("NCTSSoS"; test_args=["--solvers", "--local"])
#   Pkg.test("NCTSSoS"; test_args=["--local"])  # Full suite with Mosek
# =============================================================================

# Full test suite with Mosek (includes physics)
test:
	$(JL) -e 'using Pkg; Pkg.test(test_args=["--local"])'

# CI test suite (no physics, uses COSMO)
test-ci:
	$(JL) -e 'using Pkg; Pkg.test()'

# Individual test groups
test-polynomials:
	$(JL) -e 'using Pkg; Pkg.test(test_args=["--polynomials"])'

test-quality:
	$(JL) -e 'using Pkg; Pkg.test(test_args=["--quality"])'

test-solvers:
	$(JL) -e 'using Pkg; Pkg.test(test_args=["--solvers"])'

test-physics:
	$(JL) -e 'using Pkg; Pkg.test(test_args=["--physics", "--local"])'

# Combination targets
test-no-physics:
	$(JL) -e 'using Pkg; Pkg.test(test_args=["--polynomials", "--quality", "--solvers"])'

test-core:
	$(JL) -e 'using Pkg; Pkg.test(test_args=["--polynomials", "--solvers"])'

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

.PHONY: init init-docs update update-docs \
        test test-ci test-polynomials test-quality test-solvers test-physics \
        test-no-physics test-core \
        servedocs examples clean
