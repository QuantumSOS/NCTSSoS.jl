JL = julia --project

ifndef TARGET
override TARGET = main
endif

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
# Test categories:
#   - polynomials/  Core algebra (no solver)
#   - quality/      Code quality checks
#   - solvers/      SDP solver integration
#   - physics/      Physics models (LOCAL_TESTING only)
#
# Solver config:
#   - Default: COSMO (open-source)
#   - LOCAL_TESTING=true: Mosek (commercial, required for physics tests)
# =============================================================================

# Full test suite with Mosek (requires LOCAL_TESTING for physics tests)
test:
	LOCAL_TESTING=true $(JL) -e 'using Pkg; Pkg.test()'

# Quick test with COSMO (no physics tests, CI-compatible)
test-quick:
	$(JL) -e 'using Pkg; Pkg.test()'

# Test individual categories
test-polynomials:
	$(JL) -e 'using NCTSSoS, Test; include("test/polynomials/runtests.jl")'

test-solvers:
	$(JL) -e 'using NCTSSoS, Test; include("test/setup.jl"); \
		include("test/solvers/moment.jl"); \
		include("test/solvers/sos.jl"); \
		include("test/solvers/interface.jl"); \
		include("test/solvers/sparsity.jl"); \
		include("test/solvers/state_poly.jl"); \
		include("test/solvers/trace_poly.jl")'

test-physics:
	LOCAL_TESTING=true $(JL) -e 'using NCTSSoS, Test; include("test/setup.jl"); \
		include("test/physics/heisenberg.jl"); \
		include("test/physics/xy_model.jl"); \
		include("test/physics/bose_hubbard.jl"); \
		include("test/physics/bell_inequalities.jl"); \
		include("test/physics/fermionic.jl")'

# =============================================================================
# Documentation
# =============================================================================

servedocs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs = ["docs/src/assets", "docs/src/generated"])'

# Generate markdown files from Literate.jl examples
examples:
	@echo "Generating markdown files from Literate.jl examples..."
	$(JL) docs/generate_examples.jl

# =============================================================================
# Benchmarks
# =============================================================================

init-bench:
	$(JL) -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("AirspeedVelocity"); Pkg.build("AirspeedVelocity")'

bench: 
	benchpkg \
	--rev $(TARGET),dirty \
	--script "benchmark/benchmarks.jl"

benchtable:
	benchpkgtable \
	--rev $(TARGET),dirty \
	--ratio true

# =============================================================================
# Cleanup
# =============================================================================

clean:
	rm -rf docs/build
	find . -name "*.cov" -type f -print0 | xargs -0 /bin/rm -f

.PHONY: init init-docs update update-docs \
        test test-quick test-polynomials test-solvers test-physics \
        servedocs examples \
        init-bench bench benchtable \
        clean
