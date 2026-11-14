JULIA ?= julia
JL = $(JULIA) --project

# Comma variable for escaping in macro calls
, := ,

define pkg_test
	$(JL) -e 'using Pkg; Pkg.test($(if $(1),test_args=$(1)))'
endef

default: init test

# =============================================================================
# Setup
# =============================================================================

init:
	$(JL) -e 'using Pkg; Pkg.precompile()'

init-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.develop(path="."); Pkg.instantiate(); Pkg.precompile()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile()'

update-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.update(); Pkg.precompile()'

# =============================================================================
# Testing
# =============================================================================
# Canonical testing documentation: TESTING.md (repo root).
# =============================================================================

# Full test suite (COSMO)
test:
	$(call pkg_test)

# Standalone local-only scripts that are intentionally excluded from
# `test/runtests.jl` (for example, heavier Mosek-only reproductions).
LOCAL_ONLY_TEST_SCRIPTS = \
	test/problems/trace_polynomial/t1_broyden_banded_trace.jl \
	test/problems/trace_polynomial/t4_nc_motzkin_polynomial.jl

test-local:
	@set -e; \
	for script in $(LOCAL_ONLY_TEST_SCRIPTS); do \
		echo "==> $$script"; \
		$(JL) $$script; \
	done

test-local-one:
	@test -n "$(SCRIPT)" || (echo "Usage: make test-local-one SCRIPT=path/to/script.jl"; exit 1)
	$(JL) $(SCRIPT)

# CI-style coverage (lcov.info), matching julia-actions/julia-runtest + julia-processcoverage.
coverage-ci:
	CI=true $(JULIA) --color=yes -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path="."); Pkg.instantiate(); Pkg.test("NCTSSoS"; coverage=true)'
	$(JULIA) --color=yes -e 'using Pkg; Pkg.activate("coveragetempenv", shared=true); Pkg.add(PackageSpec(name="CoverageTools")); using CoverageTools; directories = get(ENV, "INPUT_DIRECTORIES", "src,ext"); dirs = filter!(!isempty, split(directories, ",")); for dir in dirs; if dir == "ext"; continue; elseif !isdir(dir); error("directory \\\"" * dir * "\\\" not found!"); end; end; filter!(isdir, dirs); pfs = mapreduce(process_folder, vcat, dirs); LCOV.writefile("lcov.info", pfs)'

# =============================================================================
# Oracles (requires external NCTSSOS repo)
# =============================================================================
# Regenerate oracle values from NCTSSOS reference implementation.
# Set NCTSSOS_PATH to the legacy repository.
#
# Usage: make oracle-chsh
#        make oracle-i3322
#        NCTSSOS_PATH=/custom/path make oracle-chsh
oracle-%:
	@if [ -z "$$NCTSSOS_PATH" ]; then \
		echo "Error: NCTSSOS not found. Set NCTSSOS_PATH environment variable."; \
		exit 1; \
	fi && \
	echo "Using NCTSSOS at: $$NCTSSOS_PATH" && \
	cd "$$NCTSSOS_PATH" && julia --project "$(CURDIR)/test/oracles/nctssos_$*.jl"

# =============================================================================
# Documentation
# =============================================================================

servedocs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.develop(path="."); Pkg.instantiate(); using LiveServer; servedocs(;skip_dirs=["docs/src/assets","docs/src/generated"])'

examples:
	$(JL) docs/generate_examples.jl

# =============================================================================
# Cleanup
# =============================================================================

clean:
	rm -rf docs/build
	find . -name "*.cov" -delete

.PHONY: init init-docs update update-docs \
        test test-local test-local-one coverage-ci \
        servedocs examples clean
