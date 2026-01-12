JL = julia --project

# Comma variable for escaping in macro calls
, := ,

# Macro for Pkg.test calls - usage: $(call pkg_test,["--flag1"$(,)"--flag2"])
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
#   --minimal      Fast smoke test (5 cases) for relaxation correctness
#   --polynomials  Run polynomial algebra tests
#   --quality      Run code quality checks
#   --solvers      Run SDP solver tests
#   --physics      Run physics model tests
#
# Usage:
#   make test              - Full suite with Mosek
#   make test-ci           - CI suite (no physics, uses COSMO)
#   make test-minimal      - Fast smoke test (5 cases, ~30s) for relaxation correctness
#   make test-minimal-local - Minimal with Mosek
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
	$(call pkg_test,["--local"])

# CI test suite (no physics, uses COSMO)
test-ci:
	$(call pkg_test,)

# Individual test groups
test-polynomials:
	$(call pkg_test,["--polynomials"])

test-quality:
	$(call pkg_test,["--quality"])

test-solvers:
	$(call pkg_test,["--solvers"])

test-physics:
	$(call pkg_test,["--physics"$(,)"--local"])

# Minimal correctness smoke test (5 cases covering critical algorithm paths)
# Covers: Dense, CS, TS, CS+TS, Dualization - ~30s with COSMO
test-minimal:
	$(call pkg_test,["--minimal"])

test-minimal-local:
	$(call pkg_test,["--minimal"$(,)"--local"])

# Combination targets
test-no-physics:
	$(call pkg_test,["--polynomials"$(,)"--quality"$(,)"--solvers"])

test-core:
	$(call pkg_test,["--polynomials"$(,)"--solvers"])

# Run a single test file with Mosek
# Usage: make test-file FILE=test/relaxations/sparsity.jl
test-file:
ifndef FILE
	$(error FILE is required. Usage: make test-file FILE=test/relaxations/sparsity.jl)
endif
	$(JL) -e 'include("$(FILE)")'

# =============================================================================
# Oracles (requires external NCTSSOS repo)
# =============================================================================
# Regenerate oracle values from NCTSSOS reference implementation.
# Set NCTSSOS_PATH env var or use default locations.
#
# Usage: make oracle-chsh
#        make oracle-i3322
#        NCTSSOS_PATH=/custom/path make oracle-chsh
oracle-%:
	@if [ -z "$$NCTSSOS_PATH" ]; then \
		if [ -d "/Users/yushengzhao/projects/NCTSSOS" ]; then \
			NCTSSOS_PATH="/Users/yushengzhao/projects/NCTSSOS"; \
		elif [ -d "/home/yushengzhao/NCTSSOS" ]; then \
			NCTSSOS_PATH="/home/yushengzhao/NCTSSOS"; \
		else \
			echo "Error: NCTSSOS not found. Set NCTSSOS_PATH environment variable."; \
			exit 1; \
		fi; \
	fi && \
	echo "Using NCTSSOS at: $$NCTSSOS_PATH" && \
	cd "$$NCTSSOS_PATH" && julia --project "$(CURDIR)/test/oracles/scripts/nctssos_$*.jl"

# =============================================================================
# Documentation
# =============================================================================

servedocs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs=["docs/src/assets","docs/src/generated"])'

examples:
	$(JL) docs/generate_examples.jl

# =============================================================================
# Remote Sync (a800 server via mutagen)
# =============================================================================
# Real-time bidirectional sync to a800 GPU server.
# Requires: mutagen (brew install mutagen-io/mutagen/mutagen)
#
# Usage:
#   make sync-start   - Create and start sync session
#   make sync-status  - Check sync status
#   make sync-stop    - Terminate sync session
#   make sync-pause   - Pause syncing
#   make sync-resume  - Resume syncing
#   make sync-flush   - Force immediate sync

SYNC_NAME = nctssos-a800
SYNC_REMOTE = a800:~/projects/NCTSSoS.jl-review-fastpolynomial

sync-start:
	@mutagen sync list | grep -q "$(SYNC_NAME)" && echo "Sync already running" || \
	mutagen sync create \
		--name="$(SYNC_NAME)" \
		--ignore-vcs \
		--ignore="Manifest.toml" \
		--ignore="*.jl.cov" \
		--ignore="*.jl.mem" \
		--ignore="docs/build/" \
		--ignore="docs/site/" \
		--ignore=".DS_Store" \
		--ignore=".vscode/" \
		--ignore="*.json" \
		--sync-mode="two-way-resolved" \
		$(CURDIR) $(SYNC_REMOTE)

sync-status:
	mutagen sync list

sync-stop:
	mutagen sync terminate $(SYNC_NAME)

sync-pause:
	mutagen sync pause $(SYNC_NAME)

sync-resume:
	mutagen sync resume $(SYNC_NAME)

sync-flush:
	mutagen sync flush $(SYNC_NAME)

# =============================================================================
# Cleanup
# =============================================================================

clean:
	rm -rf docs/build
	find . -name "*.cov" -delete

.PHONY: init init-docs update update-docs \
        test test-ci test-polynomials test-quality test-solvers test-physics \
        test-minimal test-minimal-local test-no-physics test-core test-file \
        sync-start sync-status sync-stop sync-pause sync-resume sync-flush \
        servedocs examples clean
