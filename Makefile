JL = julia --project

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
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.develop(path="."); Pkg.precompile()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile()'

update-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.update(); Pkg.precompile()'

# =============================================================================
# Testing
# =============================================================================
# Canonical testing documentation: TESTING.md (repo root).
# =============================================================================

# Full test suite with Mosek
test:
	$(call pkg_test,["--local"])

# CI test suite (fast feedback; COSMO). Includes TS layer-1 + layer-2 coverage.
test-ci:
	$(call pkg_test,["--polynomials"$(,)"--relaxations"$(,)"--minimal"])

# Individual test groups
test-polynomials:
	$(call pkg_test,["--polynomials"])

test-quality:
	$(call pkg_test,["--quality"])

test-relaxations:
	$(call pkg_test,["--relaxations"])

test-problems:
	$(call pkg_test,["--problems"$(,)"--local"])

# Minimal correctness smoke test (5 cases covering critical algorithm paths)
# Covers: Dense, CS, TS, CS+TS, Dualization - ~30s with COSMO
test-minimal:
	$(call pkg_test,["--minimal"])

test-minimal-local:
	$(call pkg_test,["--minimal"$(,)"--local"])

# Combination targets
test-no-problems:
	$(call pkg_test,["--polynomials"$(,)"--quality"$(,)"--relaxations"])

test-core:
	$(call pkg_test,["--polynomials"$(,)"--relaxations"])

# Run a single test file with Mosek
test-file:
ifndef FILE
	$(error FILE is required. See TESTING.md.)
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
        test test-ci test-polynomials test-quality test-relaxations test-problems \
        test-minimal test-minimal-local test-no-problems test-core test-file \
        sync-start sync-status sync-stop sync-pause sync-resume sync-flush \
        servedocs examples clean
