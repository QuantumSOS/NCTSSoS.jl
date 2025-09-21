JL = julia --project

ifndef TARGET
override TARGET = main
endif

default: init test

init:
	$(JL) -e 'using Pkg; Pkg.precompile()'

init-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.develop(path="."); Pkg.precompile()'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile()'

update-docs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); Pkg.update(); Pkg.precompile()'

test:
	LOCAL_TESTING=true $(JL) -e 'using Pkg; Pkg.test()'

servedocs:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs = ["docs/src/assets", "docs/src/generated"])'

# Generate markdown files from Literate.jl examples
examples:
	@echo "Generating markdown files from Literate.jl examples..."
	$(JL) docs/generate_examples.jl

# Build the documentation with examples
docs: examples
	@echo "Building documentation..."
	$(JL) -e 'using Pkg; Pkg.activate("docs"); include("docs/make.jl")'

clean:
	rm -rf docs/build
	find . -name "*.cov" -type f -print0 | xargs -0 /bin/rm -f

test-FastPoly:
	$(JL) -e 'using Pkg; Pkg.status(); include("test/fastpoly_test/runtests.jl")'

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

.PHONY: init test examples docs clean
