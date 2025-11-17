# Hubbard MPSKit Example Task

## Task Description
Create an example in `/Users/exaclior/projects/NCTSSoS.jl-fermionic-algebra/docs/src/examples/literate` that computes the ground state energy of a fermi-hubbard model at half-filling on a 1D chain lattice of length 4 using `MPSKit.jl` and `MPSKitModels.jl`.

## Requirements
1. Replicate code from https://quantumkithub.github.io/MPSKit.jl/dev/examples/quantum1d/6.hubbard/#hubbard (only ground state energy computation portion)
2. Use MPSKit.jl and MPSKitModels.jl packages
3. Focus on 1D chain lattice of length 4
4. Half-filling condition
5. Compute ground state energy at specific t, mu, and U parameters
6. Create as literate example in docs/src/examples/literate/

## Current Project Status
- NCTSSoS.jl project focused on fermionic algebra and quantum systems
- Existing examples in docs/src/examples/literate/
- Need to understand MPSKit integration patterns

## Key Parameters to Investigate
- t (hopping parameter): 1.0
- mu (chemical potential): 0.0 for half-filling
- U (on-site interaction): 3.0 (example value)
- Lattice size: 4 sites
- Filling: half-filling condition with mu = 0.0
- Bond dimension: 16 for computational efficiency

## MPSKit Implementation Details
- Model function: `hubbard_model(ComplexF64, Trivial, Trivial, FiniteChain(4); t=1.0, U=3.0, mu=0.0)`
- Ground state algorithm: DMRG with `find_groundstate(psi, H, DMRG())`
- Energy extraction: `expectation_value(psi, H)` / number_of_sites
- Analytic reference: Available for infinite system via `hubbard_energy(U/4) - U/4`

## Expected Output
- Literate Julia file with ground state energy computation
- Should work with MPSKit.jl and MPSKitModels.jl
- Clear documentation and explanations
- Verification against analytic solution where applicable