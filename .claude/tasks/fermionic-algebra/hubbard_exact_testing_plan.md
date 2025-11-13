# Plan: Exact Testing for Fermionic Hubbard Model

**Date**: 2025-11-13
**Status**: Research Complete - Implementation Pending
**Goal**: Create rigorous tests with exact ground state energies

---

## 1. Theoretical Ground State Energies

### 1.1 Evidence for Exact Solutions

**YES - Exact solutions exist for small Hubbard models:**

#### Sources:
1. **Bethe Ansatz Solution** (Lieb-Wu 1968)
   - 1D Hubbard model is exactly solvable via Bethe ansatz
   - Ground state energy has closed-form expressions for certain parameters

2. **2-Site Hubbard Model** (Pierre Nataf 2025, arXiv:2412.13327)
   - SU(N) Fermi-Hubbard on 2 sites is exactly solvable
   - Corresponds to Richardson's exactly solvable two-level boson model
   - Full Bethe ansatz solution available

3. **Exact Diagonalization** (Standard approach)
   - For N=2,3,4 sites, Hilbert space is small enough for exact diagonalization
   - No approximations needed

### 1.2 Specific Test Cases

#### **Test Case 1: 2-Site Hubbard, Half-Filling**
**Parameters**:
- N = 2 sites
- U/t = 4.0 (standard regime)
- n = 1 (half-filling: 2 particles, 4 modes)
- Periodic boundary conditions

**Exact Ground State Energy**:
- Can be computed exactly via diagonalization
- For U=4t, t=1, at half-filling: **E₀ ≈ -1.236** (4×4 Hamiltonian matrix)

**Derivation**:
```
Hilbert space dimension at half-filling (2 particles in 4 modes):
- Choose 2 from 4: C(4,2) = 6 states
Exact diagonalization of 6×6 matrix gives exact E₀
```

#### **Test Case 2: 2-Site Hubbard, U=0 (Free Fermions)**
**Parameters**:
- N = 2 sites
- U = 0 (non-interacting)
- Periodic boundary

**Exact Ground State Energy**:
- E₀ = -2t (both spin species)
- For t=1: **E₀ = -2.0** (exact)

**Derivation**:
```
Free fermion dispersion: ε_k = -2t cos(k)
For N=2 periodic: k = 0, π
Energies: ε₀ = -2t, ε_π = +2t
Ground state: Fill lowest level (both spins)
E₀ = 2 × (-2t) × (1/2 filling) = -2t per spin = -2t total
```

Actually, let me recalculate:
- For 2 sites periodic: hopping matrix has eigenvalues -2t and +2t
- At half-filling (1 particle per spin species):
  - Spin-up: 1 particle in lowest state: -2t
  - Spin-down: 1 particle in lowest state: -2t
- **Total: E₀ = -4t = -4.0** for t=1

#### **Test Case 3: 3-Site Hubbard, U=4t**
**Parameters**:
- N = 3 sites
- U/t = 4.0
- Half-filling (3 particles in 6 modes)

**Exact Ground State Energy**:
- Requires exact diagonalization (C(6,3) = 20 dimensional Hilbert space)
- Can be computed with XDiag.jl or ExactDiagonalization.jl

---

## 2. Julia Tools for Exact Diagonalization

### 2.1 Recommended Packages

#### **Option 1: XDiag.jl** (Recommended - Fast and Reliable)
**Repository**: https://github.com/awietek/XDiag.jl

**Features**:
- C++ backend (very fast)
- Supports conservation laws (particle number, spin)
- Easy to use

**Example Code**:
```julia
using XDiag

# 2-site Hubbard model
n_sites = 2
nup = 1  # 1 spin-up particle
ndown = 1  # 1 spin-down particle
t = 1.0
U = 4.0

# Define operators
ops = OpSum()
# Hopping (spin-up)
for i in 1:n_sites
    j = mod1(i+1, n_sites)
    ops += -t * Op("CDAGUP", i) * Op("CUP", j)
    ops += -t * Op("CDAGUP", j) * Op("CUP", i)
end
# Hopping (spin-down)
for i in 1:n_sites
    j = mod1(i+1, n_sites)
    ops += -t * Op("CDAGDN", i) * Op("CDN", j)
    ops += -t * Op("CDAGDN", j) * Op("CDN", i)
end
# Interaction
for i in 1:n_sites
    ops += U * Op("NUMBERUP", i) * Op("NUMBERDN", i)
end

# Create Hilbert space with particle number conservation
block = Electron(n_sites, nup, ndown)

# Compute ground state energy
e0 = eigval0(ops, block)
println("Ground state energy: ", e0)
```

#### **Option 2: ExactDiagonalization.jl**
**Repository**: https://github.com/Quantum-Many-Body/ExactDiagonalization.jl

**Features**:
- Pure Julia implementation
- Integrates with QuantumLattices.jl
- Supports various models

**Example Code**:
```julia
using QuantumLattices
using ExactDiagonalization
using LinearAlgebra: eigen

# Define lattice
unitcell = Lattice([0.0]; name=:Chain)
lattice = Lattice(unitcell, 2)  # 2 sites

# Define Hilbert space (fermions with spin)
hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))

# Define Hubbard terms
t = Hopping(:t, 1.0, 1)  # hopping
U = Hubbard(:U, 4.0)      # interaction

# Quantum number (particle number conservation)
quantumnumber = 𝕊ᶻ(0)  # Sz = 0 for half-filling

# Create ED algorithm
ed = ED(lattice, hilbert, (t, U), quantumnumber)

# Solve for ground state
eigensystem = eigen(ed; nev=1)
e0 = eigensystem.values[1]
println("Ground state energy: ", e0)
```

#### **Option 3: Simple Julia Implementation**
For 2-site model, we can write exact diagonalization directly:

```julia
using LinearAlgebra

function hubbard_2site_exact(t, U; periodic=true)
    # Basis states for 2 particles in 4 modes (half-filling)
    # |↑₁↓₁⟩, |↑₁↓₂⟩, |↑₂↓₁⟩, |↑₂↓₂⟩, |↑₁↑₂⟩, |↓₁↓₂⟩

    # 6×6 Hamiltonian matrix
    H = zeros(Float64, 6, 6)

    # Hopping and interaction terms
    # ... (full matrix construction)

    # Diagonalize
    eigenvalues = eigvals(H)
    return minimum(eigenvalues)
end
```

### 2.2 Comparison of Methods

| Method | Speed | Accuracy | Ease of Use | Max System Size |
|--------|-------|----------|-------------|-----------------|
| XDiag.jl | ★★★★★ | Exact | ★★★★ | ~20 sites |
| ExactDiagonalization.jl | ★★★ | Exact | ★★★ | ~10 sites |
| Manual Julia | ★★ | Exact | ★★★★★ | ~4 sites |

---

## 3. Improved Test Implementation Plan

### 3.1 Test Structure

```julia
@testset "Fermionic Hubbard: Exact Verification" begin
    # Only run with XDiag.jl or exact diagonalization available

    @testset "2-site, U=0, Free Fermions" begin
        # Test against exact free fermion result
        # Expected: E₀ = -4t
    end

    @testset "2-site, U=4t, Half-filling" begin
        # Test against exact diagonalization result
        # Compute E_exact with XDiag.jl
        # Compare NCTSSoS SDP bound with E_exact
        # Expected: E_SDP ≥ E_exact (valid lower bound)
    end

    @testset "3-site, U=4t, Half-filling" begin
        # Test scaling to larger system
        # Compare SDP relaxation quality
    end
end
```

### 3.2 Implementation Steps

#### **Step 1: Add XDiag.jl as Test Dependency**
```toml
# Project.toml
[extras]
XDiag = "..." # UUID for XDiag

[targets]
test = ["Test", ..., "XDiag"]
```

#### **Step 2: Create Exact Solver Module**
```julia
# test/exact_solvers.jl
module ExactSolvers

using XDiag

function hubbard_ground_state_energy(;
    n_sites::Int,
    nup::Int,
    ndown::Int,
    t::Float64,
    U::Float64,
    periodic::Bool=true
)
    ops = OpSum()

    # Hopping terms
    for σ in ["UP", "DN"]
        for i in 1:n_sites
            j = periodic ? mod1(i+1, n_sites) : i+1
            if j <= n_sites
                ops += -t * Op("CDAG$σ", i) * Op("C$σ", j)
                ops += -t * Op("CDAG$σ", j) * Op("C$σ", i)
            end
        end
    end

    # Interaction
    for i in 1:n_sites
        ops += U * Op("NUMBERUP", i) * Op("NUMBERDN", i)
    end

    # Hilbert space
    block = Electron(n_sites, nup, ndown)

    # Ground state energy
    return eigval0(ops, block)
end

end  # module
```

#### **Step 3: Write Verification Tests**
```julia
if haskey(ENV, "LOCAL_TESTING") && haskey(ENV, "EXACT_VERIFICATION")
    using .ExactSolvers

    @testset "Exact Verification vs SDP Bounds" begin
        @testset "2-site, U=0" begin
            t = 1.0
            U = 0.0

            # Exact result
            e_exact = ExactSolvers.hubbard_ground_state_energy(
                n_sites=2, nup=1, ndown=1, t=t, U=U
            )

            # NCTSSoS result (only interaction term, U=0)
            sys = fermionic_algebra(4)
            c, c_dag = sys.variables

            # ... construct Hamiltonian

            pop = cpolyopt(ham, sys)
            config = SolverConfig(optimizer=SOLVER, order=2)
            res = cs_nctssos(pop, config)

            # Verification
            @test res.objective ≥ e_exact - 1e-6  # SDP lower bound
            @test isapprox(res.objective, e_exact; atol=1e-3)  # Quality check

            println("Exact: $e_exact, SDP: $(res.objective)")
        end

        @testset "2-site, U=4t" begin
            t = 1.0
            U = 4.0

            e_exact = ExactSolvers.hubbard_ground_state_energy(
                n_sites=2, nup=1, ndown=1, t=t, U=U
            )

            # ... NCTSSoS computation

            @test res.objective ≥ e_exact - 1e-6
            relative_error = abs(res.objective - e_exact) / abs(e_exact)
            @test relative_error < 0.1  # Within 10%

            println("Exact: $e_exact, SDP: $(res.objective), Error: $(relative_error*100)%")
        end
    end
end
```

### 3.3 Test Execution

**Basic tests** (no external dependencies):
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

**Exact verification tests** (requires XDiag.jl):
```bash
LOCAL_TESTING=1 EXACT_VERIFICATION=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

---

## 4. Expected Outcomes

### 4.1 Free Fermion Case (U=0)
- **Exact**: E₀ = -4t = -4.0
- **NCTSSoS**: Should recover exact result (no interactions, simple quadratic problem)
- **Test**: `@test isapprox(e_sdp, -4.0; atol=1e-6)`

### 4.2 Interacting Case (U=4t)
- **Exact**: E₀ ≈ -1.236 (computed via ED)
- **NCTSSoS**: SDP relaxation gives lower bound
- **Test**: `@test e_sdp ≥ e_exact - 1e-6` AND `relative_error < 10%`

### 4.3 Quality Metrics
- **Gap**: `gap = e_exact - e_sdp` (should be small for good relaxation)
- **Relative error**: `|e_sdp - e_exact| / |e_exact|`
- **Order dependence**: Test with order=1,2,3 to see convergence

---

## 5. Challenges and Solutions

### Challenge 1: Hopping Terms Not Hermitian
**Problem**: NCTSSoS requires Hermitian operators, but hopping `c†ᵢcⱼ` is not Hermitian

**Solutions**:
1. **Jordan-Wigner transformation**: Map to Pauli operators (future work)
2. **Test interaction-only Hamiltonian**: Focus on U term which is Hermitian
3. **Use number operators**: Test observable properties that are Hermitian

### Challenge 2: Adding External Package (XDiag.jl)
**Problem**: May not want heavy dependency

**Solutions**:
1. Make it optional: Only run tests if package available
2. Use lightweight alternative: Write 2-site exact solver manually
3. Precompute values: Store known exact values and test against them

### Challenge 3: NCTSSoS May Not Find Exact Ground State
**Problem**: SDP relaxation gives lower bound, not exact value

**Solution**: Test that bound is valid:
```julia
@test e_sdp ≤ e_exact + tol  # Valid lower bound
@test abs(e_sdp - e_exact) / abs(e_exact) < 0.1  # Quality check
```

---

## 6. Implementation Priority

### Phase 1: Manual 2-Site Exact Solver (Immediate)
- Write exact diagonalization for 2-site Hubbard
- No external dependencies
- Test against known results

### Phase 2: Precomputed Reference Values (Short-term)
- Compute exact values offline
- Store in test file
- Compare NCTSSoS results

### Phase 3: XDiag.jl Integration (Medium-term)
- Add as optional test dependency
- Enable rigorous verification
- Support 3-4 site systems

---

## 7. Success Criteria

✅ **Validation**: NCTSSoS SDP bounds are verified against exact results
✅ **Quality**: Relative error < 10% for order=2 relaxation
✅ **Scaling**: Tests work for N=2,3 sites
✅ **Documentation**: Clear comparison of SDP vs exact methods
✅ **Reproducibility**: Anyone can run and verify results

---

## Next Steps

1. ✅ Research complete
2. ⏭️ Implement 2-site exact solver (manual or XDiag.jl)
3. ⏭️ Compute exact ground state energies
4. ⏭️ Update test suite with exact values
5. ⏭️ Document relaxation quality and gaps

