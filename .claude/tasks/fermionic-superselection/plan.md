# Implementation Plan: Fermionic Parity Superselection

## Approach

Implement **basis filtering** to remove odd-parity monomials from the moment relaxation basis. This is the cleanest and most efficient approach: it prevents creation of moment variables for unphysical expectations (odd-parity operators), reduces SDP size, and is physically motivated.

The implementation adds a parity check function and filters the total basis in `moment_relax()` for fermionic algebras.

## Steps

1. [x] Add parity detection helper function
   - Location: `src/FastPolynomials/src/simplification/fermionic.jl`
   - Function: `has_even_parity(m::Monomial{FermionicAlgebra,T})`
   - Logic: Return `iseven(length(m.word))`
   - Purpose: Determine if monomial has even number of operators (even parity)

2. [x] Modify moment relaxation to filter basis
   - Location: `src/moment_solver.jl`, after line 150
   - Code: Filter `total_basis` to keep only even-parity monomials
   - Condition: Only for `FermionicAlgebra`
   - Effect: Moment matrix will only have variables for even-parity expectations

3. [x] Export parity function in FastPolynomials
   - Location: `src/FastPolynomials/src/FastPolynomials.jl`
   - Add: `export has_even_parity`
   - Purpose: Make function available for testing and external use

4. [x] Create parity validation test
   - Location: `test/fermionic_parity_test.jl` (new file)
   - Purpose: Verify that moment basis contains only even-parity monomials
   - Test setup: Simple 2-mode fermionic system with Hamiltonian
   - Assertion: All monomials in `mp.total_basis` have `has_even_parity() == true`

5. [!] Validate with existing fermionic test
   - Location: `test/fermionic_chain.jl`
   - Action: Run existing test with parity filtering enabled
   - NOTE: Test was already failing BEFORE this PR (primal infeasible)
   - This is a pre-existing issue, not caused by parity filtering
   - Parity filtering did not break anything that was working

6. [ ] Document parity assumption
   - Location: `src/pop.jl`, in `PolyOpt` docstring
   - Addition: Note for `FermionicAlgebra` about parity superselection
   - Content: Explain that fermionic problems assume physical states have definite parity

7. [x] (Optional) Add objective validation
   - Location: `src/pop.jl`, in `polyopt()` constructor
   - Check: If objective has odd parity for fermionic problem
   - Action: Throw informative error suggesting reformulation
   - Purpose: Catch unphysical problem formulations early

## Testing Strategy

### Unit Test (New)
**File:** `test/fermionic_parity_test.jl`

```julia
@testset "Fermionic Parity Superselection" begin
    # 2-mode fermionic system
    registry, (a, a_dag) = create_fermionic_variables(1:2)
    H = a_dag[1] * a[2] + a_dag[2] * a[1]  # Even parity

    pop = polyopt(H, registry)
    solver_config = SolverConfig(order=2, ts_algo=MMD())
    corr_sparsity = correlative_sparsity(pop, solver_config)
    term_sparsities = [term_sparsity(pop, corr_sparsity, clique) for clique in 1:length(corr_sparsity.cliques)]
    mp = moment_relax(pop, corr_sparsity, term_sparsities)

    # All basis monomials must have even parity
    for mono in mp.total_basis
        @test has_even_parity(mono)
    end
end
```

### Integration Test (Existing)
**File:** `test/fermionic_chain.jl`

- Run with `LOCAL_TESTING=true make test`
- Should still pass: `res.objective ≈ -4.0 atol=1e-4`
- If fails, investigate which odd-parity terms were in basis before

### Edge Cases to Test

1. **Identity monomial:** Empty word (length 0) → even parity ✓
2. **Single operator:** `a₁` or `a₁†` (length 1) → odd parity → filtered out ✓
3. **Pair product:** `a₁† a₂` (length 2) → even parity ✓
4. **Long product:** `a₁† a₂† a₃ a₄` (length 4) → even parity ✓

## Implementation Notes

### Why Filtering is Better Than Constraints

**Filtering (Recommended):**
- ✅ Smaller SDP (fewer variables)
- ✅ Faster solve time
- ✅ Physically motivated (unphysical expectations never created)
- ✅ Simple implementation (~5 lines)

**Explicit Constraints (Alternative):**
- ❌ Larger SDP (variables + constraints)
- ❌ Slower solve time
- ❌ Redundant (filtering achieves same result)
- ✅ Numerically verifiable (solver confirms y[M]=0)

**Trust Simplification (Current?):**
- ❓ Unclear if current implementation accidentally handles this
- ❓ Depends on problem structure
- ❌ Not guaranteed for all inputs
- ❌ No explicit protection against unphysical results

### Potential Issues and Mitigations

**Issue 1:** Filtering removes too many terms, breaking existing tests
- **Diagnosis:** Check if odd-parity terms are actually needed
- **Mitigation:** Investigate why tests relied on odd-parity expectations
- **Likely cause:** Misunderstanding of basis construction

**Issue 2:** Hamiltonian itself has odd parity
- **Diagnosis:** User error - physical Hamiltonians must have even parity [H,P]=0
- **Mitigation:** Add validation in `polyopt()` constructor
- **Action:** Throw error with explanation

**Issue 3:** User wants to optimize ⟨a₁⟩ (odd objective)
- **Diagnosis:** Unphysical problem - result is always zero
- **Mitigation:** Detect and error in `polyopt()` constructor
- **Suggestion:** Reformulate as even-parity objective

## Open Questions for User

1. **Particle number conservation:** Do you also need fixed N (stronger than parity)?
   - Parity only: States can have N=0,2,4,... or N=1,3,5,...
   - Fixed N: State has exactly N particles
   - If yes, requires additional constraints beyond this PR

2. **Spinful fermions:** Are you working with spin-1/2 fermions?
   - Current approach works for spinful (count all operators)
   - May need separate spin-up/spin-down indexing scheme
   - Clarify variable encoding if needed

3. **Complex coefficients:** Should FermionicAlgebra use ComplexF64?
   - Currently uses Float64 (real coefficients)
   - Uses HPSD cone (Hermitian matrices)
   - Physical Hamiltonians can have complex hopping terms
   - Separate investigation needed?

## Citations

1. Siddhant Midha - Parity Superselection: https://siddhantmidha.com/blog/2025/fermions/
2. Lexin Ding (2020) - Fermionic Entanglement thesis: https://arxiv.org/pdf/2207.03848
3. Tosini et al. (2015) - Fermionic quantum theory: https://www.cs.ox.ac.uk/qpl2015/slides/21.pdf
4. NCTSSoS codebase - Implementation details: src/FastPolynomials/src/simplification/fermionic.jl
