# Implementation Plan: Simplification Scheme Review

## Approach
This is a READ-ONLY review task. The goal is to systematically examine all simplification-related code across the six algebra types to understand:
1. Algorithmic approaches and their correctness
2. Consistency vs inconsistency in design patterns
3. Performance characteristics and potential bottlenecks
4. Integration with canonicalization and other operations
5. Test coverage and edge cases

The review will be conducted in pair programming mode with the user, working through each area methodically.

## Review Structure
Each review area below includes:
- **Files**: Specific files and line ranges to examine
- **Questions**: Key questions to answer during review
- **Tests**: Corresponding test files to check coverage

## Steps

### 1. [x] [DONE] Review: Core Algebra Type Definitions
**Files:**
- `src/FastPolynomials/src/algebra_types.jl` (lines 1-290)

**Questions:**
- Are the algebraic rules documented correctly for each type?
- Is the site-based encoding scheme (encode_index, decode_site) correct?
- Are the bit allocation formulas (site_bits, max_sites, max_operators) optimal?
- Should select_uint_type be more conservative or aggressive?

**Tests:**
- `test/fastpoly_test/algebra_types.jl`

---

### 2. [x] [DONE] Review: NonCommutativeAlgebra Simplification
**Files:**
- `src/FastPolynomials/src/simplification/noncommutative.jl` (lines 1-124)

**Questions:**
- Does the site-based sorting algorithm correctly implement commutation between sites?
- Why do we have both Unsigned (site-aware) and Signed (legacy) versions?
- Is the bubble sort for site ordering acceptable? Should it be a stable sort?
- Are empty word edge cases handled correctly?

**Tests:**
- `test/fastpoly_test/simplify.jl` (NonCommutativeAlgebra sections)
- `test/fastpoly_test/monomials.jl`

---

### 3. [x] [DONE - REFACTORED] Review: PauliAlgebra Simplification
**Files:**
- `src/FastPolynomials/src/simplification/pauli.jl` (lines 1-273)

**Questions:**
- Is the Pauli product table (_pauli_product) mathematically correct?
- Does the site sorting correctly handle commutation between different sites?
- Is the cyclic product rule (XY=iZ, YZ=iX, ZX=iY) implemented correctly?
- Why does simplify! use a while loop with `changed` flag? Could this be optimized?
- Should we simplify during multiplication (see TODO at line 270)?
- Are complex phases tracked correctly throughout?

**Tests:**
- `test/fastpoly_test/simplify.jl` (PauliAlgebra sections)

---

### 4. [x] [DONE] Review: FermionicAlgebra Simplification (Wick's Theorem)
**Files:**
- `src/FastPolynomials/src/simplification/fermionic.jl` (lines 1-407)

**Questions:**
- Is the Wick's theorem implementation correct per the generalized time-independent form?
- Does _find_valid_contractions correctly identify all annihilation-creation pairs?
- Is _generate_nonoverlapping_combinations exhaustive and correct?
- Does _permutation_parity correctly compute sign from cycle counting?
- Is _contraction_sign correct in counting swaps needed?
- Does iszero(m::Monomial{FermionicAlgebra}) correctly detect nilpotent terms?
- Could the algorithm be optimized for common cases (e.g., already normal ordered)?
- Is the complexity acceptable? (Currently potentially exponential in worst case)

**Tests:**
- `test/fastpoly_test/simplify.jl` (FermionicAlgebra sections)

---

### 5. [x] [DONE] Review: BosonicAlgebra Simplification (Rook Numbers)
**Files:**
- `src/FastPolynomials/src/simplification/bosonic.jl` (lines 1-519)

**Questions:**
- Is the rook number algorithm implementation correct per arXiv:quant-ph/0507206?
- Does build_ferrers_board correctly construct the board from operator sequence?
- Is compute_rook_numbers DP algorithm correct?
- Does single_mode_normal_form apply Eq. 1.40 correctly?
- Is expand_and_construct handling the Cartesian product expansion correctly?
- Could group_by_mode! be optimized further?
- Is the claimed polynomial time complexity accurate? (TODO at line 53)
- Are there edge cases (e.g., all creations, all annihilations) handled correctly?

**Tests:**
- `test/fastpoly_test/simplify.jl` (BosonicAlgebra sections)

---

### 6. [x] [DONE] Review: ProjectorAlgebra Simplification - REFACTORED
**Files:**
- `src/FastPolynomials/src/simplification/projector.jl` (lines 1-203)

**Questions:**
- Is the idempotency rule (P² = P) correctly removing consecutive duplicates?
- Does the site grouping algorithm preserve within-site order correctly?
- Is the Dict-based grouping approach efficient? (TODO at line 86 mentions algorithm seems bad)
- Should we use a different approach (see NCTSSOS reference in TODO)?
- Are there performance issues with the current implementation?

**Tests:**
- `test/fastpoly_test/simplify.jl` (ProjectorAlgebra sections)

---

### 7. [x] [DONE] Review: UnipotentAlgebra Simplification - REFACTORED
**Files:**
- `src/FastPolynomials/src/simplification/unipotent.jl` (lines 1-253)

**Questions:**
- Is the U² = I rule correctly removing consecutive pairs?
- Does site-based sorting work correctly?
- Why do we have both Unsigned (site-aware) and Signed (legacy) versions?
- Is the bubble sort for site ordering acceptable?
- How does this differ from ProjectorAlgebra? (U²=I vs P²=P)

**Tests:**
- `test/fastpoly_test/simplify.jl` (UnipotentAlgebra sections)

---

### 8. [ ] Review: Canonicalization Algorithms
**Files:**
- `src/FastPolynomials/src/canonicalization.jl` (lines 1-344)

**Questions:**
- Is symmetric_canon correctly choosing lexicographically smaller of word and reverse?
- Is cyclic_canon finding the minimum rotation correctly? O(n²) acceptable?
- Does cyclic_symmetric_canon correctly combine both transformations?
- When should cyclic=true vs cyclic=false be used?
- Should canonicalization happen before or after simplification? (TODO at line 143)
- Is the Polynomial canonicalize correctly combining like terms?

**Tests:**
- `test/fastpoly_test/canonicalization.jl`

---

### 9. [x] [DONE] Review: Multiplication Consistency Across Algebras
**Files:**
- All simplification/*.jl files, focus on `Base.:*` definitions

**Questions:**
- Is the pattern of concatenate + simplify! consistent across all algebras?
- Are empty word edge cases handled consistently?
- Should multiplication optimize for special cases (identity, already simplified)?
- Is the return type inconsistency (Term vs Vector{Term}) necessary?

**COMPLETED:** Refactored all `*` functions to just concatenate words and return Monomial.
Simplification is now caller's responsibility. All 1060 FastPolynomials tests pass.

**Tests:**
- `test/fastpoly_test/arithmetic.jl`

---

### 10. [ ] Review: Adjoint Operations
**Files:**
- `src/FastPolynomials/src/monomial.jl` (adjoint definition)
- `src/FastPolynomials/src/simplification/unipotent.jl` (specialized adjoint for Signed)

**Questions:**
- Is the default adjoint (reverse word) correct for all algebras?
- Why does UnipotentAlgebra need a specialized Signed adjoint?
- Should other algebras have specialized adjoints?
- Is sign handling for Signed types correct (-idx for adjoint)?

**Tests:**
- `test/fastpoly_test/monomials.jl`

---

### 11. [ ] Review: Legacy Compatibility Layer
**Files:**
- `src/FastPolynomials/src/utils.jl` (lines 161-832)

**Questions:**
- Is the SimplifyAlgorithm compatibility layer complete?
- Are there cases where legacy code would break with the new API?
- Should we deprecate SimplifyAlgorithm or maintain it long-term?
- Is the Variable struct still necessary? What uses it?
- Is the mapping from comm_gps to AlgebraType correct?

**Tests:**
- Check which tests still use SimplifyAlgorithm vs direct AlgebraType dispatch

---

### 12. [ ] Review: Test Coverage Analysis
**Files:**
- `test/fastpoly_test/simplify.jl`
- `test/fastpoly_test/canonicalization.jl`
- `test/fastpoly_test/arithmetic.jl`

**Questions:**
- What percentage of simplification code is covered by tests?
- Are edge cases tested (empty words, single operators, large expressions)?
- Are the algebraic rules verified mathematically (e.g., Pauli anticommutation)?
- Are performance benchmarks in place?
- Do tests cover both Signed and Unsigned index types?

**Cross-reference:**
- `benchmark/bench_monomials.jl`
- `benchmark/bench_polynomials.jl`

---

### 13. [ ] Review: Return Type Heterogeneity Strategy
**Files:**
- All simplification/*.jl files

**Questions:**
- Why do 4 algebras return Term but 2 return Vector{Term}?
- Could we unify to always return Vector{Term} for consistency?
- What would be the performance impact of wrapping single-term results in vectors?
- Does this affect polynomial arithmetic and higher-level code?

**Impact Analysis:**
- Search codebase for all uses of simplify/simplify! results
- Check if code assumes single Term vs handles Vector{Term}

---

### 14. [ ] Review: Site-Based Encoding Applicability
**Files:**
- `src/FastPolynomials/src/algebra_types.jl` (encoding functions)
- Simplification files for NonCommutative, Projector, Unipotent

**Questions:**
- Why is site-based encoding only for Unsigned types?
- Could we extend it to Signed types (for Fermionic/Bosonic with sites)?
- Is the bit allocation scheme optimal for typical use cases?
- Should site_bits be configurable rather than fixed at sizeof(T)*2?

**Explore:**
- Are there real-world use cases that hit the site/operator limits?
- What's the typical distribution of sites and operators per site?

---

### 15. [ ] Review: Integration with Higher-Level Code
**Files:**
- `src/moment_solver.jl` (uses simplify with SimplifyAlgorithm)
- `src/complex_moment_solver.jl`
- `src/sparse.jl`
- `src/sos_solver.jl`

**Questions:**
- How does higher-level code use simplification?
- Are there performance bottlenecks in the simplification pipeline?
- Is canonicalize used appropriately in moment matrix construction?
- Could caching/memoization help?

**Profile:**
- Run examples and profile to identify hot spots
- Check if simplification is called redundantly

---

## Testing Strategy
For each review area:
1. Read the implementation code carefully
2. Read corresponding tests to understand expected behavior
3. Run tests: `julia --project=. test/fastpoly_test/runtests.jl`
4. Run benchmarks if available
5. Experiment with edge cases in REPL
6. Document findings in review notes

## Completion Criteria
- [ ] All 15 review areas examined
- [ ] Questions answered or marked as "needs investigation"
- [ ] Inconsistencies and patterns documented
- [ ] Test coverage gaps identified
- [ ] Performance bottlenecks noted
- [ ] Recommendations summary created (separate doc)

## Output Artifacts
After review, create:
1. `review-findings.md` - Detailed findings for each area
2. `inconsistencies.md` - Catalog of design inconsistencies
3. `recommendations.md` - Suggested improvements (prioritized)
4. `test-gaps.md` - Missing test coverage

## References
- arXiv:quant-ph/0507206 - "Rook numbers and the normal ordering problem" (Bosonic)
- Generalized Time-Independent Wick Theorem (Fermionic)
- NCTSSOS original implementation (for comparison)
