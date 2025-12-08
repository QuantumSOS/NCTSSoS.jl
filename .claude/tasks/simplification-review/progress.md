# Progress Log: simplification-review

## Session: 2025-12-08 20:17 - Checkpoint

**Agent:** orchestrator + polyglot-implementation-engineer
**Feature:** R009 (multiplication-consistency)

### Actions
- Delegated to lead-researcher to analyze multiplication functions across all 6 algebras
- Discovered all usage sites of `neat_dot` and `_neat_dot3` already call `simplify!` explicitly
- Delegated to polyglot-implementation-engineer to refactor multiplication
- Refactored all algebra-specific `*` functions to simply concatenate `word` fields
- Updated `neat_dot`, `_neat_dot3`, `Variable * Variable`, `Variable ^ n` in utils.jl
- Added `_add_simplified_terms!` overload for `Monomial` in polynomial.jl
- Fixed `NCStateWord` multiplication in state_word.jl
- Fixed equality/hash for zero-hash monomials in monomial.jl
- Updated 5 test files to match new behavior

### Outcome
- All `Monomial * Monomial` operations now return `Monomial` (just word concatenation)
- Simplification is caller's responsibility via explicit `simplify!` calls
- All 1060 FastPolynomials tests pass
- Clean separation of concerns achieved

### Files Modified (18 total)
**Source (10):**
- src/FastPolynomials/src/algebra_types.jl
- src/FastPolynomials/src/monomial.jl
- src/FastPolynomials/src/polynomial.jl
- src/FastPolynomials/src/state_word.jl
- src/FastPolynomials/src/utils.jl
- src/FastPolynomials/src/simplification/noncommutative.jl
- src/FastPolynomials/src/simplification/pauli.jl
- src/FastPolynomials/src/simplification/fermionic.jl
- src/FastPolynomials/src/simplification/bosonic.jl
- src/FastPolynomials/src/simplification/projector.jl
- src/FastPolynomials/src/simplification/unipotent.jl

**Tests (8):**
- test/fastpoly_test/algebra_types.jl
- test/fastpoly_test/arithmetic.jl
- test/fastpoly_test/monomials.jl
- test/fastpoly_test/polynomial.jl
- test/fastpoly_test/simplify.jl
- test/fastpoly_test/utils.jl
- test/fastpoly_test/variables.jl

### Next Steps
- Continue with remaining review items (R001-R008, R010-R015)
- R009 is now complete - multiplication is consistent across all algebras

**Commit:** `b01a6a7` - wip(fastpoly): refactor monomial multiplication to pure concatenation - checkpoint

---

## Session: 2025-12-08 21:15

**Agent:** orchestrator
**Feature:** R001 (algebra-definitions)

### Actions
- Read algebra_types.jl documentation (295 lines)
- Verified all 6 algebra type definitions and their algebraic rules
- Validated site-based encoding bit allocation formulas
- Tested `encode_index`/`decode_site`/`decode_operator_id` round-trips
- Validated `select_uint_type` boundary cases

### Findings

**6 Algebra Types Confirmed:**
| Algebra | Index Type | Key Rule | Correct |
|---------|-----------|----------|---------|
| NonCommutativeAlgebra | Unsigned/Signed | No simplification | ✓ |
| PauliAlgebra | Unsigned | σ² = I, cyclic products | ✓ |
| FermionicAlgebra | Signed | CAR: {a, a†} = δ | ✓ |
| BosonicAlgebra | Signed | CCR: [c, c†] = δ | ✓ |
| ProjectorAlgebra | Unsigned | P² = P (idempotent) | ✓ |
| UnipotentAlgebra | Unsigned | U² = I (involution) | ✓ |

**Site Encoding Validated:**
- Formula: `index = (operator_id << site_bits) | site`
- `site_bits(T) = sizeof(T) * 2` (n/4 bits)
- Round-trip encoding verified for UInt8/16/32/64
- `select_uint_type` correctly selects smallest fitting type

**Documentation Issue Found:**
- Line 51 in algebra_types.jl incorrectly calls σᵢ² = I "idempotency"
- Correct term: **involution** (P² = I), not idempotency (P² = P)
- Same issue in pauli.jl line 5
- Low priority: cosmetic only, code is correct

### Outcome
- R001 PASS: All algebra definitions and encodings are mathematically correct
- Minor doc terminology issue noted (non-blocking)

### Next Steps
- Continue with R002 (NonCommutative simplification) or other review items

**Commit:** `f60de30` - docs(fastpoly): fix Pauli algebra terminology

---

## Session: 2025-12-08 21:30

**Agent:** orchestrator
**Feature:** R002 (noncommutative-simplification)

### Actions
- Read noncommutative.jl (165 lines) - site-based simplification for Unsigned types
- Tested Signed type behavior - confirmed MethodError (intentional, no site encoding)
- Verified edge cases: empty word, single element, same site, different sites, mixed

### Findings

**Algorithm (Unsigned only):**
```julia
sort!(word, alg=InsertionSort, by=decode_site)  # Stable sort preserves within-site order
```

**Design Decisions:**
- Only `T<:Unsigned` supported - site encoding requires unsigned bit operations
- Signed types have no simplify method (MethodError) - intentional
- Stable sort ensures within-site operator order is preserved
- InsertionSort used (optimal for small arrays typical in monomials)

**Edge Cases Verified:**
| Case | Input Sites | Output Sites | Status |
|------|-------------|--------------|--------|
| Empty | [] | [] | ✓ |
| Single | [1] | [1] | ✓ |
| Same site | [1,1] | [1,1] (order preserved) | ✓ |
| Different | [3,1,2] | [1,2,3] | ✓ |
| Mixed | [2,1,2,1] | [1,1,2,2] | ✓ |

**Tests:** All 29 simplify tests pass

### Outcome
- R002 PASS: NonCommutative simplification is correct and well-designed
- Unsigned-only constraint is intentional (site encoding requirement)

### Next Steps
- Continue with R003 (Pauli simplification) or other review items

**Commit:** (no code changes - review only)

---

## Session: 2025-12-08 21:45

**Agent:** orchestrator
**Feature:** Refactoring (user-initiated)

### Actions
- User identified that all algebra-specific `Base.:*` methods were identical
- Added single generic `Base.:*(m1::Monomial{A,T}, m2::Monomial{A,T})` in monomial.jl
- Removed 7 duplicate implementations from simplification files
- Fixed hash handling: use `Monomial{A}(vcat(...))` instead of passing `zero(UInt64)`

### Impact
- **Lines removed**: 204
- **Lines added**: 27 (generic impl + docstring)
- **Net reduction**: 177 lines
- **Files changed**: 7 simplification files + monomial.jl

### Tests
All 1060 FastPolynomials tests pass.

**Commit:** `a54f899` - refactor(fastpoly): consolidate Monomial multiplication into single generic method

---

## Session: 2025-12-08 21:55

**Agent:** orchestrator
**Features:** R006 (projector), R007 (unipotent)

### Actions
- Read projector.jl (161 lines) and unipotent.jl (207 lines)
- Tested edge cases for both algebras
- Analyzed TODO about "bad algorithm"
- Compared the two implementations

### Findings

**Both use identical structure:**
1. Group by site using Dict
2. Sort sites ascending
3. Process each site group (algebra-specific)
4. Concatenate results

**Key Difference - Within-Site Processing:**

| Algebra | Rule | Processing |
|---------|------|------------|
| Projector | P² = P | Keep first of consecutive duplicates |
| Unipotent | U² = I | Stack-based pair cancellation |

**Edge Cases Verified:**
| Test | Projector | Unipotent |
|------|-----------|-----------|
| XX | → X | → empty |
| XXX | → X | → X |
| XY | → XY | → XY |
| XYXY | → XYXY | → XYXY |
| XYYX | → XY | → empty |

**TODO Analysis (line 86):**
- The Dict-based approach allocates per-site vectors
- Alternative: stable sort + single pass (like NonCommutative)
- For typical short words (<20 ops), overhead is negligible
- Low priority optimization opportunity

**Unipotent Signed types:**
- Has `adjoint` override for Signed (self-adjoint, just reverses)
- No `simplify` for Signed (MethodError, intentional)

### Outcome
- R006 PASS: Projector idempotency correct
- R007 PASS: Unipotent U²=I correct
- Both could share code (future refactor opportunity)

### Next Steps
- Continue with remaining review items

**Commit:** (no code changes - review only)

---

## Session: 2025-12-08 22:10

**Agent:** orchestrator + lead-researcher
**Features:** Refactoring (Projector + Unipotent simplification)

### Research
- Delegated to lead-researcher to analyze NCTSSOS codebase at `/Users/yushengzhao/projects/NCTSSOS`
- Found NCTSSOS uses `constraint_reduce!` with backtracking (utils.jl:114-130)
- Key insight: stable sort by site + NCTSSOS-style backtrack works for our site-encoded indices

### Implementation
Replaced Dict-based grouping with simpler sort+backtrack:

**Projector (P²=P)**: `sort! + remove consecutive duplicates`
**Unipotent (U²=I)**: `sort! + remove consecutive pairs + backtrack`

### Impact
- **Lines removed**: 61
- **Lines added**: 26
- **Net reduction**: 35 lines
- All 1060 tests pass, all edge cases verified

**Commit:** `9a66db0` - refactor(fastpoly): simplify Projector/Unipotent to sort+backtrack algorithm

---

## Session: 2025-12-08 23:00

**Agent:** orchestrator
**Feature:** R003 (pauli-simplification)

### Actions
- Researched QMBCertify codebase for Pauli simplification reference (via lead-researcher)
- Found QMBCertify uses staged pipeline: `reduce1!` (sort) → `reduce3!` (remove pairs) → `reduce2!` (Pauli algebra) → `reduce3!` (cleanup)
- Refactored pauli.jl to use `sort!` API instead of bubble sort
- Replaced while loop with single linear pass through sorted site groups
- Algorithm: sort by site → for each site group, reduce all operators to one (or identity) with coefficient

### Key Changes (pauli.jl)
- Replaced bubble sort (nested for loops) with `sort!(word, alg=InsertionSort, by=_pauli_site)`
- Replaced while-changed loop with single linear pass
- For each site group: accumulate Pauli products into single operator or identity
- Track complex phase coefficient throughout reduction

### Impact
- **Lines removed**: 45
- **Lines added**: 38
- **Net reduction**: 7 lines
- All 1060 FastPolynomials tests pass
- Cleaner, more maintainable code following same pattern as Projector/Unipotent

### Next Steps
- Continue with remaining review items (R004-R015)
- R003 is now complete - Pauli simplification uses sort!-based algorithm

**Commit:** `805efdf` - refactor(fastpoly): simplify Pauli simplification to sort!-based algorithm

---

## Session: 2025-12-08 23:30

**Agent:** orchestrator
**Feature:** R004 (fermionic-simplification)

### Actions
- Read fermionic.jl (385 lines) - Generalized Wick's Theorem implementation
- Tested all algebraic rules: CAR {aᵢ, aⱼ†} = δᵢⱼ, {aᵢ, aⱼ} = 0, {aᵢ†, aⱼ†} = 0
- Verified nilpotency detection (iszero function)
- Analyzed algorithmic complexity

### Findings

**Algorithm Structure (4-step Wick's Theorem):**
1. `_find_valid_contractions` - O(n²) scan for annihilation-creation pairs
2. `_generate_nonoverlapping_combinations` - Recursive enumeration of non-overlapping subsets
3. `_compute_normal_ordered_term` - Normal ordering with permutation parity
4. `_combine_like_terms_fermi` - Dict-based term aggregation

**Correctness Verified:**
| Test Case | Expected | Actual | Status |
|-----------|----------|--------|--------|
| a₁ a₁† | 1 - a₁† a₁ | ✓ | PASS |
| a₁ a₁ | 0 (nilpotent) | ✓ | PASS |
| a₁ a₂ + a₂ a₁ | 0 (anticommute) | ✓ | PASS |
| a₁† a₂† + a₂† a₁† | 0 (anticommute) | ✓ | PASS |
| a₁ a₂ a₃ a₃† a₂† a₁† | 8 terms | ✓ | PASS |

**Complexity Analysis:**
- For k valid contractions: O(2^k) combinations in worst case
- Typical pattern `a₁ a₁† a₂ a₂† ... aₙ aₙ†`: n contractions → 2^n terms
- Performance after JIT: ~0.02ms for n=10 operators
- Nilpotency early-exit via `iszero()` check (~0.04μs)

**Helper Functions Verified:**
- `_permutation_parity`: Cycle counting algorithm correct (sign = (-1)^(n-c))
- `_contraction_sign`: Swap counting for bringing pairs together
- `_is_creation`, `_fermi_mode`: Inline helpers for signed encoding

**Test Coverage Gap Identified:**
- No dedicated FermionicAlgebra tests in test/fastpoly_test/simplify.jl
- Only doctest examples in fermionic.jl itself
- Recommendation: Add property-based tests for anticommutation relations

### Outcome
- R004 PASS: Wick's theorem implementation is mathematically correct
- Algorithm is exponential in number of contractions (inherent to Wick's theorem)
- Performance is acceptable for typical use cases (<20 operators)

### Next Steps
- Continue with R005 (Bosonic rook number algorithm)

**Commit:** (no code changes - review only)

---

## Session: 2025-12-08 23:45

**Agent:** orchestrator
**Feature:** R005 (bosonic-simplification)

### Actions
- Read bosonic.jl (483 lines) - Rook number algorithm for normal ordering
- Verified CCR: [c, c†] = 1 (c c† = c† c + 1)
- Tested Ferrers board construction and rook number computation
- Analyzed multi-mode vs single-mode complexity

### Findings

**Algorithm Structure (4-step Rook Number Method):**
1. `group_by_mode!` - O(n log n) stable sort to partition by mode
2. `build_ferrers_board` - O(K) per mode, count creations after each annihilation
3. `compute_rook_numbers` - O(K²) DP for K = operators per mode
4. `expand_and_construct` - O(∏ᵢ |terms_i|) Cartesian product

**Correctness Verified:**
| Test Case | Expected | Actual | Status |
|-----------|----------|--------|--------|
| c c† | c† c + 1 | ✓ | PASS |
| c† c | c† c | ✓ | PASS |
| c c (not nilpotent) | c c | ✓ | PASS |
| c₁ c₂† | c₂† c₁ | ✓ | PASS |
| c c c† c† | c†² c² + 4 c† c + 2 | ✓ | PASS |
| c₁ c₂ c₁† c₂† | 4 terms | ✓ | PASS |

**Rook Number Formula Verified (Eq. 1.40):**
- For word ω with m creations, n annihilations:
- ω = Σₖ rₖ(Bω) (a†)^{m-k} a^{n-k}
- Board: cell (i,j) included if creation i AFTER annihilation j
- Example: c c c† c† → board=[2,2] → rook=[1,4,2] → c†² c² + 4 c† c + 2 ✓

**Complexity Analysis:**
| Pattern | Result Size | Time (post-JIT) |
|---------|-------------|-----------------|
| n modes × (c c†) | 2^n terms | O(ms) |
| c^n (c†)^n single mode | n+1 terms | O(μs) |

**Key Insight:** Rook numbers give polynomial-time per mode (n+1 terms max), but multi-mode expansion is exponential in number of modes. This is inherent to the math.

**Test Coverage Gap:**
- No dedicated BosonicAlgebra tests in test/fastpoly_test/simplify.jl
- Only doctest examples in bosonic.jl itself

### Outcome
- R005 PASS: Rook number algorithm correctly implements arXiv:quant-ph/0507206
- Complexity claims in TODO are accurate: polynomial per mode, exponential across modes
- Well-documented with mathematical references

### Next Steps
- Continue with R008 (Canonicalization) - next in priority

**Commit:** (no code changes - review only)

---

## Session: 2025-12-09 00:00 - Checkpoint

**Agent:** orchestrator + qa-gatekeeper
**Feature:** R012 (test-coverage) - partial

### Actions
- Completed R004 review (Fermionic Wick's theorem)
- Completed R005 review (Bosonic rook numbers)
- Identified test coverage gaps for Fermionic and Bosonic algebras
- Delegated to qa-gatekeeper to add comprehensive simplification tests
- qa-gatekeeper added 110 new tests (38 Fermionic + 45 Bosonic + 27 related)

### Outcome
- R004 PASS: Fermionic Wick's theorem mathematically correct
- R005 PASS: Bosonic rook number algorithm correct (arXiv:quant-ph/0507206)
- Test coverage gap addressed: 1031 → 1141 tests

### Current State
- Features completed: 8/15 (R001-R007, R009)
- Features pending: 7 (R008, R010-R015)
- R012 (test-coverage) partially addressed by new tests

### Next Steps
- Continue with R008 (Canonicalization algorithms)
- Complete remaining review items (R010-R015)

**Commit:** `10ceae7` - test(fastpoly): add comprehensive FermionicAlgebra and BosonicAlgebra simplification tests

---
