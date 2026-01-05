---
active: true
iteration: 1
max_iterations: 100
completion_promise: "TYPE_SYSTEM_REVAMP_COMPLETE"
started_at: "2026-01-05T16:24:29Z"
---


  ## Type System Revamp (NCTSSoS.jl v2)

  Implement type system revamp per .memory/tasks/type-system-revamp/*.md

  ### PHASE 0: DISABLE ALL TESTS FIRST
  Before ANY code changes, comment out ALL test includes in test/runtests.jl. Tests will break immediately due to fundamental type changes.

  ### Implementation Order
  1. Internal helpers (_validate_*_word!, update _simplify_*_word! return types)
  2. AbstractMonomial{A,T} hierarchy
  3. PauliMonomial{T} (new file: src/types/pauli_monomial.jl)
  4. PhysicsMonomial{A,T} (new file: src/types/physics_monomial.jl)
  5. Monomial{A,T} constructor validation (BREAKING)
  6. Multiplication dispatch changes (BREAKING)
  7. Basis construction (get_ncbasis returns wrapper types)
  8. Moment matrix integration
  9. State polynomial restriction (drop Pauli/Fermi/Boson support)
  10. Type conversion (Polynomial from wrappers)

  ### Test Strategy
  - Tests disabled in Phase 0
  - Re-enable ONE test file at a time after implementing relevant phase
  - ADAPT tests to new return types (document: # ADAPTED: old → new)
  - DROP tests explicitly with reason (document: # DROPPED: reason)
  - ADD new tests for new types/behaviors

  ### Tests to DROP (document each):
  - State poly tests for Pauli/Fermi/Boson (no longer supported)
  - Tests expecting invalid Monomial construction (now throws)

  ### Per-Phase Workflow
  1. Read plan.md for phase details
  2. Implement phase
  3. Re-enable relevant tests, adapt as needed
  4. Run: make test-polynomials
  5. Fix failures
  6. Commit with: git commit -m 'type-revamp: phase N description'
  7. Update .memory/tasks/type-system-revamp/plan.md (check off completed items)

  ### Completion Criteria
  - All 10 phases complete (checkboxes in plan.md checked)
  - make test-polynomials passes
  - make test passes
  - All dropped tests documented with reason
  - New types have tests

  Output <promise>TYPE_SYSTEM_REVAMP_COMPLETE</promise> when ALL criteria met.
  
