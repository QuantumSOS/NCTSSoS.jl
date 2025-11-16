# Annotated Source Files - Moment Matrix Extraction

This directory contains annotated versions of all files modified during the moment matrix extraction implementation. Each file includes detailed inline annotations explaining:

- **WHY** changes were made (rationale and design decisions)
- **HOW** changes integrate with existing code
- **WHAT** trade-offs were considered
- **WHERE** the change fits in the overall architecture

---

## Directory Structure

```
.claude/tasks/moment-matrix-extraction/summary/
├── SUMMARY.md                    # Comprehensive implementation summary
├── ANNOTATED_FILES.md           # This file
└── annotated/
    ├── moment_extraction.jl      # NEW: Core data structures and support building
    ├── moment_extraction_api.jl  # NEW: Public API and extraction functions
    ├── interface.jl              # MODIFIED: PolyOptResult with moment support
    ├── moment_solver.jl          # MODIFIED: Primal formulation integration
    └── sos_solver.jl             # MODIFIED: Dual formulation integration
```

---

## Annotated Files

### 1. Core Implementation (NEW)

#### `annotated/moment_extraction.jl`
**Lines**: 190 total (~140 implementation + ~50 annotations)

**Purpose**: Core data structures and support building logic

**Key Annotations**:
- Hierarchical vs flat structure design decision
- Type parameter flexibility (M vs M2) rationale
- Only processing moment matrix blocks (not localizing matrices)
- Direct indexing strategy for O(1) lookup
- Matching constraint construction pattern exactly

**Major Components**:
- `BlockSupport` struct - Maps (i,j) to dual variable index
- `CliqueSupport` struct - Contains blocks for one clique
- `MomentSupport{M}` struct - Complete hierarchical structure
- `build_moment_support()` - Constructs mapping during solving

**Critical Insights**:
- Why only store indices, not monomials
- How canonicalization affects type parameters
- Why dictionary lookup vs binary search
- How structure mirrors problem sparsity

---

#### `annotated/moment_extraction_api.jl`
**Lines**: 195 total (~150 implementation + ~45 annotations)

**Purpose**: Public API for moment matrix extraction

**Key Annotations**:
- Separate function vs embedded extraction decision
- Primal vs dual formulation branching logic
- Negative sign convention for dual variables
- Direct indexing reconstruction algorithm
- Simplification vs NCTSSOS comparison

**Major Components**:
- `get_moment_matrices()` - Main public API
- `extract_dual_variables()` - For dual SOS formulation
- `extract_primal_variables()` - For primal moment formulation
- `reconstruct_moment_matrices()` - Pure indexing reconstruction

**Critical Insights**:
- Why dual variables have negative sign
- How NCTSSOS recomputes vs our direct indexing
- Mathematical relationship: dual variables = moment values
- Performance characteristics: O(n²) pure array indexing

---

### 2. Integration Points (MODIFIED)

#### `annotated/interface.jl`
**Lines**: 151 total (~10 lines modified, rest context)

**Purpose**: Main solving interface with moment support

**Key Annotations**:
- Why M2 type parameter added (canonicalization type flexibility)
- How moment_support field enables extraction
- Why is_dual flag needed (formulation detection)
- Integration with moment_relax/sos_dualize workflow
- Linear workflow: solve → extract as separate steps

**Major Changes**:
- `PolyOptResult{T,P,M,M2}` - Added M2 parameter and two fields
- `cs_nctssos()` - Updated to pass moment_support and is_dual
- `cs_nctssos_higher()` - Same modifications

**Critical Insights**:
- No breaking API changes (signature unchanged)
- Minimal overhead (single boolean, ~1KB support structure)
- Type parameter necessary for state polynomials
- How dualize parameter flows through to result

---

#### `annotated/moment_solver.jl`
**Lines**: 116 total (~15 lines modified, rest context)

**Purpose**: Primal moment formulation with moment support

**Key Annotations**:
- Why return tuple (MomentProblem, MomentSupport)
- How global_support construction works
- Why canonicalization needed
- Matching constraint construction pattern
- Performance: reusing basis products

**Major Changes**:
- `moment_relax()` - Returns tuple instead of just MomentProblem
- Added global_support construction via canonicalization
- Added build_moment_support() call

**Critical Insights**:
- Constraint order must match global_support order
- Basis products already computed for constraints
- Why canonicalize after collecting monomap keys
- Pattern: expval(_neat_dot3(row, mono, col))

---

#### `annotated/sos_solver.jl`
**Lines**: 173 total (~20 lines modified, rest context)

**Purpose**: Dual SOS formulation with moment support

**Key Annotations**:
- Why use symmetric_basis (already canonicalized)
- Dual variables of equality constraints = moment values
- Type flexibility: M vs M2 parameters
- Complex polynomial support
- Mathematical foundation of dualization

**Major Changes**:
- `sos_dualize()` - Returns tuple (SOSProblem, MomentSupport)
- Both real and complex versions modified
- Added build_moment_support() calls

**Critical Insights**:
- This is DEFAULT formulation (dualize=true)
- symmetric_basis already canonical, no extra work needed
- Dual variables extracted from fα_constraints
- Complex case uses real/imag parts separately

---

## Annotation Format

All annotations follow this consistent format:

```julia
# CHANGED: What changed (brief description)
# WHY: Rationale - link to task context/requirements
# [Optional sections based on relevance:]
# CONTEXT: Additional background from task plan
# IMPACT: How this affects other code or functionality
# ALTERNATIVES: Other approaches considered and why rejected
# TRADE-OFF: Any compromises made
# BUG FIX: If fixing an issue - what was wrong and how fixed
# CRITICAL: Highlight particularly important points
# DESIGN DECISION: Explain choice between alternatives
# ALGORITHM: Describe algorithmic approach
# PERFORMANCE: Performance characteristics
# TYPE FLEXIBILITY: Type parameter reasoning
# INTEGRATION: How it connects with existing code
# MATHEMATICAL FOUNDATION: Math behind the implementation
```

---

## How to Read the Annotations

### For Understanding Design Decisions

Look for `DESIGN DECISION`, `ALTERNATIVES`, and `TRADE-OFF` annotations to understand why specific approaches were chosen.

**Example**:
```julia
# DESIGN DECISION: Hierarchical vs Flat Structure
# CHOSEN: Hierarchical (cliques → blocks → dual_indices)
# RATIONALE: Semantic clarity, self-documenting, mirrors problem sparsity
# ALTERNATIVE: NCTSSOS uses flat global support with algorithmic lookup
# TRADE-OFF: Slightly more complex structure, but much clearer intent
```

### For Understanding Integration

Look for `INTEGRATION`, `IMPACT`, and `CONTEXT` annotations to see how changes fit into the existing codebase.

**Example**:
```julia
# INTEGRATION: build_moment_support() matches constraint construction pattern
# IMPACT: Minimal overhead (~1ms), reuses existing computation
# CONTEXT: Called during moment_relax/sos_dualize, before optimization
```

### For Understanding Algorithms

Look for `ALGORITHM`, `PERFORMANCE`, and `CRITICAL` annotations to understand computational approaches and efficiency.

**Example**:
```julia
# ALGORITHM: Direct indexing via precomputed dual_indices
# PERFORMANCE: O(n²) pure array indexing - no recomputation!
# CRITICAL: All hard work done during support construction; extraction is trivial
```

### For Understanding Type System

Look for `TYPE FLEXIBILITY`, `TYPE PARAMETER`, and related annotations to understand Julia type system usage.

**Example**:
```julia
# TYPE FLEXIBILITY: Separate M and M2 parameters
# WHY: Canonicalization changes type (NCStateWord → StateWord)
# IMPACT: Enables state polynomial problems to work correctly
```

---

## Comparison with NCTSSOS

Throughout the annotations, comparisons with NCTSSOS implementation are provided to highlight:

- **Algorithmic differences**: Hierarchical structure vs flat support
- **Extraction approach**: Direct indexing vs recomputation + lookup
- **Code clarity**: Self-documenting structure vs implicit patterns
- **Performance**: Precomputed mappings vs runtime computation

**Key Insight**: NCTSSoS.jl achieves simpler extraction code by building the support structure during solving rather than on-demand during extraction.

---

## References to Task Documentation

Annotations reference these task documents:

- **context.md**: Overall task objectives, constraints, success criteria
- **plan.md**: Detailed implementation plan with 11 sections
- **integration-notes.md**: Integration challenges and solutions
- **primal-dual-fix.md**: Fixing formulation detection issues
- **fastpolynomials-api-research.md**: Understanding monomial operations

Cross-references help connect implementation choices back to planning decisions.

---

## Usage

### Reading the Annotations

1. Start with `SUMMARY.md` for overall picture
2. Read `annotated/moment_extraction.jl` for core abstractions
3. Read `annotated/moment_extraction_api.jl` for public API
4. Read integration files to see how it all connects

### Understanding Specific Changes

Use annotations to:
- Trace design decisions back to requirements
- Understand rationale for specific implementations
- Learn trade-offs and alternatives considered
- See how changes integrate with existing code

### Code Review

Annotations provide context for reviewers:
- Why specific approaches were chosen
- What alternatives were considered
- How changes affect system behavior
- Where potential issues might arise

---

**Total Annotated Lines**: ~900 lines across 5 files
**Annotation Density**: ~20-30% annotation comments
**Focus**: WHY over WHAT (diffs show WHAT, annotations explain WHY)
