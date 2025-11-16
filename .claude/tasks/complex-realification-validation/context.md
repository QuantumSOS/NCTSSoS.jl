# Task: Complex SDP Realification Validation

## Objective
Investigate whether the complex moment problem to real sum-of-squares translation in NCTSSoS.jl is implemented correctly by comparing it with the theoretical approach in arxiv.org/2307.11599.

## Background
The user suspects that the realification of complex SDP might be done incorrectly in this codebase. We need to:

1. Understand how NCTSSoS.jl currently translates complex moment problems into real sum of squares problems
2. Understand the correct approach from the paper: https://arxiv.org/pdf/2307.11599
3. Identify any discrepancies between the two approaches
4. If discrepancies exist, propose how to adapt the paper's approach to work in the codebase

## Project Context
- **Project**: NCTSSoS.jl - Julia package for solving sparse noncommutative polynomial optimization
- **Branch**: ys/gns-application
- **Related work**: Currently working on GNS reconstruction examples (Mermin Square)
- **Concern**: The realification process may affect correctness of complex moment problems

## Scope of Investigation

### Part 1: Codebase Analysis
Understand the current implementation:
- Locate code that handles complex polynomials and moment matrices
- Understand how complex variables are converted to real variables
- Identify how complex constraints are translated to real constraints
- Find how the moment matrix is constructed for complex problems
- Understand the relationship between complex and real basis elements

### Part 2: Paper Analysis
From https://arxiv.org/pdf/2307.11599:
- Understand the correct theoretical approach for complex to real translation
- Identify the key mathematical transformations
- Note any specific requirements or constraints
- Understand how complex positive semidefiniteness translates to real constraints

### Part 3: Comparison and Validation
- Compare the codebase implementation with the paper's approach
- Identify any mathematical discrepancies
- Assess the impact of any differences found
- Determine if the current implementation is incorrect or just different

### Part 4: Adaptation Strategy (if needed)
If discrepancies are found:
- Propose specific changes to align with the paper's approach
- Identify which functions/modules need modification
- Consider backward compatibility concerns
- Outline implementation steps

## Key Questions to Answer
1. How does the codebase currently perform realification?
2. What is the correct approach according to the paper?
3. Are there mathematical errors or just different conventions?
4. If there are errors, what is the minimal set of changes needed?
5. Will fixing this affect existing examples and tests?

## Expected Output
A detailed analysis document (.claude/tasks/complex-realification-validation/plan.md) containing:
- Summary of current implementation approach
- Summary of paper's approach
- Detailed comparison highlighting discrepancies
- Assessment of correctness/incorrectness
- If needed: adaptation plan with specific implementation steps
- References to specific code locations (file:line format)

## Notes
- This is a research and analysis task - DO NOT modify any code yet
- Focus on understanding and documenting first
- Be mathematically rigorous in the comparison
- Provide concrete examples if discrepancies are found

---

## Research Findings (Completed)

**Date**: 2025-11-10

**Result**: ✅ **No mathematical errors found - implementation is correct!**

### Key Findings

1. **Codebase Uses Efficient Formulation**: The NCTSSoS.jl implementation in `src/sos_solver.jl` (lines 92-146) correctly implements the **efficient dual-based complex-to-real SDP reformulation** (PSDP-R') from the paper.

2. **Matches State-of-the-Art**: The approach matches the 2023 paper by Jie Wang (arxiv:2307.11599) which is the most efficient known reformulation.

3. **No Extra Constraints**: The implementation correctly avoids the n(n+1) extra affine constraints that plague the standard formulation, achieving better computational efficiency.

4. **Block Structure**: The code properly uses:
   - `H_R = X1 + X2` (real part)
   - `H_I = X3 - X3^T` (imaginary part)
   from the 2n×2n real PSD matrix `X = [X1 X3^T; X3 X2]`

### Recommendations

While mathematically correct, some improvements for clarity:
1. Add citations and documentation explaining which formulation is used
2. Clarify the 4-tuple constraint assembly loop (lines 138-140)
3. Add validation tests comparing against reference complex SDP solvers

**Full analysis available in**: `.claude/tasks/complex-realification-validation/plan.md`
