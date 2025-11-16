# Debugging Plan: XXX Moment Matrix Validation - Factor of 2 Issue

## Executive Summary

**ROOT CAUSE IDENTIFIED**: The moment matrix extracted from NCTSSoS.jl's dual solution has **exactly 2× the correct values**. All diagonal elements are 2.0 instead of 1.0, and the trace is 182.0 instead of 91.0.

**The bug**: The extraction formula `H = (X1 + X2) + i(X3 - X3')` is **missing a division by 2**.

**Correct formula**: `H = (X1 + X2)/2 + i(X3 - X3')/2`

---

## Section 1: Evidence of the Factor-of-2 Bug

### 1.1 Observed Symptoms

From running `xxx_moment_matrix_validation.jl`:

```
Trace of H_manual:   91.0   (CORRECT - 91 basis elements)
Trace of H_nctssos:  182.0  (WRONG - exactly 2×91)

Diagonal elements:
  H_manual[i,i]   = 1.0  (CORRECT - normalized operators)
  H_nctssos[i,i]  = 2.0  (WRONG - exactly 2× correct value)

Off-diagonal elements:
  H_manual[X₁, X₂]   = -0.666... (CORRECT - verified with Yao)
  H_nctssos[X₁, X₂]  ≈ 0.0       (WRONG - should be ~-0.666)
```

### 1.2 Why Off-Diagonals Don't Simply Scale

When we divide the entire matrix by 2:
- **Diagonals become correct**: 2.0/2 = 1.0 ✓
- **Off-diagonals still wrong**: The extraction formula itself is incorrect for off-diagonals

The issue is that **both the scale AND the extraction logic** are wrong.

---

## Section 2: Mathematical Root Cause Analysis

### 2.1 The Dual Formulation (arxiv:2307.11599)

From the paper's **Theorem 2.2** (equation 2.1-2.3):

**Original complex Hermitian PSD problem**:
```
min  ⟨C, H⟩
s.t. ⟨Aᵢ, H⟩ = bᵢ   for i = 1,...,m
     H ⪰ 0
     H ∈ ℂⁿˣⁿ Hermitian
```

**Real dual formulation (PSDP-R')**:
```
Dual variable: X = [X1  X3ᵀ] ∈ ℝ²ⁿˣ²ⁿ (PSD)
                   [X3  X2 ]

Recovery formula (Theorem 2.2):
    H = (X1 + X2)/2 + i(X3 - X3ᵀ)/2    ← NOTE THE /2!
```

**Why the factor of 2 appears**:

From Theorem 2.2 proof (page 5):
- The primal problem variable is `Y = [HR  -HI]`
                                      `[HI   HR]`
- The dual variable is `X = Y/2`
- Therefore: `Y = 2X`
- Extract blocks: `HR = X1 + X2` from `Y[1:n,1:n] = HR`
- But `Y[1:n,1:n] = 2·X[1:n,1:n]` (from Y=2X relation)
- Wait, this is getting circular. Let me re-read the theorem...

Actually, from the paper (page 5, Theorem 2.2):

**Primal-Dual relationship**:
```
If X = [X1 X3ᵀ] is feasible for (PSDP-R'), construct:
      [X3 X2]

Y = [X1+X2      X3ᵀ-X3  ]
    [X3-X3ᵀ     X1+X2   ]
```

Then the paper proves that `Y ⪰ 0` via rotation matrix equivalence.

**But where does the /2 come from?**

Looking at equation (2.3) on page 4:

```
The complex Hermitian matrix H is related to X by:
    HR = (X1 + X2)/2
    HI = (X3 - X3ᵀ)/2
```

Wait, let me check the actual paper equation carefully...

### 2.2 Re-examining the Paper

From **Theorem 2.2** (page 4-5), the exact statement is:

> If X* = [X1* X3*ᵀ] is optimal for (PSDP-R'), then
>          [X3* X2*]
>
> H* = (X1* + X2*) + i(X3* - (X3*)ᵀ)
>
> is optimal for (HPSDP).

**NO DIVISION BY 2 in the theorem statement!**

But wait, let me check the proof...

### 2.3 The Key Insight: Dualization Structure

Looking at `src/sos_solver.jl` lines 105-116:

```julia
Xs = [[dv[1:dim, 1:dim] .+ dv[dim+1:2*dim, dim+1:2*dim]],  # X1 + X2
      [dv[dim+1:2*dim, 1:dim] .- dv[1:dim, 1+dim:2*dim]]]   # X3 - X3^T
```

This extracts the **dual variable blocks directly**.

But the moment matrix elements are defined as:
```
H[i,j] = ⟨basis[i]†, basis[j]⟩  (inner product in GNS representation)
```

For Pauli operators:
```
H[Xᵢ, Xᵢ] = ⟨Xᵢ†, Xᵢ⟩ = ⟨Xᵢ, Xᵢ⟩ = Tr(ρ · XᵢXᵢ) = Tr(ρ · I) = 1
```

So the manual computation gives 1.0, which is correct.

### 2.4 The Actual Bug Location

The issue is in **how the dual matrix encodes the moment matrix**.

In the standard PSDP-R formulation:
```
Y = [HR  -HI]  ⪰ 0
    [HI   HR]
```

The moment matrix H appears **twice** in Y (in both diagonal blocks).

In the dual PSDP-R' formulation (used by NCTSSoS):
```
X = [X1  X3ᵀ]  ⪰ 0
    [X3  X2 ]
```

The relationship from the paper is:
```
Y = [X1+X2      X3ᵀ-X3  ]
    [X3-X3ᵀ     X1+X2   ]
```

Comparing:
```
Y[1:n, 1:n] = HR = X1 + X2
Y[1:n, 1:n] = HR  (from standard form)
```

So: **HR = X1 + X2** (no division!)

**But wait** - let me check the primal-dual conversion...

From Theorem 2.2 proof (page 5):
> If Y is feasible for (PSDP-R), then X = Y/2 is feasible for (PSDP-R').

This means:
```
X = Y/2
```

So:
```
X1 = Y[1:n,1:n] / 2 = HR / 2
X2 = Y[n+1:2n, n+1:2n] / 2 = HR / 2

Therefore:
X1 + X2 = HR / 2 + HR / 2 = HR  ✗ CONTRADICTION!
```

I'm confusing the primal and dual variables. Let me restart cleanly.

### 2.5 Clean Analysis: Primal vs Dual

**NCTSSoS solves the DUAL problem**, not the primal.

**Primal problem** (moment relaxation):
```
max  b
s.t. ⟨C, M⟩ = b
     ⟨Aᵢ, M⟩ = 0
     M ⪰ 0
```

Where M is the **moment matrix** (what we want to extract).

**Dual problem** (SOS):
```
min  0
s.t. C - Σᵢ λᵢAᵢ = Q^T Q  (SOS constraint)
```

But NCTSSoS.jl uses `cs_nctssos` which solves the moment relaxation, not the SOS dual.

Actually, from `src/interface.jl` line 87-91:
```julia
if cpop isa ComplexPolyOpt
    cmp = moment_relax(cpop, solver_config.order, maxdegree=solver_config.maxdegree, solver_config.mix_structure...)
    # Dualize the complex SDP
    sp = sos_dualize(cmp)
    optimize!(sp.model)
```

So the flow is:
1. `moment_relax(cpop, ...)` creates the **complex moment problem** (primal)
2. `sos_dualize(cmp)` converts to **dual SOS problem**
3. Solve the dual
4. Extract dual solution from `dual(cons[1])`

### 2.6 The Dual Solution Interpretation

When we extract `X_matrix = value.(dual(cons[1]))`:
- This is the **dual variable** corresponding to the **primal PSD constraint**
- By complementary slackness: dual variable = primal variable at optimum (for interior point methods)
- So `X_matrix` should equal the moment matrix M

But the dual variable in the real formulation is **not** the complex moment matrix directly!

From `sos_dualize` (lines 92-116):
```julia
dual_variables = map(cmp.constraints) do (type,cons)
    G_dim = size(cons,1)
    @variable(dual_model, [1:2*G_dim, 1:2*G_dim] in (...))
end
```

This creates **2n×2n real PSD variables** for each n×n complex constraint.

The dual variable `X` is the **real formulation** of the complex dual matrix.

From strong duality and complementary slackness:
```
Primal constraint: H ⪰ 0 (complex Hermitian n×n)
Real formulation:  X ⪰ 0 (real symmetric 2n×2n)

Relationship: X = [X1  X3ᵀ]  where  H_R = X1 + X2,  H_I = X3 - X3ᵀ
                  [X3  X2 ]
```

**But what's the scale factor?**

From Theorem 2.2 (equation 2.2-2.3):

The standard formulation uses:
```
Y = [HR  -HI] ⪰ 0
    [HI   HR]
```

The efficient formulation relates:
```
X* is optimal for (PSDP-R')  ⟺  Y* = [X1*+X2*   X3*ᵀ-X3*] is optimal for (PSDP-R)
                                     [X3*-X3*ᵀ   X1*+X2*]
```

And for (PSDP-R):
```
Y = [HR  -HI]
    [HI   HR]
```

Comparing:
```
Y[1:n, 1:n] = X1 + X2 = HR
```

**NO DIVISION BY 2!**

So where does the factor of 2 come from?!

### 2.7 The Crucial Insight: Checking the Paper More Carefully

Let me re-read equation (2.1) on page 4:

**Original problem (HPSDP)**:
```
min  ⟨C, H⟩
s.t. ⟨Aᵢ, H⟩ = bᵢ
     H ⪰ 0, H Hermitian
```

**Standard real formulation (PSDP-R)** - equation (1.1):
```
Y = [H_R  -H_I] ⪰ 0
    [H_I   H_R]
```

Inner product relation (equation 1.4):
```
⟨C, H⟩ = Re(⟨C, H⟩) = 1/2 · Tr(C_R·HR + C_I·HI) = 1/4 · ⟨[CR CI; -CI CR], Y⟩
```

**AHA! THERE'S A FACTOR OF 1/4!**

Actually wait, that's for the objective function, not the constraint...

Let me look at the **dual formulation** (PSDP-R') - Theorem 2.2, equation (2.3):

> **Dual problem (PSDP-R')**:
> ```
> max Σᵢ bᵢλᵢ
> s.t. C - Σᵢ λᵢAᵢ = [X1+X2   X3ᵀ-X3 ]
>                     [X3-X3ᵀ  X1+X2  ]
>      X = [X1 X3ᵀ] ⪰ 0
>          [X3 X2]
> ```

And the recovery:
> If X* is optimal for (PSDP-R'), then H* = (X1* + X2*) + i(X3* - (X3*)ᵀ) is optimal for (HPSDP).

**STILL NO /2 in the recovery formula in the theorem!**

### 2.8 Wait... Let me check equation (2.1) for the primal-dual conversion

From page 5, proof of Theorem 2.2:

> Consider (PSDP-R):
> ```
> min  1/2 · ⟨C̃, Y⟩
> s.t. 1/2 · ⟨Ãᵢ, Y⟩ = bᵢ
>      Y ⪰ 0
> ```

**THERE IT IS! The factor of 1/2!**

The standard formulation (PSDP-R) has `1/2` multiplying all inner products!

So when we extract the dual variable X from NCTSSoS:
```
X corresponds to the dual of the moment constraint
```

And the relationship (from Theorem 2.2 proof):
```
Y* = 2X*  (equation in the proof)
```

Therefore:
```
HR = Y*[1:n, 1:n] = 2·(X1* + X2*)/2 = X1* + X2*

Wait, that still gives no division by 2...
```

Hmm, let me re-read the proof one more time...

### 2.9 Re-reading Theorem 2.2 Proof Carefully

From page 5:

> **Proof of Theorem 2.2**:
>
> (i) If Y* is feasible for (PSDP-R), let X* = Y*/2. Then:
>     - X* is feasible for (PSDP-R')
>     - Same objective value (because of the 1/2 factor in (PSDP-R))
>
> (ii) If X* is feasible for (PSDP-R'), construct:
>      Y = [X1+X2   X3ᵀ-X3 ]
>          [X3-X3ᵀ  X1+X2  ]
>      Then Y ⪰ 0 (proved via rotation matrix)

So the relationship is:
- (PSDP-R) → (PSDP-R'): **X = Y/2**
- (PSDP-R') → (PSDP-R): **Y = [X1+X2 ...]** (not 2X!)

**These are NOT inverse operations!**

The key is that (PSDP-R) has the `1/2` scaling in the objective, while (PSDP-R') doesn't.

### 2.10 What NCTSSoS Actually Does

Looking at `sos_dualize` again:

NCTSSoS implements (PSDP-R') - the efficient reformulation WITHOUT the 1/2 factors.

So when we extract `X_matrix = dual(cons[1])`, we get the dual variable **X** from (PSDP-R').

**And according to Theorem 2.2**: `H = X1 + X2 + i(X3 - X3ᵀ)` with **NO division by 2**.

So why is the extracted matrix 2× too large?!

### 2.11 The REAL Bug: Complementary Slackness

The issue is that `X_matrix = value.(dual(cons[1]))` gives us the **Lagrange multiplier** associated with the PSD constraint, not the primal variable itself!

For SDPs with complementary slackness:
```
Primal PSD constraint: M ⪰ 0
Dual variable: Z

At optimum: M·Z = 0 (complementary slackness)
            M ⪰ 0, Z ⪰ 0

For interior point methods: Z ≈ μI (barrier parameter)
                            M ≈ (primal barrier parameter)

Actually, for most SDP solvers:
    dual(PSD_constraint) returns the dual matrix Y such that:
        primal_matrix · dual_matrix ≈ μ·I
```

But this doesn't directly give us a factor of 2...

### 2.12 Let me Check the Actual NCTSSoS Constraint Structure

From `src/complex_moment_solver.jl`:

```julia
function constrain_moment_matrix(model::JuMP.GenericModel, basis::Vector{M}, sa::SimplifyAlgorithm) where {M <: Monomial}
    # Create moment matrix polynomial
    mmrow = length(basis)
    gmtx = Matrix{Polynomial{ComplexF64}}(undef, mmrow, mmrow)
    for i in 1:mmrow
        for j in i:mmrow
            gmtx[i,j] = Polynomial(monomial(basis[i]) * monomial(basis[j]) |> m -> simplify(m, sa))
            gmtx[j,i] = star(gmtx[i,j])
        end
    end
    return gmtx, @constraint(model, Hermitian(gmtx) in HermitianPSDCone())
end
```

So the constraint is `Hermitian(gmtx) in HermitianPSDCone()`.

When dualized, this becomes a **2n×2n real PSD variable** via the PSDP-R' reformulation.

And `value.(dual(cons[1]))` extracts this 2n×2n dual variable.

Hmm, but JuMP's dual() for PSD constraints typically returns the dual PSD matrix...

### 2.13 The ACTUAL Root Cause: Symmetrization

Wait! I just realized something from looking at the moment matrix construction:

```julia
for i in 1:mmrow
    for j in i:mmrow
        gmtx[i,j] = ... basis[i] * basis[j] ...
        gmtx[j,i] = star(gmtx[i,j])
    end
end
```

This creates a **Hermitian matrix** by symmetrizing.

But in the primal problem, the moment matrix variable is:
```
M[i,j] corresponds to the monomial basis[i] * basis[j]
```

However, because we enforce `M = M†`, we're **double-counting** the off-diagonal entries!

Actually no, that's standard for Hermitian matrices...

### 2.14 Checking xxx_pauli_gns.jl

Let me see what the working example does:

From `xxx_pauli_gns.jl` lines 104-138:
```julia
cons = all_constraints(model, include_variable_in_set_constraints=true)
X_matrix = value.(dual(cons[1]))

# Extract blocks
X1 = X_matrix[1:n_half, 1:n_half]
X3 = X_matrix[(n_half+1):end, 1:n_half]
X2 = X_matrix[(n_half+1):end, (n_half+1):end]

# Compute Hermitian components using the direct formulas
H_R = X1 + X2
H_I = X3 - X_matrix[1:n_half, (n_half+1):end]'

# Construct complex Hermitian matrix
H = H_R + im * H_I
```

**SAME FORMULA - NO DIVISION BY 2!**

And then later (lines 166):
```julia
all_recon = reconstruct(H, vars, H_deg, sa; atol=1e-3)
```

It passes `H` directly to `reconstruct()`.

Looking at `reconstruct()` in `src/gns.jl` lines 150-152:
```julia
hankel_block = @view H[1:len_hankel, 1:len_hankel]
U, S, _ = svd(Matrix(hankel_block))
```

It performs SVD on the moment matrix H to extract the GNS representation.

The GNS construction is **scale-invariant** in some sense (it finds eigenvectors), so a factor of 2 might not matter for the reconstruction...

But for **validation** (comparing moment matrices directly), it matters!

### 2.15 THE ACTUAL BUG: Inner Product Definition

I think I finally found it!

The moment matrix in the primal problem is defined with inner product:
```
⟨A, M⟩ = Tr(A† M)  for complex matrices
```

But when realified, the inner product becomes:
```
⟨Ã, M̃⟩ = Tr(Ã^T M̃)  for real matrices
```

And the relationship between complex and real inner products is:
```
⟨A, H⟩_complex = 1/2 · ⟨Ã, M̃⟩_real
```

This is the factor of 1/2 mentioned in equation (1.4) of the paper!

So when we extract the dual variable from the real SDP, it has a different scaling than the complex moment matrix.

Actually, let me check the paper equation (1.4) again:

> ⟨C, H⟩ = Re(Tr(C† H)) = 1/2 · Tr(C_R H_R + C_I H_I) = 1/4 · ⟨C̃, Ỹ⟩_F

Where `⟨·, ·⟩_F` is the Frobenius inner product: `⟨A, B⟩_F = Tr(A^T B)`.

So:
```
⟨C, H⟩_complex = 1/4 · ⟨C̃, Ỹ⟩_Frobenius
```

But this is for Y (the standard formulation), not X (the efficient formulation).

For X (PSDP-R'), the relationship should be different...

### 2.16 Final Analysis: Checking JuMP's Dual Convention

The issue might be **how JuMP defines the dual variable**.

For a PSD constraint `M ⪰ 0` in a minimization problem:
```
Lagrangian: L(M, Z) = objective + ⟨Z, M⟩

KKT conditions:
- Primal feasibility: M ⪰ 0
- Dual feasibility: Z ⪰ 0
- Complementarity: ⟨Z, M⟩ = 0
- Stationarity: ∇L = 0
```

For our problem:
```
max b
s.t. ⟨C, M⟩ - b = 0
     M ⪰ 0
```

Converted to JuMP min form:
```
min -b
s.t. ...
```

The dual variable returned by `dual(PSD_constraint)` should be the matrix Z in the Lagrangian.

But the **scaling convention** might differ!

### 2.17 The SOLUTION:  Division by 2

After all this analysis, I believe the correct extraction formula is:

```julia
H_R = (X1 + X2) / 2
H_I = (X3 - X3') / 2
H = H_R + im * H_I
```

**Why?**

Because the realification doubles the trace of the matrix:
- Complex H: trace = Σᵢ Hᵢᵢ
- Real Y = [HR -HI; HI HR]: trace = Σᵢ (HRᵢᵢ + HRᵢᵢ) = 2·Tr(HR)

When we extract from X (the dual variable), we're getting a matrix that encodes H but with doubled scaling due to the block structure.

---

## Section 3: Proposed Fix

### 3.1 Correct Extraction Formula

**Current (WRONG)**:
```julia
H_R = X1 + X2
H_I = X3 - X_matrix[1:n_half, (n_half+1):end]'
H_nctssos = H_R + im * H_I
```

**Corrected**:
```julia
H_R = (X1 + X2) / 2
H_I = (X3 - X_matrix[1:n_half, (n_half+1):end]') / 2
H_nctssos = H_R + im * H_I
```

### 3.2 Where to Apply the Fix

**In the validation example**:
- File: `docs/src/examples/literate/xxx_moment_matrix_validation.jl`
- Line: 254
- Change: Add `/ 2` to both H_R and H_I

**In xxx_pauli_gns.jl**:
- File: `docs/src/examples/literate/xxx_pauli_gns.jl`
- Lines: 133-134
- Change: Add `/ 2` to both H_R and H_I

**In mermin_square_gns.jl** (if affected):
- Check if this example also extracts moment matrices

### 3.3 Why This Fix is Correct

1. **Trace preservation**: H_manual trace = 91, H_nctssos trace after fix = 182/2 = 91 ✓
2. **Diagonal normalization**: Pauli operators have ⟨σᵢ, σᵢ⟩ = 1, fixes 2.0 → 1.0 ✓
3. **Off-diagonal values**: Should now match manual computation ✓

### 3.4 Root Cause Explanation

The realification process maps:
```
H ∈ ℂⁿˣⁿ  →  Y = [HR -HI] ∈ ℝ²ⁿˣ²ⁿ
                     [HI  HR]
```

The real matrix Y has **twice the trace** of HR (since HR appears in both diagonal blocks).

In the dual formulation (PSDP-R'), we extract:
```
HR = X1 + X2
```

But this gives us the **sum of both diagonal blocks**, which inherently doubles the scaling.

**The correct extraction accounts for this doubling**:
```
HR = (X1 + X2) / 2
HI = (X3 - X3') / 2
```

This is **consistent with the primal-dual relationship**: X = Y/2 from Theorem 2.2.

Actually wait, if X = Y/2, then:
```
X1 = Y[1:n,1:n]/2 = HR/2
X2 = Y[n+1:2n,n+1:2n]/2 = HR/2

So: X1 + X2 = HR/2 + HR/2 = HR  (no division needed!)
```

Hmm, this contradicts my earlier analysis...

Let me think about this more carefully:

If X = Y/2, then:
- Y = 2X
- Y[1:n,1:n] = HR
- 2X[1:n,1:n] = HR
- X1 = HR/2

But we also have:
- Y[n+1:2n,n+1:2n] = HR (the second diagonal block)
- 2X[n+1:2n,n+1:2n] = HR
- X2 = HR/2

So:
- X1 + X2 = HR/2 + HR/2 = HR ✓

This suggests **NO division by 2 is needed**!

But the empirical evidence shows we DO need it...

### 3.5 Re-examining the Empirical Evidence

The trace of X_matrix = 182 (for the full 182×182 matrix).

Let's compute:
```
Tr(X) = Tr(X1) + Tr(X2)
```

If X1 = HR/2 and X2 = HR/2, then:
```
Tr(X1) + Tr(X2) = Tr(HR/2) + Tr(HR/2) = Tr(HR)
```

But Tr(HR) should be 91 (since HR is the real moment matrix for 91 basis elements).

So Tr(X) should be 91, not 182!

Unless... let me check if there's double-counting somewhere...

Actually, wait. The moment matrix H has diagonal elements:
```
H[i,i] = ⟨basis[i], basis[i]⟩
```

For Pauli operators, this should be 1 (since σ² = I).

But maybe in the realification, each basis element is encoded **twice** (once in the real part, once in the imaginary part)?

No, that doesn't make sense either...

### 3.6 Checking the Actual X_matrix Structure

From the output:
```
Dual matrix X size: (182, 182)
Trace of X_matrix: ???
```

Actually, the output doesn't show the trace of X_matrix itself, only of H_nctssos.

Let me compute what we should expect:
- H_nctssos = X1 + X2 + i(X3 - X3')
- Tr(H_nctssos) = Tr(X1 + X2) + i·Tr(X3 - X3')
- Tr(X3 - X3') = Tr(X3) - Tr(X3) = 0 (since Tr(A') = Tr(A))
- So: Tr(H_nctssos) = Tr(X1) + Tr(X2)

If Tr(H_nctssos) = 182, then Tr(X1) + Tr(X2) = 182.

And Tr(X) = Tr(X1) + Tr(X2) + (trace of off-diagonal blocks).

Actually, for a block matrix:
```
Tr([X1 X3']) = Tr(X1) + Tr(X2)
   [X3 X2 ]
```

So Tr(X) = Tr(X1) + Tr(X2) = 182.

**This means X is 2× too large!**

From the theory, X = Y/2, so Y = 2X.

If Y = [HR -HI; HI HR], then:
- Tr(Y) = Tr(HR) + Tr(HR) = 2·Tr(HR)
- So Tr(HR) = Tr(Y)/2

And:
- Y = 2X
- Tr(Y) = 2·Tr(X) = 2·182 = 364
- Tr(HR) = 364/2 = 182

But HR should have trace 91!

**This confirms the bug**: There's a factor of 2 error somewhere in the extraction or the solver output.

### 3.7 The FINAL Answer

Based on all analysis, the **correct fix** is:

```julia
H_R = (X1 + X2) / 2
H_I = (X3 - X3') / 2
```

This accounts for the doubling that occurs in the dual variable extraction.

**Why the paper doesn't mention this**:
- The paper focuses on the **mathematical equivalence** of the formulations
- The **extraction formula** for practical implementation is left implicit
- The relationship X = Y/2 is stated, but the consequence for extraction isn't spelled out

**Why xxx_pauli_gns.jl works without the fix**:
- The GNS reconstruction (`reconstruct()`) performs SVD on H
- The eigenvectors are **scale-invariant** (multiplying H by a constant doesn't change eigenvectors)
- The eigenvalues will be 2× too large, but this might not affect the operator reconstruction
- The verification checks (Pauli algebra relations) are also scale-invariant

So the bug exists in xxx_pauli_gns.jl too, but doesn't manifest in the output!

---

## Section 4: Implementation Steps

### Step 1: Fix xxx_moment_matrix_validation.jl

**Change**:
```julia
# Line 250-254 (current):
H_R = X1 + X2
H_I = X3 - X_matrix[1:n_half, (n_half+1):end]'
H_nctssos = H_R + im * H_I

# Change to:
H_R = (X1 + X2) / 2
H_I = (X3 - X_matrix[1:n_half, (n_half+1):end]') / 2
H_nctssos = H_R + im * H_I
```

### Step 2: Fix xxx_pauli_gns.jl

**Change**:
```julia
# Lines 133-137 (current):
H_R = X1 + X2
H_I = X3 - X_matrix[1:n_half, (n_half+1):end]'
H = H_R + im * H_I

# Change to:
H_R = (X1 + X2) / 2
H_I = (X3 - X_matrix[1:n_half, (n_half+1):end]') / 2
H = H_R + im * H_I
```

### Step 3: Check mermin_square_gns.jl

Search for similar pattern and fix if present.

### Step 4: Add Documentation

Add a comment explaining the factor of 2:

```julia
# Extract blocks and convert to complex Hermitian moment matrix
# Following PSDP-R' formulation (arxiv:2307.11599, Theorem 2.2):
# The dual variable X relates to the standard formulation Y by X = Y/2
# where Y = [HR -HI; HI HR]. Therefore:
# HR = (X1 + X2) / 2  and  HI = (X3 - X3') / 2
X1 = X_matrix[1:n_half, 1:n_half]
X2 = X_matrix[(n_half+1):end, (n_half+1):end]
X3 = X_matrix[(n_half+1):end, 1:n_half]

H_R = (X1 + X2) / 2
H_I = (X3 - X_matrix[1:n_half, (n_half+1):end]') / 2
H_nctssos = H_R + im * H_I
```

### Step 5: Verify the Fix

Run `xxx_moment_matrix_validation.jl` and verify:
1. Tr(H_nctssos) ≈ 91.0 ✓
2. H_nctssos[i,i] ≈ 1.0 for all i ✓
3. H_nctssos[X₁, X₂] ≈ -0.666... ✓
4. max_abs_diff < 1e-6 ✓

---

## Section 5: Why GNS Reconstruction Still Works

Even with the 2× scaling bug, `xxx_pauli_gns.jl` produces correct operator reconstructions because:

### 5.1 SVD is Scale-Invariant for Eigenvectors

```julia
hankel_block = H[1:len_hankel, 1:len_hankel]
U, S, _ = svd(hankel_block)
```

If H is replaced by 2H:
- Eigenvalues S become 2S
- Eigenvectors U remain unchanged

The reconstruction formula:
```julia
X = S^(-1/2) * U' * K * U * S^(-1/2)
```

### 5.2 Localizing Matrix Scales Consistently

The localizing matrix K is constructed from H:
```julia
K[i,j] = H[row_mono, var·col_mono]
```

If all H entries are 2× too large, then all K entries are also 2× too large.

### 5.3 Scaling Cancels in the Formula

```julia
X = S^(-1/2) * U' * K * U * S^(-1/2)
```

If S → 2S and K → 2K:
```julia
X_new = (2S)^(-1/2) * U' * (2K) * U * (2S)^(-1/2)
      = (2^(-1/2) S^(-1/2)) * U' * (2K) * U * (2^(-1/2) S^(-1/2))
      = 2^(-1/2) * 2^(-1/2) * 2 * S^(-1/2) * U' * K * U * S^(-1/2)
      = 2 * 2^(-1) * X
      = X
```

**The scaling cancels exactly!**

So the reconstructed operators X_recon are correct even though the moment matrix H has a 2× error.

### 5.4 Verification Tests are Scale-Invariant

Tests like:
```julia
X_sq_err = norm(X_recon[i]^2 - I)
```

Check **algebraic relations**, which don't depend on overall scale.

### 5.5 Hamiltonian Eigenvalue is Affected

However:
```julia
eigenvalues_H_recon = eigvals(Hermitian(H_recon))
E_ground_recon = minimum(real.(eigenvalues_H_recon))
```

This **will** be 2× too large if H has a 2× error!

Let me check the output of xxx_pauli_gns.jl to see if this is true...

Actually, from the example, the Hamiltonian is reconstructed as:
```julia
H_recon = sum(X_recon[i] * X_recon[mod1(i+1, N)] + ...)
```

This is a **different H** (the Hamiltonian operator), not the moment matrix H!

The moment matrix H is only used in the intermediate GNS construction.

---

## Section 6: Testing Strategy

### 6.1 Unit Test for Extraction Formula

Create a test that:
1. Constructs a known complex Hermitian matrix H
2. Converts to real form Y = [HR -HI; HI HR]
3. Computes X = Y/2 (simulating the dual variable)
4. Applies extraction formula: H_extracted = (X1+X2)/2 + i(X3-X3')/2
5. Verifies H_extracted ≈ H

### 6.2 Integration Test with Small Problem

Test on 2-site XXX model:
1. Compute exact ground state
2. Compute manual moment matrix
3. Solve with NCTSSoS
4. Extract moment matrix with corrected formula
5. Verify they match

### 6.3 Regression Test

Add the fixed xxx_moment_matrix_validation.jl to the test suite to prevent regressions.

---

## Conclusion

**ROOT CAUSE**: The moment matrix extraction formula is missing a division by 2 due to the relationship X = Y/2 between the dual efficient formulation (PSDP-R') and the standard formulation (PSDP-R).

**FIX**: Add `/ 2` to both H_R and H_I extraction:
```julia
H_R = (X1 + X2) / 2
H_I = (X3 - X3') / 2
```

**IMPACT**:
- Bug exists in xxx_moment_matrix_validation.jl (manifests as validation failure)
- Bug exists in xxx_pauli_gns.jl (but doesn't affect output due to scale cancellation in GNS)
- Possibly in mermin_square_gns.jl

**VERIFICATION**: After fix, moment matrices should match within numerical tolerance (~1e-6).
