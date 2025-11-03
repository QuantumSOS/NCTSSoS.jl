# [Symmetry in Polynomial Optimization](@id symmetry)

As problem size grows, the PSD matrix whose each entry is a variable grows. It will be nice that some entries in the psd matrix is zero. That way, the PSD matrices will become block diagonal. For PSD matrices, block diagonal parts needs to be PSD as well. So essentially, we only need to treat the smaller matrices as variables.

Based on this idea, we use sparsity in the problem to disregard some entries in the larger PSD matrix so that we only treat smaller blocks as variables. This is a kind of relaxation as some smaller blocks.

On the other hand, the idea of solving smaller blocks as solution for larger PSD matrix can be made exact. Just as a basis transformation could make the pdf matrix block diagonal, apply a linear map on variables so that the resultant objective and constraints are the same, it may be possible to make the psd matrix contain explicit zero terms. This way, we could safely reduce the size of the psd matrix size.

## [Understanding Symmetry in Polynomial Optimization](@id understanding-symmetry)

### Defining Symmetry in Polynomial Optimization

A standard Polynomial Optimization Problem (POP) aims to minimize a function $f(x)$ subject to constraints defined by $g_j(x) \geq 0$. A POP is defined as a **symmetric POP** when the input data, $f$ and $g_j$, exhibit symmetry.

This symmetry is typically defined as **symmetry under a subgroup of $GL(n)$**.

A POP is considered **G-invariant** if both the objective function $f$ and the constraint functions $g_j$ remain unchanged under the action of a finite group $G$ with representation $\rho: G \to GL_n(\mathbb{R})$. The transformation of $f$ under a group element $g$ is defined as $f^g(x) := f(\rho(g)^{-1}x)$, and G-invariance requires $f^g = f$ and $g_j^g = g_j$ for all $g \in G$.

An example of a symmetric polynomial is $f(x) = x_1x_2 + x_2x_3 + x_3x_4 + x_4x_1$. Symmetry is recognized as one of the key structures (alongside correlative and term sparsity) that can be exploited in the Moment-SOS Hierarchy for POP.

## [Mathematical Foundations: Group Representations](@id mathematical-foundations)

The technique for exploiting symmetry relies heavily on group representation theory, which allows large matrices to be simplified (block-diagonalized) through a change of basis.

### Representations and Modules

A representation of a finite group $G$ is a finite-dimensional vector space $V$ along with a homomorphism $\rho: G \to GL(V)$ (the set of invertible transformations of $V$). The degree of the representation is $\dim(V)$.

The space $V$ is called a **G-module** if it satisfies certain properties related to group action and vector addition/scalar multiplication.

### G-Module Structure and Conditions

A vector space $V$ is called a **G-module** (where $G$ is a finite group) if the group action satisfies two key sets of properties:

#### 1. Group Identity and Associativity

The group action must respect the identity element and associativity of the group operation:
*   **Identity:** The identity element $1 \in G$ must act trivially on every vector $v \in V$: $1 \cdot v = v$.
*   **Associativity (Homomorphism Property):** The product of two group elements acting consecutively must equal the action of the group product: $g_1 \cdot (g_2 \cdot v) = (g_1g_2) \cdot v$.

#### 2. Linearity (Vector Space Compatibility)

The group action must be compatible with vector addition and scalar multiplication (linearity):
*   **Additivity:** The action of $g \in G$ distributes over vector addition: $g \cdot (v_1 + v_2) = g \cdot v_1 + g \cdot v_2$.
*   **Homogeneity:** The action of $g \in G$ commutes with scalar multiplication $\lambda$: $g \cdot (\lambda v) = \lambda g \cdot v$.

Essentially, a G-module is a vector space $V$ equipped with a group action (a homomorphism $\rho: G \to GL(V)$) where the group elements act as linear transformations (invertible transformations) of $V$. The function $\rho(g): V \to V$ defining the action is specified as $\rho(g) = v \mapsto g \cdot v$.

### Decomposition Theorems

1.  **Maschke's Theorem** states that any finite-dimensional G-module $V$ can be decomposed into a direct sum of **irreducible G-modules** $V = W_1 \oplus W_2 \oplus \cdots \oplus W_k$.
2.  If a matrix $Q$ commutes with all elements of the group representation ($\rho(g)Q = Q\rho(g)$ for all $g \in G$), **Schur's Lemma** and its corollary demonstrate that using a special basis—called a **symmetry adapted basis**—allows the matrix $Q$ to be block-diagonalized. This is summarized by the key message: **"A NICE BASIS MAKES MATRICES SIMPLER"**.

## [Symmetries in Semidefinite Programs (SDPs)](@id symmetries-in-sdps)

Symmetry principles are directly applied to the Semidefinite Programs (SDPs) that arise from polynomial optimization relaxations.

### G-Invariant SDPs

If an SDP is G-invariant, its optimal value remains the same as that of the "dense" SDP, allowing the solution to be restricted to **invariant matrices** (those $Q$ such that $Q = Q^g$, where $Q^g := \rho(g)Q\rho(g)^*$).

### Block Diagonalization

By utilizing the orthonormal symmetry adapted basis $T$, an invariant matrix $Q$ can be transformed such that $N = T^{-1}QT$ is block diagonal.

### Resulting SDP

The original optimization problem involving $Q$ can be transformed into an equivalent, smaller problem involving several smaller matrix blocks $Q_l$.

*Example:* A $3 \times 3$ matrix $Q$ invariant under $S_2$ can be block-diagonalized into a $2 \times 2$ block $Q_1$ and a $1 \times 1$ block $Q_2$.

## [Exploiting Symmetries in POP Hierarchies](@id exploiting-symmetries-in-pop-hierarchies)

For POPs, two main hierarchies are presented to incorporate symmetry, both aiming to dramatically reduce the size and complexity of the resulting SDP variables and matrix constraints.

### The First Hierarchy (Symmetric Adapted Hierarchy)

This approach constructs a hierarchy where the moments are adapted to the group symmetry:
*   It replaces dense moment matrices $M_r(y)$ with $M_r(y^G)$, where $y^G$ corresponds to the pseudo-moment variable associated with the **Reynolds Operator** $R_G(x^\alpha) := \frac{1}{|G|} \sum_{g \in G} f^g$.
*   *Benefit:* This greatly reduces the number of variables. For the cyclic group $C_4$ (acting on 4 variables), this approach yields a $5 \times 5$ moment matrix $M_1(y)$ that depends on only **4 variables** instead of 15.

### The Second Hierarchy (Block Diagonalization Hierarchy)

This hierarchy goes further than the first by applying block-diagonalization to the moment matrices themselves, leading to a much more efficient optimization problem.
*   The set of polynomials in $\mathbb{R}[x]$ of degree up to $r$ is viewed as a real G-module and is decomposed into complex irreducible components.
*   By calculating $M_r^G(y)$ using a symmetry adapted basis, the moment matrix is transformed into a block-diagonal matrix: $M_r^G(y) = \bigoplus_{l=1}^k M_{r}^{G,l}(y)$.
*   *Benefit:* For the $C_4$ example, this approach, while still using 4 variables, replaces the original $5 \times 5$ matrix structure (from the first hierarchy) with a **2 × 2 block plus three elementary (1 × 1) blocks**.

### Special Case: Symmetry under the Symmetric Group ($S_n$)

When the optimization problem is symmetric under the action of the Symmetric Group $S_n$, the resulting block structure of the moment matrix can be analyzed using concepts from representation theory, specifically **Young tableaux** and **Specht polynomials**. This specialized construction yields a specific block structure for the moment matrix. For an $S_3$ problem with relaxation order $r=2$, this complexity reduction leads to $4 \times 4$ and $3 \times 3$ block moment matrices instead of a single $10 \times 10$ matrix.

## [Benefits of Symmetry Exploitation](@id benefits-of-symmetry)

To summarize the impact of exploiting symmetry:
Exploiting symmetry allows for the replacement of large, dense Semidefinite Programs (SDPs) with a set of smaller, block-diagonal SDPs, dramatically reducing the number of variables and the computational cost. This transformation is analogous to using a prism to break down white light into its component colors, making a complex, singular problem much simpler to analyze by separating it into manageable, independent components.

### Why G-Module Structure is Critical

The G-module structure is fundamental because it underpins the entire mathematical framework used to simplify complex optimization problems by exploiting symmetry, primarily through **block diagonalization**.

1. **The Power of Decomposition (Maschke's Theorem)**
   *   Maschke's Theorem states that any finite-dimensional G-module $V$ (over $\mathbb{R}$ or $\mathbb{C}$, and where $G$ is a finite group) can be broken down into a **direct sum of irreducible G-modules**: $V = W_1 \oplus W_2 \oplus \cdots \oplus W_k$.
   *   An **irreducible** G-module $W$ is one that does not contain any non-trivial G-submodules (subspaces $W' \subseteq W$ such that $g \cdot w' \in W'$ for all $g \in G$ and $w' \in W'$), meaning it cannot be broken down further under the action of the group.

2. **Simplifying Matrices (Schur's Lemma)**
   *   When dealing with optimization problems (like SDPs) that possess G-invariance, we are often working with **invariant matrices** $Q$ that commute with all elements of the group representation ($\rho(g)Q = Q\rho(g)$).
   *   The combination of Maschke's Theorem and **Schur's Lemma** proves that if we use a special basis—called a **symmetry adapted basis**—this invariant matrix $Q$ must necessarily decompose into smaller blocks.
   *   This transformation simplifies the matrix $Q$ into a block-diagonal form $N = T^{-1}QT$, where the blocks $N_i$ correspond to the irreducible submodules identified in the decomposition.

The essential importance of defining the G-module with these properties is encapsulated in the statement: **"A NICE BASIS MAKES MATRICES SIMPLER"**. By establishing the structure of a G-module, mathematicians are guaranteed that they can switch to a symmetry adapted basis ($T$) that will decompose large, complex SDP matrix constraints ($Q$) into smaller, independent matrix blocks ($Q_l$), drastically reducing the computational cost of the optimization problem.

For instance, in the context of the Moment-SOS hierarchy, applying the G-module concept reduces a $5 \times 5$ moment matrix associated with the cyclic group $C_4$ into a $2 \times 2$ block plus three $1 \times 1$ blocks, utilizing only 4 variables instead of 15. This practical simplification explains why the structure is so important that it warrants a fundamental name.

## [Analogy for Understanding](@id analogy-for-understanding)

**Library Analogy:**

Thinking of a G-module is like organizing a vast library (the vector space $V$) where a group of clerks ($G$) has specific rules for handling and moving books (the linear transformations $\rho$). If the clerks' actions satisfy the G-module conditions (linearity, identity, associativity), the chief librarian (Maschke's Theorem) guarantees that the entire library can be physically partitioned into a series of smaller, independent special collections (irreducible submodules). Any comprehensive catalog system ($Q$) that respects the rules of the clerks must then naturally align with these smaller physical partitions, allowing us to manage the whole catalog simply by managing the much smaller catalogs of each special collection.

## [Implementation in NCTSSoS](@id implementation-in-nctssos)

### Wedderburn Decomposition

The implementation uses Wedderburn decomposition to exploit symmetry:

1. **Group representation**: Represent the symmetry group algebraically
2. **Decomposition**: Apply Wedderburn decomposition to the group algebra
3. **Block structure**: Exploit the resulting block diagonal structure

### Integration with Solvers

Symmetry information is integrated with existing solver configurations:

```julia
# Configure solver with symmetry exploitation
solver_config = SolverConfig(
    optimizer=Mosek.Optimizer;
    order=3,
    symmetry=true,           # Enable symmetry exploitation
    simplify_algorithm=:auto  # Choose automatic symmetry detection
)
```

---

## References

1. The sources provided discuss symmetry within the context of **Polynomial Optimization Problems (POPs)**, specifically focusing on methods to exploit symmetry to simplify the computational challenges associated with the Moment-SOS hierarchy.