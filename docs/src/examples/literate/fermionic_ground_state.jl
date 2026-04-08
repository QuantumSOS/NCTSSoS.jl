# # [Fermionic Ground State (XY Model)](@id fermionic-ground-state)
#
# This example computes **certified lower bounds** on the ground-state energy
# of a quantum spin chain — and it does so by translating the problem into the
# language of **fermions**.  If those words don't mean anything to you yet,
# good — we'll build every concept from scratch.
#
# ## What is a ground-state energy?
#
# Every quantum system has a ladder of allowed energy levels.  The **ground
# state** is the lowest rung — the state the system settles into when all
# thermal energy is removed (absolute zero).  Its energy $E_0$ is the smallest
# the system can ever have.
#
# A natural guess is that this minimum energy should be zero — a perfectly
# still, inert system.  But quantum mechanics forbids it: the Heisenberg
# uncertainty principle requires every bound particle to retain some residual
# motion ("zero-point energy"), so $E_0$ is generically nonzero.  For the
# spin chains in this example, $E_0$ turns out to be *negative*, reflecting
# how quantum correlations between neighbouring spins lower the total
# energy.  The precise value of $E_0$ encodes essential physics: phase
# transitions, stability, and the nature of those correlations.
#
# ## Bosons vs. fermions — the two families of particles
#
# Nature divides all particles into two families based on a single rule about
# **swapping identical particles**:
#
# - **Bosons** (photons, Higgs): the quantum state is *unchanged* under a
#   swap.  Bosons happily pile into the same quantum state — this is why
#   lasers work.
# - **Fermions** (electrons, quarks, neutrons): the quantum state picks up a
#   **minus sign** under a swap.  This sign has a stunning consequence called
#   the **Pauli exclusion principle**: no two fermions can occupy the same
#   state.  (If two were in the same state, swapping them should change the
#   sign, but also change nothing — both conditions require the state to be
#   zero, i.e. it cannot exist.)
#
# The Pauli exclusion principle is why atoms have shells, why matter doesn't
# collapse, and why white-dwarf stars resist gravity.
#
# ## Creation and annihilation operators
#
# Tracking many-particle quantum states by writing wavefunctions
# $\psi(x_1, x_2, \ldots)$ becomes intractable fast.  A smarter language
# asks: *how many particles are in each available state?*
#
# | Operator | What it does |
# |----------|-------------|
# | $a_i^\dagger$ (**creation**) | adds one particle to site $i$ |
# | $a_i$ (**annihilation**) | removes one particle from site $i$ |
# | $n_i = a_i^\dagger a_i$ (**number**) | counts particles at site $i$ (0 or 1 for fermions) |
#
# Think of sites as hotel rooms: $a_i^\dagger$ checks a guest into room $i$,
# $a_i$ checks them out, and $n_i$ tells the front desk whether the room is
# occupied.  For fermions, each room holds **at most one guest**.
#
# ### Anticommutation relations (CAR)
#
# The exclusion principle and the swap-sign rule are encoded in three algebraic
# identities called the **canonical anticommutation relations** (CAR), where
# $\{A,B\} \equiv AB + BA$:
#
# ```math
# \{a_i,\, a_j^\dagger\} = \delta_{ij}, \qquad
# \{a_i,\, a_j\} = 0, \qquad
# \{a_i^\dagger,\, a_j^\dagger\} = 0.
# ```
#
# The first says that at the same site ($i = j$), the two orderings of
# create-and-annihilate always add up to the identity:
# $a_i a_i^\dagger + a_i^\dagger a_i = 1$.  The second and third say
# "swapping the order of two same-type operators flips the sign."
# A direct corollary is $(a_i)^2 = 0$: removing a fermion from the same
# site twice always gives zero — the operator itself is nilpotent.
#
# In `NCTSSoS.jl`, [`create_fermionic_variables`](@ref) gives you operators
# that already obey these rules — every product is automatically simplified
# (normal-ordered) by the library.
#
# ## The XY spin chain
#
# Now for a concrete physical system.  Picture $N$ tiny magnets (spin-½
# particles) arranged on a ring.  Each spin can point "up" ($\uparrow$) or
# "down" ($\downarrow$), and neighbouring spins interact.  The **periodic XY
# model** describes the interaction:
#
# ```math
# H_{\mathrm{spin}} = \frac{j_c}{4}\sum_{i=1}^{N}
#     \bigl(\sigma^x_i \sigma^x_{i+1} + \sigma^y_i \sigma^y_{i+1}\bigr),
#     \qquad \sigma_{N+1} \equiv \sigma_1,
# ```
#
# where $\sigma^x, \sigma^y$ are Pauli matrices — $2\times 2$ matrices that
# encode the spin's quantum state along different axes.  The products
# $\sigma^x_i \sigma^x_{i+1}$ and $\sigma^y_i \sigma^y_{i+1}$ describe an
# **exchange interaction**: a spin excitation at site $i$ can hop to site
# $i{+}1$, and vice versa — like passing a ball around a ring.  The
# ground-state energy tells us the minimum energy stored in this ring of
# interacting magnets at zero temperature.
#
# The XY model is the canonical test case of many-body quantum physics:
# simple enough to solve exactly, yet rich enough to exhibit quantum phase
# transitions and long-range entanglement.  If a new numerical method can't
# reproduce it, something is wrong.
#
# ## Jordan–Wigner transformation — from spins to fermions
#
# Spins on different sites commute (they don't care about each other's
# ordering), but fermions on different sites *anticommute* (swapping costs a
# sign).  The **Jordan–Wigner transformation** (1928) bridges this gap by
# attaching a "parity string" to each site:
#
# ```math
# a_j = \Bigl(\prod_{k=1}^{j-1} \sigma^z_k\Bigr)\,\sigma^-_j,
# ```
#
# where $\sigma^-_j = (\sigma^x_j - i\,\sigma^y_j)/2$ is the spin-lowering
# operator (it flips a spin from $\uparrow$ to $\downarrow$), and each
# $\sigma^z_k$ equals $+1$ or $-1$ depending on the spin state at site
# $k$.  The product $\prod_{k<j}\sigma^z_k$ — the **Jordan–Wigner
# string** — accumulates a minus sign for every occupied site to the left
# of $j$, converting the local (commuting) spin algebra into the nonlocal
# (anticommuting) fermion algebra.
#
# **Analogy:** Imagine fermions are people walking single-file in a hallway.
# If person $j$ wants to move, they must squeeze past everyone to their left;
# each pass costs a sign.  The Jordan–Wigner string is the mathematical
# ledger that keeps track of all those sign flips.
#
# ### Why nearest-neighbour terms simplify
#
# For adjacent sites $i$ and $i{+}1$, the JW strings almost cancel.
# Site $i{+}1$'s string contains every factor of site $i$'s string plus
# one extra $\sigma^z_i$.  In the product $a_i^\dagger a_{i+1}$, the
# common part squares to identity, leaving
# $\sigma^+_i \sigma^z_i \sigma^-_{i+1}$.  A quick matrix check gives
# $\sigma^+ \sigma^z = -\sigma^+$, so:
#
# ```math
# a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i
#   = -\bigl(\sigma^+_i \sigma^-_{i+1} + \sigma^+_{i+1} \sigma^-_i\bigr).
# ```
#
# Since $\sigma^x_i \sigma^x_{i+1} + \sigma^y_i \sigma^y_{i+1} =
# 2(\sigma^+_i \sigma^-_{i+1} + \sigma^-_i \sigma^+_{i+1})$, each
# interior spin term maps to a simple fermion hop — **no leftover JW
# strings**.  The XY Hamiltonian, which looked like an interacting spin
# problem, becomes a **free-fermion** Hamiltonian — no interactions
# between fermions at all!  That makes it exactly solvable.
#
# ## The parity operator — why the boundary term is special
#
# The interior hopping terms translate cleanly, but the **boundary** term
# (site $N$ ↔ site $1$, closing the ring) picks up the full-chain parity
# string.  This string is exactly the **fermion parity operator**:
#
# ```math
# P = \prod_{i=1}^{N}(1 - 2\,a_i^\dagger a_i).
# ```
#
# What does $P$ measure?  Each factor $(1-2n_i)$ equals $+1$ if site $i$ is
# empty and $-1$ if occupied.  The product over all sites gives:
# - $P = +1$ when the total number of fermions is **even**,
# - $P = -1$ when it is **odd**.
#
# $P$ commutes with the Hamiltonian and satisfies $P^2 = I$, so the Hilbert
# space splits into two independent **parity sectors**: even ($P=+1$) and
# odd ($P=-1$).  The full fermionic Hamiltonian is
#
# ```math
# H = -\frac{j_c}{2}\sum_{i=1}^{N-1}
#     \bigl(a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i\bigr)
#   + \frac{j_c}{2}\,P\,
#     \bigl(a_N^\dagger a_1 + a_1^\dagger a_N\bigr).
# ```
#
# ### Restricting to the even sector
#
# In the even-parity sector $P = +1$, the boundary factor $P$ becomes just
# $+1$, so the Hamiltonian simplifies to a purely quadratic (free) form:
#
# ```math
# H_{\mathrm{even}} = -\frac{j_c}{2}\sum_{i=1}^{N-1}
#     \bigl(a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i\bigr)
#   + \frac{j_c}{2}
#     \bigl(a_N^\dagger a_1 + a_1^\dagger a_N\bigr).
# ```
#
# This reduced Hamiltonian already *is* the even-sector problem.  For a very
# small system we can additionally certify $P = I$ inside the SDP, but that
# parity polynomial has degree $2N$, so it stops fitting inside low-order
# relaxations once $N$ grows.
#
# !!! warning "Do not impose `P - I = 0` blindly"
#     The polynomial `P - I` has degree `2N`.  For `N = 4`, an order-4
#     relaxation can still represent it, so the explicit equality constraint
#     works.  For `N = 6` and `N = 8`, the order-3 sparse relaxations used
#     below only carry moments up to degree 6, so a model with
#     `eq_constraints = [P - I]` is under-supported.
#
#     `NCTSSoS.jl` now rejects that setup instead of silently dropping the
#     missing high-degree moments.  So for the larger cases below we solve the
#     already-reduced Hamiltonian $H_{\mathrm{even}}$ directly rather than
#     adding a separate parity polynomial.

# !!! tip "What about the odd-parity sector?"
#     In the odd sector ($P = -1$), the boundary sign flips: all hopping
#     terms get the *same* sign, giving periodic (not anti-periodic)
#     boundary conditions.  For the XY model with $j_c > 0$, the ground
#     state lives in the even sector for all chain lengths we consider —
#     you can verify this by checking that the number of filled Fermi-sea
#     modes is even.  In general, one should solve both sectors and take
#     the lower energy.

# ---
#
# With the physics in place, we can now compute.  We'll first obtain the
# **exact** answer analytically (to have something to check against), then
# solve the same problem numerically with `NCTSSoS.jl`.
#
# ## Exact energy from the free-fermion spectrum
#
# Because $H_{\mathrm{even}}$ is quadratic in the fermion operators, the
# problem reduces to diagonalising an $N \times N$ **single-particle
# hopping matrix** $h$:
#
# ```math
# h = -\frac{j_c}{2}\begin{pmatrix}
#   0 & 1 & 0 & \cdots & 0 & -1 \\
#   1 & 0 & 1 & \cdots & 0 &  0 \\
#   0 & 1 & 0 & \cdots & 0 &  0 \\
#   \vdots & & & \ddots & & \vdots \\
#   0 & 0 & 0 & \cdots & 0 &  1 \\
#  -1 & 0 & 0 & \cdots & 1 &  0
# \end{pmatrix}.
# ```
#
# The $-1$ entries in the corners (instead of $+1$) are the
# **anti-periodic boundary conditions** from the parity twist.  Each
# eigenvalue $\varepsilon_k$ of $h$ represents one available energy level.
# The ground state is obtained by filling every level with
# $\varepsilon_k < 0$ (like stacking blocks into the lowest shelves of a
# bookcase).  This filled set of levels is called the **Fermi sea**, and
# the ground-state energy is the sum of all negative eigenvalues.
#
# ### Method 1 — numerical diagonalisation

using LinearAlgebra

"""
    xy_exact_energy(N; j_c=1.0)

Exact ground-state energy of the periodic N-site XY spin chain via
numerical diagonalisation of the hopping matrix (even-parity sector).
"""
function xy_exact_energy(N::Int; j_c::Real = 1.0)
    ## Single-particle hopping matrix (anti-periodic BC from parity twist)
    h = zeros(N, N)
    for i in 1:N-1
        h[i, i+1] = h[i+1, i] = -j_c / 2
    end
    h[N, 1] = h[1, N] = j_c / 2          # anti-periodic boundary
    εs = eigvals(Symmetric(h))
    return sum(ε for ε in εs if ε < -1e-14)
end

# ### Method 2 — analytic Fourier formula
#
# The anti-periodic hopping matrix is diagonalised by a discrete Fourier
# transform.  The **half-integer shift** in the momenta comes from the
# anti-periodic boundary:
#
# ```math
# k_n = \frac{2\pi}{N}\Bigl(n + \tfrac{1}{2}\Bigr), \qquad
# \varepsilon_n = -j_c \cos(k_n), \qquad n = 0, 1, \ldots, N{-}1.
# ```
#
# **Worked example for $N = 4$.**  The momenta are
# $k_n = \frac{\pi}{4}(2n+1)$ for $n = 0,1,2,3$:
#
# | $n$ | $k_n$     | $\varepsilon_n = -\cos(k_n)$ |
# |:---:|:---------:|:----------------------------:|
# | 0   | $\pi/4$   | $-1/\sqrt{2} \approx -0.707$ |
# | 1   | $3\pi/4$  | $+1/\sqrt{2} \approx +0.707$ |
# | 2   | $5\pi/4$  | $+1/\sqrt{2} \approx +0.707$ |
# | 3   | $7\pi/4$  | $-1/\sqrt{2} \approx -0.707$ |
#
# Two eigenvalues are negative ($n = 0$ and $n = 3$).  Filling those two
# modes — an **even** number, consistent with the even-parity sector —
# gives $E_0 = -2/\sqrt{2} = -\sqrt{2} \approx -1.414$.

"""
    xy_exact_energy_analytic(N; j_c=1.0)

Same result as `xy_exact_energy`, but via the closed-form Fourier
eigenvalues — no matrix diagonalisation needed.
"""
function xy_exact_energy_analytic(N::Int; j_c::Real = 1.0)
    E = 0.0
    for n in 0:N-1
        kn = 2π / N * (n + 0.5)
        εn = -j_c * cos(kn)
        if εn < -1e-14
            E += εn
        end
    end
    return E
end

# Both methods agree:

for N in [4, 6, 8]
    E_num = xy_exact_energy(N)
    E_ana = xy_exact_energy_analytic(N)
    @assert abs(E_num - E_ana) < 1e-12   #src
    println("N = $N :  E₀ = $E_num  (analytic: $E_ana)")
end

# ## Solving with NCTSSoS.jl
#
# Now we solve the *same* problem using **semidefinite programming (SDP)**
# relaxations, without knowing the exact answer.
#
# The core idea: diagonalising the full quantum Hamiltonian scales
# exponentially in the number of sites, but we can ask a different,
# tractable question.  An SDP solver searches for the *lowest* value of
# $\langle H \rangle$ that is consistent with the CAR, the chosen
# even-sector Hamiltonian, and the requirement that expectation values come
# from a valid quantum state.  Because every valid quantum state satisfies
# those conditions, the SDP answer is guaranteed to be **at or below** the
# true ground-state energy — it is a **certified lower bound**.

# !!! note "SDP Solver"
#     These examples use [Mosek](https://www.mosek.com/) via `MosekTools`.
#     Any SDP-capable solver works: replace `Mosek.Optimizer` with
#     `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source alternatives.

using NCTSSoS, MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0,                               # quiet output
    "MSK_IPAR_NUM_THREADS" => 0)                        # use all threads

# -------------------------------------------------------------------
# ## Case 1 — N = 4 (tight at first iteration)
#
# For $N = 4$ the parity polynomial has degree $2N = 8$.  At **order 4**
# the moment matrix includes monomials up to degree 8, which covers the
# full fermionic Fock space, and the parity constraint is enforced.
# A single call to [`cs_nctssos`](@ref) is sufficient.
# -------------------------------------------------------------------

# ### Step 1 — Create fermionic variables
#
# [`create_fermionic_variables`](@ref) returns a registry and two vectors:
# annihilation operators `a[i]` and creation operators `a_dag[i]`.
# All CAR relations are encoded automatically — you never need to manually
# specify $\{a_i, a_j^\dagger\} = \delta_{ij}$.

N₁ = 4
j_c = 1.0
registry₁, (a₁, a₁_dag) = create_fermionic_variables(1:N₁)

# ### Step 2 — Build the even-parity Hamiltonian
#
# We write $H_{\mathrm{even}}$ directly in terms of the fermionic operators.
# Each term $a_i^\dagger a_{i+1}$ describes a fermion hopping from site
# $i{+}1$ to site $i$.  Notice the **sign flip** on the boundary term
# ($+j_c/2$ instead of $-j_c/2$): this is the anti-periodic boundary
# condition imposed by the even-parity sector.

ham₁ = sum(
    -ComplexF64(j_c / 2) * (a₁_dag[i] * a₁[i+1] + a₁_dag[i+1] * a₁[i])
    for i in 1:N₁-1
) + ComplexF64(j_c / 2) * (a₁_dag[N₁] * a₁[1] + a₁_dag[1] * a₁[N₁])

# ### Step 3 — Parity constraint: $P - I = 0$
#
# We build $P = \prod_{i=1}^{N}(1 - 2\,a_i^\dagger a_i)$ and impose
# $P = I$ as an equality constraint.  Each factor $(1-2n_i)$ maps site $i$
# to $+1$ (empty) or $-1$ (occupied); the product is $+1$ exactly when the
# total particle number is even.

parity₁ = prod(one(ham₁) - 2.0 * a₁_dag[i] * a₁[i] for i in 1:N₁)

# ### Step 4 — Solve the SDP relaxation
#
# [`polyopt`](@ref) packages the Hamiltonian and constraints into a
# polynomial optimisation problem.  [`cs_nctssos`](@ref) builds and solves
# the SDP relaxation, exploiting sparsity.  The `order` parameter controls
# the size of the moment matrix (higher = tighter but costlier).

pop₁ = polyopt(ham₁, registry₁; eq_constraints = [parity₁ - one(ham₁)])
config₁ = SolverConfig(optimizer = SOLVER, order = 4, ts_algo = MMD())
result₁ = cs_nctssos(pop₁, config₁)

# ### Step 5 — Verify against the exact answer

exact₁ = xy_exact_energy(N₁)
@assert abs(result₁.objective - exact₁) < 1e-5   #src
println("N = $N₁:  SDP = $(result₁.objective),  exact = $exact₁,  gap = $(abs(result₁.objective - exact₁))")

# The SDP bound matches the exact energy to many digits — the relaxation
# is **tight** for this small system.

# -------------------------------------------------------------------
# ## Scaling up: N = 6 and N = 8
#
# For larger systems, one SDP iteration may produce a loose bound.  Why?
# The solver exploits **term sparsity** to break the (exponentially large)
# moment matrix into smaller blocks.  Some cross-correlations between
# distant sites can be missed in the first decomposition.
#
# There is also a modeling detail that matters here: for `N = 6` and `N = 8`,
# the explicit parity polynomial has degree 12 and 16, so the order-3 sparse
# relaxation below cannot support it as an `eq_constraint`.  That was exactly
# what broke the old page.  We therefore keep the same reduced Hamiltonian
# $H_{\mathrm{even}}$, but we do **not** add `P - I = 0` again for these
# larger cases.
#
# A second call to [`cs_nctssos_higher`](@ref) reads the moment values from
# the first solve, discovers which cross-correlations matter, and refines
# the sparsity pattern.  The result is typically tight.
# -------------------------------------------------------------------

# ### Helper — build and solve for arbitrary N

"""
    solve_xy(N; j_c=1.0, order=3)

Build the already-reduced even-sector fermionic XY Hamiltonian for `N` sites,
run one round of `cs_nctssos` and one of `cs_nctssos_higher`,
and return `(first_obj, refined_obj, exact)`.

For `N ≥ 6` we intentionally do **not** add `P - I = 0` as an explicit
`eq_constraint`: that degree-`2N` polynomial does not fit inside the low-order
sparse relaxation used here.
"""
function solve_xy(N::Int; j_c::Real = 1.0, order::Int = 3)
    registry, (a, a_dag) = create_fermionic_variables(1:N)

    ham = sum(
        -ComplexF64(j_c / 2) * (a_dag[i] * a[i+1] + a_dag[i+1] * a[i])
        for i in 1:N-1
    ) + ComplexF64(j_c / 2) * (a_dag[N] * a[1] + a_dag[1] * a[N])

    pop = polyopt(ham, registry)
    config = SolverConfig(optimizer = SOLVER, order = order, ts_algo = MMD())

    res1 = cs_nctssos(pop, config)
    res2 = cs_nctssos_higher(pop, res1, config)

    return (first = res1.objective, refined = res2.objective,
            exact = xy_exact_energy(N; j_c))
end

# ### N = 6

r6 = solve_xy(6; order = 3)
println("N = 6:  first = $(r6.first),  refined = $(r6.refined),  exact = $(r6.exact)")
println("  gap after refinement: $(abs(r6.refined - r6.exact))")
@assert abs(r6.refined - r6.exact) < 1e-5   #src

# ### N = 8

r8 = solve_xy(8; order = 3)
println("N = 8:  first = $(r8.first),  refined = $(r8.refined),  exact = $(r8.exact)")
println("  gap after refinement: $(abs(r8.refined - r8.exact))")
@assert abs(r8.refined - r8.exact) < 1e-5   #src

# -------------------------------------------------------------------
# ## Results
#
# | ``N`` | Exact ``E_0``             | First iteration | After `cs_nctssos_higher` |
# |:-----:|:-------------------------:|:---------------:|:-------------------------:|
# | 4     | ``-\sqrt{2} \approx -1.414``  | tight       | —                         |
# | 6     | ``-\sqrt{3} \approx -1.732``  | loose       | tight                     |
# | 8     | ``\approx -2.613``            | loose       | tight                     |
#
# For $N = 4$ the first SDP iteration already matches the exact energy.
# For $N = 6$ and $N = 8$, the sparsity-refined second iteration
# ([`cs_nctssos_higher`](@ref)) tightens the bound to the exact value.
# -------------------------------------------------------------------

# ## Summary
#
# ### The physics in one paragraph
#
# We started with $N$ interacting spins on a ring (the XY model).  The
# Jordan–Wigner transformation rewrote the problem in terms of fermions —
# particles that obey anticommutation rules enforced automatically by
# `NCTSSoS.jl`.  The boundary term of the ring acquired a parity operator
# $P$ that measures whether the total fermion count is even or odd.
# Restricting to the even sector fixes the boundary sign and turns the
# Hamiltonian into a free-fermion problem whose ground-state energy we can
# check exactly.  For `N = 4` we also enforce `P = I` explicitly; for
# `N = 6` and `N = 8` we keep the already-reduced Hamiltonian and avoid the
# high-degree parity polynomial that would break the low-order sparse
# relaxation.  In every case the SDP matches the exact answer, certifying a
# rigorous lower bound.
#
# ### The code recipe
#
# | Step | What we did | API |
# |------|-------------|-----|
# | 1 | Create fermionic operators (CAR built in) | [`create_fermionic_variables`](@ref) |
# | 2 | Write the even-sector Hamiltonian | standard Julia arithmetic on operators |
# | 3 | Keep the large-`N` model low-degree by solving $H_{\mathrm{even}}$ directly (the `N = 4` sanity check also enforces `P = I`) | [`polyopt`](@ref) |
# | 4 | Solve & refine the SDP relaxation | [`cs_nctssos`](@ref), [`cs_nctssos_higher`](@ref) |
# | 5 | Compare with exact free-fermion energy | `xy_exact_energy` (above) |
#
# ### Spin–fermion dictionary
#
# The Jordan–Wigner transformation sets up a one-to-one correspondence
# between spin and fermion languages.  This table collects every identity
# used above:
#
# | Identity | Spin form | Fermion form |
# |----------|-----------|-------------|
# | Raising / lowering | $\sigma^\pm = (\sigma^x \pm i\sigma^y)/2$ | $a_j = (\prod_{k<j}\sigma^z_k)\,\sigma^-_j$ |
# | Nilpotency | $(\sigma^\pm)^2 = 0$ | $(a_j)^2 = 0,\; (a_j^\dagger)^2 = 0$ |
# | Number operator | $n_j = \sigma^+_j\sigma^-_j = (1+\sigma^z_j)/2$ | $n_j = a_j^\dagger a_j$ |
# | CAR (same site) | $\sigma^+\sigma^- + \sigma^-\sigma^+ = 1$ | $\{a_j, a_j^\dagger\} = 1$ |
# | Parity per site | $1 - 2\sigma^+_j\sigma^-_j = -\sigma^z_j$ | $1 - 2n_j$ |
# | Full parity ($N$ even) | $P = \prod_j \sigma^z_j$ | $P = \prod_j(1 - 2n_j)$ |
# | Nearest-neighbour hop | $\sigma^+_i\sigma^-_{i+1} + \text{h.c.}$ | $-(a_i^\dagger a_{i+1} + \text{h.c.})$ |
#
# ### Key physical concepts
#
# | Concept | One-liner |
# |---------|-----------|
# | **Fermion** | Particle that obeys the Pauli exclusion principle (at most one per state) |
# | **Creation / annihilation** | Operators $a_i^\dagger, a_i$ that add / remove a fermion at site $i$ |
# | **CAR** | $\{a_i, a_j^\dagger\}=\delta_{ij}$ — the algebraic encoding of exclusion + swap sign |
# | **Jordan–Wigner** | Exact mapping from spin operators to fermion operators via a parity string |
# | **Parity operator** $P$ | $\prod_i(1-2n_i)$: $+1$ for even total particle count, $-1$ for odd |
# | **Ground-state energy** | Lowest eigenvalue of $H$ — the energy at absolute zero |
#
# ### When to use the fermionic interface
#
# - Models that are naturally quadratic (or low-degree) in fermion operators
# - Jordan–Wigner reductions of spin chains
# - Problems where parity or particle-number sectors must be fixed
#
# For Pauli-level formulations of ground-state problems, see
# [Ground State Energy](@ref ground-state-energy).
# For a showcase of fermionic algebra rules (CAR, nilpotency, parity), see
# [PBW Algebra Showcase](@ref pbw-algebras-showcase).
#
# ### References
#
# - P. Jordan and E. Wigner, "Über das Paulische Äquivalenzverbot,"
#   *Zeitschrift für Physik* **47**, 631 (1928).
# - E. Lieb, T. Schultz, D. Mattis, "Two soluble models of an
#   antiferromagnetic chain," *Annals of Physics* **16**, 407 (1961) —
#   the classic paper solving the XY model via Jordan–Wigner.
# - M. Nielsen, "The fermionic canonical commutation relations and the
#   Jordan–Wigner transform," notes at
#   [futureofmatter.com](https://futureofmatter.com/assets/fermions_and_jordan_wigner.pdf).
# - S. Sachdev, *Quantum Phase Transitions*, 2nd ed., Cambridge University
#   Press (2011) — textbook treatment of parity sectors and spin-chain
#   phase transitions.
