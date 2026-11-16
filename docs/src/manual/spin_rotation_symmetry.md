# [Spin-Rotation and Spatial Symmetry](@id pauli-charge-singlet-symmetry)

Spin rotations and spatial permutations can reduce a Pauli relaxation without
changing the ground-energy problem. A fixed spin sector is a different choice:
it restricts the states over which the energy is minimized.

| Declaration | Assumption about the Hamiltonian | Assumption about the state |
|:--|:--|:--|
| `SpinRotation(registry, :z)` | Invariant under rotation about Z | No fixed magnetization; averaging over rotations preserves the energy |
| `spin_rotations(registry)` | Invariant under all spin rotations | No fixed total spin; averaging over rotations preserves the energy |
| `SitePermutation(...)` | Invariant under that site permutation | No fixed spin sector; averaging preserves the energy |
| `Sector(S², 0)` | Hermitian objective; any declared symmetry must preserve both it and S² | Supported on the singlet eigenspace S² = 0 |

If the problem has operator constraints, each declared symmetry must preserve
them too. `Sector(S², 0)` is not the operator identity `Zero(S²)`.

## [When to use it](@id pauli-charge-when)

The periodic XXX Heisenberg chain has both spin and spatial invariance:

```math
H_N = \frac14\sum_{i=1}^N\sum_{a\in\{x,y,z\}}
\sigma_i^a\sigma_{i+1\;(\mathrm{mod}\;N)}^a.
```

Use `SpinRotation(registry, :z)` for U(1) invariance, or
`spin_rotations(registry)` for the three SU(2) generators. Neither fixes a
charge or total-spin eigenvalue. Translation and reflection are
`SitePermutation` actions. Declare all generators on `Problem`.

## [API pattern](@id pauli-charge-api-pattern)

```julia
using NCTSSoS, Mosek, LinearAlgebra
import NCTSSoS

BLAS.set_num_threads(1)
solver = MosekBackend(MSK_IPAR_NUM_THREADS=1)

function chain(n)
    registry, (x, y, z) = pauli_variables(n)
    H = 0.25sum(op[i]*op[mod1(i+1,n)] for op in (x,y,z) for i in 1:n)
    S² = 0.25(sum(x)^2 + sum(y)^2 + sum(z)^2)
    spatial = (SitePermutation([2:n; 1]), SitePermutation(reverse(1:n)))
    return registry, H, S², spatial
end

registry, H, S², spatial = chain(4)
problem = Problem(H; symmetry=(spatial..., spin_rotations(registry)...))
relaxation = Relaxation(order=2)
built = NCTSSoS.build(problem, relaxation)
result = solve(built; solver)
@assert isapprox(result.bound, -2; atol=1e-6)
println(result.bound_kind)
println(result.valid_for)
NCTSSoS.block_sides(built)
```

`NCTSSoS.build` constructs the reduced program without solving it. Inspect
`built.sizes` for moment and variable counts, and `NCTSSoS.block_sides(built)` for
final PSD/HPSD sides before Hermitian realification. These counts are not the
native solver-task dimensions; use `NCTSSoS.task_sizes(NCTSSoS.prepare(built; solver))`
for those.

## [What each option does](@id pauli-charge-options)

U(1) rotations separate charge components of the row space. The local
operators Z and S± have charges 0 and ±1. Spatial symmetries reduce the row
space further. The full spin-rotation declaration imposes spin invariance and
can split the remaining cones further. The pipeline handles this projection;
the objective and observables stay in the original Pauli variables.

Spin-invariant first and second moments include zero one-site means, zero
cross-axis two-point means, and equal XX, YY and ZZ two-point means. These
equalities do **not** imply that the state is a singlet. An invariant mixture
of states with nonzero total spin can satisfy them.

## Independent matrix references

For an even ring, use bit strings with N/2 down spins to construct the sparse
Hamiltonian in the Sz = 0 space. An unlike neighboring pair contributes −1/4
on the diagonal and 1/2 when exchanged; a like pair contributes +1/4.
Every integer-spin multiplet contains an Sz = 0 state. Since this Hamiltonian
is SU(2)-invariant and N is even, the minimum in this space is also the
unrestricted minimum. This argument does not restrict the SDP to singlets.

```julia
using SparseArrays, KrylovKit, Random

function ring_reference(n; sign=1)
    @assert iseven(n)
    states = [s for s in 0:(1 << n)-1 if count_ones(s) == n÷2]
    index = Dict(s => i for (i, s) in enumerate(states))
    rows, cols, values = Int[], Int[], Float64[]
    for (column, s) in enumerate(states)
        diagonal = 0.0
        for i in 0:n-1
            j = mod(i+1, n)
            unlike = ((s >> i) & 1) != ((s >> j) & 1)
            diagonal += unlike ? -0.25 : 0.25
            if unlike
                flipped = xor(s, (1 << i) | (1 << j))
                push!(rows, index[flipped]); push!(cols, column); push!(values, 0.5)
            end
        end
        push!(rows, column); push!(cols, column); push!(values, diagonal)
    end
    matrix = sign*sparse(rows, cols, values, length(states), length(states))
    @assert issymmetric(matrix)
    start = randn(MersenneTwister(42), length(states))
    energies, vectors, info = eigsolve(matrix, start, 1, :SR;
                                      ishermitian=true, tol=1e-12)
    @assert info.converged >= 1
    residual = norm(matrix*vectors[1] - energies[1]*vectors[1])
    @assert residual < 1e-10
    return (energy=energies[1], dimension=length(states), residual=residual)
end
```

Build the Pauli matrices independently of the polynomial package. Diagonalize
the full four- and six-site Hamiltonians to check the sparse reference. For
the singlet reference, project onto the zero eigenspace of S² first. These are floating-point exact
diagonalization references, not interval certificates.

```julia
function matrix_references(n)
    paulis = ([0 1; 1 0], [0 -im; im 0], [1 0; 0 -1])
    identity2 = Matrix{ComplexF64}(I, 2, 2)
    ops = [[kron([j == i ? p : identity2 for j in 1:n]...)
            for i in 1:n] for p in paulis]
    Hmat = 0.25sum(op[i]*op[mod1(i+1,n)] for op in ops for i in 1:n)
    S2mat = 0.25sum(sum(op)^2 for op in ops)
    spin = eigen(Hermitian(S2mat))
    singlets = spin.vectors[:, abs.(spin.values) .< 1e-10]
    @assert norm(S2mat*singlets) < 1e-10
    return (antiferromagnet=eigmin(Hermitian(Hmat)),
            ferromagnet=eigmin(Hermitian(-Hmat)),
            singlet_ferromagnet=eigmin(Hermitian(-singlets'*Hmat*singlets)))
end

for n in (4, 6)
    registry, H, S², spatial = chain(n)
    reference = matrix_references(n)
    sparse_reference = ring_reference(n)
    @assert isapprox(sparse_reference.energy, reference.antiferromagnet; atol=1e-10)
    for (label, symmetry) in (
        ("Dense", nothing),
        ("U(1) + spatial", (spatial..., SpinRotation(registry, :z))),
        ("SU(2) + spatial", (spatial..., spin_rotations(registry)...)))
        result = solve(Problem(H; symmetry), Relaxation(order=2); solver)
        @assert isapprox(result.bound, reference.antiferromagnet; atol=1e-6)
        println((n=n, label=label, bound=result.bound,
                 reference=reference.antiferromagnet))
    end
end
```

All order-2 bounds below apply to unrestricted normalized states and agree
with the sparse matrix reference within `1e-6`.

| N | Matrix minimum | Dense | U(1) + spatial | SU(2) + spatial |
|--:|--:|--:|--:|--:|
| 4 | −2.0000000000 | −2.0000000003 | −2.0000000007 | −2.0000000038 |
| 6 | −2.8027756377 | −2.8027756541 | −2.8027756397 | −2.8027756408 |

## Singlet restriction changes the problem

Negate H to obtain a ferromagnetic chain. The unrestricted ground space is
not a singlet. Restricting to `Sector(S², 0)` raises the minimum:

```julia
for n in (4, 6)
    registry, H, S², spatial = chain(n)
    reference = matrix_references(n)
    @assert isapprox(ring_reference(n; sign=-1).energy, reference.ferromagnet; atol=1e-10)
    problem = Problem(-H; symmetry=(spatial..., spin_rotations(registry)...))
    unrestricted = solve(problem, Relaxation(order=2); solver)
    singlet = solve(problem, Relaxation(order=2, strengthen=Sector(S², 0)); solver)
    @assert isapprox(unrestricted.bound, reference.ferromagnet; atol=1e-6)
    @assert isapprox(singlet.bound, reference.singlet_ferromagnet; atol=1e-6)
    @assert singlet.bound > reference.ferromagnet + 0.5
    println((n=n, unrestricted=unrestricted.bound, singlet=singlet.bound))
    println(singlet.bound_kind)
    println(singlet.valid_for)
end
```

| N | Unrestricted matrix minimum | Unrestricted bound | Singlet matrix minimum | Singlet bound |
|--:|--:|--:|--:|--:|
| 4 | −1.0000000000 | −1.0000000047 | 0.0000000000 | −0.0000000001 |
| 6 | −1.5000000000 | −1.5000000043 | −0.8027756377 | −0.8027756409 |

These singlet bounds are **not** lower bounds of the unrestricted minimum.
Read `valid_for`: it states support on the eigenspace S² = 0.

## [Scope and guardrails](@id pauli-charge-scope)

The examples use full order-2 Pauli rows and no grouping or term split.
Declare only symmetries that preserve the Hamiltonian and its constraints.
When adding a sector, its Hermitian operator must use the same registry and
be invariant under all declared actions and derivations. A failed invariance
check names the generator and polynomial; fix the declaration before solving.
There is no public option to skip this check.

The continuous projection may keep a large block when its size or coefficient
memory cap prevents a further split. Inspect the closure report rather than
assuming that every requested symmetry produces smaller blocks.

## [N=16 XXX Heisenberg validation](@id pauli-charge-n16-evidence)

This is the same periodic Hamiltonian at order 2, with 1129 starting rows.
The following function separates build, backend preparation and solve time.
Mosek's optimizer time excludes package work and includes its own presolve.

```julia
function validate_chain16()
    reference = ring_reference(16)
    @assert reference.dimension == binomial(16, 8)
    registry, H, _, spatial = chain(16)
    problem = Problem(H; symmetry=(spatial..., spin_rotations(registry)...))
    build_seconds = @elapsed built = NCTSSoS.build(problem, Relaxation(order=2))
    prepare_seconds = @elapsed prepared = NCTSSoS.prepare(built; solver)
    solve_seconds = @elapsed result = solve(built; solver)
    @assert !isnothing(result.bound)
    @assert isapprox(result.bound, -7.15349290; atol=1e-6)
    @assert result.bound <= reference.energy + 1e-6
    return (sizes=built.sizes, sides=NCTSSoS.block_sides(built),
            task=NCTSSoS.task_sizes(prepared), bound=result.bound,
            reference=reference, valid_for=result.valid_for,
            build_seconds=build_seconds, prepare_seconds=prepare_seconds,
            solve_seconds=solve_seconds,
            optimizer_seconds=result.solution.diagnostics.solve_seconds)
end

validation = validate_chain16()
```

The dense column retains the earlier resource-limit run. The reduced column
was checked again at [6ec1e38](https://github.com/QuantumSOS/NCTSSoS.jl/commit/6ec1e38).
Its row, moment and cone counts are unchanged.

| Quantity | Dense | SU(2) + spatial |
|:--|--:|--:|
| Starting rows | 1129 | 1129 |
| Moments | 163669 | 235 |
| Free real moment coordinates | 163668 | 234 |
| PSD/HPSD blocks | 1 | 31 |
| Largest PSD/HPSD side | 1129 | 9 |
| Native affine rows | 163668 | 234 |
| Largest native real PSD side | 2258 | 18 |
| Build time | 0.421 s | 2.396 s |
| Preparation time | 1.980 s | 0.638 s |
| Solve wall time (prepared task cached) | Not completed | 1.011 s |
| Mosek optimizer time | Not completed | 1.011 s |
| Bound | Not available | −7.15349290243866 |
| Independent Sz = 0 minimum | −7.1422963606 | −7.1422963606 |

The sparse reference has dimension 12870 and eigenpair residual below `1e-10`.
The reduced bound applies to all
normalized states, not just singlets. It lies about 0.01120 below the independent
minimum; order 2 is not exact here. Both numbers are numerical, not interval
certificates.

The earlier dense run had a 600-second process limit. The system's memory
limit stopped it after about 42 seconds, with no bound (signal 9, exit 137).
The kernel recorded about 29.4 GiB resident memory at the kill; peak memory
was not recorded. This is not a solver infeasibility result. The reduced case
and the four- and six-site checks completed together in 32.02 seconds including
startup and compilation. The runs used Julia 1.12.7 on an Intel Xeon at
2.60 GHz, with one Julia, BLAS and Mosek thread.

## [Reading the timings](@id pauli-charge-bottleneck)

These are single validation runs, not a speed benchmark. Build and preparation
times can include compilation. They do not establish a speedup. Native task
dimensions include Hermitian realification; the reported row and cone sides
above do not.

## See also

- [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for the reduction pipeline.
- [Strengthening](@ref strengthening) for sector validity and other conditions.
- [Pauli Symmetry Reduction](@ref pauli-clifford-symmetry) for finite Clifford actions.
