# Polynomial Optimization Problems Solved in `test`

## Oracle Validation Task

### Objective
Validate NCTSSoS.jl against NCTSSOS (original package) as ground truth oracle:
1. **Optimal value match** (within `rtol=1e-6`)
2. **Moment block side lengths identical** (`sides`)
3. **Unique elements count match** (`nuniq = length(data.ksupp)`)
4. **MosekTools** used as backend for both packages (solve locally)

### Prerequisites
- **NCTSSOS** (a800): `/home/yushengzhao/NCTSSOS`
- **NCTSSOS** (local): `/Users/yushengzhao/projects/NCTSSOS`
- **NCTSSoS.jl** (local): `/Users/yushengzhao/projects/NCTSSoS-review`
- **MosekTools** licensed on both machines

### Verification Criteria

| Metric | NCTSSOS Field | NCTSSoS Equivalent | Match Required |
|--------|---------------|-------------------|----------------|
| Optimal value | `opt` | `result.objective` | ✓ (rtol=1e-6) |
| Block side lengths | see below | `block_sizes(result)` | ✓ exact |
| Unique elements | `length(data.ksupp)` | `nuniq(result)` | ✓ exact |

**Block side lengths extraction:**
```julia
# nctssos_first / nctssos_higher! (dense/TS)
opt, data = nctssos_first(...)
sides = [size(M, 1) for M in data.moment]

# cs_nctssos_first / cs_nctssos_higher! (CS, flattened over cliques)
opt, data = cs_nctssos_first(...)
sides = [size(M, 1) for clique in data.moment for M in clique]
```

**Unique elements (affine constraints in SDP):**
```julia
nuniq = length(data.ksupp)
```

### Sign Convention (Critical)
| Package | Default | To maximize f |
|---------|---------|---------------|
| NCTSSOS | maximizes | solve directly |
| NCTSSoS | minimizes | solve `min(-f)`, negate result |

**Relation**: `NCTSSOS_max(f) = -NCTSSoS_min(-f)`

### Algebra Type Mapping
| NCTSSoS | NCTSSOS | Notes |
|---------|---------|-------|
| `NonCommutativeAlgebra` | (no constraint) | Default |
| `UnipotentAlgebra` | `constraint="unipotent"` | U² = I |
| `ProjectorAlgebra` | `constraint="projector"` | P² = P |
| `PauliAlgebra` | N/A | NCTSSoS-specific |
| `FermionicAlgebra` | N/A | NCTSSoS-specific |
| `BosonicAlgebra` | N/A | NCTSSoS-specific |

### Sparsity Mapping
| NCTSSoS | NCTSSOS |
|---------|---------|
| `NoElimination` | `TS=false, CS=false` |
| `MaximalFusion` (CS) | `CS="MF"` |
| `MMD` (TS) | `TS="MD"` |
| Combined | `CS="MF", TS="MD"` |

### Workflow

1. **Generate oracles** (on a800 with MosekTools):
   ```bash
   cd /home/yushengzhao/NCTSSOS
   julia --project test/oracles/<category>_oracle.jl  # e.g., moment_oracle.jl
   ```

2. **Update oracle file**: Copy output to corresponding `test/oracles/<category>_oracles.jl`

3. **Run NCTSSoS tests locally** (with MosekTools):
   ```bash
   julia --project -e 'using Pkg; Pkg.test(test_args=["--local"])'
   ```

4. **Validate**:
   ```julia
   using NCTSSoS.Test.Oracles: validate_result
   passed, msg = validate_result("moment.CHSH_Unipotent", result)
   @test passed
   ```

### Oracle Infrastructure
```
test/oracles/
├── scripts/                       # NCTSSOS oracle generation scripts (run on a800)
│   ├── oracle_utils.jl            # Shared helpers for oracle extraction
│   ├── nctssos_chsh.jl            # CHSH with all sparsity variants
│   ├── nctssos_i3322.jl           # I_3322 Bell inequality
│   ├── nctssos_heisenberg_star.jl # Heisenberg star graph
│   ├── nctssos_cs_ts_n10.jl       # Large CS+TS example
│   ├── nctssos_example1.jl        # Unconstrained 3-var NC
│   ├── nctssos_example2.jl        # Constrained 2-var NC
│   ├── nctssos_corr_sparsity.jl   # Correlative sparsity example
│   ├── nctssos_rosenbrock.jl      # Generalized Rosenbrock
│   └── nctssos_state_poly.jl      # State polynomial (7.2.x)
└── results/                       # Generated oracle values
    ├── chsh_oracles.jl            # CHSH ground truth (opt, sides, nuniq)
    ├── i3322_oracles.jl           # I3322 oracles
    ├── heisenberg_star_oracles.jl # Heisenberg star oracles
    ├── cs_ts_n10_oracles.jl       # CS+TS n=10 oracles
    └── ...                        # Other oracle results
```

### Oracle Entry Format
```julia
# Each oracle entry contains:
"TestName" => (
    opt = 2.8284,                    # NCTSSOS optimal value
    sides = [5],                     # block side lengths (see extraction above)
    nuniq = 11,                      # length(data.ksupp) = unique elements
    notes = "Dense, MosekTools",
)
```

---

## Physics Tests (`test/physics/`)

*Note: Pauli/Fermionic/Bosonic algebras have no NCTSSOS oracle; validate against exact diagonalization.*

- [ ] Bell inequalities — all instances via `tester` — `bell_inequalities.jl:88`
- [ ] Bilocal networks — Example 8.1.1 — `bilocal_networks.jl:113`
- [ ] Bose-Hubbard ground state — vacuum reference — `bose_hubbard.jl:204`
- [ ] Bose-Hubbard chain ground state — (broken/unbounded note) — `bose_hubbard.jl:255`
- [ ] Fermionic chain ground state — N=4, anti-periodic BC — `fermionic_chain.jl:48`
- [ ] Fermionic parity superselection — constraint injection in moment relaxation — `fermionic.jl:89`
- [ ] Heisenberg XXX model — N=6 — `heisenberg.jl:23`
- [ ] Heisenberg J1-J2 model — N=4, loop over `J2s` — `heisenberg.jl:45`
- [ ] Heisenberg J1-J2 model — N=6 — `heisenberg.jl:69`
- [ ] Heisenberg 2D model — 3x3 — `heisenberg.jl:94`
- [ ] XY model ground state — N=4, periodic BC — `xy_model.jl:27`
- [ ] PXP model bounds — `cs_nctssos_with_blockade` (lower bound) — `pxp.jl:92`
- [ ] PXP model bounds — `cs_nctssos_with_blockade` (upper bound) — `pxp.jl:96`
- [ ] PXP model bounds — JuMP solve call inside helper (`optimize!`) — `pxp.jl:60`

## Trace Polynomial (`test/solvers/trace_poly.jl`)

| Test | Oracle Key | Algebra | Order | Expected |
|------|------------|---------|-------|----------|
| Example 6.1 (order=2) | `Trace_6_1_d2` | Projector | 2 | -0.0467 |
| Example 6.1 (order=3, `USE_LOCAL`) | `Trace_6_1_d3` | Projector | 3 | -0.0312 |
| Example 6.2.0 (CHSH trace) | `Trace_6_2_0` | Unipotent | 1 | -2.8284 |
| Example 6.2.1 (squared trace, `USE_LOCAL`) | `Trace_6_2_1_d2` | Unipotent | 2 | -4.0 |
| Example 6.2.2 (covariance trace) | `Trace_6_2_2` | Unipotent | 2 | -5.0 |

## State Polynomial (`test/solvers/state_poly.jl`)

| Test | Oracle Key | Algebra | Order | Expected |
|------|------------|---------|-------|----------|
| 7.2.0 | `State_7_2_0` | Unipotent | 1 | -2.8284 |
| 7.2.0 (Sparse) | `State_7_2_0_TS` | Unipotent | 1 | -2.8284 |
| 7.2.1 | `State_7_2_1_d3` | Unipotent | 3 | -4.0 |
| 7.2.2 | `State_7_2_2` | Unipotent | 2 | -5.0 |
| 7.2.3 | `State_7_2_3` | Unipotent | 2 | -3.5115 |
| 7.2.3 (Sparse) | `State_7_2_3_TS` | Unipotent | 2 | -3.5115 |

- [ ] 7.2.0 direct moment (`dualize=false`) — `state_poly.jl:128`
- [ ] 7.2.0 SOS dual (`dualize=true`) — `state_poly.jl:129`
- [ ] 7.2.1 order-3 tight bound — `state_poly.jl:152`
- [ ] 7.2.2 covariance (duplicate section) — `state_poly.jl:169`
- [ ] 7.2.3 direct moment (`dualize=false`) — `state_poly.jl:188`
- [ ] 7.2.3 SOS dual (`dualize=true`) — `state_poly.jl:189`
- [ ] direct moment solve (moment, `dualize=false`) — `state_poly.jl:206`
- [ ] direct moment solve (SOS, `dualize=true`) — `state_poly.jl:207`

## Moment Solver (`test/solvers/moment.jl`)

| Test | Oracle Key | Algebra | Sparsity | Expected |
|------|------------|---------|----------|----------|
| CHSH inequality | `CHSH_Unipotent` | Unipotent | Dense | -2.8284 |
| CS/TS example (n=10) | `CS_TS_n10_d3` | NC | MF+MMD | 3.0113 |
| Heisenberg star graph | `Heisenberg_Star_n10` | Unipotent | MF | -1.0 |
| Example 1 (Dense) | `Example1_Dense` | NC | Dense | 0.0 |
| Example 1 (Sparse) | `Example1_TS_MMD` | NC | MMD | -0.0036 |
| Example 2 (Dense) | `Example2_Dense` | NC | Dense | -1.0 |
| Example 2 (Term Sparse) | `Example2_CS_TS` | NC | MF+MMD | -1.0 |
| Correlative sparsity (CS) | `CorrSparsity_CS` | NC | MF | 0.9975 |
| Correlative sparsity (TS) | `CorrSparsity_TS` | NC | MMD | 0.9975 |

## SOS Dualization (`test/solvers/sos.jl`)

| Test | Oracle Key | Algebra | Notes |
|------|------------|---------|-------|
| I_3322 inequality | `Projector_I3322_Dense_d3` | Projector | order=3 |
| CS/TS example | `CS_TS_n10_d3` | NC | Same as moment |
| trivial example | `NC_Dualization_Trivial` | NC | true_min=3.0 |
| example 1 | `NC_Dualization_Example1` | NC | 3-var |
| example 2 (Dense) | `NC_Dualization_Example2_Dense` | NC | Constrained |
| example 2 (Term Sparse) | `NC_Dualization_Example2_TS` | NC | MMD |
| Heisenberg star graph | `Unipotent_Heisenberg_Star_n8` | Unipotent | 8 sites |

- [ ] trivial example 2 (`USE_LOCAL`, `dualize=false`) — `sos.jl:146`
- [ ] trivial example 2 (`USE_LOCAL`, `dualize=true`) — `sos.jl:147`
- [ ] correlative sparsity (CS, `dualize=true`) — `sos.jl:271`
- [ ] term sparsity (TS, `dualize=true`) — `sos.jl:283`

## API Interface (`test/solvers/interface.jl`)

| Test | Oracle Key | Algebra | Sparsity |
|------|------------|---------|----------|
| naive example | N/A | Pauli | Dense |
| 1D transverse field Ising | N/A | Pauli | - |
| 1D Heisenberg chain | N/A | Pauli | - |
| I_3322 with sparsity (Dense) | `Projector_I3322_Dense_d3` | Projector | Dense |
| I_3322 with sparsity (MF+MMD) | `Projector_I3322_CS_TS_d3` | Projector | MF+MMD |
| Majumdar Gosh model | `Projector_Majumdar_Ghosh` | Projector | Dense |
| README unconstrained (dense) | `NC_README_Unconstrained_Dense` | NC | Dense |
| README unconstrained (CS) | `NC_README_Unconstrained_CS` | NC | MF |
| README unconstrained (CS+TS) | `NC_README_Unconstrained_CS_TS` | NC | MF+MMD |
| README constrained (dense) | `NC_README_Constrained_Dense` | NC | Dense |
| README constrained (CS) | `NC_README_Constrained_CS` | NC | MF |
| README constrained (CS+TS) | `NC_README_Constrained_CS_TS` | NC | MF+MMD |

- [ ] problem creation interface example — `interface.jl:347`

## NC Benchmarks (`test/solvers/ncpop_benchmarks.jl`) — `USE_LOCAL`

| Problem | Variants | Oracle Keys |
|---------|----------|-------------|
| Rosenbrock | NoElimination, MF, MF+MMD | `NC_Rosenbrock_*` |
| Broyden banded | NoElimination, MF, MF+MMD | - |
| Broyden tridiagonal | NoElimination, MF, MF+MMD | - |
| Chained singular | NoElimination, MF, MF+MMD | - |
| Chained wood | NoElimination, MF, MF+MMD | - |

- [ ] sparsity method consistency (loop) — `ncpop_benchmarks.jl:396`

## Sparsity Verification (`test/solvers/sparsity.jl`)

### Bell Inequalities
| Test | Oracle Key | Sparsity | Expected |
|------|------------|----------|----------|
| CHSH (dense) | `CHSH_Dense` | Dense | -2.8284 |
| CHSH (CS only) | `CHSH_CS` | MF | -2.8284 |
| CHSH (TS only) | `CHSH_TS` | MMD | -2.8284 |
| CHSH (CS+TS) | `CHSH_CS_TS` | MF+MMD | -2.8284 |
| I_3322 (dense) | `I3322_Dense_d2` | Dense | -0.25 |
| I_3322 (CS) | `I3322_CS_d2` | MF | -0.25 |
| I_3322 (TS) | `I3322_TS_d2` | MMD | -0.25 |
| I_3322 (CS+TS) | `I3322_CS_TS_d2` | MF+MMD | -0.25 |

### Trace Polynomial
| Test | Oracle Key | Sparsity |
|------|------------|----------|
| CHSH trace (dense) | `Trace_6_2_0` | Dense |
| CHSH trace (MaximalElimination TS) | - | MaxElim |
| covariance trace (dense) | `Trace_6_2_2` | Dense |
| covariance trace (MF+MMD) | - | MF+MMD |

### State Polynomial
| Test | Oracle Key | Sparsity |
|------|------------|----------|
| CHSH state (dense) | `State_7_2_0` | Dense |
| CHSH state (MMD TS) | `State_7_2_0_TS` | MMD |
| covariance state (dense) | `State_7_2_2` | Dense |
| covariance state (MF+MMD) | - | MF+MMD |

### Constrained Problems
| Test | Oracle Key | Sparsity |
|------|------------|----------|
| ball constraint POP (dense) | `NC_Ball_Dense` | Dense |
| ball constraint POP (MMD TS) | `NC_Ball_TS_MMD` | MMD |
| Rosenbrock (dense) | `NC_Rosenbrock_Dense` | Dense |
| Rosenbrock (MF CS) | `NC_Rosenbrock_CS_MF` | MF |
| Rosenbrock (MF+MMD) | `NC_Rosenbrock_CS_TS` | MF+MMD |

## Sparsity Algorithm Variants (`test/solvers/sparsity.jl`)

| Test | Elimination | CS | TS |
|------|-------------|----|----|
| I_3322 (NoElimination) | None | - | - |
| I_3322 (MF) | - | MF | - |
| I_3322 (AsIsElimination) | AsIs | - | - |
| CHSH trace term-sparsity loop | Various | - | Various |

---

## Oracle Problem Consolidation Status

Tracking migration of oracle problems to `test/oracles/problems/` and solver tests to `test/solvers/problems/`.

| Problem | Oracle Definition | Solver Tests | Status |
|---------|-------------------|--------------|--------|
| CHSH | `oracles/problems/chsh.jl` | `solvers/problems/chsh.jl` | ✅ Verified |
| I3322 | `oracles/problems/i3322.jl` | `solvers/problems/i3322.jl` | ✅ Verified |
| Heisenberg Star | `oracles/problems/heisenberg_star.jl` | `solvers/problems/heisenberg_star.jl` | ✅ Verified |
| CS/TS n=10 | `oracles/problems/cs_ts_n10.jl` | `solvers/problems/nc_examples.jl` | ✅ Verified |
| Example 1 | `oracles/problems/example1.jl` | `solvers/problems/nc_examples.jl` | ⬜ Pending |
| Example 2 | `oracles/problems/example2.jl` | `solvers/problems/nc_examples.jl` | ⬜ Pending |
| Correlative Sparsity | `oracles/problems/corr_sparsity.jl` | - | ⬜ Pending |
| Rosenbrock | `oracles/problems/rosenbrock.jl` | - | ⬜ Pending |
| State Polynomial | `oracles/problems/state_poly.jl` | `solvers/problems/state_polynomial.jl` | ⬜ Pending |
| Trace Polynomial | - | `solvers/problems/trace_polynomial.jl` | ⬜ Pending |

---

## Summary Statistics

| Category | Total | Has Oracle | NCTSSoS-only |
|----------|-------|------------|--------------|
| Physics | 14 | 2 | 12 (Pauli/Fermionic) |
| Trace Polynomial | 5 | 5 | 0 |
| State Polynomial | 14 | 6 | 8 |
| Moment Solver | 9 | 9 | 0 |
| SOS Dualization | 11 | 7 | 4 |
| API Interface | 17 | 9 | 8 (Pauli) |
| NC Benchmarks | 16 | 3 | 13 |
| Sparsity | 31 | 15 | 16 |
| **Total** | **117** | **56** | **61** |
