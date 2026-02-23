# Oracle System

Oracles are reference values from the original NCTSSOS package for validation.

## Directory Structure

```
test/oracles/
├── oracle_utils.jl              # Shared utilities for oracle generation
├── nctssos_bell_template.jl     # Bell inequality template
├── nctssos_benchmarks.jl        # NCPOP benchmarks
├── nctssos_bilocal.jl           # Bilocal network
├── nctssos_chsh.jl              # CHSH Bell inequality
├── nctssos_corr_sparsity.jl     # Correlative sparsity examples
├── nctssos_cs_ts_n10.jl         # Large-scale n=10 problems
├── nctssos_example1.jl          # NC polynomial example 1
├── nctssos_example2.jl          # NC polynomial example 2
├── nctssos_heisenberg_star.jl   # Heisenberg star model
├── nctssos_i3322.jl             # I3322 Bell inequality
├── nctssos_rosenbrock.jl        # Rosenbrock optimization
├── nctssos_state_poly.jl        # State polynomial examples
├── nctssos_state_poly_extended.jl # Extended state polynomial
├── nctssos_trace_poly.jl        # Trace polynomial examples
└── README.md
```

## Oracle Script → Test File Mapping

| Oracle Script | Updates Test File(s) |
|---------------|---------------------|
| `nctssos_chsh.jl` | `test/problems/bell_inequalities/chsh_simple.jl`, `chsh_state.jl`, `chsh_trace.jl`, `chsh_high_order.jl` |
| `nctssos_i3322.jl` | `test/problems/bell_inequalities/i3322.jl` |
| `nctssos_bell_template.jl` | `test/problems/bell_inequalities/bell_inequalities.jl` |
| `nctssos_bilocal.jl` | `test/problems/quantum_networks/bilocal_networks.jl` |
| `nctssos_example1.jl` | `test/problems/nc_polynomial/nc_example1.jl` |
| `nctssos_example2.jl` | `test/problems/nc_polynomial/nc_example2.jl` |
| `nctssos_corr_sparsity.jl` | `test/correlated_sparsity/core_pipeline_numeric.jl` |
| `nctssos_cs_ts_n10.jl` | `test/problems/nc_polynomial/nc_large_scale.jl` |
| `nctssos_rosenbrock.jl` | `test/problems/benchmarks/ncpop_benchmarks.jl` |
| `nctssos_benchmarks.jl` | `test/problems/benchmarks/ncpop_benchmarks.jl` |
| `nctssos_heisenberg_star.jl` | `test/problems/condensed_matter/heisenberg_star.jl` |
| `nctssos_state_poly.jl` | `test/problems/state_polynomial/state_polynomial.jl` |
| `nctssos_state_poly_extended.jl` | `test/problems/state_polynomial/state_polynomial.jl` |
| `nctssos_trace_poly.jl` | `test/problems/trace_polynomial/trace_polynomial.jl` |

## Oracle Format

Each oracle entry is a named tuple:
```julia
(opt, sides, nuniq)
```

| Field | Type | Description |
|-------|------|-------------|
| `opt` | Float64 | Optimal objective value (minimization) |
| `sides` | Vector{Int} | Moment matrix block sizes |
| `nuniq` | Int | Unique moment indices (affine constraints) |

Example:
```julia
const CHSH_ORACLES = Dict(
    "CHSH_Dense_d1" => (opt=-2.8284271321623193, sides=[5], nuniq=11),
    "CHSH_CS_d1" => (opt=-2.8284271247170496, sides=[4, 4], nuniq=10),
)
```

## Regenerating Oracles

### Prerequisites

1. NCTSSOS repository
2. Solver installed:
   - **Mosek** (default): MosekTools + valid license
   - **COSMO** (open-source): COSMO.jl package

### NCTSSOS Repository Paths

| Machine | Path |
|---------|------|
| Local (macOS) | `/Users/yushengzhao/projects/NCTSSOS` |
| a800 server | `/home/yushengzhao/NCTSSOS` |

Set `NCTSSOS_PATH` environment variable to override.

### Using Makefile (Recommended)

From NCTSSoS.jl repository:
```bash
# Run with Mosek (default)
make oracle-chsh
make oracle-i3322

# Run with COSMO solver
make oracle-chsh SOLVER=cosmo
make oracle-i3322 SOLVER=cosmo

# With custom NCTSSOS path
NCTSSOS_PATH=/custom/path make oracle-chsh SOLVER=mosek
```

After running, paste the output to Claude to update the test file.

### Manual Execution

```bash
# Navigate to NCTSSOS repo
cd /path/to/NCTSSOS

# Run with Mosek (default)
julia --project /path/to/NCTSSoS.jl/test/oracles/nctssos_chsh.jl

# Run with COSMO
ORACLE_SOLVER=cosmo julia --project /path/to/NCTSSoS.jl/test/oracles/nctssos_chsh.jl
```
