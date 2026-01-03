# Oracle System

Oracles are reference values from the original NCTSSOS package for validation.

## Directory Structure

```
test/oracles/
├── scripts/           # Julia scripts to generate oracles (run in NCTSSOS repo)
│   ├── oracle_utils.jl        # Shared utilities
│   ├── nctssos_chsh.jl        # CHSH Bell inequality
│   ├── nctssos_i3322.jl       # I3322 Bell inequality
│   ├── nctssos_example1.jl    # NC polynomial example 1
│   ├── nctssos_example2.jl    # NC polynomial example 2
│   ├── nctssos_rosenbrock.jl  # Rosenbrock benchmark
│   ├── nctssos_corr_sparsity.jl
│   ├── nctssos_cs_ts_n10.jl
│   ├── nctssos_heisenberg_star.jl
│   ├── nctssos_state_poly.jl
│   ├── nctssos_state_poly_extended.jl
│   └── nctssos_trace_poly.jl
└── results/           # Generated oracle dictionaries
    ├── chsh_oracles.jl
    ├── i3322_oracles.jl
    ├── example1_oracles.jl
    ├── example2_oracles.jl
    ├── rosenbrock_oracles.jl
    ├── corr_sparsity_oracles.jl
    ├── cs_ts_n10_oracles.jl
    ├── heisenberg_star_oracles.jl
    ├── state_poly_oracles.jl
    ├── state_poly_extended_oracles.jl  # Has placeholders
    └── trace_poly_oracles.jl           # Has placeholders
```

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

1. NCTSSOS repository with MosekTools installed
2. Mosek license active

### NCTSSOS Repository Paths

| Machine | Path |
|---------|------|
| Local (macOS) | `/Users/yushengzhao/projects/NCTSSOS` |
| a800 server | `/home/yushengzhao/NCTSSOS` |

Set `NCTSSOS_PATH` environment variable to override.

### Using Makefile (Recommended)

From NCTSSoS.jl repository:
```bash
# Run specific oracle script
make oracle-chsh
make oracle-i3322
make oracle-state_poly

# With custom NCTSSOS path
NCTSSOS_PATH=/custom/path make oracle-chsh
```

### Manual Execution

```bash
# Navigate to NCTSSOS repo
cd /path/to/NCTSSOS

# Run oracle script
julia --project /path/to/NCTSSoS.jl/test/oracles/scripts/nctssos_chsh.jl

# Copy output to results file
# Output format: const CHSH_ORACLES = Dict(...)
```

## Placeholder Values

The following files contain placeholder values (marked with `# PLACEHOLDER`):

### `state_poly_extended_oracles.jl`
- `State_7_2_1_Dense_d3`, `State_7_2_1_TS_d3`
- `State_7_2_2_Dense_d2`, `State_7_2_2_TS_d2`
- `State_7_2_3_Dense_d2`, `State_7_2_3_TS_d2`

### `trace_poly_oracles.jl`
- `Trace_6_1_Dense_d2`, `Trace_6_1_Dense_d3`, `Trace_6_1_TS_d2`, `Trace_6_1_TS_d3`
- `Trace_6_2_0_Dense_d1`, `Trace_6_2_0_TS_d1`
- `Trace_6_2_1_Dense_d2`, `Trace_6_2_1_TS_d2`
- `Trace_6_2_2_Dense_d2`, `Trace_6_2_2_TS_d2`

Placeholder format: `sides=[Int[]], nuniq=0`

Run the corresponding oracle scripts to generate actual values.

## Available Oracle Scripts

| Script | Problem | Notes |
|--------|---------|-------|
| `nctssos_chsh.jl` | CHSH Bell inequality | 3 formulations (NC, State, Trace) |
| `nctssos_i3322.jl` | I3322 Bell inequality | Multiple orders |
| `nctssos_example1.jl` | 3-variable NC polynomial | Dense, TS |
| `nctssos_example2.jl` | Rosenbrock-like | Dense, TS |
| `nctssos_rosenbrock.jl` | Rosenbrock optimization | Dense, TS |
| `nctssos_corr_sparsity.jl` | Correlative sparsity | CS, TS |
| `nctssos_cs_ts_n10.jl` | Large-scale n=10 | Multiple combinations |
| `nctssos_heisenberg_star.jl` | Condensed matter | Multiple sizes |
| `nctssos_state_poly.jl` | State polynomial 7.2.0 | Pure state expectations |
| `nctssos_state_poly_extended.jl` | State polynomial 7.2.1-7.2.3 | Products |
| `nctssos_trace_poly.jl` | Trace polynomial 6.x | Examples 6.1, 6.2.0-6.2.2 |
