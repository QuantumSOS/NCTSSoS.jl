using JuMP, COSMO

const SOLVER = optimizer_with_attributes(
    COSMO.Optimizer,
    "verbose" => false,
    "eps_abs" => 1e-8,
    "eps_rel" => 1e-8,
    "eps_prim_inf" => 1e-6,
    "eps_dual_inf" => 1e-6,
    "max_iter" => 50_000,
    "rho" => 1.0,
    "adaptive_rho" => true,
    "alpha" => 1.0,
    "scaling" => 10,
)

"""
    flatten_sizes(sizes)

Flatten nested moment matrix sizes for comparison with oracle values.

# Example
```julia
flatten_sizes([[3, 3], [2]]) == [3, 3, 2]
```
"""
flatten_sizes(sizes) = reduce(vcat, sizes)
