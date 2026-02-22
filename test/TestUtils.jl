using JuMP, COSMO

const SOLVER = optimizer_with_attributes(
    COSMO.Optimizer,
    "verbose" => false,
    "eps_abs" => 1e-7,
    "eps_rel" => 1e-7,
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
