# NCTSSoS

[NCTSSoS.jl](https://github.com/wangjie212/NCTSSoS.jl) aims to provide a *user-friendly* and *efficient* tool for solving [noncommutative polynomial optimization problems](@ref polynomial-optimization) which is based on the structured [moment-SOHS hierarchy](@ref moment-sohs-hierarchy).

## Features

- Ergonomic API: Easy to use and intuitive interface for defining polynomial optimization problems!
- General Objectives: [Eigenvalue](@ref eigenvalue-polynomial), [State](@ref state-polynomial), [Tracial](@ref tracial-polynomial) Polynomial optimizations are supported!
- [Correlative and Term Sparsities](@ref sparsities): Plug-and-Play API for utilizing sparsities in reducing cost of solving optimization problems.

!!! warning "Performance Issue"
	Version 0.1.0 has performance issues! We plan to release a new version with better performance in the near future. Stay tuned!

## Installation

[NCTSSoS.jl](https://github.com/wangjie212/NCTSSoS.jl) could be installed by running

```julia
using Pkg
Pkg.add("NCTSSoS")
```
To make sure everything works, you can execute the test suite.
See `TESTING.md` in the repository root.

Can't wait to try NCTSSoS.jl in action? Hop to [Quick Start](@ref quick-start).

## Credits

The following people are involved in the development of NCTSSoS.jl:

- [Jie Wang](https://wangjie212.github.io/jiewang), Academy of Mathematics and Systems Science, Chinese Academy of Sciences.
- [Jin-Guo Liu](https://giggleliu.github.io/), Advanced Materials Thrust , The Hong Kong University of Science and Technology(Guangzhou).
- [Yusheng Zhao](https://exaclior.github.io/), Advanced Materials Thrust , The Hong Kong University of Science and Technology(Guangzhou).
- [Huanhai Zhou](https://github.com/fliingelephant), Advanced Materials Thrust , The Hong Kong University of Science and Technology(Guangzhou).

If this project is useful for your work please consider

- Citing NCTSSoS.jl
- Leave us a star on [GitHub!](https://github.com/wangjie212/NCTSSoS.jl)!

## License

NCTSSoS.jl is licensed under the MIT License.

# Related packages

- [TSSOS](https://github.com/wangjie212/TSSOS): Commutative polynomial optimization
- [ChordalGraph](https://github.com/wangjie212/ChordalGraph): Chordal graphs and chordal extentions
