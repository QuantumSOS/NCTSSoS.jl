<div align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/dark_logo.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/src/assets/light_logo.png">
  <img alt="NCTSSoS Logo">
</picture>
</div>

---
[![][docs-dev-img]][docs-dev-url]
[![CI][main-ci-img]][main-ci-url]
[![codecov][codecov-img]][codecov-url]
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

[NCTSSoS.jl](https://github.com/wangjie212/NCTSSoS) aims to provide a *efficient* tool for solving sparse noncommutative polynomial optimization problems which is based on the structured moment-SOHS hierarchy.

It is a successor to [NCTSSOS](https://github.com/wangjie212/NCTSSOS).

## Installation

<p>
NCTSSoS is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install NCTSSoS,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, and then type the following command:
</p>

For stable release:

```julia
pkg> add NCTSSoS 
```

For current master:

```julia
pkg> add NCTSSoS#master
```
[main-ci-img]: https://github.com/wangjie212/NCTSSoS.jl/actions/workflows/CI.yml/badge.svg
[main-ci-url]: https://github.com/wangjie212/NCTSSoS.jl/actions/workflows/CI.yml

[codecov-img]: https://codecov.io/gh/wangjie212/NCTSSoS.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/wangjie212/NCTSSoS.jl

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://wangjie212.github.io/NCTSSoS.jl/stable

## Supporting and Citing

Much of the software in this ecosystem was developed as part of academic research. If you would like to help support it, please star the repository as such metrics may help us secure funding in the future. If you use our software as part of your research, teaching, or other activities, we would be grateful if you could cite our work. The [CITATION.bib](CITATION.bib) file in the root of this repository lists the relevant papers.

## Contact
[Jie Wang](https://wangjie212.github.io/jiewang/): wangjie212@amss.ac.cn
