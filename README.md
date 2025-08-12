# The Coevolution of Cooperation and Competition
*Authors: Angelo Romano, Martin Smit, Carsten K. W. De Dreu, Fernando P. Santos*

This repository contains the code and instructions to reproduce the results and figures from the paper of the same name.

## System requirements

This code was written in [Julia](https://julialang.org/) v1.11.6 and run on a Macbook Air M2 using macOS 14.4.
The instructions below should additionally work with no alterations for newer Julia versions and [other common operating systems and architectures](https://julialang.org/downloads/#supported_platforms).

The full list of dependencies and their versions can be found in the `compat` section of the `Project.toml` but for thoroughness we list them below:
```
julia = "1.11"
CairoMakie = "0.14"
JLD2 = "0.5"
StaticArrays = "1.9"
StatsBase = "0.34"
ThreadsX = "0.1.12"
```

The instructions cover the installation of the correct versions of Julia and the required packages, and details which figures can be reproduced by which scripts.

## Installation guide
The recommended way to install Julia is via [juliaup](https://github.com/JuliaLang/juliaup).
[Installation instructions for juliaup](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#installation) are detailed for both [Windows](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#windows) and [Unix-like systems](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#mac-linux-and-freebsd).
Installation of Julia and compiling the project's dependencies takes 15 minutes on my Macbook Air M2.

Once juliaup is installed, run the following command to install the correct version of Julia and optionally set it as the default Julia version in case another version is already installed
```
juliaup add 1.11
juliaup default 1.11
```

If you don't want to change your defaults, simply replace `julia` with `julia +1.11` in the code snippets below.

The recommended way to use Julia is with a single persistent Julia instance to which lines of code are sent to be compiled and executed.
The most popular way of doing this is with [VSCode](https://code.visualstudio.com/) and the [Julia extension](https://marketplace.visualstudio.com/items?itemName=julialang.language-julia), but support for [Vim](https://github.com/JuliaEditorSupport/julia-vim), [Emacs](https://github.com/JuliaEditorSupport/julia-emacs) and [other popular editors](https://github.com/JuliaEditorSupport) is also available.
That said, the our instructions will *not* assume that these tools are installed and will detail how to run each script from the command line, a slower but more universal method.

With `julia` installed and the corresponding command available from the command line, clone the repository with 
```
git clone https://github.com/jacobusmmsmit/nastyneighbours
```
and navigate to the downloaded folder.
Then run `julia` and once the Julia REPL is open, go into `Pkg` mode by pressing `]`.

Once in `Pkg` mode run ```activate .``` followed by `instantiate` and to install the dependencies at their required versions.
If not using Julia 1.11 and macOS 14 or if encountering any other error, you should delete the `Manifest.toml` file and rebuild it for your Julia version and operating system by, still in `Pkg` mode, running `resolve`.

## Reproducing results

To run a script from the terminal (not the Julia REPL), navigate to the projects home directory and run
```
julia --project "." path/to/script
```

The scripts in `scripts` each reproduce a plot from the main text or supplimentary information:
- `unstructured_results.jl` produces Figure 2
- `structured_results.jl` produces Figure 4
- `ESS_regions.jl` produces Appendix Figure 1
- `analytical_new.jl` produces Appendix Figure 2
- `analytical.jl` produces Appendix Figure 3

To avoid re-running simulations, one can use `@save` and `@load` from the [JLD2.jl package](https://github.com/JuliaIO/JLD2.jl).
An example dataset resulting from a simulation can be found in `data/`, with how to produce and use it in `scripts/analytical.jl`.

To run every simulation and produce every plot from the paper takes around 5 hours on my Macbook. This can be substantially decreased for a loss in heatmap resolution by changing the various `_range` variables such that the argument `length` is smaller.
