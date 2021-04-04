# Module analysis in cancer diagnosis

Project realized as part of the course *Discrete Optimization* given by **Quentin Louveaux** to graduate computer engineering students at the [University of Liège](https://www.uliege.be/) during the academic year 2020-2021.

## Dependencies

The scripts are implemented using [Julia 1.5](https://julialang.org/) and require several dependencies to be executed. A [`Project.toml`](Project.toml) file is provided to install these dependencies.

```bash
julia --project=.
julia> using Pkg
julia> Pkg.instantiate()
```

The functions from these dependencies (*e.g.* `CSV.File`) take a long time to compile. Creating a *system image* allows to save the
the compiled binaries and reload them later. To do so, create the system image with

```bash
julia --project=.
julia> using PackageCompiler
julia> PackageCompiler.create_sysimage([:CSV, :DataFrames]; sysimage_path="image.so", precompile_execution_file="src/precompile.jl")
```

and load it with the `-J` (or `--sysimage`) flag.

```bash
julia --project=. -J image.so path/to/script.jl
```

## Authors

* **François Rozet** - [francois-rozet](https://github.com/francois-rozet)
* **Maxime Meurisse** - [meurissemax](https://github.com/meurissemax)
