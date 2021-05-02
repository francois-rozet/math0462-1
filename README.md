# Module analysis in cancer diagnosis

Project realized as part of the course *Discrete Optimization* given by **Quentin Louveaux** to graduate computer engineering students at the [University of Liège](https://www.uliege.be/) during the academic year 2020-2021.

## Dependencies

The scripts are implemented using [Julia 1.6](https://julialang.org/) and require several dependencies to be executed. A [`Project.toml`](Project.toml) file is provided to install these dependencies.

```bash
julia --project=.
julia> using Pkg
julia> Pkg.instantiate()
julia> Pkg.build()
```

The functions from these dependencies (*e.g.* `CSV.File`) take a long time to compile. Creating a *system image* allows to save the compiled binaries and reload them later. To do so, create the system image with

```bash
julia --project=.
julia> using PackageCompiler
julia> deps = [:ArgParse, :CSV, :DataFrames, :Gurobi, :JuMP]
julia> PackageCompiler.create_sysimage(deps; sysimage_path="image.so", precompile_execution_file="src/precompile.jl")
```

and load it with the `-J` (or `--sysimage`) flag.

```bash
julia --project=. -J image.so path/to/script.jl
```

## Usage

The central file of this project is [`script.jl`](src/script.jl). It can be executed on any genomic co-expression file from [`resources/hsa/`](resources/hsa/).

```bash
julia --project=. -J image.so src/script.jl resources/hsa/BGSE1456.txt
```

The execution will display (on standard output) a few interesting metrics evaluated on the graph represented by the provided file. These metrics include, among others, the number of independent complete sub-graphs and, according to several algorithms, the size of the largest module in the largest sub-graph. It will also create a file (by default `out.json`) containing a covering of the genes.

### Parameters

A few parameters can be passed on to `script.jl` using flags. For example, to use a constant (*e.g.* `4.2069`) tolerance function,

```bash
julia --project=. -J image.so src/script.jl resources/hsa/BGSE1456.txt --tol 4.2069
```

## Authors

* **François Rozet** - [francois-rozet](https://github.com/francois-rozet)
* **Maxime Meurisse** - [meurissemax](https://github.com/meurissemax)
