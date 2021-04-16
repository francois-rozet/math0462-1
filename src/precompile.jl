# ArgParse
using ArgParse

parser = ArgParseSettings()

@add_arg_table parser begin
	"--input", "-i"
	help = "input file"
end

parse_args(parser)

# CSV, DataFrames
using CSV, DataFrames

DataFrame(
	CSV.File(
		first(readdir("resources/hsa/", join=true)),
		header=["i"; "j"; "v"],
		delim='\t'
	)
)

# Gurobi, JuMP
using Gurobi, JuMP

model = Model(Gurobi.Optimizer)

A = [0 1 0; 1 0 0; 0 0 0]

@variable(model, x[1:3], Bin)
@constraint(model, A * x .<= 3 * (1 .- x))
@objective(model, Max, sum(x))

optimize!(model)
