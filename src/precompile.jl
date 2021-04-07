using ArgParse, CSV, DataFrames

# ArgParse
parser = ArgParseSettings()

@add_arg_table parser begin
	"--input", "-i"
	help = "input file"
end

parse_args(parser)

# CSV, DataFrames
DataFrame(
	CSV.File(
		first(readdir("resources/hsa/", join=true)),
		header=["i"; "j"; "v"],
		delim='\t'
	)
)
