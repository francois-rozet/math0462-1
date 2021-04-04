import CSV, DataFrames

DataFrames.DataFrame(
	CSV.File(
		first(readdir("resources/hsa/", join=true)),
		header=["i"; "j"; "v"],
		delim='\t'
	)
)
