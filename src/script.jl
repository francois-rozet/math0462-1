include("utils.jl")

using ArgParse
using Statistics

# Arguments
function parse()
	parser = ArgParseSettings()

	@add_arg_table parser begin
		"input"
			help = "input file (.txt)"
			required = true
		"--output", "-o"
			help = "output file (.json)"
			default = nothing
		"--tol"
			help = "tolerance function"
			default = "0"
		"--overlap"
			help = "overlapping coverage"
			default = "none"
	end

	return parse_args(ARGS, parser)
end

args = parse()

## Output
if isnothing(args["output"])
	output = replace(args["input"], ".txt" => ".json")
else
	output = args["output"]
end

## Tolerance function
C = tryparse(Float64, args["tol"])

alpha = 0.
beta = 0.
gamma = 0.

if isnothing(C)
	if args["tol"] == "lin"
		tol = n -> n / 2
		beta = 1.
	elseif args["tol"] == "log"
		tol = n -> log1p(n) * (n - 1) / 2
	elseif args["tol"] == "root"
		tol = n -> sqrt(n) * (n - 1) / 2
	else  # half
		tol = n -> n * (n - 1) / 4
		alpha = .5
		beta = -alpha
	end
else
	tol = n -> C
	gamma = 2 * C
end

## Overlap
if args["overlap"] == "limited"
	over = n -> ceil(Int, log1p(n))
elseif args["overlap"] == "full"
	over = n -> -1
else  # == "none"
	over = n -> 0
end

# Co-Expression Matrix
S, genes = load(args["input"])

println("Stats")
println("≡≡≡≡≡")

println("Number of genes: $(vertices(S))")
println("Number of co-expressions: $(edges(S))")

# Adjacency Matrix
A = adjacency(S)
K = Vector(degree(A))

println("Mean degree: $(mean(K))")
println("Median degree: $(median(K))")
println("Max degree: $(maximum(K))")
println()

# Blocks

println("Blocks")
println("≡≡≡≡≡≡")

list = @time blocks(A)

sort!(list, by=length, rev=true)

println("Number of blocks: $(length(list))")

lengths = map(length, list)

println("Mean block size: $(mean(lengths))")
println("Median block size: $(median(lengths))")
println("Max block size: $(maximum(lengths))")
println()

for i in 1:1
	## Analysis
	println("Block $i")
	println("=======")

	idx = list[i]
	B = A[idx, idx]

	println("|V| = $(vertices(B))")
	println("|E| = $(edges(B))")
	println("√G = $(root(B))")
	println("pivot(G) = $(pivot(B))")

	@time pivot(B)
	println()

	for (name, routine) in [
		("Solver", gurobisolver),
		("Best-In", bestin),
		("Worst-Out", worstout),
		("Simulated annealing", annealing),
		("Warm-Start annealing", warmstart)
	]
		println(name)
		println(repeat('-', length(name)))

		if name == "Solver"
			@time x = routine(B, alpha=alpha, beta=beta, gamma=gamma)
		else
			@time x = routine(B, f=tol)
		end

		println("|M| = $(nnz(x))")
		println("𝛿(M) = $(delta(B, x))")
		println()
	end
end

# Covering
println("Covering")
println("≡≡≡≡≡≡≡≡")

modules = Vector{Int}[]

@time for idx in list
	B = A[idx, idx]

	temp = ice(B, f=tol, overlap=over)
	temp = map(x -> idx[x], temp)

	append!(modules, temp)
end

sort!(modules, by=length, rev=true)

println("Number of modules: $(length(modules))")

## Sizes
lengths = map(length, modules)

println("Mean module size: $(mean(lengths))")
println("Median module size: $(median(lengths))")
println("Max module size: $(maximum(lengths))")
println()

## Overlap
cnt = zeros(Int, size(A, 1))

for M in modules
	cnt[M] .+= 1
end

println("Mean vertex coverage: $(mean(cnt))")
println("Median vertex coverage: $(median(cnt))")
println("Max vertex coverage: $(maximum(cnt))")

## Export
open(output, "w") do io
	for M in modules
		println(io, sort(genes[M]))
	end
end
