include("utils.jl")

# Arguments
filename = ARGS[1]

if length(ARGS) > 1
	relax = tryparse(Int, ARGS[2])

	if isnothing(relax)
		if ARGS[2] == "log"
			f = n -> floor(Int, log1p(n))
		else  # == "linear"
			f = n -> n
		end
	else
		f = n -> relax
	end
else
	f = n -> 0
end

# Co-Expression Matrix
S, genes = load(filename)

println("Stats")
println("≡≡≡≡≡")

println("Number of genes: $(length(genes))")
println("Number of co-expressions: $(nnz(S))")

# Adjacency Matrix
A = adjacency(S)
list = blocks(A)
list = sort(list, by=length, rev=true)

println("Number of blocks: $(length(list))")

@time blocks(A)
println()

for i in 1:2
	## Analysis
	println("Block $i")
	println("=======")

	idx = list[i]
	x = sparsevec(idx, 1, size(A, 1))
	B = x' .* A .* x

	println("|V| = $(nnz(diag(B)))")
	println("|E| = $(nnz(B))")
	println("√|E| = $(floor(Int, sqrt(nnz(B))))")
	println("pivot(G) = $(pivot(B))")

	@time pivot(B)
	println()

	## Heuristics
	for (name, routine) in [
		("Best-In", bestin),
		("Worst-Out", worstout),
		("Simulated-Annealing", annealing)
	]
		x = routine(B, tolerance=f)

		println(name)
		println(repeat('-', length(name)))
		println("|M| = $(nnz(x))")
		println("|E ∩ M×M| = $(x' * B * x)")

		@time routine(B, tolerance=f)
		println()
	end
end
