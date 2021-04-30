using CSV, DataFrames, Gurobi, JuMP, LinearAlgebra, SparseArrays

"""Extract the gene co-expression matrix from `file`."""
function load(file::String)::Tuple{SparseMatrixCSC, Vector{Int}}
	# Read file
	df = DataFrame(CSV.File(file, header=["i"; "j"; "v"], delim='\t'))

	# Set of represented genes
	genes = unique([df.i; df.j])
	n = length(genes)

	# Map to lower dimensionality graph
	mapping = Dict(genes .=> 1:n)

	df.i = map(x -> mapping[x], df.i)
	df.j = map(x -> mapping[x], df.j)

	# Co-expression matrix
	S = sparse(df.i, df.j, df.v, n, n)
	S += S'
	S += sparse(Inf * I, n, n)

	return S, genes
end

"""Compute the adjacency matrix by hard-thresholding `S`."""
function adjacency(S::SparseMatrixCSC; tau::Float64 = 0.)::SparseMatrixCSC
	return convert.(Int, S .> tau)
end

"""Compute the degree vector of `A`."""
function degree(A::SparseMatrixCSC)::SparseVector
	return A * diag(A)
end

"""Compute the number of vertices of `A`, i.e. |V|."""
function vertices(A::SparseMatrixCSC)::Int
	return nnz(diag(A))
end

"""Compute the number of edges of `A`, i.e. |E|."""
function edges(A::SparseMatrixCSC)::Int
	return (nnz(A) + vertices(A)) ÷ 2
end

"""Compute the number of missing edges in `A` for `x` to be a module, i.e. 𝛿(x)."""
function delta(A::SparseMatrixCSC, x::SparseVector)::Int
	return (nnz(x)^2 - x' * A * x) ÷ 2
end

"""Compute the root of `A`, i.e.
	(1 + √(1 + 8 (|E| - |V|))) / 2 .
"""
function root(A::SparseMatrixCSC)::Int
	rt = (1 + sqrt(1 + 8 * (edges(A) - vertices(A)))) / 2
	return floor(Int, rt)
end

"""Compute the pivot of `A`.

The pivot is the largest number `p` such that `p`
is smaller than the degree of `p` vertices.
"""
function pivot(A::SparseMatrixCSC)::Int
	K = degree(A)
	K = sort(K.nzval, rev=true)

	# Dichotomic search
	low, up = 1, root(A)

	while low < up
		p = (low + up + 1) ÷ 2

		if p > K[p]
			up = p - 1
		else
			low = p
		end
	end

	return low
end

"""Find all block sub-matrices of `A`.

A block sub-matrix corresponds to a connected sub-graph.
"""
function blocks(A::SparseMatrixCSC)::Vector{Vector{Int}}
	V = Set(A.rowval)
	list = Vector{Int}[]

	while !isempty(V)
		i = first(V)

		x = sparsevec([i], 1, size(A, 1))

		while true
			y = A * x

			if nnz(y) == nnz(x)
				break
			end

			x = y
		end

		setdiff!(V, x.nzind)
		append!(list, [x.nzind])
	end

	return list
end

"""Gurobi solver"""
function gurobisolver(
	A::SparseMatrixCSC;
	alpha::Float64 = 0.,
	beta::Float64 = 0.,
	gamma::Float64 = 0.,
	maxtime::Float64 = 60.
)::SparseVector
	model = Model(Gurobi.Optimizer)
	set_silent(model)
	set_time_limit_sec(model, maxtime)

	n = size(A, 1)
	K = Vector(degree(A))
	A = Array(A)

	@variable(model, x[1:n], Bin)

	if alpha == beta == gamma == 0.
		@constraint(model, [i=1:n, j=1:i-1], (1 - A[i, j]) * (x[i] + x[j] - 1) <= 0)
	else
		@variable(model, y[1:n])
		@constraint(model, sum(y) <= beta * sum(x) + gamma)
		@constraint(model, (1 - alpha .- A) * x .- (1 - alpha) * (n .- K) .* (1 .- x) .<= y)
		@constraint(model, -alpha * K .* x .<= y)
	end

	@objective(model, Max, sum(x))

	optimize!(model)

	println(raw_status(model))

	return sparse(round.(Int, value.(x)))
end

"""Return the element of `itr` whose value in `f` is optimal."""
function findopt(f::Function, itr; comp::Function = (a, b) -> a < b)
	x = first(itr)
	fx = f(x)

	for y in itr
		fy = f(y)
		if comp(fy, fx)
			x, fx = y, fy
		end
	end

	return x, fx
end

"""Best-In heuristic"""
function bestin(
	A::SparseMatrixCSC,
	x::Union{SparseVector, Nothing} = nothing;
	f::Function = n -> 0,
	maxiter::Int = -1
)::SparseVector
	if isnothing(x)
		x = spzeros(Int, size(A, 1))
	end

	M = Set(x.nzind)

	V = setdiff(Set(diag(A).nzind), M)
	K = Vector(degree(A))

	d = Vector(A * x)
	vertices = nnz(x)
	edges = x' * d

	while !isempty(V) && maxiter != 0
		i, (d_i, _) = findopt(i -> (d[i], K[i]), V; comp=(a, b) -> a > b)

		n = vertices + 1
		m = edges + 2 * d_i + 1

		if n^2 - m <= 2 * f(n)
			push!(M, i)
			d .+= A[:, i]
			vertices = n
			edges = m
		else
			break
		end

		delete!(V, i)
		maxiter -= 1
	end

	return sparsevec(collect(M), 1, size(A, 1))
end

"""Worst-Out heuristic"""
function worstout(
	A::SparseMatrixCSC,
	x::Union{SparseVector, Nothing} = nothing;
	f::Function = n -> 0
)::SparseVector
	if isnothing(x)
		x = diag(A)
	end

	M = Set(x.nzind)

	d = Vector(A * x)
	vertices = nnz(x)
	edges = x' * d

	while vertices^2 - edges > 2 * f(vertices)
		i, d_i = findopt(i -> d[i], M)

		delete!(M, i)
		d .-= A[:, i]
		vertices -= 1
		edges -= 2 * d_i - 1
	end

	return sparsevec(collect(M), 1, size(A, 1))
end

"""Simulated annealing meta-heuristic"""
function annealing(
	A::SparseMatrixCSC,
	x::Union{SparseVector, Nothing} = nothing;
	f::Function = n -> 0,
	alpha::Function = t -> 4 + 0.5 * log10(t),
	steps::Int = 1000000
)::SparseVector
	if isnothing(x)
		x = spzeros(Int, size(A, 1))
	end

	V = diag(A).nzind
	K = Vector(degree(A))
	delta = 1 .+ K / maximum(K)

	d = A * x
	vertices = nnz(x)
	edges = x' * d

	best = copy(x)

	for t in 1:steps
		i = rand(V)

		if Bool(x[i])
			n = vertices - 1
			m = edges - 2 * d[i] + 1
		else
			n = vertices + 1
			m = edges + 2 * d[i] + 1
		end

		if n^2 - m <= 2 * f(n)
			if (n < vertices ? rand() < alpha(t)^-delta[i] : true)
				x[i] = 1 - x[i]

				if n < vertices
					d -= A[:, i]
				else
					d += A[:, i]
				end

				vertices = n
				edges = m

				if vertices > nnz(best)
					best = dropzeros(x)
				end
			end
		end
	end

	return best
end

"""Warm-Start simulated annealing"""
function warmstart(
	A::SparseMatrixCSC;
	f::Function = n -> 0,
)::SparseVector
	x = worstout(A, f=f)
	return annealing(A, x, f=f)
end

"""Iterative clique enumeration heuristic"""
function ice(A::SparseMatrixCSC; f::Function = n -> 0, overlap::Function = n -> 0)::Vector{Vector{Int}}
	x = diag(A)

	list = Vector{Int}[]

	while nnz(x) > 0
		y = worstout(A, x, f=f)
		x -= y
		dropzeros!(x)
		y = bestin(A, y, f=f, maxiter=overlap(nnz(y)))
		append!(list, [y.nzind])
	end

	return list
end
