using CSV, DataFrames, LinearAlgebra, SparseArrays

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

"""Compute the connectivity vector of `A`."""
function connectivity(A::SparseMatrixCSC)::SparseVector
	return A * diag(A)
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

"""Compute the pivot of `A`.

The pivot is the largest number `p` such that `p`
is smaller than the connectivity of `p` vertices.
"""
function pivot(A::SparseMatrixCSC)::Int
	K = connectivity(A)
	K = sort(K.nzval, rev=true)

	# Dichotomic search
	low, up = 1, floor(Int, sqrt(nnz(A)))

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
	V = Set(diag(A).nzind)
	setdiff!(V, M)

	K = A * x
	vertices = nnz(x)
	edges = x' * K

	K = Vector(K)
	C = Vector(connectivity(A))

	while !isempty(V) && maxiter != 0
		i, (k, _) = findopt(i -> (K[i], C[i]), V; comp=(a, b) -> a > b)

		n = vertices + 1
		m = edges + 2 * k + 1

		if n^2 - m <= 2 * f(n)
			push!(M, i)
			K .+= A[:, i]
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

	K = A * x
	vertices = nnz(x)
	edges = x' * K

	K = Vector(K)

	while vertices^2 - edges > 2 * f(vertices)
		i, k = findopt(i -> K[i], M)

		delete!(M, i)
		K .-= A[:, i]
		vertices -= 1
		edges -= 2 * k - 1
	end

	return sparsevec(collect(M), 1, size(A, 1))
end

"""Simulated-Annealing meta-heuristic"""
function annealing(
	A::SparseMatrixCSC,
	x::Union{SparseVector, Nothing} = nothing;
	f::Function = n -> 0,
	alpha::Float64 = 0.5,
	steps::Int = 1000000
)::SparseVector
	V = diag(A).nzind

	if isnothing(x)
		x = spzeros(Int, size(A, 1))
	end

	K = A * x
	vertices = nnz(x)
	edges = x' * K

	best = copy(x)

	for _ in 1:steps
		i = rand(V)

		remove = Bool(x[i])
		if remove
			n = vertices - 1
			m = edges - 2 * K[i] + 1
		else
			n = vertices + 1
			m = edges + 2 * K[i] + 1
		end

		if n^2 - m <= 2 * f(n)
			if (remove ? rand() < alpha : true)
				x[i] = 1 - x[i]

				if remove
					K -= A[:, i]
				else
					K += A[:, i]
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
