import CSV
import DataFrames

using LinearAlgebra
using SparseArrays

"""Extract the gene co-expression matrix from `file`."""
function load(file::String)::Tuple{SparseMatrixCSC, Vector{Int}}
	# Read file
	df = DataFrames.DataFrame(CSV.File(file, header=["i"; "j"; "value"], delim='\t'))

	# Set of represented genes
	genes = unique([df.i; df.j])
	n = length(genes)

	# Map to lower dimensionality graph
	mapping = Dict(genes .=> 1:n)

	df.i = map(x -> mapping[x], df.i)
	df.j = map(x -> mapping[x], df.j)

	# Co-expression matrix
	S = sparse(df.i, df.j, df.value, n, n)
	S += S'
	S += sparse(Inf * I, n, n)

	return S, genes
end

"""Compute the adjacency matrix by hard-thresholding `S`."""
function adjacency(S::SparseMatrixCSC, tau::Float64 = 0.)::SparseMatrixCSC
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
	n = size(A, 1)
	used = spzeros(Int, n)
	i = 1

	list = Vector{Int}[]

	while !isnothing(i)
		x = sparsevec([i], [1], n)

		while true
			y = A * x
			y.nzval .= 1

			if all(y .== x)
				break
			else
				x = y
			end
		end

		used += x
		i = findnext(iszero, used, i)

		append!(list, [x.nzind])
	end

	return list
end

"""Compute the pivot of `A`.

The pivot is the largest number `p` such that `p`
is smaller than the connectivity of `p` vertices.
"""
function pivot(A::SparseMatrixCSC)::Int
	k = connectivity(A)
	k = sort(k.nzval, rev=true)

	# Dichotomic search
	low, up = 1, floor(Int, sqrt(sum(k)))

	while low < up
		p = (low + up + 1) ÷ 2

		if p > k[p]
			up = p - 1
		else
			low = p
		end
	end

	return low
end

"""Apply the Worst-Out heuristic to `A`."""
function worstout(A::SparseMatrixCSC)::SparseVector
	x = diag(A)
	k = connectivity(A)

	while sum(k) < nnz(x)^2
		i = k.nzind[argmin(k.nzval)]

		x[i] = 0
		dropzeros!(x)

		k = (k .- A[:, i]) .* x
	end

	return x
end
