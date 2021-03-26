include("utils.jl")

# Arguments

filename = ARGS[1]

S, genes = load(filename)
A = adjacency(S)  # tau = 0

list = blocks(A)
list = sort(list, by=length, rev=true)

subset = list[1]
A = A[subset, subset]

println(size(A))

println(pivot(A))

x = worstout(A)
x = x.nzind
x = genes[subset[x]]

println(x)
