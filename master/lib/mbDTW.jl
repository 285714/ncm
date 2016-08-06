module mbDTW
export dtwDist, dtwPath

using mbUtil

function dtwMatrix(A,B,d=(x,y)->norm(x-y))
	m,n = length(A),length(B)
	M = Array{Float64}(m,n)

	M[1,1] = d(A[1],B[1])
	for i in 2:m M[i,1] = M[i-1,1] + d(A[i],B[1]) end
	for j in 2:n M[1,j] = M[1,j-1] + d(A[1],B[j]) end

	for i in 2:m, j in 2:n
		M[i,j] = d(A[i],B[j]) + min(M[i-1, j-1], M[i-1, j], M[i, j-1])
	end

	return M
end

dtwDist = last ∘ dtwMatrix

function dtwPath(A,B)
	M = dtwMatrix(A,B)
	m,n = size(M)

	P = [(m,n)]

	function locmin(M,i,j)
		i≤1 && return (i,j-1)
		j≤1 && return (i-1,j)
		d = [(i-1, j-1), (i, j-1), (i-1, j)]
		rtn = d[indmin(map(i->M[i...], d))]
		return rtn
	end

	while P[1] ≠ (1,1)
		unshift!(P, locmin(M, P[1]...))
	end

	return P
end



end #module
